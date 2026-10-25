//! `GeomFill_CircularBlendFunc` — the circular fillet section law that
//! `Approx_SweepApproximation` drives.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_CircularBlendFunc.cxx:28-677` and
//! `GeomFill_CircularBlendFunc.hxx:36-155` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Adaptor3d_Curve` is [`Curve`]; `A->Trim(First, Last, Eps)` becomes
//!   [`GeomTrimmedCurveBasis`] over the same curve (`GeomAdaptor_Curve::Trim`
//!   keeps the basis parameterisation, `GeomAdaptor_Curve.cxx:567-572`).
//! - `SetInterval` mutates `myTPath` / `myTCurve1` / `myTCurve2` through a
//!   `const` method in OCCT; the port stores the trimmed trio behind a `Mutex`
//!   like the other `ApproxSweepFunction` ports in this crate.
//! - `GeomFill::GetShape` / `Knots` / `Mults` / `GetTolerance` /
//!   `GetMinimalWeights` / `GetCircle` are [`crate::geom_fill::section_convert`].
//! - `GeomFillFusInt` is [`crate::geom_lib::fuse_intervals`] with
//!   `Precision::PConfusion() * 0.99` and no re-anchoring to the first interval.
//! - `GeomFill::GetTolerance` can raise out of `GeomConvert`; the trait method
//!   cannot, so the failure falls back to `SurfTol` (see [`Self::get_tolerance`]).

use std::sync::{Arc, Mutex};

use occt_core::approx::{ApproxSweepFunction, SweepD1, SweepD2};
use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpXyz};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use crate::curve::Curve;
use crate::gcpnts::QuasiUniformDeflection;
use crate::geom_fill::section_convert::{
    get_circle, get_circle_d1, get_circle_d2, get_minimal_weights, get_shape, get_tolerance, knots,
    mults,
};
use crate::geom_lib::fuse_intervals;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `static const double TolAng = 1.e-6` (`cxx:29`).
const TOL_ANG: f64 = 1.0e-6;

/// `GeomFillNextShape(S)` (`cxx:31-48`): only the plain `C*` chain climbs, and
/// `CN` is its fixed point. `GeomAbs_Shape` keeps `G1`/`G2` between the `C*`
/// values (`GeomAbs_Shape` = 0..6), so the mapping goes through the enum.
fn next_shape(s: u8) -> u8 {
    match s {
        0 => Shape::C1 as u8,
        2 => Shape::C2 as u8,
        4 => Shape::C3 as u8,
        _ => Shape::CN as u8,
    }
}

/// The trio `SetInterval` swaps in (`cxx:635-642`).
type TrimmedCurves = (Arc<dyn Curve>, Arc<dyn Curve>, Arc<dyn Curve>);

/// `GeomFill_CircularBlendFunc` (`GeomFill_CircularBlendFunc.hxx:36-155`).
pub struct CircularBlendFunc {
    /// `myBary` (`hxx:141`).
    bary: GpPnt,
    /// `myRadius` (`hxx:142`).
    radius: f64,
    /// `maxang` (`hxx:143`).
    max_ang: f64,
    /// `minang` (`hxx:144`).
    min_ang: f64,
    /// `distmin` (`hxx:145`). Computed by `Discret` but, like OCCT, never read
    /// afterwards.
    #[allow(dead_code)]
    dist_min: f64,
    /// `myPath` (`hxx:146`).
    path: Arc<dyn Curve>,
    /// `myCurve1` (`hxx:147`).
    curve1: Arc<dyn Curve>,
    /// `myCurve2` (`hxx:148`).
    curve2: Arc<dyn Curve>,
    /// `myTPath` / `myTCurve1` / `myTCurve2` (`hxx:149-151`).
    trimmed: Mutex<TrimmedCurves>,
    /// `myDegree` (`hxx:152`).
    degree: i32,
    /// `myNbKnots` (`hxx:153`).
    nb_knots: i32,
    /// `myNbPoles` (`hxx:154`).
    nb_poles: i32,
    /// `myTConv` (`hxx:155`).
    tconv: ParameterisationType,
    /// `myreverse` (`hxx:156`). `Discret` always clears it (`cxx:264-266`).
    reverse: bool,
}

impl CircularBlendFunc {
    /// `GeomFill_CircularBlendFunc(Path, Curve1, Curve2, Radius, Polynomial)`
    /// (`cxx:108-141`).
    pub fn new(
        path: Arc<dyn Curve>,
        curve1: Arc<dyn Curve>,
        curve2: Arc<dyn Curve>,
        radius: f64,
        polynomial: bool,
    ) -> Self {
        let dis = Self::discret(&path, &curve1, &curve2);
        let mut tconv = if polynomial {
            ParameterisationType::Polynomial
        } else if dis.max_ang > 0.65 * std::f64::consts::PI {
            // Continuous is preferable (`cxx:133`).
            ParameterisationType::QuasiAngular
        } else {
            ParameterisationType::TgtThetaOver2
        };
        let (mut nb_poles, mut nb_knots, mut degree) = (0, 0, 0);
        get_shape(dis.max_ang, &mut nb_poles, &mut nb_knots, &mut degree, &mut tconv);
        Self {
            bary: dis.bary,
            radius,
            max_ang: dis.max_ang,
            min_ang: dis.min_ang,
            dist_min: dis.dist_min,
            trimmed: Mutex::new((path.clone(), curve1.clone(), curve2.clone())),
            path,
            curve1,
            curve2,
            degree,
            nb_knots,
            nb_poles,
            tconv,
            // `myreverse = false` (`cxx:266`).
            reverse: false,
        }
    }

    /// `Discret` (`cxx:146-266`): numeric estimates of the sweep angle
    /// (`minang` / `maxang`), of the minimal curve distance and of the surface
    /// barycentre.
    fn discret(path: &Arc<dyn Curve>, curve1: &Arc<dyn Curve>, curve2: &Arc<dyn Curve>) -> Discret {
        let t_first = path.first_parameter();
        let t_last = path.last_parameter();
        let mid = (t_first + t_last) / 2.0;

        let p1 = curve1.value(t_first);
        let p2 = curve1.value(mid);
        let p3 = curve1.value(t_last);
        let l1 = p1.distance(&p2) + p2.distance(&p3);

        let p1 = curve2.value(t_first);
        let p2 = curve2.value(mid);
        let p3 = curve2.value(t_last);
        let l2 = p1.distance(&p2) + p2.distance(&p3);

        let (l, c) = if l1 > l2 {
            (l1, curve1.clone())
        } else {
            (l2, curve2.clone())
        };

        let fleche = 1.0e-2 * l;
        let mut samp = QuasiUniformDeflection::new();
        samp.initialize_whole(c.as_ref(), fleche, Shape::C0 as u8);

        // `myBary` starts at the origin (`cxx:179`).
        let mut dis = Discret {
            max_ang: f64::NEG_INFINITY,
            min_ang: f64::INFINITY,
            dist_min: f64::INFINITY,
            bary: GpPnt::zero(),
        };
        let mut bary = GpXyz::zero();

        if samp.is_done() {
            let percent = 1.0 / (2 * samp.nb_points()) as f64;
            for ii in 1..=samp.nb_points() {
                let Some(t) = samp.parameter(ii) else {
                    continue;
                };
                let p1 = curve1.value(t);
                let p2 = curve2.value(t);
                let center = path.value(t);
                let mut ns1 = GpVec::from_pnts(&p1, &center);
                let mut ns2 = GpVec::from_pnts(&p2, &center);
                ns1.normalize();
                ns2.normalize();
                dis.accumulate(&ns1, &ns2, &p1, &p2, &mut bary);
            }
            bary = bary.multiplied(percent);
        } else {
            let delta = (t_last - t_first) / 20.0;
            let percent = 1.0 / 42.0;
            let mut t = t_first;
            for _ in 0..=20 {
                let p1 = curve1.value(t);
                let p2 = curve2.value(t);
                let center = path.value(t);
                let mut ns1 = GpVec::from_pnts(&p1, &center);
                let mut ns2 = GpVec::from_pnts(&p2, &center);
                ns1.normalize();
                ns2.normalize();
                dis.accumulate(&ns1, &ns2, &p1, &p2, &mut bary);
                t += delta;
            }
            bary = bary.multiplied(percent);
        }
        dis.bary = GpPnt::from_xyz(&bary);
        dis
    }

    /// The trimmed curves `SetInterval` installed (`cxx:635-642`).
    fn trimmed(&self) -> (Arc<dyn Curve>, Arc<dyn Curve>, Arc<dyn Curve>) {
        let g = self.trimmed.lock().unwrap_or_else(|e| e.into_inner());
        g.clone()
    }

    /// The (centre, normal1, normal2) frame `D0` needs (`cxx:278-297`).
    fn frame(&self, param: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (tpath, tc1, tc2) = self.trimmed();
        let center = tpath.value(param);
        let p1 = tc1.value(param);
        let p2 = tc2.value(param);
        let mut ns1 = GpVec::from_pnts(&p1, &center);
        let mut ns2 = GpVec::from_pnts(&p2, &center);
        let nplan = if !ns1.is_parallel_ang(&ns2, TOL_ANG) {
            ns1.crossed(&ns2)
        } else {
            let (_, mut n) = tpath.d1(param);
            if self.reverse {
                n = n.reversed();
            }
            n
        };
        ns1.normalize();
        ns2.normalize();
        let mut nplan = nplan;
        nplan.normalize();
        (center, ns1, ns2, nplan)
    }

    /// `temp.SetLinearForm(myRadius, ns1, myRadius, ns2, 1, P1, P2)` followed by
    /// `Center = 0.5 * temp` (`cxx:303-304`, `:371-372`, `:482-483`).
    fn blend_center(&self, ns1: &GpVec, ns2: &GpVec, p1: &GpPnt, p2: &GpPnt) -> GpPnt {
        let mut temp = GpXyz::zero();
        temp.set_linear_form_4(
            self.radius,
            ns1.xyz(),
            self.radius,
            ns2.xyz(),
            1.0,
            p1.xyz(),
            1.0,
            p2.xyz(),
        );
        GpPnt::from_xyz(&temp.multiplied(0.5))
    }
}

/// The `Discret` outputs (`cxx:147-266`).
struct Discret {
    max_ang: f64,
    min_ang: f64,
    dist_min: f64,
    bary: GpPnt,
}

impl Discret {
    /// One sample of the `Discret` loop (`cxx:191-218`): angle between the two
    /// normals, minimal pair distance and the barycentre accumulation.
    fn accumulate(&mut self, ns1: &GpVec, ns2: &GpVec, p1: &GpPnt, p2: &GpPnt, bary: &mut GpXyz) {
        let mut cosa = ns1.dot(ns2);
        if cosa > 1.0 {
            cosa = 1.0;
        }
        let angle = cosa.acos().abs();
        if angle > self.max_ang {
            self.max_ang = angle;
        }
        if angle < self.min_ang {
            self.min_ang = angle;
        }
        self.dist_min = self.dist_min.min(p1.distance(p2));
        *bary = bary.added(&p1.xyz().add(&p2.xyz()));
    }
}

impl ApproxSweepFunction for CircularBlendFunc {
    /// `D0(Param, ...)` (`cxx:267-310`).
    fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
    ) -> Result<(Vec<GpPnt>, Vec<GpPnt2d>, Vec<f64>), String> {
        let (_, ns1, ns2, nplan) = self.frame(param);
        let (_, tc1, tc2) = self.trimmed();
        let p1 = tc1.value(param);
        let p2 = tc2.value(param);
        let center = self.blend_center(&ns1, &ns2, &p1, &p2);
        let n = self.nb_poles.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); n];
        let mut weights = vec![1.0; n];
        get_circle(
            self.tconv, &ns1, &ns2, &nplan, &p1, &p2, self.radius, &center, &mut poles,
            &mut weights,
        );
        Ok((poles, Vec::new(), weights))
    }

    /// `D1(Param, ...)` (`cxx:311-401`).
    fn d1(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD1, String> {
        let (tpath, tc1, tc2) = self.trimmed();
        let (center0, dcenter) = tpath.d1(param);
        let (p1, dp1) = tc1.d1(param);
        let (p2, dp2) = tc2.d1(param);
        let mut ns1 = GpVec::from_pnts(&p1, &center0);
        let mut ns2 = GpVec::from_pnts(&p2, &center0);
        let mut dns1 = dcenter.subtracted(&dp1);
        let mut dns2 = dcenter.subtracted(&dp2);

        let (nplan, dnplan) = if !ns1.is_parallel_ang(&ns2, TOL_ANG) {
            (
                ns1.crossed(&ns2),
                dns1.crossed(&ns2).added(&ns1.crossed(&dns2)),
            )
        } else {
            let (_, mut n, mut dn) = tpath.d2(param);
            if self.reverse {
                n = n.reversed();
                dn = dn.reversed();
            }
            (n, dn)
        };

        let invnorm1 = 1.0 / ns1.magnitude();
        let invnorm2 = 1.0 / ns2.magnitude();
        ns1 = ns1.multiplied_scalar(invnorm1);
        dns1 = ns1.multiplied_scalar(-dns1.dot(&ns1)).added(&dns1);
        dns1 = dns1.multiplied_scalar(invnorm1);
        ns2 = ns2.multiplied_scalar(invnorm2);
        dns2 = ns2.multiplied_scalar(-dns2.dot(&ns2)).added(&dns2);
        dns2 = dns2.multiplied_scalar(invnorm2);

        let center = self.blend_center(&ns1, &ns2, &p1, &p2);
        let mut dcenter_out = GpXyz::zero();
        dcenter_out.set_linear_form_4(
            self.radius,
            dns1.xyz(),
            self.radius,
            dns2.xyz(),
            1.0,
            dp1.xyz(),
            1.0,
            dp2.xyz(),
        );
        let dcenter_out = GpVec::from_xyz(&dcenter_out).multiplied_scalar(0.5);

        let invnormp = 1.0 / nplan.magnitude();
        let nplan = nplan.multiplied_scalar(invnormp);
        let dnplan = nplan
            .multiplied_scalar(-dnplan.dot(&nplan))
            .added(&dnplan)
            .multiplied_scalar(invnormp);

        let n = self.nb_poles.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); n];
        let mut dpoles = vec![GpVec::zero(); n];
        let mut weights = vec![1.0; n];
        let mut dweights = vec![0.0; n];
        get_circle_d1(
            self.tconv,
            &ns1,
            &ns2,
            &dns1,
            &dns2,
            &nplan,
            &dnplan,
            &p1,
            &p2,
            &dp1,
            &dp2,
            self.radius,
            0.0,
            &center,
            &dcenter_out,
            &mut poles,
            &mut dpoles,
            &mut weights,
            &mut dweights,
        );
        Ok(SweepD1 {
            poles,
            dpoles,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
            weights,
            dweights,
        })
    }

    /// `D2(Param, ...)` (`cxx:404-533`).
    fn d2(&self, param: f64, _first: f64, _last: f64) -> Result<SweepD2, String> {
        let (tpath, tc1, tc2) = self.trimmed();
        let (center0, dcenter, d2center) = tpath.d2(param);
        let (p1, dp1, d2p1) = tc1.d2(param);
        let (p2, dp2, d2p2) = tc2.d2(param);
        let mut ns1 = GpVec::from_pnts(&p1, &center0);
        let mut dns1 = dcenter.subtracted(&dp1);
        let mut d2ns1 = d2center.subtracted(&d2p1);
        let mut ns2 = GpVec::from_pnts(&p2, &center0);
        let mut dns2 = dcenter.subtracted(&dp2);
        let mut d2ns2 = d2center.subtracted(&d2p2);

        let (nplan, dnplan, d2nplan) = if !ns1.is_parallel_ang(&ns2, TOL_ANG) {
            let nplan = ns1.crossed(&ns2);
            let dnplan = dns1.crossed(&ns2).added(&ns1.crossed(&dns2));
            let mut d2nplan = d2ns1.crossed(&ns2);
            d2nplan = d2nplan.added(&dns1.crossed(&dns2).multiplied_scalar(2.0));
            d2nplan = d2nplan.added(&ns1.crossed(&d2ns2));
            (nplan, dnplan, d2nplan)
        } else {
            let (_, mut n, mut dn, mut d2n) = tpath.d3(param);
            if self.reverse {
                n = n.reversed();
                dn = dn.reversed();
                d2n = d2n.reversed();
            }
            (n, dn, d2n)
        };

        let invnorm1 = 1.0 / ns1.magnitude();
        let invnorm2 = 1.0 / ns2.magnitude();
        ns1 = ns1.multiplied_scalar(invnorm1);
        let mut sc = dns1.dot(&ns1);
        let old = d2ns1.clone();
        d2ns1.set_linear_form_3(
            3.0 * sc * sc * invnorm1 - d2ns1.dot(&ns1) - invnorm1 * dns1.square_magnitude(),
            &ns1,
            -2.0 * sc * invnorm1,
            &dns1,
            1.0,
            &old,
        );
        dns1 = ns1.multiplied_scalar(-dns1.dot(&ns1)).added(&dns1);
        d2ns1 = d2ns1.multiplied_scalar(invnorm1);
        dns1 = dns1.multiplied_scalar(invnorm1);

        ns2 = ns2.multiplied_scalar(invnorm2);
        sc = dns2.dot(&ns2);
        let old = d2ns2.clone();
        d2ns2.set_linear_form_3(
            3.0 * sc * sc * invnorm2 - d2ns2.dot(&ns2) - invnorm2 * dns2.square_magnitude(),
            &ns2,
            -2.0 * sc * invnorm2,
            &dns2,
            1.0,
            &old,
        );
        dns2 = ns2.multiplied_scalar(-sc).added(&dns2);
        d2ns2 = d2ns2.multiplied_scalar(invnorm2);
        dns2 = dns2.multiplied_scalar(invnorm2);

        let center = self.blend_center(&ns1, &ns2, &p1, &p2);
        let mut dc = GpXyz::zero();
        dc.set_linear_form_4(
            self.radius,
            dns1.xyz(),
            self.radius,
            dns2.xyz(),
            1.0,
            dp1.xyz(),
            1.0,
            dp2.xyz(),
        );
        let dcenter_out = GpVec::from_xyz(&dc).multiplied_scalar(0.5);
        let mut d2c = GpXyz::zero();
        d2c.set_linear_form_4(
            self.radius,
            d2ns1.xyz(),
            self.radius,
            d2ns2.xyz(),
            1.0,
            d2p1.xyz(),
            1.0,
            d2p2.xyz(),
        );
        let d2center_out = GpVec::from_xyz(&d2c).multiplied_scalar(0.5);

        let invnormp = 1.0 / nplan.magnitude();
        let nplan = nplan.multiplied_scalar(invnormp);
        let mut dnplan = dnplan;
        let mut d2nplan = d2nplan;
        let sc = dnplan.dot(&nplan);
        let old = d2nplan.clone();
        d2nplan.set_linear_form_3(
            3.0 * sc * sc * invnormp - d2nplan.dot(&nplan) - invnormp * dnplan.square_magnitude(),
            &nplan,
            -2.0 * sc * invnormp,
            &dnplan,
            1.0,
            &old,
        );
        dnplan = nplan.multiplied_scalar(-sc).added(&dnplan);
        dnplan = dnplan.multiplied_scalar(invnormp);
        d2nplan = d2nplan.multiplied_scalar(invnormp);

        let n = self.nb_poles.max(0) as usize;
        let mut poles = vec![GpPnt::zero(); n];
        let mut dpoles = vec![GpVec::zero(); n];
        let mut d2poles = vec![GpVec::zero(); n];
        let mut weights = vec![1.0; n];
        let mut dweights = vec![0.0; n];
        let mut d2weights = vec![0.0; n];
        get_circle_d2(
            self.tconv,
            &ns1,
            &ns2,
            &dns1,
            &dns2,
            &d2ns1,
            &d2ns2,
            &nplan,
            &dnplan,
            &d2nplan,
            &p1,
            &p2,
            &dp1,
            &dp2,
            &d2p1,
            &d2p2,
            self.radius,
            0.0,
            0.0,
            &center,
            &dcenter_out,
            &d2center_out,
            &mut poles,
            &mut dpoles,
            &mut d2poles,
            &mut weights,
            &mut dweights,
            &mut d2weights,
        );
        Ok(SweepD2 {
            poles,
            dpoles,
            d2poles,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
            d2poles2d: Vec::new(),
            weights,
            dweights,
            d2weights,
        })
    }

    /// `Nb2dCurves` (`cxx:536-539`).
    fn nb_2d_curves(&self) -> i32 {
        0
    }

    /// `SectionShape` (`cxx:538-543`).
    fn section_shape(&self) -> (i32, i32, i32) {
        (self.nb_poles, self.nb_knots, self.degree)
    }

    /// `Knots` (`cxx:545-549`).
    fn knots(&self) -> Vec<f64> {
        let mut tk = vec![0.0; self.nb_knots.max(0) as usize];
        knots(self.tconv, &mut tk);
        tk
    }

    /// `Mults` (`cxx:551-555`).
    fn mults(&self) -> Vec<i32> {
        let mut tm = vec![0; self.nb_knots.max(0) as usize];
        mults(self.tconv, &mut tm);
        tm
    }

    /// `IsRational` (`cxx:557-560`).
    fn is_rational(&self) -> bool {
        self.tconv != ParameterisationType::Polynomial
    }

    /// `NbIntervals` (`cxx:561-598`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let icenter = self.path.parameter_intervals(next_shape(s));
        let icb1 = self.curve1.parameter_intervals(s);
        let icb2 = self.curve2.parameter_intervals(s);
        let epspar = PCONFUSION * 0.99;
        let icbs = fuse_intervals(&icb1, &icb2, epspar, false);
        let inter = fuse_intervals(&icenter, &icbs, epspar, false);
        inter.len() as i32 - 1
    }

    /// `Intervals` (`cxx:601-633`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let icenter = self.path.parameter_intervals(next_shape(s));
        let icb1 = self.curve1.parameter_intervals(s);
        let icb2 = self.curve2.parameter_intervals(s);
        let epspar = PCONFUSION * 0.99;
        let icbs = fuse_intervals(&icb1, &icb2, epspar, false);
        fuse_intervals(&icenter, &icbs, epspar, false)
    }

    /// `SetInterval` (`cxx:635-642`): trim the path and both rails.
    fn set_interval(&self, first: f64, last: f64) {
        let mut g = self.trimmed.lock().unwrap_or_else(|e| e.into_inner());
        *g = (
            Arc::new(GeomTrimmedCurveBasis::new(self.path.clone(), first, last)),
            Arc::new(GeomTrimmedCurveBasis::new(self.curve1.clone(), first, last)),
            Arc::new(GeomTrimmedCurveBasis::new(self.curve2.clone(), first, last)),
        );
    }

    /// `GetTolerance` (`cxx:644-657`).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64) -> Vec<f64> {
        let n = self.nb_poles.max(0) as usize;
        let mut tol3d = vec![surf_tol; n];
        if n < 2 {
            return tol3d;
        }
        // OCCT raises out of `GeomConvert`; the trait cannot, so the failure
        // keeps `SurfTol` everywhere.
        let Ok(tol) = get_tolerance(self.tconv, self.min_ang, self.radius, angle_tol, surf_tol)
        else {
            return tol3d;
        };
        tol3d[1] = tol.min(surf_tol);
        tol3d[n - 2] = tol.min(surf_tol);
        tol3d[0] = tol.min(bound_tol);
        tol3d[n - 1] = tol.min(bound_tol);
        tol3d
    }

    /// `SetTolerance` (`cxx:659-661`): nothing to do.
    fn set_tolerance(&self, _tol3d: f64, _tol2d: f64) {}

    /// `BarycentreOfSurf` (`cxx:663-666`).
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        Ok(self.bary.clone())
    }

    /// `MaximalSection` (`cxx:668-671`).
    fn maximal_section(&self) -> Result<f64, String> {
        Ok(self.max_ang * self.radius)
    }

    /// `GetMinimalWeight` (`cxx:673-677`).
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        let mut w = vec![1.0; self.nb_poles.max(0) as usize];
        get_minimal_weights(self.tconv, self.min_ang, self.max_ang, &mut w)
            .map_err(|e| format!("GeomFill_CircularBlendFunc::GetMinimalWeight: {e:?}"))?;
        Ok(w)
    }
}
