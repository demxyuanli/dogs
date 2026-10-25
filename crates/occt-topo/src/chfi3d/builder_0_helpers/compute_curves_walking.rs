//! The `IntWalk_PWalking` + `GeomInt_WLApprox` arm of `ChFi3d_ComputeCurves`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4072-4294`
//! (OCCT 8.0.0).
//!
//! The block starts its own tolerance/step loop (`Step` halves while the walk
//! answers at most 30 points); every other outcome leaves through a `return`,
//! so the loop runs at most five times before `Step` drops under `1.e-4`.
//!
//! Translation notes:
//! - OCCT truncates the `IntSurf_LineOn2S` in place (`Reverse`,
//!   `RemovePoint`, `InsertBefore`, `Value(i, p)`). The port's [`WLine`]
//!   owns its points vector, so the same edits are made on
//!   [`WLine::points`],
//! - the approximation result is rebuilt from the poles because OCCT resets
//!   `pol3d(1)` / `pol3d(nbpol)` after the fit (`cxx:4267-4281`). The port's
//!   [`WlApprox`] answers a degree-1 curve over the sampled points
//!   (`int_ss_make_bspline`), so the rebuilt curves use degree 1 and the
//!   same knot vector (`cxx:4283-4288`),
//! - `Approx.TolReached3d()` / `TolReached2d()` are commented out in OCCT
//!   (`cxx:4252-4255`); the port keeps `tolreached` from
//!   [`ch_fi3d_eval_tol_reached`] plus the two extremity distances.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::INFINITE;
use occt_geom::curve::Curve;
use occt_geom::{GeomBSplineCurve};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::chfi3d::builder_0_helpers::compute_curves::ComputeCurves;
use crate::chfi3d::builder_0_helpers::eval_tol_reached::ch_fi3d_eval_tol_reached;
use crate::geom_int::{ParametrizationType, WlApprox};
use crate::int_tools_wline::{PntOn2S, WLine};
use crate::intwalk::PWalking;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `Step` of the `IntWalk_PWalking` construction (`cxx:4075`).
const STEP_INIT: f64 = 0.1;

/// The smallest `Step` the loop accepts (`cxx:4115`).
const STEP_MIN: f64 = 0.0001;

/// The fewest points a usable walk must answer (`cxx:4111`).
const MIN_POINTS: i32 = 30;

/// `cxx:4072-4294`.
#[allow(clippy::too_many_arguments)]
pub fn walking_arm(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    pdeb: &GpPnt,
    pfin: &GpPnt,
    pardeb: [f64; 4],
    parfin: [f64; 4],
    tol3d: f64,
    tol2d: f64,
    out: &mut ComputeCurves,
) -> bool {
    let (Some(surf1), Some(surf2)) = (s1.surface().cloned(), s2.surface().cloned()) else {
        return false;
    };
    let mut step = STEP_INIT;
    loop {
        // `double fleche = 1.e-3 * pdeb.Distance(pfin); double tolap =
        // 1.e-7; IntWalk_PWalking IntKK(S1, S2, tol3d, tol3d, fleche, Step);`
        // (`cxx:4082-4084`).
        let fleche = 1.0e-3 * pdeb.distance(pfin);
        let tolap = 1.0e-7;
        let mut int_kk = PWalking::new(surf1.as_ref(), surf2.as_ref(), tol3d, tol3d, fleche, step);

        // `cxx:4088-4102`: walk out from a point between the two known
        // extremities.
        let mut found: Option<PntOn2S> = None;
        for ipdep in 2..=7 {
            let alpha = 0.1 * ipdep as f64;
            let unmoinsalpha = 1.0 - alpha;
            let depart = [
                alpha * pardeb[0] + unmoinsalpha * parfin[0],
                alpha * pardeb[1] + unmoinsalpha * parfin[1],
                alpha * pardeb[2] + unmoinsalpha * parfin[2],
                alpha * pardeb[3] + unmoinsalpha * parfin[3],
            ];
            if let Some(pintdep) = int_kk.perform_first_point(depart) {
                found = Some(pintdep);
                break;
            }
        }
        // `if (!depok) { return false; }` (`cxx:4101-4104`).
        let Some(pintdep) = found else {
            return false;
        };

        // `pintdep.Parameters(depart(1), depart(2), depart(3), depart(4));
        // IntKK.Perform(depart);` (`cxx:4105-4106`).
        let depart = [pintdep.u1, pintdep.v1, pintdep.u2, pintdep.v2];
        int_kk.perform(depart);
        if !int_kk.is_done() {
            return false;
        }
        // `cxx:4110-4120`.
        if int_kk.nb_points() <= MIN_POINTS {
            step *= 0.5;
            if step <= STEP_MIN {
                return false;
            }
            continue;
        }

        // `const occ::handle<IntSurf_LineOn2S>& L2S = IntKK.Line();`
        // (`cxx:4125`).
        let mut l2s = WLine::new();
        l2s.points = int_kk.line().to_vec();

        // `cxx:4127-4137`: where the walk ends relative to the known
        // extremities.
        let value = |s: &SurfaceAdaptor, u: f64, v: f64| s.value(u, v).unwrap_or_default();
        let codeb1 = value(s1, pardeb[0], pardeb[1]);
        let codeb2 = value(s2, pardeb[2], pardeb[3]);
        let tol1 = codeb1.distance(&codeb2).max(tol3d);
        let bondeb = tol1 == tol3d;
        let pntd = midpoint(&codeb1, &codeb2);

        let cofin1 = value(s1, parfin[0], parfin[1]);
        let cofin2 = value(s2, parfin[2], parfin[3]);
        let tol2 = cofin1.distance(&cofin2).max(tol3d);
        let bonfin = tol2 == tol3d;
        let pntf = midpoint(&cofin1, &cofin2);

        // `cxx:4139-4158`: the walk points closest to each extremity.
        let mut nbp = l2s.nb_pnts();
        let mut ddeb = INFINITE;
        let mut dfin = INFINITE;
        let mut indd = 0;
        let mut indf = 0;
        for i in 1..=nbp {
            let dd = l2s.point(i).value().distance(&pntd);
            if dd <= ddeb {
                ddeb = dd;
                indd = i;
            }
            let dd = l2s.point(i).value().distance(&pntf);
            if dd < dfin {
                dfin = dd;
                indf = i;
            }
        }
        // `cxx:4159-4164`.
        if indd > indf {
            l2s.points.reverse();
            indd = nbp - indd + 1;
            indf = nbp - indf + 1;
        }
        // `cxx:4165-4178`: cut the walk down to `indd .. indf`.
        for _ in 1..indd {
            l2s.remove_point(1);
            nbp -= 1;
            indf -= 1;
        }
        for _ in indf + 1..=nbp {
            l2s.remove_point(indf + 1);
        }
        nbp = indf;
        if nbp == 1 {
            return false;
        }

        // `cxx:4180-4211`: insert `pntd` at the start when the walk stops far
        // from it while the corresponding input point is good.
        if ddeb >= tol3d && bondeb {
            let p1 = *l2s.point(1);
            let p2 = *l2s.point(2);
            let v1 = GpVec::from_pnts(&pntd, &p1.value());
            let v2 = GpVec::from_pnts(&p1.value(), &p2.value());
            let v3 = GpVec::from_pnts(&pntd, &p2.value());
            let mut p1new = p1;
            p1new.set_value(pntd, pardeb[0], pardeb[1], pardeb[2], pardeb[3]);
            if v1.dot(&v3) < 0.0 {
                if v3.magnitude() < 0.2 * v2.magnitude() {
                    l2s.remove_point(1);
                    nbp -= 1;
                }
                *l2s.point_mut(1) = p1new;
            } else if v1.magnitude() > 0.2 * v2.magnitude() {
                l2s.insert_before(1, p1new);
                nbp += 1;
            } else {
                *l2s.point_mut(1) = p1new;
            }
            ddeb = 0.0;
        }
        // `cxx:4212-4239`: the same at the end.
        if dfin >= tol3d && bonfin {
            let p1 = *l2s.point(nbp);
            let p2 = *l2s.point(nbp - 1);
            let v1 = GpVec::from_pnts(&pntf, &p1.value());
            let v2 = GpVec::from_pnts(&p1.value(), &p2.value());
            let v3 = GpVec::from_pnts(&pntf, &p2.value());
            let mut p1new = p1;
            p1new.set_value(pntf, parfin[0], parfin[1], parfin[2], parfin[3]);
            if v1.dot(&v3) < 0.0 {
                if v3.magnitude() < 0.2 * v2.magnitude() {
                    l2s.remove_point(nbp);
                    nbp -= 1;
                }
                *l2s.point_mut(nbp) = p1new;
            } else if v1.magnitude() > 0.2 * v2.magnitude() {
                l2s.add(p1new);
                nbp += 1;
            } else {
                *l2s.point_mut(nbp) = p1new;
            }
            dfin = 0.0;
        }

        // `GeomInt_WLApprox approx; approx.SetParameters(tolap, tol2d, 4, 8,
        // 0, 30, true); approx.Perform(S1, S2, WL, true, true, true, 1, nbp);`
        // (`cxx:4246-4249`).
        let mut approx = WlApprox::new();
        approx.set_parameters(
            tolap,
            tol2d,
            4,
            8,
            0,
            30,
            true,
            ParametrizationType::ChordLength,
        );
        approx.perform(surf1.as_ref(), surf2.as_ref(), &l2s, true, true, true, 1, nbp);
        if !approx.is_done() {
            return false;
        }
        // `const AppParCurves_MultiBSpCurve& mbs = approx.Value(1);`
        // (`cxx:4258`).
        let Some(mbs) = approx.value(1).cloned() else {
            return false;
        };
        if mbs.poles3d.len() < 2 {
            return false;
        }
        let mut pol3d = mbs.poles3d.clone();
        let mut pol2d1 = mbs.poles2d_s1.clone();
        let mut pol2d2 = mbs.poles2d_s2.clone();
        let nbpol = pol3d.len();

        // `cxx:4267-4281`: reset the extremities on the known input points.
        if ddeb >= tol1 {
            pol3d[0] = pntd;
            pol2d1[0] = GpPnt2d::new(pardeb[0], pardeb[1]);
            pol2d2[0] = GpPnt2d::new(pardeb[2], pardeb[3]);
        }
        if dfin >= tol2 {
            pol3d[nbpol - 1] = pntf;
            pol2d1[nbpol - 1] = GpPnt2d::new(parfin[0], parfin[1]);
            pol2d2[nbpol - 1] = GpPnt2d::new(parfin[2], parfin[3]);
        }

        // `cxx:4283-4288`: `Geom_BSplineCurve(pol3d, knots, mults, deg)` and
        // the two `Geom2d_BSplineCurve` twins.
        let Ok(c3d) = GeomBSplineCurve::new(pol3d, mbs.knots.clone(), 1) else {
            return false;
        };
        let c3d: Arc<dyn Curve> = Arc::new(c3d);
        let Ok(pc1) = bspline2d(&pol2d1, &mbs.knots) else {
            return false;
        };
        let Ok(pc2) = bspline2d(&pol2d2, &mbs.knots) else {
            return false;
        };
        let pc1: Arc<dyn Curve2d> = Arc::new(pc1);
        let pc2: Arc<dyn Curve2d> = Arc::new(pc2);

        // `cxx:4289-4291`.
        let mut tolreached = ch_fi3d_eval_tol_reached(
            surf1.as_ref(),
            pc1.as_ref(),
            surf2.as_ref(),
            pc2.as_ref(),
            c3d.as_ref(),
        );
        tolreached = tolreached.max(ddeb).max(dfin);

        out.c3d = Some(c3d.clone());
        out.pc1 = Some(pc1);
        out.pc2 = Some(pc2);
        out.pardeb = c3d.first_parameter();
        out.parfin = c3d.last_parameter();
        out.tolreached = tolreached;
        return true;
    }
}

/// `gp_Pnt(0.5 * (A.Coord() + B.Coord()))`.
fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (a.x() + b.x()),
        0.5 * (a.y() + b.y()),
        0.5 * (a.z() + b.z()),
    )
}

/// A degree-1 `Geom2d_BSplineCurve` through the poles.
fn bspline2d(poles: &[GpPnt2d], knots: &[f64]) -> Result<Geom2dBSplineCurve, &'static str> {
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots.to_vec(), 1)
}
