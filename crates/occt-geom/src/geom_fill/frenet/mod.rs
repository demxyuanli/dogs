//! `GeomFill_Frenet` — the Frenet trihedron law, including the singular-point
//! search of `Init`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Frenet.hxx` and
//! `GeomFill_Frenet.cxx:30-1108` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `mySngl` / `mySnglLen` / `isSngl` (`hxx:86-88`) are written by `Init`,
//!   which every other law calls from the non-const `SetCurve` through a
//!   shared handle. The port takes `&self` there (see
//!   [`crate::geom_fill::trihedron_law`]), so they live behind a `Mutex` as
//!   [`SingularPoints`]; `None` is `isSngl == false`.
//! - `Extrema_ExtPC` (`cxx:203`, `cxx:261`) is
//!   [`crate::extrema_pc::extrema_ext_pc_range_done`].
//! - `P` is an OCCT member (`hxx:80`) used as scratch by `D0`; the port keeps
//!   it a local, since no other call reads it.
//! - `Standard_OutOfRange` from the `NbIntervals` / `Intervals` switch
//!   (`cxx:781`, `cxx:820`) is a panic, OCCT's exception counterpart.

pub(crate) mod deriv;
mod singular;

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpVec, GpXyz};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{PCONFUSION, RESOLUTION};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::snglr_func::SnglrFunc;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};

/// `NullTol` (`cxx:38`).
pub(crate) const NULL_TOL: f64 = 1.0e-10;
/// `MaxSingular` (`cxx:39`).
pub(crate) const MAX_SINGULAR: f64 = 1.0e-5;
/// `maxDerivOrder` (`cxx:41`).
pub(crate) const MAX_DERIV_ORDER: i32 = 3;

/// `mySngl` / `mySnglLen` (`hxx:86-87`) once `isSngl` is true.
#[derive(Clone, Debug, Default)]
pub(crate) struct SingularPoints {
    /// `mySngl` (`hxx:86`): the singular parameters.
    pub sngl: Vec<f64>,
    /// `mySnglLen` (`hxx:87`): the half-width of each singular interval.
    pub len: Vec<f64>,
}

/// `GeomFill_Frenet` (`GeomFill_Frenet.hxx:37-95`).
pub struct Frenet {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `isSngl`, `mySngl`, `mySnglLen` (`hxx:85-87`); `None` is `!isSngl`.
    singular: Mutex<Option<SingularPoints>>,
}

impl Default for Frenet {
    fn default() -> Self {
        Self::new()
    }
}

impl Frenet {
    /// `GeomFill_Frenet::GeomFill_Frenet` (`cxx:99-102`): `isSngl = false`.
    pub fn new() -> Self {
        Self {
            data: TrihedronLawCurve::new(),
            singular: Mutex::new(None),
        }
    }

    /// `isSngl` (`hxx:85`) together with the two arrays.
    pub(crate) fn singular_points(&self) -> Option<SingularPoints> {
        self.singular
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .clone()
    }

    pub(crate) fn set_singular_points(&self, value: Option<SingularPoints>) {
        *self.singular.lock().unwrap_or_else(|e| e.into_inner()) = value;
    }

    /// `myTrimmed` (`hxx:79`).
    pub(crate) fn trimmed(&self) -> Result<Arc<dyn Curve>, String> {
        self.data.curve_or_err("GeomFill_Frenet")
    }

    /// `myCurve->GetType()` for the `IsConstant` / `SetCurve` switches.
    pub(crate) fn curve_type(&self) -> Option<GeomAbsCurveType> {
        self.data.curve().map(|c| extrema_curve_type(c.as_ref()))
    }

    /// `GeomFill_Frenet::Init` (`cxx:147-392`): the singular-point search.
    ///
    /// The OCCT loop keeps `Average` running across the C2 intervals
    /// (`cxx:168`, `cxx:185`, `cxx:193`) -- `AveFunc(i)` is the *cumulative*
    /// mean divided by `NbControl`. That is reproduced verbatim.
    fn init(&self) -> Result<(), String> {
        let curve = self.trimmed()?;
        let nb_int_c2 = curve.nb_intervals(Shape::C2 as u8);
        if nb_int_c2 <= 0 {
            self.set_singular_points(None);
            return Ok(());
        }
        let c2_disc = curve.parameter_intervals(Shape::C2 as u8);

        const TOL_F: f64 = 1.0e-10;
        const TOL: f64 = 10.0 * TOL_F;
        const TOL2: f64 = TOL * TOL;
        const PTOL: f64 = PCONFUSION;
        const NB_CONTROL: usize = 10;

        let mut func = SnglrFunc::new(curve.clone());
        let origin = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);

        let n = nb_int_c2 as usize;
        let mut is_lin = vec![true; n];
        let mut is_const = vec![true; n];
        let mut ave_func = vec![0.0f64; n];
        let mut average = 0.0f64;

        for i in 0..n {
            let step = (c2_disc[i + 1] - c2_disc[i]) / NB_CONTROL as f64;
            let mut c1 = GpXyz::zero();
            for j in 0..NB_CONTROL {
                let c = func.d0(c2_disc[i] + j as f64 * step);
                if j == 0 {
                    c1 = *c.xyz();
                }
                let modulus = c.xyz().modulus();
                if modulus > TOL {
                    is_lin[i] = false;
                }
                average += modulus;
                if is_const[i]
                    && ((c.x() - c1.x).abs() > TOL
                        || (c.y() - c1.y).abs() > TOL
                        || (c.z() - c1.z).abs() > TOL)
                {
                    is_const[i] = false;
                }
            }
            ave_func[i] = average / NB_CONTROL as f64;
        }

        // First sets of roots (`cxx:199-241`).
        let mut seq: Vec<Vec<f64>> = vec![Vec::new(); n];
        for i in 0..n {
            if is_lin[i] || is_const[i] {
                continue;
            }
            func.set_ratio(1.0 / ave_func[i]);
            let sols = crate::extrema_pc::extrema_ext_pc_range_done(
                &func,
                &origin,
                c2_disc[i],
                c2_disc[i + 1],
            )
            .0;
            for s in sols {
                if s.sq_dist < TOL2 {
                    seq[i].push(s.u);
                }
            }
            seq[i].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        }

        let mut sngl_seq: Vec<f64> = Vec::new();
        for s in &seq {
            sngl_seq.extend_from_slice(s);
        }

        // Extrema works bad, need to pass second time (`cxx:249-283`).
        for i in 0..n {
            if seq[i].is_empty() {
                continue;
            }
            let mut bounds = Vec::with_capacity(seq[i].len() + 2);
            bounds.push(c2_disc[i + 1]);
            bounds.extend_from_slice(&seq[i]);
            bounds.push(c2_disc[i + 2]);
            func.set_ratio(1.0 / ave_func[i]);
            for j in 0..bounds.len() - 1 {
                if bounds[j + 1] - bounds[j] <= PTOL {
                    continue;
                }
                let sols =
                    crate::extrema_pc::extrema_ext_pc_range_done(&func, &origin, bounds[j], bounds[j + 1])
                        .0;
                for s in sols {
                    if s.sq_dist < TOL2
                        && s.u - bounds[j] > PTOL
                        && bounds[j + 1] - s.u > PTOL
                    {
                        sngl_seq.push(s.u);
                    }
                }
            }
        }

        if sngl_seq.is_empty() {
            self.set_singular_points(None);
            return Ok(());
        }

        // Sorting then discarding repeating elements (`cxx:288-315`).
        sngl_seq.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut j = 0usize;
        loop {
            let mut found = false;
            let mut i = j;
            while i + 1 < sngl_seq.len() {
                if sngl_seq[i + 1] - sngl_seq[i] <= PTOL {
                    sngl_seq.remove(i + 1);
                    j = i;
                    found = true;
                    break;
                }
                i += 1;
            }
            if !found {
                break;
            }
        }

        // Computation of the length of the singular interval (`cxx:322-343`).
        let mut lens = Vec::with_capacity(sngl_seq.len());
        for &u in &sngl_seq {
            let (_, d1, d2) = func.d2(u);
            let n1 = d1.magnitude();
            let mut norm = n1;
            let mut value;
            if norm > RESOLUTION {
                value = NULL_TOL / norm;
            } else {
                norm = d2.magnitude();
                if norm > RESOLUTION {
                    value = (2.0 * NULL_TOL / norm).sqrt();
                } else {
                    value = MAX_SINGULAR;
                }
            }
            if value > MAX_SINGULAR {
                value = MAX_SINGULAR;
            }
            lens.push(value);
        }

        // Merge singular points that share parts of their intervals
        // (`cxx:351-379`).
        if sngl_seq.len() > 1 {
            let mut tmp: Vec<(f64, f64)> = vec![(sngl_seq[0], lens[0])];
            for i in 1..sngl_seq.len() {
                let (lx, ly) = *tmp.last().unwrap();
                let u12 = lx + ly;
                let u21 = sngl_seq[i] - lens[i];
                if u12 >= u21 {
                    let u11 = lx - ly;
                    let u22 = sngl_seq[i] + lens[i];
                    let last = tmp.len() - 1;
                    tmp[last] = ((u11 + u22) / 2.0, (u22 - u11) / 2.0);
                } else {
                    tmp.push((sngl_seq[i], lens[i]));
                }
            }
            sngl_seq = tmp.iter().map(|p| p.0).collect();
            lens = tmp.iter().map(|p| p.1).collect();
        }

        self.set_singular_points(Some(SingularPoints {
            sngl: sngl_seq,
            len: lens,
        }));
        Ok(())
    }
}

impl TrihedronLaw for Frenet {
    /// `GeomFill_Frenet::SetCurve` (`cxx:119-145`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve);
        let Some(curve) = self.data.curve() else {
            return true;
        };
        match extrema_curve_type(curve.as_ref()) {
            GeomAbsCurveType::Line
            | GeomAbsCurveType::Circle
            | GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::Hyperbola
            | GeomAbsCurveType::Parabola => {
                self.set_singular_points(None);
            }
            _ => {
                let _ = self.init();
            }
        }
        true
    }

    /// `GeomFill_Frenet::Copy` (`cxx:107-114`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let copy = Frenet::new();
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&copy, curve);
        }
        Arc::new(copy)
    }

    /// `GeomFill_Frenet::D0` (`cxx:472-626`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let tol = RESOLUTION;
        let mut delta = 0.0;
        if let Some(index) = self.is_singular(param) {
            let (frame, d) = self.singular_d0(param, index)?;
            if let Some(frame) = frame {
                return Ok(frame);
            }
            delta = d;
        }
        let a_param = param + delta;
        let curve = self.trimmed()?;
        let (_p, mut tangent, mut bi_normal) = curve.d2(a_param);

        const DIVISION_FACTOR: f64 = 1.0e-3;
        let (u_inf, u_sup) = self.data.interval();
        let a_delta = (u_sup - u_inf) * DIVISION_FACTOR;
        let ndu = tangent.magnitude();

        if ndu <= tol {
            // Derivative approximated by Taylor series (`cxx:501-539`), then by
            // three points (`cxx:540-575`).
            let mut index = 1;
            let mut found = false;
            let mut a_tn = GpVec::zero();
            while index < MAX_DERIV_ORDER {
                index += 1;
                a_tn = curve.eval_dn(param, index);
                if a_tn.magnitude() > tol {
                    found = true;
                    break;
                }
            }
            if found {
                let u = if param - u_inf < a_delta {
                    param + a_delta
                } else {
                    param - a_delta
                };
                let p1 = curve.d0(param.min(u));
                let p2 = curve.d0(param.max(u));
                let v1 = occt_core::gp::GpVec::from_pnts(&p1, &p2);
                if a_tn.dot(&v1) < 0.0 {
                    a_tn = a_tn.reversed();
                }
            } else {
                let p1;
                let p2;
                let p3;
                let grown;
                if param - u_inf < 2.0 * a_delta {
                    p1 = curve.d0(param);
                    p2 = curve.d0(param + a_delta);
                    p3 = curve.d0(param + 2.0 * a_delta);
                    grown = true;
                } else {
                    p1 = curve.d0(param - 2.0 * a_delta);
                    p2 = curve.d0(param - a_delta);
                    p3 = curve.d0(param);
                    grown = false;
                }
                let origin = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);
                let v1 = occt_core::gp::GpVec::from_pnts(&origin, &p1);
                let v2 = occt_core::gp::GpVec::from_pnts(&origin, &p2);
                let v3 = occt_core::gp::GpVec::from_pnts(&origin, &p3);
                a_tn = if grown {
                    v1.multiplied_scalar(-3.0)
                        .added(&v2.multiplied_scalar(4.0))
                        .subtracted(&v3)
                } else {
                    v1.subtracted(&v2.multiplied_scalar(4.0))
                        .added(&v3.multiplied_scalar(3.0))
                };
            }
            let d_par = 10.0 * a_delta;
            let (mut tangent, mut normal, mut bi_normal) = if param - u_inf < d_par {
                self.d0(a_param + d_par)?
            } else {
                self.d0(a_param - d_par)?
            };
            if !deriv::rotate_trihedron(&mut tangent, &mut normal, &mut bi_normal, &a_tn) {
                return Err("GeomFill_Frenet::D0 : cannot coincide two tangents".into());
            }
            return Ok((tangent, normal, bi_normal));
        }

        tangent = tangent.divided(ndu);
        bi_normal = tangent.crossed(&bi_normal);
        let norm = bi_normal.magnitude();
        if norm <= RESOLUTION {
            let axe = occt_core::gp::GpAx2::from_axis(
                occt_core::gp::GpPnt::new(0.0, 0.0, 0.0),
                occt_core::gp::GpDir::from_vec(&tangent)
                    .map_err(|e| format!("GeomFill_Frenet::D0: {e}"))?,
            );
            bi_normal = GpVec::from_xyz(axe.y_direction().xyz());
        } else {
            bi_normal.normalize();
        }
        let normal = bi_normal.crossed(&tangent);
        Ok((tangent, normal, bi_normal))
    }

    /// `GeomFill_Frenet::D1` (`cxx:628-679`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let mut delta = 0.0;
        if let Some(index) = self.is_singular(param) {
            let (frame, d) = self.singular_d1(param, index)?;
            if let Some(frame) = frame {
                return Ok(frame);
            }
            delta = d;
        }
        let the_param = param + delta;
        let curve = self.trimmed()?;
        let (_p, dc1, dc2, dc3) = curve.d3(the_param);
        let tangent = dc1.normalized();

        if tangent.crossed(&dc2).magnitude() <= RESOLUTION {
            let axe = occt_core::gp::GpAx2::from_axis(
                occt_core::gp::GpPnt::new(0.0, 0.0, 0.0),
                occt_core::gp::GpDir::from_vec(&tangent)
                    .map_err(|e| format!("GeomFill_Frenet::D1: {e}"))?,
            );
            return Ok(TrihedronD1 {
                tangent,
                dtangent: GpVec::zero(),
                normal: GpVec::from_xyz(axe.x_direction().xyz()),
                dnormal: GpVec::zero(),
                binormal: GpVec::from_xyz(axe.y_direction().xyz()),
                dbinormal: GpVec::zero(),
            });
        }
        let bi_normal = tangent.crossed(&dc2).normalized();
        let normal = bi_normal.crossed(&tangent);
        let d_tangent = deriv::fderiv(&dc1, &dc2);

        let instead_dc1 = tangent.crossed(&dc2);
        let instead_dc2 = d_tangent.crossed(&dc2).added(&tangent.crossed(&dc3));
        let d_bi_normal = deriv::fderiv(&instead_dc1, &instead_dc2);
        let d_normal = d_bi_normal.crossed(&tangent).added(&bi_normal.crossed(&d_tangent));
        Ok(TrihedronD1 {
            tangent,
            dtangent: d_tangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_Frenet::D2` (`cxx:682-760`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let mut delta = 0.0;
        if let Some(index) = self.is_singular(param) {
            let (frame, d) = self.singular_d2(param, index)?;
            if let Some(frame) = frame {
                return Ok(frame);
            }
            delta = d;
        }
        let the_param = param + delta;
        let curve = self.trimmed()?;
        let (_p, dc1, dc2, dc3) = curve.d3(the_param);
        let dc4 = curve.eval_dn(the_param, 4);
        let tangent = dc1.normalized();

        if tangent.crossed(&dc2).magnitude() <= RESOLUTION {
            let axe = occt_core::gp::GpAx2::from_axis(
                occt_core::gp::GpPnt::new(0.0, 0.0, 0.0),
                occt_core::gp::GpDir::from_vec(&tangent)
                    .map_err(|e| format!("GeomFill_Frenet::D2: {e}"))?,
            );
            return Ok(TrihedronD2 {
                tangent,
                dtangent: GpVec::zero(),
                d2tangent: GpVec::zero(),
                normal: GpVec::from_xyz(axe.x_direction().xyz()),
                dnormal: GpVec::zero(),
                d2normal: GpVec::zero(),
                binormal: GpVec::from_xyz(axe.y_direction().xyz()),
                dbinormal: GpVec::zero(),
                d2binormal: GpVec::zero(),
            });
        }
        let bi_normal = tangent.crossed(&dc2).normalized();
        let normal = bi_normal.crossed(&tangent);
        let d_tangent = deriv::fderiv(&dc1, &dc2);
        let d2_tangent = deriv::dderiv(&dc1, &dc2, &dc3);

        let instead_dc1 = tangent.crossed(&dc2);
        let instead_dc2 = d_tangent.crossed(&dc2).added(&tangent.crossed(&dc3));
        let instead_dc3 = d2_tangent
            .crossed(&dc2)
            .added(&d_tangent.crossed(&dc3).multiplied_scalar(2.0))
            .added(&tangent.crossed(&dc4));
        let d_bi_normal = deriv::fderiv(&instead_dc1, &instead_dc2);
        let d2_bi_normal = deriv::dderiv(&instead_dc1, &instead_dc2, &instead_dc3);
        let d_normal = d_bi_normal.crossed(&tangent).added(&bi_normal.crossed(&d_tangent));
        let d2_normal = d2_bi_normal
            .crossed(&tangent)
            .added(&d_bi_normal.crossed(&d_tangent).multiplied_scalar(2.0))
            .added(&bi_normal.crossed(&d2_tangent));
        Ok(TrihedronD2 {
            tangent,
            dtangent: d_tangent,
            d2tangent: d2_tangent,
            normal,
            dnormal: d_normal,
            d2normal: d2_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
            d2binormal: d2_bi_normal,
        })
    }

    /// `GeomFill_Frenet::NbIntervals` (`cxx:764-798`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let tmp_s = frenet_shape(s);
        let curve = match self.data.curve() {
            Some(c) => c,
            None => return 1,
        };
        let nb_trimmed = curve.nb_intervals(tmp_s as u8);
        let Some(sngl) = self.singular_points() else {
            return nb_trimmed;
        };
        let trim_int = curve.parameter_intervals(tmp_s as u8);
        let fusion =
            crate::geom_lib::fuse_intervals(&trim_int, &sngl.sngl, PCONFUSION, true);
        fusion.len() as i32 - 1
    }

    /// `GeomFill_Frenet::Intervals` (`cxx:802-843`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let tmp_s = frenet_shape(s);
        let curve = match self.data.curve() {
            Some(c) => c,
            None => return vec![0.0, 0.0],
        };
        let Some(sngl) = self.singular_points() else {
            return curve.parameter_intervals(tmp_s as u8);
        };
        let trim_int = curve.parameter_intervals(tmp_s as u8);
        crate::geom_lib::fuse_intervals(&trim_int, &sngl.sngl, PCONFUSION, true)
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_Frenet::GetAverageLaw` (`cxx:845-870`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        const NUM: usize = 20;
        let mut a_tangent = GpVec::zero();
        let mut a_normal = GpVec::zero();
        let mut a_bi_normal = GpVec::zero();
        let (first, last) = self.data.interval();
        let step = (last - first) / NUM as f64;
        for i in 0..=NUM {
            let mut param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            if let Ok((t, n, bn)) = self.d0(param) {
                a_tangent = a_tangent.added(&t);
                a_normal = a_normal.added(&n);
                a_bi_normal = a_bi_normal.added(&bn);
            }
        }
        a_tangent = a_tangent.divided((NUM + 1) as f64);
        a_normal = a_normal.divided((NUM + 1) as f64);
        a_tangent.normalize();
        a_bi_normal = a_tangent.crossed(&a_normal).normalized();
        a_normal = a_bi_normal.crossed(&a_tangent);
        (a_tangent, a_normal, a_bi_normal)
    }

    /// `GeomFill_Frenet::IsConstant` (`cxx:874-878`): `GeomAbs_Line`.
    fn is_constant(&self) -> bool {
        matches!(self.curve_type(), Some(GeomAbsCurveType::Line))
    }

    /// `GeomFill_Frenet::IsOnlyBy3dCurve` (`cxx:882-885`).
    fn is_only_by_3d_curve(&self) -> bool {
        true
    }
}

/// The `NbIntervals` / `Intervals` switch (`cxx:767-782`).
fn frenet_shape(s: u8) -> Shape {
    match s {
        0 => Shape::C2,
        2 => Shape::C3,
        4 | 5 | 6 => Shape::CN,
        _ => panic!("Standard_OutOfRange: GeomFill_Frenet::Intervals"),
    }
}
