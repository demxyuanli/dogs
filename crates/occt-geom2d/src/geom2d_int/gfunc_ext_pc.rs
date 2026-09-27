//! Port of `Extrema_GFuncExtPC` (`Extrema_GFuncExtPC.hxx:32-482`)
//! monomorphised onto the `Geom2dInt` instantiation
//! `Extrema_GFuncExtPC<Adaptor2d_Curve2d, Geom2dInt_Geom2dCurveTool,
//! Extrema_POnCurv2d, gp_Pnt2d, gp_Vec2d, NCollection_Sequence<Extrema_POnCurv2d>>`
//! (`Geom2dInt_PCLocFOfTheLocateExtPCOfTheProjPCurOfGInter.hxx:29-35`).
//!
//! This is the `math_FunctionWithDerivative` whose zero is searched by
//! `Extrema_GenLocateExtPC` (and hence by `Geom2dInt_TheProjPCurOfGInter`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::MathFunctionWithDerivative;
use occt_core::precision::{INFINITE, Precision};

use crate::curve::Curve2d;
use super::curve_tool::{self, GeomAbsCurveType};

const TOL_FACTOR: f64 = 1.0e-12;
const MIN_TOL: f64 = 1.0e-20;
const MIN_STEP: f64 = 1.0e-7;
const MAX_ORDER: i32 = 3;

/// `Extrema_POnCurv2d` (`Extrema_POnCurv2d.hxx`): a parameter and the point.
#[derive(Clone, Copy, Debug)]
pub struct POnCurv2d {
    pub parameter: f64,
    pub point: GpPnt2d,
}

impl POnCurv2d {
    pub fn new(parameter: f64, point: GpPnt2d) -> Self {
        Self { parameter, point }
    }

    pub fn parameter(&self) -> f64 {
        self.parameter
    }

    pub fn value(&self) -> GpPnt2d {
        self.point
    }
}

/// `Extrema_GFuncExtPC` (`hxx:53-480`), 2D instantiation.
pub struct GFuncExtPC<'a> {
    my_p: GpPnt2d,
    my_c: Option<&'a dyn Curve2d>,
    my_u: f64,
    my_pc: GpPnt2d,
    my_d1f: f64,
    my_sq_dist: Vec<f64>,
    my_is_min: Vec<i32>,
    my_point: Vec<POnCurv2d>,
    my_p_init: bool,
    my_c_init: bool,
    my_d1_init: bool,
    my_tol: f64,
    my_max_deriv_order: i32,
    my_u_infium: f64,
    my_u_supremum: f64,
}

impl<'a> Default for GFuncExtPC<'a> {
    /// `Extrema_GFuncExtPC()` (`hxx:59-71`).
    fn default() -> Self {
        Self {
            my_p: GpPnt2d::new(0.0, 0.0),
            my_c: None,
            my_u: 0.0,
            my_pc: GpPnt2d::new(0.0, 0.0),
            my_d1f: 0.0,
            my_sq_dist: Vec::new(),
            my_is_min: Vec::new(),
            my_point: Vec::new(),
            my_p_init: false,
            my_c_init: false,
            my_d1_init: false,
            my_tol: MIN_TOL,
            my_max_deriv_order: 0,
            my_u_infium: 0.0,
            my_u_supremum: 0.0,
        }
    }
}

impl<'a> GFuncExtPC<'a> {
    /// `Initialize(theC)` (`hxx:105-129`).
    pub fn initialize(&mut self, c: &'a dyn Curve2d) {
        self.my_c = Some(c);
        self.my_c_init = true;
        self.my_point.clear();
        self.my_sq_dist.clear();
        self.my_is_min.clear();
        self.sub_interval_initialize(curve_tool::first_parameter(c), curve_tool::last_parameter(c));
        match curve_tool::get_type(c) {
            GeomAbsCurveType::BezierCurve
            | GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::OtherCurve => {
                self.my_max_deriv_order = MAX_ORDER;
                self.my_tol = self.search_of_tolerance();
            }
            _ => {
                self.my_max_deriv_order = 0;
                self.my_tol = MIN_TOL;
            }
        }
    }

    /// `SetPoint(theP)` (`hxx:133-140`).
    pub fn set_point(&mut self, p: &GpPnt2d) {
        self.my_p = *p;
        self.my_p_init = true;
        self.my_point.clear();
        self.my_sq_dist.clear();
        self.my_is_min.clear();
    }

    /// `SubIntervalInitialize` (`hxx:420-424`).
    pub fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.my_u_infium = u_first;
        self.my_u_supremum = u_last;
    }

    /// `SearchOfTolerance` (`hxx:428-457`).
    fn search_of_tolerance(&self) -> f64 {
        const N_POINT: i32 = 10;
        let a_step = (self.my_u_supremum - self.my_u_infium) / f64::from(N_POINT);
        let c = self.my_c.expect("curve not initialized");
        let mut a_num = 0;
        let mut a_max = -INFINITE;
        loop {
            let mut u = self.my_u_infium + f64::from(a_num) * a_step;
            if u > self.my_u_supremum {
                u = self.my_u_supremum;
            }
            let (_p, v_der) = curve_tool::d1(c, u);
            if !Precision::is_infinite(v_der.x()) && !Precision::is_infinite(v_der.y()) {
                let vm = v_der.magnitude();
                if vm > a_max {
                    a_max = vm;
                }
            }
            a_num += 1;
            if a_num >= N_POINT + 1 {
                break;
            }
        }
        (a_max * TOL_FACTOR).max(MIN_TOL)
    }

    /// `Value` (`hxx:146-253`).
    fn value_impl(&mut self, the_u: f64, the_f: &mut f64) -> bool {
        let c = match self.my_c {
            Some(c) => c,
            None => return false,
        };
        self.my_u = the_u;
        let (pc, mut d1c) = curve_tool::d1(c, self.my_u);
        self.my_pc = pc;
        if Precision::is_infinite(d1c.x()) || Precision::is_infinite(d1c.y()) {
            *the_f = INFINITE;
            return false;
        }
        let mut ndu = d1c.magnitude();
        if self.my_max_deriv_order != 0 && ndu <= self.my_tol {
            const DIVISION_FACTOR: f64 = 1.0e-3;
            let du = if self.my_u_supremum >= f64::MAX || self.my_u_infium <= -f64::MAX {
                0.0
            } else {
                self.my_u_supremum - self.my_u_infium
            };
            let a_delta = (du * DIVISION_FACTOR).max(MIN_STEP);
            let mut n = 1i32;
            let mut v = GpVec2d::new(0.0, 0.0);
            let mut is_derive_found = false;
            loop {
                n += 1;
                v = curve_tool::dn(c, self.my_u, n);
                ndu = v.magnitude();
                is_derive_found = ndu > self.my_tol;
                if is_derive_found || n >= self.my_max_deriv_order {
                    break;
                }
            }
            if is_derive_found {
                let u = if self.my_u - self.my_u_infium < a_delta {
                    self.my_u + a_delta
                } else {
                    self.my_u - a_delta
                };
                let p1 = curve_tool::d0(c, self.my_u.min(u));
                let p2 = curve_tool::d0(c, self.my_u.max(u));
                let v1 = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
                let a_dir_factor = v.dot(&v1);
                d1c = if a_dir_factor < 0.0 {
                    GpVec2d::new(-v.x(), -v.y())
                } else {
                    v
                };
            } else {
                let (p1, p2, p3, grown);
                if self.my_u - self.my_u_infium < 2.0 * a_delta {
                    p1 = curve_tool::d0(c, self.my_u);
                    p2 = curve_tool::d0(c, self.my_u + a_delta);
                    p3 = curve_tool::d0(c, self.my_u + 2.0 * a_delta);
                    grown = true;
                } else {
                    p1 = curve_tool::d0(c, self.my_u - 2.0 * a_delta);
                    p2 = curve_tool::d0(c, self.my_u - a_delta);
                    p3 = curve_tool::d0(c, self.my_u);
                    grown = false;
                }
                let (v1, v2, v3) = (
                    GpVec2d::new(p1.x(), p1.y()),
                    GpVec2d::new(p2.x(), p2.y()),
                    GpVec2d::new(p3.x(), p3.y()),
                );
                d1c = if grown {
                    GpVec2d::new(-3.0 * v1.x() + 4.0 * v2.x() - v3.x(), -3.0 * v1.y() + 4.0 * v2.y() - v3.y())
                } else {
                    GpVec2d::new(v1.x() - 4.0 * v2.x() + 3.0 * v3.x(), v1.y() - 4.0 * v2.y() + 3.0 * v3.y())
                };
            }
            ndu = d1c.magnitude();
        }
        if ndu <= MIN_TOL {
            return false;
        }
        let ppc = GpVec2d::new(self.my_pc.x() - self.my_p.x(), self.my_pc.y() - self.my_p.y());
        *the_f = ppc.dot(&d1c) / ndu;
        true
    }

    /// `Values` (`hxx:274-353`).
    fn values_impl(&mut self, the_u: f64, the_f: &mut f64, the_df: &mut f64) -> bool {
        let c = match self.my_c {
            Some(c) => c,
            None => return false,
        };
        let pc_old = self.my_pc;
        let p_old = self.my_p;
        if !self.value_impl(the_u, the_f) {
            self.my_d1_init = false;
            return false;
        }
        self.my_u = the_u;
        self.my_pc = pc_old;
        self.my_p = p_old;
        let (pc, d1c, d2c) = curve_tool::d2(c, self.my_u);
        self.my_pc = pc;
        let ndu = d1c.magnitude();
        if ndu <= self.my_tol {
            const DIVISION_FACTOR: f64 = 0.01;
            let du = if self.my_u_supremum >= f64::MAX || self.my_u_infium <= -f64::MAX {
                0.0
            } else {
                self.my_u_supremum - self.my_u_infium
            };
            let a_delta = (du * DIVISION_FACTOR).max(MIN_STEP);
            if self.my_u - self.my_u_infium < 2.0 * a_delta {
                let f1 = *the_f;
                let u2 = self.my_u + a_delta;
                let u3 = self.my_u + a_delta * 2.0;
                let (mut f2, mut f3) = (0.0, 0.0);
                if !self.value_impl(u2, &mut f2) || !self.value_impl(u3, &mut f3) {
                    self.my_d1_init = false;
                    return false;
                }
                *the_df = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * a_delta);
            } else {
                let f3 = *the_f;
                let u1 = self.my_u - a_delta * 2.0;
                let u2 = self.my_u - a_delta;
                let (mut f1, mut f2) = (0.0, 0.0);
                if !self.value_impl(u2, &mut f2) || !self.value_impl(u1, &mut f1) {
                    self.my_d1_init = false;
                    return false;
                }
                *the_df = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * a_delta);
            }
            self.my_u = the_u;
            self.my_pc = pc_old;
            self.my_p = p_old;
        } else {
            let ppc = GpVec2d::new(self.my_p.x() - self.my_pc.x(), self.my_p.y() - self.my_pc.y());
            *the_df = ndu + (ppc.dot(&d2c) / ndu) - *the_f * (d1c.dot(&d2c)) / (ndu * ndu);
        }
        self.my_d1f = *the_df;
        self.my_d1_init = true;
        true
    }

    /// `GetStateNumber` (`hxx:357-379`).
    fn get_state_number_impl(&mut self) -> i32 {
        self.my_sq_dist.push(self.my_pc.square_distance(&self.my_p));
        self.my_d1_init = true;
        let (mut ff, mut dd) = (0.0, 0.0);
        self.values_impl(self.my_u, &mut ff, &mut dd);
        let int_val = if self.my_d1f > 0.0 { 1 } else { 0 };
        self.my_is_min.push(int_val);
        self.my_point.push(POnCurv2d::new(self.my_u, self.my_pc));
        0
    }

    /// `NbExt` (`hxx:382`).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// `SquareDistance` (`hxx:386-393`), 1-based.
    pub fn square_distance(&self, the_n: usize) -> f64 {
        self.my_sq_dist[the_n - 1]
    }

    /// `IsMin` (`hxx:397-404`), 1-based.
    pub fn is_min(&self, the_n: usize) -> bool {
        self.my_is_min[the_n - 1] == 1
    }

    /// `Point` (`hxx:408-415`), 1-based.
    pub fn point(&self, the_n: usize) -> POnCurv2d {
        self.my_point[the_n - 1]
    }
}

impl MathFunctionWithDerivative for GFuncExtPC<'_> {
    /// `Value` (`hxx:146`).
    fn value(&mut self, param: f64, approx_distance: &mut f64) -> bool {
        self.value_impl(param, approx_distance)
    }

    /// `Derivative` (`hxx:259-267`).
    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let mut f = 0.0;
        self.values_impl(param, &mut f, d)
    }

    /// `Values` (`hxx:274`).
    fn values(&mut self, param: f64, approx_distance: &mut f64, deriv: &mut f64) -> bool {
        self.values_impl(param, approx_distance, deriv)
    }

    /// `GetStateNumber` (`hxx:357-379`).
    fn get_state_number(&mut self) -> i32 {
        self.get_state_number_impl()
    }
}
