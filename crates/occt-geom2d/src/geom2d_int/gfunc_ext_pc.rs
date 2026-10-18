//! Port of `Extrema_GFuncExtPC` (`Extrema_GFuncExtPC.hxx:32-482`)
//! monomorphised onto the `Geom2dInt` instantiation
//! `Extrema_GFuncExtPC<Adaptor2d_Curve2d, Geom2dInt_Geom2dCurveTool,
//! Extrema_POnCurv2d, gp_Pnt2d, gp_Vec2d, NCollection_Sequence<Extrema_POnCurv2d>>`
//! (`Geom2dInt_PCLocFOfTheLocateExtPCOfTheProjPCurOfGInter.hxx:29-35`).
//!
//! This is the `math_FunctionWithDerivative` whose zero is searched by
//! `Extrema_GenLocateExtPC` (and hence by `Geom2dInt_TheProjPCurOfGInter`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::FunctionWithDerivative;
use occt_core::precision::{INFINITE, Precision};

use crate::curve::Curve2d;
use super::curve_tool::{self, GeomAbsCurveType2d};

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
pub struct GFuncExtPC2d<'a> {
    p: GpPnt2d,
    c: Option<&'a dyn Curve2d>,
    u: f64,
    pc: GpPnt2d,
    d1f: f64,
    sq_dist: Vec<f64>,
    is_min: Vec<i32>,
    point: Vec<POnCurv2d>,
    p_init: bool,
    c_init: bool,
    d1_init: bool,
    tol: f64,
    max_deriv_order: i32,
    u_infium: f64,
    u_supremum: f64,
}

impl<'a> Default for GFuncExtPC2d<'a> {
    /// `Extrema_GFuncExtPC()` (`hxx:59-71`).
    fn default() -> Self {
        Self {
            p: GpPnt2d::new(0.0, 0.0),
            c: None,
            u: 0.0,
            pc: GpPnt2d::new(0.0, 0.0),
            d1f: 0.0,
            sq_dist: Vec::new(),
            is_min: Vec::new(),
            point: Vec::new(),
            p_init: false,
            c_init: false,
            d1_init: false,
            tol: MIN_TOL,
            max_deriv_order: 0,
            u_infium: 0.0,
            u_supremum: 0.0,
        }
    }
}

impl<'a> GFuncExtPC2d<'a> {
    /// `Initialize(theC)` (`hxx:105-129`).
    pub fn initialize(&mut self, c: &'a dyn Curve2d) {
        self.c = Some(c);
        self.c_init = true;
        self.point.clear();
        self.sq_dist.clear();
        self.is_min.clear();
        self.sub_interval_initialize(curve_tool::geom2d_first_parameter(c), curve_tool::geom2d_last_parameter(c));
        match curve_tool::geom2d_curve_type(c) {
            GeomAbsCurveType2d::BezierCurve
            | GeomAbsCurveType2d::BSplineCurve
            | GeomAbsCurveType2d::OffsetCurve
            | GeomAbsCurveType2d::OtherCurve => {
                self.max_deriv_order = MAX_ORDER;
                self.tol = self.search_of_tolerance();
            }
            _ => {
                self.max_deriv_order = 0;
                self.tol = MIN_TOL;
            }
        }
    }

    /// `SetPoint(theP)` (`hxx:133-140`).
    pub fn set_point(&mut self, p: &GpPnt2d) {
        self.p = *p;
        self.p_init = true;
        self.point.clear();
        self.sq_dist.clear();
        self.is_min.clear();
    }

    /// `SubIntervalInitialize` (`hxx:420-424`).
    pub fn sub_interval_initialize(&mut self, u_first: f64, u_last: f64) {
        self.u_infium = u_first;
        self.u_supremum = u_last;
    }

    /// `SearchOfTolerance` (`hxx:428-457`).
    fn search_of_tolerance(&self) -> f64 {
        const N_POINT: i32 = 10;
        let step = (self.u_supremum - self.u_infium) / f64::from(N_POINT);
        let c = self.c.expect("curve not initialized");
        let mut num = 0;
        let mut a_max = -INFINITE;
        loop {
            let mut u = self.u_infium + f64::from(num) * step;
            if u > self.u_supremum {
                u = self.u_supremum;
            }
            let (_p, v_der) = curve_tool::geom2d_d1(c, u);
            if !Precision::is_infinite(v_der.x()) && !Precision::is_infinite(v_der.y()) {
                let vm = v_der.magnitude();
                if vm > a_max {
                    a_max = vm;
                }
            }
            num += 1;
            if num >= N_POINT + 1 {
                break;
            }
        }
        (a_max * TOL_FACTOR).max(MIN_TOL)
    }

    /// `Value` (`hxx:146-253`).
    fn value_impl(&mut self, the_u: f64, f: &mut f64) -> bool {
        let c = match self.c {
            Some(c) => c,
            None => return false,
        };
        self.u = the_u;
        let (pc, mut d1c) = curve_tool::geom2d_d1(c, self.u);
        self.pc = pc;
        if Precision::is_infinite(d1c.x()) || Precision::is_infinite(d1c.y()) {
            *f = INFINITE;
            return false;
        }
        let mut ndu = d1c.magnitude();
        if self.max_deriv_order != 0 && ndu <= self.tol {
            const DIVISION_FACTOR: f64 = 1.0e-3;
            let du = if self.u_supremum >= f64::MAX || self.u_infium <= -f64::MAX {
                0.0
            } else {
                self.u_supremum - self.u_infium
            };
            let delta = (du * DIVISION_FACTOR).max(MIN_STEP);
            let mut n = 1i32;
            let mut v = GpVec2d::new(0.0, 0.0);
            let mut is_derive_found = false;
            loop {
                n += 1;
                v = curve_tool::geom2d_dn(c, self.u, n);
                ndu = v.magnitude();
                is_derive_found = ndu > self.tol;
                if is_derive_found || n >= self.max_deriv_order {
                    break;
                }
            }
            if is_derive_found {
                let u = if self.u - self.u_infium < delta {
                    self.u + delta
                } else {
                    self.u - delta
                };
                let p1 = curve_tool::geom2d_d0(c, self.u.min(u));
                let p2 = curve_tool::geom2d_d0(c, self.u.max(u));
                let v1 = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
                let dir_factor = v.dot(&v1);
                d1c = if dir_factor < 0.0 {
                    GpVec2d::new(-v.x(), -v.y())
                } else {
                    v
                };
            } else {
                let (p1, p2, p3, grown);
                if self.u - self.u_infium < 2.0 * delta {
                    p1 = curve_tool::geom2d_d0(c, self.u);
                    p2 = curve_tool::geom2d_d0(c, self.u + delta);
                    p3 = curve_tool::geom2d_d0(c, self.u + 2.0 * delta);
                    grown = true;
                } else {
                    p1 = curve_tool::geom2d_d0(c, self.u - 2.0 * delta);
                    p2 = curve_tool::geom2d_d0(c, self.u - delta);
                    p3 = curve_tool::geom2d_d0(c, self.u);
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
        let ppc = GpVec2d::new(self.pc.x() - self.p.x(), self.pc.y() - self.p.y());
        *f = ppc.dot(&d1c) / ndu;
        true
    }

    /// `Values` (`hxx:274-353`).
    fn values_impl(&mut self, the_u: f64, f: &mut f64, df: &mut f64) -> bool {
        let c = match self.c {
            Some(c) => c,
            None => return false,
        };
        let pc_old = self.pc;
        let p_old = self.p;
        if !self.value_impl(the_u, f) {
            self.d1_init = false;
            return false;
        }
        self.u = the_u;
        self.pc = pc_old;
        self.p = p_old;
        let (pc, d1c, d2c) = curve_tool::geom2d_d2(c, self.u);
        self.pc = pc;
        let ndu = d1c.magnitude();
        if ndu <= self.tol {
            const DIVISION_FACTOR: f64 = 0.01;
            let du = if self.u_supremum >= f64::MAX || self.u_infium <= -f64::MAX {
                0.0
            } else {
                self.u_supremum - self.u_infium
            };
            let delta = (du * DIVISION_FACTOR).max(MIN_STEP);
            if self.u - self.u_infium < 2.0 * delta {
                let f1 = *f;
                let u2 = self.u + delta;
                let u3 = self.u + delta * 2.0;
                let (mut f2, mut f3) = (0.0, 0.0);
                if !self.value_impl(u2, &mut f2) || !self.value_impl(u3, &mut f3) {
                    self.d1_init = false;
                    return false;
                }
                *df = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * delta);
            } else {
                let f3 = *f;
                let u1 = self.u - delta * 2.0;
                let u2 = self.u - delta;
                let (mut f1, mut f2) = (0.0, 0.0);
                if !self.value_impl(u2, &mut f2) || !self.value_impl(u1, &mut f1) {
                    self.d1_init = false;
                    return false;
                }
                *df = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * delta);
            }
            self.u = the_u;
            self.pc = pc_old;
            self.p = p_old;
        } else {
            let ppc = GpVec2d::new(self.p.x() - self.pc.x(), self.p.y() - self.pc.y());
            *df = ndu + (ppc.dot(&d2c) / ndu) - *f * (d1c.dot(&d2c)) / (ndu * ndu);
        }
        self.d1f = *df;
        self.d1_init = true;
        true
    }

    /// `GetStateNumber` (`hxx:357-379`).
    fn get_state_number_impl(&mut self) -> i32 {
        self.sq_dist.push(self.pc.square_distance(&self.p));
        self.d1_init = true;
        let (mut ff, mut dd) = (0.0, 0.0);
        self.values_impl(self.u, &mut ff, &mut dd);
        let int_val = if self.d1f > 0.0 { 1 } else { 0 };
        self.is_min.push(int_val);
        self.point.push(POnCurv2d::new(self.u, self.pc));
        0
    }

    /// `NbExt` (`hxx:382`).
    pub fn nb_ext(&self) -> usize {
        self.sq_dist.len()
    }

    /// `SquareDistance` (`hxx:386-393`), 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_dist[n - 1]
    }

    /// `IsMin` (`hxx:397-404`), 1-based.
    pub fn is_min(&self, n: usize) -> bool {
        self.is_min[n - 1] == 1
    }

    /// `Point` (`hxx:408-415`), 1-based.
    pub fn point(&self, n: usize) -> POnCurv2d {
        self.point[n - 1]
    }
}

impl FunctionWithDerivative for GFuncExtPC2d<'_> {
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
