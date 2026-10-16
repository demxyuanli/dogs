//! Port of `Extrema_GFuncExtPC` (`Extrema_GFuncExtPC.hxx:32-482`)
//! monomorphised onto the 3D instantiation
//! `Extrema_GFuncExtPC<Adaptor3d_Curve, Extrema_CurveTool, Extrema_POnCurv,
//! gp_Pnt, gp_Vec, NCollection_Sequence<Extrema_POnCurv>>`
//! (`Extrema_PCLocFOfLocEPCOfLocateExtPC.hxx:29-35`).
//!
//! This is the `math_FunctionWithDerivative` whose zero is searched by
//! `Extrema_GenLocateExtPC` (`Extrema_GenLocateExtPC.hxx:108-126`), and hence
//! by `Extrema_LocateExtPC` (`Extrema_LocateExtPC.hxx:29-34`).
//!
//! The 2D twin is `occt_geom2d::geom2d_int::gfunc_ext_pc`.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{INFINITE, Precision};
use occt_math::MathFunctionWithDerivative;

use crate::curve::Curve;
use crate::extrema_cc::POnCurv;

use super::curve_tool::{self, GeomAbsCurveType};

const TOL_FACTOR: f64 = 1.0e-12;
const MIN_TOL: f64 = 1.0e-20;
const MIN_STEP: f64 = 1.0e-7;
const MAX_ORDER: i32 = 3;

/// `RealLast()` (`Standard_Real.hxx:166-169`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx:175-178`).
const REAL_FIRST: f64 = -f64::MAX;

/// `Extrema_GFuncExtPC` (`hxx:53-480`), 3D instantiation.
pub struct GFuncExtPC<'a> {
    my_p: GpPnt,
    my_c: Option<&'a dyn Curve>,
    my_u: f64,
    my_pc: GpPnt,
    my_d1f: f64,
    my_sq_dist: Vec<f64>,
    my_is_min: Vec<i32>,
    my_point: Vec<POnCurv>,
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
            my_p: GpPnt::zero(),
            my_c: None,
            my_u: 0.0,
            my_pc: GpPnt::zero(),
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
    pub fn initialize(&mut self, c: &'a dyn Curve) {
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
    pub fn set_point(&mut self, p: &GpPnt) {
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
            let du = if self.my_u_supremum >= REAL_LAST || self.my_u_infium <= REAL_FIRST {
                0.0
            } else {
                self.my_u_supremum - self.my_u_infium
            };
            let a_delta = (du * DIVISION_FACTOR).max(MIN_STEP);
            let mut n = 1i32;
            // `TheVector V; bool IsDeriveFound;` (`hxx:180-181`): both are
            // declared uninitialised and assigned by the first `do` iteration.
            let (mut v, mut is_derive_found);
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
                let v1 = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                let a_dir_factor = v.dot(&v1);
                d1c = if a_dir_factor < 0.0 {
                    GpVec::new(-v.x(), -v.y(), -v.z())
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
                // `gp_Vec(Ptemp, P1)` with the default-constructed `Ptemp` is
                // the point taken as a vector from the origin.
                let (v1, v2, v3) = (
                    GpVec::new(p1.x(), p1.y(), p1.z()),
                    GpVec::new(p2.x(), p2.y(), p2.z()),
                    GpVec::new(p3.x(), p3.y(), p3.z()),
                );
                d1c = if grown {
                    GpVec::new(
                        -3.0 * v1.x() + 4.0 * v2.x() - v3.x(),
                        -3.0 * v1.y() + 4.0 * v2.y() - v3.y(),
                        -3.0 * v1.z() + 4.0 * v2.z() - v3.z(),
                    )
                } else {
                    GpVec::new(
                        v1.x() - 4.0 * v2.x() + 3.0 * v3.x(),
                        v1.y() - 4.0 * v2.y() + 3.0 * v3.y(),
                        v1.z() - 4.0 * v2.z() + 3.0 * v3.z(),
                    )
                };
            }
            ndu = d1c.magnitude();
        }
        if ndu <= MIN_TOL {
            return false;
        }
        let ppc = GpVec::new(
            self.my_pc.x() - self.my_p.x(),
            self.my_pc.y() - self.my_p.y(),
            self.my_pc.z() - self.my_p.z(),
        );
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
            let du = if self.my_u_supremum >= REAL_LAST || self.my_u_infium <= REAL_FIRST {
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
            let ppc = GpVec::new(
                self.my_p.x() - self.my_pc.x(),
                self.my_p.y() - self.my_pc.y(),
                self.my_p.z() - self.my_pc.z(),
            );
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
        self.my_point.push(POnCurv::new(self.my_u, self.my_pc));
        0
    }

    /// `NbExt` (`hxx:382`). The only cxx reader is
    /// `Extrema_GLocateExtPC.hxx:287`, which calls `NbExt()` on `myExtremPC`,
    /// i.e. `Extrema_ExtPC`, not on this class; the ported callers of this port
    /// use [`Self::point`] / [`Self::square_distance`] directly, so the accessor
    /// is kept for the member set and marked unused.
    #[allow(dead_code)]
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
    pub fn point(&self, the_n: usize) -> POnCurv {
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
