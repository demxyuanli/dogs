//! `Approx_CurvilinearParameter` case 2
//! (`Approx_CurvilinearParameter.cxx:285-397`): reparametrize a curve-on-surface
//! by arc length through `Approx_CurvlinFunc` + `AdvApprox_ApproxAFunction`.

use crate::approx_curvlin_func::CurvlinFunc;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_core::adv_approx::{ApproxAFunction1d3d, PrefAndRec};
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use std::cell::RefCell;
use std::sync::Arc;

/// Result of the case-2 constructor used by `BRepLib::SameParameter`.
pub struct CurvilinearParameter {
    pub done: bool,
    pub has_result: bool,
    pub curve2d1: Option<Geom2dBSplineCurve>,
    pub max_error_2d1: f64,
    pub max_error_3d: f64,
}

impl CurvilinearParameter {
    /// `Approx_CurvilinearParameter(C2D, Surf, Tol, Order, MaxDegree, MaxSegments)`
    /// (`cxx:285-397`).
    ///
    /// `cos` is the matching `Adaptor3d_CurveOnSurface(C2D, Surf)`.
    /// `rebuild_cos` rebuilds that Cos after `CurvlinFunc::Trim`.
    pub fn from_curve_on_surface<F>(
        c2d: Arc<dyn Curve2d>,
        surf: Arc<dyn Surface>,
        cos: Arc<dyn Curve>,
        tol: f64,
        order: Shape,
        max_degree: i32,
        max_segments: i32,
        rebuild_cos: F,
    ) -> Self
    where
        F: Fn(Arc<dyn Curve2d>, Arc<dyn Surface>) -> Arc<dyn Curve>,
    {
        let fonct = RefCell::new(CurvlinFunc::from_curve_on_surface(
            c2d,
            surf,
            cos,
            tol / 20.0,
        ));
        let first_s = fonct.borrow().first_parameter();
        let last_s = fonct.borrow().last_parameter();
        let cut_c2 = fonct.borrow().intervals(Shape::C2);
        let cut_c3 = fonct.borrow().intervals(Shape::C3);
        let cut_tool = match PrefAndRec::new(cut_c2, cut_c3, 5.0) {
            Ok(c) => c,
            Err(()) => return Self::empty(),
        };
        let rebuild = RefCell::new(rebuild_cos);
        let eval = |start: f64, end: f64, param: f64, der_order: i32, result: &mut [f64]| -> i32 {
            if result.len() < 5 {
                return 1;
            }
            let mut f = fonct.borrow_mut();
            if (start - f.first_parameter()).abs() > CONFUSION
                || (end - f.last_parameter()).abs() > CONFUSION
            {
                f.trim(start, end, CONFUSION, |pc, s| rebuild.borrow()(pc, s));
            }
            if param < start || param > end {
                return 2;
            }
            if !f.eval_case2(param, der_order, result) {
                return 3;
            }
            0
        };
        let approx = match ApproxAFunction1d3d::approx(
            first_s,
            last_s,
            order,
            max_degree,
            max_segments,
            tol,
            tol / 2.0,
            &cut_tool,
            &eval,
        ) {
            Ok(a) => a,
            Err(()) => return Self::empty(),
        };
        let curve2d1 = if approx.has_result && !approx.poles_1d_u.is_empty() {
            let degree = approx.degree as usize;
            let flat = knot_sequence(&approx.knots, &approx.mults, degree as i32);
            Geom2dBSplineCurve::from_flat(
                approx.poles_1d_u,
                approx.poles_1d_v,
                None,
                flat,
                degree,
                false,
            )
            .ok()
        } else {
            None
        };
        Self {
            done: approx.done,
            has_result: approx.has_result,
            curve2d1,
            max_error_2d1: approx.max_error_1d,
            max_error_3d: approx.max_error_3d,
        }
    }

    fn empty() -> Self {
        Self {
            done: false,
            has_result: false,
            curve2d1: None,
            max_error_2d1: 0.0,
            max_error_3d: 0.0,
        }
    }
}
