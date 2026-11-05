//! The null `Geom2d_Curve` handle.
//!
//! `Bisector_BisecPC` (`cxx:36-50`) and `Bisector_BisecCC` leave their curve
//! members as null handles after the default constructor. OCCT would
//! dereference those null handles; the port stores this placeholder instead so
//! that the default-constructed objects answer every query with a zero value.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};

use crate::curve::Curve2d;

/// Placeholder for a null `occ::handle<Geom2d_Curve>`.
#[derive(Debug, Clone, Copy)]
pub struct NullCurve;

impl Curve2d for NullCurve {
    fn d0(&self, _u: f64) -> GpPnt2d {
        GpPnt2d::zero()
    }

    fn d1(&self, _u: f64) -> (GpPnt2d, GpVec2d) {
        (GpPnt2d::zero(), GpVec2d::zero())
    }

    fn d2(&self, _u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        (GpPnt2d::zero(), GpVec2d::zero(), GpVec2d::zero())
    }

    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        0.0
    }

    fn continuity(&self) -> u8 {
        0
    }

    fn transform(&mut self, _t: &GpTrsf2d) {}

    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(*self)
    }
}

/// A fresh null handle.
pub fn null_curve() -> Arc<dyn Curve2d> {
    Arc::new(NullCurve)
}
