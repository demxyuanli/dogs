//! `ProjLib_Function` — the `AppCont_Function` adapter handed to
//! `Approx_FitAndDivide2d`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_ComputeApprox.cxx:1026-1117`
//! (OCCT 8.0.0).
//!
//! The C++ class holds the 3d curve, the surface, the parameter window and the
//! seam flags computed by `Function_SetUVBounds`, and answers `Value` / `D1`
//! through the file statics [`function_value`] / [`function_d1`]. `myNbPnt` is
//! `0` and `myNbPnt2d` is `1` (`cxx:1048-1049`): the function carries one 2d
//! point and no 3d one.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use super::function::{function_d1, function_value};
use super::uv_bounds::{function_set_uv_bounds, UvBounds};
use crate::app_cont::function::Function;
use crate::curve::Curve;
use crate::surface::Surface;

/// `class ProjLib_Function` (`cxx:1026-1117`).
pub struct ProjLibFunction<'a> {
    /// `myCurve` (`cxx:1028`).
    curve: &'a dyn Curve,
    /// `mySurface` (`cxx:1029`).
    surface: &'a dyn Surface,
    /// `myIsPeriodic[2]` (`cxx:1030`).
    is_periodic: [bool; 2],
    /// `myPeriod[2]` (`cxx:1031`).
    period: [f64; 2],
    /// The four public members `myU1` / `myU2` / `myV1` / `myV2` plus
    /// `UCouture` / `VCouture` (`cxx:1034-1035`).
    bounds: UvBounds,
}

impl<'a> ProjLibFunction<'a> {
    /// `ProjLib_Function(C, S)` (`cxx:1037-1070`).
    pub fn new(curve: &'a dyn Curve, surface: &'a dyn Surface) -> Self {
        let bounds = function_set_uv_bounds(curve, surface);
        let is_periodic = [surface.is_u_periodic(), surface.is_v_periodic()];
        let mut period = [0.0f64; 2];
        if is_periodic[0] {
            period[0] = surface.u_period();
        }
        if is_periodic[1] {
            period[1] = surface.v_period();
        }
        Self {
            curve,
            surface,
            is_periodic,
            period,
            bounds,
        }
    }

    /// The `myU1` / `myU2` / `myV1` / `myV2` / `UCouture` / `VCouture` members
    /// (`cxx:1034-1035`), read by `Perform`'s "return curve home" tail.
    pub fn bounds(&self) -> &UvBounds {
        &self.bounds
    }

    /// `gp_Pnt2d Value(const double theT)` (`cxx:1091-1094`).
    pub fn value_2d(&self, the_t: f64) -> GpPnt2d {
        let b = &self.bounds;
        function_value(
            the_t,
            self.curve,
            self.surface,
            b.u1,
            b.u2,
            b.v1,
            b.v2,
            b.u_couture,
            b.v_couture,
        )
    }
}

impl Function for ProjLibFunction<'_> {
    /// `myNbPnt = 0`, `myNbPnt2d = 1` (`cxx:1048-1049`).
    fn get_number_of_points(&self) -> (i32, i32) {
        (0, 1)
    }

    /// `FirstParameter()` (`cxx:1079`).
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    /// `LastParameter()` (`cxx:1081`).
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// `Value(theT, thePnt2d, thePnt)` (`cxx:1083-1089`).
    fn value(&self, the_t: f64, the_pnt2d: &mut [GpPnt2d], _the_pnt: &mut [GpPnt]) -> bool {
        if let Some(slot) = the_pnt2d.first_mut() {
            *slot = self.value_2d(the_t);
        }
        true
    }

    /// `D1(theT, theVec2d, theVec)` (`cxx:1096-1115`): `theVec2d(1)` is written
    /// unconditionally and carries the zero vector when `Function_D1` fails.
    fn d1(&self, the_t: f64, the_vec2d: &mut [GpVec2d], _the_vec: &mut [GpVec]) -> bool {
        let b = &self.bounds;
        let result = function_d1(
            the_t,
            self.curve,
            self.surface,
            b.u1,
            b.u2,
            b.v1,
            b.v2,
            b.u_couture,
            b.v_couture,
        );
        let (is_ok, vec2d) = match result {
            Some((_, v)) => (true, v),
            None => (false, GpVec2d::new(0.0, 0.0)),
        };
        if let Some(slot) = the_vec2d.first_mut() {
            *slot = vec2d;
        }
        is_ok
    }

    /// `PeriodInformation(theDimIdx, IsPeriodic, thePeriod)` (`cxx:1072-1076`).
    fn period_information(&self, the_dim_idx: i32) -> (bool, f64) {
        let i = (the_dim_idx - 1) as usize;
        (self.is_periodic[i], self.period[i])
    }
}
