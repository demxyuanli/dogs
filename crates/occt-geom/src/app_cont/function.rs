//! `AppCont_Function` - the continuous 3d/2d vector function handed to the
//! continuous approximation algorithms.
//!
//! Source: `ModelingData/TKGeomBase/AppCont/AppCont_Function.hxx:26-81`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Array1<gp_Pnt2d>&` / `NCollection_Array1<gp_Pnt>&` output
//!   parameters become `&mut [GpPnt2d]` / `&mut [GpPnt]`; OCCT's 1-based
//!   `thePnt2d(1)` is the port's index `0`.
//! - `myNbPnt` / `myNbPnt2d` default to `-1` (`hxx:32-35`) and stay negative
//!   until the concrete function sets them, so the default accessors answer
//!   `-1` exactly like OCCT.
//! - `PeriodInformation` has the OCCT default body (`hxx:70-75`): not periodic,
//!   period `0.0`.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

/// `AppCont_Function` (`AppCont_Function.hxx:26-81`).
pub trait Function {
    /// `GetNumberOfPoints` (`hxx:37-41`): `(myNbPnt, myNbPnt2d)`.
    fn get_number_of_points(&self) -> (i32, i32);

    /// `GetNbOf3dPoints` (`hxx:43-44`).
    fn nb_of_3d_points(&self) -> i32 {
        self.get_number_of_points().0
    }

    /// `GetNbOf2dPoints` (`hxx:46-47`).
    fn nb_of_2d_points(&self) -> i32 {
        self.get_number_of_points().1
    }

    /// `FirstParameter` (`hxx:53-54`).
    fn first_parameter(&self) -> f64;

    /// `LastParameter` (`hxx:57-58`).
    fn last_parameter(&self) -> f64;

    /// `Value(theU, thePnt2d, thePnt)` (`hxx:59-61`).
    ///
    /// `the_pnt2d` has `nb_of_2d_points()` slots, `the_pnt` has
    /// `nb_of_3d_points()` slots.
    fn value(&self, the_u: f64, the_pnt2d: &mut [GpPnt2d], the_pnt: &mut [GpPnt]) -> bool;

    /// `D1(theU, theVec2d, theVec)` (`hxx:63-65`).
    fn d1(&self, the_u: f64, the_vec2d: &mut [GpVec2d], the_vec: &mut [GpVec]) -> bool;

    /// `PeriodInformation(theDimIdx, IsPeriodic, thePeriod)` (`hxx:68-75`).
    /// `the_dim_idx` is 1-based over `3 * nb_of_3d_points() + 2 * nb_2d`.
    fn period_information(&self, _the_dim_idx: i32) -> (bool, f64) {
        (false, 0.0)
    }
}
