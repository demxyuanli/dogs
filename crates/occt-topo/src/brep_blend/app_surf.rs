//! `BRepBlend_AppSurf` — `AppBlend_AppSurf` instantiated with
//! `TheSectionGenerator = Blend_AppFunction` and `TheLine = BRepBlend_Line`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_AppSurf.hxx:44-193`
//! and `BRepBlend_AppSurf_0.cxx:17-30` (OCCT 8.0.0); `BRepBlend_AppSurf.hxx`
//! only sets the template macros, so the instantiated body is
//! `TKGeomAlgo/AppBlend/AppBlend_AppSurf.gxx`.
//!
//! Translation notes:
//! - `AppBlend_AppSurf.gxx` has no member that depends on `TheSectionGenerator`
//!   and no method stores one; it only uses `Lin->Point(i)` to feed
//!   `F.Section(...)` / `F.Parameter(...)`. The port therefore shares
//!   [`occt_geom::geom_fill::AppSurf`] (the `GeomFill_AppSurf` /
//!   `GeomFill_AppSweep` struct) and selects this pairing through the
//!   generator passed to [`perform`] / [`perform_nb_max_p`], exactly as
//!   `AppSweep` does.
//! - `TheLine = BRepBlend_Line` reaches the shared driver as the index list
//!   `GeomFill_Line::with_nb_points(Lin->NbPoints())`, whose `Point(Index)` is
//!   the identity; [`LineSectionGenerator`] then maps that index back to
//!   `Lin->Point(Index)` for every `Section` / `Parameter` call, so the driver
//!   sees the same sequence OCCT hands it.
//! - `PerformSmoothing` keeps the `AppDef_Variational` gap of the shared
//!   driver (`AppBlend_AppSurf.gxx:472-541`).
//! - The `Blend_AppFunction&` that `Perform` receives is non-`const` in OCCT
//!   (`Blend_AppFunction.hxx:158-185`) while `AppBlendSectionGenerator` takes
//!   `&self`; the adapter holds it as `RefCell<&mut dyn BlendAppFunction>`.

use std::cell::RefCell;

use crate::blend::BlendAppFunction;
use crate::brep_blend::Line;

use super::app_func::line_point;

pub use occt_geom::geom_fill::{AppBlendSectionGenerator, AppSurf, SectionValues};

use occt_geom::geom_fill::Line as IndexList;

/// `TheSectionGenerator = Blend_AppFunction`, `TheLine = BRepBlend_Line`
/// (`BRepBlend_AppSurf.hxx:179-187`): the section of one stored `Blend_Point`
/// of the walking line.
struct LineSectionGenerator<'a> {
    line: &'a Line,
    func: RefCell<&'a mut dyn BlendAppFunction>,
}

impl AppBlendSectionGenerator for LineSectionGenerator<'_> {
    /// `Blend_AppFunction::GetShape(NbPoles, NbKnots, Degree, NbPoles2d)`
    /// (`Blend_AppFunction.hxx:135-139`).
    fn get_shape(&self) -> Result<(i32, i32, i32, i32), String> {
        let (mut nb_poles, mut nb_knots, mut degree, mut nb_poles2d) = (0, 0, 0, 0);
        self.func.borrow_mut().get_shape(
            &mut nb_poles,
            &mut nb_knots,
            &mut degree,
            &mut nb_poles2d,
        );
        Ok((nb_poles, nb_knots, degree, nb_poles2d))
    }

    /// `Blend_AppFunction::Knots(TKnots)` (`Blend_AppFunction.hxx:152`).
    fn knots(&self) -> Result<Vec<f64>, String> {
        let (_nb_poles, nb_knots, _degree, _nb_poles2d) = self.get_shape()?;
        let mut tknots = vec![0.0; nb_knots.max(0) as usize];
        self.func.borrow_mut().knots(&mut tknots);
        Ok(tknots)
    }

    /// `Blend_AppFunction::Mults(TMults)` (`Blend_AppFunction.hxx:154`).
    fn mults(&self) -> Result<Vec<i32>, String> {
        let (_nb_poles, nb_knots, _degree, _nb_poles2d) = self.get_shape()?;
        let mut tmults = vec![0; nb_knots.max(0) as usize];
        self.func.borrow_mut().mults(&mut tmults);
        Ok(tmults)
    }

    /// `Blend_AppFunction::Section(P, Poles, DPoles, Poles2d, DPoles2d,
    /// Weigths, DWeigths)` (`Blend_AppFunction.hxx:158-164`); the `bool` is
    /// OCCT's `withderiv`.
    fn section_with_derivatives(&self, p: i32) -> Result<(SectionValues, bool), String> {
        let pnt = line_point(self.line, p as usize);
        let (nb_poles, _nb_knots, _degree, nb_poles2d) = self.get_shape()?;
        let nb_poles = nb_poles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut poles = vec![occt_core::gp::GpPnt::zero(); nb_poles];
        let mut dpoles = vec![occt_core::gp::GpVec::zero(); nb_poles];
        let mut poles2d = vec![occt_core::gp::GpPnt2d::zero(); nb_poles2d];
        let mut dpoles2d = vec![occt_core::gp::GpVec2d::zero(); nb_poles2d];
        let mut weights = vec![0.0; nb_poles];
        let mut dweights = vec![0.0; nb_poles];
        let with_deriv = self.func.borrow_mut().section(
            &pnt,
            &mut poles,
            &mut dpoles,
            &mut poles2d,
            &mut dpoles2d,
            &mut weights,
            &mut dweights,
        );
        Ok((
            SectionValues {
                poles,
                d_poles: dpoles,
                poles2d,
                d_poles2d: dpoles2d,
                weights,
                d_weights: dweights,
            },
            with_deriv,
        ))
    }

    /// `Blend_AppFunction::Section(P, Poles, Poles2d, Weigths)`
    /// (`Blend_AppFunction.hxx:166-169`): the `SpApprox` overload.
    fn section(&self, p: i32) -> Result<SectionValues, String> {
        let pnt = line_point(self.line, p as usize);
        let (nb_poles, _nb_knots, _degree, nb_poles2d) = self.get_shape()?;
        let nb_poles = nb_poles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut poles = vec![occt_core::gp::GpPnt::zero(); nb_poles];
        let mut poles2d = vec![occt_core::gp::GpPnt2d::zero(); nb_poles2d];
        let mut weights = vec![0.0; nb_poles];
        self.func
            .borrow_mut()
            .section_poles(&pnt, &mut poles, &mut poles2d, &mut weights);
        Ok(SectionValues {
            poles,
            poles2d,
            weights,
            ..SectionValues::default()
        })
    }

    /// `Blend_AppFunction::Parameter(P)` (`Blend_AppFunction.cxx:20-23`).
    fn parameter(&self, p: i32) -> Result<f64, String> {
        let pnt = line_point(self.line, p as usize);
        Ok(self.func.borrow().parameter(&pnt))
    }
}

/// `Perform(Lin, SecGen, SpApprox)` (`BRepBlend_AppSurf.hxx:86-88`).
pub fn perform(
    app: &mut AppSurf,
    line: &Line,
    func: &mut dyn BlendAppFunction,
    sp_approx: bool,
) -> Result<(), String> {
    let lin = IndexList::with_nb_points(line.nb_points() as i32);
    let generator = LineSectionGenerator {
        line,
        func: RefCell::new(func),
    };
    app.perform(&lin, &generator, sp_approx)
}

/// `PerformSmoothing(Lin, SecGen)` (`BRepBlend_AppSurf.hxx:91-92`).
pub fn perform_smoothing(
    app: &mut AppSurf,
    line: &Line,
    func: &mut dyn BlendAppFunction,
) -> Result<(), String> {
    let lin = IndexList::with_nb_points(line.nb_points() as i32);
    let generator = LineSectionGenerator {
        line,
        func: RefCell::new(func),
    };
    app.perform_smoothing(&lin, &generator)
}

/// `Perform(Lin, SecGen, NbMaxP)` (`BRepBlend_AppSurf.hxx:94-96`).
pub fn perform_nb_max_p(
    app: &mut AppSurf,
    line: &Line,
    func: &mut dyn BlendAppFunction,
    nb_max_p: i32,
) -> Result<(), String> {
    let lin = IndexList::with_nb_points(line.nb_points() as i32);
    let generator = LineSectionGenerator {
        line,
        func: RefCell::new(func),
    };
    app.perform_nb_max_p(&lin, &generator, nb_max_p)
}
