//! `LocOpe_LinearForm` (`LocOpe_LinearForm.hxx:25-84`, `.cxx:35-172`).
//!
//! The feature front end of a translated prism. `IntPerf`
//! (`LocOpe_LinearForm.cxx:110-229`) is a verbatim copy of
//! `LocOpe_Prism::IntPerf` (`LocOpe_Prism.cxx:107-226`) - same sweep, same
//! `Face` / non-`Face` split, same `LocOpe_BuildShape` fallback, same
//! descendant remap after a translation and the same
//! `BRepTools_Modifier` round trip. The port therefore drives one
//! [`super::prism::Prism`] body instead of duplicating that control flow, and
//! keeps the two points the reference stores but never reads (`myPnt1`,
//! `myPnt2`, `hxx:80-81`).
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | `Perform` (`cxx:35-88`) | this file |
//! | `IntPerf` (`cxx:110-229`) | [`super::prism::Prism::perform`] / [`super::prism::Prism::perform_with_translation`] |
//! | `Shape` / `FirstShape` / `LastShape` / `Shapes` (`cxx:...`) | accessors |
//!
//! Deviations: `myRes` is `Option<TopoShape>` and `Shape()` answers `None`
//! where OCCT raises `StdFail_NotDone`; `Shapes(S)` answers `None` where the
//! reference's `myMap(S)` raises on an unbound key.

use occt_core::gp::{GpPnt, GpVec};

use crate::locope::prism::Prism;
use crate::shape::TopoShape;

/// `LocOpe_LinearForm` (`hxx:38-84`).
#[derive(Debug, Clone, Default)]
pub struct LinearForm {
    /// `TopoDS_Shape myBase`
    my_base: Option<TopoShape>,
    /// `gp_Vec myVec`
    my_vec: GpVec,
    /// `gp_Vec myTra`
    my_tra: GpVec,
    /// `bool myDone`
    my_done: bool,
    /// `bool myIsTrans`
    my_is_trans: bool,
    /// `gp_Pnt myPnt1` - stored, never read (`hxx:80`).
    my_pnt1: GpPnt,
    /// `gp_Pnt myPnt2` - stored, never read (`hxx:81`).
    my_pnt2: GpPnt,
    /// The shared `IntPerf` body (`LocOpe_Prism.cxx:107-226`).
    prism: Prism,
}

impl LinearForm {
    /// `LocOpe_LinearForm()` (`hxx:42`).
    pub fn new() -> Self {
        Self {
            my_pnt1: GpPnt::new(0.0, 0.0, 0.0),
            my_pnt2: GpPnt::new(0.0, 0.0, 0.0),
            ..Default::default()
        }
    }

    /// `Perform(Base, V, Pnt1, Pnt2)` (`cxx:35-58`).
    pub fn perform(&mut self, base: &TopoShape, v: &GpVec, pnt1: &GpPnt, pnt2: &GpPnt) {
        self.my_is_trans = false;
        self.my_base = Some(base.clone());
        self.my_vec = *v;
        self.my_pnt1 = *pnt1;
        self.my_pnt2 = *pnt2;
        self.prism.perform(base, v);
        self.my_done = self.prism.is_done();
    }

    /// `Perform(Base, V, Vectra, Pnt1, Pnt2)` (`cxx:64-88`).
    pub fn perform_with_translation(
        &mut self,
        base: &TopoShape,
        v: &GpVec,
        vectra: &GpVec,
        pnt1: &GpPnt,
        pnt2: &GpPnt,
    ) {
        self.my_is_trans = true;
        self.my_tra = *vectra;
        self.my_base = Some(base.clone());
        self.my_vec = *v;
        self.my_pnt1 = *pnt1;
        self.my_pnt2 = *pnt2;
        self.prism.perform_with_translation(base, v, vectra);
        self.my_done = self.prism.is_done();
    }

    /// `Shape()`; `None` is the `StdFail_NotDone` raise (`cxx:...`).
    pub fn shape(&self) -> Option<&TopoShape> {
        if self.my_done {
            self.prism.shape()
        } else {
            None
        }
    }

    /// `FirstShape()`.
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.prism.first_shape()
    }

    /// `LastShape()`.
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.prism.last_shape()
    }

    /// `Shapes(S)`.
    pub fn shapes(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        self.prism.shapes(s)
    }

    /// The two stored points (`hxx:80-81`); kept so the port does not silently
    /// drop reference state.
    pub fn points(&self) -> (&GpPnt, &GpPnt) {
        (&self.my_pnt1, &self.my_pnt2)
    }
}
