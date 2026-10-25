//! The face / face half of the `ChFi3d_Builder` protected virtual interface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.hxx` (OCCT 8.0.0):
//!
//! - `PerformSurf` (`hxx:445-466`) — pure virtual (`= 0`), the *construction*
//!   of the fillet / chamfer surface between two faces,
//! - `SimulSurf` (`hxx:203-221`) — pure virtual (`= 0`), the *simulation* of
//!   the same band, answering a single `ChFiDS_SurfData`.
//!
//! OCCT leaves both to `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`, so the base
//! class cannot define them. The port states the same thing with a trait: the
//! subclasses are not ported yet, so [`ChFi3dSurfOps`] has no implementor in
//! the crate, and the callers that need it ([`super::call_perform_surf`]) take
//! `&mut dyn ChFi3dSurfOps`.
//!
//! Divergences that are not heuristics:
//! - `occ::handle<ChFiDS_SurfData>&` is `&mut SurfData` (`SimulSurf`) and
//!   `NCollection_Sequence<occ::handle<ChFiDS_SurfData>>&` is
//!   `&mut Vec<Arc<SurfData>>` (`PerformSurf`), matching
//!   [`crate::chfi_ds::Stripe::set_of_surf_data`],
//! - `occ::handle<ChFiDS_ElSpine>&` is `&mut Arc<ElSpine>`: the two virtuals
//!   onward call `ComputeData`, which writes the guideline through
//!   `Arc::make_mut` (see `builder/compute_data_face_face/mod.rs`),
//! - `occ::handle<Adaptor3d_TopolTool>&` is `Rc<RefCell<TopolTool>>`, the
//!   crate's shared-mutable topology tool,
//! - `BRepAdaptor_Surface&` is `&SurfaceAdaptor`; the port's adaptor carries
//!   its `TopoDS_Face`, which is what `S1->Face()` answers in OCCT.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_math::MathVector;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{ElSpine, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// The protected virtuals of `ChFi3d_Builder` that `ChFi3d_FilBuilder` /
/// `ChFi3d_ChBuilder` override (`hxx:203-221`, `:445-466`).
///
/// `data` / `guide` / `spine` / `soldep` mirror the OCCT references; `first` /
/// `last` / `intf` / `intl` / `pos1` / `pos2`-style outputs are `&mut` the same
/// way, so an implementor cannot silently drop an output.
pub trait ChFi3dSurfOps {
    /// `PerformSurf(Data, Guide, Spine, Choix, S1, I1, S2, I2, MaxStep, Fleche,
    /// TolGuide, First, Last, Inside, Appro, Forward, RecOnS1, RecOnS2,
    /// Soldep, Intf, Intl)` (`hxx:445-466`).
    ///
    /// `data` is `NCollection_Sequence<handle(ChFiDS_SurfData)>`, `false` is
    /// OCCT's "did not compute".
    #[allow(clippy::too_many_arguments)]
    fn perform_surf(
        &mut self,
        data: &mut Vec<Arc<SurfData>>,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool;

    /// `SimulSurf(Data, Guide, Spine, Choix, S1, I1, S2, I2, TolGuide, First,
    /// Last, Inside, Appro, Forward, RecOnS1, RecOnS2, Soldep, Intf, Intl)`
    /// (`hxx:203-221`).
    ///
    /// `false` is OCCT's "did not simulate". `SimulKPart` and the three
    /// edge-contact overloads of both functions (`hxx:223-302`, `:399-540`) are
    /// not part of this trait yet; see the module UNPORTED list.
    #[allow(clippy::too_many_arguments)]
    fn simul_surf(
        &mut self,
        data: &mut SurfData,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool;
}
