//! The face / face half of the `ChFi3d_Builder` protected virtual interface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.hxx` (OCCT 8.0.0):
//!
//! - `PerformSurf` (`hxx:445-466`) — pure virtual (`= 0`), the *construction*
//!   of the fillet / chamfer surface between two faces,
//! - `SimulSurf` (`hxx:203-221`) — pure virtual (`= 0`), the *simulation* of
//!   the same band, answering a single `ChFiDS_SurfData`,
//! - the six edge-contact overloads of the two (`hxx:223-302`, `hxx:470-556`)
//!   — `virtual` with a base body that throws
//!   (`ChFi3d_Builder_NotImp.cxx:29-220`); they are trait defaults returning
//!   those `Err`s here.
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

use crate::abs::{Orientation, State};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
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
    /// `false` is OCCT's "did not simulate". The six edge-contact overloads of
    /// both functions are the throwing defaults further down; `SimulKPart`
    /// (`hxx:201`) is still missing.
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

    /// `PerformFirstSection(S, HGuide, Choix, S1, S2, I1, I2, Par, SolDep, Pos1,
    /// Pos2)` (`hxx:577-588`) — pure virtual (`= 0`), the initial-section
    /// construction `ChFi3d_Builder::StartSol` dispatches to.
    ///
    /// `par` is the guideline parameter the caller starts from; the subclass
    /// answers the four `SolDep` coordinates and the two `TopAbs_State`s.
    /// OCCT declares the method `const`; the port keeps `&mut self` so a caller
    /// can hold one `&mut dyn ChFi3dSurfOps` for both this and `PerformSurf` /
    /// `SimulSurf`.
    #[allow(clippy::too_many_arguments)]
    fn perform_first_section(
        &mut self,
        spine: &SpineHandle,
        guide: &Arc<ElSpine>,
        choix: i32,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        i2: &Rc<RefCell<TopolTool>>,
        par: f64,
        sol_dep: &mut MathVector,
        pos1: &mut State,
        pos2: &mut State,
    ) -> bool;

    // ------------------------------------------------------------------
    // The six edge-contact overloads: virtual with a base body that throws.
    // ------------------------------------------------------------------
    //
    // `ChFi3d_Builder` declares all six `virtual` without `= 0`
    // (`hxx:223-302`, `hxx:470-556`) and defines them in
    // `ChFi3d_Builder_NotImp.cxx:29-220` as `throw Standard_Failure(...)`. Both
    // subclasses repeat exactly the same throwing bodies
    // (`ChFi3d_ChBuilder.cxx:1223-1310`, `:1854-1940`;
    // `ChFi3d_FilBuilder.cxx:1020-1110` overrides two of them with real work
    // for the fillet case, which is a separate port). Until those are ported,
    // the faithful answer of the *base* class — and of the chamfer builder — is
    // the `Err`, which is what these defaults return.

    /// `PerformSurf(Data, Guide, Spine, Choix, S1, I1, PC1, Sref1, PCref1,
    /// Decroch1, S2, I2, Or2, MaxStep, Fleche, TolGuide, First, Last, Inside,
    /// Appro, Forward, RecP, RecS, RecRst, Soldep)` (`hxx:470-492`): edge / face
    /// contact on `S1`.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:125-152`.
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_on_s1(
        &mut self,
        _data: &mut Vec<Arc<SurfData>>,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _pc1: Option<&Curve2dAdaptor>,
        _sref1: &SurfaceAdaptor,
        _pcref1: &Curve2dAdaptor,
        _decroch1: &mut bool,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _or2: Orientation,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p: bool,
        _rec_s: bool,
        _rec_rst: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("PerformSurf Not Implemented".to_string())
    }

    /// `PerformSurf(Data, Guide, Spine, Choix, S1, I1, Or1, S2, I2, PC2, Sref2,
    /// PCref2, Decroch2, MaxStep, Fleche, TolGuide, First, Last, Inside, Appro,
    /// Forward, RecP, RecS, RecRst, Soldep)` (`hxx:498-522`): edge / face
    /// contact on `S2`.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:156-182`.
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_on_s2(
        &mut self,
        _data: &mut Vec<Arc<SurfData>>,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _or1: Orientation,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _pc2: Option<&Curve2dAdaptor>,
        _sref2: &SurfaceAdaptor,
        _pcref2: &Curve2dAdaptor,
        _decroch2: &mut bool,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p: bool,
        _rec_s: bool,
        _rec_rst: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("PerformSurf Not Implemented".to_string())
    }

    /// `PerformSurf(Data, Guide, Spine, Choix, S1, I1, PC1, Sref1, PCref1,
    /// Decroch1, Or1, S2, I2, PC2, Sref2, PCref2, Decroch2, Or2, MaxStep,
    /// Fleche, TolGuide, First, Last, Inside, Appro, Forward, RecP1, RecRst1,
    /// RecP2, RecRst2, Soldep)` (`hxx:527-556`): edge / edge contact.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:187-220`.
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_contact(
        &mut self,
        _data: &mut Vec<Arc<SurfData>>,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _pc1: Option<&Curve2dAdaptor>,
        _sref1: &SurfaceAdaptor,
        _pcref1: &Curve2dAdaptor,
        _decroch1: &mut bool,
        _or1: Orientation,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _pc2: Option<&Curve2dAdaptor>,
        _sref2: &SurfaceAdaptor,
        _pcref2: &Curve2dAdaptor,
        _decroch2: &mut bool,
        _or2: Orientation,
        _max_step: f64,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p1: bool,
        _rec_rst1: bool,
        _rec_p2: bool,
        _rec_rst2: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("PerformSurf Not Implemented".to_string())
    }

    /// `SimulSurf(Data, Guide, Spine, Choix, S1, I1, PC1, Sref1, PCref1,
    /// Decroch1, S2, I2, Or2, Fleche, TolGuide, First, Last, Inside, Appro,
    /// Forward, RecP, RecS, RecRst, Soldep)` (`hxx:223-247`): edge / face
    /// contact on `S1`.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:29-55`.
    #[allow(clippy::too_many_arguments)]
    fn simul_surf_on_s1(
        &mut self,
        _data: &mut SurfData,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _pc1: Option<&Curve2dAdaptor>,
        _sref1: &SurfaceAdaptor,
        _pcref1: &Curve2dAdaptor,
        _decroch1: &mut bool,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _or2: Orientation,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p: bool,
        _rec_s: bool,
        _rec_rst: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("SimulSurf Not Implemented".to_string())
    }

    /// `SimulSurf(Data, Guide, Spine, Choix, S1, I1, Or1, S2, I2, PC2, Sref2,
    /// PCref2, Decroch2, Fleche, TolGuide, First, Last, Inside, Appro, Forward,
    /// RecP, RecS, RecRst, Soldep)` (`hxx:248-272`): edge / face contact on `S2`.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:59-85`.
    #[allow(clippy::too_many_arguments)]
    fn simul_surf_on_s2(
        &mut self,
        _data: &mut SurfData,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _or1: Orientation,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _pc2: Option<&Curve2dAdaptor>,
        _sref2: &SurfaceAdaptor,
        _pcref2: &Curve2dAdaptor,
        _decroch2: &mut bool,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p: bool,
        _rec_s: bool,
        _rec_rst: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("SimulSurf Not Implemented".to_string())
    }

    /// `SimulSurf(Data, Guide, Spine, Choix, S1, I1, PC1, Sref1, PCref1,
    /// Decroch1, Or1, S2, I2, PC2, Sref2, PCref2, Decroch2, Or2, Fleche,
    /// TolGuide, First, Last, Inside, Appro, Forward, RecP1, RecRst1, RecP2,
    /// RecRst2, Soldep)` (`hxx:277-302`): edge / edge contact.
    ///
    /// Base body: `ChFi3d_Builder_NotImp.cxx:89-121`.
    #[allow(clippy::too_many_arguments)]
    fn simul_surf_contact(
        &mut self,
        _data: &mut SurfData,
        _guide: &Arc<ElSpine>,
        _spine: Option<&SpineHandle>,
        _choix: i32,
        _s1: &SurfaceAdaptor,
        _i1: &Rc<RefCell<TopolTool>>,
        _pc1: Option<&Curve2dAdaptor>,
        _sref1: &SurfaceAdaptor,
        _pcref1: &Curve2dAdaptor,
        _decroch1: &mut bool,
        _or1: Orientation,
        _s2: &SurfaceAdaptor,
        _i2: &Rc<RefCell<TopolTool>>,
        _pc2: Option<&Curve2dAdaptor>,
        _sref2: &SurfaceAdaptor,
        _pcref2: &Curve2dAdaptor,
        _decroch2: &mut bool,
        _or2: Orientation,
        _fleche: f64,
        _tol_guide: f64,
        _first: &mut f64,
        _last: &mut f64,
        _inside: bool,
        _approx: bool,
        _forward: bool,
        _rec_p1: bool,
        _rec_rst1: bool,
        _rec_p2: bool,
        _rec_rst2: bool,
        _soldep: &MathVector,
    ) -> Result<(), String> {
        Err("SimulSurf Not Implemented".to_string())
    }
}
