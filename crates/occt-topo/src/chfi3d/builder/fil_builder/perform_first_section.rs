//! `ChFi3d_FilBuilder::PerformFirstSection`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:1500-1537`
//! (OCCT 8.0.0). The initial section of the fillet band along one guideline:
//! the blend function of the section shape, then
//! `BRepBlend_Walking::PerformFirstSection` on it.
//!
//! The body is split in two by `fsp->IsConstant()` (`cxx:1519`): the constant
//! radius arm runs `BRepBlend_ConstRad` / `Func.Set(fsp->Radius(), Choix)` /
//! `Func.Set(myShape)` (`cxx:1520-1522`), the variable one
//! `BRepBlend_EvolRad` / `Func.Set(Choix)` / `Func.Set(myShape)` with
//! `fsp->Law(HGuide)` (`cxx:1528-1530`). Both then build
//! `BRepBlend_Walking(S1, S2, I1, I2, HGuide)` and call the same
//! `PerformFirstSection(Func, Par, SolDep, tolapp3d, TolGuide, Pos1, Pos2)`.
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("PerformSurf : this is not the spine of
//!   a fillet")` (`cxx:1514-1517`) is the `Err` of [`perform_first_section`];
//!   the trait entry point answers `false`,
//! - `TopAbs_State& Pos1 / Pos2` is [`State`] at the trait boundary (the
//!   `ChFi3d_Builder::StartSol` signature), while
//!   `BRepBlend_Walking::PerformFirstSection`'s `Pos1` / `Pos2` are the
//!   `TopAbs_State`-shaped [`FaceState`] of its own port; [`state_face`] and
//!   [`state_of`] are the one-to-one rename between the two four-valued
//!   enumerations,
//! - `HGuide->Resolution(tolapp3d)` (`cxx:1518`) raises in OCCT when the
//!   guideline carries no curve; that is the `Err` here.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::law::{LawFunction, LawHandle};
use occt_math::MathVector;

use crate::abs::State;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{ConstRad, EvolRad};
use crate::brep_blend::walking::Walking;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::ChFi3dCtx;
use crate::chfi_ds::{ElSpine, SpineClass, SpineHandle};
use crate::fclass2d::FaceState;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::simul_surf::guide_as_curve;

/// `ChFi3d_FilBuilder::PerformFirstSection(Spine, HGuide, Choix, S1, S2, I1, I2,
/// Par, SolDep, Pos1, Pos2)` (`cxx:1500-1537`).
///
/// `shape` is the subclass member `myShape` (`hxx:377`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_first_section(
    shape: BlendFuncSectionShape,
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
    ctx: &ChFi3dCtx<'_>,
) -> Result<bool, String> {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of a
    // fillet");` (`cxx:1512-1517`).
    let is_fil = {
        let guard = spine
            .read()
            .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?;
        guard.class() == SpineClass::Fil
    };
    if !is_fil {
        return Err("PerformSurf : this is not the spine of a fillet".to_string());
    }

    // `double TolGuide = HGuide->Resolution(tolapp3d);` (`cxx:1518`).
    let tol_guide = guide
        .resolution(ctx.tolapp3d)
        .ok_or_else(|| "PerformFirstSection: HGuide has no curve".to_string())?;

    // The two adaptor surfaces, shared by the blend function and the walker;
    // OCCT passes the `BRepAdaptor_Surface` handles themselves (`cxx:1503-1504`).
    let s1_surf = s1
        .surface()
        .ok_or_else(|| "PerformFirstSection: S1 has no surface".to_string())?
        .clone();
    let s2_surf = s2
        .surface()
        .ok_or_else(|| "PerformFirstSection: S2 has no surface".to_string())?
        .clone();
    let guide_curve: Arc<dyn Curve> = guide_as_curve(guide);

    // `if (fsp->IsConstant())` (`cxx:1519`).
    let constant = spine
        .read()
        .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?
        .is_fil_constant();
    if constant {
        // `BRepBlend_ConstRad Func(S1, S2, HGuide); Func.Set(fsp->Radius(),
        // Choix); Func.Set(myShape);` (`cxx:1520-1522`).
        let radius = spine
            .read()
            .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?
            .fil_radius()
            .ok_or_else(|| {
                "PerformFirstSection: ChFiDS_FilSpine::Radius raised Standard_DomainError"
                    .to_string()
            })?;
        let mut func = ConstRad::new(Arc::clone(&s1_surf), Arc::clone(&s2_surf), guide_curve);
        func.set_radius(radius, choix);
        func.set_section_shape(shape);

        // `BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide);` (`cxx:1523`).
        let mut walk = Walking::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Rc::clone(i1),
            Rc::clone(i2),
            Arc::clone(guide),
        );
        // `return TheWalk.PerformFirstSection(Func, Par, SolDep, tolapp3d,
        // TolGuide, Pos1, Pos2);` (`cxx:1524-1525`).
        let (mut f1, mut f2) = (state_face(*pos1), state_face(*pos2));
        let done = walk.perform_first_section(
            &mut func,
            par,
            sol_dep,
            ctx.tolapp3d,
            tol_guide,
            &mut f1,
            &mut f2,
        );
        *pos1 = state_of(f1);
        *pos2 = state_of(f2);
        Ok(done)
    } else {
        // `BRepBlend_EvolRad Func(S1, S2, HGuide, fsp->Law(HGuide));
        // Func.Set(Choix); Func.Set(myShape);` (`cxx:1528-1530`).
        let law = spine
            .read()
            .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?
            .fil_law(guide)
            .ok_or_else(|| {
                "PerformFirstSection: ChFiDS_FilSpine::Law raised Standard_DomainError".to_string()
            })?;
        let law_f: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
        let mut func = EvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            guide_curve,
            law_f,
        );
        func.set_choix(choix);
        func.set_section_shape(shape);

        // `BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide);` (`cxx:1531`).
        let mut walk = Walking::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Rc::clone(i1),
            Rc::clone(i2),
            Arc::clone(guide),
        );
        // `return TheWalk.PerformFirstSection(...)` (`cxx:1532-1534`).
        let (mut f1, mut f2) = (state_face(*pos1), state_face(*pos2));
        let done = walk.perform_first_section(
            &mut func,
            par,
            sol_dep,
            ctx.tolapp3d,
            tol_guide,
            &mut f1,
            &mut f2,
        );
        *pos1 = state_of(f1);
        *pos2 = state_of(f2);
        Ok(done)
    }
}

/// `TopAbs_State` as the four-valued [`FaceState`] of the classification.
///
/// The two enumerations have the same four members
/// (`TopAbs_IN` / `TopAbs_OUT` / `TopAbs_ON` / `TopAbs_UNKNOWN`), so this is a
/// rename, not a mapping between different sets.
pub(crate) fn state_face(s: State) -> FaceState {
    match s {
        State::In => FaceState::In,
        State::Out => FaceState::Out,
        State::On => FaceState::On,
        State::Unknown => FaceState::Unknown,
    }
}

/// The inverse of [`state_face`].
pub(crate) fn state_of(s: FaceState) -> State {
    match s {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}
