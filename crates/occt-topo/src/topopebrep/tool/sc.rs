//! `TopOpeBRepTool_SC` — the point-on-face classification helpers the 2d
//! transitions of `TopOpeBRepDS_EXPORT` read.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_SC.cxx:27-85` (OCCT 8.0.0).
//!
//! OCCT keeps one `thread_local TopOpeBRepTool_ShapeClassifier` behind
//! `FSC_GetPSC` (`cxx:27-45`) and hands it to the callers by reference. The
//! port takes the classifier as an explicit argument instead: the singleton's
//! only field that survives a call is `myRef`, and it is re-set on both entry
//! paths (`cxx:44`, `:64`), while `mySameDomain` — the one field that would
//! carry over — is never written by `FSC_*` (`TopOpeBRepTool_SC.cxx:27-69`).
//!
//! `FUN_tool_projPonS(P, BRep_Tool::Surface(F), UV, dist)` (`cxx:55-58`) is
//! [`fun_tool_proj_pon_f`], the port's `Extrema_ExtPS` projection on the same
//! surface.

use occt_core::gp::GpPnt;

use crate::abs::State;
use crate::shape::{Edge, Face, TopoShape};

use super::adaptor::CurveAdaptor;
use super::project::fun_tool_proj_pon_f;
use super::shape_classifier::ShapeClassifier;
use super::fun_tool_bounds;

/// `FSC_GetPSC()` / `FSC_GetPSC(S)` (`cxx:29-45`): a classifier whose
/// reference is `reference`.
pub fn fsc_get_psc(reference: &TopoShape) -> ShapeClassifier {
    ShapeClassifier::with_reference(reference)
}

/// `FSC_StatePonFace(P, F, PSC)` (`cxx:48-68`): the point `P` is projected on
/// the surface of `F` and classified in the face.
///
/// `false` from `FUN_tool_projPonS` (`cxx:57-61`) yields `UNKNOWN`.
pub fn fsc_state_pon_face(p: &GpPnt, f: &TopoShape, psc: &mut ShapeClassifier) -> State {
    let face = Face(f.clone());
    let Some((u, v, _dist)) = fun_tool_proj_pon_f(p, &face) else {
        return State::Unknown;
    };
    psc.set_reference(f);
    if psc
        .state_p2d_reference(&occt_core::gp::GpPnt2d::new(u, v))
        .is_err()
    {
        return State::Unknown;
    }
    psc.state()
}

/// `FSC_StateEonFace(E, t, F, PSC)` (`cxx:71-85`): the point of `E` at
/// `par = (1 - t) * f + t * l` classified in `F`.
pub fn fsc_state_eon_face(
    e: &TopoShape,
    t: f64,
    f: &TopoShape,
    psc: &mut ShapeClassifier,
) -> State {
    let (fpar, lpar) = fun_tool_bounds(e);
    let par = (1.0 - t) * fpar + t * lpar;
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let Some(p) = bc.d0(par) else {
        return State::Unknown;
    };
    fsc_state_pon_face(&p, f, psc)
}
