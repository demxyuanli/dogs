//! `HLRBRep_EdgeIList`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_EdgeIList.{hxx,cxx}`,
//! `hxx:29-44`, `cxx:28-182`): the two static list operations
//! `HLRBRep_Hider` uses on its interference lists.
//!
//! `AddInterference` (`cxx:28-45`) is an insertion sort on the interference
//! parameter: the new element goes before the first strictly greater one and is
//! appended when the list is empty or already ordered.
//!
//! `ProcessComplex` (`cxx:69-127`) collapses every run of consecutive
//! interferences that share one geometric locus into the first of them, with
//! the cumulated transitions of [`EdgeFaceTransition`] over the run.
//!
//! Deviations:
//! * `NCollection_List<HLRAlgo_Interference>` is a `Vec`; `IL.InsertBefore` /
//!   `IL.Remove(It)` become `insert(idx, ..)` / `remove(idx)`, which keeps the
//!   OCCT cursor semantics because the removals always happen strictly after
//!   the outer cursor.
//! * `EdgeGeometry` / `InterferenceBoundaryGeometry` raise `Standard_Failure`
//!   when `LocalLEGeometry2D` / `LocalFEGeometry2D` cannot evaluate; the port
//!   returns `Err`, and `HLRBRep_Hider::hide` maps it to the `catch` arm of
//!   `HLRBRep_Hider.cxx:842-850` (skip the rest of the edge).
//! * `SimilarInterference` (`cxx:50-68`) is behind `OCCT_DEBUG_SI`, not
//!   defined in the release build; not ported.

use crate::hlr_algo::HlrAlgoInterference;
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrbrep_edge_interference_tool::HlrBRepEdgeInterferenceTool;
use crate::topcnx_edge_face_transition::EdgeFaceTransition;

/// `const double TolAng = 0.0001` (`cxx:73`).
const TOL_ANG: f64 = 0.0001;

/// `HLRBRep_EdgeIList::AddInterference(IL, I, T)` (`cxx:28-45`): insert `i`
/// before the first interference of greater parameter, else append.
pub fn add_interference(
    il: &mut Vec<HlrAlgoInterference>,
    i: &HlrAlgoInterference,
    t: &HlrBRepEdgeInterferenceTool,
) {
    let p = t.parameter_of_interference(i);
    for (idx, it) in il.iter().enumerate() {
        if p < t.parameter_of_interference(it) {
            il.insert(idx, *i);
            return;
        }
    }
    il.push(*i);
}

/// `HLRBRep_EdgeIList::ProcessComplex(IL, T)` (`cxx:69-127`): replace every run
/// of same-locus interferences by its first element carrying the cumulated
/// transition and boundary transition.
pub fn process_complex(
    il: &mut Vec<HlrAlgoInterference>,
    ds: &HlrBRepData,
    t: &HlrBRepEdgeInterferenceTool,
) -> Result<(), &'static str> {
    let mut trans_tool = EdgeFaceTransition::new();
    let mut i = 0;
    while i < il.len() {
        if i + 1 < il.len() && t.same_interferences(&il[i], &il[i + 1]) {
            let first = il[i];
            let param = t.parameter_of_interference(&first);
            let (tgt_e, norm_e, curv_e) = t
                .edge_geometry(ds, param)
                .ok_or("HLRBRep_EdgeIList::ProcessComplex: EdgeGeometry failed")?;
            trans_tool.reset(&tgt_e, &norm_e, curv_e);
            let (tgt_i, norm_i, curv_i) = t
                .interference_boundary_geometry(ds, &first)
                .ok_or("HLRBRep_EdgeIList::ProcessComplex: InterferenceBoundaryGeometry failed")?;
            trans_tool.add_interference(
                TOL_ANG,
                &tgt_i,
                &norm_i,
                curv_i,
                first.orientation(),
                first.transition(),
                first.boundary_transition(),
            );

            while i + 1 < il.len() && t.same_interferences(&il[i], &il[i + 1]) {
                let next = il[i + 1];
                let (tgt_i, norm_i, curv_i) =
                    t.interference_boundary_geometry(ds, &next).ok_or(
                        "HLRBRep_EdgeIList::ProcessComplex: InterferenceBoundaryGeometry failed",
                    )?;
                trans_tool.add_interference(
                    TOL_ANG,
                    &tgt_i,
                    &norm_i,
                    curv_i,
                    next.orientation(),
                    next.transition(),
                    next.boundary_transition(),
                );
                il.remove(i + 1);
            }

            // `It1.ChangeValue().Transition(transTool.Transition())` and the
            // boundary transition (`cxx:123-124`).
            il[i].set_transition(trans_tool.transition());
            il[i].set_boundary_transition(trans_tool.boundary_transition());
        }
        i += 1;
    }
    Ok(())
}
