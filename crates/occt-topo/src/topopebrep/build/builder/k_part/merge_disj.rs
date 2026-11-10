//! `TopOpeBRepBuild_Builder` — the `isdisj` KPart case.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! (OCCT 8.0.0):
//! - `MergeKPartisdisj` (`:840-1030`),
//! - `LocalKPisdisjanalyse` (`:2388-2484`) and its
//!   `TopOpeBRepBuild_KPart_Operation` (`:66-74`),
//! - `BuildNewSolid` (`:2486-2655`),
//! - `disjPerformFuse` (`:2660-2775`), `disjPerformCommon` (`:2781-2867`),
//!   `disjPerformCut` (`:2873-2978`).
//!
//! `LocalKPisdisjanalyse` duplicates `KPisdisjanalyse` (`:2103-2194`) with the
//! operation passed in; the latter is the `Ope*()`-dispatched form and lives in
//! `predicates.rs`'s sibling `lists.rs` only through `KPisdisjanalyse`, which
//! `MergeKPartisdisj` itself does not call.

use std::collections::HashSet;

use crate::abs::{ShapeType, State};
use crate::brep_class3d::outer_shell;
use crate::shape::{Solid, TopoShape};
use crate::topo_tools::{explorer_shapes, is_equal, shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::build_tool::BuildTool;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;

use super::super::k_part::{
    KP_ISDISJ, RESNEWSHA1, RESNEWSHA2, RESNULL, RESSHAPE1, RESSHAPE12, RESSHAPE2, RESUNDEF,
    SHEAUCU, SHEGARDCOUR, SHEGARDAUTR, SHEGARDTOUS, SHECLASAUTR, SHEUNDEF,
};
use super::Builder;

/// `enum TopOpeBRepBuild_KPart_Operation` (`KPart.cxx:66-74`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Operation {
    Fuse,
    Common,
    Cut12,
    Cut21,
}

/// `LocalKPisdisjanalyse(Stsol1, Stsol2, theOperation, ires, icla1, icla2)`
/// (`KPart.cxx:2388-2484`).
fn local_kp_isdisj_analyse(
    stsol1: State,
    stsol2: State,
    op: Operation,
) -> (i32, i32, i32) {
    let mut ires = RESUNDEF;
    let mut icla1 = SHEUNDEF;
    let mut icla2 = SHEUNDEF;

    match op {
        Operation::Fuse => {
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE12;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNEWSHA1;
                icla1 = SHECLASAUTR;
                icla2 = SHECLASAUTR;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNEWSHA2;
                icla1 = SHECLASAUTR;
                icla2 = SHECLASAUTR;
            }
        }
        Operation::Cut12 => {
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE1;
                icla1 = SHEGARDTOUS;
                icla2 = SHEAUCU;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNEWSHA1;
                icla1 = SHECLASAUTR;
                icla2 = SHEGARDCOUR;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNULL;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            }
        }
        Operation::Cut21 => {
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE2;
                icla1 = SHEAUCU;
                icla2 = SHEGARDTOUS;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNULL;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNEWSHA2;
                icla1 = SHEGARDCOUR;
                icla2 = SHECLASAUTR;
            }
        }
        Operation::Common => {
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESNULL;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNEWSHA2;
                icla1 = SHECLASAUTR;
                icla2 = SHEGARDAUTR;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNEWSHA1;
                icla1 = SHEGARDAUTR;
                icla2 = SHECLASAUTR;
            }
        }
    }
    (ires, icla1, icla2)
}

/// The `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` of
/// `BuildNewSolid` (`KPart.cxx:2502`): `IsSame` keys, orientation ignored.
fn map_has(set: &HashSet<ShapeMapKey>, s: &TopoShape) -> bool {
    set.contains(&shape_map_key(s))
}

/// `ChangeMerged`'s `NCollection_List` append of one shape per solid, in order,
/// without duplicates removed — OCCT's list keeps duplicates.
fn add_to_indexed_map(map: &mut Vec<TopoShape>, s: &TopoShape) {
    if map.iter().any(|m| shape_map_key(m) == shape_map_key(s)) {
        return;
    }
    map.push(s.clone());
}

/// `TopoDS_Shape::Complement()`.
fn complemented(s: &TopoShape) -> TopoShape {
    let mut r = s.clone();
    r.set_orientation(s.orientation().complement());
    r
}

impl Builder {
    /// `MergeKPartisdisj()` (`KPart.cxx:840-1030`).
    ///
    /// OCCT holds `lmergesha1 = ChangeMerged(myShape1, myState1)` for the whole
    /// body; the port fills the same list at the end so the borrow of
    /// `myMerged1` does not outlive one statement.
    pub fn merge_k_part_isdisj(&mut self) {
        if self.my_is_k_part != KP_ISDISJ {
            return;
        }
        let (Some(shape1), Some(shape2)) = (self.my_shape1.clone(), self.my_shape2.clone()) else {
            // OCCT dereferences both operands; the port's null is `None`.
            return;
        };
        let st1 = self.my_state1;
        let st2 = self.my_state2;

        // cxx:853-854
        self.change_merged(&shape1, st1);
        self.change_merged(&shape2, st2);

        let merged = self.isdisj_merged(&shape1, &shape2);
        self.change_merged(&shape1, st1).extend(merged);
    }

    /// The body of `MergeKPartisdisj` (`cxx:856-1028`); the returned list is
    /// what the body appended to `lmergesha1`.
    fn isdisj_merged(&mut self, shape1: &TopoShape, shape2: &TopoShape) -> Vec<TopoShape> {
        let mut lmergesha1: Vec<TopoShape> = Vec::new();

        // cxx:857-862
        let exsol1 = explorer_shapes(shape1, ShapeType::Solid);
        let exsol2 = explorer_shapes(shape2, ShapeType::Solid);
        let hassol1 = !exsol1.is_empty();
        let hassol2 = !exsol2.is_empty();

        // cxx:864-880: a free shell / face / wire / edge / vertex on either
        // operand sends the case to the general path.
        let mut hasnotsol1 = false;
        let mut hasnotsol2 = false;
        for t in [
            ShapeType::Shell,
            ShapeType::Face,
            ShapeType::Wire,
            ShapeType::Edge,
            ShapeType::Vertex,
        ] {
            if hasnotsol1 || hasnotsol2 {
                break;
            }
            if !explorer_shapes(shape1, t).is_empty() {
                hasnotsol1 = true;
            }
            if !explorer_shapes(shape2, t).is_empty() {
                hasnotsol2 = true;
            }
        }
        let soldisj = (hassol1 && hassol2) && !(hasnotsol1 || hasnotsol2);

        // cxx:882-886
        let sol1;
        let outsha1;
        let sol2;
        let outsha2;

        if soldisj {
            // cxx:891-895
            let map_of_solid1 = explorer_shapes(shape1, ShapeType::Solid);
            let map_of_solid2 = explorer_shapes(shape2, ShapeType::Solid);

            // cxx:896-1001
            if map_of_solid1.len() > 1 || map_of_solid2.len() > 1 {
                let result = if self.opefus() {
                    disj_perform_fuse(&map_of_solid1, &map_of_solid2)
                } else if self.opec12() {
                    disj_perform_cut(&map_of_solid1, &map_of_solid2)
                } else if self.opec21() {
                    disj_perform_cut(&map_of_solid2, &map_of_solid1)
                } else if self.opecom() {
                    disj_perform_common(&map_of_solid1, &map_of_solid2)
                } else {
                    None
                };
                let Some(result) = result else {
                    return lmergesha1;
                };
                // cxx:926-929
                lmergesha1.extend(result);
                return lmergesha1;
            }

            // cxx:933-1000
            sol1 = exsol1[0].clone();
            self.change_merged(&sol1, self.my_state1);
            outsha1 = outer_shell(&sol1).unwrap_or_else(|| TopoShape::new(ShapeType::Shell));

            sol2 = exsol2[0].clone();
            self.change_merged(&sol2, self.my_state2);
            outsha2 = outer_shell(&sol2).unwrap_or_else(|| TopoShape::new(ShapeType::Shell));

            let stsol1 = self.kp_clas_ss(&outsha1, &sol2);
            let stsol2 = self.kp_clas_ss(&outsha2, &sol1);

            let (ires, icla1, icla2) = self.kp_isdisj_analyse(stsol1, stsol2);

            if ires == RESUNDEF
                || icla1 == SHEUNDEF
                || icla2 == SHEUNDEF
                || ires == RESNULL
            {
                return lmergesha1;
            } else if ires == RESSHAPE12 {
                lmergesha1.push(shape1.clone());
                lmergesha1.push(shape2.clone());
                return lmergesha1;
            } else if ires == RESSHAPE1 {
                lmergesha1.push(shape1.clone());
                return lmergesha1;
            } else if ires == RESSHAPE2 {
                lmergesha1.push(shape2.clone());
                return lmergesha1;
            } else if ires == RESNEWSHA1 || ires == RESNEWSHA2 {
                // cxx:984-990
                if let Some(newsol) = build_new_solid(
                    &sol1,
                    &sol2,
                    stsol1,
                    stsol2,
                    ires,
                    icla1,
                    icla2,
                    self.my_state1,
                    self.my_state2,
                ) {
                    lmergesha1.push(newsol);
                }
                return lmergesha1;
            } else {
                return lmergesha1;
            }
        }

        // cxx:1004-1028: not a disjoint-solid case.
        if self.opec12() {
            lmergesha1.push(shape1.clone());
        } else if self.opec21() {
            lmergesha1.push(shape2.clone());
        } else if self.opecom() {
            lmergesha1.clear();
        } else if self.opefus() {
            lmergesha1.push(shape1.clone());
            lmergesha1.push(shape2.clone());
        } else {
            return lmergesha1;
        }

        lmergesha1
    }

    /// `KPisdisjanalyse(Stsol1, Stsol2, ires, icla1, icla2)` (`KPart.cxx:2103-2194`).
    ///
    /// The `Ope*()` dispatch of `LocalKPisdisjanalyse`; `MergeKPartisdisj` uses
    /// the latter, `KPiskoleanalyse` the former.
    pub fn kp_isdisj_analyse(&self, stsol1: State, stsol2: State) -> (i32, i32, i32) {
        let op = if self.opefus() {
            Operation::Fuse
        } else if self.opec12() {
            Operation::Cut12
        } else if self.opec21() {
            Operation::Cut21
        } else if self.opecom() {
            Operation::Common
        } else {
            // cxx:2107-2110: none of the four matches leaves `RESUNDEF`.
            return (RESUNDEF, SHEUNDEF, SHEUNDEF);
        };
        local_kp_isdisj_analyse(stsol1, stsol2, op)
    }
}

/// `BuildNewSolid(sol1, sol2, stsol1, stsol2, ires, icla1, icla2, theState1,
/// theState2)` (`KPart.cxx:2486-2655`).
#[allow(clippy::too_many_arguments)]
fn build_new_solid(
    sol1: &TopoShape,
    sol2: &TopoShape,
    stsol1: State,
    stsol2: State,
    ires: i32,
    icla1: i32,
    icla2: i32,
    the_state1: State,
    the_state2: State,
) -> Option<TopoShape> {
    let a_build_tool = BuildTool::new();
    let outsha1 = outer_shell(sol1);
    let outsha2 = outer_shell(sol2);
    // `Snull` for `StateShapeShape(shacur, Snull, othersol)`.
    let _null_shape = TopoShape::new(ShapeType::Shape);

    let mut isdisjmap: HashSet<ShapeMapKey> = HashSet::new();
    let mut newsol = TopoShape::new(ShapeType::Solid);
    a_build_tool.make_solid(&mut newsol);

    // cxx:2514-2529
    if ires == RESNEWSHA1 {
        if let Some(outsha1) = &outsha1 {
            if !map_has(&isdisjmap, outsha1) {
                isdisjmap.insert(shape_map_key(outsha1));
                a_build_tool.add_solid_shell(&mut newsol, outsha1);
            }
        }
    } else if ires == RESNEWSHA2 {
        if let Some(outsha2) = &outsha2 {
            if !map_has(&isdisjmap, outsha2) {
                isdisjmap.insert(shape_map_key(outsha2));
                a_build_tool.add_solid_shell(&mut newsol, outsha2);
            }
        }
    }

    let mut classifier = ShapeClassifier::new();

    // cxx:2531-2594 (block1: sol1)
    add_new_solid_shells(
        &a_build_tool,
        &mut newsol,
        &mut isdisjmap,
        &mut classifier,
        sol1,
        stsol1,
        the_state1,
        outsha1.as_ref(),
        icla1,
        sol2,
        true,
    );
    // cxx:2596-2651 (block2: sol2)
    add_new_solid_shells(
        &a_build_tool,
        &mut newsol,
        &mut isdisjmap,
        &mut classifier,
        sol2,
        stsol2,
        the_state2,
        outsha2.as_ref(),
        icla2,
        sol1,
        false,
    );

    // cxx:2655
    Some(newsol)
}

/// One block of `BuildNewSolid` (`KPart.cxx:2537-2594` for `sol1`,
/// `:2601-2651` for `sol2`).
#[allow(clippy::too_many_arguments)]
fn add_new_solid_shells(
    a_build_tool: &BuildTool,
    newsol: &mut TopoShape,
    isdisjmap: &mut HashSet<ShapeMapKey>,
    classifier: &mut ShapeClassifier,
    sol: &TopoShape,
    solstate: State,
    shastatetoadd: State,
    outsha: Option<&TopoShape>,
    icla: i32,
    othersol: &TopoShape,
    guard_with_map: bool,
) {
    let empty_shell = TopoShape::new(ShapeType::Shell);
    let outsha = outsha.cloned().unwrap_or(empty_shell);

    for shacur in explorer_shapes(sol, ShapeType::Shell) {
        let isoutsha = is_equal(&shacur, &outsha);

        // cxx:2547-2557
        let mut garde = true;
        if icla == SHEAUCU {
            garde = false;
        } else if icla == SHEGARDAUTR || icla == SHECLASAUTR {
            garde = !isoutsha;
        }
        if !garde {
            continue;
        }

        // cxx:2559-2574
        let mut add = false;
        if icla == SHEGARDCOUR || icla == SHEGARDAUTR || icla == SHEGARDTOUS {
            add = true;
        } else if icla == SHECLASAUTR {
            let state = classifier
                .state_shape_shape_null_avoid(&shacur, othersol)
                .unwrap_or(State::Unknown);
            add = state == shastatetoadd;
        }

        if add {
            let mut shaori = shacur.clone();
            // cxx:2580-2584
            if solstate == State::In {
                shaori = complemented(&shaori);
            }
            // cxx:2585-2589 / cxx:2643-2648: block2 has no map guard.
            if guard_with_map {
                if !map_has(isdisjmap, &shaori) {
                    isdisjmap.insert(shape_map_key(&shaori));
                    a_build_tool.add_solid_shell(newsol, &shaori);
                }
            } else {
                a_build_tool.add_solid_shell(newsol, &shaori);
            }
        }
    }
}

/// `disjPerformFuse(theMapOfSolid1, theMapOfSolid2, theMapOfResult)`
/// (`KPart.cxx:2660-2775`).
fn disj_perform_fuse(
    the_map_of_solid1: &[TopoShape],
    the_map_of_solid2: &[TopoShape],
) -> Option<Vec<TopoShape>> {
    let mut a_map_of_solid: Vec<TopoShape> = Vec::new();
    for s in the_map_of_solid1.iter().chain(the_map_of_solid2.iter()) {
        add_to_indexed_map(&mut a_map_of_solid, s);
    }

    let a_max_number_of_iterations = a_map_of_solid.len() * a_map_of_solid.len();
    let mut a_map_of_used_solids: HashSet<ShapeMapKey> = HashSet::new();
    let mut classifier = ShapeClassifier::new();
    let _null_shape = TopoShape::new(ShapeType::Shape);
    let mut the_map_of_result: Vec<TopoShape> = Vec::new();

    for i in 0..a_map_of_solid.len() {
        let localshape1 = a_map_of_solid[i].clone();
        // cxx:2689-2692
        if localshape1.shape_type() != ShapeType::Solid {
            return None;
        }
        let sol1 = Solid::wrap(localshape1.clone())?;
        let mut acurrentsolid = sol1.0;
        if map_has(&a_map_of_used_solids, &localshape1) {
            continue;
        }

        let mut j = 0usize;
        let mut acheckiterator = 0usize;
        while j < a_map_of_solid.len() && acheckiterator <= a_max_number_of_iterations {
            acheckiterator += 1;
            if j == i {
                j += 1;
                continue;
            }
            let localshape2 = a_map_of_solid[j].clone();
            // cxx:2711-2714
            if localshape2.shape_type() != ShapeType::Solid {
                return None;
            }
            j += 1;

            if map_has(&a_map_of_used_solids, &localshape2) {
                continue;
            }
            let sol2 = Solid::wrap(localshape2.clone())?.0;
            let outsha2 = outer_shell(&sol2)?;
            let outsha1 = outer_shell(&acurrentsolid)?;

            let stsol1 = classifier
                .state_shape_shape_null_avoid(&outsha1, &sol2)
                .unwrap_or(State::Unknown);
            let stsol2 = classifier
                .state_shape_shape_null_avoid(&outsha2, &acurrentsolid)
                .unwrap_or(State::Unknown);

            let (ires, icla1, icla2) = local_kp_isdisj_analyse(stsol1, stsol2, Operation::Fuse);
            // cxx:2733-2738
            if ires == RESUNDEF || icla1 == SHEUNDEF || icla2 == SHEUNDEF || ires == RESNULL {
                eprintln!("Warning: disjPerformFuse: can not determine solid's states");
                continue;
            }
            // cxx:2740-2742
            if ires == RESSHAPE12 {
                continue;
            }

            // cxx:2744-2767
            if ires == RESNEWSHA1 || ires == RESNEWSHA2 {
                let Some(newsol) = build_new_solid(
                    &acurrentsolid,
                    &sol2,
                    stsol1,
                    stsol2,
                    ires,
                    icla1,
                    icla2,
                    State::Out,
                    State::Out,
                ) else {
                    return None;
                };
                j = 0;
                acurrentsolid = newsol;
                a_map_of_used_solids.insert(shape_map_key(&localshape2));
            }
        }
        // cxx:2768-2772
        if acheckiterator > a_max_number_of_iterations {
            eprintln!("disjPerformFuse: programming error");
            return None;
        }
        the_map_of_result.push(acurrentsolid);
    }

    Some(the_map_of_result)
}

/// `disjPerformCommon(theMapOfSolid1, theMapOfSolid2, theMapOfResult)`
/// (`KPart.cxx:2781-2867`).
fn disj_perform_common(
    the_map_of_solid1: &[TopoShape],
    the_map_of_solid2: &[TopoShape],
) -> Option<Vec<TopoShape>> {
    let mut classifier = ShapeClassifier::new();

    // cxx:2797-2798
    let a_map_of_separated_solid1 = disj_perform_fuse(the_map_of_solid1, the_map_of_solid1)?;
    let a_map_of_separated_solid2 = disj_perform_fuse(the_map_of_solid2, the_map_of_solid2)?;

    let mut a_map_of_common_of_couple: Vec<TopoShape> = Vec::new();

    // cxx:2800-2865
    for localshape1 in &a_map_of_separated_solid1 {
        if localshape1.shape_type() != ShapeType::Solid {
            return None;
        }
        let sol1 = Solid::wrap(localshape1.clone())?.0;
        let outsha1 = outer_shell(&sol1)?;

        for localshape2 in &a_map_of_separated_solid2 {
            if localshape2.shape_type() != ShapeType::Solid {
                return None;
            }
            let sol2 = Solid::wrap(localshape2.clone())?.0;
            let outsha2 = outer_shell(&sol2)?;

            let stsol1 = classifier
                .state_shape_shape_null_avoid(&outsha1, &sol2)
                .unwrap_or(State::Unknown);
            let stsol2 = classifier
                .state_shape_shape_null_avoid(&outsha2, &sol1)
                .unwrap_or(State::Unknown);

            let (ires, icla1, icla2) =
                local_kp_isdisj_analyse(stsol1, stsol2, Operation::Common);
            // cxx:2830-2834
            if ires == RESUNDEF || icla1 == SHEUNDEF || icla2 == SHEUNDEF {
                eprintln!("Warning: disjPerformCommon: can not determine solid's states");
                continue;
            }
            // cxx:2835-2862
            match ires {
                RESNULL => continue,
                RESSHAPE12 => {
                    add_to_indexed_map(&mut a_map_of_common_of_couple, &sol1);
                    add_to_indexed_map(&mut a_map_of_common_of_couple, &sol2);
                    continue;
                }
                RESSHAPE1 => {
                    add_to_indexed_map(&mut a_map_of_common_of_couple, &sol1);
                    continue;
                }
                RESSHAPE2 => {
                    add_to_indexed_map(&mut a_map_of_common_of_couple, &sol2);
                }
                RESNEWSHA1 | RESNEWSHA2 => {
                    let Some(newsol) = build_new_solid(
                        &sol1, &sol2, stsol1, stsol2, ires, icla1, icla2, State::In, State::In,
                    ) else {
                        return None;
                    };
                    add_to_indexed_map(&mut a_map_of_common_of_couple, &newsol);
                }
                _ => continue,
            }
        }
    }

    // cxx:2866
    disj_perform_fuse(&a_map_of_common_of_couple, &a_map_of_common_of_couple)
}

/// `disjPerformCut(theMapOfSolid1, theMapOfSolid2, theMapOfResult)`
/// (`KPart.cxx:2873-2978`).
fn disj_perform_cut(
    the_map_of_solid1: &[TopoShape],
    the_map_of_solid2: &[TopoShape],
) -> Option<Vec<TopoShape>> {
    let mut classifier = ShapeClassifier::new();
    let mut the_map_of_result: Vec<TopoShape> = Vec::new();

    // cxx:2892-2893
    let a_map_of_separated_solid1 = disj_perform_fuse(the_map_of_solid1, the_map_of_solid1)?;
    let a_map_of_separated_solid2 = disj_perform_fuse(the_map_of_solid2, the_map_of_solid2)?;

    // cxx:2895-2976
    for localshape1 in &a_map_of_separated_solid1 {
        if localshape1.shape_type() != ShapeType::Solid {
            return None;
        }
        let sol1 = Solid::wrap(localshape1.clone())?.0;
        let mut acurrentsolid = sol1;
        let mut null_result = false;

        for localshape2 in &a_map_of_separated_solid2 {
            if null_result {
                break;
            }
            if localshape2.shape_type() != ShapeType::Solid {
                return None;
            }
            let sol2 = Solid::wrap(localshape2.clone())?.0;
            let outsha2 = outer_shell(&sol2)?;
            let outsha1 = outer_shell(&acurrentsolid)?;

            let stsol1 = classifier
                .state_shape_shape_null_avoid(&outsha1, &sol2)
                .unwrap_or(State::Unknown);
            let stsol2 = classifier
                .state_shape_shape_null_avoid(&outsha2, &acurrentsolid)
                .unwrap_or(State::Unknown);

            let (ires, icla1, icla2) =
                local_kp_isdisj_analyse(stsol1, stsol2, Operation::Cut12);
            // cxx:2927-2932
            if ires == RESUNDEF || icla1 == SHEUNDEF || icla2 == SHEUNDEF {
                eprintln!("Warning: disjPerformCut: can not determine solid's states");
                continue;
            }
            // cxx:2933-2966
            match ires {
                RESNULL => null_result = true,
                RESSHAPE12 => null_result = true,
                RESSHAPE1 => null_result = false,
                RESSHAPE2 => null_result = true,
                RESNEWSHA1 | RESNEWSHA2 => {
                    let Some(newsol) = build_new_solid(
                        &acurrentsolid,
                        &sol2,
                        stsol1,
                        stsol2,
                        ires,
                        icla1,
                        icla2,
                        State::Out,
                        State::In,
                    ) else {
                        return None;
                    };
                    acurrentsolid = newsol;
                }
                _ => continue,
            }
        }

        // cxx:2968-2975
        if !null_result {
            the_map_of_result.push(acurrentsolid);
        }
    }
    Some(the_map_of_result)
}
