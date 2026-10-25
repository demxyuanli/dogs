//! `TopOpeBRepTool::RegularizeShells` (`TopOpeBRepTool_RegularizeS.cxx:39-74`).
//!
//! Drives one `TopOpeBRepTool_REGUS` over every shell of a solid: `Init`,
//! `MapS`, `SplitFaces`, then `REGU`. The instance is created once so
//! `myOshNsh` and `myFsplits` accumulate over the shells; the returned maps
//! answer which old shell produced which new shells, and which ancestor face
//! was split into what.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, ShapeMapKey};
use crate::topopebrep::tool::regus::Regus;

/// `TopOpeBRepTool::RegularizeShells(theSolid, OldSheNewShe, FSplits)`
/// (`cxx:39-74`).
///
/// `cxx:70-72`: `REGU`'s answer is ignored on purpose; the two maps are the
/// result. `false` only comes out of `MapS` / `SplitFaces` (`cxx:60-68`).
pub fn regularize_shells(
    the_solid: &TopoShape,
    old_she_new_she: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    f_splits: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    old_she_new_she.clear();
    f_splits.clear();

    let mut regus = Regus::new();
    regus.set_osh_nsh(old_she_new_she);
    regus.set_f_splits(f_splits);

    for sh in explorer_shapes(the_solid, ShapeType::Shell) {
        regus.init(&sh);
        if !regus.map_s() {
            return false;
        }
        if !regus.split_faces() {
            return false;
        }
        regus.regu();
    }

    *old_she_new_she = regus.osh_nsh().clone();
    *f_splits = regus.f_splits().clone();
    true
}
