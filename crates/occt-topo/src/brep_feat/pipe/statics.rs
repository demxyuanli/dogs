//! File-local helpers of `BRepFeat_MakePipe.cxx`.

use crate::abs::ShapeType;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::Pipe;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:354-397`).
///
/// Records the first/last wire of the pipe and, for every profile edge not yet
/// keyed, the shapes swept from it.
pub(super) fn maj_map(
    the_b: &TopoShape,
    the_p: &mut Pipe,
    the_map: &mut ShapeValueMap<Vec<TopoShape>>,
    the_f_shape: &mut Option<TopoShape>,
    the_l_shape: &mut Option<TopoShape>,
) -> Result<(), String> {
    // `cxx:363-374`.
    if let Some(first) = the_p.first_shape().cloned() {
        if let Some(w) = explorer_shapes(&first, ShapeType::Wire)
            .into_iter()
            .next()
        {
            *the_f_shape = Some(w.clone());
            the_map.bind(&w, Vec::new());
            let faces = explorer_shapes(&first, ShapeType::Face);
            if let Some(slot) = the_map.change(&w) {
                slot.extend(faces);
            }
        }
    }

    // `cxx:376-387`.
    if let Some(last) = the_p.last_shape().cloned() {
        if let Some(w) = explorer_shapes(&last, ShapeType::Wire)
            .into_iter()
            .next()
        {
            *the_l_shape = Some(w.clone());
            the_map.bind(&w, Vec::new());
            let faces = explorer_shapes(&last, ShapeType::Face);
            if let Some(slot) = the_map.change(&w) {
                slot.extend(faces);
            }
        }
    }

    // `cxx:389-396`. `LocOpe_Pipe::Shapes` clears its own result buffer, so the
    // answer is copied out before the map slot is written.
    for e in explorer_shapes(the_b, ShapeType::Edge) {
        if !the_map.is_bound(&e) {
            the_map.bind(&e, Vec::new());
            let shapes = the_p.shapes(&e)?.to_vec();
            if let Some(slot) = the_map.change(&e) {
                *slot = shapes;
            }
        }
    }
    Ok(())
}

/// The `Shape()` half of a `BRepAlgoAPI_*` operation, which raises
/// `StdFail_NotDone` when the run failed.
///
/// `BRepFeat_MakePipe` keeps the operation alive because
/// `UpdateDescendants` queries its history (`cxx:201-211`).
pub(super) fn boolean_result(
    mut bop: BooleanOperation,
) -> Result<(BooleanOperation, TopoShape), String> {
    bop.build()?;
    let shape = bop
        .shape()
        .cloned()
        .ok_or_else(|| "StdFail_NotDone".to_string())?;
    Ok((bop, shape))
}
