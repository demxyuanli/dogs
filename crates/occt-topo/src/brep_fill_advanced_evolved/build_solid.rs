//! `BuildSolid` (`BRepFill_AdvancedEvolved.cxx:704-831`), `UnifyShape`
//! (`cxx:835-848`), `ExtractOuterSolid` (`cxx:852-906`) and
//! `RemoveExcessSolids` (`cxx:910-1119`).
//!
//! `BuildSolid` collects the pipe-shell faces (after `ReduceVertexTolerance`
//! and `CheckSingularityAndAdd`), fuses them together with the splits the
//! singularity pass produced, appends the lids, and lets
//! [`BopMakerVolume`](crate::bop_maker_volume::BopMakerVolume) answer the solid
//! the (not necessarily closed) face set bounds. `RemoveExcessSolids` then
//! keeps the solids that carry the lids and drop the ones MakerVolume built
//! out of flipped faces.
//!
//! Translation notes:
//! - `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Add` answers
//!   `false` for an already present shape and the hasher ignores orientation;
//!   [`shape_map_key`] models exactly that (`TopoDS_Shape.hxx:262-267`).
//! - `ShapeUpgrade_UnifySameDomain` has no port (see
//!   [`AdvancedEvolved::unify_shape`]) so `UnifyShape` is UNPORTED.
//! - `TopoDS_Shape::IsEqual` (same `TShape`, same location, **same**
//!   orientation, `TopoDS_Shape.hxx:262-267`) is
//!   [`is_equal_oriented`]; `IsSame` is [`crate::topo_tools::is_same`].

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_maker_volume::BopMakerVolume;
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, map_shapes, shape_map_key, ShapeMapKey};

use super::maps::ordered_map_shapes_and_ancestors;
use super::{AdvancedEvolved, STATUS_NOT_VOLUME, STATUS_OK};

/// `aF.IsEqual(aF1)` (`cxx:973`, `:1080`): same `TShape`, same location **and**
/// the same orientation.
fn is_equal_oriented(a: &TopoShape, b: &TopoShape) -> bool {
    is_same(a, b) && a.orientation() == b.orientation()
}

/// `BRep_Builder::MakeCompSolid` + one `Add` per solid (`cxx:1005-1013`,
/// `:1107-1118`).
fn make_comp_solid(solids: &[TopoShape]) -> TopoShape {
    let builder = TopoBuilder::new();
    let mut comp_sol = builder.make_shape(ShapeType::CompSolid);
    for s in solids {
        builder.add(&mut comp_sol, s);
    }
    comp_sol
}

/// `TopExp::MapShapes(aSh, TopAbs_FACE, aMapF)` reduced to its two lookups:
/// `FindIndex` (is the face present at all) and `FindKey` (the stored view,
/// which keeps the orientation met by the traversal).
fn find_mapped_face(
    map: &HashMap<ShapeType, Vec<TopoShape>>,
    key: &ShapeMapKey,
) -> Option<TopoShape> {
    map.get(&ShapeType::Face)?
        .iter()
        .find(|f| shape_map_key(f) == *key)
        .cloned()
}

impl AdvancedEvolved {
    /// `BuildSolid()` (`cxx:704-831`).
    pub(super) fn build_solid(&mut self) {
        // `cxx:706-709`.
        if self.status != STATUS_NOT_VOLUME {
            return;
        }

        // `cxx:711`.
        self.status = STATUS_NOT_VOLUME;

        // `cxx:713-715`: `aMapF`, `aLF`, `aLSplits`.
        let mut a_map_f: Vec<ShapeMapKey> = Vec::new();
        let mut a_lf: Vec<TopoShape> = Vec::new();
        let mut a_l_splits: Vec<TopoShape> = Vec::new();

        // `cxx:724-734`.
        let Some(pipe_shell) = self.pipe_shell.clone() else {
            return;
        };
        for face_shape in explorer_shapes(&pipe_shell, ShapeType::Face) {
            let key = shape_map_key(&face_shape);
            if a_map_f.contains(&key) {
                continue;
            }
            a_map_f.push(key);

            let a_f = crate::shape::Face(face_shape);
            self.reduce_vertex_tolerance(&a_f.0);
            self.check_singularity_and_add(&a_f.0, self.fuzzy_value, &mut a_lf, &mut a_l_splits);
        }

        // `cxx:736-766`: the split faces join the argument list and everything
        // is fused ("split interfered faces").
        a_lf.extend(a_l_splits.iter().cloned());
        if let Some(fused) = self.perform_boolean(&a_lf) {
            self.pipe_shell = Some(fused);
        }

        // `cxx:768-778`: the faces of the fused shell, deduplicated again.
        a_lf.clear();
        a_map_f.clear();
        let Some(pipe_shell) = self.pipe_shell.clone() else {
            return;
        };
        for face_shape in explorer_shapes(&pipe_shell, ShapeType::Face) {
            let key = shape_map_key(&face_shape);
            if a_map_f.contains(&key) {
                continue;
            }
            a_map_f.push(key);
            a_lf.push(face_shape);
        }

        // `cxx:780-788`: `TopoDS_Iterator(myTopBottom)` -- the direct children
        // of the lid compound.
        if let Some(top_bottom) = self.top_bottom.clone() {
            for lid in crate::iterator::ShapeIterator::of_shape(&top_bottom) {
                a_lf.push(lid);
            }
        }

        // `cxx:806-818`.
        let mut a_mv = BopMakerVolume::new();
        a_mv.set_arguments(&a_lf);
        a_mv.set_fuzzy_value(self.fuzzy_value);
        a_mv.set_intersect(true);
        a_mv.set_run_parallel(self.is_parallel);
        a_mv.set_avoid_internal_shapes(true);
        let _ = a_mv.perform();

        if a_mv.has_errors() {
            return;
        }

        // `cxx:819`.
        let a_result = a_mv.shape().clone();
        self.result = Some(a_result.clone());

        // `cxx:824-828`.
        self.remove_excess_solids(&a_l_splits, &a_result, &mut a_lf, &a_mv);
        self.unify_shape();
        self.remove_internal_wires();

        // `cxx:829`.
        self.status = STATUS_OK;
    }

    /// `UnifyShape()` (`cxx:835-848`).
    ///
    /// UNPORTED: `ShapeUpgrade_UnifySameDomain` (`ShapeUpgrade_UnifySameDomain.cxx:1-1500`,
    /// `ShapeUpgrade_UnifySameDomain.hxx:1-120`) has no port, so the
    /// `Initialize(myResult, true, true, false)` / `SetSafeInputMode(true)` /
    /// `AllowInternalEdges(false)` / `SetLinearTolerance(1e-4)` /
    /// `SetAngularTolerance(1e-2)` / `Build()` sequence (`cxx:836-846`) has no
    /// equivalent. `myResult` therefore keeps the faces MakerVolume answered,
    /// without merging the same-domain ones.
    pub(super) fn unify_shape(&mut self) {
        let _ = &self.result;
    }

    /// `ExtractOuterSolid(theShape, theArgsList)` (`cxx:852-906`).
    fn extract_outer_solid(&self, the_shape: &mut TopoShape, args: &mut Vec<TopoShape>) {
        // `cxx:854-857`: `TopExp::MapShapesAndAncestors(theShape, FACE, SOLID)`.
        let a_map_s =
            ordered_map_shapes_and_ancestors(the_shape, ShapeType::Face, ShapeType::Solid);

        // `cxx:860-868`: the faces belonging to a single solid.
        let mut a_new_list: Vec<TopoShape> = Vec::new();
        for entry in &a_map_s {
            if entry.ancestors.len() == 1 {
                a_new_list.push(entry.key.clone());
            }
        }

        // `cxx:869-872`.
        if a_new_list.is_empty() {
            return;
        }

        // `cxx:874-891`: drop from the arguments every face that is not in
        // `aNewList`.
        let mut is_removed = true;
        while is_removed {
            is_removed = false;
            let mut to_remove = None;
            for (i, f) in args.iter().enumerate() {
                if !a_new_list.iter().any(|n| is_same(n, f)) {
                    to_remove = Some(i);
                    break;
                }
            }
            if let Some(i) = to_remove {
                args.remove(i);
                is_removed = true;
            }
        }

        // `cxx:893-904`.
        let mut a_mv = BopMakerVolume::new();
        a_mv.set_arguments(&a_new_list);
        a_mv.set_intersect(true);
        a_mv.set_run_parallel(self.is_parallel);
        a_mv.set_avoid_internal_shapes(true);
        let _ = a_mv.perform();

        if a_mv.has_errors() {
            return;
        }

        *the_shape = a_mv.shape().clone();
    }

    /// `RemoveExcessSolids(theLSplits, theShape, theArgsList, theMV)`
    /// (`cxx:910-1119`).
    fn remove_excess_solids(
        &mut self,
        the_l_splits: &[TopoShape],
        the_shape: &TopoShape,
        the_args_list: &mut Vec<TopoShape>,
        the_mv: &BopMakerVolume,
    ) {
        // `cxx:916-919`.
        if self.status != STATUS_NOT_VOLUME {
            return;
        }

        let mut a_res_shape = the_shape.clone();

        // `cxx:921-947`.
        for i in 0..2 {
            let solids = explorer_shapes(&a_res_shape, ShapeType::Solid);
            // `cxx:925-931`: "No any solids".
            if solids.is_empty() {
                self.result = Some(a_res_shape);
                return;
            }
            // `cxx:933-938`: "Only one solid has been generated".
            if solids.len() == 1 {
                self.result = Some(solids[0].clone());
                return;
            }

            if i != 0 {
                break;
            }

            self.extract_outer_solid(&mut a_res_shape, the_args_list);
        }

        // `cxx:949-950`.
        let mut a_solid_list: Vec<TopoShape> = Vec::new();

        // `cxx:952-998`: "look for all solids containing lids".
        {
            // `cxx:962`: `TopExp_Explorer anExpLids(myTopBottom, TopAbs_FACE)`;
            // a null compound visits nothing.
            let lids = match &self.top_bottom {
                Some(tb) => explorer_shapes(tb, ShapeType::Face),
                None => Vec::new(),
            };
            for sol in explorer_shapes(&a_res_shape, ShapeType::Solid) {
                let a_map_f = map_shapes(&sol, &[ShapeType::Face]);

                let mut are_there_lids = false;
                // `cxx:961-983`.
                for lid in &lids {
                    are_there_lids = true;
                    let Some(stored) = find_mapped_face(&a_map_f, &shape_map_key(lid)) else {
                        continue;
                    };
                    // `cxx:973-977`.
                    if is_equal_oriented(&stored, lid) {
                        a_solid_list.push(sol.clone());
                    }
                    break;
                }

                // `cxx:981-984`.
                if !are_there_lids {
                    break;
                }
            }

            // `cxx:986-998`.
            if a_solid_list.is_empty() {
                self.result = Some(a_res_shape);
                return;
            }
            if a_solid_list.len() == 1 {
                self.result = Some(a_solid_list[0].clone());
                return;
            }
            if !a_solid_list.is_empty() {
                a_res_shape = make_comp_solid(&a_solid_list);
                a_solid_list.clear();
            }
        }

        // `cxx:1017-1094`.
        {
            // `cxx:1018-1024`: remove the split faces from the arguments.
            for f in the_l_splits {
                the_args_list.retain(|a| !is_same(a, f));
            }

            // `cxx:1026-1042`: "the face is invalid if BOPAlgo_MakerVolume
            // changes its orientation while creating solids".
            let mut a_list_inv_faces: Vec<TopoShape> = Vec::new();
            for f in the_args_list.iter() {
                for f_m in the_mv.modified(f) {
                    if f_m.orientation() != f.orientation() {
                        a_list_inv_faces.push(f_m);
                    }
                }
            }

            // `cxx:1045-1093`.
            for so in explorer_shapes(&a_res_shape, ShapeType::Solid) {
                let a_map_f = map_shapes(&so, &[ShapeType::Face]);
                let mut is_to_delete = false;

                for f in &a_list_inv_faces {
                    if find_mapped_face(&a_map_f, &shape_map_key(f)).is_some() {
                        is_to_delete = true;
                        break;
                    }
                }

                if is_to_delete {
                    continue;
                }

                for f in the_args_list.iter() {
                    // `cxx:1070-1075`: `anIdx == 0` -- the face is not in the
                    // solid.
                    let Some(f1) = find_mapped_face(&a_map_f, &shape_map_key(f)) else {
                        continue;
                    };
                    // `cxx:1080-1085`: `aF` and `aF1` are the same shape; the
                    // face is dropped when the orientations differ.
                    if !is_equal_oriented(f, &f1) {
                        is_to_delete = true;
                        break;
                    }
                }

                if is_to_delete {
                    continue;
                }

                a_solid_list.push(so);
            }
        }

        // `cxx:1097-1105`.
        if a_solid_list.is_empty() {
            self.result = Some(a_res_shape);
            return;
        }
        if a_solid_list.len() == 1 {
            self.result = Some(a_solid_list[0].clone());
            return;
        }

        // `cxx:1107-1118`.
        self.result = Some(make_comp_solid(&a_solid_list));
    }
}
