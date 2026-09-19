//! 3-D region construction — a Rust port of `BOPAlgo_BuilderSolid`.
//!
//! The algorithm that builds closed volumes ("areas") from a set of faces and,
//! given a source solid, splits it into its regions — the `BOPAlgo_SplitSolid`
//! that `BOPAlgo_Builder::BuildSplitSolids` runs per interfered solid
//! (`BOPAlgo_Builder_3.cxx`). The pipeline is the four-phase sequence of
//! `BOPAlgo_BuilderArea`:
//!
//! ```text
//!   ShapesToAvoid → Loops → Areas → InternalShapes
//! ```
//!
//! * `PerformShapesToAvoid` strips the faces that cannot belong to a closed
//!   boundary (a face whose boundary edge is used by a single face, or whose
//!   coincident copies wrap an edge without closing on the face).
//! * `PerformLoops` closes the surviving faces into shells
//!   ([`crate::shell_splitter::ShellSplitter`], OCCT `BOPAlgo_ShellSplitter`),
//!   then sets aside the faces that reached no loop and joins the set-aside
//!   faces through shared edges into internal shells.
//! * `PerformAreas` classifies each shell as a *growth* (the outer boundary of
//!   a region) or a *hole* (a cavity lying inside a growth) and assembles the
//!   regions into solids: every growth becomes one solid and receives the holes
//!   it contains.
//! * `PerformInternalShapes` classifies the internal faces against the regions
//!   with `BOPAlgo_Tools::ClassifyFaces` and adds IN faces as INTERNAL shells.
//!
//! ## Translation boundaries vs OCCT
//!
//! * `IsHole` is `BRepClass3d_SolidClassifier::PerformInfinitePoint` on the
//!   shell wrapped as a solid (`BuilderSolid.cxx:823`). `IsGrowthShell`
//!   (`:864`) treats a later shell as growth when it shares a face with an
//!   already-identified hole. Hole-to-growth assignment still uses
//!   `BOPTools_BoxTree` (`:460-502`) and `IsInside` (`:520`).
//! * `IntTools_Context::IsInfiniteFace` (open bounding box) is ported through
//!   the `BndBox` open flags.
//! * `BOPAlgo_ShellSplitter` is the existing
//!   [`crate::shell_splitter::ShellSplitter`], and the shared-section-face
//!   closing (`crate::bop_build_solids::close_open_shells`) is applied inside
//!   `PerformLoops` — the pipeline shares the section-face `TShape` between
//!   adjacent pieces where OCCT gives each piece its own face.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::abs::Orientation;
use crate::algo_tools::AlgoTools;
use crate::bbox_from_geometry::shape_bbox;
use crate::bop_aabb_faces::AabbTree;
use crate::bop_classify_occt::classify_faces_occt;
use crate::bop_occt_util::{iter_children, shape_key as occt_shape_key};
use crate::brep_extrema::{closest_point_on_face, is_inside};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Face, Shell, Solid, TopoShape};
use crate::shell_splitter::{ShellSplitter, edge_key, EKey};
use crate::topo_tools_full::{edges_of, faces_of};

/// Stable identity key of a shape (the address of its shared `TShape`).
fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Undirected identity keys of the boundary edges of `face`.
fn face_edge_keys(face: &TopoShape) -> HashSet<EKey> {
    edges_of(face).into_iter().map(|e| edge_key(&e)).collect()
}

/// Whether `edge_key` closes on `face` — the edge is used twice by the face
/// boundary (OCCT `BRep_Tool::IsClosed(aE, aF)`).
fn edge_closed_on_face(face: &TopoShape, key: &EKey) -> bool {
    edges_of(face).iter().filter(|e| edge_key(e) == *key).count() >= 2
}

/// Whether a face is unbounded (its bounding box is open in some direction).
/// OCCT `IntTools_Context::IsInfiniteFace`.
fn is_infinite_face(face: &TopoShape) -> bool {
    let b = shape_bbox(face);
    b.is_open_xmin()
        || b.is_open_xmax()
        || b.is_open_ymin()
        || b.is_open_ymax()
        || b.is_open_zmin()
        || b.is_open_zmax()
}

/// A 3-D point on `face`: the surface centre when its UV range is bounded,
/// else the midpoint of the first boundary edge.
fn face_sample_point(face: &Face) -> Option<GpPnt> {
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && v1.is_finite() {
        let s = BRepTool::face_surface(face)?;
        return Some(s.d0(0.5 * (u1 + u2), 0.5 * (v1 + v2)));
    }
    let e = edges_of(&face.0).into_iter().next()?;
    let (a, b) = BRepTool::edge_parameters(&e);
    if a.is_finite() && b.is_finite() {
        let c = BRepTool::edge_curve(&e)?;
        Some(c.d0(0.5 * (a + b)))
    } else {
        None
    }
}

/// State of a 3-D point relative to a solid: `On` within `tol` of a boundary
/// face, else `In`/`Out` by the even-odd parity test. Mirrors
/// `BOPTools_AlgoTools::ComputeState` reduced to a point.
fn classify_point(solid: &TopoShape, p: &GpPnt, tol: f64) -> FaceState {
    let mut d = f64::INFINITY;
    for f in faces_of(solid) {
        let (_, q) = closest_point_on_face(&f, p, 16, 16);
        d = d.min(q.distance(p));
    }
    if d.is_finite() && d <= tol {
        return FaceState::On;
    }
    if is_inside(solid, p) {
        FaceState::In
    } else {
        FaceState::Out
    }
}

/// Wrap a shell into a solid.
fn make_solid_from_shell(shell: &TopoShape) -> TopoShape {
    let bld = TopoBuilder::new();
    let mut solid = Solid::new();
    bld.add(&mut solid.0, shell);
    solid.0
}

/// Whether solid `inner` lies inside solid `outer` (`IsInside` of two solids,
/// `BOPAlgo_BuilderSolid.cxx:520`): a face point of `inner` classified IN
/// `outer`. Used only to pick the innermost hole owner, not as `IsHole`.
fn solid_inside_solid(inner: &TopoShape, outer: &TopoShape, tol: f64) -> bool {
    let Some(f) = faces_of(inner).into_iter().next() else {
        return false;
    };
    let Some(p) = face_sample_point(&f) else {
        return false;
    };
    classify_point(outer, &p, tol) == FaceState::In
}

/// `IsGrowthShell` (`BuilderSolid.cxx:864`): a later shell is a growth when
/// one of its direct children is a face of an already-identified hole.
fn is_growth_shell(shell: &TopoShape, hole_faces: &HashSet<usize>) -> bool {
    if hole_faces.is_empty() {
        return false;
    }
    iter_children(shell)
        .iter()
        .any(|c| hole_faces.contains(&shape_key(c)))
}

/// `IsHole` (`BuilderSolid.cxx:823`): infinite point of the shell-as-solid
/// classifies `IN`.
fn is_hole_shell(shell: &TopoShape) -> bool {
    AlgoTools::is_inverted_solid(&make_solid_from_shell(shell))
}

/// Whether the whole shell lies strictly inside the solid `solid`: a
/// representative face of the shell is sampled and classified `In`
/// (`IsInside`, OCCT `BOPTools_AlgoTools::ComputeState(Face, Solid)`).
fn shell_inside_solid(shell: &TopoShape, solid: &TopoShape, tol: f64) -> bool {
    let Some(f) = faces_of(shell).into_iter().next() else { return false };
    let Some(p) = face_sample_point(&f) else { return false };
    classify_point(solid, &p, tol) == FaceState::In
}

/// Join `faces` into shells through shared edges, each shell a connected group.
///
/// Mirrors the internal-shell builders of `BOPAlgo_BuilderSolid`: the
/// `PerformLoops` post-treatment (faces added as-is) and `MakeInternalShells`
/// (faces forced `INTERNAL`, as used by `PerformInternalShapes`).
fn connect_faces_into_shells(faces: &[TopoShape], force_internal: bool) -> Vec<TopoShape> {
    let n = faces.len();
    if n == 0 {
        return Vec::new();
    }
    // aEFMap: edge -> faces containing it.
    let face_keys: Vec<HashSet<EKey>> =
        faces.iter().map(|f| face_edge_keys(f)).collect();
    let mut edge_faces: HashMap<EKey, Vec<usize>> = HashMap::new();
    for (i, keys) in face_keys.iter().enumerate() {
        for &k in keys {
            edge_faces.entry(k).or_default().push(i);
        }
    }
    let bld = TopoBuilder::new();
    let mut added = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if added[start] {
            continue;
        }
        added[start] = true;
        let mut shell = Shell::new();
        let mut f = faces[start].clone();
        if force_internal {
            f.set_orientation(Orientation::Internal);
        }
        bld.add(&mut shell.0, &f);
        let mut queue = vec![start];
        while let Some(i) = queue.pop() {
            for &k in &face_keys[i] {
                if let Some(neigh) = edge_faces.get(&k) {
                    for &j in neigh {
                        if !added[j] {
                            added[j] = true;
                            let mut f = faces[j].clone();
                            if force_internal {
                                f.set_orientation(Orientation::Internal);
                            }
                            bld.add(&mut shell.0, &f);
                            queue.push(j);
                        }
                    }
                }
            }
        }
        shell.0.set_closed(true);
        out.push(shell.0);
    }
    out
}

/// The `BOPAlgo_BuilderSolid` state machine over a set of faces.
#[derive(Debug, Default)]
pub struct BuilderSolid {
    /// The input faces (`myShapes`).
    shapes: Vec<TopoShape>,
    /// Faces excluded from the boundary by `PerformShapesToAvoid`
    /// (`myShapesToAvoid`).
    shapes_to_avoid: Vec<TopoShape>,
    /// Closed shells of the boundary (`myLoops`).
    loops: Vec<TopoShape>,
    /// Internal shells reconstructed from the avoided faces (`myLoopsInternal`).
    loops_internal: Vec<TopoShape>,
    /// The region solids (`myAreas`).
    areas: Vec<TopoShape>,
    /// `myBoxes` — bounding box of each area / unbound hole solid.
    boxes: HashMap<usize, BndBox>,
    /// When true, internal parts are not merged back into the result
    /// (`myAvoidInternalShapes`).
    avoid_internal_shapes: bool,
    /// The tolerance of the classification (`myFuzzyValue`).
    fuzzy: f64,
    /// Non-fatal diagnostics.
    warnings: Vec<String>,
}

impl BuilderSolid {
    /// An empty builder with default options (fuzzy `1e-7`, internal shapes on).
    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            shapes_to_avoid: Vec::new(),
            loops: Vec::new(),
            loops_internal: Vec::new(),
            areas: Vec::new(),
            boxes: HashMap::new(),
            avoid_internal_shapes: false,
            fuzzy: 1e-7,
            warnings: Vec::new(),
        }
    }

    /// Sets the input faces (`SetShapes`).
    pub fn set_shapes(&mut self, shapes: Vec<TopoShape>) {
        self.shapes = shapes;
    }

    /// Sets the classification tolerance (`SetFuzzyValue`).
    pub fn set_fuzzy(&mut self, v: f64) {
        self.fuzzy = if v.is_nan() { 1e-7 } else { v.max(1e-7) };
    }

    /// Enables/disables the internal-shapes reconstruction
    /// (`SetAvoidInternalShapes`).
    pub fn set_avoid_internal_shapes(&mut self, v: bool) {
        self.avoid_internal_shapes = v;
    }

    /// The region solids built by the pipeline (`Areas`).
    pub fn areas(&self) -> &[TopoShape] {
        &self.areas
    }

    /// The closed boundary shells found by `PerformLoops`.
    pub fn loops(&self) -> &[TopoShape] {
        &self.loops
    }

    /// The faces set aside by `PerformShapesToAvoid`.
    pub fn shapes_to_avoid(&self) -> &[TopoShape] {
        &self.shapes_to_avoid
    }

    /// The internal shells reconstructed from the avoided faces.
    pub fn loops_internal(&self) -> &[TopoShape] {
        &self.loops_internal
    }

    /// The non-fatal diagnostics of the run.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Runs the `ShapesToAvoid → Loops → Areas → InternalShapes` sequence
    /// (`BOPAlgo_BuilderSolid::Perform`).
    pub fn perform(&mut self) -> Result<(), String> {
        self.perform_shapes_to_avoid();
        self.perform_loops()?;
        self.perform_areas();
        self.perform_internal_shapes();
        Ok(())
    }

    /// Whether `s` is currently marked to be avoided (same `TShape` identity).
    fn is_avoided(&self, s: &TopoShape) -> bool {
        self.shapes_to_avoid.iter().any(|a| a.same_tshape(s))
    }

    /// Phase 1 — strip the faces that cannot belong to a closed boundary.
    ///
    /// Mirrors `BOPAlgo_BuilderSolid::PerformShapesToAvoid`: repeatedly build
    /// the edge→faces map (`TopExp::MapShapesAndAncestors`) of the faces not
    /// yet avoided and mark a face when one of its boundary edges is used by a
    /// single face (a dangling leaf), or by its own two coincident copies that
    /// do not close on the face. Iterates until no new face is marked.
    fn perform_shapes_to_avoid(&mut self) {
        self.shapes_to_avoid.clear();
        let faces: Vec<TopoShape> =
            self.shapes.iter().filter(|s| s.is_face()).cloned().collect();
        loop {
            // MEF: edge key -> (face indices, any degenerated, any INTERNAL).
            let mut mef: HashMap<EKey, (Vec<usize>, bool, bool)> = HashMap::new();
            for (i, f) in faces.iter().enumerate() {
                if self.is_avoided(f) {
                    continue;
                }
                for e in edges_of(f) {
                    let key = edge_key(&e);
                    let entry = mef.entry(key).or_default();
                    entry.0.push(i);
                    entry.1 |= BRepTool::is_degenerated(&e);
                    entry.2 |= e.orientation() == Orientation::Internal;
                }
            }
            let mut found = false;
            for (_key, (lf, degenerated, has_internal)) in mef {
                if degenerated || lf.is_empty() {
                    continue;
                }
                if lf.len() == 1 {
                    if has_internal {
                        continue;
                    }
                    let f1 = faces[lf[0]].clone();
                    if !self.is_avoided(&f1) {
                        self.shapes_to_avoid.push(f1);
                        found = true;
                    }
                } else if lf.len() == 2 {
                    let f1 = faces[lf[0]].clone();
                    let f2 = faces[lf[1]].clone();
                    if f1.same_tshape(&f2) {
                        if edge_closed_on_face(&f1, &_key) {
                            continue;
                        }
                        if has_internal {
                            continue;
                        }
                        if !self.is_avoided(&f1) {
                            self.shapes_to_avoid.push(f1);
                            found = true;
                        }
                        if !self.is_avoided(&f2) {
                            self.shapes_to_avoid.push(f2);
                            found = true;
                        }
                    }
                }
            }
            if !found {
                break;
            }
        }
    }

    /// Phase 2 — close the surviving faces into shells.
    ///
    /// Mirrors `BOPAlgo_BuilderSolid::PerformLoops`: unbounded faces become
    /// single-face shells, the rest go to the shell splitter; the splitter's
    /// shells are closed over the shared-section-face artifact
    /// ([`crate::bop_build_solids::close_open_shells`]); faces that reached no
    /// loop are set aside; and the set-aside faces are joined through shared
    /// edges into the internal shells.
    fn perform_loops(&mut self) -> Result<(), String> {
        self.loops.clear();
        self.loops_internal.clear();
        let bld = TopoBuilder::new();
        // 1. Shells usual.
        let mut splitter = ShellSplitter::new();
        for f in &self.shapes {
            if !f.is_face() {
                continue;
            }
            if is_infinite_face(f) {
                let mut shell = Shell::new();
                bld.add(&mut shell.0, f);
                shell.0.set_closed(true);
                self.loops.push(shell.0);
                continue;
            }
            if !self.is_avoided(f) {
                splitter.add_start_element(f.clone());
            }
        }
        splitter.perform()?;
        let shells = crate::bop_build_solids::close_open_shells(
            splitter.shells(),
            &self.shapes,
        );
        for sh in shells {
            if crate::bop_build_solids::geometrically_open(&sh) {
                continue;
            }
            self.loops.push(sh);
        }
        // 2. Post treatment: faces of the input that reached no loop are set
        //    aside (OCCT appends them to myShapesToAvoid).
        let in_loops: Vec<TopoShape> = self
            .loops
            .iter()
            .flat_map(|l| faces_of(l).into_iter().map(|f| f.0))
            .collect();
        for f in &self.shapes {
            if f.is_face()
                && !is_infinite_face(f)
                && !self.is_avoided(f)
                && !in_loops.iter().any(|x| x.same_tshape(f))
            {
                self.shapes_to_avoid.push(f.clone());
            }
        }
        // 3. Internal shells from the faces set aside.
        let avoided: Vec<TopoShape> = self
            .shapes_to_avoid
            .iter()
            .filter(|s| s.is_face())
            .cloned()
            .collect();
        self.loops_internal = connect_faces_into_shells(&avoided, false);
        Ok(())
    }

    /// Phase 3 — turn the shells into region solids.
    ///
    /// Mirrors `BOPAlgo_BuilderSolid::PerformAreas` (`:397`).
    /// `IsGrowthShell` then `IsHole` (`:422-427`); hole-to-growth assignment
    /// uses `BOPTools_BoxTree` (`:460-502`) and `IsInside` (`:520`). Unbound
    /// holes get a `SetWhole` box (`:591`).
    fn perform_areas(&mut self) {
        self.areas.clear();
        self.boxes.clear();
        let n = self.loops.len();
        if n == 0 {
            return;
        }
        let tol = self.fuzzy.max(1e-7);
        let solids: Vec<TopoShape> = self.loops.iter().map(make_solid_from_shell).collect();
        let mut hole_faces: HashSet<usize> = HashSet::new();
        let mut growth: Vec<usize> = Vec::new();
        let mut hole_idx: Vec<usize> = Vec::new();
        for i in 0..n {
            let mut b_growth = is_growth_shell(&self.loops[i], &hole_faces);
            if !b_growth {
                b_growth = !is_hole_shell(&self.loops[i]);
            }
            if b_growth {
                growth.push(i);
            } else {
                hole_idx.push(i);
                for c in iter_children(&self.loops[i]) {
                    hole_faces.insert(shape_key(&c));
                }
            }
        }
        let bld = TopoBuilder::new();

        if hole_idx.is_empty() {
            for &j in &growth {
                let s = solids[j].clone();
                self.boxes.insert(occt_shape_key(&s), shape_bbox(&s));
                self.areas.push(s);
            }
            return;
        }

        let mut tree = AabbTree::new();
        tree.set_size(hole_idx.len());
        for (k, &i) in hole_idx.iter().enumerate() {
            tree.add(k, shape_bbox(&self.loops[i]));
        }
        tree.build();

        let mut hole_owner: Vec<Option<usize>> = vec![None; n];
        for &j in &growth {
            let box_s = shape_bbox(&solids[j]);
            for k in tree.select(&box_s) {
                if k >= hole_idx.len() {
                    continue;
                }
                let i = hole_idx[k];
                if !shell_inside_solid(&self.loops[i], &solids[j], tol) {
                    continue;
                }
                hole_owner[i] = match hole_owner[i] {
                    None => Some(j),
                    Some(o) => {
                        if solid_inside_solid(&solids[j], &solids[o], tol) {
                            Some(j)
                        } else {
                            Some(o)
                        }
                    }
                };
            }
        }

        for &j in &growth {
            let mut solid = solids[j].clone();
            for &i in &hole_idx {
                if hole_owner[i] == Some(j) {
                    bld.add(&mut solid, &self.loops[i]);
                }
            }
            self.boxes.insert(occt_shape_key(&solid), shape_bbox(&solid));
            self.areas.push(solid);
        }
        for &i in &hole_idx {
            if hole_owner[i].is_none() {
                let s = solids[i].clone();
                let mut whole = BndBox::new();
                whole.set_whole();
                self.boxes.insert(occt_shape_key(&s), whole);
                self.areas.push(s);
            }
        }
    }

    /// Phase 4 — classify the internal shells against the regions.
    ///
    /// Mirrors `BOPAlgo_BuilderSolid::PerformInternalShapes` (`:602`): with no
    /// regions the internal faces alone form the solid (`MakeInternalShells`);
    /// otherwise `BOPAlgo_Tools::ClassifyFaces` (`:673`) maps each area to its
    /// IN faces, which become INTERNAL shells. Faces that fit no region raise
    /// `BOPAlgo_AlertSolidBuilderUnusedFaces` and are not added to the result.
    fn perform_internal_shapes(&mut self) {
        if self.avoid_internal_shapes {
            return;
        }
        if self.loops_internal.is_empty() {
            return;
        }
        // Shell children (`TopoDS_Iterator`), not a uniquifying face explorer.
        let mut faces: Vec<TopoShape> = Vec::new();
        for sh in &self.loops_internal {
            for child in iter_children(sh) {
                if !faces.iter().any(|x| x.same_tshape(&child)) {
                    faces.push(child);
                }
            }
        }
        if faces.is_empty() {
            return;
        }
        let bld = TopoBuilder::new();
        if self.areas.is_empty() {
            let mut solid = Solid::new();
            for sh in connect_faces_into_shells(&faces, true) {
                bld.add(&mut solid.0, &sh);
            }
            self.areas.push(solid.0);
            return;
        }

        let areas = self.areas.clone();
        let mut box_map = HashMap::new();
        for a in &areas {
            let k = occt_shape_key(a);
            let box_ = self.boxes.get(&k).copied().unwrap_or_else(|| shape_bbox(a));
            box_map.insert(k, box_);
        }
        let ctx = IntToolsContext::new();
        let empty_if: HashMap<usize, Vec<TopoShape>> = HashMap::new();
        let in_parts = classify_faces_occt(&faces, &areas, &ctx, &box_map, &empty_if);

        let mut done: HashSet<usize> = HashSet::new();
        for (k, area) in areas.iter().enumerate() {
            let Some(list) = in_parts.get(&occt_shape_key(area)) else {
                continue;
            };
            if list.is_empty() {
                continue;
            }
            for f in list {
                done.insert(occt_shape_key(f));
            }
            let mut solid = self.areas[k].clone();
            for sh in connect_faces_into_shells(list, true) {
                bld.add(&mut solid, &sh);
            }
            self.areas[k] = solid;
        }
        let unused = faces
            .iter()
            .filter(|f| !done.contains(&occt_shape_key(f)))
            .count();
        if unused > 0 {
            self.warnings.push(format!(
                "BuilderSolid: {unused} internal face(s) fit no region and were dropped"
            ));
        }
    }
}

/// Builds the region solids from `faces` — the drop-in equivalent of the
/// `BOPAlgo_BuilderSolid::Perform` on a set of (already split) faces.
pub fn build_solids_from_faces(faces: &[TopoShape]) -> Result<Vec<TopoShape>, String> {
    let mut bs = BuilderSolid::new();
    bs.set_shapes(faces.to_vec());
    bs.perform()?;
    Ok(bs.areas().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;
    use crate::brep_extrema::test_box::unit_box;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt};

    /// The six faces of a unit box (from the shared fixture).
    fn box_faces() -> Vec<TopoShape> {
        faces_of(&unit_box().solid.0).into_iter().map(|f| f.0).collect()
    }

    #[test]
    fn closed_box_faces_produce_one_solid() {
        // The six connected faces of a box close into one growth shell and one
        // region solid; nothing is avoided, no internal shells.
        let faces = box_faces();
        let mut bs = BuilderSolid::new();
        bs.set_shapes(faces);
        bs.perform().expect("builder solid perform");
        assert!(bs.shapes_to_avoid().is_empty(), "no dangling face in a closed box");
        assert!(bs.loops_internal().is_empty());
        assert_eq!(bs.areas().len(), 1, "six connected faces -> one region solid");
        assert!(bs.areas()[0].is_solid());
        assert_eq!(faces_of(&bs.areas()[0]).len(), 6);
    }

    #[test]
    fn dangling_face_is_avoided_and_ignored() {
        // A box plus a small quad floating above it: the quad shares no edge
        // with the box, so PerformShapesToAvoid strips it and only the box is
        // rebuilt as the single region.
        let mut faces = box_faces();
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(0.0, 0.0, 2.0),
            GpPnt::new(0.5, 0.0, 2.0),
            GpPnt::new(0.5, 0.5, 2.0),
            GpPnt::new(0.0, 0.5, 2.0),
        ];
        let wire = b.make_wire(&[
            b.make_edge_segment(&corners[0], &corners[1]),
            b.make_edge_segment(&corners[1], &corners[2]),
            b.make_edge_segment(&corners[2], &corners[3]),
            b.make_edge_segment(&corners[3], &corners[0]),
        ]);
        let ax = GpAx3::new(
            GpPnt::new(0.0, 0.0, 2.0),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ).unwrap();
        let mut quad = b.make_face_plane(&GpPln::new(ax));
        b.add_wire(&mut quad, &wire);
        faces.push(quad.0);
        let mut bs = BuilderSolid::new();
        bs.set_shapes(faces);
        bs.perform().expect("perform");
        assert!(
            bs.shapes_to_avoid().iter().any(|a| a.shape_type() == ShapeType::Face),
            "a dangling face is avoided"
        );
        assert_eq!(bs.areas().len(), 1, "the box remains the single region");
    }

    #[test]
    fn inner_cavity_is_absorbed_as_a_hole() {
        // A closed outer box shell and a smaller closed box shell strictly
        // inside it: the inner shell is a hole of the outer one, so the assembly
        // yields ONE solid containing both shells (not two separate solids).
        let outer = unit_box(); // [0,1]^3
        let inner = BRepPrimBox::make_box_corner(
            &GpPnt::new(0.25, 0.25, 0.25),
            &GpPnt::new(0.75, 0.75, 0.75),
        );
        let mut faces: Vec<TopoShape> = faces_of(&outer.solid.0).into_iter().map(|f| f.0).collect();
        // A cavity shell is stored INVERTED: `IsHole`
        // (`BOPAlgo_BuilderSolid.cxx:823-831`) asks whether the *infinite* point
        // classifies `IN` on the shell wrapped as a solid, which only holds when
        // the shell's faces point into its own void. `BRepPrimBox::make_box_corner`
        // builds a normal outward box, so the cavity faces are reversed here —
        // exactly how the split of a Cut operand hands its inner shell to the
        // builder. (With outward faces OCCT sees a second *growth* and emits two
        // areas, which is what the un-reversed fixture was really asserting.)
        faces.extend(
            faces_of(&inner.solid.0)
                .into_iter()
                .map(|f| f.0.oriented(crate::abs::Orientation::Reversed)),
        );
        let mut bs = BuilderSolid::new();
        bs.set_shapes(faces);
        bs.perform().expect("perform");
        assert_eq!(bs.areas().len(), 1, "the cavity is absorbed into the outer solid");
        let area = &bs.areas()[0];
        assert_eq!(faces_of(area).len(), 12, "outer 6 + inner 6 shells in one solid");
    }

    #[test]
    fn build_solids_from_faces_two_disjoint_boxes() {
        // Two fully disjoint box shells: neither lies inside the other, so both
        // become their own region solid.
        let a = unit_box();
        let b = BRepPrimBox::make_box_corner(
            &GpPnt::new(3.0, 0.0, 0.0),
            &GpPnt::new(4.0, 1.0, 1.0),
        );
        let mut faces: Vec<TopoShape> = faces_of(&a.solid.0).into_iter().map(|f| f.0).collect();
        faces.extend(faces_of(&b.solid.0).into_iter().map(|f| f.0));
        let solids = build_solids_from_faces(&faces).expect("build");
        assert_eq!(solids.len(), 2, "two disjoint shells -> two region solids");
    }
}
