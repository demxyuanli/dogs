//! Pave-block building and edge splitting — Phase 19 wave C2b.
//!
//! Port of the `BOPAlgo_PaveFiller` block-building methods:
//!
//! | OCCT method                   | Rust function                          |
//! |-------------------------------|----------------------------------------|
//! | `MakeBlocks`                  | [`make_blocks`]                        |
//! | `MakeBlocks` (full)          | [`make_blocks_full`]                   |
//! | `FilterPavesOnCurves`         | [`filter_paves_on_curves`]             |
//! | `MakeSplitEdges`              | [`make_split_edges`]                   |
//! | `MakePCurves`                 | [`make_pcurves`]                       |
//! | `FindPaveBlocks`              | [`find_pave_blocks`]                   |
//! | `FillPaves`                   | [`fill_paves`]                         |
//! | `SplitEdge`                   | [`make_split_edge`]                    |
//! | `MakeSplitEdge` (ProcessDE)   | [`make_split_edge_de`]                 |
//! | `ProcessDE`                   | [`process_de`]                         |
//!
//! `MakeBlocks` in OCCT builds the face/face section edges (the `PerformFF`
//! step). That step is **not** part of this phase, so [`make_blocks`] is the
//! reduced form: it makes sure every source edge has a default pave block,
//! orders each edge's blocks by their parametric range, then splits any block
//! carrying extra paves (the OCCT `BOPDS_DS::UpdatePaveBlocks`). The FF-only
//! helpers (`CorrectToleranceOfSE`, `PutPavesOnCurve`, …) are omitted.
//!
//! Every function operates on a [`PaveFillerLike`] — the minimal contract the
//! concrete `crate::pave_filler::PaveFiller` (filled by a parallel agent)
//! satisfies. The port keeps the logic decoupled from that type so the two
//! can be developed and verified independently.
mod prelude {

pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION};
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use occt_geom2d::geom2d_api::{intersect_curves, project_point_on_curve};
pub(crate) use occt_geom2d::Geom2dBSplineCurve;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::algo_tools::AlgoTools;
pub(crate) use crate::bopds::{BopdsDS, BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
pub(crate) use crate::boptools_2d;
pub(crate) use crate::brep_projection::project_point_on_face;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::pcurve_full;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex};
pub(crate) use crate::tgeometry::GeometryRegistry;

}

use prelude::*;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;

    /// Minimal filler satisfying [`PaveFillerLike`] for isolated tests.
    #[derive(Debug, Default)]
    struct StubFiller {
        ds: BopdsDS,
        errors: Vec<String>,
        warnings: Vec<String>,
    }

    impl PaveFillerLike for StubFiller {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn ds_mut(&mut self) -> &mut BopdsDS {
            &mut self.ds
        }
        fn add_error(&mut self, msg: String) {
            self.errors.push(msg);
        }
        fn add_warning(&mut self, msg: String) {
            self.warnings.push(msg);
        }
    }

    /// A stub DS holding a unit box whose edge 0 runs (0,0,0) → (1,0,0).
    /// Uses `init` so the source/new-shape split (`nb_source_shapes`) is set.
    fn box_filler() -> (StubFiller, usize) {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        (f, e0)
    }

    /// Append a fresh vertex to the DS and return its index.
    fn push_vertex(f: &mut StubFiller, p: GpPnt) -> usize {
        let v = AlgoTools::make_new_vertex(&p, 1e-7).unwrap();
        f.ds.append(v).unwrap()
    }

    /// Build a pave block on `edge` bounded by `v1`/`v2` over `(t1, t2)`.
    fn mk_block(edge: usize, v1: usize, t1: f64, v2: usize, t2: f64) -> BopdsPaveBlock {
        let mut pb = BopdsPaveBlock::new();
        pb.set_edge(edge);
        pb.set_original_edge(edge);
        pb.set_pave1(BopdsPave::new(v1, t1));
        pb.set_pave2(BopdsPave::new(v2, t2));
        pb
    }

    /// A filler whose DS holds a unit box, with default pave blocks created on
    /// the box edges at the given indices. Returns the filler and the DS
    /// indices of the requested edges.
    fn box_edges(indices: &[usize]) -> (StubFiller, Vec<usize>) {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let inds: Vec<usize> = indices.iter().map(|&i| f.ds.index(&b.edges[i].0).unwrap()).collect();
        for &e in &inds {
            f.ds.init_pave_blocks_for_edge(e);
        }
        (f, inds)
    }

    #[test]
    fn make_blocks_splits_edge_with_two_extra_paves() {
        let (mut f, e0) = box_filler();
        // Edge 0 has one default block [0, 1]. Add two intersection paves.
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            assert_eq!(blocks.len(), 1);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
            assert!(blocks[0].is_to_update());
        }

        make_blocks(&mut f).unwrap();

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3, "expected 3 elementary blocks");
        assert_eq!(blocks[0].range(), (0.0, 0.25));
        assert_eq!(blocks[1].range(), (0.25, 0.75));
        assert_eq!(blocks[2].range(), (0.75, 1.0));
        // The middle block is bounded by the two new vertices.
        assert_eq!(blocks[1].indices(), (i1, i2));
    }

    #[test]
    fn make_split_edge_creates_sub_edge_with_correct_range() {
        let (mut f, e0) = box_filler();
        let v1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let v2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        let n_sp = make_split_edge(&mut f, e0, v1, 0.25, v2, 0.75).expect("split edge");
        assert!(n_sp >= f.ds.nb_source_shapes(), "new edge is a new shape");
        let new_shape = f.ds.shape(n_sp).cloned().expect("new edge shape");
        assert_eq!(new_shape.shape_type(), ShapeType::Edge);
        assert_eq!(BRepTool::edge_parameters(&Edge(new_shape)), (0.25, 0.75));
        // Sub-shapes of the new edge are the two bound vertices.
        let si = f.ds.shape_info(n_sp).unwrap();
        assert_eq!(si.sub_shapes(), &[v1, v2]);
        // 3D endpoints land on the original curve.
        let p0 = AlgoTools::point_on_edge(&Edge(f.ds.shape(n_sp).cloned().unwrap()), 0.25).unwrap();
        let p1 = AlgoTools::point_on_edge(&Edge(f.ds.shape(n_sp).cloned().unwrap()), 0.75).unwrap();
        assert!(p0.distance(&GpPnt::new(0.25, 0.0, 0.0)) < 1e-9, "start {p0:?}");
        assert!(p1.distance(&GpPnt::new(0.75, 0.0, 0.0)) < 1e-9, "end {p1:?}");
    }

    #[test]
    fn make_split_edges_splits_every_block_onto_new_edge() {
        let (mut f, e0) = box_filler();
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
        }
        make_blocks(&mut f).unwrap();
        let before = f.ds.nb_shapes();
        make_split_edges(&mut f).unwrap();
        assert_eq!(f.ds.nb_shapes(), before + 3, "3 sub-edges appended");

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        let expected = [(0.0, 0.25), (0.25, 0.75), (0.75, 1.0)];
        for (pb, (lo, hi)) in blocks.iter().zip(expected.iter()) {
            assert_eq!(pb.range(), (*lo, *hi));
            let n_sp = pb.edge();
            assert!(n_sp != e0, "block re-pointed to a new edge");
            let shape = f.ds.shape(n_sp).cloned().expect("split edge shape");
            assert_eq!(BRepTool::edge_parameters(&Edge(shape)), (*lo, *hi));
            assert_eq!(pb.original_edge(), e0, "original edge preserved");
        }
    }

    #[test]
    fn make_pcurves_stores_pcurve_in_registry() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Register face 0 with edge 0 lying on it.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        assert!(!boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());

        make_pcurves(&mut f).unwrap();

        let pc = boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).expect("pcurve stored");
        // Face 0 (bottom, z = 0) parameterizes (u, v) → (v, u, 0): its u axis
        // runs +Y, so edge 0 (0,0,0)→(1,0,0) is a v-line at u = 0.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        let q1 = pc.d0(1.0);
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
    }

    #[test]
    fn make_pcurves_builds_pcurve_for_in_and_on_edges() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Edge 0 has no pcurve on face 0 yet.
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_none());
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave_on(e0, 0.0, 1.0);
            fi.add_pave_in(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }

        make_pcurves(&mut f).unwrap();

        // The IN/ON paves got a pcurve (deduplicated: the same edge appears in
        // both lists, the pcurve is built once).
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some(), "IN/ON pcurve built");
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_builds_on_pcurve_once_when_already_present() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Edge 0 carries a pcurve on face 0 already.
        let pc = pcurve_full::make_pcurve_full(&b.edges[0], &b.faces[0]).unwrap();
        let key = GeometryRegistry::shape_key(&b.faces[0].0);
        GeometryRegistry::global().set_edge_pcurve(&b.edges[0].0, key, pc);
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave_on(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }

        make_pcurves(&mut f).unwrap();

        // The ON block whose edge already has a pcurve is skipped (OCCT bHasPC).
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_updates_vertices_of_section_edge_with_existing_pcurve() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Store a deliberately wrong pcurve for edge 0 on face 0 (shifted by
        // (0.1, 0.1) in UV) — simulating a section pcurve that deviates from
        // the true surface mapping.
        let true_pc = pcurve_full::make_pcurve_full(&b.edges[0], &b.faces[0]).unwrap();
        let mut tr = occt_core::gp::GpTrsf2d::identity();
        tr.set_translation_vec(&occt_core::gp::GpVec2d::new(0.1, 0.1));
        let shifted: Arc<dyn Curve2d> = Arc::from(true_pc.transformed(&tr));
        let face_key = GeometryRegistry::shape_key(&b.faces[0].0);
        GeometryRegistry::global().set_edge_pcurve(&b.edges[0].0, face_key, shifted);
        // Register the edge as a section pave on the face.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        let v0 = edge_vertex_shapes(&b.edges[0])[0].clone();
        let v1 = edge_vertex_shapes(&b.edges[0])[1].clone();
        let tol_before = BRepTool::vertex_tolerance(&Vertex(v0.clone()));

        make_pcurves(&mut f).unwrap();

        // The section branch calls UpdateVertices even though the pcurve already
        // exists: the boundary vertex tolerances grow to cover the 3D/2D gap.
        let tol_after = BRepTool::vertex_tolerance(&Vertex(v0.clone()));
        assert!(tol_after > tol_before, "vertex tolerance grew from {tol_before} to {tol_after}");
        let tol_after1 = BRepTool::vertex_tolerance(&Vertex(v1.clone()));
        assert!(tol_after1 > tol_before, "second vertex tolerance grew: {tol_after1}");
    }

    #[test]
    fn filter_paves_dedupes_bound_coincident_extra_pave() {
        let (mut f, e0) = box_filler();
        let i1 = push_vertex(&mut f, GpPnt::new(0.25, 0.0, 0.0));
        let i2 = push_vertex(&mut f, GpPnt::new(0.75, 0.0, 0.0));
        // Split into three blocks.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i1, 0.25));
            blocks[0].append_ext_pave(BopdsPave::new(i2, 0.75));
        }
        make_blocks(&mut f).unwrap();
        assert_eq!(f.ds.pave_blocks(e0).len(), 3);

        // Add a spurious extra pave at the exact upper bound of the first
        // block (a different vertex, so only the parameter collides).
        let i3 = push_vertex(&mut f, GpPnt::new(0.5, 0.5, 0.5));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(i3, 0.25));
            assert!(blocks[0].is_to_update());
        }
        filter_paves_on_curves(&mut f).unwrap();

        // The bound-coincident pave is dropped; the 3 elementary blocks remain
        // (their shared vertices are new, so no spurious merge happens).
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        for pb in blocks {
            assert!(pb.ext_paves().is_empty());
        }
    }

    #[test]
    fn filter_paves_merges_adjacent_touching_blocks() {
        let (mut f, e0) = box_filler();
        let v0 = f.ds.pave_blocks(e0)[0].indices().0;
        let v1 = f.ds.pave_blocks(e0)[0].indices().1;
        // Two adjacent blocks sharing the middle vertex/parameter.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            let mut a = BopdsPaveBlock::new();
            a.set_edge(e0);
            a.set_original_edge(e0);
            a.set_pave1(BopdsPave::new(v0, 0.0));
            a.set_pave2(BopdsPave::new(v1, 0.5));
            let mut b = BopdsPaveBlock::new();
            b.set_edge(e0);
            b.set_original_edge(e0);
            b.set_pave1(BopdsPave::new(v1, 0.5));
            b.set_pave2(BopdsPave::new(v0, 1.0));
            blocks.push(a);
            blocks.push(b);
        }
        filter_paves_on_curves(&mut f).unwrap();
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 1, "adjacent blocks merged");
        assert_eq!(blocks[0].range(), (0.0, 1.0));
        assert_eq!(blocks[0].indices(), (v0, v0));
    }

    #[test]
    fn make_blocks_full_groups_parallel_edges_with_same_range() {
        let (mut f, inds) = box_edges(&[0, 4]);
        let (e0, e4) = (inds[0], inds[1]);
        assert_ne!(e0, e4);
        // Edge 0 (bottom front, along X) and edge 4 (top front, along X) are
        // parallel and equal-length. Artificially force the same pave-block
        // range on both — the CommonBlock detection must group them.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].set_range(0.25, 0.75);
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks[0].set_range(0.25, 0.75);
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "expected one common block, got {}", cbs.len());
        let cb = &cbs[0];
        assert!(cb.contains_index(e0) && cb.contains_index(e4), "indices {:?}", cb.indices());
        assert!(cb.contains_range(0.25, 0.75, 1e-9), "ranges {:?}", cb.ranges());
    }

    #[test]
    fn make_blocks_full_does_not_group_different_ranges() {
        let (mut f, inds) = box_edges(&[0, 4]);
        let (e0, e4) = (inds[0], inds[1]);
        {
            // e4's block covers only the second half of its curve.
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks[0].set_range(0.5, 1.0);
        }

        make_blocks_full(&mut f).unwrap();

        assert!(
            f.ds.common_blocks().is_empty(),
            "no group expected, got {}",
            f.ds.common_blocks().len()
        );
        // The two edges still have their own (unmerged) block lists.
        assert_eq!(f.ds.pave_blocks(e0).len(), 1);
        assert_eq!(f.ds.pave_blocks(e4).len(), 1);
    }

    #[test]
    fn make_blocks_full_groups_by_shared_bound_vertices() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        f.ds.init_pave_blocks_for_edge(e4);
        // Give e0's block the SAME bound vertices as e4's block, but a
        // different parametric range — the vertex match alone must group them.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].set_indices(v4, v5);
            blocks[0].set_range(0.2, 0.8);
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "expected one common block, got {}", cbs.len());
        assert!(cbs[0].contains_index(e0) && cbs[0].contains_index(e4));
    }

    #[test]
    fn make_blocks_full_groups_three_edges_sharing_interval() {
        let (mut f, inds) = box_edges(&[0, 4, 6]);
        let (e0, e4, e6) = (inds[0], inds[1], inds[2]);
        // Edges 0, 4 and 6 are three parallel, equal-length box edges; their
        // default blocks all span [0, 1], so all three join one common block.
        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 1, "one common block for three coincident intervals");
        let cb = &cbs[0];
        assert!(
            cb.contains_index(e0) && cb.contains_index(e4) && cb.contains_index(e6),
            "indices {:?}",
            cb.indices()
        );
        assert!(cb.contains_range(0.0, 1.0, 1e-9));
    }

    #[test]
    fn make_blocks_full_records_two_intervals_on_same_edge() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v0 = f.ds.index(&b.vertices[0].0).unwrap();
        let v1 = f.ds.index(&b.vertices[1].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        f.ds.init_pave_blocks_for_edge(e4);
        // A new (intersection) vertex splits both edges at t = 0.5, so the two
        // halves are not mergeable and form two separate common intervals.
        let mid = push_vertex(&mut f, GpPnt::new(0.5, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, mid, 0.5));
            blocks.push(mk_block(e0, mid, 0.5, v1, 1.0));
            let blocks = f.ds_mut().change_pave_blocks_mut(e4);
            blocks.clear();
            blocks.push(mk_block(e4, v4, 0.0, mid, 0.5));
            blocks.push(mk_block(e4, mid, 0.5, v5, 1.0));
        }

        make_blocks_full(&mut f).unwrap();

        let cbs = f.ds.common_blocks();
        assert_eq!(cbs.len(), 2, "two intervals on the same edges, got {}", cbs.len());
        for cb in cbs {
            assert!(cb.contains_index(e0) && cb.contains_index(e4), "indices {:?}", cb.indices());
        }
        let has_first = cbs.iter().any(|cb| cb.contains_range(0.0, 0.5, 1e-9));
        let has_second = cbs.iter().any(|cb| cb.contains_range(0.5, 1.0, 1e-9));
        assert!(has_first && has_second, "ranges: {:?}", cbs.iter().map(|c| c.ranges().to_vec()).collect::<Vec<_>>());
    }

    #[test]
    fn make_blocks_full_merges_touching_blocks() {
        let (mut f, e0) = box_filler();
        let (v0, v1) = {
            let blocks = f.ds.pave_blocks(e0);
            (blocks[0].indices().0, blocks[0].indices().1)
        };
        // Two clean adjacent blocks sharing source vertex v1 at t = 0.5.
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, v1, 0.5));
            blocks.push(mk_block(e0, v1, 0.5, v0, 1.0));
        }

        make_blocks_full(&mut f).unwrap();

        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 1, "touching blocks merged into one");
        assert_eq!(blocks[0].range(), (0.0, 1.0));
    }

    #[test]
    fn make_blocks_full_does_not_merge_to_update_blocks() {
        let (mut f, e0) = box_filler();
        let (v0, v1) = {
            let blocks = f.ds.pave_blocks(e0);
            (blocks[0].indices().0, blocks[0].indices().1)
        };
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks.clear();
            blocks.push(mk_block(e0, v0, 0.0, v1, 0.5));
            // The second block carries an intersection split point, so it is
            // not merged away even though it touches the first one.
            let mut b = mk_block(e0, v1, 0.5, v0, 1.0);
            b.append_ext_pave(BopdsPave::new(v0, 0.75));
            blocks.push(b);
        }

        make_blocks_full(&mut f).unwrap();

        // No merge; the second block is split by update_pave_blocks.
        let blocks = f.ds.pave_blocks(e0);
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[0].range(), (0.0, 0.5));
        assert_eq!(blocks[1].range(), (0.5, 0.75));
        assert_eq!(blocks[2].range(), (0.75, 1.0));
    }

    #[test]
    fn blocks_match_handles_edges_intervals_and_orientation() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let e4 = f.ds.index(&b.edges[4].0).unwrap();
        let v0 = f.ds.index(&b.vertices[0].0).unwrap();
        let v1 = f.ds.index(&b.vertices[1].0).unwrap();
        let v4 = f.ds.index(&b.vertices[4].0).unwrap();
        let v5 = f.ds.index(&b.vertices[5].0).unwrap();
        let a = mk_block(e0, v0, 0.0, v1, 1.0);
        // Blocks on the same edge never match (a self-overlap is not a common
        // block between distinct edges).
        assert!(!blocks_match(&a, &mk_block(e0, v0, 0.0, v1, 1.0), PCONFUSION), "same edge");
        // Same interval on a different edge matches.
        assert!(blocks_match(&a, &mk_block(e4, v4, 0.0, v5, 1.0), PCONFUSION), "same interval");
        // Disjoint intervals do not match.
        assert!(!blocks_match(&a, &mk_block(e4, v4, 0.5, v5, 1.0), PCONFUSION), "different interval");
        // A reversed orientation still covers the same interval.
        assert!(blocks_match(&a, &mk_block(e4, v5, 1.0, v4, 0.0), PCONFUSION), "reversed interval");
        // The same pair of bound vertices groups regardless of range.
        assert!(blocks_match(&a, &mk_block(e4, v0, 0.2, v1, 0.8), PCONFUSION), "shared bounds");
    }

    #[test]
    fn make_pcurves_full_stores_trimmed_aligned_pcurve_on_box_edge() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_none());

        make_pcurves_full(&mut f).unwrap();

        let pc = boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).expect("pcurve stored");
        // Edge 0 runs (0,0,0)→(1,0,0) on the bottom face, whose (u,v)→(v,u,0)
        // parameterization maps the edge to the v-line u = 0, v ∈ [0, 1]. The
        // endpoints hit the face UV quad boundary.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        let q1 = pc.d0(1.0);
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
        // The endpoints equal the projection of the edge's 3D endpoints.
        let (u0, v0, _) = project_point_on_face(&b.faces[0], &GpPnt::new(0.0, 0.0, 0.0)).unwrap();
        assert!((q0.x() - u0).abs() < 1e-6 && (q0.y() - v0).abs() < 1e-6, "proj start");
        let (u1, v1, _) = project_point_on_face(&b.faces[0], &GpPnt::new(1.0, 0.0, 0.0)).unwrap();
        assert!((q1.x() - u1).abs() < 1e-6 && (q1.y() - v1).abs() < 1e-6, "proj end");
        // The trimmed pcurve lies within the face UV domain.
        let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(&b.faces[0]);
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            let q = pc.d0(t);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {q:?} at {t}");
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {q:?} at {t}");
        }
        assert!(f.warnings.is_empty(), "warnings: {:?}", f.warnings);
    }

    #[test]
    fn make_pcurves_full_covers_split_edges() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.init(&[b.solid.0.clone()]);
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        f.ds.init_pave_blocks_for_edge(e0);
        // Split edge 0 at t = 0.5 with an intersection vertex.
        let mid = push_vertex(&mut f, GpPnt::new(0.5, 0.0, 0.0));
        {
            let blocks = f.ds_mut().change_pave_blocks_mut(e0);
            blocks[0].append_ext_pave(BopdsPave::new(mid, 0.5));
        }
        make_blocks(&mut f).unwrap();
        make_split_edges(&mut f).unwrap();
        let split: Vec<usize> = f.ds.pave_blocks(e0).iter().map(|pb| pb.edge()).collect();
        assert!(split.iter().any(|&e| e != e0), "split edges created: {split:?}");
        // The face info lists only the source edge.
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        make_pcurves_full(&mut f).unwrap();

        // The source edge got a pcurve.
        assert!(boptools_2d::curve_on_surface(&b.edges[0], &b.faces[0]).is_some());
        // Every split edge got a pcurve on the same face, with endpoints on the
        // face UV quad boundary (u = 0, v within [0, 1]).
        for &sp in &split {
            let Some(shape) = f.ds.shape(sp).cloned() else { continue };
            let e = Edge(shape);
            let pc = boptools_2d::curve_on_surface(&e, &b.faces[0]).expect("split edge pcurve");
            let (a, b) = BRepTool::edge_parameters(&e);
            let q0 = pc.d0(a);
            let q1 = pc.d0(b);
            assert!(q0.x().abs() < 1e-6, "split start {q0:?}");
            assert!(q1.x().abs() < 1e-6, "split end {q1:?}");
            assert!(q0.y() >= -1e-6 && q0.y() <= 1.0 + 1e-6, "split start y {q0:?}");
            assert!(q1.y() >= -1e-6 && q1.y() <= 1.0 + 1e-6, "split end y {q1:?}");
        }
    }

    #[test]
    fn make_pcurves_full_updates_vertex_2d_tolerance() {
        let b = unit_box();
        let mut f = StubFiller::default();
        f.ds.append(b.solid.0.clone()).unwrap();
        let e0 = f.ds.index(&b.edges[0].0).unwrap();
        let f0 = f.ds.index(&b.faces[0].0).unwrap();
        // Give the face a non-trivial tolerance.
        let reg = GeometryRegistry::global();
        let mut fg = reg.face_geom(&b.faces[0].0).expect("face geom");
        fg.tolerance = 0.05;
        reg.set_face(&b.faces[0].0, fg);
        {
            let mut fi = crate::bopds::BopdsFaceInfo::new(f0);
            fi.add_pave(e0, 0.0, 1.0);
            f.ds.change_face_info_pool().push(fi);
        }
        make_pcurves_full(&mut f).unwrap();

        // The boundary vertices of edge 0 carry a 2D tolerance ≥ face tolerance.
        let verts = edge_vertex_shapes(&b.edges[0]);
        for v in &verts {
            let tol = BRepTool::vertex_tolerance(&Vertex(v.clone()));
            assert!(tol >= 0.05 - 1e-12, "vertex tolerance {tol}");
        }
    }

    #[test]
    fn align_pcurve_endpoints_repairs_off_endpoints() {
        let b = unit_box();
        let edge = b.edges[0].clone();
        let face = b.faces[0].clone();
        let pc = pcurve_full::make_pcurve_full(&edge, &face).unwrap();
        // Shift the pcurve in UV so both endpoints leave the projected
        // positions by more than the alignment tolerance.
        let mut tr = occt_core::gp::GpTrsf2d::identity();
        tr.set_translation_vec(&occt_core::gp::GpVec2d::new(0.25, 0.0));
        let shifted: Arc<dyn Curve2d> = Arc::from(pc.transformed(&tr));
        let aligned = align_pcurve_endpoints(&shifted, &edge, &face, 1e-9).expect("aligned");
        // The endpoints snap back onto the projections of the edge's 3D ends.
        let q0 = aligned.d0(0.0);
        let q1 = aligned.d0(1.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {q0:?}");
        assert!((q1.x() - 0.0).abs() < 1e-6 && (q1.y() - 1.0).abs() < 1e-6, "end {q1:?}");
        // The interior is resampled from the shifted curve (still on the face).
        let qm = aligned.d0(0.5);
        assert!((qm.x() - 0.25).abs() < 1e-6, "mid {qm:?}");
    }

    #[test]
    fn align_pcurve_endpoints_leaves_trimmed_subrange_alone() {
        let b = unit_box();
        let edge = b.edges[0].clone();
        let face = b.faces[0].clone();
        let pc = pcurve_full::make_pcurve_full(&edge, &face).unwrap();
        // Simulate a trim that clipped the pcurve to a sub-range: the endpoints
        // are no longer at the edge boundary parameters, so alignment is skipped.
        let n = 17;
        let pts: Vec<GpPnt2d> = (0..n)
            .map(|i| pc.d0(0.25 + 0.5 * i as f64 / (n - 1) as f64))
            .collect();
        let bs = bspline_from_pts_2d(&pts, 0.25, 0.75, 1).unwrap();
        let trimmed: Arc<dyn Curve2d> = Arc::new(bs);
        let out = align_pcurve_endpoints(&trimmed, &edge, &face, 1e-9).expect("aligned");
        assert!(Arc::ptr_eq(&trimmed, &out), "trimmed pcurve returned unchanged");
    }
}

mod p01;
mod p02;
pub use p01::*;
pub use p02::*;
