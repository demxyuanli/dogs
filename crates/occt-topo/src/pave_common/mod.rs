//! Shrunk-range data and pave-block updates — Phase 19 wave C2b.
//!
//! Port of the `BOPAlgo_PaveFiller` shrunk-data/update/self-interference
//! methods:
//!
//! | OCCT method                        | Rust function                              |
//! |------------------------------------|--------------------------------------------|
//! | `FillShrunkData(PaveBlock)`        | [`fill_shrunk_data_for_block`]             |
//! | `FillShrunkData(TopAbs, TopAbs)`   | [`fill_shrunk_data`]                       |
//! | `AnalyzeShrunkData`                | [`analyze_shrunk_data`]                    |
//! | `SetNonDestructive`                | [`set_non_destructive`]                    |
//! | `UpdatePaveBlocksWithSDVertices`   | [`update_pave_blocks_with_sd_vertices`]    |
//! | `UpdateEdgeTolerance`              | [`update_edge_tolerance`]                  |
//! | `UpdateVertex`                     | [`update_vertex_sd`]                       |
//! | `UpdateInterfsWithSDVertices`      | [`update_interfs_with_sd_vertices`]        |
//! | `UpdateCommonBlocksWithSDVertices` | [`update_common_blocks_with_sd_vertices`]  |
//! | `CheckSelfInterference`            | [`check_self_interference`]                |
//! | `RemoveMicroEdges`                 | [`remove_micro_edges`]                     |
//!
//! Source: `BOPAlgo_PaveFiller_3.cxx` (FillShrunkData/AnalyzeShrunkData),
//! `BOPAlgo_PaveFiller_9.cxx` (the edge-wide FillShrunkData),
//! `BOPAlgo_PaveFiller_10.cxx` (SetNonDestructive/Update*),
//! `BOPAlgo_PaveFiller_11.cxx` (CheckSelfInterference).
//!
//! Every function operates directly on the concrete
//! [`crate::pave_filler::PaveFiller`] (landed by the sibling Phase-19 agent),
//! accessing the data structure through `ds()`/`ds_mut()` and the alert report
//! through `add_error()`/`add_warning()`.
//!
//! ## Port simplifications
//!
//! - [`update_interfs_with_sd_vertices`] rewrites `index_new` on the typed
//!   `InterfVV/VE/VF/EE/EF` arrays through [`BopdsDS::update_interfs_with_sd_vertices`].
//! - The Rust `BopdsCommonBlock` records only shared edge indices and ranges
//!   (no bound vertices), so [`update_common_blocks_with_sd_vertices`] reduces
//!   to the pave-block SD pass.
//! - [`check_self_interference`] uses the [`crate::bopds::BopdsIteratorSI`]
//!   self-intersection candidates filtered by a bbox-overlap measure rather
//!   than the OCCT face-info connection map (the Rust `BopdsFaceInfo` does not
//!   retain the IN/Section vertex sets). It detects overlapping solids, faces
//!   and crossing edges within one argument; boundary-touching sub-shapes are
//!   excluded.
mod prelude {

pub(crate) use std::collections::HashSet;

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::algo_tools::AlgoTools;
pub(crate) use crate::bopds::{BopdsDS, BopdsInterf, BopdsPaveBlock, BopdsShapeInfo};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::inttools_range::ShrunkRange;
pub(crate) use crate::pave_filler::PaveFiller;
pub(crate) use crate::shape::{Edge, Vertex};
pub(crate) use crate::tgeometry::GeometryRegistry;

}

use prelude::*;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpDir, GpLin, GpPnt};
    use occt_geom::GeomLine;

    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;

    /// Release the geometry-side-table entries of a shape tree so tests do not
    /// leave stale registry entries keyed by a freed TShape address.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("unit direction")
    }

    /// A filler whose DS holds a single box with a default pave block on every
    /// source edge (the input `fill_shrunk_data` expects).
    ///
    /// Uses `init` (not a bare `append`) so the per-argument ranges and the
    /// source-shape count are set up, which the `is_new_shape` checks in the
    /// update helpers rely on.
    fn filler_with_box_edges() -> (PaveFiller, crate::brep_extrema::test_box::UnitBox) {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        let n = pf.ds().nb_source_shapes();
        for i in 0..n {
            if pf.ds().shape_info(i).map(|s| s.kind) == Some(ShapeType::Edge) {
                pf.ds_mut().init_pave_blocks_for_edge(i);
            }
        }
        (pf, b)
    }

    /// A filler whose DS holds a single box initialized through `init` (source
    /// shapes, ranges and ranks set up).
    fn filler_with_box() -> (PaveFiller, crate::brep_extrema::test_box::UnitBox) {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        (pf, b)
    }

    // -----------------------------------------------------------------------
    // fill_shrunk_data / analyze_shrunk_data
    // -----------------------------------------------------------------------

    #[test]
    fn fill_shrunk_data_sets_valid_ranges_for_box_edges() {
        let (mut pf, b) = filler_with_box_edges();
        fill_shrunk_data(&mut pf).unwrap();
        assert!(!pf.has_errors());

        let e0 = pf.ds().index(&b.edges[0].0).expect("edge indexed");
        let blocks = pf.ds().pave_blocks(e0);
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].has_shrunk_data(), "pave block must carry shrunk data");
        let (ts1, ts2, splittable) = blocks[0].shrunk_data();
        assert!(ts1 > 0.0, "shrunk first {ts1} must be > 0 (edge range is [0, 1])");
        assert!(ts2 < 1.0, "shrunk last {ts2} must be < 1 (edge range is [0, 1])");
        assert!(ts2 > ts1, "shrunk range must be non-empty: {ts1}..{ts2}");
        assert!(splittable, "a unit-length edge must be splittable");

        // Every box edge got a valid shrunk range.
        for edge in &b.edges {
            let e = pf.ds().index(&edge.0).expect("edge indexed");
            let blocks = pf.ds().pave_blocks(e);
            assert_eq!(blocks.len(), 1);
            let (s1, s2, _) = blocks[0].shrunk_data();
            assert!(s1 > 0.0 && s2 < 1.0 && s2 > s1, "edge {e}: shrunk {s1}..{s2}");
        }
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fill_shrunk_data_is_idempotent() {
        let (mut pf, b) = filler_with_box_edges();
        fill_shrunk_data(&mut pf).unwrap();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let before = pf.ds().pave_blocks(e0)[0].shrunk_data();
        let n_warnings = pf.warnings().len();

        // A second pass must not recompute / re-warn for blocks that already
        // carry shrunk data.
        fill_shrunk_data(&mut pf).unwrap();
        let after = pf.ds().pave_blocks(e0)[0].shrunk_data();
        assert_eq!(before, after);
        assert_eq!(pf.warnings().len(), n_warnings);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fill_shrunk_data_on_edge_without_curve_warns() {
        let mut pf = PaveFiller::new();
        let e = Edge::new();
        let mut pb = BopdsPaveBlock::new();
        pb.set_range(0.0, 1.0);
        fill_shrunk_data_for_block(&mut pf, &e, 1e-7, &mut pb);
        assert!(
            !pb.has_shrunk_data(),
            "OCCT HasShrunkData is false when AnalyzeShrunkData stores a void box"
        );
        assert!(!pb.is_splittable());
        assert!(
            pf.warnings().iter().any(|w| w.contains("bad positioning")),
            "expected a bad-positioning warning, got: {:?}",
            pf.warnings()
        );
    }

    // -----------------------------------------------------------------------
    // set_non_destructive
    // -----------------------------------------------------------------------

    #[test]
    fn set_non_destructive_toggles_flag() {
        let mut pf = PaveFiller::new();
        assert!(!pf.non_destructive());
        set_non_destructive(&mut pf, true);
        assert!(pf.non_destructive());
        set_non_destructive(&mut pf, false);
        assert!(!pf.non_destructive());
    }

    // -----------------------------------------------------------------------
    // update_pave_blocks_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_pave_blocks_with_sd_vertices_rewrites_bounds() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v1 = pf.ds().index(&b.vertices[1].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();

        pf.ds_mut().init_pave_blocks_for_edge(e0);
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v0, v1));

        // v0 is merged into v5 (same-domain chain).
        pf.ds_mut().add_shape_sd(v0, v5);
        update_pave_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v5, v1));

        // A block with no SD vertices is left untouched.
        let e1 = pf.ds().index(&b.edges[1].0).unwrap();
        pf.ds_mut().init_pave_blocks_for_edge(e1);
        let before = pf.ds().pave_blocks(e1)[0].indices();
        update_pave_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e1)[0].indices(), before);
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_edge_tolerance / update_vertex_sd
    // -----------------------------------------------------------------------

    #[test]
    fn update_edge_tolerance_raises_edge_and_vertex_tolerances() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();

        update_edge_tolerance(&mut pf, e0, 0.01).unwrap();
        let edge = Edge(pf.ds().shape(e0).unwrap().clone());
        assert!(
            (BRepTool::edge_tolerance(&edge) - 0.01).abs() < 1e-12,
            "edge tolerance raised to 0.01, got {}",
            BRepTool::edge_tolerance(&edge)
        );
        let vtx = Vertex(pf.ds().shape(v0).unwrap().clone());
        assert!(
            (BRepTool::vertex_tolerance(&vtx) - 0.01).abs() < 1e-12,
            "bound vertex tolerance raised to 0.01, got {}",
            BRepTool::vertex_tolerance(&vtx)
        );
        clear_tree(&b.solid.0);
    }

    #[test]
    fn update_edge_tolerance_non_destructive_skips_original_edge() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        set_non_destructive(&mut pf, true);

        update_edge_tolerance(&mut pf, e0, 0.01).unwrap();
        // The original edge is not modified in safe-input mode.
        let edge = Edge(pf.ds().shape(e0).unwrap().clone());
        assert!(
            (BRepTool::edge_tolerance(&edge) - 0.0).abs() < 1e-12,
            "original edge tolerance must stay 0.0 in non-destructive mode, got {}",
            BRepTool::edge_tolerance(&edge)
        );
        clear_tree(&b.solid.0);
    }

    #[test]
    fn update_vertex_sd_creates_replacement_in_non_destructive_mode() {
        let (mut pf, b) = filler_with_box();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let n_before = pf.ds().nb_shapes();
        set_non_destructive(&mut pf, true);

        let n_v_new = update_vertex_sd(&mut pf, v0, 0.02).unwrap();
        assert_ne!(n_v_new, v0, "a replacement vertex must be allocated");
        assert_eq!(pf.ds().nb_shapes(), n_before + 1);
        assert_eq!(pf.ds().has_shape_sd(v0), Some(n_v_new));
        let vtx = Vertex(pf.ds().shape(n_v_new).unwrap().clone());
        assert!(
            (BRepTool::vertex_tolerance(&vtx) - 0.02).abs() < 1e-12,
            "replacement tolerance raised to 0.02, got {}",
            BRepTool::vertex_tolerance(&vtx)
        );
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_interfs_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_intfs_with_sd_vertices_rewrites_index_new() {
        let (mut pf, b) = filler_with_box();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();
        pf.ds_mut().add_shape_sd(v0, v5);

        let mut interfs = [
            BopdsInterf::new(1, 2),
            BopdsInterf::new(3, 4),
        ];
        interfs[0].set_index_new(v0);
        update_intfs_with_sd_vertices(pf.ds(), &mut interfs);
        assert_eq!(interfs[0].get_index_new(), Some(v5), "index_new rewritten through the SD map");
        assert_eq!(interfs[1].get_index_new(), None, "an unset index_new stays unset");

        update_interfs_with_sd_vertices(&mut pf).unwrap();
        assert!(!pf.has_errors());
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // update_common_blocks_with_sd_vertices
    // -----------------------------------------------------------------------

    #[test]
    fn update_common_blocks_with_sd_vertices_updates_pave_bounds() {
        let (mut pf, b) = filler_with_box();
        let e0 = pf.ds().index(&b.edges[0].0).unwrap();
        let v0 = pf.ds().index(&b.vertices[0].0).unwrap();
        let v5 = pf.ds().index(&b.vertices[5].0).unwrap();

        pf.ds_mut().init_pave_blocks_for_edge(e0);
        pf.ds_mut().add_shape_sd(v0, v5);
        update_common_blocks_with_sd_vertices(&mut pf).unwrap();
        assert_eq!(pf.ds().pave_blocks(e0)[0].indices(), (v5, pf.ds().index(&b.vertices[1].0).unwrap()));
        clear_tree(&b.solid.0);
    }

    // -----------------------------------------------------------------------
    // check_self_interference
    // -----------------------------------------------------------------------

    #[test]
    fn check_self_interference_clean_box_is_false() {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "a single closed box is not self-interfering");
        assert!(!pf.has_warnings(), "no warnings expected, got {:?}", pf.warnings());
        clear_tree(&b.solid.0);
    }

    #[test]
    fn check_self_interference_overlapping_boxes_is_true() {
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        let b = TopoBuilder::new();
        // Two identically-positioned boxes as one argument: their solids
        // overlap with positive volume.
        let comp = b.make_compound_of(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.set_arguments(&[comp.into()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "SI mode (one argument) skips CheckSelfInterference");
        assert!(!pf.has_warnings());
        clear_tree(&a.solid.0);
        clear_tree(&c.solid.0);
    }

    #[test]
    fn check_self_interference_crossing_edges_is_true() {
        let mut pf = PaveFiller::new();
        let b = TopoBuilder::new();
        // Edge (0,0,0)->(1,1,0) and edge (0,1,0)->(1,0,0) cross at (0.5, 0.5).
        let e1 = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), dir(1.0, 1.0, 0.0)))),
            0.0,
            2.0f64.sqrt(),
        );
        let e2 = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(GpPnt::new(0.0, 1.0, 0.0), dir(1.0, -1.0, 0.0)))),
            0.0,
            2.0f64.sqrt(),
        );
        let comp = b.make_compound_of(&[e1.0.clone(), e2.0.clone()]);
        pf.set_arguments(&[comp.into()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "SI mode (one argument) skips CheckSelfInterference");
        assert!(!pf.has_warnings());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn check_self_interference_separate_arguments_is_false() {
        // Two overlapping boxes passed as TWO arguments are not acquired
        // self-interference — the overlap is the intended operation.
        let mut pf = PaveFiller::new();
        let a = unit_box();
        let c = unit_box();
        pf.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        pf.init().unwrap();

        let res = check_self_interference(&mut pf).unwrap();
        assert!(!res, "inter-argument overlap is not self-interference");
        assert!(!pf.has_warnings());
        clear_tree(&a.solid.0);
        clear_tree(&c.solid.0);
    }

    // -----------------------------------------------------------------------
    // remove_micro_edges
    // -----------------------------------------------------------------------

    #[test]
    fn remove_micro_edges_removes_degenerate_blocks() {
        let mut pf = PaveFiller::new();
        let b = unit_box();
        pf.set_arguments(&[b.solid.0.clone()]);
        pf.init().unwrap();
        let e = pf.ds().index(&b.edges[0].0).expect("edge indexed");
        // A curve-less edge appended to the DS: any block on it fails the
        // shrunk-range computation, so it is a micro edge.
        let n_empty = pf.ds_mut().append(Edge::new().0).unwrap();
        let subs = pf.ds().shape_info(e).unwrap().sub_shapes().to_vec();
        let (v1, v2) = (subs[0], subs[1]);
        {
            let blocks = pf.ds_mut().change_pave_blocks_mut(e);
            // Block 1: normal (distinct bounds on a real edge).
            let mut pb1 = BopdsPaveBlock::new();
            pb1.set_edge(e);
            pb1.set_original_edge(e);
            pb1.set_range(0.0, 0.4);
            pb1.set_indices(v1, v2);
            // Block 2: degenerate (coincident bounds) on a curve-less edge.
            let mut pb2 = BopdsPaveBlock::new();
            pb2.set_edge(n_empty);
            pb2.set_original_edge(e);
            pb2.set_range(0.4, 0.6);
            pb2.set_indices(n_empty, n_empty);
            blocks.clear();
            blocks.push(pb1);
            blocks.push(pb2);
        }

        remove_micro_edges(&mut pf);

        let blocks = pf.ds().pave_blocks(e);
        assert!(
            !blocks.iter().any(|pb| pb.edge() == n_empty),
            "the degenerate block must be removed, got: {:?}",
            blocks
        );
        assert!(blocks.iter().any(|pb| pb.edge() == e), "the normal block is kept");
        clear_tree(&b.solid.0);
    }
}

mod p01;
pub use p01::*;
