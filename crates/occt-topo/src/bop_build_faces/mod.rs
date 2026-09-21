//! Face-image reconstruction — Phase 20 wave C2b-2a.
//!
//! Port of `BOPAlgo_Builder::FillImagesFaces` (with its two stages
//! `BuildSplitFaces` and `FillSameDomainFaces`) from
//! `BOPAlgo_Builder_2.cxx`:
//!
//! | OCCT method               | Rust function                         |
//! |---------------------------|---------------------------------------|
//! | `FillImagesFaces`         | [`fill_images_faces`]                 |
//! | `BuildSplitFaces`         | [`build_split_faces`]                 |
//! | `FillSameDomainFaces`     | [`fill_same_domain_faces`]            |
//!
//! The three entry points are generic over a [`BopBuilderLike`] host — the
//! minimal contract `crate::bop_builder2::BopBuilder` (General Fuse builder,
//! Phase 20) satisfies. The host exposes the BOPDS, the images history, the
//! same-domain map and the option flags the OCCT builder stores as members
//! (`myDS`, `myImages`, `myShapesSD`, `myPaveFiller`, …). Keeping the logic
//! decoupled from the concrete type follows the [`crate::pave_blocks::PaveFillerLike`]
//! / [`crate::bop_build_common::BopBuildOps`] pattern: the module is verifiable
//! in isolation with a local test stub, and the builder main class adapts by
//! implementing the trait (the only addition it needs is a mutable accessor to
//! its same-domain map).
//!
//! Algorithm overview (matching the C++ flow):
//!
//! 1. [`build_split_faces`] — for every source face carrying on-face
//!    (IN/section) edges in the BOPDS face-info pool, collect the face's
//!    bounding edges (substituting split images) plus the on-face edges, close
//!    the set into wires with [`crate::wire_splitter::WireSplitter`], then
//!    rebuild one closed face per wire on the original surface with
//!    [`crate::builder_face::FaceBuilder`]. The resulting faces are recorded as
//!    images of the source face.
//! 2. [`fill_same_domain_faces`] — group the result faces by (boundary edge
//!    signature, geometric surface); each group of coincident faces collapses
//!    to one representative recorded in the same-domain map
//!    ([`BopBuilderLike::bind_shapes_sd`]).
//!
//! Dependencies: [`crate::builder_face`], [`crate::wire_splitter`],
//! [`crate::bopds`], [`crate::pave_filler`], [`crate::algo_tools`],
//! [`crate::bop_builder2`].
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_geom::Surface;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::algo_tools::AlgoTools;
pub(crate) use crate::bop_hist::BopHistory;
pub(crate) use crate::bopds::BopdsDS;
pub(crate) use crate::boptools_2d;
pub(crate) use crate::brep_surface::surface_closest_params;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder_area::AreaBuilder;
pub(crate) use crate::builder_face::{make_face_from_wire, FaceBuilder};
pub(crate) use crate::fclass2d::{FClass2d, FaceState};
pub(crate) use crate::int_tools_full::IntToolsContext;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, vertex_position};
pub(crate) use crate::wire_splitter::{WireEdgeSet, WireSplitter};

}

use prelude::*;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bopds::BopdsFaceInfo;
    use crate::brep_surface::face_plane;
    use crate::builder::TopoBuilder;
    use crate::builder_face::loop_signed_area;
    use crate::brep_extrema::test_box::unit_box;
    use crate::topo_tools_full::wires_of_face;
    use occt_core::gp::GpPnt;

    /// A minimal `BopBuilderLike` host for the isolated tests.
    #[derive(Debug)]
    struct StubBopBuilder {
        ds: BopdsDS,
        history: BopHistory,
        errors: Vec<String>,
        warnings: Vec<String>,
        shapes_sd: Vec<(TopoShape, TopoShape)>,
        origins: HashMap<usize, Vec<TopoShape>>,
        fuzzy: f64,
    }

    impl StubBopBuilder {
        fn new() -> Self {
            Self {
                ds: BopdsDS::new(),
                history: BopHistory::new(),
                errors: Vec::new(),
                warnings: Vec::new(),
                shapes_sd: Vec::new(),
                origins: HashMap::new(),
                fuzzy: 1e-7,
            }
        }
    }

    impl BopBuilderLike for StubBopBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn ds_mut(&mut self) -> &mut BopdsDS {
            &mut self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn has_errors(&self) -> bool {
            !self.errors.is_empty()
        }
        fn add_error(&mut self, msg: String) {
            self.errors.push(msg);
        }
        fn add_warning(&mut self, msg: String) {
            self.warnings.push(msg);
        }
        fn non_destructive(&self) -> bool {
            false
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn bind_shapes_sd(&mut self, shape: TopoShape, sd: TopoShape) {
            if let Some(e) = self.shapes_sd.iter_mut().find(|(s, _)| s.same_tshape(&shape)) {
                e.1 = sd;
            } else {
                self.shapes_sd.push((shape, sd));
            }
        }
        fn seek_shapes_sd(&self, shape: &TopoShape) -> Option<TopoShape> {
            self.shapes_sd
                .iter()
                .find(|(s, _)| s.same_tshape(shape))
                .map(|(_, sd)| sd.clone())
        }
        fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
            &mut self.origins
        }
    }

    /// True polygon area of a planar face: the sum of |signed loop areas| of
    /// its wires projected on the supporting plane. Unlike `face_area_exact`
    /// (which triangulates the whole UV rectangle), this gives the trimmed
    /// area of an unbounded-plane face.
    fn polygon_area_on_plane(face: &Face) -> f64 {
        let Some(pln) = face_plane(face) else { return f64::NAN };
        wires_of_face(face)
            .iter()
            .map(|w| loop_signed_area(&edges_of_wire(w), &pln).abs())
            .sum()
    }

    /// Record the diagonal `(0,0,0)-(1,1,0)` as an on-face edge of the bottom
    /// face (`z = 0`) of `ub`, in the stub's BOPDS.
    fn stub_with_split_bottom_face(ub: &crate::brep_extrema::test_box::UnitBox) -> (StubBopBuilder, usize) {
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        let bottom = &ub.faces[0];
        let face_idx = stub.ds().index(&bottom.0).expect("bottom face indexed");

        let b = TopoBuilder::new();
        let diag = b.make_edge_segment_with_vertices(
            &GpPnt::new(0.0, 0.0, 0.0),
            &GpPnt::new(1.0, 1.0, 0.0),
            &ub.vertices[0],
            &ub.vertices[2],
        );
        let diag_idx = stub.ds_mut().append(diag.0.clone()).expect("append diagonal edge");

        let pool = stub.ds_mut().change_face_info_pool();
        let fi = match pool.iter_mut().find(|fi| fi.face_index == face_idx) {
            Some(fi) => fi,
            None => {
                pool.push(BopdsFaceInfo::new(face_idx));
                pool.last_mut().expect("just pushed")
            }
        };
        fi.add_pave(diag_idx, 0.0, 1.0);
        (stub, face_idx)
    }

    #[test]
    fn box_face_split_by_diagonal_builds_two_closed_faces() {
        let ub = unit_box();
        let (mut stub, face_idx) = stub_with_split_bottom_face(&ub);

        build_split_faces(&mut stub).expect("build split faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        assert!(stub.warnings.is_empty(), "warnings: {:?}", stub.warnings);

        // The bottom face must have exactly two images (the two triangles).
        let orig = stub.ds().shape(face_idx).cloned().expect("face shape");
        let imgs = stub.history().image(&orig).expect("images recorded");
        assert_eq!(imgs.len(), 2, "two split faces, got {}", imgs.len());

        let a1 = polygon_area_on_plane(&Face(imgs[0].clone()));
        let a2 = polygon_area_on_plane(&Face(imgs[1].clone()));
        assert!(
            (a1 + a2 - 1.0).abs() < 1e-6,
            "split areas must sum to the original face area: {a1} + {a2}"
        );
        assert!(
            (a1 - 0.5).abs() < 1e-6 && (a2 - 0.5).abs() < 1e-6,
            "each diagonal half is a triangle of area 0.5, got {a1}, {a2}"
        );

        // Both split faces are closed (single wire, closed loop).
        for im in imgs {
            let f = Face(im.clone());
            assert_eq!(wires_of_face(&f).len(), 1, "one wire per split face");
        }
    }

    #[test]
    fn fill_same_domain_merges_coincident_keeps_parallel() {
        // OCCT FillSameDomainFaces groups by BOPTools_Set of EDGE TShapes
        // (after intersection common-blocks unify coincident edges). A second
        // face rebuilt on the same wires shares those TShapes; InterfFF +
        // FaceInfo are the other two filters (`_2.cxx:586-680`).
        let ub1 = unit_box();
        let bottom = &ub1.faces[0];
        let wires = wires_of_face(bottom);
        let surf = BRepTool::face_surface_world(bottom).expect("bottom surface");
        let bld = TopoBuilder::new();
        let f2 = bld.make_face(surf, &wires);
        let sh2 = bld.make_shell(&[f2.clone()]);
        let s2 = bld.make_solid(&[sh2]);

        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub1.solid.0.clone(), s2.0.clone()]);

        let i1 = stub.ds().index(&bottom.0).expect("bottom index");
        let i2 = stub.ds().index(&f2.0).expect("copy index");
        stub.ds_mut().ensure_face_info(i1);
        stub.ds_mut().ensure_face_info(i2);
        stub.ds_mut()
            .append_interf_ff(crate::bopds_ff::BopdsInterfFf::new(i1, i2));

        fill_same_domain_faces(&mut stub).expect("fill same domain");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);

        let expected = if i1 < i2 { bottom.0.clone() } else { f2.0.clone() };
        let sd2 = stub.seek_shapes_sd(&f2.0);
        assert!(sd2.is_some(), "coincident faces that share edges must merge");
        assert!(
            sd2.unwrap().same_tshape(&expected),
            "representative is the smaller DS index"
        );
        let sd1 = stub.seek_shapes_sd(&bottom.0);
        assert!(sd1.is_some());
        assert!(sd1.unwrap().same_tshape(&expected));

        // Parallel but distinct faces of one box do not merge: the bottom and
        // the top of the same box are never bound together (same parent solid,
        // `_2.cxx:776-779`, and no shared EDGE set).
        let f_top = &ub1.faces[1];
        let sd_bottom = stub.seek_shapes_sd(&bottom.0);
        let sd_top = stub.seek_shapes_sd(&f_top.0);
        assert!(!sd_bottom.map_or(false, |s| s.same_tshape(&f_top.0)));
        assert!(!sd_top.map_or(false, |s| s.same_tshape(&bottom.0)));

        // Original faces that acquired a same-domain twin got themselves as
        // images (they are unchanged but now have a representative).
        assert!(stub.history().has_image(&bottom.0));
        assert!(stub.history().has_image(&f2.0));
    }

    #[test]
    fn fill_same_domain_single_box_produces_no_merge() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        fill_same_domain_faces(&mut stub).expect("fill same domain");
        assert!(!stub.has_errors());

        // None of the six distinct box faces is coincident with another.
        assert!(stub.shapes_sd.is_empty(), "no same-domain bindings expected");
        assert!(!stub.history().has_any_images(), "no face should become its own image");
    }

    #[test]
    fn fill_images_faces_runs_clean_without_face_info() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        // No face carries face-info: nothing to split, nothing to merge.
        fill_images_faces(&mut stub).expect("fill images faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        assert!(!stub.history().has_any_images(), "no split faces recorded");
        assert!(stub.shapes_sd.is_empty(), "no same-domain bindings");
    }

    #[test]
    fn face_without_face_info_is_skipped() {
        let ub = unit_box();
        let mut stub = StubBopBuilder::new();
        stub.ds_mut().init(&[ub.solid.0.clone()]);

        // Only the top face carries a face-info entry — and it is empty.
        let top = &ub.faces[1];
        let top_idx = stub.ds().index(&top.0).expect("top face indexed");
        stub.ds_mut().change_face_info_pool().push(BopdsFaceInfo::new(top_idx));

        build_split_faces(&mut stub).expect("build split faces");
        assert!(!stub.has_errors(), "errors: {:?}", stub.errors);
        // An empty face-info record means the face was not split.
        assert!(!stub.history().has_any_images(), "no split faces recorded");
    }

    /// A hole wire passed *before* its growth wire must still be attached to
    /// the growth, not emitted as a standalone overlapping face.
    ///
    /// The old greedy grouping scanned `j > i` and took the first enclosing
    /// wire, so a hole at index 0 had no growth to attach to. The
    /// `PerformAreas` port classifies all wires first, then attaches every hole
    /// to its (nearest) containing growth — order-independent.
    #[test]
    fn group_wires_as_areas_attaches_hole_regardless_of_order() {
        use occt_core::gp::{GpAx3, GpPln};
        use occt_geom::GeomPlane;

        let b = TopoBuilder::new();
        // Outer square and an inner diamond, both CCW in the +Z plane.
        let sq = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let dia = [
            GpPnt::new(0.2, 0.5, 0.0),
            GpPnt::new(0.5, 0.2, 0.0),
            GpPnt::new(0.8, 0.5, 0.0),
            GpPnt::new(0.5, 0.8, 0.0),
        ];
        let sq_es: Vec<Edge> = (0..4).map(|i| b.make_edge_segment(&sq[i], &sq[(i + 1) % 4])).collect();
        let dia_es: Vec<Edge> = (0..4).map(|i| b.make_edge_segment(&dia[i], &dia[(i + 1) % 4])).collect();
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));

        // Hole wire first: `[diamond, square]`.
        let groups = group_wires_as_areas(&[dia_es, sq_es], &Some(surf));
        assert_eq!(groups.len(), 1, "exactly one growth group");
        assert_eq!(groups[0].0, 1, "the growth is the square (wire 1)");
        assert_eq!(groups[0].1, vec![0], "the diamond (wire 0) is its hole");
    }

    /// Two disjoint growth pieces (a diagonal split) never become each other's
    /// holes: they share no containing relation, so both come out as growths.
    #[test]
    fn group_wires_as_areas_disjoint_pieces_are_both_growths() {
        use occt_core::gp::{GpAx3, GpPln};
        use occt_geom::GeomPlane;

        let b = TopoBuilder::new();
        // Two disjoint triangles inside a unit square (a diagonal split).
        let t1 = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
        ];
        let t2 = [
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(0.0, 0.0, 0.0),
        ];
        let e1: Vec<Edge> = (0..3).map(|i| b.make_edge_segment(&t1[i], &t1[(i + 1) % 3])).collect();
        let e2: Vec<Edge> = (0..3).map(|i| b.make_edge_segment(&t2[i], &t2[(i + 1) % 3])).collect();
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));

        let groups = group_wires_as_areas(&[e1, e2], &Some(surf));
        assert_eq!(groups.len(), 2, "two disjoint pieces -> two growths");
        assert!(groups.iter().all(|(_, hs)| hs.is_empty()), "no holes");
    }
}

mod builder_like;
pub use builder_like::*;
