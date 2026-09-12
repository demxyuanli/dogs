use super::*;
use crate::brep_gprop::volume as shape_volume;
use crate::primitives::BRepPrimBox;
use crate::shell_check::shell_invariants;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{faces_of, shapes_of};
use occt_core::gp::GpPnt;

    pub(super) fn box_vol(s: &TopoShape) -> f64 {
        shape_volume(s, 0.02)
    }

    pub(super) fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    /// Faceted cylinder solid (axis +Z, base at z=1, centered at (1,1)) from a
    /// triangle mesh — all faces planar, like the boss tool.
    pub(super) fn faceted_cylinder() -> TopoShape {
        let (r, h, slices) = (0.25, 0.8, 24usize);
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();
        for i in 0..=slices {
            let th = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
            vertices.push(GpPnt::new(1.0 + r * th.cos(), 1.0 + r * th.sin(), 1.0));
            vertices.push(GpPnt::new(1.0 + r * th.cos(), 1.0 + r * th.sin(), 1.0 + h));
        }
        for i in 0..slices {
            let (a, b, c, d) = (2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3);
            triangles.push(occt_core::poly::triangulation::Triangle::new(a, b, c));
            triangles.push(occt_core::poly::triangulation::Triangle::new(b, d, c));
        }
        let (top_idx, bot_idx) = (vertices.len(), vertices.len() + 1);
        vertices.push(GpPnt::new(1.0, 1.0, 1.0 + h));
        vertices.push(GpPnt::new(1.0, 1.0, 1.0));
        for i in 0..slices {
            let (tb, tt) = (2 * i + 1, 2 * i + 3);
            let (bb, bt) = (2 * i, 2 * i + 2);
            triangles.push(occt_core::poly::triangulation::Triangle::new(top_idx, tb, tt));
            triangles.push(occt_core::poly::triangulation::Triangle::new(bot_idx, bt, bb));
        }
        let mesh = crate::mesh::ShapeMesh { vertices, triangles, source_shape: ShapeType::Solid };
        crate::mesh_to_brep::shape_mesh_to_brep(&mesh).solid.expect("closed").0
    }

    /// Wave 3/4 gate: a planar fuse of a box with a protruding faceted cylinder
    /// must yield a closed, genus-0 solid (manifold + Euler characteristic 2).
    #[test]
    fn fuse_box_cylinder_is_closed_solid() {
        let a = BRepPrimBox::make_box(2.0, 2.0, 1.0).solid.0;
        let tool = faceted_cylinder();
        let r = boolean(&a, &tool, BoolOp::Fuse, 1e-4).expect("fuse");
        let n_solids = shapes_of(&r.shape, ShapeType::Solid).len();
        let n_faces = faces_of(&r.shape).len();
        let shape = r.solid.clone().map(|s| s.0).unwrap_or(r.shape.clone());
        let shells = shapes_of(&shape, ShapeType::Shell);
        assert!(
            !shells.is_empty(),
            "fuse produced no shell; nsolids={n_solids} nfaces={n_faces} type={:?} warnings {:?}",
            r.shape.shape_type(),
            r.warnings
        );
        for s in &shells {
            let inv = shell_invariants(&Shell(s.clone()));
            assert!(
                inv.is_valid_solid(),
                "invariants {inv:?} nsolids={n_solids} nfaces={n_faces} nshells={} type={:?} warnings {:?}",
                shells.len(),
                r.shape.shape_type(),
                r.warnings
            );
        }
    }

    pub(super) fn overlapping_boxes() -> (Solid, Solid) {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        (a.solid, b.solid)
    }

    #[test]
    fn fuse_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(
            r.solid.is_some(),
            "fuse produces a solid; type={:?} nsolids={} nshells={} nfaces={} warnings {:?}",
            r.shape.shape_type(),
            shapes_of(&r.shape, ShapeType::Solid).len(),
            r.shells.len(),
            faces_of(&r.shape).len(),
            r.warnings
        );
        assert!(shell_is_closed(&r.shells[0]), "fuse shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "fuse volume {v} (expected 1.5)");
        assert!(faces_of(&r.shape).len() > 12, "faces {} (expected > 12)", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn cut_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(r.solid.is_some(), "cut produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "cut shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 0.5).abs() < 0.05, "cut volume {v} (expected 0.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn common_overlapping_boxes() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!(r.solid.is_some(), "common produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "common shell is closed");
        let v = box_vol(&r.shape);
        assert!((v - 0.5).abs() < 0.05, "common volume {v} (expected 0.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn disjoint_boxes() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(2.0, 2.0, 2.0), &GpPnt::new(3.0, 3.0, 3.0));

        let f = boolean(&a.solid.0, &b.solid.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!((box_vol(&f.shape) - 2.0).abs() < 0.05, "disjoint fuse volume {}", box_vol(&f.shape));

        let c = boolean(&a.solid.0, &b.solid.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!((box_vol(&c.shape) - 1.0).abs() < 0.05, "disjoint cut volume {}", box_vol(&c.shape));

        let m = boolean(&a.solid.0, &b.solid.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!((box_vol(&m.shape)).abs() < 1e-9, "disjoint common volume {}", box_vol(&m.shape));

        clear_tree(&f.shape);
        clear_tree(&c.shape);
        clear_tree(&m.shape);
        clear_tree(&a.solid.0);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn cube_minus_internal_box() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.2, 0.2, 0.2), &GpPnt::new(0.8, 0.8, 0.8));
        let r = boolean(&a.solid.0, &b.solid.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(r.solid.is_some(), "hollow cube produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "hollow cube shell is closed");
        let v = box_vol(&r.shape);
        let expected = 1.0 - 0.6 * 0.6 * 0.6;
        assert!((v - expected).abs() < 0.05, "hollow cube volume {v} (expected {expected})");
        clear_tree(&r.shape);
        clear_tree(&a.solid.0);
        clear_tree(&b.solid.0);
    }

    #[test]
    fn fuse_matches_voxel_union() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let exact = box_vol(&r.shape);
        let voxel = crate::solid_union::union_volume(&a.0, &b.0, 32);
        assert!(voxel > 1e-9, "voxel volume should be non-zero");
        let rel = (exact - voxel).abs() / voxel;
        assert!(rel < 0.15, "fuse {exact} vs voxel {voxel} (rel {rel:.3})");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // ------------------------------------------------------------------
    // Full-topology tests
    // ------------------------------------------------------------------

    /// Axis-aligned box at `[lo, hi]` via `BRepPrim_GWedge` (`make_box_corner`).
    /// The previous handmade 8-vertex builder is not used: `make_box_corner`
    /// already shares TShape with `make_box`.
    pub(super) fn test_box_at(lo: &GpPnt, hi: &GpPnt) -> Solid {
        BRepPrimBox::make_box_corner(lo, hi).solid
    }

    pub(super) fn disjoint_box_shapes() -> Vec<TopoShape> {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let b3 = test_box_at(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(1.0, 3.0, 1.0));
        vec![b1.0, b2.0, b3.0]
    }

    #[test]
    fn fuse_three_boxes_compound_or_solid() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        assert!(r.shape.is_compound(), "three disjoint boxes fuse to a compound");
        let subs = decompose_compound(&r.shape);
        assert_eq!(subs.len(), 3, "compound has 3 sub-shapes, got {}", subs.len());
        for s in &subs {
            assert!(s.is_solid(), "each sub-shape is a solid");
        }
        let v = box_vol(&r.shape);
        assert!((v - 3.0).abs() < 0.1, "disjoint multi fuse volume {v} (expected 3.0)");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

