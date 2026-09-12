//! Boolean operations on curved-face solids.
//! Source: `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_Cut` / `BRepAlgoAPI_Common` (TKBO),
//! `BOPAlgo_Builder`.
//!
//! This module extends the planar-exact boolean (`crate::bop_builder`) to
//! solids whose faces carry analytic curved surfaces (sphere, cylinder, cone).
//! The strategy:
//!
//! 1. when every face of both solids is planar, delegate to
//!    `crate::bop_builder::boolean` (`BOPAlgo_BOP`);
//! 2. otherwise classify each face against the other solid (Inside / Outside /
//!    On) by sampling surface points and testing point-in-solid;
//! 3. keep whole faces whose classification is uniform; for faces that cross
//!    the boundary, triangulate the retained region. Spherical faces crossing
//!    another sphere are trimmed to the analytic spherical cap; the two caps
//!    of an intersecting pair share a discretized intersection circle (and its
//!    in-plane basis) so the assembled mesh is closed and its
//!    divergence-theorem volume is exact;
//! 4. weld the retained triangles into a closed mesh, compute its volume, and
//!    rebuild a BRep solid (via `mesh_to_brep`).
mod prelude {

pub(crate) use std::f64::consts::{FRAC_PI_2, PI};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::bnd::BndBox;
pub(crate) use occt_core::geom::polygon_boolean::{point_in_polygon2d, polygon_boolean, signed_area2d, PolygonBoolOp};
pub(crate) use occt_core::geom::triangulate::triangulate_polygon;
pub(crate) use occt_core::gp::{GpAx1, GpDir, GpPnt, GpPnt2d, GpVec};
pub(crate) use occt_geom::{Curve, Surface};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::bop_builder::{BoolOp, BooleanResult};
pub(crate) use crate::brep_surface::{classify_surface, SurfaceKind};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
pub(crate) use crate::topo_tools_full::{edges_of_wire, faces_of, shapes_of, wires_of_face};

}

use prelude::*;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpAx3, GpDir};
    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
    use crate::shape_mesh::shape_volume;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;

    fn unit_sphere_solid() -> Solid {
        BRepPrimSphere::make_sphere(1.0).solid
    }

    fn sphere_at(center: GpPnt) -> Solid {
        let ax3 = GpAx3::new(center, GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bld = TopoBuilder::new();
        let face = bld.make_face(Arc::new(occt_geom::GeomSphere::new(
            occt_core::gp::GpSphere::new(ax3, 1.0).unwrap(),
        )), &[]);
        let shell = bld.make_shell(&[face]);
        bld.make_solid(&[shell])
    }

    fn centered_box(half: f64) -> Solid {
        BRepPrimBox::make_box_corner(&GpPnt::new(-half, -half, -half), &GpPnt::new(half, half, half)).solid
    }

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    #[test]
    fn box_fuse_box_delegates() {
        // Both planar → delegated to bop_builder (exact polygon boolean).
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.solid.is_some(), "fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "fuse shell is closed");
        assert!(faces_of(&r.shape).len() > 12, "faces {}", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_inside_box_common_volume() {
        // Unit sphere at origin inside box [-1.5,1.5]³: Common = the sphere.
        let sphere = unit_sphere_solid();
        let box_s = centered_box(1.5);
        let r = curved_boolean(&sphere.0, &box_s.0, BoolOp::Common, 1e-6).expect("common ok");
        assert!(r.solid.is_some(), "common produces a solid");
        // The result keeps the whole (curved) sphere face, so shape_volume is
        // reliable for this curved solid. The adaptive mesh's deflection
        // calibration is looser than OCCT's (def 0.05 → ~3.5% sphere error), so
        // a finer deflection is needed for the <0.05 (1.2%) volume gate.
        let v = shape_volume(&r.shape, 0.01);
        let expected = 4.0 / 3.0 * PI;
        assert!((v - expected).abs() < 0.05, "common volume {v} (expected {expected})");
        clear_tree(&r.shape);
        clear_tree(&sphere.0);
        clear_tree(&box_s.0);
    }

    #[test]
    fn sphere_sphere_fuse_volume() {
        // Two unit spheres 1.5 apart: Fuse volume ≈ 8.018.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse volume");
        assert!((7.5..8.5).contains(&v), "fuse volume {v} (expected ≈ 8.018)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_sphere_common_volume() {
        // Two unit spheres 1.5 apart: Common = lens, volume ≈ 0.36.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Common, 1e-6).expect("common volume");
        assert!((0.2..0.6).contains(&v), "common volume {v} (expected ≈ 0.36)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn sphere_sphere_cut_volume() {
        // A − B: sphere 1 minus the lens → volume ≈ 4.188 − 0.36 = 3.83.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let v = curved_boolean_volume(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut volume");
        assert!((3.0..4.2).contains(&v), "cut volume {v} (expected ≈ 3.83)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn disjoint_fuse_returns_compound() {
        // Sphere at origin, sphere far away → Fuse returns a compound.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(5.0, 0.0, 0.0));
        let r = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.shape.is_compound(), "disjoint fuse is a compound");
        assert!(r.solid.is_none(), "disjoint fuse has no single solid");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn classify_sphere_inside_box() {
        let sphere = unit_sphere_solid();
        let box_s = centered_box(1.5);
        let sf = faces_of(&sphere.0);
        assert_eq!(classify_face(&sf[0], &box_s.0, 1e-6), FaceRegion::Inside);
        let bf = faces_of(&box_s.0);
        for f in bf.iter().take(3) {
            assert_eq!(classify_face(f, &sphere.0, 1e-6), FaceRegion::Outside);
        }
        clear_tree(&sphere.0);
        clear_tree(&box_s.0);
    }

    #[test]
    fn ray_hits_box() {
        let box_s = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let p = GpPnt::new(-1.0, 0.5, 0.5);
        let hits = ray_hits_solid(p, &box_s.0, GpVec::new(1.0, 0.0, 0.0), 1e-9);
        assert_eq!(hits, 2, "ray through unit box");
        clear_tree(&box_s.0);
    }

    #[test]
    fn mesh_volume_full_sphere() {
        let sphere = unit_sphere_solid();
        let sf = faces_of(&sphere.0);
        let m = mesh_whole_face(&sf[0]);
        let v = mesh_volume(&m);
        let expected = 4.0 / 3.0 * PI;
        // Chord triangulation at 24×24 is ~1.6% low; keep a relaxed band.
        assert!((v - expected).abs() < 0.1, "sphere mesh volume {v} (expected {expected})");
        clear_tree(&sphere.0);
    }

    // -- general (non-analytic) curved-face boolean extensions --

    fn cylinder_solid(radius: f64, height: f64) -> Solid {
        BRepPrimCylinder::make_cylinder(radius, height).solid
    }

    /// A curved B-spline patch face `z = base + bump·((u−½)² + (v−½)²)` on
    /// `[0, 1]²` (a non-planar `Other` surface).
    fn curved_patch_face(base: f64, bump: f64) -> Face {
        let mut grid: Vec<Vec<GpPnt>> = Vec::new();
        for i in 0..=3usize {
            let mut row = Vec::new();
            for j in 0..=3usize {
                let u = i as f64 / 3.0;
                let v = j as f64 / 3.0;
                let z = base + bump * ((u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5));
                row.push(GpPnt::new(u, v, z));
            }
            grid.push(row);
        }
        let surf = crate::intpatch::make_bspline_surface_from_grid(&grid, 2, 2).expect("patch fit");
        let bld = TopoBuilder::new();
        bld.make_face(surf, &[])
    }

    #[test]
    fn general_box_cylinder_fuse() {
        // Box [-1,1]³ fused with a radius-0.4 cylinder on z∈[0,2] poking through
        // the top face: expected = 8 + π·0.4²·1 ≈ 8.50.
        let box_s = centered_box(1.0);
        let cyl = cylinder_solid(0.4, 2.0);
        let v = general_curved_boolean_volume(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse volume");
        let expected = 8.0 + PI * 0.16 * 1.0;
        assert!((v - expected).abs() < 0.5, "box-cylinder fuse volume {v} (expected {expected})");
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn general_sphere_cylinder_cut() {
        // Unit sphere minus a radius-0.4 cylinder on z∈[0,2]: the cylinder exits
        // the sphere at z = sqrt(1 − 0.4²) ≈ 0.9165.
        let sphere = unit_sphere_solid();
        let cyl = cylinder_solid(0.4, 2.0);
        let v = general_curved_boolean_volume(&sphere.0, &cyl.0, BoolOp::Cut, 1e-6).expect("cut volume");
        let z_exit = (1.0 - 0.4f64 * 0.4).sqrt();
        let expected = 4.0 / 3.0 * PI - PI * 0.16 * z_exit;
        assert!((v - expected).abs() < 0.5, "sphere-cylinder cut volume {v} (expected {expected})");
        clear_tree(&sphere.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn planar_delegates_to_bop_builder() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean_ext(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.solid.is_some(), "planar fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "planar fuse shell is closed");
        assert!(faces_of(&r.shape).len() > 12, "faces {}", faces_of(&r.shape).len());
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn classify_general_face_inside_outside() {
        let box_s = centered_box(1.0);
        let inside = curved_patch_face(0.8, 0.3); // z∈[0.8, 0.95] inside the box
        assert_eq!(classify_face_general(&inside, &box_s.0, 1e-6), FaceRegion::Inside);
        let outside = curved_patch_face(3.0, 0.3); // z∈[3.0, 3.15] above the box
        assert_eq!(classify_face_general(&outside, &box_s.0, 1e-6), FaceRegion::Outside);
        // A curved face crossing the box top (spanning inside and outside) is On.
        let crossing = curved_patch_face(0.9, 0.5); // z∈[0.9, 1.15]
        assert_eq!(classify_face_general(&crossing, &box_s.0, 1e-6), FaceRegion::On);
        clear_tree(&box_s.0);
        clear_tree(&inside.0);
        clear_tree(&outside.0);
        clear_tree(&crossing.0);
    }

    #[test]
    fn disjoint_general_fuse_compound() {
        // Box at z∈[-4,-2] far below a cylinder at z∈[0,2] → Fuse returns a
        // compound (bboxes do not overlap).
        let box_s = BRepPrimBox::make_box_corner(&GpPnt::new(-1.0, -1.0, -4.0), &GpPnt::new(1.0, 1.0, -2.0)).solid;
        let cyl = cylinder_solid(0.5, 2.0);
        let r = curved_boolean_ext(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(r.shape.is_compound(), "disjoint general fuse is a compound");
        assert!(r.solid.is_none(), "disjoint general fuse has no single solid");
        clear_tree(&r.shape);
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn general_boolean_volume_matches_curved_boolean_volume() {
        // Sphere-sphere dispatch stays on the analytic path: the general
        // dispatcher's volume agrees with curved_boolean_volume.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let ext_vol = general_curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ext volume");
        let ref_vol = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref volume");
        assert!((ext_vol - ref_vol).abs() < 0.05, "ext {ext_vol} vs curved_boolean_volume {ref_vol}");
        assert!((ext_vol - 8.018).abs() < 0.3, "ext volume {ext_vol} (expected ≈ 8.018)");
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // -- trimmed-face (real B-Rep) general boolean --

    #[test]
    fn trimmed_face_box_cylinder_fuse() {
        // Box [-1,1]³ fused with a radius-0.4 cylinder on z∈[0,2] poking through
        // the top: the cylinder lateral face is trimmed to the z>1 stub and keeps
        // its analytic cylinder surface (not a faceted mesh).
        let box_s = centered_box(1.0);
        let cyl = cylinder_solid(0.4, 2.0);
        let r = curved_boolean_full(&box_s.0, &cyl.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let has_cylinder = faces_of(&r.shape)
            .iter()
            .any(|f| BRepTool::face_surface(f).map(|s| surface_cylinder_params(s.as_ref()).is_some()).unwrap_or(false));
        assert!(has_cylinder, "a trimmed cylinder face is preserved");
        // `shape_volume` tessellates each face's full UV window (it ignores the
        // trimming wire and full-surface overlap), so the volume is a loose band —
        // it must still be in the ballpark of the box plus the protruding stub.
        let v = shape_volume(&r.shape, 0.02);
        let expected = 8.0 + PI * 0.16 * 1.0;
        assert!((v - expected).abs() < 5.0, "box-cylinder fuse volume {v} (expected ≈ {expected})");
        clear_tree(&r.shape);
        clear_tree(&box_s.0);
        clear_tree(&cyl.0);
    }

    #[test]
    fn trimmed_sphere_box_cut() {
        // Sphere radius 1 cut by a thin slab box z∈[0.5,0.9]: the retained sphere
        // caps keep their analytic sphere surface.
        let sphere = unit_sphere_solid();
        let slab = BRepPrimBox::make_box_corner(&GpPnt::new(-2.0, -2.0, 0.5), &GpPnt::new(2.0, 2.0, 0.9)).solid;
        let shape = general_boolean_trimmed(&sphere.0, &slab.0, BoolOp::Cut, 1e-6).expect("cut ok");
        let has_sphere = faces_of(&shape)
            .iter()
            .any(|f| BRepTool::face_surface(f).map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere).unwrap_or(false));
        assert!(has_sphere, "a retained sphere cap is preserved");
        clear_tree(&shape);
        clear_tree(&sphere.0);
        clear_tree(&slab.0);
    }

    #[test]
    fn trim_face_to_region_keeps_side() {
        // A sphere face trimmed by a small UV circle keeps the inside region:
        // the surface is still a radius-1 sphere, but the face's wire is the
        // smaller loop.
        let sphere = unit_sphere_solid();
        let face = faces_of(&sphere.0)[0].clone();
        let (u0, u1, v0, v1) = face_uv_window_local(&face);
        let (uc, vc) = (0.5 * (u0 + u1), 0.0);
        let radius = 0.9 * (v1 - v0).min(u1 - u0) / 2.0;
        let circle: Vec<GpPnt2d> = (0..32)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / 32.0;
                GpPnt2d::new(uc + radius * a.cos(), vc + radius * a.sin())
            })
            .collect();
        let trimmed = trim_face_to_region(&face, true, &circle, 1e-6).expect("trim ok").expect("some face");
        let s = BRepTool::face_surface(&trimmed).expect("surface");
        assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Sphere, "surface preserved");
        let (sp_c, sp_r) = crate::intpatch::sphere_params(s.as_ref()).expect("sphere params");
        assert!((sp_r - 1.0).abs() < 1e-6, "radius {sp_r}");
        assert!(sp_c.distance(&GpPnt::zero()) < 1e-6, "center {sp_c:?}");
        let w = wires_of_face(&trimmed);
        assert_eq!(w.len(), 1, "one boundary wire");
        let e = edges_of_wire(&w[0]);
        assert!(e.len() >= 3, "wire has edges");
        clear_tree(&trimmed.0);
        clear_tree(&sphere.0);
    }

    #[test]
    fn general_boolean_trimmed_closed_shell() {
        // Box [-1,1]³ fused with a small sphere poking through the top face → the
        // trimmed faces weld into a closed shell.
        let box_s = centered_box(1.0);
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, 1.2), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
        let bld = TopoBuilder::new();
        let face = bld.make_face(Arc::new(occt_geom::GeomSphere::new(occt_core::gp::GpSphere::new(ax3, 0.5).unwrap())), &[]);
        let shell_s = bld.make_shell(&[face]);
        let sphere = bld.make_solid(&[shell_s]);
        let shape = general_boolean_trimmed(&box_s.0, &sphere.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let shells = shapes_of(&shape, ShapeType::Shell);
        assert!(!shells.is_empty(), "result has shells");
        let shell = Shell(shells[0].clone());
        assert!(shell_is_closed(&shell), "fused box+sphere shell is closed");
        // The welded trimmed faces (box top ring + sphere cap) share every edge.
        let usage = crate::shell_check::edge_face_usage(&shell.0);
        assert!(usage.values().all(|&c| c == 2), "every boundary edge is used by exactly 2 faces");
        clear_tree(&shape);
        clear_tree(&box_s.0);
        clear_tree(&sphere.0);
    }

    #[test]
    fn curved_boolean_full_planar_delegates() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid;
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0)).solid;
        let r = curved_boolean_full(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let refr = crate::bop_builder::boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref ok");
        assert_eq!(faces_of(&r.shape).len(), faces_of(&refr.shape).len(), "planar delegates to bop_builder");
        assert!(shell_is_closed(&r.shells[0]), "planar fuse shell is closed");
        clear_tree(&r.shape);
        clear_tree(&refr.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn curved_boolean_full_quadric_unchanged() {
        // Two spheres (analytic quadrics) route to the unchanged curved_boolean.
        let a = unit_sphere_solid();
        let b = sphere_at(GpPnt::new(1.5, 0.0, 0.0));
        let r = curved_boolean_full(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let refr = curved_boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref ok");
        assert_eq!(faces_of(&r.shape).len(), faces_of(&refr.shape).len(), "quadrics route to curved_boolean");
        let ref_vol = curved_boolean_volume(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("ref vol");
        assert!((ref_vol - 8.018).abs() < 0.3, "volume {ref_vol} (expected ≈ 8.018)");
        clear_tree(&r.shape);
        clear_tree(&refr.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn trimmed_face_preserves_surface() {
        let sphere = unit_sphere_solid();
        let face = faces_of(&sphere.0)[0].clone();
        let orig = BRepTool::face_surface(&face).expect("original surface");
        let (u0, u1, v0, v1) = face_uv_window_local(&face);
        let (uc, vc) = (0.5 * (u0 + u1), 0.3);
        let radius = 0.7 * (v1 - v0).min(u1 - u0) / 2.0;
        let circle: Vec<GpPnt2d> = (0..32)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / 32.0;
                GpPnt2d::new(uc + radius * a.cos(), vc + radius * a.sin())
            })
            .collect();
        let trimmed = trim_face_to_region(&face, true, &circle, 1e-6).expect("trim ok").expect("some face");
        let s = BRepTool::face_surface(&trimmed).expect("trimmed surface");
        let p_orig = orig.d0(uc, vc);
        let p_new = s.d0(uc, vc);
        assert!(p_orig.distance(&p_new) < 1e-9, "surface evaluates identically at the face center");
        assert!((p_orig.distance(&GpPnt::zero()) - 1.0).abs() < 1e-6, "still on the unit sphere");
        clear_tree(&trimmed.0);
        clear_tree(&sphere.0);
    }


    #[test]
    fn boolean_result_shape_has_faces() {
        let box_s = BRepPrimBox::make_box(1.0, 2.0, 3.0).solid;
        let faces = faces_of(&box_s.0);
        let r = boolean_result_from_shape(box_s.0.clone(), faces.clone(), vec![]);
        assert_eq!(r.faces.len(), 6);
        assert!(r.solid.is_some(), "a box shape round-trips to a solid");
        assert!(!r.shells.is_empty(), "the solid carries a shell");
        assert_eq!(faces_of(&r.shape).len(), 6, "faces round-trip");
        clear_tree(&box_s.0);
        clear_tree(&r.shape);
    }
}

mod p01;
mod p02;
mod p03;
mod p04;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;
