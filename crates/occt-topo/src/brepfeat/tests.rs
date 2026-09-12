use super::prelude::*;
use super::*;
    use crate::brep_surface::{classify_surface, surface_normal, SurfaceKind};
    use crate::brep_tool::BRepTool;
    use crate::primitives::{BRepPrimBox, BRepPrimCylinder};
    use crate::topo_tools_full::{faces_of, vertex_position, vertices_of};
    use occt_core::gp::{GpAx3, GpXyz};

    fn z_axis() -> GpAx1 {
        GpAx1::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap())
    }

    fn xz_pivot() -> GpPln {
        GpPln::new(
            GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .expect("xz pivot"),
        )
    }

    fn xy_plane_at(z: f64) -> GpPln {
        GpPln::new(
            GpAx3::new(
                GpPnt::new(0.0, 0.0, z),
                GpDir::new(0.0, 0.0, 1.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .expect("xy plane"),
        )
    }

    /// A cylinder target built as a faceted mesh (planar faces), because the
    /// analytic `BRepPrimCylinder` lateral surface cannot be meshed by the
    /// wireframe mesher (`GeomCylinder::d1` returns zero partials).
    fn faceted_cyl(r: f64, h: f64) -> Solid {
        let tool = mesh_cylinder(r, h, 24);
        Solid::wrap(tool).expect("faceted cylinder solid")
    }

    #[test]
    fn draft_zero_angle_noop() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let n = faces_of(&b.solid.0).len();
        let r = draft(&b.solid, &[GpPnt::new(0.5, 0.5, 1.0)], 0.0, &xz_pivot(), 1e-6).unwrap();
        assert_eq!(faces_of(&r.shape).len(), n, "draft at angle 0 leaves the solid unchanged");
        // The top face is still the horizontal z = 1 plane.
        let faces = faces_of(&r.shape);
        let top = faces
            .iter()
            .find(|f| face_point_distance(f, &GpPnt::new(0.5, 0.5, 1.0)) < 1e-6)
            .expect("top face");
        let s = BRepTool::face_surface(top).unwrap();
        assert!((surface_normal(s.as_ref(), 0.5, 0.5).z() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn draft_rotates_face_geometry() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let r = draft(&b.solid, &[GpPnt::new(0.5, 0.5, 1.0)], 30.0, &xz_pivot(), 1e-6).unwrap();
        // The drafted face is still planar but leans: a plane face whose normal
        // keeps a strong +Z component yet is tilted away from +Z.
        let faces = faces_of(&r.shape);
        let drafted = faces
            .iter()
            .find(|f| {
                let Some(s) = BRepTool::face_surface(f) else { return false };
                if classify_surface(s.as_ref()) != SurfaceKind::Plane {
                    return false;
                }
                let n = surface_normal(s.as_ref(), 0.5, 0.5);
                n.z() > 0.5 && (n.x().abs() + n.y().abs()) > 0.1
            })
            .expect("a drafted (tilted) planar face exists");
        let s = BRepTool::face_surface(drafted).unwrap();
        let n = surface_normal(s.as_ref(), 0.5, 0.5);
        assert!(classify_surface(s.as_ref()) == SurfaceKind::Plane, "drafted face stays planar");
        assert!(
            (n.xyz().z - 1.0).abs() > 0.05,
            "drafted normal {n:?} differs from the original +Z"
        );
        assert_eq!(faces_of(&r.shape).len(), 6, "draft keeps the box's face count");
    }

    #[test]
    fn groove_cuts_cylinder() {
        let cyl = faceted_cyl(1.0, 3.0);
        let before = solid_volume(&cyl.0, 40, 40);
        // Annular groove ring: cut a 0.4-deep band out of the wall (r 0.6..1)
        // between z 1.0 and 1.8.
        let profile = [
            GpPnt2d::new(0.6, 1.0),
            GpPnt2d::new(1.3, 1.0),
            GpPnt2d::new(1.3, 1.8),
            GpPnt2d::new(0.6, 1.8),
            GpPnt2d::new(0.6, 1.0),
        ];
        let after = groove(&cyl, &profile, &z_axis(), 16, 1e-4).unwrap();
        let groove_vol = std::f64::consts::PI * (1.0 - 0.36) * 0.8;
        assert!(
            after.volume < before,
            "groove must remove material: after {} < before {}",
            after.volume,
            before
        );
        assert!(
            (after.volume - (before - groove_vol)).abs() < 0.5,
            "grooved volume {} vs expected {}",
            after.volume,
            before - groove_vol
        );
    }

    #[test]
    fn neck_adds_material() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let before = crate::shape_mesh::shape_volume(&box_s.solid.0, 0.05);
        // Neck: cylinder r 0.4 standing from the box top (z 1) to z 2.0,
        // centered on the box.
        let profile = [
            GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(0.4, 1.0),
            GpPnt2d::new(0.4, 2.0),
            GpPnt2d::new(0.0, 2.0),
        ];
        let axis = GpAx1::new(GpPnt::new(1.0, 1.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let after = neck(&box_s.solid, &profile, &axis, 16, 1e-4).unwrap();
        let added = std::f64::consts::PI * 0.16 * (2.0 - 1.0);
        assert!(after.volume > before, "neck must add material: {} > {}", after.volume, before);
        assert!(
            (after.volume - before - added).abs() < 0.6,
            "neck added {} vs expected {}",
            after.volume - before,
            added
        );
    }

    #[test]
    fn rib_fuses_to_box() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let before = crate::shape_mesh::shape_volume(&box_s.solid.0, 0.05);
        let before_faces = faces_of(&box_s.solid.0).len();
        let plane = xy_plane_at(1.0);
        let profile = [
            GpPnt2d::new(0.5, 0.5),
            GpPnt2d::new(0.8, 0.5),
            GpPnt2d::new(0.8, 0.8),
            GpPnt2d::new(0.5, 0.8),
        ];
        let after = rib(&box_s.solid, &profile, &plane, 0.3, 1.0, 0.05).unwrap();
        let added = 0.3 * 0.3 * 1.0;
        assert!(after.volume > before, "rib must add material: {} > {}", after.volume, before);
        assert!(
            (after.volume - (before + added)).abs() < 0.08,
            "rib added {} vs expected {}",
            after.volume - before,
            added
        );
        assert!(
            faces_of(&after.shape).len() > before_faces,
            "rib grows the face count ({} > {})",
            faces_of(&after.shape).len(),
            before_faces
        );
    }

    #[test]
    fn boss_thru_all_pierces() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 0.5);
        let before = crate::shape_mesh::shape_volume(&box_s.solid.0, 0.05);
        let before_faces = faces_of(&box_s.solid.0).len();
        let center = GpPnt::new(1.0, 1.0, 0.25);
        let after = boss_thru_all(&box_s.solid, &center, 0.3, 0.15).unwrap();
        assert!(after.volume > before, "boss must add material: {} > {}", after.volume, before);
        // The boss pierces the box: the result's bounding box is taller than the
        // box, and it carries a cylinder-like wall — a face whose vertices all
        // sit a boss-radius away from the boss axis (a voxel-faceted boss wall).
        let (_, _, _, _, bz0, bz1) = crate::bbox_from_geometry::shape_bbox(&after.shape).get().unwrap();
        assert!(
            bz1 - bz0 > 0.5 + 1.0,
            "boss protrudes through the box (bbox z-span {})",
            bz1 - bz0
        );
        let has_wall = faces_of(&after.shape).iter().any(|f| {
            let vs = vertices_of(&f.0);
            vs.iter().all(|v| {
                let p = vertex_position(v);
                let r = ((p.x() - 1.0).powi(2) + (p.y() - 1.0).powi(2)).sqrt();
                r > 0.1 && r < 0.6
            })
        });
        assert!(has_wall, "boss adds a cylinder-like wall face");
        assert!(faces_of(&after.shape).len() > before_faces, "boss grows the face count");
    }

    #[test]
    fn revolve_profile_creates_tool() {
        let profile = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(0.5, 0.0),
            GpPnt2d::new(0.5, 2.0),
            GpPnt2d::new(0.0, 2.0),
        ];
        let rev = crate::sweep_revolve::revolve_polyline_around_z(&profile, 16).unwrap();
        let v = crate::sweep_revolve::revolved_volume(&rev);
        assert!(
            (v - std::f64::consts::PI * 0.25 * 2.0).abs() < 1e-9,
            "revolved volume {v}"
        );
    }

    #[test]
    fn groove_negative_volume_delta() {
        let cyl = faceted_cyl(1.0, 3.0);
        let profile = [
            GpPnt2d::new(0.6, 1.0),
            GpPnt2d::new(1.3, 1.0),
            GpPnt2d::new(1.3, 1.8),
            GpPnt2d::new(0.6, 1.8),
            GpPnt2d::new(0.6, 1.0),
        ];
        let after = groove(&cyl, &profile, &z_axis(), 16, 1e-4).unwrap();
        let delta = feature_before_after_delta(&cyl, &after);
        assert!(delta > 1e-6, "groove changes the volume (delta {delta})");
    }

    #[test]
    fn neck_positive_volume_delta() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let profile = [
            GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(0.4, 1.0),
            GpPnt2d::new(0.4, 2.0),
            GpPnt2d::new(0.0, 2.0),
        ];
        let axis = GpAx1::new(GpPnt::new(1.0, 1.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let after = neck(&box_s.solid, &profile, &axis, 16, 1e-4).unwrap();
        let delta = feature_before_after_delta(&box_s.solid, &after);
        assert!(delta > 1e-6, "neck changes the volume (delta {delta})");
    }

    #[test]
    fn draft_face_found_by_point() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&b.solid.0);
        let p = GpPnt::new(0.5, 0.5, 1.0);
        let found = faces
            .iter()
            .min_by(|a, b| face_point_distance(a, &p).partial_cmp(&face_point_distance(b, &p)).unwrap())
            .expect("a nearest face");
        assert!(face_point_distance(found, &p) < 1e-6, "face found by point is on the surface");
    }

    #[test]
    fn invalid_profile_errors() {
        let cyl = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let axis = z_axis();
        let plane = xy_plane_at(1.0);
        assert!(groove(&cyl.solid, &[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.5, 0.0)], &axis, 16, 0.05).is_err());
        assert!(neck(&cyl.solid, &[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.5, 0.0)], &axis, 16, 0.05).is_err());
        assert!(rib(&cyl.solid, &[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.5, 0.0)], &plane, 0.3, 1.0, 0.05).is_err());
    }

    #[test]
    fn rib_through_wire_ok() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let before = crate::shape_mesh::shape_volume(&box_s.solid.0, 0.05);
        let plane = xy_plane_at(1.0);
        let profile = [
            GpPnt2d::new(0.5, 0.5),
            GpPnt2d::new(0.8, 0.5),
            GpPnt2d::new(0.8, 0.8),
            GpPnt2d::new(0.5, 0.8),
        ];
        let after = rib(&box_s.solid, &profile, &plane, 0.3, 1.0, 0.05).unwrap();
        assert!(after.volume > before, "rib through the profile succeeds and adds volume");
    }

    #[test]
    fn boss_adds_material() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let before = crate::shape_mesh::shape_volume(&box_s.solid.0, 0.05);
        let after = boss(&box_s.solid, &GpPnt::new(1.0, 1.0, 0.0), 0.25, 0.8, 1e-4).unwrap();
        let added = std::f64::consts::PI * 0.25 * 0.25 * 0.8;
        assert!(after.volume > before, "boss must add material: {} > {}", after.volume, before);
        assert!(
            (after.volume - before - added).abs() < 0.6,
            "boss added {} vs expected {}",
            after.volume - before,
            added
        );
    }

    #[test]
    fn rib_volume_analytic_matches_tool() {
        let profile = [
            GpPnt2d::new(0.5, 0.5),
            GpPnt2d::new(0.8, 0.5),
            GpPnt2d::new(0.8, 0.8),
            GpPnt2d::new(0.5, 0.8),
        ];
        let v = rib_volume(&profile, 0.3, 1.0).unwrap();
        assert!((v - 0.09).abs() < 1e-9, "rib volume {v}");
    }

    #[test]
    fn revolve_profile_validation() {
        assert!(validate_revolve_profile(&[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.5, 0.0)], 16).is_err());
        assert!(validate_revolve_profile(
            &[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.5, 0.0), GpPnt2d::new(0.5, 1.0)],
            3,
        )
        .is_err());
        assert!(validate_revolve_profile(
            &[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(-0.5, 0.0), GpPnt2d::new(0.5, 1.0)],
            16,
        )
        .is_err());
        let v = revolved_profile_volume(
            &[
                GpPnt2d::new(0.0, 0.0),
                GpPnt2d::new(0.5, 0.0),
                GpPnt2d::new(0.5, 2.0),
                GpPnt2d::new(0.0, 2.0),
            ],
            16,
        )
        .unwrap();
        assert!((v - std::f64::consts::PI * 0.25 * 2.0).abs() < 1e-9, "revolved volume {v}");
    }

    #[test]
    fn draft_hinge_axis_helper() {
        let top = GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap());
        let pivot = xz_pivot();
        let hinge = draft_hinge_axis(&top, &pivot).unwrap();
        // The hinge lies in the top face plane (z = 1) and runs along X.
        let p = hinge.location();
        assert!((p.z() - 1.0).abs() < 1e-9, "hinge on the top face z = 1");
        let x = GpXyz::new(1.0, 0.0, 0.0);
        let mx = GpXyz::new(-1.0, 0.0, 0.0);
        assert!(
            hinge.direction().xyz().crossed(&x).modulus() < 1e-9
                || hinge.direction().xyz().crossed(&mx).modulus() < 1e-9
        );
    }

    #[test]
    fn feature_kind_classifies_add_subtract() {
        let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
        let profile = [
            GpPnt2d::new(0.0, 1.0),
            GpPnt2d::new(0.4, 1.0),
            GpPnt2d::new(0.4, 2.0),
            GpPnt2d::new(0.0, 2.0),
        ];
        let axis = GpAx1::new(GpPnt::new(1.0, 1.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let after = neck(&box_s.solid, &profile, &axis, 16, 1e-4).unwrap();
        assert_eq!(after.kind(&box_s.solid), FeatKind::Additive);
        assert!(after.volume_delta(&box_s.solid) > 0.0);
        assert!(after.volume_ratio(&box_s.solid) > 1.0);
    }
