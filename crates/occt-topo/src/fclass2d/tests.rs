use super::prelude::*;
use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_surface::{face_is_planar, surface_closest_params};
    use crate::builder::TopoBuilder;
    use crate::builder_face::build_face_with_holes;
    use crate::primitives::BRepPrimCylinder;
    use crate::topo_tools_full::faces_of;
    use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpDir, GpPln, GpPnt, GpPnt2d};
    use occt_geom::{GeomCylinder, GeomPlane, Surface};

    fn p2(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    /// UV parameters of a 3D point on a face surface (for tests).
    fn project_uv(face: &Face, p: &GpPnt) -> GpPnt2d {
        let surf = GeometryRegistry::global().face_surface(&face.0).expect("face surface");
        let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
        p2(u, v)
    }

    fn square_edges(b: &TopoBuilder, pts: &[GpPnt; 4]) -> Vec<Edge> {
        (0..4).map(|i| b.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
    }

    #[test]
    fn unit_box_bottom_face_classifies_square() {
        let ub = unit_box();
        let face = &ub.faces[0]; // bottom (z = 0), UV square [0,1]²
        let cl = FClass2d::new(face, 1e-6).expect("classifier");
        assert_eq!(cl.perform(p2(0.5, 0.5)), FaceState::In);
        assert_eq!(cl.perform(p2(-0.5, 0.5)), FaceState::Out);
        assert_eq!(cl.perform(p2(0.5, 0.0)), FaceState::On);
        // The open bottom face has no closed periodic boundary: the infinite
        // point is Out.
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn face_with_square_hole() {
        let b = TopoBuilder::new();
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let hole = [
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let face = build_face_with_holes(&square_edges(&b, &outer), &[square_edges(&b, &hole)])
            .expect("face with hole");
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");

        // Hole centre -> Out; hole left-edge midpoint -> On; annulus -> In;
        // in the annulus near the outer corner -> In; outside the outer square
        // (projection clamps the grid to [-1,1], so use a direct UV point) -> Out.
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.5, 0.5, 0.0))), FaceState::Out);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.25, 0.5, 0.0))), FaceState::On);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.9, 0.9, 0.0))), FaceState::In);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.05, 0.05, 0.0))), FaceState::In);
        assert_eq!(cl.perform(p2(1.5, 0.5)), FaceState::Out);
    }

    /// A clean cylinder lateral face: two cap circles plus two *distinct* seam
    /// edges whose stored pcurves sit at u = 0 (up) and u = 2π (down). The
    /// port's flat wire model cannot store the reversed seam orientation, so
    /// the pcurves are attached explicitly to make the UV ring a clean
    /// `[0, 2π] × [0, h]` rectangle.
    fn clean_cylinder_lateral(radius: f64, height: f64) -> Face {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("cylinder axis");
        let surface: Arc<dyn Surface> = Arc::new(GeomCylinder::new(
            GpCylinder::new(ax, radius).expect("cylinder radius"),
        ));

        let bottom = GpPnt::new(radius, 0.0, 0.0);
        let top = GpPnt::new(radius, 0.0, height);
        let ax2_bot = GpAx2::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("bottom circle axis");
        let ax2_top = GpAx2::new(GpPnt::new(0.0, 0.0, height), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("top circle axis");
        let bottom_circle = b.make_edge_circle(&ax2_bot, radius, 0.0, 2.0 * PI);
        let top_circle = b.make_edge_circle(&ax2_top, radius, 0.0, 2.0 * PI);
        let seam_up = b.make_edge_segment(&bottom, &top);
        let seam_down = b.make_edge_segment(&bottom, &top);

        let lateral = b.make_face(surface, &[]);
        let face_key = GeometryRegistry::shape_key(&lateral.0);
        let pc_up: Arc<dyn Curve2d> = Arc::new(occt_geom2d::Geom2dLine::from_pnt_dir(
            p2(0.0, 0.0),
            occt_core::gp::GpDir2d::new(0.0, 1.0).unwrap(),
        ));
        let pc_down: Arc<dyn Curve2d> = Arc::new(occt_geom2d::Geom2dLine::from_pnt_dir(
            p2(2.0 * PI, height),
            occt_core::gp::GpDir2d::new(0.0, -1.0).unwrap(),
        ));
        GeometryRegistry::global().set_edge_pcurve(&seam_up.0, face_key, pc_up);
        GeometryRegistry::global().set_edge_pcurve(&seam_down.0, face_key, pc_down);

        let wire = b.make_wire(&[bottom_circle, seam_up, top_circle, seam_down]);
        let mut lateral = lateral;
        b.add_wire(&mut lateral, &wire);
        lateral
    }

    #[test]
    fn cylinder_lateral_periodic_folding() {
        let radius = 1.0;
        let height = 2.0;
        let face = clean_cylinder_lateral(radius, height);
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");
        let h = height;

        // Interior point is In, and its (u + 2π) periodic image classifies the
        // same.
        let u = 0.5;
        let v = 0.25 * h;
        assert_eq!(cl.perform(p2(u, v)), FaceState::In);
        assert_eq!(cl.perform(p2(u + 2.0 * PI, v)), FaceState::In);
        // A domain-outside left point folds by +2π into the face -> In.
        assert_eq!(cl.perform(p2(-PI, v)), FaceState::In);
        assert_eq!(cl.perform(p2(-PI, v)), cl.perform(p2(PI, v)));
        // Outside above the top edge (v not periodic) stays Out in both images.
        assert_eq!(cl.perform(p2(u, 1.5 * h)), FaceState::Out);
        assert_eq!(cl.perform(p2(u + 2.0 * PI, 1.5 * h)), FaceState::Out);

        // The lateral face is open (bounded caps): the infinite point is Out.
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn cylinder_primitive_lateral_periodic_equality() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let lateral = faces_of(&c.solid.0).into_iter().find(|f| !face_is_planar(f)).expect("lateral");
        let cl = FClass2d::new(&lateral, 1e-6).expect("classifier");
        let u = 0.5;
        let v = 0.5;
        // The periodic-image equality holds even for the primitive's (seam-
        // duplicated) lateral wire.
        assert_eq!(cl.perform(p2(u, v)), cl.perform(p2(u + 2.0 * PI, v)));
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn is_hole_reflects_winding() {
        let b = TopoBuilder::new();
        // CCW ring in the standard-plane UV -> not a hole.
        let ccw = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let pln = GpPln::new(GpAx3::standard());
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let f_ccw = crate::builder_face::make_face_from_wire(&square_edges(&b, &ccw), Some(surf))
            .expect("ccw face");
        let cl = FClass2d::new(&f_ccw, 1e-6).expect("classifier");
        assert!(!cl.is_hole(), "CCW ring is a bounded face, not a hole");

        // CW ring in the same UV -> a hole.
        let cw = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln));
        let f_cw = crate::builder_face::make_face_from_wire(&square_edges(&b, &cw), Some(surf))
            .expect("cw face");
        let cl = FClass2d::new(&f_cw, 1e-6).expect("classifier");
        assert!(cl.is_hole(), "CW ring reports a hole face");
    }

    #[test]
    fn infinite_point_closed_periodic_face_is_in() {
        // A full sphere face has no boundary wires -> the classifier sees an
        // empty boundary and the infinite point is In (closed periodic face).
        let sph = crate::primitives::BRepPrimSphere::make_sphere(2.0);
        let face = faces_of(&sph.solid.0).into_iter().next().expect("sphere face");
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");
        assert_eq!(cl.perform_infinite_point(), FaceState::In);
        // Any 2D point is In (no restriction).
        assert_eq!(cl.perform(p2(1.0, 0.5)), FaceState::In);
    }

    #[test]
    fn test_on_restriction_detects_boundary() {
        let ub = unit_box();
        let face = &ub.faces[0];
        let cl = FClass2d::new(face, 1e-6).expect("classifier");
        // A point exactly on the boundary is On; an interior point is In.
        assert_eq!(cl.test_on_restriction(p2(0.5, 0.0), 1e-6), FaceState::On);
        assert_eq!(cl.test_on_restriction(p2(0.5, 0.5), 1e-6), FaceState::In);
        assert_eq!(cl.test_on_restriction(p2(1.5, 0.5), 1e-6), FaceState::Out);
    }
