use super::prelude::*;
use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use occt_core::gp::{GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus};
    use occt_geom::{GeomBSplineCurve, GeomCone, GeomCylinder, GeomSphere, GeomTorus};

    fn cyl_face(radius: f64) -> Face {
        let b = TopoBuilder::new();
        let cyl = GpCylinder::new(GpAx3::standard(), radius).unwrap();
        b.make_face(Arc::new(GeomCylinder::new(cyl)), &[])
    }

    #[test]
    fn classify_surface_kinds() {
        let b = TopoBuilder::new();
        let pln = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let s = GeometryRegistry::global().face_surface(&pln.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Plane);

        let fc = cyl_face(1.0);
        let s = GeometryRegistry::global().face_surface(&fc.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Cylinder);

        let sph = GpSphere::new(GpAx3::standard(), 2.0).unwrap();
        let fs = b.make_face(Arc::new(GeomSphere::new(sph)), &[]);
        let s = GeometryRegistry::global().face_surface(&fs.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Sphere);

        let cone = GpCone::new(GpAx3::standard(), 1.0, 0.5f64.atan()).unwrap();
        let fcn = b.make_face(Arc::new(GeomCone::new(cone)), &[]);
        let s = GeometryRegistry::global().face_surface(&fcn.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Cone);

        let tor = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let ft = b.make_face(Arc::new(GeomTorus::new(tor)), &[]);
        let s = GeometryRegistry::global().face_surface(&ft.0).unwrap();
        assert_eq!(classify_surface_kind(s.as_ref()), SurfaceKind::Torus);
    }

    #[test]
    fn unit_box_line_edge_is_line_pcurve_with_boundary_endpoints() {
        let b = unit_box();
        let edge = &b.edges[4]; // (0,0,1) → (1,0,1)
        let face = &b.faces[1]; // top (z = 1)
        let pc = make_pcurve_full(edge, face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn cylinder_generatrix_is_v_line() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        // Generatrix at angle 0: u = atan2(0, 1) = 0, v = axial coordinate.
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 0.0).abs() < 1e-6 && (p1.y() - 1.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn cylinder_top_circle_is_u_line_across_the_seam() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(
            GpPnt::new(0.0, 0.0, 1.0),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let edge = b.make_edge_circle(&ax2, 1.0, 0.0, 2.0 * PI);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        // The cap circle is a v = 1 isoline in the U-periodic direction.
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::Line);
        let p0 = pc.d0(0.0);
        let pm = pc.d0(PI);
        let p2 = pc.d0(2.0 * PI);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 1.0).abs() < 1e-6, "start {:?}", p0);
        assert!((pm.x() - PI).abs() < 1e-6 && (pm.y() - 1.0).abs() < 1e-6, "mid {:?}", pm);
        assert!((p2.x() - 2.0 * PI).abs() < 1e-6 && (p2.y() - 1.0).abs() < 1e-6, "end {:?}", p2);
    }

    #[test]
    fn general_bspline_edge_on_plane_face_samples_boundary() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 0.4, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap());
        let edge = b.make_edge(curve, 0.0, 1.0);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(crate::pcurve::pc_curve_kind(pc.as_ref()), crate::pcurve::CurveKind::BSpline);
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(1.0);
        assert!((p0.x() - 0.0).abs() < 1e-6 && (p0.y() - 0.0).abs() < 1e-6, "start {:?}", p0);
        assert!((p1.x() - 1.0).abs() < 1e-6 && (p1.y() - 0.0).abs() < 1e-6, "end {:?}", p1);
    }

    #[test]
    fn agrees_with_make_pcurve_on_face_for_plane_line() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let full = make_pcurve_full(edge, face).expect("full pcurve");
        let base = crate::pcurve::make_pcurve_on_face(edge, face).expect("base pcurve");
        let (a, z) = GeometryRegistry::global().edge_parameters(&edge.0);
        for i in 0..=8 {
            let t = a + (z - a) * i as f64 / 8.0;
            let q1 = full.d0(t);
            let q2 = base.d0(t);
            assert!((q1.x() - q2.x()).abs() < 1e-9 && (q1.y() - q2.y()).abs() < 1e-9, "t {t}");
        }
    }

    #[test]
    fn trim_shifts_periodic_pcurve_into_face_bounds() {
        let face = cyl_face(1.0);
        // A degree-1 B-spline with u in [6.5, 7.5] at v = 0: u is outside the
        // face's [0, 2π] range.
        let pts = [
            GpPnt2d::new(6.5, 0.0),
            GpPnt2d::new(7.0, 0.0),
            GpPnt2d::new(7.5, 0.0),
        ];
        let bs = bspline_from_samples(&pts, 6.5, 7.5, 1).unwrap();
        let pc: Arc<dyn Curve2d> = Arc::new(bs);
        let trimmed = trim_pcurve_to_face(&pc, &face, 1e-7).expect("trimmed");
        let (umin, umax, vmin, vmax) = face_uv_bounds(&face);
        for i in 0..=16 {
            let t = 6.5 + (7.5 - 6.5) * i as f64 / 16.0;
            let q = trimmed.d0(t);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {} at t {}", q.x(), t);
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {} at t {}", q.y(), t);
        }
    }

    #[test]
    fn orientation_matches_edge_direction() {
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let edge = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        assert_eq!(pc_curve_orientation(&edge, &face, pc.as_ref()), 1.0);
        // A reversed edge runs against the stored curve direction.
        let mut rev = edge.clone();
        rev.0.set_orientation(crate::abs::Orientation::Reversed);
        assert_eq!(pc_curve_orientation(&rev, &face, pc.as_ref()), -1.0);
    }

    #[test]
    fn pcurve_point_evaluation() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        let q = pcurve_point_on_face(pc.as_ref(), 0.5).expect("point");
        // The unit-speed Geom2dLine pcurve maps edge parameter → (u, 0).
        assert!((q.x() - 0.5).abs() < 1e-6 && (q.y() - 0.0).abs() < 1e-6, "q {:?}", q);

        // A finite-range B-spline pcurve rejects out-of-range parameters.
        let poles = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.5, 0.4, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap());
        let be = b.make_edge(curve, 0.0, 1.0);
        let bpc = make_pcurve_full(&be, &face).expect("bspline pcurve");
        assert!(pcurve_point_on_face(bpc.as_ref(), 5.0).is_err(), "out of range must fail");
    }

    #[test]
    fn non_isoparametric_circle_on_cylinder_falls_back_to_sampling() {
        // A circle in the X-Z plane (not a cap circle, and only touching the
        // cylinder at the two axis points) is not an isoparametric curve, so
        // the pcurve is built by the sampling path. Its endpoints still land on
        // the projected surface parameters.
        let face = cyl_face(1.0);
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(
            GpPnt::zero(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let edge = b.make_edge_circle(&ax2, 1.0, 0.0, PI);
        let pc = make_pcurve_full(&edge, &face).expect("pcurve");
        let p0 = pc.d0(0.0);
        let p1 = pc.d0(PI);
        // Endpoint (1, 0, 0) → (u=0, v=0); endpoint (−1, 0, 0) → (u=π, v=0).
        assert!((p0.x() - 0.0).abs() < 1e-5 && (p0.y() - 0.0).abs() < 1e-5, "start {:?}", p0);
        assert!((p1.x() - PI).abs() < 1e-5 && (p1.y() - 0.0).abs() < 1e-5, "end {:?}", p1);
    }
