use super::prelude::*;
use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpPln, GpSphere, GpTorus};
    use occt_geom::{GeomBSplineCurve, GeomCircle, GeomCone, GeomCylinder, GeomLine, GeomPlane, GeomSphere, GeomTorus};

    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_tool::BRepTool;

    const PI: f64 = std::f64::consts::PI;

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GeomSphere {
        let ax3 = GpAx3::new(origin, dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        GeomSphere::new(GpSphere::new(ax3, r).unwrap())
    }

    #[test]
    fn line_through_unit_sphere_two_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -2.0), dir(0.0, 0.0, 1.0));
        // z = −1 at t = 1, z = +1 at t = 3.
        let res = perform_curve_surface(&line, &sphere, (-2.0, 4.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        let p0 = res.point(0).pnt();
        let p1 = res.point(1).pnt();
        assert!((p0.coord.modulus() - 1.0).abs() < 1e-9, "p0 {:?}", p0);
        assert!((p1.coord.modulus() - 1.0).abs() < 1e-9, "p1 {:?}", p1);
        assert!((p0.z() + 1.0).abs() < 1e-9 || (p1.z() + 1.0).abs() < 1e-9);
        assert!(res.point(0).param() < res.point(1).param());
    }

    #[test]
    fn line_outside_sphere_no_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(3.0, 0.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &sphere, (-1.0, 1.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
    }

    #[test]
    fn line_through_plane_single_point() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &plane, (-1.0, 1.0), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9, "p {:?}", p);
    }

    #[test]
    fn line_in_plane_is_on_segment() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        // Line lying in z = 0.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &plane, (-2.0, 2.0), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_in_plane_is_on_segment() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let res = perform_curve_surface(&circle, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_parallel_to_plane_no_intersection() {
        // Circle in z = 1, plane z = 0.
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 1.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 0);
    }

    #[test]
    fn circle_crosses_sphere_two_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Circle radius 1 centered (1.5, 0, 0) in the XY plane.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!((p.coord.modulus() - 1.0).abs() < 1e-9, "p {:?}", p);
        }
    }

    #[test]
    fn circle_inside_sphere_no_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 0.5));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 0);
    }

    #[test]
    fn bspline_through_sphere_general_path() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // A B-spline through (−2,0,0) and (2,0,0) with a bow toward z = 0.5;
        // it crosses the unit sphere near x = ±1.
        let poles = vec![
            GpPnt::new(-2.0, 0.0, 0.0),
            GpPnt::new(-0.5, 0.0, 0.6),
            GpPnt::new(0.5, 0.0, 0.6),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        // Degree-2 clamped B-spline: 4 poles → 4 + 2 + 1 = 7 knots.
        let curve = Arc::new(GeomBSplineCurve::new(poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap());
        let res = perform_curve_surface(curve.as_ref(), &sphere, (0.0, 1.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!((p.coord.modulus() - 1.0).abs() < 1e-5, "p {:?}", p);
        }
    }

    #[test]
    fn box_face_surface_crossed_by_line() {
        let b = unit_box();
        // Bottom face (z = 0).
        let face = &b.faces[0];
        let surf = BRepTool::face_surface(face).expect("face surface");
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.5, 0.5, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, surf.as_ref(), (-1.0, 1.0), (0.0, 1.0, 0.0, 1.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6, "p {:?}", p);
    }

    #[test]
    fn box_face_surface_missed_by_line() {
        let b = unit_box();
        let face = &b.faces[0]; // bottom, z = 0
        let surf = BRepTool::face_surface(face).expect("face surface");
        // Line outside the face's UV domain (x beyond the box).
        let line = GeomLine::from_pnt_dir(GpPnt::new(5.0, 5.0, -1.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, surf.as_ref(), (-1.0, 1.0), (0.0, 1.0, 0.0, 1.0)).unwrap();
        assert_eq!(res.nb_points(), 0, "points: {:?}", res.points());
    }

    #[test]
    fn tangent_line_touches_sphere_single_on_point() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Line at y = 1, tangent to the unit sphere at (0, 1, 0).
        let line = GeomLine::from_pnt_dir(GpPnt::new(-2.0, 1.0, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &sphere, (-2.0, 2.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!(p.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6, "p {:?}", p);
        assert_eq!(res.point(0).state(), State::On);
    }

    #[test]
    fn ellipse_crosses_plane_two_points() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        // Ellipse centered (0,0,0.5), semi-axes 2 (X) and 1 (Z), lying in the
        // XZ plane: d0(u) = (2·cos u, 0, 0.5 + sin u). It crosses the plane
        // z = 0 where sin u = −0.5 → u = 7π/6, 11π/6, giving two points.
        let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.5), dir(0.0, 1.0, 0.0), dir(1.0, 0.0, 0.0)).unwrap();
        let elips = occt_core::gp::GpElips::new(ax2, 2.0, 1.0);
        let ellipse = occt_geom::GeomEllipse::new(elips);
        let res = perform_curve_surface(&ellipse, &plane, (0.0, 2.0 * PI), (-5.0, 5.0, -5.0, 5.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            assert!(p.z().abs() < 1e-9, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_cylinder_two_points() {
        // Cylinder radius 1 along Z; line at y = 0.5, z = 0 crosses at
        // x = ±√(1 − 0.25) = ±0.866.
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let surface = GeomCylinder::new(cyl);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-2.0, 0.5, 0.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &surface, (-2.0, 4.0), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - 1.0).abs() < 1e-6, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_cone_two_points() {
        // Cone (geometric tip at (0,0,−1), radius 1 at z = 0, half-angle 45°):
        // radius at z = 1 is 2, so the line along X at z = 1 crosses at x = ±2.
        // Placement at the origin (radius 1 at z = 0); the tip is then at
        // z = −RefRadius/tan(45°) = −1.
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0)).unwrap();
        let cone = GpCone::new(ax3, 1.0, PI / 4.0).unwrap();
        let surface = GeomCone::new(cone);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-3.0, 0.0, 1.0), dir(1.0, 0.0, 0.0));
        let res = perform_curve_surface(&line, &surface, (-2.0, 6.0), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            // On the cone: radius = z + 1 (tip at z = −1, slope 1).
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - (p.z() + 1.0)).abs() < 1e-6, "p {:?}", p);
        }
    }

    #[test]
    fn line_through_torus_four_points() {
        // Torus at origin R = 3, r = 1; line along X through the centre cuts at
        // x = ±2, ±4 (port of IntAna_IntLinTorus).
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let surface = GeomTorus::new(torus);
        let line = GeomLine::from_pnt_dir(GpPnt::new(-6.0, 0.0, 0.0), dir(1.0, 0.0, 0.0));
        // x = −6 + t cuts the torus at t = 2, 4, 8, 10 (x = −4, −2, 2, 4).
        let res = perform_curve_surface(&line, &surface, (0.0, 12.0), (0.0, 2.0 * PI, 0.0, 2.0 * PI)).unwrap();
        assert_eq!(res.nb_points(), 4, "points: {:?}", res.points());
        let mut xs: Vec<f64> = res.points().iter().map(|p| p.pnt().x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (x, expected) in xs.iter().zip([-4.0, -2.0, 2.0, 4.0].iter()) {
            assert!((x - expected).abs() < 1e-5, "x {x} expected {expected}");
        }
    }

    #[test]
    fn line_entering_sphere_state_in_out() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, -2.0), dir(0.0, 0.0, 1.0));
        let res = perform_curve_surface(&line, &sphere, (-2.0, 4.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 2);
        assert_eq!(res.point(0).state(), State::In, "first crossing enters");
        assert_eq!(res.point(1).state(), State::Out, "second crossing exits");
    }

    #[test]
    fn circle_on_sphere_is_on_segment() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Unit circle in the XY plane = the unit sphere's equator.
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let res = perform_curve_surface(&circle, &sphere, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 0);
        assert_eq!(res.nb_segments(), 1, "segments: {:?}", res.segments());
        assert_eq!(res.segment(0).first_point().state(), State::On);
    }

    #[test]
    fn circle_crosses_cylinder_general_path() {
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let surface = GeomCylinder::new(cyl);
        // Circle radius 1 centred (1.5, 0, 0) in the XY plane: points satisfy
        // x = 1.5 + cos u, y = sin u; distance to the Z axis is √(3.25 + 3 cos u),
        // equal to 1 when cos u = −0.75 → two points.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &surface, (0.0, 2.0 * PI), (0.0, 2.0 * PI, -10.0, 10.0)).unwrap();
        assert_eq!(res.nb_points(), 2, "points: {:?}", res.points());
        for i in 0..2 {
            let p = res.point(i).pnt();
            let rho = (p.x() * p.x() + p.y() * p.y()).sqrt();
            assert!((rho - 1.0).abs() < 1e-5, "p {:?}", p);
        }
    }

    #[test]
    fn trimmed_circle_range_filters_points() {
        let sphere = sphere_at(GpPnt::zero(), 1.0);
        // Circle radius 1 centred (1.5, 0, 0) in the XY plane crosses the unit
        // sphere at u = acos(−0.75) ≈ 2.42 and u = 2π − 2.42 ≈ 3.86. Restricting
        // to [2.0, 3.0] keeps only the first.
        let ax2 = GpAx2::new(GpPnt::new(1.5, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let circle = GeomCircle::new(GpCirc::new(ax2, 1.0));
        let res = perform_curve_surface(&circle, &sphere, (2.0, 3.0), (0.0, 2.0 * PI, -PI / 2.0, PI / 2.0)).unwrap();
        assert_eq!(res.nb_points(), 1, "points: {:?}", res.points());
        let p = res.point(0).pnt();
        assert!((p.coord.modulus() - 1.0).abs() < 1e-9, "p {:?}", p);
    }
