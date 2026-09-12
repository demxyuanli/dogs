use super::prelude::*;
use super::*;
    use occt_core::gp::GpAx3;
    use occt_core::precision::ANGULAR;

    const PI: f64 = std::f64::consts::PI;

    fn plane(origin: GpPnt, normal: GpDir) -> GpPln {
        GpPln::new(
            GpAx3::new(origin, normal, &perp_x_dir(&normal)).unwrap_or_default(),
        )
    }

    fn sphere_at(origin: GpPnt, r: f64) -> GpSphere {
        let ax3 = GpAx3::new(origin, GpDir::from_axis(DirAxis::Z), &GpDir::from_axis(DirAxis::X))
            .unwrap_or_default();
        GpSphere::new(ax3, r).unwrap()
    }

    fn on_plane(p: &GpPnt, pl: &GpPln, tol: f64) -> bool {
        let (a, b, c, d) = plane_coeffs(pl);
        (a * p.x() + b * p.y() + c * p.z() + d).abs() < tol
    }

    fn on_cone(p: &GpPnt, cone: &GpCone, tol: f64) -> bool {
        let apex = cone.apex();
        let axis = *cone.axis().direction();
        let v = GpVec::from_pnts(&apex, p);
        let h = v.dot(&GpVec::from_xyz(axis.xyz()));
        let perp = v.coord.subtracted(&axis.xyz().multiplied(h));
        let rho = perp.modulus();
        (rho - h.abs() * cone.semi_angle().tan()).abs() < tol
    }

    fn on_cylinder(p: &GpPnt, cyl: &GpCylinder, tol: f64) -> bool {
        let axis = *cyl.axis().direction();
        let v = GpVec::from_pnts(cyl.axis().location(), p);
        let h = v.dot(&GpVec::from_xyz(axis.xyz()));
        let perp = v.coord.subtracted(&axis.xyz().multiplied(h));
        (perp.modulus() - cyl.radius()).abs() < tol
    }

    fn on_sphere(p: &GpPnt, sph: &GpSphere, tol: f64) -> bool {
        (p.distance(&sph.location()) - sph.radius()).abs() < tol
    }

    fn cylinder_axis(center: GpPnt, axis: GpDir, r: f64) -> GpCylinder {
        let x = perp_x_dir(&axis);
        GpCylinder::new(GpAx3::new(center, axis, &x).unwrap_or_default(), r).unwrap()
    }

    /// A cone whose geometric apex is `apex` (radius 0 anchor).
    fn cone_apex_at(apex: GpPnt, axis: GpDir, semi: f64) -> GpCone {
        let x = perp_x_dir(&axis);
        GpCone::new(GpAx3::new(apex, axis, &x).unwrap_or_default(), 0.0, semi).unwrap()
    }

    #[test]
    fn three_planes_known_point() {
        // x=1, y=2, z=3 meet at (1,2,3).
        let p1 = plane(GpPnt::new(1.0, 0.0, 0.0), GpDir::from_axis(DirAxis::X));
        let p2 = plane(GpPnt::new(0.0, 2.0, 0.0), GpDir::from_axis(DirAxis::Y));
        let p3 = plane(GpPnt::new(0.0, 0.0, 3.0), GpDir::from_axis(DirAxis::Z));
        match three_planes_intersect(&p1, &p2, &p3).unwrap() {
            Intersect3Pln::Point(p) => {
                assert!((p.x() - 1.0).abs() < 1e-9);
                assert!((p.y() - 2.0).abs() < 1e-9);
                assert!((p.z() - 3.0).abs() < 1e-9, "got {p:?}");
            }
            other => panic!("expected point, got {other:?}"),
        }
    }

    #[test]
    fn three_planes_parallel_no_result() {
        let p1 = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let p2 = plane(GpPnt::new(0.0, 0.0, 1.0), GpDir::from_axis(DirAxis::Z));
        let p3 = plane(GpPnt::new(0.0, 0.0, 2.0), GpDir::from_axis(DirAxis::Z));
        assert_eq!(three_planes_intersect(&p1, &p2, &p3).unwrap(), Intersect3Pln::NoResult);
    }

    #[test]
    fn plane_sphere_circle_z0() {
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let s = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        match quadric_quadric_plane_sphere(&pl, &s) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - 1.0).abs() < 1e-9, "radius {}", c.radius());
                // Every sampled point lies in z=0 at distance 1 from the origin.
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9), "{p:?} not in plane");
                    assert!((p.distance(&GpPnt::new(0., 0., 0.)) - 1.0).abs() < 1e-9);
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_sphere_no_intersection() {
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let s = sphere_at(GpPnt::new(0.0, 0.0, 5.0), 1.0);
        assert_eq!(quadric_quadric_plane_sphere(&pl, &s), QuadricIntersection::None);
    }

    #[test]
    fn sphere_sphere_circle() {
        // Two unit spheres, centers 1 apart: circle radius √(1 − (d/2)²) = √0.75
        // at the midpoint (0.5, 0, 0).
        let s1 = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let s2 = sphere_at(GpPnt::new(1.0, 0.0, 0.0), 1.0);
        match quadric_quadric_sphere_sphere(&s1, &s2, 1e-7) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0.5, 0.0, 0.0)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - (0.75f64).sqrt()).abs() < 1e-9, "radius {}", c.radius());
                // Every point on the circle is on both spheres.
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!((p.distance(&GpPnt::new(0., 0., 0.)) - 1.0).abs() < 1e-9);
                    assert!((p.distance(&GpPnt::new(1., 0., 0.)) - 1.0).abs() < 1e-9);
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn sphere_sphere_disjoint() {
        let s1 = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        let s2 = sphere_at(GpPnt::new(5.0, 0.0, 0.0), 1.0);
        assert_eq!(quadric_quadric_sphere_sphere(&s1, &s2, 1e-7), QuadricIntersection::None);
    }

    #[test]
    fn plane_cylinder_circle_perpendicular_axis() {
        // Cylinder along Z at origin r=1, plane z=0 → circle radius 1 at origin.
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(DirAxis::Z));
        let cyl = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        match quadric_quadric_plane_cylinder(&pl, &cyl, ANGULAR, CONFUSION) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                assert!((c.radius() - 1.0).abs() < 1e-9);
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9));
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_cylinder_two_lines_parallel_axis() {
        // Cylinder along X, base at (0,1,0.5), r=1; plane z=0 → two lines at
        // y = 1 ± √(1 − 0.5²) in z=0.
        let xdir = GpDir::from_axis(DirAxis::X);
        let zdir = GpDir::from_axis(DirAxis::Z);
        let ax3 = GpAx3::new(GpPnt::new(0.0, 1.0, 0.5), xdir, &zdir).unwrap();
        let cyl = GpCylinder::new(ax3, 1.0).unwrap();
        let pl = plane(GpPnt::new(0.0, 0.0, 0.0), zdir);
        match quadric_quadric_plane_cylinder(&pl, &cyl, ANGULAR, CONFUSION) {
            QuadricIntersection::TwoLines(l1, l2) => {
                let h = (0.75f64).sqrt();
                let p1 = clib::line_value(&l1, 0.0);
                let p2 = clib::line_value(&l2, 0.0);
                for p in [p1, p2] {
                    assert!(on_plane(&p, &pl, 1e-9), "{p:?} not in plane");
                    assert!((p.y() - (1.0 + h)).abs() < 1e-9 || (p.y() - (1.0 - h)).abs() < 1e-9, "y {}", p.y());
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn plane_cone_circle_perpendicular_axis() {
        // Cone apex at origin, axis +Z, semi-angle 30°; plane z=1 → circle of
        // radius tan(30°) at (0,0,1).
        let zdir = GpDir::from_axis(DirAxis::Z);
        let h = 1.0;
        let radius = h * (PI / 6.0).tan();
        let loc = GpPnt::new(0.0, 0.0, -h);
        let ax3 = GpAx3::new(loc, zdir, &GpDir::from_axis(DirAxis::X)).unwrap();
        let cone = GpCone::new(ax3, radius, PI / 6.0).unwrap();
        assert!(cone.apex().distance(&GpPnt::new(0., 0., 0.)) < 1e-9, "apex {:?}", cone.apex());
        let pl = plane(GpPnt::new(0.0, 0.0, 1.0), zdir);
        match quadric_quadric_plane_cone(&pl, &cone, ANGULAR, CONFUSION) {
            QuadricIntersection::Circle(c) => {
                assert!(c.location().distance(&GpPnt::new(0., 0., 1.)) < 1e-9, "center {:?}", c.location());
                assert!((c.radius() - radius).abs() < 1e-9, "radius {} vs {}", c.radius(), radius);
                for k in 0..8 {
                    let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                    assert!(on_plane(&p, &pl, 1e-9));
                    assert!(on_cone(&p, &cone, 1e-6), "{p:?} not on cone");
                }
            }
            other => panic!("expected circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_cone_ellipse_lies_on_quadrics() {
        // Oblique plane (not through the apex, not perpendicular, not parallel
        // to a generatrix) cutting a cone → ellipse; verify membership by
        // sampling.
        let cone = GpCone::new(GpAx3::standard(), 1.0, PI / 6.0).unwrap();
        let pl = plane(GpPnt::new(0.0, 0.0, 2.0), GpDir::new(0.0, 0.3, 0.954).unwrap());
        match quadric_quadric_plane_cone(&pl, &cone, ANGULAR, CONFUSION) {
            QuadricIntersection::Ellipse(e) => {
                for k in 0..16 {
                    let p = clib::ellipse_value(&e, 2.0 * PI * k as f64 / 16.0);
                    assert!(on_plane(&p, &pl, 1e-6), "{p:?} not in plane");
                    assert!(on_cone(&p, &cone, 1e-5), "{p:?} not on cone");
                }
            }
            other => panic!("expected ellipse, got {other:?}"),
        }
    }

    #[test]
    fn line_torus_x_axis_four_points() {
        // Torus at origin R=3 r=1; line along X through origin cuts at x=±2, ±4.
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let line = GpLin::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::from_axis(DirAxis::X));
        let pts = line_torus_intersect(&line, &torus);
        assert_eq!(pts.len(), 4, "got {pts:?}");
        let expected = [-4.0, -2.0, 2.0, 4.0];
        let mut got: Vec<f64> = pts.iter().map(|p| p.x()).collect();
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (g, e) in got.iter().zip(expected.iter()) {
            assert!((g - e).abs() < 1e-6, "got {g} expected {e}: {pts:?}");
        }
    }

    #[test]
    fn line_torus_miss() {
        // Line along Z through the center misses (the minor circle is off-axis).
        let torus = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let line = GpLin::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::from_axis(DirAxis::Z));
        assert!(line_torus_intersect(&line, &torus).is_empty());
    }

    #[test]
    fn cylinder_cylinder_parallel_two_lines() {
        // Two unit cylinders, parallel Z axes through (0,0,0) and (0.5,0,0):
        // base circles intersect at x=0.25, y=±√(1−0.25²).
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.5, 0.0, 0.0), z, 1.0);
        match quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7) {
            QuadricIntersection::TwoLines(l1, l2) => {
                let h = (0.9375f64).sqrt();
                for l in [l1, l2] {
                    let p = clib::line_value(&l, 0.0);
                    assert!((p.y().abs() - h).abs() < 1e-9, "y {}", p.y());
                    assert!((p.x() - 0.25).abs() < 1e-9, "x {}", p.x());
                    assert!(on_cylinder(&p, &c1, 1e-9), "{p:?} not on cyl1");
                    assert!(on_cylinder(&p, &c2, 1e-9), "{p:?} not on cyl2");
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_cylinder_nested_empty() {
        // Concentric cylinders of different radii never meet.
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 2.0);
        assert_eq!(
            quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7),
            QuadricIntersection::None
        );
    }

    #[test]
    fn cylinder_cylinder_intersecting_ellipses() {
        // Equal unit cylinders, perpendicular axes through the origin: two
        // bisector-plane ellipses (x²+z²=1 ∧ x²+y²=1 → y=±z).
        let z = GpDir::from_axis(DirAxis::Z);
        let x = GpDir::from_axis(DirAxis::X);
        let c1 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 1.0);
        let c2 = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), x, 1.0);
        match quadric_quadric_cylinder_cylinder(&c1, &c2, 1e-7) {
            QuadricIntersection::TwoEllipses(e1, e2) => {
                assert!(e1.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                assert!(e2.location().distance(&GpPnt::new(0., 0., 0.)) < 1e-9);
                for e in [e1, e2] {
                    for k in 0..16 {
                        let p = clib::ellipse_value(&e, 2.0 * PI * k as f64 / 16.0);
                        assert!(on_cylinder(&p, &c1, 1e-6), "{p:?} not on cyl1");
                        assert!(on_cylinder(&p, &c2, 1e-6), "{p:?} not on cyl2");
                    }
                }
            }
            other => panic!("expected two ellipses, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_sphere_two_circles() {
        // Unit sphere at the origin, Z-axis cylinder r=0.5 through it: circles
        // at z=±√(1−0.25)=±0.866, radius 0.5.
        let z = GpDir::from_axis(DirAxis::Z);
        let cyl = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 0.5);
        let sph = sphere_at(GpPnt::new(0.0, 0.0, 0.0), 1.0);
        match quadric_quadric_cylinder_sphere(&cyl, &sph, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    assert!((c.radius() - 0.5).abs() < 1e-9, "radius {}", c.radius());
                    assert!(c.location().x().abs() < 1e-9 && c.location().y().abs() < 1e-9);
                    assert!((c.location().z().abs() - (0.75f64).sqrt()).abs() < 1e-9);
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cylinder(&p, &cyl, 1e-6), "{p:?} not on cylinder");
                        assert!(on_sphere(&p, &sph, 1e-6), "{p:?} not on sphere");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn sphere_cone_two_circles() {
        // Cone apex at the origin, axis +Z, 30°; sphere center (0,0,3), r=2:
        // two circles at z≈3.40 (r≈1.96) and z≈1.10 (r≈0.64).
        let z = GpDir::from_axis(DirAxis::Z);
        let cone = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let sph = sphere_at(GpPnt::new(0.0, 0.0, 3.0), 2.0);
        match quadric_quadric_sphere_cone(&sph, &cone, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cone(&p, &cone, 1e-5), "{p:?} not on cone");
                        assert!(on_sphere(&p, &sph, 1e-5), "{p:?} not on sphere");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn cone_cone_same_axis_two_circles() {
        // Same axis +Z; cone1 apex origin 30°, cone2 apex (0,0,1) 45°.
        let z = GpDir::from_axis(DirAxis::Z);
        let c1 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let c2 = cone_apex_at(GpPnt::new(0.0, 0.0, 1.0), z, PI / 4.0);
        match quadric_quadric_cone_cone(&c1, &c2, ANGULAR, 1e-7) {
            QuadricIntersection::TwoCircles(a, b) => {
                for c in [a, b] {
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cone(&p, &c1, 1e-5), "{p:?} not on cone1");
                        assert!(on_cone(&p, &c2, 1e-5), "{p:?} not on cone2");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }

    #[test]
    fn cone_cone_common_apex_two_lines() {
        // Two 30° cones sharing the apex, axes differing by 40°: two generatrix
        // lines through the apex.
        let z = GpDir::from_axis(DirAxis::Z);
        let axis2 = GpDir::new((40.0f64).to_radians().sin(), 0.0, (40.0f64).to_radians().cos()).unwrap();
        let c1 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 6.0);
        let c2 = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), axis2, PI / 6.0);
        match quadric_quadric_cone_cone(&c1, &c2, ANGULAR, 1e-7) {
            QuadricIntersection::TwoLines(l1, l2) => {
                for l in [l1, l2] {
                    for t in [-2.0, -1.0, 1.0, 2.0] {
                        let p = clib::line_value(&l, t);
                        assert!(on_cone(&p, &c1, 1e-5), "{p:?} not on cone1");
                        assert!(on_cone(&p, &c2, 1e-5), "{p:?} not on cone2");
                    }
                }
            }
            other => panic!("expected two lines, got {other:?}"),
        }
    }

    #[test]
    fn cylinder_cone_same_axis_two_circles() {
        // Cylinder r=0.5 axis Z through origin; cone apex origin 45°: circles at
        // z=±0.5/tan(45°)=±0.5, radius 0.5.
        let z = GpDir::from_axis(DirAxis::Z);
        let cyl = cylinder_axis(GpPnt::new(0.0, 0.0, 0.0), z, 0.5);
        let cone = cone_apex_at(GpPnt::new(0.0, 0.0, 0.0), z, PI / 4.0);
        match quadric_quadric_cylinder_cone(&cyl, &cone, 1e-7) {
            QuadricIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    assert!((c.radius() - 0.5).abs() < 1e-9);
                    assert!(c.location().x().abs() < 1e-9 && c.location().y().abs() < 1e-9);
                    assert!((c.location().z().abs() - 0.5).abs() < 1e-9);
                    for k in 0..8 {
                        let p = clib::circle_value(&c, 2.0 * PI * k as f64 / 8.0);
                        assert!(on_cylinder(&p, &cyl, 1e-6), "{p:?} not on cylinder");
                        assert!(on_cone(&p, &cone, 1e-6), "{p:?} not on cone");
                    }
                }
            }
            other => panic!("expected two circles, got {other:?}"),
        }
    }
