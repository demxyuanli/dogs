use super::prelude::*;
use super::*;
    use crate::gp::dir::DirAxis;
    use std::f64::consts::PI;

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} != {b} within {tol}");
    }

    fn z_axis() -> GpAx1 {
        GpAx1::new(GpPnt::zero(), GpDir::from_axis(DirAxis::Z))
    }

    #[test]
    fn h_operator_point_above_origin() {
        let g = GpPnt::new(0.0, 0.0, 2.0);
        let m = h_operator(&g, &GpPnt::zero(), 1.0);
        assert_approx(m.m[0][0], 4.0, 1e-12); // y²+z²
        assert_approx(m.m[1][1], 4.0, 1e-12); // x²+z²
        assert_approx(m.m[2][2], 0.0, 1e-12); // x²+y²
        assert_approx(m.m[0][1], 0.0, 1e-12);
        assert_approx(m.m[0][2], 0.0, 1e-12);
    }

    #[test]
    fn empty_system_is_zero() {
        let s = GProps::new();
        assert_approx(s.mass(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().x(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().y(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().z(), 0.0, 1e-30);
        for r in 0..3 {
            for c in 0..3 {
                assert_approx(s.matrix_of_inertia().m[r][c], 0.0, 1e-30);
            }
        }
    }

    #[test]
    fn density_must_be_positive() {
        let mut s = GProps::new();
        let p = GProps::new();
        assert!(s.add(&p, 0.0).is_err());
        assert!(s.add(&p, -1.0).is_err());
        assert!(s.add(&p, 1e-8).is_err());
        assert!(s.add(&p, f64::NAN).is_err());
        assert!(s.add(&p, 1.0).is_ok());
    }

    #[test]
    fn point_mass_props() {
        let mut s = GProps::new();
        s.add_point_mass(&GpPnt::new(1.0, 0.0, 0.0), 2.0).unwrap();
        assert_approx(s.mass(), 2.0, 1e-12);
        assert_approx(s.centre_of_mass().x(), 1.0, 1e-12);
        // Moment about the z axis through the origin: m·1² = 2.
        assert_approx(s.moment_of_inertia(&z_axis()), 2.0, 1e-12);
        assert_approx(s.radius_of_gyration(&z_axis()), 1.0, 1e-12);
        // Static moments: COM·mass = (2, 0, 0).
        let (ix, iy, iz) = s.static_moments();
        assert_approx(ix, 2.0, 1e-12);
        assert_approx(iy, 0.0, 1e-12);
        assert_approx(iz, 0.0, 1e-12);
    }

    #[test]
    fn point_cloud_cube_corners() {
        // 8 unit masses at the corners of [0,1]³: mass 8, COM at the centre,
        // each principal moment about the COM = Σ(y²+z²) = 8·(0.25+0.25) = 4.
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.),
            GpPnt::new(1., 0., 1.),
            GpPnt::new(1., 1., 1.),
            GpPnt::new(0., 1., 1.),
        ];
        let mut pg = PGProps::new();
        for p in &pts {
            pg.add_point(p).unwrap();
        }
        assert_approx(pg.mass(), 8.0, 1e-12);
        let c = pg.centre_of_mass();
        assert_approx(c.x(), 0.5, 1e-12);
        assert_approx(c.y(), 0.5, 1e-12);
        assert_approx(c.z(), 0.5, 1e-12);
        let pp = pg.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 4.0, 1e-9);
        assert_approx(i2, 4.0, 1e-9);
        assert_approx(i3, 4.0, 1e-9);
        assert!(pp.has_symmetry_point());
        assert!(pp.has_symmetry_axis());
        let (r1, _, _) = pp.radius_of_gyration();
        assert_approx(r1, 0.5f64.sqrt(), 1e-9);
    }

    #[test]
    fn solid_unit_cube() {
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 1., 1.), 1.0).unwrap();
        assert_approx(v.mass(), 1.0, 1e-12);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 0.5, 1e-9);
        assert_approx(c.y(), 0.5, 1e-9);
        assert_approx(c.z(), 0.5, 1e-9);
        // Solid cube: I = m/12·(1+1) = 1/6 along each axis.
        let pp = v.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 1.0 / 6.0, 1e-9);
        assert_approx(i2, 1.0 / 6.0, 1e-9);
        assert_approx(i3, 1.0 / 6.0, 1e-9);
        assert!(pp.has_symmetry_point());
        let (r1, _, _) = pp.radius_of_gyration();
        assert_approx(r1, (1.0f64 / 6.0).sqrt(), 1e-9);
    }

    #[test]
    fn solid_box_2x3x4() {
        // Density 1, box 2×3×4: volume 24, COM at the centre.
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 3., 4.), 1.0).unwrap();
        assert_approx(v.mass(), 24.0, 1e-9);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 1.0, 1e-9);
        assert_approx(c.y(), 1.5, 1e-9);
        assert_approx(c.z(), 2.0, 1e-9);
        // Analytic box inertia about COM:
        //   Ixx = m/12·(dy²+dz²) = 24/12·(9+16) = 50
        //   Iyy = m/12·(dx²+dz²) = 24/12·(4+16) = 40
        //   Izz = m/12·(dx²+dy²) = 24/12·(4+9)  = 26
        let m = v.matrix_of_inertia();
        assert_approx(m.m[0][0], 50.0, 1e-9);
        assert_approx(m.m[1][1], 40.0, 1e-9);
        assert_approx(m.m[2][2], 26.0, 1e-9);
        assert_approx(m.m[0][1], 0.0, 1e-9);
        assert_approx(m.m[0][2], 0.0, 1e-9);
        assert_approx(m.m[1][2], 0.0, 1e-9);
        let pp = v.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 50.0, 1e-9);
        assert_approx(i2, 40.0, 1e-9);
        assert_approx(i3, 26.0, 1e-9);
        assert!(!pp.has_symmetry_point());
        assert!(!pp.has_symmetry_axis());
        // Principal axes of the axis-aligned box are the coordinate axes.
        let a = pp.first_axis_of_inertia();
        assert!(a.x().abs() > 0.9);
    }

    #[test]
    fn segment_properties() {
        let mut c = CelGProps::new();
        c.add_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 0., 0.), 1.0).unwrap();
        assert_approx(c.mass(), 2.0, 1e-12);
        assert_approx(c.centre_of_mass().x(), 1.0, 1e-12);
        // About COM: Ixx = 0, Iyy = Izz = m·L²/12 = 2·4/12 = 2/3.
        let m = c.matrix_of_inertia();
        assert_approx(m.m[0][0], 0.0, 1e-12);
        assert_approx(m.m[1][1], 2.0 / 3.0, 1e-12);
        assert_approx(m.m[2][2], 2.0 / 3.0, 1e-12);
    }

    #[test]
    fn surface_disk() {
        // Triangulate a radius-1 disk from its centre: area → π, Izz → π/2.
        let r = 1.0f64;
        let n = 256usize;
        let mut s = SelGProps::new();
        let center = GpPnt::new(0., 0., 0.);
        for i in 0..n {
            let t1 = 2.0 * PI * i as f64 / n as f64;
            let t2 = 2.0 * PI * (i + 1) as f64 / n as f64;
            let p1 = GpPnt::new(r * t1.cos(), r * t1.sin(), 0.);
            let p2 = GpPnt::new(r * t2.cos(), r * t2.sin(), 0.);
            s.add_triangle(&center, &p1, &p2, 1.0).unwrap();
        }
        assert_approx(s.mass(), PI * r * r, 0.01);
        let i = s.matrix_of_inertia();
        assert_approx(i.m[2][2], PI * r * r * r * r / 2.0, 0.02);
        assert_approx(s.centre_of_mass().x(), 0.0, 1e-9);
        assert_approx(s.centre_of_mass().y(), 0.0, 1e-9);
    }

    #[test]
    fn surface_sphere_shell() {
        // UV-grid triangulation of a unit sphere surface: area → 4π,
        // Izz → ∫(x²+y²)dA = (8/3)π (isotropic shell, about the centre).
        let r = 1.0f64;
        let (nu, nv) = (48usize, 32usize);
        let mut s = SelGProps::new();
        let pt = |u: f64, v: f64| GpPnt::new(r * u.cos() * v.sin(), r * u.sin() * v.sin(), r * v.cos());
        for i in 0..nu {
            let u0 = 2.0 * PI * i as f64 / nu as f64;
            let u1 = 2.0 * PI * (i + 1) as f64 / nu as f64;
            for j in 0..nv {
                let v0 = PI * j as f64 / nv as f64;
                let v1 = PI * (j + 1) as f64 / nv as f64;
                s.add_triangle(&pt(u0, v0), &pt(u1, v0), &pt(u0, v1), 1.0).unwrap();
                s.add_triangle(&pt(u1, v0), &pt(u1, v1), &pt(u0, v1), 1.0).unwrap();
            }
        }
        assert_approx(s.mass(), 4.0 * PI * r * r, 0.5);
        let i = s.matrix_of_inertia();
        assert_approx(i.m[2][2], (8.0 / 3.0) * PI * r.powi(4), 0.5);
    }

    #[test]
    fn solid_cylinder() {
        // Decompose a radius-1, height-2 cylinder into 48 triangular prisms,
        // each split into 3 tetrahedra. Volume → πr²h, Izz → ½mr².
        let r = 1.0f64;
        let h = 2.0f64;
        let n = 48usize;
        let z0 = -h / 2.0;
        let z1 = h / 2.0;
        let mut v = VelGProps::new();
        for i in 0..n {
            let t1 = 2.0 * PI * i as f64 / n as f64;
            let t2 = 2.0 * PI * (i + 1) as f64 / n as f64;
            let ob = GpPnt::new(0.0, 0.0, z0);
            let ot = GpPnt::new(0.0, 0.0, z1);
            let b1 = GpPnt::new(r * t1.cos(), r * t1.sin(), z0);
            let b2 = GpPnt::new(r * t2.cos(), r * t2.sin(), z0);
            let tp1 = GpPnt::new(r * t1.cos(), r * t1.sin(), z1);
            let tp2 = GpPnt::new(r * t2.cos(), r * t2.sin(), z1);
            v.add_tetrahedron(&ob, &b1, &b2, &tp2, 1.0).unwrap();
            v.add_tetrahedron(&ob, &b1, &tp2, &tp1, 1.0).unwrap();
            v.add_tetrahedron(&ob, &ot, &tp1, &tp2, 1.0).unwrap();
        }
        let m = v.mass();
        assert_approx(m, PI * r * r * h, 0.05 * PI * r * r * h);
        let i = v.matrix_of_inertia();
        assert_approx(i.m[2][2], 0.5 * m * r * r, 0.05 * (0.5 * m * r * r));
        // A cylinder has an axis of symmetry.
        let pp = v.principal_properties().unwrap();
        assert!(pp.has_symmetry_axis());
        assert!(!pp.has_symmetry_point());
    }

    #[test]
    fn compose_with_density() {
        let mut pg = PGProps::new();
        pg.add_point(&GpPnt::new(1., 0., 0.)).unwrap();
        let mut sys = GProps::new();
        sys.add(&pg, 3.0).unwrap(); // density 3
        assert_approx(sys.mass(), 3.0, 1e-12);
        assert_approx(sys.centre_of_mass().x(), 1.0, 1e-12);
    }

    #[test]
    fn composite_two_locations() {
        // A unit mass at the origin and another at (10,0,0): COM at (5,0,0),
        // inertia about COM: each contributes m·5² to Iyy and Izz → 50 each.
        let mut a = GProps::new();
        a.add_point_mass(&GpPnt::new(0., 0., 0.), 1.0).unwrap();
        let mut b = GProps::new_at(GpPnt::new(10., 0., 0.));
        b.add_point_mass(&GpPnt::new(10., 0., 0.), 1.0).unwrap();
        let mut sys = GProps::new();
        sys.add(&a, 1.0).unwrap();
        sys.add(&b, 1.0).unwrap();
        assert_approx(sys.mass(), 2.0, 1e-12);
        assert_approx(sys.centre_of_mass().x(), 5.0, 1e-9);
        assert_approx(sys.centre_of_mass().y(), 0.0, 1e-9);
        let m = sys.matrix_of_inertia();
        assert_approx(m.m[1][1], 50.0, 1e-9);
        assert_approx(m.m[2][2], 50.0, 1e-9);
        assert_approx(m.m[0][0], 0.0, 1e-9);
    }

    #[test]
    fn principal_axes_orthonormal() {
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 3., 4.), 1.0).unwrap();
        let pp = v.principal_properties().unwrap();
        let (a, b, c) = (
            pp.first_axis_of_inertia(),
            pp.second_axis_of_inertia(),
            pp.third_axis_of_inertia(),
        );
        assert_approx(a.dot(&b), 0.0, 1e-9);
        assert_approx(a.dot(&c), 0.0, 1e-9);
        assert_approx(b.dot(&c), 0.0, 1e-9);
        assert_approx(a.magnitude(), 1.0, 1e-9);
        assert_approx(b.magnitude(), 1.0, 1e-9);
        assert_approx(c.magnitude(), 1.0, 1e-9);
    }

    #[test]
    fn barycentre_helpers() {
        let pts = [GpPnt::new(0., 0., 0.), GpPnt::new(2., 0., 0.), GpPnt::new(0., 2., 0.)];
        let b = PGProps::barycentre(&pts);
        assert_approx(b.x(), 2.0 / 3.0, 1e-12);
        assert_approx(b.y(), 2.0 / 3.0, 1e-12);
        let (mass, g) = PGProps::weighted_barycentre(&pts, &[1.0, 2.0, 1.0]).unwrap();
        assert_approx(mass, 4.0, 1e-12);
        assert_approx(g.x(), 1.0, 1e-12); // (0·1 + 2·2 + 0·1)/4
        assert_approx(g.y(), 0.5, 1e-12); // (0·1 + 0·2 + 2·1)/4
        assert!(PGProps::weighted_barycentre(&pts, &[1.0, 0.0, 1.0]).is_err());
    }

    #[test]
    fn plane_equation_from_points() {
        let eq = PEquation::from_points(
            &GpPnt::new(0., 0., 0.),
            &GpPnt::new(1., 0., 0.),
            &GpPnt::new(0., 1., 0.),
        )
        .unwrap();
        assert_approx(eq.a, 0.0, 1e-12);
        assert_approx(eq.b, 0.0, 1e-12);
        assert!(eq.c.abs() > 0.9);
        let p = GpPnt::new(1.0, 2.0, 5.0);
        assert_approx(eq.distance_to(&p), 5.0, 1e-12);
        let proj = eq.project(&p);
        assert_approx(proj.x(), 1.0, 1e-9);
        assert_approx(proj.y(), 2.0, 1e-9);
        assert_approx(proj.z(), 0.0, 1e-9);
        // Collinear points are rejected.
        assert!(
            PEquation::from_points(
                &GpPnt::new(0., 0., 0.),
                &GpPnt::new(1., 0., 0.),
                &GpPnt::new(2., 0., 0.),
            )
            .is_err()
        );
    }

    #[test]
    fn plane_equation_from_point_normal() {
        let eq =
            PEquation::from_point_normal(&GpPnt::new(0., 0., 3.), &GpVec::new(0., 0., 2.)).unwrap();
        assert_approx(eq.signed_distance(&GpPnt::new(0., 0., 5.)), 2.0, 1e-12);
        assert_approx(eq.distance_to(&GpPnt::new(1., 1., 3.)), 0.0, 1e-12);
        let n = eq.normal_unit();
        assert_approx(n.z(), 1.0, 1e-12);
        // Zero normal is rejected.
        assert!(PEquation::from_point_normal(&GpPnt::zero(), &GpVec::zero()).is_err());
    }

    #[test]
    fn plane_equation_to_plane() {
        let eq =
            PEquation::from_point_normal(&GpPnt::new(1., 2., 3.), &GpVec::new(1., 1., 1.)).unwrap();
        let pln = eq.to_plane();
        // Every point of the plane satisfies the equation.
        let loc = pln.location();
        assert_approx(eq.distance_to(&loc), 0.0, 1e-9);
        let axis = pln.axis();
        let nu = eq.normal_unit();
        assert!(axis.direction().xyz().dot(&nu.coord) > 0.999);
    }

    #[test]
    fn value_type_enum() {
        // Sanity: the enum carries the OCCT GProp_ValueType variants.
        let _ = [
            ValueType::Mass,
            ValueType::CenterMassX,
            ValueType::CenterMassY,
            ValueType::CenterMassZ,
            ValueType::InertiaXx,
            ValueType::InertiaYy,
            ValueType::InertiaZz,
            ValueType::InertiaXy,
            ValueType::InertiaXz,
            ValueType::InertiaYz,
            ValueType::Unknown,
        ];
        assert_eq!(ValueType::Mass as u8, 0);
        assert_eq!(ValueType::Unknown as u8, 10);
    }

    #[test]
    fn octahedron_volume() {
        // Regular octahedron with vertices (±1,0,0),(0,±1,0),(0,0,±1):
        // volume = 4/3, centroid at origin.
        let top = GpPnt::new(0., 0., 1.);
        let bottom = GpPnt::new(0., 0., -1.);
        let equator = [
            GpPnt::new(1., 0., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(-1., 0., 0.),
            GpPnt::new(0., -1., 0.),
        ];
        let mut v = VelGProps::new();
        v.add_octahedron(&top, &bottom, &equator, 1.0).unwrap();
        assert_approx(v.mass(), 4.0 / 3.0, 1e-9);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 0.0, 1e-9);
        assert_approx(c.y(), 0.0, 1e-9);
        assert_approx(c.z(), 0.0, 1e-9);
        // The regular octahedron has a point of symmetry.
        assert!(v.principal_properties().unwrap().has_symmetry_point());
    }
