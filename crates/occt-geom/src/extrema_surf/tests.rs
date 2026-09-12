use super::prelude::*;
use super::*;
    use crate::bspline_surface::{bspline_surface_uniform_knots, fit_surface_grid, GeomBSplineSurface};
    use crate::{GeomBSplineCurve, GeomLine, GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPnt, GpSphere as GpSphereT, GpTorus};

    const PI: f64 = std::f64::consts::PI;

    fn unit_sphere() -> GpSphere {
        GpSphereT::new(GpAx3::standard(), 1.0).unwrap()
    }

    #[test]
    fn point_sphere_extrema_min_and_max() {
        let sp = unit_sphere();
        let all = point_sphere_extrema(&sp, &GpPnt::new(3.0, 0.0, 0.0));
        assert_eq!(all.len(), 2, "min+max: {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 1.0).abs() < 1e-9, "closest {min:?}");
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.p2.x() + 1.0).abs() < 1e-9, "farthest {max:?}");
    }

    #[test]
    fn point_plane_extrema_distance() {
        let mut pl = GpPln::new(GpAx3::standard());
        pl.set_location(&GpPnt::new(0.0, 0.0, 5.0));
        let e = point_plane_extrema(&pl, &GpPnt::new(0.0, 0.0, 0.0));
        assert!((e.distance - 5.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.p2.z() - 5.0).abs() < 1e-9, "closest z {}", e.p2.z());
    }

    #[test]
    fn point_cylinder_extrema_min() {
        let cy = GpCylinder::new(GpAx3::standard(), 1.0).unwrap();
        let all = point_cylinder_extrema(&cy, &GpPnt::new(3.0, 0.0, 0.0));
        assert_eq!(all.len(), 2);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 1.0).abs() < 1e-9, "closest {min:?}");
    }

    #[test]
    fn point_cone_extrema_min() {
        // Cone placement at origin, RefRadius 1, semi-angle 45° → vertex at
        // (0,0,-1); in the XZ plane the generatrix is ρ = 1 + z. Point
        // (0.5, 0, 0.1): closest point at ρ = 0.8, z = -0.2, distance² = 0.18.
        let co = GpCone::new(GpAx3::standard(), 1.0, PI / 4.0).unwrap();
        let all = point_cone_extrema(&co, &GpPnt::new(0.5, 0.0, 0.1));
        assert_eq!(all.len(), 2, "cone extrema {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - (0.18f64).sqrt()).abs() < 1e-7, "min {}", min.distance);
        assert!((min.p2.x() - 0.8).abs() < 1e-6 && (min.p2.z() - (-0.2)).abs() < 1e-6, "closest {min:?}");
    }

    #[test]
    fn point_torus_extrema_min() {
        // Major 3, minor 1; point (6,0,0): torus crosses the X axis at
        // x = 4, 2, -2, -4 → min distance 2.
        let to = GpTorus::new(GpAx3::standard(), 3.0, 1.0).unwrap();
        let all = point_torus_extrema(&to, &GpPnt::new(6.0, 0.0, 0.0));
        assert_eq!(all.len(), 4, "torus extrema {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p2.x() - 4.0).abs() < 1e-9, "closest {min:?}");
    }

    #[test]
    fn sphere_parameters_roundtrip() {
        let sp = unit_sphere();
        for (u, v) in [(0.0, 0.0), (PI / 2.0, PI / 4.0), (PI, -PI / 3.0), (3.0 * PI / 2.0, 0.2)] {
            let q = slib::sphere_value(&sp, u, v);
            let (u2, v2) = sphere_parameters(&sp, &q);
            let q2 = slib::sphere_value(&sp, u2, v2);
            assert!(q.distance(&q2) < 1e-9, "params roundtrip (u={u},v={v})");
        }
    }

    #[test]
    fn curve_surface_line_sphere_min() {
        let line = GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 3.0, 0.0),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ));
        let sphere = GeomSphere::new(unit_sphere());
        let e = curve_surface_extrema(&line, &sphere);
        assert!((e.distance - 2.0).abs() < 1e-9, "min {}", e.distance);
        assert!((e.p1.y() - 3.0).abs() < 1e-9 && e.p1.x().abs() < 1e-9, "curve point {e:?}");
        assert!((e.p2.x() - 0.0).abs() < 1e-9 && (e.p2.y() - 1.0).abs() < 1e-9, "surf point {e:?}");
    }

    #[test]
    fn curve_surface_line_plane_intersect() {
        // Line (0,0,3) + t·(0,0,-1) pierces the plane z=0 at t=3 → distance 0.
        let line = GeomLine::new(GpLin::from_pnt_dir(
            GpPnt::new(0.0, 0.0, 3.0),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
        ));
        let pl = GpPln::new(GpAx3::standard());
        let e = curve_surface_extrema(&line, &GeomPlane::new(pl));
        assert!(e.distance.abs() < 1e-9, "line meets plane: {}", e.distance);
    }

    fn paraboloid() -> GeomBSplineSurface {
        let (nu, nv) = (3, 3);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v)
                    })
                    .collect()
            })
            .collect();
        fit_surface_grid(&points, 2, 2).unwrap()
    }

    #[test]
    fn newton_path_bspline_paraboloid_min() {
        // S(u,v) = (u, v, u²+v²) exactly (degree-2 Bernstein reproduces the
        // quadratic). Point P = (0.5, 0.5, -1). The closest point solves
        // F = 0; by symmetry u = v = t with 4t³ + 3t - 0.5 = 0 → t ≈ 0.16071.
        let s = paraboloid();
        let p = GpPnt::new(0.5, 0.5, -1.0);
        let all = point_surface_newton_all(&s, &p);
        assert!(!all.is_empty(), "no extrema found");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        // Independent reference: dense scan over the same domain.
        let mut ref_min = f64::INFINITY;
        for i in 0..=200 {
            for j in 0..=200 {
                let u = i as f64 / 200.0;
                let v = j as f64 / 200.0;
                let d = s.d0(u, v).distance(&p);
                if d < ref_min {
                    ref_min = d;
                }
            }
        }
        assert!(
            (min.distance - ref_min).abs() < 1e-4,
            "newton {} vs reference {}",
            min.distance,
            ref_min
        );
        // Orthogonality condition holds at the closest point.
        let f = ps_f(&s, &p, min.u1, min.v2.unwrap());
        assert!(
            f[0].abs() + f[1].abs() < 1e-6,
            "orthogonality F={f:?} at u={},v={}",
            min.u1,
            min.v2.unwrap()
        );
        // Local min: neighbours are farther.
        let (u, v) = (min.u1, min.v2.unwrap());
        let eps = 1e-3;
        let du = p.distance(&s.d0(u + eps, v));
        let dv = p.distance(&s.d0(u, v + eps));
        assert!(min.distance <= du + 1e-9 && min.distance <= dv + 1e-9, "not a local min");
    }

    #[test]
    fn point_surface_extrema_dispatch_sphere() {
        let sphere = GeomSphere::new(unit_sphere());
        let e = point_surface_extrema(&sphere, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((e.distance - 2.0).abs() < 1e-9, "min {}", e.distance);
    }

    #[test]
    fn classification_plane_and_sphere() {
        let sphere = GeomSphere::new(unit_sphere());
        assert!(classify_sphere(&sphere).is_some());
        assert!(classify_plane(&sphere).is_none());
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        assert!(classify_plane(&plane).is_some());
        assert!(classify_sphere(&plane).is_none());
    }

    #[test]
    fn curve_surface_newton_parabola_plane_min() {
        // Curve c(t) = (t, 0, 1+t²) over [0,1] (degree-2 B-spline), plane
        // z = 0. Minimum distance 1 at t = 0, with orthogonality F ≈ 0.
        // Degree-2 B-spline through control points (0,0,1), (0.5,0,1),
        // (1,0,2) reproduces c(t) = (t, 0, 1+t²) exactly.
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 1.0),
                GpPnt::new(0.5, 0.0, 1.0),
                GpPnt::new(1.0, 0.0, 2.0),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let s = GeomBSplineSurface::new(
            vec![
                vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)],
                vec![GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)],
            ],
            bspline_surface_uniform_knots(2, 2, 1, 1).0,
            bspline_surface_uniform_knots(2, 2, 1, 1).1,
            1,
            1,
        )
        .unwrap();
        let all = curve_surface_newton_all(&c, &s);
        assert!(!all.is_empty(), "no curve-surface extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 1.0).abs() < 1e-6, "min {}", min.distance);
        let f = cs_f(&c, &s, min.u1, min.u2, min.v2.unwrap()).unwrap();
        assert!(
            f[0].abs() + f[1].abs() + f[2].abs() < 1e-6,
            "orthogonality F={f:?}"
        );
    }
