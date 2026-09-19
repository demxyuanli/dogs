use super::prelude::*;
use super::*;
    use crate::curve::Curve;
    use crate::{GeomBSplineCurve, GeomCircle, GeomLine};
    use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPnt, GpVec};

    const PI: f64 = std::f64::consts::PI;

    fn line_on_x_axis() -> GpLin {
        GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap())
    }

    #[test]
    fn point_line_extrema_closest() {
        let e = point_line_extrema(&line_on_x_axis(), &GpPnt::new(3.0, 4.0, 0.0));
        assert!((e.distance - 4.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1 - 3.0).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 3.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9 && e.p2.z().abs() < 1e-9);
    }

    #[test]
    fn point_circle_min_and_max_analytic() {
        let c = GpCirc::new(GpAx2::standard(), 1.0);
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let all = circle_all(&c, &p, 0.0, 2.0 * PI);
        assert_eq!(all.len(), 2, "min+max: {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.u1).abs() < 1e-9, "min u {}", min.u1);
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.u1 - PI).abs() < 1e-9, "max u {}", max.u1);
    }

    #[test]
    fn point_ellipse_x_axis() {
        // a=2, b=1, point (3,0,0): closest at u=0 (point (2,0,0)), farthest at u=π.
        let e = GpElips::new(GpAx2::standard(), 2.0, 1.0);
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let all = ellipse_all(&e, &p, 0.0, 2.0 * PI);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 1.0).abs() < 1e-7, "min {}", min.distance);
        assert!((min.p2.x() - 2.0).abs() < 1e-6, "closest {min:?}");
        assert!((max.distance - 5.0).abs() < 1e-7, "max {}", max.distance);
        // The extremum condition must hold on the actual curve.
        for s in all {
            let q = clib::ellipse_value(&e, s.u1);
            let (_, d1) = clib::ellipse_d1(&e, s.u1);
            let v = GpVec::from_pnts(&p, &q);
            assert!(v.dot(&d1).abs() < 1e-6, "dF≠0 at u={}", s.u1);
        }
    }

    #[test]
    fn point_ellipse_y_axis() {
        // Point (0,4,0), a=2,b=1: closest at u=π/2 → (0,1,0), dist 3.
        // Parameterization is `ElCLib::EllipseValue` (`ElCLib.cxx:176-189`):
        // P = Loc + Major*cos(U)*XDir + Minor*sin(U)*YDir, so +YDir sits at
        // U = π/2 (the previous expectation 3π/2 encoded the mirrored
        // `-Minor*sin(U)` that `clib::ellipse_value` no longer uses).
        let e = GpElips::new(GpAx2::standard(), 2.0, 1.0);
        let p = GpPnt::new(0.0, 4.0, 0.0);
        let em = point_ellipse_extrema(&e, &p);
        assert!((em.distance - 3.0).abs() < 1e-7, "min {}", em.distance);
        assert!((em.p2.x()).abs() < 1e-6 && (em.p2.y() - 1.0).abs() < 1e-6, "closest {em:?}");
        assert!((em.u1 - PI / 2.0).abs() < 1e-6, "u {}", em.u1);
    }

    #[test]
    fn point_hyperbola_origin() {
        // a=b=1, point at origin: closest at vertex (1,0,0), dist 1, u=0.
        let h = GpHypr::new(GpAx2::standard(), 1.0, 1.0);
        let p = GpPnt::new(0.0, 0.0, 0.0);
        let e = point_hyperbola_extrema(&h, &p);
        assert!((e.distance - 1.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 1.0).abs() < 1e-9, "closest {e:?}");
    }

    #[test]
    fn point_parabola_side() {
        // F=1, C(u) = (u²/4, u). Point (4,0): cubic u³/4 − 2u = 0 → u=0, ±2√2.
        // Min distance √12 at u=±2√2.
        let pa = GpParab::new(GpAx2::standard(), 1.0);
        let p = GpPnt::new(4.0, 0.0, 0.0);
        let e = point_parabola_extrema(&pa, &p);
        assert!((e.distance - (12.0f64).sqrt()).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.u1.abs() - 2.0 * (2.0f64).sqrt()).abs() < 1e-6, "u {}", e.u1);
    }

    #[test]
    fn point_curve_extrema_line_dispatch() {
        let line = GeomLine::new(line_on_x_axis());
        let e = point_curve_extrema(&line, &GpPnt::new(3.0, 4.0, 0.0));
        assert!((e.distance - 4.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1 - 3.0).abs() < 1e-9, "u {}", e.u1);
    }

    #[test]
    fn point_curve_extrema_circle_dispatch_min_and_max() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let all = point_curve_extrema_all(&circle, &p);
        assert!(!all.is_empty(), "no extrema found");
        assert!((all[0].distance - 2.0).abs() < 1e-6, "min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 4.0).abs() < 1e-6, "max {}", all[all.len() - 1].distance);
        let e = point_curve_extrema(&circle, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min via wrapper {}", e.distance);
    }

    #[test]
    fn newton_path_bspline_min() {
        // A convex degree-2 B-spline arc (parabola-like). The Newton path must
        // satisfy the derivative condition and return a genuine local minimum.
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(2.0, 2.0, 0.0),
                GpPnt::new(4.0, 0.0, 0.0),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let p = GpPnt::new(2.0, 4.0, 0.0);
        let all = newton_point_curve_all(&c, &p);
        assert!(!all.is_empty());
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        // Derivative condition (C-P)·C′ ≈ 0 at the closest parameter.
        let (_, d1) = c.d1(min.u1);
        let v = GpVec::from_pnts(&p, &min.p2);
        assert!(v.dot(&d1).abs() < 1e-6, "dF={}", v.dot(&d1));
        // The returned distance must be a local minimum: any neighbour is farther.
        let eps = 1e-3;
        let d0 = p.distance(&c.d0(min.u1 - eps));
        let d2 = p.distance(&c.d0(min.u1 + eps));
        assert!(min.distance <= d0 + 1e-9 && min.distance <= d2 + 1e-9, "not a local min");
        // And point_curve_extrema delegates to it.
        let e = point_curve_extrema(&c, &p);
        assert!((e.distance - min.distance).abs() < 1e-6, "wrapper {}", e.distance);
    }

    #[test]
    fn polynomial_helpers_ground_truth() {
        assert_eq!(quadratic_roots(1.0, -5.0, 6.0), vec![2.0, 3.0]);
        let mut c = cubic_roots(1.0, -6.0, 11.0, -6.0);
        c.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(c.len(), 3);
        assert!((c[0] - 1.0).abs() < 1e-10 && (c[1] - 2.0).abs() < 1e-10 && (c[2] - 3.0).abs() < 1e-10);
        // quartic (x²-1)(x²-4) = x⁴ - 5x² + 4.
        let mut q = quartic_roots(1.0, 0.0, -5.0, 0.0, 4.0);
        q.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(q.len(), 4, "quartic roots {q:?}");
        for (got, want) in q.iter().zip([-2.0, -1.0, 1.0, 2.0]) {
            assert!((got - want).abs() < 1e-6, "root {got} vs {want}");
        }
    }

    #[test]
    fn trig_roots_sincos_ground_truth() {
        // sin(u) = 0 → d=0,e=0,f=1: roots 0, π, 2π in [0, 2π].
        let roots = trig_roots_sincos(0.0, 0.0, 1.0, 0.0, 2.0 * PI);
        assert_eq!(roots.len(), 3, "got {roots:?}");
        assert!((roots[0]).abs() < 1e-8);
        assert!((roots[1] - PI).abs() < 1e-8);
        assert!((roots[2] - 2.0 * PI).abs() < 1e-8);
        // cos(u) = 0 → u = π/2, 3π/2 (d=0, e=1, f=0)
        let roots = trig_roots_sincos(0.0, 1.0, 0.0, 0.0, 2.0 * PI);
        assert_eq!(roots.len(), 2, "got {roots:?}");
        assert!((roots[0] - PI / 2.0).abs() < 1e-8);
        assert!((roots[1] - 3.0 * PI / 2.0).abs() < 1e-8);
    }
