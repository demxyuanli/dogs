use super::prelude::*;
use super::*;
    use crate::{Geom2dCircle, Geom2dEllipse, Geom2dLine};
    use occt_core::gp::{GpAx22d, GpAx2d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d};

    fn line2d() -> Geom2dLine {
        Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()))
    }

    fn circle2d() -> Geom2dCircle {
        Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0))
    }

    #[test]
    fn point_line_min_distance() {
        let c = line2d();
        let e = point_curve_extrema2d(&c, &GpPnt2d::new(3.0, 4.0));
        assert!((e.distance - 4.0).abs() < 1e-6, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-5, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let c = circle2d();
        let p = GpPnt2d::new(3.0, 0.0);
        let e = point_curve_extrema2d(&c, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min {}", e.distance);
        let m = point_curve_max_extrema2d(&c, &p, 64);
        assert!((m.distance - 4.0).abs() < 1e-6, "max {}", m.distance);
    }

    #[test]
    fn curve_curve_circles_min() {
        // Two unit circles centers 3 apart → min distance 1.
        let c1 = circle2d();
        let c2 = Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0)));
        let es = curve_curve_extrema2d(&c1, &c2, 24);
        assert!(!es.is_empty());
        assert!((es[0].distance - 1.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    #[test]
    fn curve_curve_intersections_found() {
        // Circle at origin and horizontal line y=0 → 2 intersection points.
        let c1 = circle2d();
        let l = line2d();
        let pts = curve_curve_intersections2d(&c1, &l, 1e-3);
        assert_eq!(pts.len(), 2, "circle∩x-axis points {pts:?}");
        for p in &pts {
            assert!((p.x().abs() - 1.0).abs() < 0.05, "on circle x {}", p.x());
        }
    }

    #[test]
    fn point_polyline_closest() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(1.0, 3.0));
        assert!((e.distance - 3.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.p2.x() - 1.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9);
    }

    #[test]
    fn point_polyline_endpoint_t() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        // Point beyond the segment end clamps to the endpoint.
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(5.0, 1.0));
        assert!((e.p2.x() - 2.0).abs() < 1e-9, "clamp x {}", e.p2.x());
    }

    #[test]
    fn tangent_horizontal_line() {
        let c = line2d();
        let t = tangent2d(&c, 0.5);
        assert!((t.x() - 1.0).abs() < 1e-9 && t.y().abs() < 1e-9, "tangent {t:?}");
    }

    #[test]
    fn skew_lines_min_distance() {
        // Horizontal line and vertical line offset — skew in the plane means
        // they cross; min distance 0. Use parallel lines instead: y=0 and y=3.
        let l1 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let l2 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let es = curve_curve_extrema2d(&l1, &l2, 8);
        assert!((es[0].distance - 3.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    const PI: f64 = std::f64::consts::PI;

    #[test]
    fn analytic_line_extrema2d() {
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let e = point_line_extrema2d(&l, &GpPnt2d::new(3.0, 4.0));
        assert!((e.distance - 4.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1 - 3.0).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 3.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9);
    }

    #[test]
    fn analytic_circle_extrema2d_min_max() {
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let all = circle_all2d(&c, &p, 0.0, 2.0 * PI);
        assert_eq!(all.len(), 2, "min+max {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.u1).abs() < 1e-9, "min u {}", min.u1);
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.u1 - PI).abs() < 1e-9, "max u {}", max.u1);
    }

    #[test]
    fn analytic_ellipse_extrema2d() {
        // a=2, b=1, point (3,0): closest (2,0) at u=0 (dist 1), farthest (-2,0) at u=π.
        let e = GpElips2d::new(GpAx22d::standard(), 2.0, 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let em = point_ellipse_extrema2d(&e, &p);
        assert!((em.distance - 1.0).abs() < 1e-7, "min {}", em.distance);
        assert!((em.p2.x() - 2.0).abs() < 1e-6 && em.p2.y().abs() < 1e-6, "closest {em:?}");
        let all = ellipse_all2d(&e, &p, 0.0, 2.0 * PI);
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((max.distance - 5.0).abs() < 1e-7, "max {}", max.distance);
    }

    #[test]
    fn analytic_hyperbola_extrema2d_origin() {
        let h = GpHypr2d::new(GpAx22d::standard(), 1.0, 1.0);
        let e = point_hyperbola_extrema2d(&h, &GpPnt2d::new(0.0, 0.0));
        assert!((e.distance - 1.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 1.0).abs() < 1e-9, "closest {e:?}");
    }

    #[test]
    fn analytic_parabola_extrema2d() {
        // F=1, C(u) = (u²/4, u). Point (4,0): min √12 at u = ±2√2.
        let pa = GpParab2d::new(GpAx22d::standard(), 1.0);
        let e = point_parabola_extrema2d(&pa, &GpPnt2d::new(4.0, 0.0));
        assert!((e.distance - (12.0f64).sqrt()).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.u1.abs() - 2.0 * (2.0f64).sqrt()).abs() < 1e-6, "u {}", e.u1);
    }

    #[test]
    fn point_curve_extrema2d_all_line_and_circle() {
        let line = line2d();
        let all = point_curve_extrema2d_all(&line, &GpPnt2d::new(3.0, 4.0));
        assert!(!all.is_empty());
        assert!((all[0].distance - 4.0).abs() < 1e-9, "line min {}", all[0].distance);

        let circle = circle2d();
        let p = GpPnt2d::new(3.0, 0.0);
        let all = point_curve_extrema2d_all(&circle, &p);
        assert!((all[0].distance - 2.0).abs() < 1e-6, "circle min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 4.0).abs() < 1e-6, "circle max {}", all[all.len() - 1].distance);
    }

    #[test]
    fn newton_path_2d_bspline_min() {
        // Convex degree-2 arc: the Newton path must return a genuine minimum.
        let c = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let p = GpPnt2d::new(2.0, 4.0);
        let all = point_curve_extrema2d_all(&c, &p);
        assert!(!all.is_empty());
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let (_, d1) = c.d1(min.u1);
        let v = GpVec2d::new(min.p2.x() - p.x(), min.p2.y() - p.y());
        assert!(v.dot(&d1).abs() < 1e-6, "dF={}", v.dot(&d1));
        // The wrapper must agree.
        let e = point_curve_extrema2d(&c, &p);
        assert!((e.distance - min.distance).abs() < 1e-6, "wrapper {}", e.distance);
    }

    #[test]
    fn conic_via_dyn_curve_matches_analytic() {
        // A Geom2dEllipse driven through the dyn path: classification falls to
        // Newton (ellipse is not classified) and must still find the min.
        let e = Geom2dEllipse::from_axes(GpAx22d::standard(), 2.0, 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let em = point_curve_extrema2d(&e, &p);
        assert!((em.distance - 1.0).abs() < 1e-6, "ellipse min {}", em.distance);
        assert!((em.p2.x() - 2.0).abs() < 1e-5, "closest {em:?}");
    }

    // --- Curve-curve extrema (Phase 13 Wave 2) ---

    #[test]
    fn line_line_extrema2d_parallel_constant_distance() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let all = line_line_extrema2d(&l1, &l2);
        assert_eq!(all.len(), 1);
        assert!((all[0].distance - 3.0).abs() < 1e-9, "dist {}", all[0].distance);
    }

    #[test]
    fn line_line_extrema2d_intersecting_zero() {
        // x-axis and y-axis cross at the origin.
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(0.0, 1.0).unwrap());
        let all = line_line_extrema2d(&l1, &l2);
        assert_eq!(all.len(), 1);
        assert!(all[0].distance < 1e-9, "dist {}", all[0].distance);
        assert!((all[0].p1.x()).abs() < 1e-9 && (all[0].p1.y()).abs() < 1e-9, "{:?}", all[0]);
    }

    #[test]
    fn line_circle_extrema2d_min_max() {
        // Line y=3 and unit circle at origin: min 2, max 4.
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let all = line_circle_extrema2d(&l, &c);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
    }

    #[test]
    fn line_circle_extrema2d_intersecting_zero() {
        // Line y=0 through the unit circle: two zero-distance intersection pairs.
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let all = line_circle_extrema2d(&l, &c);
        assert!(all.iter().any(|e| e.distance < 1e-9), "no intersection {all:?}");
    }

    #[test]
    fn circle_circle_extrema2d_min_and_max() {
        // Unit circles at (0,0) and (3,0): min 1, max 5.
        let c1 = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let c2 = GpCirc2d::new(GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0));
        let all = circle_circle_extrema2d(&c1, &c2);
        assert!(all.iter().any(|e| (e.distance - 1.0).abs() < 1e-7), "min missing {all:?}");
        assert!(all.iter().any(|e| (e.distance - 5.0).abs() < 1e-7), "max missing {all:?}");
    }

    #[test]
    fn curve_curve_extrema2d_all_circles_min_max() {
        let c1 = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let c2 = Geom2dCircle::new(
            GpCirc2d::new(GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0)),
        );
        let all = curve_curve_extrema2d_all(&c1, &c2);
        assert!(!all.is_empty());
        assert!((all[0].distance - 1.0).abs() < 1e-7, "min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 5.0).abs() < 1e-7, "max {}", all[all.len() - 1].distance);
        // The samples-taking wrapper delegates to the same set.
        let es = curve_curve_extrema2d(&c1, &c2, 16);
        assert!((es[0].distance - 1.0).abs() < 1e-7, "wrapper min {}", es[0].distance);
    }

    #[test]
    fn newton_path_2d_bspline_curve_curve() {
        // Convex degree-2 arc (max y = 1 at u = 0.5) vs line y=3: min distance 2.
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let bs = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let all = curve_curve_extrema2d_all(&line, &bs);
        assert!(!all.is_empty(), "no extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-5, "min {}", min.distance);
        assert!((min.u2 - 0.5).abs() < 1e-3, "bspline param {}", min.u2);
        // Extremum condition: (C1−C2)·C′ ≈ 0 on both curves.
        let v = GpVec2d::new(min.p1.x() - min.p2.x(), min.p1.y() - min.p2.y());
        let (_, d1l) = line.d1(min.u1);
        let (_, d1b) = bs.d1(min.u2);
        assert!(v.dot(&d1l).abs() < 1e-4, "dF1 {}", v.dot(&d1l));
        assert!(v.dot(&d1b).abs() < 1e-4, "dF2 {}", v.dot(&d1b));
    }

    #[test]
    fn locate_extcc2d_refines_seed() {
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let bs = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let e = locate_extcc2d(&line, &bs, 2.0, 0.5).expect("locate failed");
        assert!((e.distance - 2.0).abs() < 1e-5, "dist {}", e.distance);
        assert!((e.u2 - 0.5).abs() < 1e-3, "param {}", e.u2);
    }
