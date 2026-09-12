use super::prelude::*;
use super::*;
    use crate::curve::Curve;
    use crate::{GeomBSplineCurve, GeomCircle, GeomLine};
    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPnt, GpVec};

    fn line_on_x_axis() -> GpLin {
        GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap())
    }

    fn line_along_y_at_z3() -> GpLin {
        GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(0.0, 1.0, 0.0).unwrap())
    }

    #[test]
    fn line_line_skew_min_distance() {
        let l1 = line_on_x_axis();
        let l2 = line_along_y_at_z3();
        let all = line_line_extrema(&l1, &l2);
        assert_eq!(all.len(), 1);
        let e = all[0];
        assert!((e.distance - 3.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1).abs() < 1e-9 && (e.u2).abs() < 1e-9, "params {} {}", e.u1, e.u2);
        assert!((e.p1.x()).abs() < 1e-9 && (e.p2.z() - 3.0).abs() < 1e-9, "{e:?}");
    }

    #[test]
    fn line_line_parallel_constant_distance() {
        let l1 = line_on_x_axis();
        let l2 = GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let all = line_line_extrema(&l1, &l2);
        assert_eq!(all.len(), 1);
        assert!((all[0].distance - 3.0).abs() < 1e-9, "dist {}", all[0].distance);
    }

    #[test]
    fn line_circle_planar_min_and_max() {
        // Line along x through (0,0,3), unit circle at origin in the xy-plane:
        // min distance 3 (line point above the nearest circle point), max √10.
        let l = GpLin::from_pnt_dir(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let c = GpCirc::new(GpAx2::standard(), 1.0);
        let all = line_circle_extrema(&l, &c);
        assert!(!all.is_empty(), "no solutions");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 3.0).abs() < 1e-7, "min {}", min.distance);
        assert!((max.distance - (10.0f64).sqrt()).abs() < 1e-7, "max {}", max.distance);
    }

    #[test]
    fn line_circle_perpendicular_min_and_max() {
        // Line through (2,0,0) along z, unit circle at origin: min 1, max 3.
        let l = GpLin::from_pnt_dir(GpPnt::new(2.0, 0.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap());
        let c = GpCirc::new(GpAx2::standard(), 1.0);
        let all = line_circle_extrema(&l, &c);
        assert!(!all.is_empty(), "no solutions");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 1.0).abs() < 1e-7, "min {}", min.distance);
        assert!((max.distance - 3.0).abs() < 1e-7, "max {}", max.distance);
    }

    #[test]
    fn circle_circle_min_and_max_both_returned() {
        // Unit circles at origin and (3,0,0): min 1, max 5.
        let c1 = GpCirc::new(GpAx2::standard(), 1.0);
        let c2 = GpCirc::new(GpAx2::standard(), 1.0).translated_vec(&GpVec::new(3.0, 0.0, 0.0));
        let all = circle_circle_extrema(&c1, &c2);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 1.0).abs() < 1e-7, "min {}", min.distance);
        assert!((max.distance - 5.0).abs() < 1e-7, "max {}", max.distance);
        // Both distances must be present in the returned set.
        assert!(all.iter().any(|e| (e.distance - 1.0).abs() < 1e-7), "min missing {all:?}");
        assert!(all.iter().any(|e| (e.distance - 5.0).abs() < 1e-7), "max missing {all:?}");
    }

    #[test]
    fn circle_circle_intersecting_zero_distance() {
        // Unit circles at origin and (1,0,0) cross; min distance 0, max 3.
        let c1 = GpCirc::new(GpAx2::standard(), 1.0);
        let c2 = GpCirc::new(GpAx2::standard(), 1.0).translated_vec(&GpVec::new(1.0, 0.0, 0.0));
        let all = circle_circle_extrema(&c1, &c2);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!(min.distance < 1e-6, "min {}", min.distance);
        assert!((max.distance - 3.0).abs() < 1e-6, "max {}", max.distance);
    }

    #[test]
    fn curve_curve_all_skew_lines() {
        let c1 = GeomLine::new(line_on_x_axis());
        let c2 = GeomLine::new(line_along_y_at_z3());
        let all = curve_curve_extrema_all(&c1, &c2);
        assert!(!all.is_empty());
        assert!((all[0].distance - 3.0).abs() < 1e-9, "min {}", all[0].distance);
    }

    #[test]
    fn curve_curve_all_circles_min_and_max() {
        let c1 = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let c2 = GeomCircle::new(
            GpCirc::new(GpAx2::standard(), 1.0).translated_vec(&GpVec::new(3.0, 0.0, 0.0)),
        );
        let all = curve_curve_extrema_all(&c1, &c2);
        assert!(!all.is_empty());
        assert!((all[0].distance - 1.0).abs() < 1e-7, "min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 5.0).abs() < 1e-7, "max {}", all[all.len() - 1].distance);
        // The convenience minimum must agree.
        let e = curve_curve_extrema(&c1, &c2);
        assert!((e.distance - 1.0).abs() < 1e-7, "min wrapper {}", e.distance);
    }

    #[test]
    fn newton_path_line_vs_bspline() {
        // Convex degree-2 arc with max y = 1 at u = 0.5 (point (2,1,0)); a line
        // along x through (0,3,0) gives min distance 2 at u = 0.5.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let bs = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(2.0, 2.0, 0.0),
                GpPnt::new(4.0, 0.0, 0.0),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let all = curve_curve_extrema_all(&line, &bs);
        assert!(!all.is_empty(), "no extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-5, "min {}", min.distance);
        assert!((min.u2 - 0.5).abs() < 1e-3, "bspline param {}", min.u2);
        // Extremum condition: (C1−C2)·C′ ≈ 0 on both curves.
        let v = GpVec::from_pnts(&min.p2, &min.p1);
        let (_, d1l) = line.d1(min.u1);
        let (_, d1b) = bs.d1(min.u2);
        assert!(v.dot(&d1l).abs() < 1e-4, "dF1 {}", v.dot(&d1l));
        assert!(v.dot(&d1b).abs() < 1e-4, "dF2 {}", v.dot(&d1b));
    }

    #[test]
    fn locate_extcc_refines_seed() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let bs = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(2.0, 2.0, 0.0),
                GpPnt::new(4.0, 0.0, 0.0),
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let e = locate_extcc(&line, &bs, 2.0, 0.5).expect("locate failed");
        assert!((e.distance - 2.0).abs() < 1e-5, "dist {}", e.distance);
        assert!((e.u2 - 0.5).abs() < 1e-3, "param {}", e.u2);
    }

    #[test]
    fn curve_curve_all_bspline_min_matches_analytic() {
        // A straight degree-1 B-spline classified as a line vs a circle: the
        // analytic dispatch must find the same min as the geometric line-circle.
        let line_bs = GeomBSplineCurve::new(
            vec![GpPnt::new(-2.0, 0.0, 3.0), GpPnt::new(4.0, 0.0, 3.0)],
            vec![0.0, 0.0, 1.0, 1.0],
            1,
        )
        .unwrap();
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let all = curve_curve_extrema_all(&line_bs, &circ);
        assert!(!all.is_empty());
        // Line at height 3 above the circle plane → min distance 3.
        assert!((all[0].distance - 3.0).abs() < 1e-5, "min {}", all[0].distance);
    }
