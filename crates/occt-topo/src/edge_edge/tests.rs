use super::prelude::*;
use super::*;
    use std::cmp::Ordering;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpCirc, GpDir};
    use occt_geom::{GeomBSplineCurve, GeomCircle};

    use crate::builder::TopoBuilder;
    use crate::inttools::edge_edge_intersections;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn run(e1: &Edge, e2: &Edge, fuzzy: f64) -> EdgeEdge {
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_fuzzy_value(fuzzy);
        ee.perform().unwrap();
        ee
    }

    fn sorted_x(ee: &EdgeEdge) -> Vec<f64> {
        let mut xs: Vec<f64> = ee.points().iter().map(|p| p.pnt1.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        xs
    }

    #[test]
    fn line_line_crossing_single_hit() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.is_done());
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        let pts = ee.points();
        assert_eq!(pts.len(), 1, "points: {pts:?}");
        assert!(pts[0].pnt1.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        let d = 2.0f64.sqrt();
        assert!((pts[0].uv1.0 - d).abs() < 1e-6, "u1={}", pts[0].uv1.0);
        assert!((pts[0].uv2.0 - d).abs() < 1e-6, "u2={}", pts[0].uv2.0);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_coincident_overlap_common_part() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(4.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 1e-7);
        assert_eq!(ee.common_parts().len(), 1, "common: {:?}", ee.common_parts());
        let cp = &ee.common_parts()[0];
        assert_eq!(cp.part_type, CommonPartType::Edge);
        assert!((cp.range.first - 1.0).abs() < 1e-6, "first={}", cp.range.first);
        assert!((cp.range.last - 3.0).abs() < 1e-6, "last={}", cp.range.last);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_collinear_non_overlap_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 1e-7);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_parallel_distinct_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(2.0, 1.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.points().is_empty());
        assert!(ee.common_parts().is_empty());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn circle_circle_two_hits() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ee = run(&c1, &c2, 1e-7);
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        let pts = ee.points();
        assert_eq!(pts.len(), 2, "points: {pts:?}");
        for p in pts {
            assert!((p.pnt1.x() - 0.5).abs() < 1e-6, "x={}", p.pnt1.x());
            assert!((p.pnt1.y().abs() - 0.75f64.sqrt()).abs() < 1e-4, "y={}", p.pnt1.y());
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn separated_circles_empty() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(5.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ee = run(&c1, &c2, 1e-7);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn line_bspline_crossing_two_hits() {
        let b = TopoBuilder::new();
        let bs = Arc::new(
            GeomBSplineCurve::new(
                vec![
                    GpPnt::new(0.0, 0.0, 0.0),
                    GpPnt::new(2.0, 2.0, 0.0),
                    GpPnt::new(4.0, 0.0, 0.0),
                ],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                2,
            )
            .unwrap(),
        );
        let bs_edge = b.make_edge(bs, 0.0, 1.0);
        let line = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(4.0, 0.5, 0.0));
        let ee = run(&line, &bs_edge, 1e-7);
        let pts = ee.points();
        assert_eq!(pts.len(), 2, "points: {pts:?}");
        let xs = sorted_x(&ee);
        // Quadratic Bezier (0,0)-(2,2)-(4,0): x=4t, y=4t(1-t); the line y=0.5
        // meets it at t=(1±1/√2)/2 → x=4t = 2∓√2.
        let xa = 2.0f64 - std::f64::consts::SQRT_2;
        let xb = 2.0f64 + std::f64::consts::SQRT_2;
        assert!((xs[0] - xa).abs() < 1e-3, "x0={} expected {xa}", xs[0]);
        assert!((xs[1] - xb).abs() < 1e-3, "x1={} expected {xb}", xs[1]);
        // The line parameter equals the x coordinate; the BSpline parameter
        // equals x/4.
        let mut u1s: Vec<f64> = pts.iter().map(|p| p.uv1.0).collect();
        u1s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        assert!((u1s[0] - xa).abs() < 1e-3, "u1={}", u1s[0]);
        let mut u2s: Vec<f64> = pts.iter().map(|p| p.uv2.0).collect();
        u2s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        assert!((u2s[0] - xa / 4.0).abs() < 1e-3, "u2={}", u2s[0]);
        assert!((u2s[1] - xb / 4.0).abs() < 1e-3, "u2={}", u2s[1]);
        clear_tree(&line.0);
        clear_tree(&bs_edge.0);
    }

    #[test]
    fn separated_edges_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(1.0, 2.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        assert!(ee.common_parts().is_empty());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn matches_inttools_line_line() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let ref_hits = edge_edge_intersections(&e1, &e2, 1e-9);
        let ee = run(&e1, &e2, 0.0);
        assert_eq!(ref_hits.len(), ee.points().len(), "ref={ref_hits:?} ee={:?}", ee.points());
        for (h, p) in ref_hits.iter().zip(ee.points()) {
            assert!((h.u1 - p.uv1.0).abs() < 1e-6, "u1 {} vs {}", h.u1, p.uv1.0);
            assert!((h.u2 - p.uv2.0).abs() < 1e-6, "u2 {} vs {}", h.u2, p.uv2.0);
            assert!(h.point.distance(&p.pnt1) < 1e-6);
        }
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn matches_inttools_circle_circle() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ref_hits = edge_edge_intersections(&c1, &c2, 1e-6);
        let ee = run(&c1, &c2, 1e-7);
        assert_eq!(ref_hits.len(), ee.points().len(), "ref={ref_hits:?} ee={:?}", ee.points());
        let mut ref_xs: Vec<f64> = ref_hits.iter().map(|h| h.point.x()).collect();
        let mut ee_xs: Vec<f64> = ee.points().iter().map(|p| p.pnt1.x()).collect();
        ref_xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        ee_xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        for (r, e) in ref_xs.iter().zip(&ee_xs) {
            assert!((r - e).abs() < 1e-6, "x {r} vs {e}");
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn range_restriction_filters_hits() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        // Restrict edge 1 to [1.9, 2.0]-ish parameters: the crossing at u=√2
        // is excluded.
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_range1(IntRange::new(2.0, 2.8).unwrap());
        ee.perform().unwrap();
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }
