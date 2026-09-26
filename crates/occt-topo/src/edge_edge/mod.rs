//! Edge–edge intersection — port of `IntTools_EdgeEdge`.
//!
//! Orchestrates the intersection of two edges into two kinds of results:
//!
//! - **discrete vertex hits** — [`PntOn2Faces`] entries carrying the parameter
//!   on each edge and the 3D point (`points()`);
//! - **coincident common parts** — [`CommonPrt`] entries with [`CommonPartType::Edge`]
//!   describing an overlapping sub-range of two (nearly) coincident edges
//!   (`common_parts()`).
//!
//! Line–line is dispatched to `compute_line_line` (port of
//! `IntTools_EdgeEdge::ComputeLineLine`); **every other combination** —
//! circles, BSpline/Bezier and general curves alike — runs the faithful
//! `IntTools_EdgeEdge::FindSolutions` parameter-box recursion in
//! [`find_solutions`](EdgeEdge::find_solutions) (module `find_solutions`), whose
//! common parts and vertex parameters come from `MergeSolutions` /
//! `AddSolution` / `FindBestSolution`. OCCT has no circle/circle branch, so the
//! port-local fast path was removed (board task R2-17). Coincidence is detected
//! by sampling one curve and projecting the samples onto the other (port of
//! `IntTools_EdgeEdge::IsCoincident`).
//!
//! The task spec suggested `geom2d_api::project_point_on_curve` for parameter
//! refinement; that is a 2-D API, so the 3-D analog `geom_api::project_point_on_curve`
//! is used instead.
mod prelude {

pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpLin, GpPnt, GpVec};
pub(crate) use occt_core::precision::ANGULAR;
pub(crate) use occt_geom::geom_api;
pub(crate) use occt_geom::Curve;

pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::inttools::{edge_edge_intersections, EdgeEdgeHit};
pub(crate) use crate::inttools_data::{CommonPartType, CommonPrt, IntRange, PntOn2Faces};
pub(crate) use crate::shape::Edge;

}

use prelude::*;


// ---------------------------------------------------------------------------
// Phase 18c tests — the appended full solvers and best-solution selection.
// (Deliberately a separate test module so the original tests stay untouched.)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_full {
    use super::*;

    use crate::builder::TopoBuilder;
    use crate::shape::{Edge, TopoShape};
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    /// Set two edges and run `prepare`, so the full solvers can run directly.
    fn prepared_ee(e1: &Edge, e2: &Edge, fuzzy: f64) -> EdgeEdge {
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_fuzzy_value(fuzzy);
        ee.prepare().unwrap();
        ee
    }

    #[test]
    fn line_line_full_coincident_segment() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(4.0, 0.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 1e-7);
        let found = ee.compute_line_line_full().unwrap();
        assert!(found, "coincident lines intersect");
        assert_eq!(ee.common_parts().len(), 1, "common: {:?}", ee.common_parts());
        // coincident_range covers the overlapping domain [1, 3].
        let cr = ee.coincident_range().unwrap();
        assert!((cr.first - 1.0).abs() < 1e-6, "first={}", cr.first);
        assert!((cr.last - 3.0).abs() < 1e-6, "last={}", cr.last);
        // intersection_points exposes the segment endpoints.
        let pts = ee.intersection_points();
        assert_eq!(pts.len(), 2, "span endpoints: {pts:?}");
        assert!((pts[0].0 - 1.0).abs() < 1e-6, "t1={}", pts[0].0);
        assert!((pts[1].0 - 3.0).abs() < 1e-6, "t1={}", pts[1].0);
        assert!(pts[0].2.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-6);
        assert!(pts[1].2.distance(&GpPnt::new(3.0, 0.0, 0.0)) < 1e-6);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_full_crossing_single_point() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 0.0);
        let found = ee.compute_line_line_full().unwrap();
        assert!(found);
        assert_eq!(ee.points().len(), 1, "points: {:?}", ee.points());
        let p = &ee.points()[0];
        assert!(p.pnt1.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6, "pnt={:?}", p.pnt1);
        let d = 2.0f64.sqrt();
        assert!((p.uv1.0 - d).abs() < 1e-6, "u1={}", p.uv1.0);
        assert!((p.uv2.0 - d).abs() < 1e-6, "u2={}", p.uv2.0);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_full_parallel_distinct_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(2.0, 1.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 0.0);
        let found = ee.compute_line_line_full().unwrap();
        assert!(!found);
        assert!(ee.points().is_empty());
        assert!(ee.common_parts().is_empty());
        assert!(ee.coincident_range().is_none());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn find_best_solution_selects_minimum_residual() {
        let b = TopoBuilder::new();
        // X-line crossing a Y-line at (1, 0, 0): c2 runs through (1,0,0) at
        // its parameter 1 (segment (1,-1,0)→(1,1,0)).
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let ee = prepared_ee(&e1, &e2, 0.0);

        let exact = (1.0, 1.0, GpPnt::new(1.0, 0.0, 0.0));
        let far = (0.0, 0.0, GpPnt::new(0.0, 0.0, 0.0)); // residual ≈ √2 → dropped
        let best = ee.find_best_solution(&[far, exact], 1e-6).unwrap();
        assert_eq!(best.0, 1.0);
        assert_eq!(best.1, 1.0);
        assert!(best.2.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-12);

        // All-degenerate (or empty) candidate lists yield None.
        assert!(ee.find_best_solution(&[far, far], 1e-6).is_none());
        assert!(ee.find_best_solution(&[], 1e-6).is_none());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }
}

mod edge_edge;
mod find_solutions;
mod solvers;
pub use edge_edge::*;


#[cfg(test)]
#[path = "tests.rs"]
mod tests;
