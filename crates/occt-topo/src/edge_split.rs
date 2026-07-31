//! Edge splitting — subdivide an edge's parameter range into several edges.
//!
//! Source: `GeomLib::SplitCurve` / `BRepAlgoAPI_Splitter` (TKTopAlgo). Each
//! sub-edge wraps the original curve in a `GeomTrimmedCurve` over the sub-range
//! and registers an `EdgeGeom` whose `first`/`last` are the sub-range expressed
//! in the *original* curve's parameter space (matching `BRep_TEdge` semantics).
//!
//! Note: `GeomTrimmedCurve` reparametrizes to `[0, 1]`, so to evaluate the
//! sub-edge at its stored parameter `t ∈ [first, last]` use
//! `curve.d0((t - first) / (last - first))`.

use std::sync::Arc;

use occt_geom::trimmed::GeomTrimmedCurve;

use crate::builder::TopoBuilder;
use crate::shape::Edge;
use crate::tgeometry::GeometryRegistry;

/// Split `e` at the given interior parameters. `params` are sorted, deduplicated
/// and clamped to the edge's `[first, last]` range; the result has
/// `params.len() + 1` edges covering `[first, last]` in order. Edges without a
/// registered curve (or with an unbounded range) are returned unchanged.
pub fn split_edge(e: &Edge, params: &[f64]) -> Vec<Edge> {
    let reg = GeometryRegistry::global();
    let Some(curve) = reg.edge_curve(&e.0) else {
        return vec![e.clone()];
    };
    let (first, last) = reg.edge_parameters(&e.0);
    if params.is_empty() || !first.is_finite() || !last.is_finite() || (first - last).abs() <= 1e-15 {
        return vec![e.clone()];
    }
    let (a0, a1) = (first.min(last), first.max(last));

    let mut ps: Vec<f64> = params
        .iter()
        .copied()
        .filter(|p| p.is_finite())
        .map(|p| p.clamp(a0, a1))
        .collect();
    ps.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    ps.dedup_by(|x, y| (*x - *y).abs() <= 1e-12);

    let mut bounds = vec![a0];
    bounds.extend(ps.into_iter().filter(|p| *p > a0 && *p < a1));
    bounds.push(a1);
    bounds.dedup_by(|x, y| (*x - *y).abs() <= 1e-12);

    let b = TopoBuilder::new();
    let mut out = Vec::with_capacity(bounds.len() - 1);
    for w in bounds.windows(2) {
        let (u0, u1) = (w[0], w[1]);
        let trimmed = GeomTrimmedCurve::new(curve.clone(), u0, u1);
        out.push(b.make_edge(Arc::new(trimmed), u0, u1));
    }
    out
}

/// Alias for [`split_edge`].
pub fn split_edge_at_curve_params(e: &Edge, u: &[f64]) -> Vec<Edge> {
    split_edge(e, u)
}

/// Batch split: apply [`split_edge`] to every edge with its own parameter list.
/// Edges with no entry in `params_by_edge` are returned unsplit.
pub fn split_edges_at(edges: &[Edge], params_by_edge: &[Vec<f64>]) -> Vec<Edge> {
    let mut out = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        let params = params_by_edge.get(i).map(Vec::as_slice).unwrap_or(&[]);
        out.extend(split_edge(e, params));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_geom::Curve;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;
    use crate::brep_tool::BRepTool;
    use occt_core::gp::GpPnt;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    #[test]
    fn split_line_edge_into_three() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        let parts = split_edge(&e, &[3.0, 7.0]);
        assert_eq!(parts.len(), 3);

        let ranges: Vec<(f64, f64)> = parts
            .iter()
            .map(|p| BRepTool::edge_parameters(p))
            .collect();
        assert_eq!(ranges, vec![(0.0, 3.0), (3.0, 7.0), (7.0, 10.0)]);

        // Geometrical continuity at the splices: the trimmed curve endpoint of
        // part i equals the trimmed curve start of part i+1.
        let c0 = BRepTool::edge_curve(&parts[0]).expect("curve");
        let c1 = BRepTool::edge_curve(&parts[1]).expect("curve");
        let c2 = BRepTool::edge_curve(&parts[2]).expect("curve");
        let p0 = c0.d0(0.0); // start of part 0 → curve at 0
        let p1 = c0.d0(1.0); // end of part 0 → curve at 3
        let p2 = c1.d0(0.0); // start of part 1 → curve at 3
        let p3 = c1.d0(1.0); // end of part 1 → curve at 7
        let p4 = c2.d0(0.0); // start of part 2 → curve at 7
        let p5 = c2.d0(1.0); // end of part 2 → curve at 10
        assert!(p0.is_equal(&GpPnt::new(0.0, 0.0, 0.0)));
        assert!(p1.is_equal(&GpPnt::new(3.0, 0.0, 0.0)));
        assert!(p2.is_equal(&GpPnt::new(3.0, 0.0, 0.0)));
        assert!(p3.is_equal(&GpPnt::new(7.0, 0.0, 0.0)));
        assert!(p4.is_equal(&GpPnt::new(7.0, 0.0, 0.0)));
        assert!(p5.is_equal(&GpPnt::new(10.0, 0.0, 0.0)));

        for p in &parts {
            clear_tree(&p.0);
        }
        clear_tree(&e.0);
    }

    #[test]
    fn split_handles_dedup_clamp_and_out_of_order() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        // Out of order + duplicates + outside-range + interior.
        let parts = split_edge(&e, &[9.0, 2.0, 2.0, -5.0, 5.0, 100.0]);
        let ranges: Vec<(f64, f64)> = parts
            .iter()
            .map(|p| BRepTool::edge_parameters(p))
            .collect();
        // -5 clamps to 0, 100 clamps to 10; 2 dup deduped; interior 2, 5, 9.
        assert_eq!(ranges, vec![(0.0, 2.0), (2.0, 5.0), (5.0, 9.0), (9.0, 10.0)]);
        for p in &parts {
            clear_tree(&p.0);
        }
        clear_tree(&e.0);
    }

    #[test]
    fn empty_params_returns_single_edge() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(5.0, 0.0, 0.0));
        let parts = split_edge(&e, &[]);
        assert_eq!(parts.len(), 1);
        assert!(parts[0].same_tshape(&e.0));
        clear_tree(&e.0);
    }

    #[test]
    fn split_edges_at_batches() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(4.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(4.0, 1.0, 0.0));
        let out = split_edges_at(&[e1.clone(), e2.clone()], &[vec![2.0], vec![]]);
        assert_eq!(out.len(), 3); // e1 → 2 parts, e2 unsplit → 1
        assert_eq!(BRepTool::edge_parameters(&out[0]), (0.0, 2.0));
        assert_eq!(BRepTool::edge_parameters(&out[1]), (2.0, 4.0));
        assert_eq!(BRepTool::edge_parameters(&out[2]), (0.0, 4.0));
        for p in &out {
            clear_tree(&p.0);
        }
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }
}
