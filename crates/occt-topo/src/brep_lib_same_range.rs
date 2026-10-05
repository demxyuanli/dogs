//! `BRepLib::CheckSameRange` and `BRepLib::SameRange`.
//! Source: `BRepLib.cxx:149-183` and `BRepLib.cxx:187-263`.
//!
//! Both walk `BRep_TEdge::ChangeCurves()` and compare the `[First, Last]` of
//! every `BRep_GCurve` representation (the 3D curve representation is one of
//! them). `SameRange` reparameterises the pcurves that disagree with the
//! reference range through `GeomLib::SameRange` and then writes the reference
//! range to every representation.
//!
//! The port stores the 3D curve range in the edge core (`EdgeGeom::first` /
//! `last`) and the CurveOnSurface ranges in `EdgePcurves::ranges`, keyed by
//! surface pointer (`GeometryRegistry::repr_key`).

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::shape::TopoShape;
use crate::shhealing::geom_lib_same_range;
use crate::tgeometry::GeometryRegistry;

/// One `BRep_GCurve` representation of an edge with its `[First, Last]`.
pub struct RepRange {
    /// `None` for the 3D curve representation (whose range is the edge core
    /// range), `Some(repr_key)` for a CurveOnSurface representation.
    pub key: Option<usize>,
    pub first: f64,
    pub last: f64,
}

/// Every `BRep_GCurve` representation of `s` with its `[First, Last]`, the 3D
/// curve representation first when the edge carries one.
///
/// UNPORTED (ordering): OCCT walks an insertion-ordered list, where the 3D
/// curve representation is appended by `UpdateCurves`
/// (`BRep_Builder.cxx:58-95`), i.e. never the head once the edge carries
/// pcurves. The port's CurveOnSurface slots live in a `HashMap`, so they are
/// emitted in key order and the 3D range is placed first. Both orders agree
/// whenever the ranges are already consistent, which is the state
/// `BuildCurve3d` requires.
pub fn rep_ranges(s: &TopoShape) -> Vec<RepRange> {
    let ts = s.tshape.read().expect("poisoned TShape lock");
    let Some(core) = ts.edge_core() else {
        return Vec::new();
    };
    let mut out = Vec::new();
    if core.curve.is_some() {
        out.push(RepRange {
            key: None,
            first: core.first,
            last: core.last,
        });
    }
    if let Some(g) = ts.edge_pcurves() {
        let mut keys: Vec<usize> = g.ranges.keys().copied().collect();
        keys.sort_unstable();
        for k in keys {
            let (first, last) = g.ranges[&k];
            out.push(RepRange {
                key: Some(k),
                first,
                last,
            });
        }
    }
    out
}

/// `BRep_Builder::Range(E, First, Last)` (`BRep_Builder.cxx:1091-1118`) with
/// `Only3d = false`: every representation gets the range.
pub(crate) fn set_rep_ranges(s: &TopoShape, first: f64, last: f64) {
    let mut ts = s.tshape.write().expect("poisoned TShape lock");
    let c = ts.edge_core_mut();
    c.first = first;
    c.last = last;
    if ts.edge_pcurves().is_some() {
        for r in ts.edge_pcurves_mut().ranges.values_mut() {
            *r = (first, last);
        }
    }
}

/// `BRepLib::CheckSameRange(E, Tolerance)` (`BRepLib.cxx:149-183`): are the
/// `[First, Last]` of every representation within `tolerance` of each other?
pub fn check_same_range(s: &TopoShape, tolerance: f64) -> bool {
    let mut current: Option<(f64, f64)> = None;
    for r in rep_ranges(s) {
        match current {
            None => current = Some((r.first, r.last)),
            Some((cf, cl)) => {
                if (cf - r.first).abs() > tolerance || (cl - r.last).abs() > tolerance {
                    return false;
                }
            }
        }
    }
    true
}

/// `BRepLib::SameRange(E, Tolerance)` (`BRepLib.cxx:187-263`).
///
/// UNPORTED: `Tolerance` is the first argument of `GeomLib::SameRange`
/// (`cxx:237-254`); the port's `shhealing::geom_lib_same_range` does not carry
/// it (pre-existing simplification of that function).
pub fn same_range(s: &TopoShape, _tolerance: f64) {
    let reg = GeometryRegistry::global();
    let reps = rep_ranges(s);
    let Some(head) = reps.first() else {
        return;
    };
    // `first_time_in` (`cxx:224-229`).
    let (cf, cl) = (head.first, head.last);
    for r in &reps[1..] {
        // `cxx:231-234`.
        if (r.first - cf).abs() <= CONFUSION && (r.last - cl).abs() <= CONFUSION {
            continue;
        }
        let Some(key) = r.key else {
            continue;
        };
        // `GeomLib::SameRange(Tolerance, Curve2dPtr, GC->First(), GC->Last(),
        // current_first, current_last, NewCurve2dPtr)` (`cxx:236-254`) for the
        // PCurve and, on a closed-surface representation, the PCurve2.
        let pcs = reg.edge_pcurves(s, key);
        if pcs.is_empty() {
            continue;
        }
        let new_pcs: Vec<Arc<dyn Curve2d>> = pcs
            .into_iter()
            .map(|pc| geom_lib_same_range(pc, r.first, r.last, cf, cl))
            .collect();
        if new_pcs.len() == 1 {
            reg.set_edge_pcurve(s, key, new_pcs.into_iter().next().unwrap());
        } else {
            reg.set_edge_pcurves(s, key, new_pcs);
        }
    }
    // `B.Range(TopoDS::Edge(AnEdge), current_first, current_last)` and
    // `B.SameRange(AnEdge, true)` (`cxx:260-263`).
    set_rep_ranges(s, cf, cl);
    reg.set_same_range(s, true);
}
