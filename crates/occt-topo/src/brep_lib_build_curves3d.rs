//! `BRepLib::BuildCurve3d` / `BuildCurves3d` (`BRepLib.cxx:301-489`).
//!
//! `TShape::Locked` (`BRep_Builder.cxx:641-643`) is not stored, so `UpdateEdge`
//! always writes. `UpdateEdge` stores `L.Predivided(E.Location())`
//! (`cxx:645`); the port keeps the curve in the edge frame and applies the
//! edge location at evaluation.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom::geom_lib;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_lib_same_range::{check_same_range, same_range, set_rep_ranges};
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::topo_tools_full::edges_of;

/// `BRepLib::BuildCurves3d(S)` (`BRepLib.cxx:460-463`): tolerance `1.0e-5`,
/// `GeomAbs_C1`, degree 14, `MaxSegment` 0 (`BRepLib.hxx:99-103`).
pub fn build_curves3d(shape: &TopoShape) -> bool {
    build_curves3d_tol(shape, 1.0e-5, Shape::C1, 14, 0)
}

/// `BRepLib::BuildCurves3d(S, Tolerance, Continuity, MaxDegree, MaxSegment)`
/// (`BRepLib.cxx:467-489`).
pub fn build_curves3d_tol(
    shape: &TopoShape,
    tolerance: f64,
    continuity: Shape,
    max_degree: i32,
    max_segment: i32,
) -> bool {
    let mut ok = true;
    let mut seen: Vec<TopoShape> = Vec::new();
    for edge in edges_of(shape) {
        if seen.iter().any(|k| same_located(k, &edge.0)) {
            continue;
        }
        seen.push(edge.0.clone());
        let built = build_curve3d(&edge, tolerance, continuity, max_degree, max_segment);
        ok = ok && built;
    }
    ok
}

/// `BRepLib::BuildCurve3d` (`BRepLib.cxx:301-455`).
pub fn build_curve3d(
    edge: &Edge,
    tolerance: f64,
    continuity: Shape,
    max_degree: i32,
    max_segment: i32,
) -> bool {
    let reg = GeometryRegistry::global();
    // `cxx:321-325`: an edge that already has a 3D curve returns true.
    if reg.edge_curve(&edge.0).is_some() {
        return true;
    }
    // `cxx:330-333`.
    if !check_same_range(&edge.0, CONFUSION) {
        same_range(&edge.0, tolerance);
    }

    // `cxx:335-357`: search a curve on a plane. The loop stops on the first
    // plane, including a trimmed plane's basis.
    let reps = reg.edge_pcurve_reps(&edge.0);
    let mut plane = None;
    for (surf, pc, f, l) in &reps {
        let basis = surf.rectangular_trimmed_basis().unwrap_or_else(|| surf.clone());
        if let Some(pln) = basis.gp_pln() {
            plane = Some((pc.clone(), *f, *l, pln.position().ax2()));
            break;
        }
    }

    if let Some((pc, f, l, axes)) = plane {
        // `cxx:361-365`: `GeomLib::To3d` null handle returns false.
        let Some(c3d) = geom_lib::to_3d(&axes, pc.as_ref()) else {
            return false;
        };
        // `cxx:371`: `UpdateEdge(..., 0)` raises tolerance only, so the stored
        // tolerance stays. `cxx:372-373` copies the plane pcurve range onto the
        // 3D curve.
        let mut g = EdgeGeom::new(c3d, f, l);
        g.tolerance = reg.edge_tolerance(&edge.0);
        g.same_parameter = reg.same_parameter(&edge.0);
        g.same_range = reg.same_range(&edge.0);
        g.degenerated = reg.is_degenerated_edge(&edge.0);
        reg.set_edge(&edge.0, g);
        set_rep_ranges(&edge.0, f, l);
        return true;
    }

    // `cxx:382-452`.
    if reg.is_degenerated_edge(&edge.0) {
        return false;
    }
    let jj = reps.len().min(2);
    let Some((surf, pc, f, l)) = reps.first().cloned() else {
        return false;
    };
    let cos: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(pc.clone(), surf.clone(), f, l));
    let max_segment = evaluate_max_segment(max_segment, &surf, &pc);
    let result = geom_lib::geom_lib_build_curve3d(
        tolerance,
        cos.as_ref(),
        &surf,
        &pc,
        f,
        l,
        continuity,
        max_degree,
        max_segment,
    );
    // `cxx:435`: `max(BRep_Tool::Tolerance(E), Tolerance)`.
    let max_deviation = reg.edge_tolerance(&edge.0).max(tolerance);
    let Some(c3d) = result.curve else {
        return false;
    };
    let mut g = EdgeGeom::new(c3d, f, l);
    g.tolerance = max_deviation;
    g.same_parameter = reg.same_parameter(&edge.0);
    g.same_range = reg.same_range(&edge.0);
    g.degenerated = reg.is_degenerated_edge(&edge.0);
    reg.set_edge(&edge.0, g);
    // `cxx:441-447`.
    if jj == 1 {
        reg.set_same_parameter(&edge.0, true);
    }
    true
}

/// `evaluateMaxSegment` (`BRepLib.cxx:273-296`).
fn evaluate_max_segment(
    max_segment: i32,
    surface: &Arc<dyn Surface>,
    pcurve: &Arc<dyn Curve2d>,
) -> i32 {
    if max_segment != 0 {
        return max_segment;
    }
    let mut nb_s_knots = 0usize;
    if surface.is_bspline_surface() {
        let nu = surface.bspline_surface_uknots().map_or(0, |k| k.len());
        let nv = surface.bspline_surface_vknots().map_or(0, |k| k.len());
        nb_s_knots = nu.max(nv);
    }
    let nb_c2d_knots = pcurve.bspline_nb_knots().unwrap_or(0);
    (30 + nb_s_knots.max(nb_c2d_knots)) as i32
}

/// `TopTools_ShapeMapHasher` / `IsSame`: same `TShape` and the same location.
fn same_located(a: &TopoShape, b: &TopoShape) -> bool {
    crate::topo_tools_full::is_same(a, b)
}
