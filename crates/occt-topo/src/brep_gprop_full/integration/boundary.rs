use super::*;

/// Build the boundary arc for an edge of `face`, honouring the edge's
/// orientation.
///
/// The face's stored **pcurve** wins (`BRepGProp_Face` trims on
/// `BRep_Tool::CurveOnSurface`): it is the only representation that keeps the
/// two sides of a seam apart (`u = 0` vs `u = 2π`) and it exists for degenerate
/// pole edges, which have no 3D curve. Without a pcurve the previous 3D-curve +
/// UV-inverse shortcut is used.
pub(in crate::brep_gprop_full) fn build_arc(e: &Edge, face: &Face, map: &UVMap) -> Option<BoundaryArc> {
    let face_key = crate::tgeometry::GeometryRegistry::shape_key(&face.0);
    let reg = crate::tgeometry::GeometryRegistry::global();
    // `BRep_Tool::CurveOnSurface` (`BRep_Tool.cxx:327-373`): on a *closed*
    // surface (`IsCurveOnClosedSurface`, i.e. a seam edge carrying two pcurves)
    // a REVERSED edge integrates `PCurve2`, any other edge `PCurve`. Picking the
    // first pcurve for both occurrences of a seam makes their boundary terms
    // cancel (which is what `BRepGProp` must not do).
    let pcurves = reg.edge_pcurves(&e.0, face_key);
    let pc = if pcurves.len() >= 2 && e.orientation().is_reversed() {
        pcurves.into_iter().nth(1)
    } else {
        pcurves.into_iter().next()
    };
    if let Some(pc) = pc {
        let (a0, b0) = reg
            .pcurve_range(&e.0, face_key)
            .unwrap_or_else(|| (pc.first_parameter(), pc.last_parameter()));
        if a0.is_finite() && b0.is_finite() && b0 > a0 {
            // `BRepGProp_Face::Load(const TopoDS_Edge&)` (`BRepGProp_Face.cxx:173-179`):
            // `C = C->Reversed(); a = C_old->ReversedParameter(b); b =
            // C_old->ReversedParameter(a);` — i.e. the arc uses the *reversed*
            // curve over the *mapped* range (per class: `-U` for a line,
            // `2*pi - U` for a circle, `first + last - U` for a BSpline).
            let (pc, a, b) = if e.orientation().is_reversed() {
                let (na, nb) = (pc.reversed_parameter(b0), pc.reversed_parameter(a0));
                (Arc::from(pc.reversed()), na, nb)
            } else {
                (pc, a0, b0)
            };
            let kind = classify_arc_kind2d(pc.as_ref(), a, b);
            return Some(BoundaryArc {
                geom: ArcGeom::Pcurve(pc),
                a,
                b,
                kind,
                // The curve itself is already reversed above.
                reversed: false,
            });
        }
    }
    let curve = BRepTool::edge_curve_world(e)?;
    let (a0, b0) = BRepTool::edge_parameters(e);
    if !(a0.is_finite() && b0.is_finite() && b0 > a0) {
        return None;
    }
    let (curve, a, b) = if e.orientation().is_reversed() {
        // `BRepGProp_Face::Load(const TopoDS_Edge&)` (`BRepGProp_Face.cxx:173-179`):
        // `C = C->Reversed(); a = C_old->ReversedParameter(b); b =
        // C_old->ReversedParameter(a);`. `ReversedParameter` is per-class
        // (`-U` for `Geom_Line`, `2*pi - U` for `Geom_Circle`, curve
        // `FirstParameter() + LastParameter() - U` for `Geom_BSplineCurve`), so
        // the mirrored range must come from the curve itself. The former
        // hardcoded `(-b0, -a0)` is only the line law: on a BSpline arc whose
        // range does not straddle 0 (T0M face 416, edge `(-0.592492, 0)`) it
        // evaluated the curve far outside the edge window and blew the face
        // area up by three orders of magnitude.
        let (na, nb) = (curve.reversed_parameter(b0), curve.reversed_parameter(a0));
        (Arc::from(curve.reversed()), na, nb)
    } else {
        (curve, a0, b0)
    };
    let kind = classify_arc_kind(curve.as_ref(), a, b);
    Some(BoundaryArc {
        geom: ArcGeom::Curve3d(curve, map.clone()),
        a,
        b,
        kind,
        // The 3D fallback already reversed the curve and its range above.
        reversed: false,
    })
}

/// [`classify_arc_kind`] for a 2D pcurve (`ArcKind::Other` is the safe answer
/// for anything that is not a straight UV line).
pub(in crate::brep_gprop_full) fn classify_arc_kind2d(c: &dyn Curve2d, a: f64, b: f64) -> ArcKind {
    let mut is_line = true;
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = c.d2(t);
        if d2.square_magnitude() > 1e-18 {
            is_line = false;
            break;
        }
    }
    if is_line {
        return ArcKind::Line;
    }
    ArcKind::Other
}

/// Classify the edge curve for integration-order purposes.
pub(in crate::brep_gprop_full) fn classify_arc_kind(c: &dyn Curve, a: f64, b: f64) -> ArcKind {
    // A line has vanishing second derivative.
    let mut is_line = true;
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = c.d2(t);
        if d2.square_magnitude() > 1e-18 {
            is_line = false;
            break;
        }
    }
    if is_line {
        return ArcKind::Line;
    }
    // A circular arc has constant tangent speed and curvature.
    let mut speeds = Vec::with_capacity(6);
    let mut accs = Vec::with_capacity(6);
    for i in 0..=5 {
        let t = a + (b - a) * i as f64 / 5.0;
        let (_, d1, d2) = c.d2(t);
        speeds.push(d1.magnitude());
        accs.push(d2.magnitude());
    }
    let s0 = speeds[0];
    let a0 = accs[0];
    if s0 > 1e-12
        && a0 > 1e-12
        && speeds.iter().all(|&s| (s - s0).abs() < 1e-6 * s0.max(1.0))
        && accs.iter().all(|&a| (a - a0).abs() < 1e-6 * a0.max(1.0))
    {
        return ArcKind::Circle;
    }
    ArcKind::Other
}

/// UV bounds of a trimmed face from its boundary pcurves (`BRepTools::UVBounds`).
pub(in crate::brep_gprop_full) fn uv_bounds(s: &dyn Surface, arcs: &[BoundaryArc]) -> (f64, f64, f64, f64) {
    let mut u1 = f64::INFINITY;
    let mut u2 = f64::NEG_INFINITY;
    let mut v1 = f64::INFINITY;
    let mut v2 = f64::NEG_INFINITY;
    for arc in arcs {
        for k in 0..=16 {
            let t = arc.a + (arc.b - arc.a) * k as f64 / 16.0;
            let p = arc.value(s, t);
            u1 = u1.min(p.x());
            u2 = u2.max(p.x());
            v1 = v1.min(p.y());
            v2 = v2.max(p.y());
        }
    }
    if u1.is_finite() && u2.is_finite() {
        (u1, u2, v1, v2)
    } else {
        let (a, b) = s.u_range();
        let (c, d) = s.v_range();
        (a, b, c, d)
    }
}
