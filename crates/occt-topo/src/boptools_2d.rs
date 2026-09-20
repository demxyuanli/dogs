//! Edge→face 2D (p-curve) utilities — a complete port of
//! `BOPTools_AlgoTools2D` (TKBO).
//!
//! The class groups the static helpers the boolean pipeline needs to build and
//! interrogate the 2D (UV-domain) representation of an edge on a face:
//!
//! * **Construction** — [`make_2d`] builds the edge→face p-curve (delegating to
//!   [`crate::pcurve_full::make_pcurve_full`]); [`build_pcurve_for_edge_on_face`]
//!   builds it *and* attaches it to the edge; [`attach_existing_pcurve`] registers
//!   an already-built curve.
//! * **Storage access** — [`curve_on_surface`] reads the stored p-curve of an
//!   edge on a face from the [`crate::tgeometry::GeometryRegistry`] side-table;
//!   [`has_curve_on_surface`] reports whether one exists.
//! * **Interrogation** — [`is_edge_isoline`] decides whether the p-curve is a
//!   `u`- or `v`-isoline of the face surface; [`edge_tangent`] gives the tangent
//!   of the edge's 3D curve at a parameter; [`point_on_surface`] evaluates the
//!   face surface at `(u, v)`.
//! * **Sampling** — [`intermediate_point`] returns the OCCT ~43.2% division
//!   parameter between two values.
//!
//! The `Arc<dyn Curve2d>` p-curves are stored per-face in `EdgeGeom::pcurves`,
//! keyed by face pointer identity (`GeometryRegistry::shape_key`), mirroring
//! `BRep_TEdge`'s `(face → Geom2d_Curve)` map.
//!
//! Source: `BOPTools_AlgoTools2D.hxx/.cxx` and `BOPTools_AlgoTools2D_1.cxx`
//! (TKBO).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec2d};
use occt_core::precision::{ANGULAR, PCONFUSION, RESOLUTION};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::SurfaceKind;
use crate::pcurve_full::classify_surface_kind;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// OCCT `BOPTools_AlgoTools2D::IntermediatePoint` division constant:
/// `10·e^(−π) = 0.43213918`. Deliberately *not* the geometric midpoint — the
/// asymmetry keeps sample/intersection parameters away from symmetric
/// (ambiguous) positions.
const PAR_T: f64 = 0.43213918;

/// The p-curve of `edge` on `face`, returning the stored one when it already
/// exists (see [`curve_on_surface`]) and otherwise building it fresh
/// (`BOPTools_AlgoTools2D::Make2D`).
///
/// The build is delegated to [`crate::pcurve_full::make_pcurve_full`], which
/// gives analytic isoparametric projections on plane/cylinder/cone/sphere/torus
/// faces and a sampled degree-1 B-spline on general faces. The returned curve
/// is *not* attached to the edge — call
/// [`build_pcurve_for_edge_on_face`] for the build-and-attach behaviour.
pub fn make_2d(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
    if let Some(pc) = curve_on_surface(edge, face) {
        return Ok(pc);
    }
    crate::pcurve_full::make_pcurve_full(edge, face)
}

/// The stored p-curve of `edge` on `face`, if one has been attached to the
/// edge for that face (`BRep_Tool::CurveOnSurface`).
///
/// Unlike OCCT's `CurveOnSurface` this never builds — it is a pure read of the
/// [`GeometryRegistry`] side-table. Use [`make_2d`] or
/// [`build_pcurve_for_edge_on_face`] to obtain a p-curve when none is stored.
pub fn curve_on_surface(edge: &Edge, face: &Face) -> Option<Arc<dyn Curve2d>> {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcs = reg.edge_pcurves(&edge.0, face_key);
    // `BRep_Tool::CurveOnSurface` (`BRep_Tool.cxx:347-357`): for a
    // representation on a *closed* surface (a seam's two pcurves) a REVERSED
    // edge reads `PCurve2`, otherwise `PCurve1`.
    if edge.0.orientation() == crate::abs::Orientation::Reversed {
        if let Some(second) = pcs.get(1) {
            return Some(second.clone());
        }
    }
    pcs.into_iter().next()
}

/// `BRep_Tool::CurveOnSurface(edge, face, first, last)`: stored p-curve plus
/// the CurveOnSurface representation range (`BRep_GCurve::First/Last`).
/// Falls back to the 3D edge range, then the 2d curve domain.
///
/// When no COS is stored, `BRep_Tool.cxx:367-372` projects onto a plane via
/// `CurveOnPlane` and sets `First/Last` to the 3D edge `f,l` (`cxx:418-419`).
pub fn curve_on_surface_range(edge: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    if let Some(pc) = curve_on_surface(edge, face) {
        let face_key = GeometryRegistry::shape_key(&face.0);
        let reg = GeometryRegistry::global();
        if let Some((t1, t2)) = reg.pcurve_range(&edge.0, face_key) {
            if t1.is_finite() && t2.is_finite() {
                return Some((pc, t1, t2));
            }
        }
        let (t1, t2) = reg.edge_parameters(&edge.0);
        if t1.is_finite() && t2.is_finite() {
            return Some((pc, t1, t2));
        }
        let (a, b) = (pc.first_parameter(), pc.last_parameter());
        if a.is_finite() && b.is_finite() {
            return Some((pc, a, b));
        }
        return None;
    }
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    if classify_surface_kind(surf.as_ref()) != SurfaceKind::Plane {
        return None;
    }
    let pc = crate::pcurve::make_pcurve_on_face(edge, face).ok()?;
    let (t1, t2) = GeometryRegistry::global().edge_parameters(&edge.0);
    if t1.is_finite() && t2.is_finite() {
        Some((pc, t1, t2))
    } else {
        None
    }
}

/// `BRep_Tool::CurveOnSurface(edge, surface)`: a seam on a closed surface
/// returns `PCurve2` when the edge is `REVERSED` (`BRep_Tool.cxx:354-361`).
/// `ShapeAnalysis_Edge::PCurve(..., orient=true)` then swaps the range
/// (`ShapeAnalysis_Edge.cxx:201-206`).
pub fn curve_on_surface_oriented(
    edge: &Edge,
    face: &Face,
    orient: bool,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcs = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
    let reversed = edge.0.orientation().is_reversed();
    let pc = if pcs.len() >= 2 && reversed {
        pcs[1].clone()
    } else {
        pcs.first()?.clone()
    };
    let reg = GeometryRegistry::global();
    let (mut t1, mut t2) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or_else(|| reg.edge_parameters(&edge.0));
    if !(t1.is_finite() && t2.is_finite()) {
        t1 = pc.first_parameter();
        t2 = pc.last_parameter();
    }
    if !(t1.is_finite() && t2.is_finite()) {
        return None;
    }
    if orient && reversed {
        std::mem::swap(&mut t1, &mut t2);
    }
    Some((pc, t1, t2))
}

/// `ShapeBuild_Edge::ReplacePCurve` (`ShapeBuild_Edge.cxx:474-503`).
pub fn replace_pcurve(edge: &Edge, face: &Face, pcurve: Arc<dyn Curve2d>) {
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    // `ShapeBuild_Edge::ReplacePCurve` (`cxx:479-502`) reads COS `f,l`
    // before UpdateEdge, then `B.Range(edge, face, f, l)` restores it.
    let saved = reg.pcurve_range(&edge.0, face_key);
    let pcs = reg.edge_pcurves(&edge.0, face_key);
    if pcs.len() < 2 || Arc::ptr_eq(&pcs[0], &pcs[1]) {
        reg.set_edge_pcurve(&edge.0, face_key, pcurve);
    } else {
        let pair = if edge.0.orientation().is_reversed() {
            vec![pcs[0].clone(), pcurve]
        } else {
            vec![pcurve, pcs[1].clone()]
        };
        reg.set_edge_pcurves(&edge.0, face_key, pair);
    }
    if let Some((f, l)) = saved {
        reg.set_pcurve_range(&edge.0, face_key, f, l);
    }
}

/// Unit tangent of the edge's 3D curve at parameter `t`
/// (`BOPTools_AlgoTools2D::EdgeTangent`).
///
/// The tangent is normalized, and for a `REVERSED` edge the direction is
/// flipped to run with the edge orientation. Because the port's return type is
/// a `GpVec2d`, only the `(x, y)` components of the 3D tangent are returned
/// (the projection onto the XY plane); for edges lying in the XY plane this is
/// the exact direction. Errors when the edge has no 3D curve or the tangent is
/// degenerate (a zero-length derivative at `t`).
pub fn edge_tangent(edge: &Edge, t: f64) -> Result<GpVec2d, String> {
    let curve = GeometryRegistry::global()
        .edge_curve(&edge.0)
        .ok_or("edge_tangent: edge has no 3D curve")?;
    let (_, mut tau) = curve.d1(t);
    let mag = tau.magnitude();
    if mag <= RESOLUTION {
        return Err(format!("edge_tangent: degenerate tangent at t={t}"));
    }
    tau = tau.divided(mag);
    if edge.orientation().is_reversed() {
        tau = tau.reversed();
    }
    Ok(GpVec2d::new(tau.x(), tau.y()))
}

/// Whether the edge's p-curve on `face` is a surface isoline — a `u = const`
/// or `v = const` line in the face's parameter domain
/// (`BOPTools_AlgoTools2D::IsEdgeIsoline`).
///
/// The p-curve tangent at the edge-range midpoint is compared against the two
/// axis directions: a tangent parallel to `v` means `u` is constant (a
/// `u`-isoline), a tangent parallel to `u` means `v` is constant (a
/// `v`-isoline). The single boolean result is true when *either* holds. When no
/// p-curve is stored the curve is built on the fly via [`make_2d`] (not
/// attached).
pub fn is_edge_isoline(edge: &Edge, face: &Face) -> Result<bool, String> {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Ok(false);
    }
    let pc = match curve_on_surface(edge, face) {
        Some(pc) => pc,
        None => make_2d(edge, face)?,
    };
    let tm = 0.5 * (a + b);
    let (_, t) = pc.d1(tm);
    let sq = t.square_magnitude();
    if sq <= RESOLUTION {
        return Ok(false);
    }
    let n = t.divided(sq.sqrt());
    // sin(da) ~ da, when da -> 0.  A tangent parallel to (0,1) leaves u
    // constant; parallel to (1,0) leaves v constant.
    let dp_v = n.crossed(&GpVec2d::new(0.0, 1.0)).abs();
    let dp_u = n.crossed(&GpVec2d::new(1.0, 0.0)).abs();
    Ok(dp_v <= ANGULAR || dp_u <= ANGULAR)
}

/// Intermediate value in between `[t1, t2]`
/// (`BOPTools_AlgoTools2D::IntermediatePoint`): a point ~43.2% of the way from
/// `t1` to `t2`, matching `IntTools_Tools::IntermediatePoint`.
pub fn intermediate_point(t1: f64, t2: f64) -> f64 {
    (1.0 - PAR_T) * t1 + PAR_T * t2
}

/// Attach the already-built p-curve `curve` of `edge` on `face` to the edge's
/// geometry entry (`BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`).
///
/// The curve is stored in the [`GeometryRegistry`] side-table keyed by the
/// face's pointer identity, so a later [`curve_on_surface`] /
/// [`has_curve_on_surface`] round-trips it back. Errors when the edge has no
/// geometry entry to attach to.
pub fn attach_existing_pcurve(
    edge: &Edge,
    face: &Face,
    curve: Arc<dyn Curve2d>,
) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    if reg.edge_geom(&edge.0).is_none() {
        return Err("attach_existing_pcurve: edge has no geometry entry".into());
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    reg.set_edge_pcurve(&edge.0, face_key, curve);
    Ok(())
}

/// Build the p-curve of `edge` on `face` and attach it to the edge
/// (`BOPTools_AlgoTools2D::BuildPCurveForEdgeOnFace`).
///
/// When a p-curve is already stored it is returned unchanged; otherwise the
/// curve is built with [`make_2d`] and registered with
/// [`attach_existing_pcurve`], so subsequent [`curve_on_surface`] /
/// [`has_curve_on_surface`] calls see it.
pub fn build_pcurve_for_edge_on_face(
    edge: &Edge,
    face: &Face,
) -> Result<Arc<dyn Curve2d>, String> {
    if let Some(pc) = curve_on_surface(edge, face) {
        return Ok(pc);
    }
    let pc = make_2d(edge, face)?;
    attach_existing_pcurve(edge, face, pc.clone())?;
    Ok(pc)
}

/// The 3D point of the face's surface at parameter `(u, v)`
/// (`BRep_Tool::Surface(aF).Value(u, v)`).
///
/// Errors when the face has no registered surface.
pub fn point_on_surface(face: &Face, u: f64, v: f64) -> Result<GpPnt, String> {
    let surf = GeometryRegistry::global()
        .face_surface(&face.0)
        .ok_or("point_on_surface: face has no surface")?;
    Ok(surf.d0(u, v))
}

/// Whether the edge has a stored p-curve on the face
/// (`BOPTools_AlgoTools2D::HasCurveOnSurface`).
///
/// Follows the OCCT predicate: an edge whose parameter range is degenerate
/// (shorter than `Precision::PConfusion()`) never has a p-curve; otherwise the
/// decision is whether the [`GeometryRegistry`] holds an entry for
/// `(edge, face)`.
pub fn has_curve_on_surface(edge: &Edge, face: &Face) -> bool {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || (b - a).abs() < PCONFUSION {
        return false;
    }
    curve_on_surface(edge, face).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::Orientation;
    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::pcurve;
    use occt_core::gp::{GpAx3, GpPln, GpPnt};

    #[test]
    fn make_2d_agrees_with_pcurve_full_for_box_edge() {
        let b = unit_box();
        let edge = &b.edges[0]; // (0,0,0) → (1,0,0)
        let face = &b.faces[0]; // bottom (z = 0)
        let a = make_2d(edge, face).expect("make_2d");
        let full = crate::pcurve_full::make_pcurve_full(edge, face).expect("pcurve_full");
        assert_eq!(pcurve::pc_curve_kind(a.as_ref()), pcurve::CurveKind::Line);
        for i in 0..=8 {
            let t = i as f64 / 8.0;
            let pa = a.d0(t);
            let pf = full.d0(t);
            assert!(
                (pa.x() - pf.x()).abs() < 1e-9 && (pa.y() - pf.y()).abs() < 1e-9,
                "sample t={t}"
            );
        }
    }

    #[test]
    fn edge_tangent_constant_direction_on_box_straight_edge() {
        let b = unit_box();
        let edge = &b.edges[0]; // (0,0,0) → (1,0,0), tangent (1,0,0)
        let t0 = edge_tangent(edge, 0.25).expect("tangent t=0.25");
        let t1 = edge_tangent(edge, 0.75).expect("tangent t=0.75");
        // Constant direction along +X.
        assert!((t0.x() - 1.0).abs() < 1e-9 && (t0.y() - 0.0).abs() < 1e-9, "t0 {:?}", t0);
        assert!((t1.x() - 1.0).abs() < 1e-9 && (t1.y() - 0.0).abs() < 1e-9, "t1 {:?}", t1);
        assert_eq!(t0, t1);

        // A REVERSED edge runs against the stored curve direction.
        let mut rev = edge.clone();
        rev.0.set_orientation(Orientation::Reversed);
        let tr = edge_tangent(&rev, 0.5).expect("reversed tangent");
        assert!((tr.x() + 1.0).abs() < 1e-9 && (tr.y() - 0.0).abs() < 1e-9, "rev {:?}", tr);
    }

    #[test]
    fn is_edge_isoline_true_for_box_straight_edge() {
        let b = unit_box();
        let edge = &b.edges[0]; // (0,0,0) → (1,0,0) on the bottom face
        let face = &b.faces[0];
        assert!(is_edge_isoline(edge, face).expect("isoline check"));
    }

    #[test]
    fn is_edge_isoline_false_for_diagonal_edge_on_plane() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let edge = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        // The p-curve is the diagonal (t, t): its tangent is not axis-aligned,
        // so the edge is not a u- or v-isoline.
        assert!(!is_edge_isoline(&edge, &face).expect("isoline check"));
    }

    #[test]
    fn intermediate_point_matches_occt_constant() {
        let m = intermediate_point(0.0, 1.0);
        assert!((m - 0.43213918).abs() < 1e-12, "m {m}");
        let m2 = intermediate_point(10.0, 20.0);
        assert!((m2 - (10.0 + 0.43213918 * 10.0)).abs() < 1e-12, "m2 {m2}");
        // Agrees with IntTools_Tools::middle_point.
        let ref_ = crate::inttools_range::IntToolsTools::middle_point(0.0, 1.0);
        assert!((m - ref_).abs() < 1e-12);
    }

    #[test]
    fn attach_existing_pcurve_roundtrips_through_registry() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        assert!(!has_curve_on_surface(edge, face));
        assert!(curve_on_surface(edge, face).is_none());

        let pc = make_2d(edge, face).expect("pcurve");
        attach_existing_pcurve(edge, face, pc.clone()).expect("attach");

        assert!(has_curve_on_surface(edge, face));
        let back = curve_on_surface(edge, face).expect("curve on surface");
        let q1 = back.d0(0.5);
        let q2 = pc.d0(0.5);
        assert!(
            (q1.x() - q2.x()).abs() < 1e-12 && (q1.y() - q2.y()).abs() < 1e-12,
            "roundtrip {:?} vs {:?}",
            q1,
            q2
        );

        // make_2d now returns the stored curve instead of rebuilding.
        let again = make_2d(edge, face).expect("make_2d after attach");
        let q3 = again.d0(0.5);
        assert!((q3.x() - q1.x()).abs() < 1e-12 && (q3.y() - q1.y()).abs() < 1e-12);
    }

    #[test]
    fn build_pcurve_for_edge_on_face_attaches() {
        let b = unit_box();
        let edge = &b.edges[1]; // (1,0,0) → (1,1,0)
        let face = &b.faces[0]; // bottom (z = 0)
        assert!(!has_curve_on_surface(edge, face));
        let pc = build_pcurve_for_edge_on_face(edge, face).expect("built");
        assert!(has_curve_on_surface(edge, face));
        let back = curve_on_surface(edge, face).expect("stored");
        let q1 = pc.d0(0.0);
        let q2 = back.d0(0.0);
        assert!(
            (q1.x() - q2.x()).abs() < 1e-12 && (q1.y() - q2.y()).abs() < 1e-12,
            "attach mismatch {:?} vs {:?}",
            q1,
            q2
        );
    }

    #[test]
    fn point_on_surface_evaluates_face_d0() {
        let b = unit_box();
        let face = &b.faces[1]; // top (z = 1), (u,v) → (u, v, 1)
        let p = point_on_surface(face, 0.3, 0.7).expect("d0");
        assert!((p.x() - 0.3).abs() < 1e-9, "x {p:?}");
        assert!((p.y() - 0.7).abs() < 1e-9, "y {p:?}");
        assert!((p.z() - 1.0).abs() < 1e-9, "z {p:?}");
    }

    #[test]
    fn has_curve_on_surface_false_for_fresh_edges() {
        let b = unit_box();
        assert!(!has_curve_on_surface(&b.edges[2], &b.faces[2]));
        assert!(!has_curve_on_surface(&b.edges[2], &b.faces[1]));
    }
}
