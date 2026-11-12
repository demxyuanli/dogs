//! `BRep_Tool::CurveOnSurface` and its parameter-range overload
//! (`BRep_Tool.cxx:300-420`), plus `BRep_Tool::CurveOnPlane`
//! (`BRep_Tool.cxx:379-449`).
//!
//! TKBRep-level: every arm reads the edge's `TShape` p-curve table (the
//! `BRep_GCurve` list, here `GeometryRegistry::edge_pcurves`) and, when nothing
//! is stored, the plane projection `CurveOnPlane` OCCT falls back to
//! (`BRep_Tool.cxx:367-372`). The geometry half of `CurveOnPlane`
//! (`GeomProjLib::ProjectOnPlane` + `ProjLib_ProjectedCurve`) lives in
//! `occt-geom::projlib_project_on_plane`.
//!
//! `BOPTools_AlgoTools2D` (`occt_topo::boptools_2d`) re-exports these so its
//! own call sites keep the `BOPTools_AlgoTools2D` surface.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// The stored p-curve of `edge` on `face`, if one has been attached to the
/// edge for that face (`BRep_Tool::CurveOnSurface`).
///
/// Unlike OCCT's `CurveOnSurface` this never builds — it is a pure read of the
/// [`GeometryRegistry`] side-table.
pub fn curve_on_surface(edge: &Edge, face: &Face) -> Option<Arc<dyn Curve2d>> {
    let reg = GeometryRegistry::global();
    let pcs = reg.edge_pcurves_on_face(&edge.0, &face.0);
    // `BRep_Tool::CurveOnSurface(E, F)` (`BRep_Tool.cxx:310-314`) reverses a
    // copy of the edge when the face is REVERSED, then the closed-surface arm
    // (`cxx:354-357`) returns `PCurve2` for that copy and `PCurve` otherwise.
    if seam_reads_pcurve2(edge, face) {
        if let Some(second) = pcs.get(1) {
            return Some(second.clone());
        }
    }
    pcs.into_iter().next()
}

/// `BRep_Tool::CurveOnSurface(E, F)` (`BRep_Tool.cxx:310-314`, `:354-357`):
/// PCurve2 is selected when exactly one of the edge and the face is REVERSED.
pub fn seam_reads_pcurve2(edge: &Edge, face: &Face) -> bool {
    let edge_rev = edge.0.orientation() == Orientation::Reversed;
    let face_rev = face.0.orientation() == Orientation::Reversed;
    edge_rev != face_rev
}

/// `BRep_Tool::CurveOnSurface(edge, face, first, last)`: stored p-curve plus
/// the CurveOnSurface representation range (`BRep_GCurve::First/Last`).
/// Falls back to the 3D edge range, then the 2d curve domain.
///
/// When no COS is stored, `BRep_Tool.cxx:367-372` projects onto a plane via
/// `CurveOnPlane` and sets `First/Last` to the 3D edge `f,l` (`cxx:418-419`).
pub fn curve_on_surface_range(edge: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    if let Some(pc) = curve_on_surface(edge, face) {
        let reg = GeometryRegistry::global();
        if let Some((t1, t2)) = reg.pcurve_range_on_face(&edge.0, &face.0) {
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
    let pc = pcurve_on_plane_face(edge, face)?;
    let (t1, t2) = GeometryRegistry::global().edge_parameters(&edge.0);
    if t1.is_finite() && t2.is_finite() {
        Some((pc, t1, t2))
    } else {
        None
    }
}

/// `BRep_Tool.cxx:367-372`: with no COS stored, `BRep_Tool::CurveOnSurface`
/// builds one with `CurveOnPlane`, whose plane arm is
/// `ProjLib::MakePCurveOnFace`'s `GeomAbs_Plane` branch
/// (`ProjLib.cxx:217-232`). A face whose adaptor type is not `GeomAbs_Plane`
/// has no implicit p-curve.
///
/// The plane test is the adaptor type (`GeomAdaptor_Surface.cxx:422-513`), not
/// a normal-sampling heuristic: a trimmed `Geom_Plane` answers true and a
/// geometrically flat B-spline patch (`GeomAbs_BSplineSurface`) answers false.
fn pcurve_on_plane_face(edge: &Edge, face: &Face) -> Option<Arc<dyn Curve2d>> {
    let reg = GeometryRegistry::global();
    let surf = reg.face_surface(&face.0)?;
    if !occt_geom::surface_kind::is_plane_type(surf.as_ref()) {
        return None;
    }
    let curve = reg.edge_curve(&edge.0)?;
    let relative = edge.0.location().predivided(face.0.location());
    let curve = if relative.is_identity() {
        curve
    } else {
        Arc::from(curve.transformed(&relative.transformation()))
    };
    occt_geom::projlib_project_on_plane::curve_on_plane(curve.as_ref(), surf.as_ref())
}
