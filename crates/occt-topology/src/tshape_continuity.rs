//! Face-pair continuity storage — the `BRep_CurveOn2Surfaces` entries on
//! `BRep_TEdge::myCurves` that `BRep_Builder::Continuity` writes and
//! `BRep_Tool::Continuity` / `HasContinuity` / `MaxContinuity` read.
//!
//! Matched by `(surface, location)` relative to the edge
//! (`L.Predivided(E.Location())`), not by face shape ids. The public face-pair
//! API resolves each face to its surface handle and shape location
//! (`BRep_Tool::Surface(F, L)` ≈ face location; the port keeps no separate
//! `BRep_TFace` location).
//!
//! Source: `BRep_Builder.cxx:1012-1043` (`Continuity`),
//! `BRep_Builder.cxx:376-405` (`UpdateCurves` regularity arm),
//! `BRep_Tool.cxx:1170-1290` (`HasContinuity` / `Continuity` / `MaxContinuity`),
//! `BRep_CurveOn2Surfaces.cxx:58-65` (`IsRegularity`).

use std::sync::Arc;

use occt_core::toploc::TopLocLocation;
use occt_geom::Surface;

use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::tshape::CurveOn2SurfacesRepr;

/// Face surface + shape location used as the Continuity key
/// (`BRep_Tool::Surface(F, L)`, `BRep_Tool.cxx:76-81`).
fn face_surface_location(face: &Face) -> Option<(Arc<dyn Surface>, TopLocLocation)> {
    let surf = GeometryRegistry::global().face_surface(&face.0)?;
    Some((surf, face.0.location().clone()))
}

/// Locations relative to the edge (`L.Predivided(E.Location())`).
fn relative_locations(
    edge: &Edge,
    l1: &TopLocLocation,
    l2: &TopLocLocation,
) -> (TopLocLocation, TopLocLocation) {
    let eloc = edge.0.location();
    (l1.predivided(eloc), l2.predivided(eloc))
}

/// `BRep_Tool::Continuity(E, F1, F2)` (`BRep_Tool.cxx:1180-1188` / `:1223-1246`):
/// stored regularity, or `GeomAbs_C0` when no matching representation exists.
pub fn stored_continuity(edge: &Edge, f1: &Face, f2: &Face) -> u8 {
    let Some((s1, l1)) = face_surface_location(f1) else {
        return 0;
    };
    let Some((s2, l2)) = face_surface_location(f2) else {
        return 0;
    };
    continuity_on_surfaces(edge, &s1, &s2, &l1, &l2)
}

/// `BRep_Tool::Continuity(E, S1, S2, L1, L2)` (`BRep_Tool.cxx:1223-1246`).
pub fn continuity_on_surfaces(
    edge: &Edge,
    s1: &Arc<dyn Surface>,
    s2: &Arc<dyn Surface>,
    l1: &TopLocLocation,
    l2: &TopLocLocation,
) -> u8 {
    let (rl1, rl2) = relative_locations(edge, l1, l2);
    let ts = edge.0.tshape.read().expect("poisoned TShape lock");
    let Some(core) = ts.edge_core() else {
        return 0;
    };
    for r in &core.regularities {
        if r.is_regularity(s1, s2, &rl1, &rl2) {
            return r.continuity;
        }
    }
    0
}

/// `BRep_Tool::HasContinuity(E, F1, F2)` (`BRep_Tool.cxx:1170-1176`).
pub fn has_continuity_faces(edge: &Edge, f1: &Face, f2: &Face) -> bool {
    let Some((s1, l1)) = face_surface_location(f1) else {
        return false;
    };
    let Some((s2, l2)) = face_surface_location(f2) else {
        return false;
    };
    has_continuity_on_surfaces(edge, &s1, &s2, &l1, &l2)
}

/// `BRep_Tool::HasContinuity(E, S1, S2, L1, L2)` (`BRep_Tool.cxx:1195-1218`).
pub fn has_continuity_on_surfaces(
    edge: &Edge,
    s1: &Arc<dyn Surface>,
    s2: &Arc<dyn Surface>,
    l1: &TopLocLocation,
    l2: &TopLocLocation,
) -> bool {
    let (rl1, rl2) = relative_locations(edge, l1, l2);
    let ts = edge.0.tshape.read().expect("poisoned TShape lock");
    let Some(core) = ts.edge_core() else {
        return false;
    };
    core.regularities
        .iter()
        .any(|r| r.is_regularity(s1, s2, &rl1, &rl2))
}

/// `BRep_Tool::HasContinuity(E)` (`BRep_Tool.cxx:1253-1267`): any
/// `IsRegularity()` representation on the edge. The port's list is only
/// `BRep_CurveOn2Surfaces`; `BRep_CurveOnClosedSurface` regularity is UNPORTED.
pub fn has_continuity(edge: &Edge) -> bool {
    let ts = edge.0.tshape.read().expect("poisoned TShape lock");
    ts.edge_core()
        .map(|c| !c.regularities.is_empty())
        .unwrap_or(false)
}

/// `BRep_Tool::MaxContinuity(E)` (`BRep_Tool.cxx:1271-1290`).
pub fn max_continuity(edge: &Edge) -> u8 {
    let ts = edge.0.tshape.read().expect("poisoned TShape lock");
    let Some(core) = ts.edge_core() else {
        return 0;
    };
    core.regularities
        .iter()
        .map(|r| r.continuity)
        .max()
        .unwrap_or(0)
}

/// `BRep_Builder::Continuity(E, F1, F2, C)` (`BRep_Builder.cxx:1012-1021`).
/// `TShape::Locked` (`cxx:1033-1036`) is not stored.
pub fn set_stored_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: u8) {
    write_continuity(edge, f1, f2, cont);
}

/// Store the `GeomAbs_Shape` code directly (`GeomAbs_C0 = 0` .. `GeomAbs_CN = 6`).
pub fn write_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: u8) {
    let Some((s1, l1)) = face_surface_location(f1) else {
        return;
    };
    let Some((s2, l2)) = face_surface_location(f2) else {
        return;
    };
    write_continuity_on_surfaces(edge, s1, s2, l1, l2, cont);
}

/// `BRep_Builder::Continuity(E, S1, S2, L1, L2, C)` + `UpdateCurves`
/// (`BRep_Builder.cxx:1025-1043`, `:376-405`).
pub fn write_continuity_on_surfaces(
    edge: &Edge,
    s1: Arc<dyn Surface>,
    s2: Arc<dyn Surface>,
    l1: TopLocLocation,
    l2: TopLocLocation,
    cont: u8,
) {
    let (rl1, rl2) = relative_locations(edge, &l1, &l2);
    let mut ts = edge.0.tshape.write().expect("poisoned TShape lock");
    let core = ts.edge_core_mut();
    if let Some(existing) = core
        .regularities
        .iter_mut()
        .find(|r| r.is_regularity(&s1, &s2, &rl1, &rl2))
    {
        existing.continuity = cont;
        return;
    }
    core.regularities.push(CurveOn2SurfacesRepr {
        surface1: s1,
        surface2: s2,
        location1: rl1,
        location2: rl2,
        continuity: cont,
    });
}
