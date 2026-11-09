//! Face-pair continuity storage — the `BRep_TEdge` continuity list that
//! `BRep_Builder::Continuity` writes and `BRep_Tool::Continuity` reads.
//!
//! OCCT keeps the list on the edge's `TShape` (`BRep_TEdge::ChangeCurves`'
//! `BRep_GCurve` pairs, matched by `(surface, location)`); the port keys a
//! side-table by `(edge, face_a, face_b)` shape ids with `face_a <= face_b`,
//! so either face order answers the same pair.
//!
//! Source: `BRep_Builder.cxx:1012-1036` (`Continuity`), `BRep_Tool.cxx:155-160`
//! (`HasContinuity`), `BRep_Tool.cxx` (`Continuity`).

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `(edge, face_a, face_b)` with `face_a <= face_b`. `BRep_Tool::Continuity`
/// matches either face order (`BRep_Tool.cxx` continuity query).
fn pair_key(edge: &Edge, f1: &Face, f2: &Face) -> (usize, usize, usize) {
    let e = GeometryRegistry::shape_key(&edge.0);
    let mut a = GeometryRegistry::shape_key(&f1.0);
    let mut b = GeometryRegistry::shape_key(&f2.0);
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    (e, a, b)
}

fn continuity_store() -> &'static Mutex<HashMap<(usize, usize, usize), u8>> {
    static STORE: OnceLock<Mutex<HashMap<(usize, usize, usize), u8>>> = OnceLock::new();
    STORE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Stored face-pair continuity, or `GeomAbs_C0` when the edge was never encoded.
pub fn stored_continuity(edge: &Edge, f1: &Face, f2: &Face) -> u8 {
    let guard = continuity_store()
        .lock()
        .unwrap_or_else(|err| err.into_inner());
    guard.get(&pair_key(edge, f1, f2)).copied().unwrap_or(0)
}

/// `BRep_Tool::HasContinuity(E)` (`BRep_Tool.cxx:155-160`):
/// `E.TShape()->HasContinuity()`, i.e. whether `BRep_Builder::Continuity` ever
/// stored a pair for this edge. The port's store is keyed by
/// `(edge, face_a, face_b)`, so the flag is "any pair bound to this edge".
pub fn has_continuity(edge: &Edge) -> bool {
    let e = GeometryRegistry::shape_key(&edge.0);
    let guard = continuity_store()
        .lock()
        .unwrap_or_else(|err| err.into_inner());
    guard.keys().any(|k| k.0 == e)
}

/// `BRep_Builder::Continuity(E, F1, F2, C)` (`BRep_Builder.cxx:1012-1021`).
/// `TShape::Locked` (`cxx:1034-1036`) is not stored.
pub fn set_stored_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: u8) {
    write_continuity(edge, f1, f2, cont);
}

/// Store the `GeomAbs_Shape` code directly (`GeomAbs_C0 = 0` .. `GeomAbs_CN = 6`).
pub fn write_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: u8) {
    let mut guard = continuity_store()
        .lock()
        .unwrap_or_else(|err| err.into_inner());
    guard.insert(pair_key(edge, f1, f2), cont);
}
