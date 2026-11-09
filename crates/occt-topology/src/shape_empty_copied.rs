//! `TopoDS_Shape::EmptyCopied` for an edge and for a face.
//!
//! Source: `TopoDS_Shape.hxx:294-302`, `BRep_TEdge::EmptyCopy`
//! (`BRep_TEdge.cxx:104-129`), `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:36-43`)
//! and `ShapeBuild_Edge::CopyRanges` (`ShapeBuild_Edge.cxx:125`).
//!
//! TKBRep-level: the copy carries the `TShape` state (tolerance, the
//! `SameParameter` / `SameRange` / `Degenerated` flags, the 3D range and every
//! curve representation, pcurves included) and no components.

use crate::abs::ShapeType;
use crate::shape::TopoShape;

/// `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) for an EDGE: the
/// copy carries `BRep_TEdge::EmptyCopy`'s state (tolerance, `SameParameter` /
/// `SameRange` / `Degenerated`, the 3D range and every curve representation,
/// pcurves included) and no components.
///
/// `BRepFill_TrimShellCorner::MakeFacesNonSec` re-cuts a bound edge with
/// `TopoDS::Edge(aE1.EmptyCopied())` + `B.Begin/Add(E, V)` +
/// `Range(E, aparam1, aparam2)` (`BRepFill_TrimShellCorner.cxx:660-688`); the
/// pcurve of the source edge has to survive that copy, which is why this is
/// the shared `EmptyCopied` helper instead of a rebuild on the 3D curve.
pub fn empty_copied_edge(src: &TopoShape) -> TopoShape {
    let mut copy = TopoShape::new(ShapeType::Edge);
    copy.set_orientation(src.orientation());
    copy.set_location(src.location());
    copy_edge_geom(&copy, src);
    copy
}

/// `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) plus
/// `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:36-43`) for a face: the surface,
/// the location and the tolerance are carried over, and the surface is
/// re-registered so the rebuilt face still answers `BRep_Tool::Surface`.
pub fn empty_copied_face(src: &TopoShape) -> TopoShape {
    let mut copy = TopoShape::new(ShapeType::Face);
    copy.set_orientation(src.orientation());
    copy.set_location(src.location());
    copy_face_geom(&copy, src);
    copy
}

/// `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.hxx:294-302`) plus
/// `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) and `CopyRanges`
/// (`ShapeBuild_Edge.cxx:125`) for an edge: the tolerance, the
/// `SameParameter` / `SameRange` / `Degenerated` flags, the 3D range and a copy
/// of every curve representation (3D curve and pcurves), with no components.
pub fn copy_edge_geom(dst: &TopoShape, src: &TopoShape) {
    let (core, pcs) = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        (ts.edge_core().cloned(), ts.edge_pcurves().cloned())
    };
    let mut ts = dst.tshape.write().expect("poisoned TShape lock");
    if let Some(c) = core {
        let d = ts.edge_core_mut();
        d.curve = c.curve;
        // `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:104-129`) copies the whole
        // curve representation list, so a `BRep_Curve3D` holding a *null* curve
        // (a pcurve-only edge) stays null on the copy. Dropping the flag here
        // would make `BRep_Tool::Curve` answer the COS adaptor and take the
        // wrong mesh branch (`discretize_pcurves.rs`).
        d.curve3d_null = c.curve3d_null;
        d.first = c.first;
        d.last = c.last;
        d.tolerance = c.tolerance;
        d.same_parameter = c.same_parameter;
        d.same_range = c.same_range;
        d.degenerated = c.degenerated;
    }
    if let Some(p) = pcs {
        *ts.edge_pcurves_mut() = p;
    }
}

/// `BRep_TFace::EmptyCopy` (`BRep_TFace.cxx:37-43`) for a face: the surface and
/// the tolerance are carried over (the location is restored by the caller).
/// `myNaturalRestriction` deliberately stays at its constructor default
/// (`false`, `BRep_TFace.cxx:28-31`): `EmptyCopy` does not copy it either.
pub fn copy_face_geom(dst: &TopoShape, src: &TopoShape) {
    let core = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        ts.face_core().cloned()
    };
    if let Some(c) = &core {
        let mut ts = dst.tshape.write().expect("poisoned TShape lock");
        let d = ts.face_core_mut();
        d.surface = c.surface.clone();
        d.tolerance = c.tolerance;
    }
    // `BRep_Tool::CurveOnSurface` resolves the edge's representation through
    // the face's surface, so the rebuilt face must be registered against that
    // very surface; otherwise `repr_key` falls back to the face pointer and
    // every pcurve (keyed by the surface) becomes invisible on the copy.
    if let Some(surf) = core.as_ref().and_then(|c| c.surface.clone()) {
        crate::tgeometry::GeometryRegistry::global().register_face_surface(dst, &surf);
    }
}
