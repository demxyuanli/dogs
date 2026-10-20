//! Geometry attachment — side-table storing analytic geometry for TShapes.
//!
//! OCCT stores geometry *inside* the TShape sub-classes (`BRep_TVertex`,
//! `BRep_TEdge`, `BRep_TFace`). This port keeps `TShape` itself light
//! (Clone + Debug, no trait objects) and attaches geometry via a process-wide
//! registry keyed by TShape identity (pointer address). `BRepBuilder` writes
//! geometry here; `BRepTool` reads it back. This mirrors `BRep_Tool`'s role of
//! separating topological structure from geometric content.
//!
//! The `TShape` address is stable for the lifetime of the `Arc`, so pointer
//! keys are valid until the shape is dropped. `clear_shape` (called on Drop
//! hooks, or manually) releases the entry.
//!
//! Source: `BRep_TVertex` / `BRep_TEdge` / `BRep_TFace` (TKBRep).

use std::collections::HashMap;
use std::sync::{Arc, RwLock, OnceLock};

use occt_core::gp::GpPnt;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::shape::TopoShape;
use crate::tshape::{EdgePcurves, OrderedSlots, VertexGeomCore};

/// Monotonic sequence number used by the in-tree pcurve/ReShape diagnostics to
/// order an attach/remove against a recorded `Replace`.
pub fn next_pcseq() -> usize {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static SEQ: AtomicUsize = AtomicUsize::new(0);
    SEQ.fetch_add(1, Ordering::Relaxed)
}

/// Vertex geometry — a 3D point and a tolerance. (BRep_TVertex)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexGeom {
    pub point: GpPnt,
    pub tolerance: f64,
}

/// PARKED (t313) `BRep_TEdge` flags for a freshly built edge.
///
/// `BRep_TEdge::BRep_TEdge()` starts with `SameParameter(true)`/`SameRange(true)`
/// (`BRep_TEdge.cxx:31-38`) and `EmptyCopy` carries both flags with the edge
/// (`cxx:124-127`), so seeding `true` matches OCCT for a freshly built edge. A
/// fix tool clears the flags only when it rewrites the curve/pcurve
/// (`ShapeFix_Edge.cxx:781-782` FixReversed2d, `ShapeBuild_Edge.cxx:326-327`
/// CopyRanges, `ShapeFix.cxx:129-133` for `enforce`, `BRepLib.cxx:945`,
/// `ShapeConstruct.cxx:449-450`), and the STEP importer clears them only for
/// COMPOSITE_CURVE edges (`StepToTopoDS_TranslateCompositeCurve.cxx:267`).
///
/// t314 ported the `!wasSP` arm of `ShapeFix_Edge::FixSameParameter`
/// (`shhealing/wire_fix.rs`), so this flag no longer hides a missing branch. Measured
/// with a temporary probe over all 16 `data/*.step` inputs of
/// `export_data_obj`: no edge is cleared, the `!wasSP` arm never runs, and no
/// input contains a COMPOSITE_CURVE. It stays dormant under FromSTEP.FixShape
/// (`FixEdgeSameParameterMode: 0` -> `ShapeFix_Wire.cxx:953` is off;
/// `FixSameParameterMode: -1` -> `ShapeFix_Shape.cxx:259-261` passes
/// `enforce = false`, so `ShapeFix.cxx:129-133` does not clear the flag;
/// `read.stdsameparameter.mode` defaults to 0, `XSAlgo.cxx:46`).
///
/// The remaining `Shape-1` face `f34` box (`v = 5.346039` / `Zmax = 75.346039`
/// against OCCT `5.283071`/`75.283071`) is not this arm: our cylinder pcurve of
/// that edge is the 10 pole B-spline built by
/// `fix_add_pcurve` -> `pcurve_full::project_curve_on_surface_perform`
/// (`0.196` same-parameter deviation), while OCCT's comes from
/// `XSAlgo_ShapeProcessor::CheckPCurve`'s `FixAddPCurve` tail
/// (`XSAlgo_ShapeProcessor.cxx:448-505`), parked here as
/// `CHECK_PCURVE_REPROJECT`.
const EDGE_NEW_IS_SAME_PARAMETER: bool = true;

/// Edge geometry — an underlying curve and a parameter range. (BRep_TEdge)
pub struct EdgeGeom {
    pub curve: Arc<dyn Curve>,
    /// `true` when the edge's `BRep_Curve3D` representation holds a null curve
    /// (see [`crate::tshape::EdgeGeomCore::curve3d_null`]); `curve` then carries
    /// what `BRepAdaptor_Curve(E, F)` builds from the pcurve, and
    /// [`GeometryRegistry::edge_curve`] answers `None`, as `BRep_Tool::Curve`
    /// does.
    pub curve3d_null: bool,
    pub first: f64,
    pub last: f64,
    pub tolerance: f64,
    pub same_parameter: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// `BRep_GCurve` CurveOnSurface representations, keyed by the surface data
    /// pointer (`GeometryRegistry::repr_key`): `BRep_Tool::CurveOnSurface`
    /// finds a representation by the `Geom_Surface` handle, so faces that share
    /// a surface share one representation, as in OCCT. A seam edge of a
    /// periodic surface carries *two* pcurves on the same representation (one
    /// per side of the seam), stored in forward-then-reversed order.
    ///
    /// Traversal follows `BRep_TEdge::myCurves`' `NCollection_List` order
    /// (`OrderedSlots`), not the hash order.
    pub pcurves: OrderedSlots<Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve::First/Last` of each CurveOnSurface representation.
    /// `UpdateEdge` copies the 3D range; `BRep_Builder::Range(edge, face)`
    /// can then clamp it (`StepToTopoDS_TranslateEdgeLoop::CheckPCurves`).
    pub pcurve_ranges: OrderedSlots<(f64, f64)>,
}

impl EdgeGeom {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve,
            curve3d_null: false,
            first,
            last,
            tolerance: 0.0,
            same_parameter: EDGE_NEW_IS_SAME_PARAMETER,
            same_range: true,
            degenerated: false,
            pcurves: OrderedSlots::default(),
            pcurve_ranges: OrderedSlots::default(),
        }
    }
    pub fn curve(&self) -> Arc<dyn Curve> { self.curve.clone() }
    pub fn parameters(&self) -> (f64, f64) { (self.first, self.last) }

    /// Attach the (single) pcurve of this edge on the face identified by
    /// `face_key`, replacing any previously attached pcurves.
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>) {
        self.pcurves.set(face_key, vec![c.clone()]);
        self.init_pcurve_range(face_key, Some(c.as_ref()));
    }

    /// Replace the pcurves of this edge on the face (one for a normal edge, two
    /// in forward-then-reversed order for a seam edge).
    pub fn set_pcurves(&mut self, face_key: usize, cs: Vec<Arc<dyn Curve2d>>) {
        let first = cs.first().cloned();
        self.pcurves.set(face_key, cs);
        self.init_pcurve_range(face_key, first.as_deref());
    }

    /// `UpdateCurves` (`BRep_Builder.cxx:149-164`): new COS range is the 3D
    /// range when finite, else the pcurve's own `[First, Last]`.
    fn init_pcurve_range(&mut self, face_key: usize, pc: Option<&dyn Curve2d>) {
        let (f, l) = if self.first.is_finite() && self.last.is_finite() {
            (self.first, self.last)
        } else if let Some(c) = pc {
            (c.first_parameter(), c.last_parameter())
        } else {
            return;
        };
        self.pcurve_ranges.insert(face_key, (f, l));
    }

    /// The first pcurve on the face (the single pcurve of a normal edge).
    pub fn pcurve(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.pcurves.get(&face_key).and_then(|v| v.first().cloned())
    }

    /// All pcurves on the face (two, in forward-then-reversed order, for a seam
    /// edge).
    pub fn pcurves(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        self.pcurves.get(&face_key).cloned().unwrap_or_default()
    }
}

/// Face geometry — an underlying surface and a tolerance. (BRep_TFace)
pub struct FaceGeom {
    pub surface: Arc<dyn Surface>,
    pub tolerance: f64,
    pub natural_restriction: bool,
}

impl FaceGeom {
    /// `BRep_TFace` default: `myNaturalRestriction = false`
    /// (`BRep_TFace.cxx:30`; neither `BRep_Builder::MakeFace` overload sets it,
    /// `BRep_Builder.cxx:500-560`). Only an *unbounded* face turns it on
    /// (`BRepPrim_FaceBuilder.cxx:168`, `BOPAlgo_BuilderFace.cxx:407-409`), and
    /// `BRepToIGES_BRShell.cxx:381` reads it as `isWholeSurface`.
    pub fn new(surface: Arc<dyn Surface>) -> Self {
        // `BRep_TFace` default tolerance is `RealEpsilon()` (`BRep_TFace.cxx:29`),
        // not 0: `BRepAdaptor_Surface::Tolerance()` (`BRepAdaptor_Surface.cxx:92`)
        // and `BOPAlgo_PaveFiller::ToleranceFF` read it straight off the face.
        Self {
            surface,
            tolerance: occt_core::precision::REAL_EPSILON,
            natural_restriction: false,
        }
    }
    pub fn surface(&self) -> Arc<dyn Surface> { self.surface.clone() }
}

/// Process-wide geometry side-table keyed by `TShape` pointer identity.
pub struct GeometryRegistry {
    // T-25 end phase: the three geometry maps are gone — pcurves, edge,
    // vertex and face geometry all live on their own `TShape` now. Only the
    // shape ids and the face->surface table remain.
    /// `TShape` address -> the `id` of the shape that registered it. A drop
    /// removes the geometry only when the ids match, so a stale drop cannot
    /// erase a later shape's entries at a reused address.
    ids: RwLock<HashMap<usize, u64>>,
    /// T-25 end phase: `face_key` -> the face's `Geom_Surface`, kept for the
    /// `BRep_Tool::CurveOnSurface` identity fallback while the full `faces`
    /// map is on its way out. Surfaces never change, so this never needs
    /// resynchronising (unlike the tolerances in `FaceGeom`).
    face_surfaces: RwLock<HashMap<usize, Arc<dyn Surface>>>,
    /// Surface data pointer -> the surface handle. `BRep_GCurve`
    /// representations are keyed by the surface handle in OCCT, so the port
    /// keys an edge's pcurves by this pointer; this map recovers the handle
    /// when a representation has to be returned.
    surface_by_ptr: RwLock<HashMap<usize, Arc<dyn Surface>>>,
}

/// Registry key: the address of the `TShape` stored inside the shared
/// `RwLock`. This is stable for the lifetime of the `Arc` (the `TShape` never
/// moves) and matches what `TShape::drop` can recover from `&self`, so entries
/// are reliably removed when the last handle to a shape is dropped.
///
/// `TopoShape::shape_addr` caches the same value at construction, so this is
/// lock-free.
fn key(s: &TopoShape) -> usize {
    s.shape_addr()
}

/// The shape's process-unique id (see `TShape::id`), cached on the view.
fn shape_id(s: &TopoShape) -> u64 {
    s.shape_id()
}


mod registry_core;
mod registry_faces;
mod registry_lifecycle;
mod registry_pcurves;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
