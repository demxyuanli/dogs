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
use crate::tshape::{EdgePcurves, FaceGeomCore, VertexGeomCore};

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
    pub pcurves: HashMap<usize, Vec<Arc<dyn Curve2d>>>,
    /// `BRep_GCurve::First/Last` of each CurveOnSurface representation.
    /// `UpdateEdge` copies the 3D range; `BRep_Builder::Range(edge, face)`
    /// can then clamp it (`StepToTopoDS_TranslateEdgeLoop::CheckPCurves`).
    pub pcurve_ranges: HashMap<usize, (f64, f64)>,
}

impl EdgeGeom {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve,
            first,
            last,
            tolerance: 0.0,
            same_parameter: EDGE_NEW_IS_SAME_PARAMETER,
            same_range: true,
            degenerated: false,
            pcurves: HashMap::new(),
            pcurve_ranges: HashMap::new(),
        }
    }
    pub fn curve(&self) -> Arc<dyn Curve> { self.curve.clone() }
    pub fn parameters(&self) -> (f64, f64) { (self.first, self.last) }

    /// Attach the (single) pcurve of this edge on the face identified by
    /// `face_key`, replacing any previously attached pcurves.
    pub fn set_pcurve(&mut self, face_key: usize, c: Arc<dyn Curve2d>) {
        self.pcurves.insert(face_key, vec![c.clone()]);
        self.init_pcurve_range(face_key, Some(c.as_ref()));
    }

    /// Replace the pcurves of this edge on the face (one for a normal edge, two
    /// in forward-then-reversed order for a seam edge).
    pub fn set_pcurves(&mut self, face_key: usize, cs: Vec<Arc<dyn Curve2d>>) {
        let first = cs.first().cloned();
        self.pcurves.insert(face_key, cs);
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
    pub fn get_pcurve(&self, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.pcurves.get(&face_key).and_then(|v| v.first().cloned())
    }

    /// All pcurves on the face (two, in forward-then-reversed order, for a seam
    /// edge).
    pub fn get_pcurves(&self, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
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
        Self { surface, tolerance: 0.0, natural_restriction: false }
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
fn key(s: &TopoShape) -> usize {
    let lock = s.tshape.read().expect("poisoned TShape lock");
    std::ptr::addr_of!(*lock) as usize
}

/// The shape's process-unique id (see `TShape::id`).
fn shape_id(s: &TopoShape) -> u64 {
    s.tshape.read().expect("poisoned TShape lock").id
}

impl GeometryRegistry {
    /// The shared registry. Shapes and geometry live for the whole process,
    /// matching OCCT's reference-counted Handle model.
    pub fn global() -> &'static GeometryRegistry {
        static REGISTRY: OnceLock<GeometryRegistry> = OnceLock::new();
        REGISTRY.get_or_init(|| GeometryRegistry {
            ids: RwLock::new(HashMap::new()),
            face_surfaces: RwLock::new(HashMap::new()),
            surface_by_ptr: RwLock::new(HashMap::new()),
        })
    }

    // ---- vertices ----

    pub fn set_vertex(&self, s: &TopoShape, geom: VertexGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25 batch 3: the vertex geometry lives on the shape itself now.
        let mut ts = s.tshape.write().unwrap();
        let c = ts.vertex_core_mut();
        c.point = geom.point;
        c.tolerance = geom.tolerance;
        drop(ts);
    }

    pub fn vertex_geom(&self, s: &TopoShape) -> Option<VertexGeom> {
        // T-25 batch 3: read off the vertex's own `TShape`; the side table is
        // only a fallback for shapes registered before this batch.
        // T-25 end phase: the vertex geometry comes off its own `TShape`.
        let c = s.tshape.read().unwrap().vertex_core().copied()?;
        Some(VertexGeom { point: c.point, tolerance: c.tolerance })
    }

    /// The vertex point, or the origin when unregistered (placeholder fallback).
    pub fn vertex_point(&self, s: &TopoShape) -> GpPnt {
        self.vertex_geom(s).map(|g| g.point).unwrap_or_else(GpPnt::zero)
    }

    pub fn vertex_tolerance(&self, s: &TopoShape) -> f64 {
        self.vertex_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    // ---- edges ----

    pub fn set_edge(&self, s: &TopoShape, geom: EdgeGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25: drain the carrier onto the edge's own `TShape`. A geom that
        // carries no pcurve must NOT clear the slot (clearing is the explicit
        // `remove_pcurves_on_surface`, `ShapeBuild_Edge::RemovePCurve`).
        if !(geom.pcurves.is_empty() && geom.pcurve_ranges.is_empty()) {
            s.tshape.write().unwrap().edge_pcurves = Some(EdgePcurves {
                curves: geom.pcurves.clone(),
                ranges: geom.pcurve_ranges.clone(),
            });
        }
        {
            let mut ts = s.tshape.write().unwrap();
            let c = ts.edge_core_mut();
            c.curve = Some(geom.curve.clone());
            c.first = geom.first;
            c.last = geom.last;
            c.tolerance = geom.tolerance;
            c.same_parameter = geom.same_parameter;
            c.same_range = geom.same_range;
            c.degenerated = geom.degenerated;
        }
    }

    /// `BRep_Builder::Range(E, First, Last)`: set the edge's parameter range.
    /// `BRepPrim_OneAxis` relies on it for the degenerate pole edges
    /// (`SetParameters(ETOP/EBOTTOM, …, 0., myAngle)`,
    /// `BRepPrim_OneAxis.cxx:407`, `:418`): the 3D curve is absent but the range
    /// still spans the full period so the edge's pcurve is a full u-isoline.
    pub fn set_edge_range(&self, s: &TopoShape, first: f64, last: f64) {
        let mut ts = s.tshape.write().unwrap();
        let c = ts.edge_core_mut();
        c.first = first;
        c.last = last;
    }

    pub fn edge_geom(&self, s: &TopoShape) -> Option<EdgeGeom> {
        // T-25 end phase: everything comes off the edge's own `TShape`; the
        // side table no longer stores edge geometry.
        let ts = s.tshape.read().unwrap();
        let c = ts.edge_core()?;
        Some(EdgeGeom {
            // An `EdgeGeom` always carries a 3D curve; an edge registered
            // without one is not representable and is reported as absent.
            curve: c.curve.clone()?,
            first: c.first,
            last: c.last,
            tolerance: c.tolerance,
            same_parameter: c.same_parameter,
            same_range: c.same_range,
            degenerated: c.degenerated,
            pcurves: ts.edge_pcurves().map(|p| p.curves.clone()).unwrap_or_default(),
            pcurve_ranges: ts.edge_pcurves().map(|p| p.ranges.clone()).unwrap_or_default(),
        })
    }

    /// The underlying curve handle (clone of the Arc).
    pub fn edge_curve(&self, s: &TopoShape) -> Option<Arc<dyn Curve>> {
        self.edge_geom(s).map(|g| g.curve)
    }

    pub fn edge_parameters(&self, s: &TopoShape) -> (f64, f64) {
        // `BRep_TEdge` always carries `myRange`; the 3D curve is a *separate*
        // CurveRepresentation (`BRep_Curve3D`) and may be absent.
        // `ShapeFix_ComposeShell::SplitByLine` builds exactly such an edge
        // (`cxx:2061-2070`: `MakeEdge` + two pcurves + `Range`), so the range
        // must be read from the edge core even when `edge_geom` reports no
        // curve. Returning `(-inf, inf)` here made `set_pcurves` /
        // `set_pcurve_range` consumers store an infinite COS window and the
        // seam edge never meshed.
        let ts = s.tshape.read().expect("poisoned TShape lock");
        match ts.edge_core() {
            Some(c) => (c.first, c.last),
            None => (f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    pub fn edge_tolerance(&self, s: &TopoShape) -> f64 {
        self.edge_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    pub fn same_parameter(&self, s: &TopoShape) -> bool {
        self.edge_geom(s).map(|g| g.same_parameter).unwrap_or(false)
    }

    pub fn is_degenerated_edge(&self, s: &TopoShape) -> bool {
        self.edge_geom(s).map(|g| g.degenerated).unwrap_or(false)
    }

    /// `BRep_Builder::Degenerated(E, D)` (`BRep_Builder.cxx:1073-1085`): set the
    /// degenerated flag. OCCT also drops the 3D curve when `D` is true
    /// (`UpdateCurves(TE->ChangeCurves(), occ::handle<Geom_Curve>(), ...)`,
    /// `cxx:1082-1084`). `EdgeGeom::curve` here is a non-nullable
    /// `Arc<dyn Curve>`, so a caller that needs OCCT's "no 3D curve" state
    /// stores what `BRepAdaptor_Curve(edge, face)` would build from the pcurve
    /// instead (`Adaptor3d_CurveOnSurface`).
    pub fn set_degenerated(&self, s: &TopoShape, v: bool) {
        let mut ts = s.tshape.write().unwrap();
        ts.edge_core_mut().degenerated = v;
    }

    // ---- edge p-curves ----

    /// The pcurve of edge `s` on the face identified by `face_key` (see
    /// `shape_key`), if one has been attached. For a seam edge this is the
    /// forward pcurve; use [`GeometryRegistry::edge_pcurves`] to get both.
    ///
    /// When no entry exists for `face_key`, a pcurve attached for another face
    /// that shares the same surface handle is returned (`BRep_Tool::CurveOnSurface`
    /// keys representations by `Geom_Surface`, not by face TShape).
    pub fn edge_pcurve(&self, s: &TopoShape, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.edge_pcurves(s, face_key).into_iter().next()
    }

    /// All pcurves of edge `s` on the face identified by `face_key`
    /// (forward-then-reversed for a seam edge, one for a normal edge).
    pub fn edge_pcurves(&self, s: &TopoShape, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        let key = self.repr_key(face_key);
        let ts = s.tshape.read().unwrap();
        ts.edge_pcurves()
            .map(|p| p.get_pcurves(key))
            .unwrap_or_default()
    }

    /// Attach a pcurve to edge `s` for the face identified by `face_key`.
    /// Mirrors `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`.
    pub fn set_edge_pcurve(&self, s: &TopoShape, face_key: usize, curve: Arc<dyn Curve2d>) {
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurve(key, curve, first, last);
    }

    /// Replace the pcurves of edge `s` on the face identified by `face_key`.
    /// Mirrors the seam overload `BRep_Builder::UpdateEdge(edge, c1, c2, face)`.
    pub fn set_edge_pcurves(&self, s: &TopoShape, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurves(key, curves, first, last);
    }

    /// `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:360-413`) at the
    /// `BRep_TEdge` level: copy every CurveOnSurface representation (pcurve
    /// list + range) of `from` onto `to`, keyed by the surface pointer.
    ///
    /// Deliberately independent of `edge_geom`: OCCT's CurveRepresentation
    /// list exists for edges that carry **only** pcurves.
    /// `ShapeFix_ComposeShell::SplitByLine` (`ShapeFix_ComposeShell.cxx:2061-2070`)
    /// creates exactly such an edge (no 3D curve, two pcurves + range), and
    /// `DispatchWires` copies it with `ShapeBuild_Edge::Copy`; routing that
    /// through `edge_geom` dropped the pcurves and the resulting face had
    /// mesh-less edges.
    pub fn copy_edge_pcurve_slots(&self, to: &TopoShape, from: &TopoShape) {
        if Arc::ptr_eq(&from.tshape, &to.tshape) {
            return;
        }
        let (curves, ranges) = {
            let fts = from.tshape.read().expect("poisoned TShape lock");
            match fts.edge_pcurves() {
                Some(g) => (g.curves.clone(), g.ranges.clone()),
                None => return,
            }
        };
        let mut tts = to.tshape.write().expect("poisoned TShape lock");
        let slot = tts.edge_pcurves_mut();
        for (k, v) in curves {
            slot.curves.insert(k, v);
        }
        for (k, v) in ranges {
            slot.ranges.insert(k, v);
        }
    }

    /// `BRep_Builder::Range(edge, face, first, last)` (`BRep_Builder.cxx:1121`).
    pub fn set_pcurve_range(&self, s: &TopoShape, face_key: usize, first: f64, last: f64) {
        let key = self.repr_key(face_key);
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .ranges
            .insert(key, (first, last));
    }

    /// COS representation `[First, Last]` for `face_key`, with the same
    /// same-surface fallback as [`GeometryRegistry::edge_pcurve`].
    pub fn pcurve_range(&self, s: &TopoShape, face_key: usize) -> Option<(f64, f64)> {
        let key = self.repr_key(face_key);
        let ts = s.tshape.read().unwrap();
        ts.edge_pcurves()?.ranges.get(&key).copied()
    }

    /// `BRep_Builder::SameRange`.
    pub fn set_same_range(&self, s: &TopoShape, value: bool) {
        s.tshape.write().unwrap().edge_core_mut().same_range = value;
    }

    /// `BRep_Builder::SameParameter`.
    pub fn set_same_parameter(&self, s: &TopoShape, value: bool) {
        s.tshape.write().unwrap().edge_core_mut().same_parameter = value;
    }

    /// `BRep_Builder::UpdateEdge` tolerance write.
    pub fn set_edge_tolerance(&self, s: &TopoShape, tol: f64) {
        s.tshape.write().unwrap().edge_core_mut().tolerance = tol;
    }

    /// CurveOnSurface representations: surface, pcurve, COS `[first, last]`.
    pub fn edge_pcurve_reps(
        &self,
        s: &TopoShape,
    ) -> Vec<(Arc<dyn Surface>, Arc<dyn Curve2d>, f64, f64)> {
        let by_ptr = self.surface_by_ptr.read().unwrap();
        // T-25 end phase: pcurves, their COS ranges and the 3D range all come
        // off the edge's own `TShape`.
        let (slot, three_d) = {
            let ts = s.tshape.read().unwrap();
            let Some(core) = ts.edge_core() else {
                return Vec::new();
            };
            (ts.edge_pcurves().cloned(), (core.first, core.last))
        };
        let Some(g) = slot else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (&fk, pcs) in &g.curves {
            let Some(surf) = by_ptr.get(&fk) else {
                continue;
            };
            let (a, b) = g.ranges.get(&fk).copied().unwrap_or(three_d);
            for pc in pcs {
                out.push((surf.clone(), pc.clone(), a, b));
            }
        }
        out
    }

    /// `ShapeBuild_Edge::RemovePCurve` — drop every pcurve of `edge` on the
    /// same `Geom_Surface` handle as `face`.
    pub fn remove_pcurves_on_surface(&self, edge: &TopoShape, face: &TopoShape) {
        let Some(want) = self.face_surface(face) else {
            return;
        };
        let want_key = Arc::as_ptr(&want) as *const () as usize;
        let mut ts = edge.tshape.write().unwrap();
        if let Some(g) = ts.edge_pcurves.as_mut() {
            g.curves.remove(&want_key);
            g.ranges.remove(&want_key);
        }
    }

    /// Stable pointer-identity key for a shape, usable as a `HashMap` key.
    /// Same value as the internal `key` used to index the side-table.
    pub fn shape_key(s: &TopoShape) -> usize {
        key(s)
    }

    // ---- faces ----

    pub fn set_face(&self, s: &TopoShape, geom: FaceGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25 batch 4: the face geometry lives on the shape itself now.
        {
            let mut ts = s.tshape.write().unwrap();
            let c = ts.face_core_mut();
            c.surface = Some(geom.surface.clone());
            c.tolerance = geom.tolerance;
            c.natural_restriction = geom.natural_restriction;
        }
        self.face_surfaces.write().unwrap().insert(k, geom.surface.clone());
        self.surface_by_ptr
            .write()
            .unwrap()
            .insert(Arc::as_ptr(&geom.surface) as *const () as usize, geom.surface);
    }

    /// Canonical `BRep_GCurve` representation key: OCCT keys a CurveOnSurface
    /// representation by the `Geom_Surface` handle (plus location); the port
    /// has identity locations, so the registered surface's data pointer is the
    /// key. Faces that share a surface therefore share one representation, as
    /// in OCCT.
    fn repr_key(&self, face_key: usize) -> usize {
        self.face_surfaces
            .read()
            .unwrap()
            .get(&face_key)
            .map(|s| Arc::as_ptr(s) as *const () as usize)
            .unwrap_or(face_key)
    }

    pub fn face_geom(&self, s: &TopoShape) -> Option<FaceGeom> {
        // T-25 batch 4: read off the face's own `TShape`; the side table is the
        // safety net while other slots are still being lifted.
        // T-25 end phase: the face geometry comes off its own `TShape`.
        let ts = s.tshape.read().unwrap();
        let c = ts.face_core()?;
        let surf = c.surface.clone()?;
        Some(FaceGeom {
            surface: surf,
            tolerance: c.tolerance,
            natural_restriction: c.natural_restriction,
        })
    }

    /// The underlying surface handle (clone of the Arc).
    pub fn face_surface(&self, s: &TopoShape) -> Option<Arc<dyn Surface>> {
        self.face_geom(s).map(|g| g.surface)
    }

    pub fn face_tolerance(&self, s: &TopoShape) -> f64 {
        self.face_geom(s).map(|g| g.tolerance).unwrap_or(0.0)
    }

    pub fn natural_restriction(&self, s: &TopoShape) -> bool {
        self.face_geom(s).map(|g| g.natural_restriction).unwrap_or(true)
    }

    /// `BRep_Builder::NaturalRestriction(face, flag)`.
    pub fn set_natural_restriction(&self, s: &TopoShape, flag: bool) {
        s.tshape.write().unwrap().face_core_mut().natural_restriction = flag;
    }

    /// `BRep_Builder::UpdateFace` tolerance write.
    pub fn set_face_tolerance(&self, s: &TopoShape, tol: f64) {
        s.tshape.write().unwrap().face_core_mut().tolerance = tol;
    }

    // ---- lifecycle ----

    /// Drop all geometry entries owned by `s`. Call when a shape is discarded
    /// to keep the side-table from growing without bound.
    pub fn clear_shape(&self, s: &TopoShape) {
        let k = key(s);

        self.face_surfaces.write().unwrap().remove(&k);
        // T-25: the edge geometry lives on the shape itself now.
        {
            let mut ts = s.tshape.write().unwrap();
            ts.edge_pcurves = None;
            ts.edge_core = None;
            ts.vertex_core = None;
            ts.face_core = None;
        }
    }

    /// Remove every geometry entry keyed by the raw `TShape` address. Called
    /// from `TShape::drop` so entries die with their shape — this prevents a
    /// stale entry from leaking into a future shape that reuses the address.
    pub fn remove_by_ptr(&self, ptr: usize, id: u64) {
        {
            let mut ids = self.ids.write().unwrap();
            if ids.get(&ptr).copied() != Some(id) {
                // The address has been reused by a later shape: its geometry is
                // not ours to erase.
                return;
            }
            ids.remove(&ptr);
        }

        self.face_surfaces.write().unwrap().remove(&ptr);
    }

    /// Number of live entries (vertices + edges + faces).
    pub fn len(&self) -> usize {
        self.ids.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    /// Remove every entry (for tests / teardown).
    pub fn clear_all(&self) {
        self.ids.write().unwrap().clear();

        self.face_surfaces.write().unwrap().clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;
    use occt_geom::GeomLine;

    #[test]
    fn vertex_roundtrip() {
        let reg = GeometryRegistry::global();
        let v = TopoShape::new(ShapeType::Vertex);
        reg.set_vertex(&v, VertexGeom { point: GpPnt::new(1.0, 2.0, 3.0), tolerance: 1e-3 });
        let g = reg.vertex_geom(&v).expect("registered");
        assert!(g.point.is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        assert_eq!(g.tolerance, 1e-3);
        assert!(reg.vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
        reg.clear_shape(&v);
        assert!(reg.vertex_geom(&v).is_none());
    }

    #[test]
    fn edge_roundtrip_with_curve() {
        use occt_core::gp::{GpAx1, GpDir, GpLin};
        let reg = GeometryRegistry::global();
        let e = TopoShape::new(ShapeType::Edge);
        let lin = GpLin::new(GpAx1::new(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap()));
        let curve = Arc::new(GeomLine::new(lin));
        reg.set_edge(&e, EdgeGeom::new(curve, 0.0, 5.0));
        assert_eq!(reg.edge_parameters(&e), (0.0, 5.0));
        let c = reg.edge_curve(&e).expect("curve present");
        let p = c.d0(2.5);
        assert!((p.x() - 2.5).abs() < 1e-12);
        reg.clear_shape(&e);
        assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
    }

    #[test]
    fn unregistered_shapes_fall_back() {
        let reg = GeometryRegistry::global();
        let e = TopoShape::new(ShapeType::Edge);
        // Default: unbounded range, no curve, origin point for vertex.
        assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
        assert!(reg.edge_curve(&e).is_none());
        let v = TopoShape::new(ShapeType::Vertex);
        assert!(reg.vertex_point(&v).is_equal(&GpPnt::zero()));
    }
}
