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
use crate::tshape::EdgePcurves;

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
    /// Per-face 2D pcurves, keyed by face pointer identity (see `shape_key`).
    /// `BRep_Tool::CurveOnSurface` finds a representation by `Geom_Surface`
    /// handle, so faces that share a surface share these pcurves: lookup
    /// falls back to another face key on the same `Arc<dyn Surface>`.
    /// A seam edge of a periodic surface carries *two* pcurves on the same
    /// face (one per side of the seam), stored in forward-then-reversed order.
    pub pcurves: HashMap<usize, Vec<Arc<dyn Curve2d>>>,
    /// Per-face `BRep_GCurve` First/Last of the CurveOnSurface representation.
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
    vertices: RwLock<HashMap<usize, VertexGeom>>,
    edges: RwLock<HashMap<usize, EdgeGeom>>,
    faces: RwLock<HashMap<usize, FaceGeom>>,
    /// `TShape` address -> the `id` of the shape that registered it. A drop
    /// removes the geometry only when the ids match, so a stale drop cannot
    /// erase a later shape's entries at a reused address.
    ids: RwLock<HashMap<usize, u64>>,
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
            vertices: RwLock::new(HashMap::new()),
            edges: RwLock::new(HashMap::new()),
            faces: RwLock::new(HashMap::new()),
            ids: RwLock::new(HashMap::new()),
        })
    }

    // ---- vertices ----

    pub fn set_vertex(&self, s: &TopoShape, geom: VertexGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        self.vertices.write().unwrap().insert(k, geom);
    }

    pub fn vertex_geom(&self, s: &TopoShape) -> Option<VertexGeom> {
        self.vertices.read().unwrap().get(&key(s)).copied()
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
        self.edges.write().unwrap().insert(k, geom);
    }

    /// `BRep_Builder::Range(E, First, Last)`: set the edge's parameter range.
    /// `BRepPrim_OneAxis` relies on it for the degenerate pole edges
    /// (`SetParameters(ETOP/EBOTTOM, …, 0., myAngle)`,
    /// `BRepPrim_OneAxis.cxx:407`, `:418`): the 3D curve is absent but the range
    /// still spans the full period so the edge's pcurve is a full u-isoline.
    pub fn set_edge_range(&self, s: &TopoShape, first: f64, last: f64) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.first = first;
            g.last = last;
        }
    }

    pub fn edge_geom(&self, s: &TopoShape) -> Option<EdgeGeom> {
        self.edges.read().unwrap().get(&key(s)).map(|g| {
            // EdgeGeom is not Clone (Arc<dyn Curve> is Clone, but we rebuild a
            // fresh struct to avoid needing Clone on the whole thing).
            EdgeGeom {
                curve: g.curve.clone(),
                first: g.first,
                last: g.last,
                tolerance: g.tolerance,
                same_parameter: g.same_parameter,
                same_range: g.same_range,
                degenerated: g.degenerated,
                // T-25: pcurves come off the edge's own `TShape`.
                pcurves: s.tshape.read().unwrap().edge_pcurves().map(|p| p.curves.clone()).unwrap_or_default(),
                pcurve_ranges: s.tshape.read().unwrap().edge_pcurves().map(|p| p.ranges.clone()).unwrap_or_default(),
            }
        })
    }

    /// The underlying curve handle (clone of the Arc).
    pub fn edge_curve(&self, s: &TopoShape) -> Option<Arc<dyn Curve>> {
        self.edge_geom(s).map(|g| g.curve)
    }

    pub fn edge_parameters(&self, s: &TopoShape) -> (f64, f64) {
        match self.edge_geom(s) {
            Some(g) => (g.first, g.last),
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
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.degenerated = v;
        }
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
        // T-25: read off the edge's own `TShape`; the same-surface fallback
        // still consults this registry's face map.
        let candidates: Vec<(usize, Vec<Arc<dyn Curve2d>>)> = {
            let ts = s.tshape.read().unwrap();
            let Some(p) = ts.edge_pcurves() else {
                return Vec::new();
            };
            let direct = p.get_pcurves(face_key);
            if !direct.is_empty() {
                return direct;
            }
            p.curves
                .iter()
                .filter(|(&fk, _)| fk != face_key)
                .map(|(&fk, cs)| (fk, cs.clone()))
                .collect()
        };
        let want = self.faces.read().unwrap().get(&face_key).map(|g| g.surface.clone());
        let Some(want) = want else {
            return Vec::new();
        };
        let faces = self.faces.read().unwrap();
        for (fk, cs) in candidates {
            if faces.get(&fk).is_some_and(|fg| Arc::ptr_eq(&fg.surface, &want)) {
                return cs;
            }
        }
        Vec::new()
    }

    /// Attach a pcurve to edge `s` for the face identified by `face_key`.
    /// Mirrors `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`.
    pub fn set_edge_pcurve(&self, s: &TopoShape, face_key: usize, curve: Arc<dyn Curve2d>) {
        let (first, last) = self.edge_parameters(s);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurve(face_key, curve, first, last);
    }

    /// Replace the pcurves of edge `s` on the face identified by `face_key`.
    /// Mirrors the seam overload `BRep_Builder::UpdateEdge(edge, c1, c2, face)`.
    pub fn set_edge_pcurves(&self, s: &TopoShape, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        let (first, last) = self.edge_parameters(s);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurves(face_key, curves, first, last);
    }

    /// `BRep_Builder::Range(edge, face, first, last)` (`BRep_Builder.cxx:1121`).
    pub fn set_pcurve_range(&self, s: &TopoShape, face_key: usize, first: f64, last: f64) {
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .ranges
            .insert(face_key, (first, last));
        if false {
            let _ = (first, last);
        }
    }

    /// COS representation `[First, Last]` for `face_key`, with the same
    /// same-surface fallback as [`GeometryRegistry::edge_pcurve`].
    pub fn pcurve_range(&self, s: &TopoShape, face_key: usize) -> Option<(f64, f64)> {
        // T-25: read off the edge's own `TShape`; the same-surface fallback
        // still consults this registry's face map.
        let candidates: Vec<(usize, (f64, f64))> = {
            let ts = s.tshape.read().unwrap();
            let p = ts.edge_pcurves()?;
            if let Some(&r) = p.ranges.get(&face_key) {
                return Some(r);
            }
            p.ranges.iter().filter(|(&fk, _)| fk != face_key).map(|(&fk, &r)| (fk, r)).collect()
        };
        let want = self.faces.read().unwrap().get(&face_key).map(|fg| fg.surface.clone());
        let Some(want) = want else {
            return None;
        };
        let faces = self.faces.read().unwrap();
        for (fk, r) in candidates {
            if faces.get(&fk).is_some_and(|fg| Arc::ptr_eq(&fg.surface, &want)) {
                return Some(r);
            }
        }
        None
    }

    /// `BRep_Builder::SameRange`.
    pub fn set_same_range(&self, s: &TopoShape, value: bool) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.same_range = value;
        }
    }

    /// `BRep_Builder::SameParameter`.
    pub fn set_same_parameter(&self, s: &TopoShape, value: bool) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.same_parameter = value;
        }
    }

    /// `BRep_Builder::UpdateEdge` tolerance write.
    pub fn set_edge_tolerance(&self, s: &TopoShape, tol: f64) {
        if let Some(g) = self.edges.write().unwrap().get_mut(&key(s)) {
            g.tolerance = tol;
        }
    }

    /// CurveOnSurface representations: surface, pcurve, COS `[first, last]`.
    pub fn edge_pcurve_reps(
        &self,
        s: &TopoShape,
    ) -> Vec<(Arc<dyn Surface>, Arc<dyn Curve2d>, f64, f64)> {
        let faces = self.faces.read().unwrap();
        // T-25: pcurves and their COS ranges come off the edge's own `TShape`;
        // the 3D range still lives on `EdgeGeom` in this batch.
        let (slot, three_d) = {
            let edges = self.edges.read().unwrap();
            let Some(g) = edges.get(&key(s)) else {
                return Vec::new();
            };
            (s.tshape.read().unwrap().edge_pcurves().cloned(), (g.first, g.last))
        };
        let Some(g) = slot else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (&fk, pcs) in &g.curves {
            let Some(fg) = faces.get(&fk) else {
                continue;
            };
            let (a, b) = g.ranges.get(&fk).copied().unwrap_or(three_d);
            for pc in pcs {
                out.push((fg.surface.clone(), pc.clone(), a, b));
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
        let drop_keys: Vec<usize> = {
            let faces = self.faces.read().unwrap();
            let ts = edge.tshape.read().unwrap();
            let Some(g) = ts.edge_pcurves() else {
                return;
            };
            g.curves
                .keys()
                .copied()
                .filter(|&fk| {
                    faces
                        .get(&fk)
                        .map(|fg| Arc::ptr_eq(&fg.surface, &want))
                        .unwrap_or(fk == key(face))
                })
                .collect()
        };
        if drop_keys.is_empty() {
            return;
        }
        let mut ts = edge.tshape.write().unwrap();
        if let Some(g) = ts.edge_pcurves.as_mut() {
            for fk in drop_keys {
                g.curves.remove(&fk);
                g.ranges.remove(&fk);
            }
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
        self.faces.write().unwrap().insert(k, geom);
    }

    pub fn face_geom(&self, s: &TopoShape) -> Option<FaceGeom> {
        self.faces.read().unwrap().get(&key(s)).map(|g| FaceGeom {
            surface: g.surface.clone(),
            tolerance: g.tolerance,
            natural_restriction: g.natural_restriction,
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
        if let Some(g) = self.faces.write().unwrap().get_mut(&key(s)) {
            g.natural_restriction = flag;
        }
    }

    /// `BRep_Builder::UpdateFace` tolerance write.
    pub fn set_face_tolerance(&self, s: &TopoShape, tol: f64) {
        if let Some(g) = self.faces.write().unwrap().get_mut(&key(s)) {
            g.tolerance = tol;
        }
    }

    // ---- lifecycle ----

    /// Drop all geometry entries owned by `s`. Call when a shape is discarded
    /// to keep the side-table from growing without bound.
    pub fn clear_shape(&self, s: &TopoShape) {
        let k = key(s);
        self.vertices.write().unwrap().remove(&k);
        self.edges.write().unwrap().remove(&k);
        self.faces.write().unwrap().remove(&k);
        // T-25: the pcurves live on the shape itself now.
        s.tshape.write().unwrap().edge_pcurves = None;
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
        self.vertices.write().unwrap().remove(&ptr);
        self.edges.write().unwrap().remove(&ptr);
        self.faces.write().unwrap().remove(&ptr);
    }

    /// Number of live entries (vertices + edges + faces).
    pub fn len(&self) -> usize {
        self.vertices.read().unwrap().len()
            + self.edges.read().unwrap().len()
            + self.faces.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    /// Remove every entry (for tests / teardown).
    pub fn clear_all(&self) {
        self.ids.write().unwrap().clear();
        self.vertices.write().unwrap().clear();
        self.edges.write().unwrap().clear();
        self.faces.write().unwrap().clear();
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
