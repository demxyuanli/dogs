//! `ChFiDS_CommonPoint` — a start/end point of a fillet, shared by two
//! adjacent fillets and by one edge on one of the two faces taking part in the
//! construction.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_CommonPoint.hxx:33-151`,
//! `ChFiDS_CommonPoint.cxx:24-111` (OCCT 8.0.0).

use occt_core::gp::{GpPnt, GpVec};

use crate::abs::Orientation;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;

/// Source: `ChFiDS_CommonPoint.hxx:35-151`.
///
/// Field order and every initial value come from the constructor at
/// `cxx:24-32`: `tol`, `prmarc` and `prmtg` are `0.0`, `traarc` is `FORWARD`,
/// `isonarc` / `isvtx` / `hasvector` are `false`.
#[derive(Debug, Clone)]
pub struct CommonPoint {
    /// `arc`
    arc: Option<TopoShape>,
    /// `vtx`
    vtx: Option<TopoShape>,
    /// `point`
    point: GpPnt,
    /// `vector`
    vector: GpVec,
    /// `tol`
    tol: f64,
    /// `prmarc`
    prmarc: f64,
    /// `prmtg`
    prmtg: f64,
    /// `traarc`
    traarc: Orientation,
    /// `isonarc`
    isonarc: bool,
    /// `isvtx`
    isvtx: bool,
    /// `hasvector`
    hasvector: bool,
}

impl Default for CommonPoint {
    /// `ChFiDS_CommonPoint()` (`cxx:24-32`).
    fn default() -> Self {
        Self {
            arc: None,
            vtx: None,
            point: GpPnt::default(),
            vector: GpVec::default(),
            tol: 0.0,
            prmarc: 0.0,
            prmtg: 0.0,
            traarc: Orientation::Forward,
            isonarc: false,
            isvtx: false,
            hasvector: false,
        }
    }
}

impl CommonPoint {
    /// `ChFiDS_CommonPoint()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Reset()` (`cxx:36-44`): clears `tol`, `isvtx`, `isonarc` and
    /// `hasvector` only. OCCT leaves `arc` / `vtx` / `point` / `vector` /
    /// `prmarc` / `prmtg` / `traarc` untouched, and the port does the same.
    pub fn reset(&mut self) {
        self.tol = 0.0;
        self.isvtx = false;
        self.isonarc = false;
        self.hasvector = false;
    }

    /// `SetVertex(V)` (`hxx:49-54`).
    pub fn set_vertex(&mut self, vertex: &TopoShape) {
        self.isvtx = true;
        self.vtx = Some(vertex.clone());
    }

    /// `SetArc(Tol, A, Param, TArc)` (`cxx:51-66`).
    ///
    /// The `if (Tol > tol)` guard at `cxx:59-62` is the "do not squash the
    /// tolerance" fix dated 30/09/1997 in the comment above the function: the
    /// stored tolerance is a running maximum, never a replacement.
    pub fn set_arc(&mut self, tol: f64, arc: &TopoShape, param: f64, t_arc: Orientation) {
        self.isonarc = true;
        if tol > self.tol {
            self.tol = tol;
        }
        self.arc = Some(arc.clone());
        self.prmarc = param;
        self.traarc = t_arc;
    }

    /// `SetParameter(Param)` (`cxx:68-71`).
    pub fn set_parameter(&mut self, param: f64) {
        self.prmtg = param;
    }

    /// `SetPoint(thePoint)` (`hxx:57`).
    pub fn set_point(&mut self, point: &GpPnt) {
        self.point = *point;
    }

    /// `SetVector(theVector)` (`hxx:69-73`).
    pub fn set_vector(&mut self, vector: &GpVec) {
        self.hasvector = true;
        self.vector = *vector;
    }

    /// `SetTolerance(Tol)` (`hxx:76-82`): a running maximum, like
    /// [`Self::set_arc`].
    pub fn set_tolerance(&mut self, tol: f64) {
        if tol > self.tol {
            self.tol = tol;
        }
    }

    /// `Tolerance()` (`hxx:85`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `IsVertex()` (`hxx:88`).
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// `Vertex()` (`hxx:96-101`): OCCT raises `Standard_DomainError` when the
    /// point is not a vertex; the port answers `None`.
    pub fn vertex(&self) -> Option<&TopoShape> {
        if !self.isvtx {
            return None;
        }
        self.vtx.as_ref()
    }

    /// `IsOnArc()` (`hxx:106`).
    pub fn is_on_arc(&self) -> bool {
        self.isonarc
    }

    /// `Arc()` (`cxx:75-84`): OCCT raises when the point is not on an arc.
    pub fn arc(&self) -> Option<&TopoShape> {
        if !self.isonarc {
            return None;
        }
        self.arc.as_ref()
    }

    /// `TransitionOnArc()` (`cxx:86-95`).
    pub fn transition_on_arc(&self) -> Option<Orientation> {
        if !self.isonarc {
            return None;
        }
        Some(self.traarc)
    }

    /// `ParameterOnArc()` (`cxx:97-106`).
    pub fn parameter_on_arc(&self) -> Option<f64> {
        if !self.isonarc {
            return None;
        }
        Some(self.prmarc)
    }

    /// `Parameter()` (`cxx:108-111`): the parameter on the spine. Unlike the
    /// three accessors above this one has no guard in OCCT.
    pub fn parameter(&self) -> f64 {
        self.prmtg
    }

    /// `Point()` (`hxx:123`).
    pub fn point(&self) -> &GpPnt {
        &self.point
    }

    /// `HasVector()` (`hxx:127`).
    pub fn has_vector(&self) -> bool {
        self.hasvector
    }

    /// `Vector()` (`hxx:130-136`): OCCT raises `Standard_DomainError` when no
    /// vector was stored.
    pub fn vector(&self) -> Option<&GpVec> {
        if !self.hasvector {
            return None;
        }
        Some(&self.vector)
    }

    /// `TopoDS_Vertex::IsSame` against [`Self::vertex`] — the comparison
    /// `ChFi3d` uses to match a common point to a spine vertex.
    pub fn is_vertex_same(&self, other: &TopoShape) -> bool {
        match self.vtx.as_ref() {
            Some(v) => is_same(v, other),
            None => false,
        }
    }
}
