//! Discrete mesh data model — port of `BRepMeshData::{Model,Edge,Face,Wire,Curve,PCurve}`
//! and `IMeshData::{Status,StatusOwner,TessellatedShape}`.
//!
//! Mirrors OCCT's `IMeshData`/`BRepMeshData` layer. A [`MeshModel`] holds
//! indexed collections of discrete edges, faces and wires:
//!
//! * [`MeshEdge`] carries the topological [`Edge`], the 3D curve geometry
//!   (`Arc<dyn Curve>`), its parameter range, a [`MeshCurve`] discretization
//!   and per-face [`MeshPCurve`]s.
//! * [`MeshFace`] carries the topological [`Face`], the surface geometry
//!   (`Arc<dyn Surface>`), its wires (referenced by model index) and internal
//!   face points.
//! * [`MeshWire`] is an ordered chain of edge indices + orientations.
//! * Every entity is a [`MeshStatus`] status owner; [`MeshModel::status_mask`]
//!   aggregates the statuses of all contained entities.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::{Curve, Surface};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;

/// OCCT `RealLast()` — largest finite double, used as the "unset" deflection
/// sentinel (matches `IMeshData_TessellatedShape` default).
const REAL_LAST: f64 = f64::MAX;

/// Out-of-range error message, mirroring OCCT's `Standard_OutOfRange`.
fn out_of_range(op: &str, index: usize, len: usize) -> String {
    format!("{op}: index {index} out of range (len {len})")
}

/// Bit flag status of a discrete model entity. Source: `IMeshData_Status.hxx`
/// + `IMeshData_StatusOwner.hxx`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct MeshStatus(u32);

impl MeshStatus {
    /// Mesh generation is successful.
    pub const NO_ERROR: MeshStatus = MeshStatus(0x0);
    /// Open wire problem — can potentially lead to incorrect results.
    pub const OPEN_WIRE: MeshStatus = MeshStatus(0x1);
    /// Self-intersections on a discretized wire.
    pub const SELF_INTERSECTING_WIRE: MeshStatus = MeshStatus(0x2);
    /// Failed to generate mesh for some faces.
    pub const FAILURE: MeshStatus = MeshStatus(0x4);
    /// Deflection of some edges was decreased due to interference.
    pub const REMESH: MeshStatus = MeshStatus(0x8);
    /// Bad orientation of a wire.
    pub const UNORIENTED_WIRE: MeshStatus = MeshStatus(0x10);
    /// Discrete model contains too few boundary points.
    pub const TOO_FEW_POINTS: MeshStatus = MeshStatus(0x20);
    /// Existing triangulation corresponds to a greater deflection than requested.
    pub const OUTDATED: MeshStatus = MeshStatus(0x40);
    /// Existing triangulation was reused as it fits the requested deflection.
    pub const REUSED: MeshStatus = MeshStatus(0x80);
    /// User break.
    pub const USER_BREAK: MeshStatus = MeshStatus(0x100);

    /// Raw bit mask.
    pub const fn bits(self) -> u32 {
        self.0
    }

    /// Build from a raw bit mask.
    pub const fn from_bits(bits: u32) -> Self {
        Self(bits)
    }

    /// True when no flag is set.
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// True when all flags of `other` are set in `self`.
    pub const fn contains(self, other: MeshStatus) -> bool {
        self.0 & other.0 == other.0
    }

    /// Union of two masks.
    pub const fn union(self, other: MeshStatus) -> MeshStatus {
        Self(self.0 | other.0)
    }

    /// Adds flags.
    pub const fn insert(&mut self, other: MeshStatus) {
        self.0 |= other.0;
    }

    /// Clears flags.
    pub const fn remove(&mut self, other: MeshStatus) {
        self.0 &= !other.0;
    }

    // ---- OCCT `IMeshData_StatusOwner` API ----

    /// True if status is strictly equal to `other`.
    pub const fn is_equal(self, other: MeshStatus) -> bool {
        self.0 == other.0
    }

    /// True if `other` flag is set.
    pub const fn is_set(self, other: MeshStatus) -> bool {
        self.0 & other.0 != 0
    }

    /// Adds `other` to the status flags.
    pub fn set_status(&mut self, other: MeshStatus) {
        self.0 |= other.0;
    }

    /// Removes `other` from the status flags.
    pub fn unset_status(&mut self, other: MeshStatus) {
        self.0 &= !other.0;
    }

    /// Complete status mask.
    pub const fn status_mask(self) -> u32 {
        self.0
    }
}

/// Discrete 3D curve of an edge — point/parameter discretization.
/// Source: `BRepMeshData_Curve.hxx` + `IMeshData_Curve.hxx`.
#[derive(Debug, Clone, Default)]
pub struct MeshCurve {
    points: Vec<GpPnt>,
    parameters: Vec<f64>,
    deflection: f64,
}

impl MeshCurve {
    /// Empty discretization; deflection defaults to `RealLast()` (unset).
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts a discretization point at the given position.
    pub fn insert_point(&mut self, position: usize, point: GpPnt, param: f64) -> Result<(), String> {
        if position > self.points.len() {
            return Err(out_of_range("MeshCurve::insert_point", position, self.points.len()));
        }
        self.points.insert(position, point);
        self.parameters.insert(position, param);
        Ok(())
    }

    /// Appends a discretization point.
    pub fn add_point(&mut self, point: GpPnt, param: f64) {
        self.points.push(point);
        self.parameters.push(param);
    }

    /// Discretization point at the given index.
    pub fn get_point(&self, index: usize) -> Result<GpPnt, String> {
        self.points
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshCurve::get_point", index, self.points.len()))
    }

    /// Mutable access to the discretization point at the given index.
    pub fn get_point_mut(&mut self, index: usize) -> Result<&mut GpPnt, String> {
        if index >= self.points.len() {
            return Err(out_of_range("MeshCurve::get_point_mut", index, self.points.len()));
        }
        Ok(&mut self.points[index])
    }

    /// Removes the point (and its parameter) at the given index.
    pub fn remove_point(&mut self, index: usize) -> Result<(), String> {
        if index >= self.points.len() {
            return Err(out_of_range("MeshCurve::remove_point", index, self.points.len()));
        }
        self.points.remove(index);
        self.parameters.remove(index);
        Ok(())
    }

    /// Curve parameter at the given index.
    pub fn get_parameter(&self, index: usize) -> Result<f64, String> {
        self.parameters
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshCurve::get_parameter", index, self.parameters.len()))
    }

    /// Mutable access to the parameter at the given index.
    pub fn get_parameter_mut(&mut self, index: usize) -> Result<&mut f64, String> {
        if index >= self.parameters.len() {
            return Err(out_of_range("MeshCurve::get_parameter_mut", index, self.parameters.len()));
        }
        Ok(&mut self.parameters[index])
    }

    /// Number of stored parameters (== number of points).
    pub fn parameters_nb(&self) -> usize {
        self.parameters.len()
    }

    /// Clears the discretization. When `keep_end_points` is set, only the
    /// first and last entries are retained.
    pub fn clear(&mut self, keep_end_points: bool) {
        if !keep_end_points {
            self.points.clear();
            self.parameters.clear();
        } else if self.parameters.len() > 2 {
            let n = self.parameters.len();
            self.points.drain(1..n - 1);
            self.parameters.drain(1..n - 1);
        }
    }

    /// All discretization points.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }

    /// All parameters.
    pub fn parameters(&self) -> &[f64] {
        &self.parameters
    }

    /// Deflection value for the discrete curve.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Sets the deflection value for the discrete curve.
    pub fn set_deflection(&mut self, value: f64) {
        self.deflection = value;
    }
}

/// Discrete pcurve of an edge associated with a discrete face.
/// Source: `BRepMeshData_PCurve.hxx` + `IMeshData_PCurve.hxx`.
#[derive(Debug, Clone)]
pub struct MeshPCurve {
    face: usize,
    orientation: Orientation,
    points: Vec<GpPnt2d>,
    parameters: Vec<f64>,
    indices: Vec<i32>,
    deflection: f64,
}

impl MeshPCurve {
    /// Creates a pcurve for the face with the given model index and edge
    /// orientation. Deflection defaults to `RealLast()` (unset).
    pub fn new(face: usize, orientation: Orientation) -> Self {
        Self {
            face,
            orientation,
            points: Vec::new(),
            parameters: Vec::new(),
            indices: Vec::new(),
            deflection: REAL_LAST,
        }
    }

    /// Inserts a discretization point at the given position.
    pub fn insert_point(
        &mut self,
        position: usize,
        point: GpPnt2d,
        param: f64,
    ) -> Result<(), String> {
        if position > self.points.len() {
            return Err(out_of_range("MeshPCurve::insert_point", position, self.points.len()));
        }
        self.points.insert(position, point);
        self.parameters.insert(position, param);
        self.indices.insert(position, 0);
        Ok(())
    }

    /// Appends a discretization point (mesh index initialized to 0).
    pub fn add_point(&mut self, point: GpPnt2d, param: f64) {
        self.points.push(point);
        self.parameters.push(param);
        self.indices.push(0);
    }

    /// Discretization point at the given index.
    pub fn get_point(&self, index: usize) -> Result<GpPnt2d, String> {
        self.points
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshPCurve::get_point", index, self.points.len()))
    }

    /// Mutable access to the discretization point at the given index.
    pub fn get_point_mut(&mut self, index: usize) -> Result<&mut GpPnt2d, String> {
        if index >= self.points.len() {
            return Err(out_of_range("MeshPCurve::get_point_mut", index, self.points.len()));
        }
        Ok(&mut self.points[index])
    }

    /// Mesh node index corresponding to the discretization point at the index.
    pub fn get_index(&self, index: usize) -> Result<i32, String> {
        self.indices
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshPCurve::get_index", index, self.indices.len()))
    }

    /// Mutable access to the mesh node index at the given index.
    pub fn get_index_mut(&mut self, index: usize) -> Result<&mut i32, String> {
        if index >= self.indices.len() {
            return Err(out_of_range("MeshPCurve::get_index_mut", index, self.indices.len()));
        }
        Ok(&mut self.indices[index])
    }

    /// Removes the point (and its parameter/index) at the given index.
    pub fn remove_point(&mut self, index: usize) -> Result<(), String> {
        if index >= self.points.len() {
            return Err(out_of_range("MeshPCurve::remove_point", index, self.points.len()));
        }
        self.points.remove(index);
        self.parameters.remove(index);
        self.indices.remove(index);
        Ok(())
    }

    /// Curve parameter at the given index.
    pub fn get_parameter(&self, index: usize) -> Result<f64, String> {
        self.parameters
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshPCurve::get_parameter", index, self.parameters.len()))
    }

    /// Mutable access to the parameter at the given index.
    pub fn get_parameter_mut(&mut self, index: usize) -> Result<&mut f64, String> {
        if index >= self.parameters.len() {
            return Err(out_of_range("MeshPCurve::get_parameter_mut", index, self.parameters.len()));
        }
        Ok(&mut self.parameters[index])
    }

    /// Number of stored parameters (== number of points).
    pub fn parameters_nb(&self) -> usize {
        self.parameters.len()
    }

    /// Clears the discretization. When `keep_end_points` is set, only the
    /// first and last entries are retained.
    pub fn clear(&mut self, keep_end_points: bool) {
        if !keep_end_points {
            self.points.clear();
            self.parameters.clear();
            self.indices.clear();
        } else if self.parameters.len() > 2 {
            let n = self.parameters.len();
            self.points.drain(1..n - 1);
            self.parameters.drain(1..n - 1);
            self.indices.drain(1..n - 1);
        }
    }

    /// True if the associated edge orientation is forward.
    pub fn is_forward(&self) -> bool {
        self.orientation == Orientation::Forward
    }

    /// True if the associated edge orientation is internal.
    pub fn is_internal(&self) -> bool {
        self.orientation == Orientation::Internal
    }

    /// Orientation of the edge associated with this pcurve.
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Model index of the discrete face this pcurve is associated to.
    pub fn face(&self) -> usize {
        self.face
    }

    /// All 2D discretization points.
    pub fn points(&self) -> &[GpPnt2d] {
        &self.points
    }

    /// All parameters.
    pub fn parameters(&self) -> &[f64] {
        &self.parameters
    }

    /// All mesh node indices.
    pub fn indices(&self) -> &[i32] {
        &self.indices
    }

    /// Deflection value for the discrete pcurve.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Sets the deflection value for the discrete pcurve.
    pub fn set_deflection(&mut self, value: f64) {
        self.deflection = value;
    }
}

/// Discrete model of an edge. Source: `BRepMeshData_Edge.hxx` + `IMeshData_Edge.hxx`.
#[derive(Clone)]
pub struct MeshEdge {
    edge: Edge,
    curve: Option<Arc<dyn Curve>>,
    first: f64,
    last: f64,
    discretization: MeshCurve,
    pcurves: Vec<MeshPCurve>,
    same_param: bool,
    same_range: bool,
    degenerated: bool,
    angular_deflection: f64,
    deflection: f64,
    status: MeshStatus,
}

impl MeshEdge {
    /// Builds a discrete edge from the topological edge, reading the 3D curve,
    /// parameter range and `BRep_TEdge` flags from the geometry registry.
    pub fn new(edge: Edge) -> Self {
        let geom = GeometryRegistry::global().edge_geom(&edge.0);
        let curve = BRepTool::edge_curve(&edge);
        let (first, last) = BRepTool::edge_parameters(&edge);
        Self {
            edge,
            curve,
            first,
            last,
            discretization: MeshCurve::new(),
            pcurves: Vec::new(),
            same_param: geom.as_ref().map(|g| g.same_parameter).unwrap_or(false),
            same_range: geom.as_ref().map(|g| g.same_range).unwrap_or(false),
            degenerated: geom.as_ref().map(|g| g.degenerated).unwrap_or(false),
            angular_deflection: REAL_LAST,
            deflection: REAL_LAST,
            status: MeshStatus::NO_ERROR,
        }
    }

    /// Topological edge attached to this discrete edge.
    pub fn edge(&self) -> &Edge {
        &self.edge
    }

    /// 3D curve geometry of the edge.
    pub fn curve(&self) -> Option<Arc<dyn Curve>> {
        self.curve.clone()
    }

    /// Replaces the 3D curve geometry.
    pub fn set_curve(&mut self, curve: Option<Arc<dyn Curve>>) {
        self.curve = curve;
    }

    /// First parameter of the edge range.
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// Last parameter of the edge range.
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// Updates the edge parameter range.
    pub fn set_parameters(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;
    }

    /// The 3D curve discretization (points + parameters).
    pub fn discretization(&self) -> &MeshCurve {
        &self.discretization
    }

    /// Mutable access to the 3D curve discretization.
    pub fn discretization_mut(&mut self) -> &mut MeshCurve {
        &mut self.discretization
    }

    /// Number of pcurves assigned to this edge.
    pub fn pcurves_nb(&self) -> usize {
        self.pcurves.len()
    }

    /// Pcurve with the given index.
    pub fn pcurve(&self, index: usize) -> Result<&MeshPCurve, String> {
        self.pcurves
            .get(index)
            .ok_or_else(|| out_of_range("MeshEdge::pcurve", index, self.pcurves.len()))
    }

    /// Adds a pcurve for the given discrete face and orientation; returns the
    /// index of the added pcurve.
    pub fn add_pcurve(&mut self, face: usize, orientation: Orientation) -> usize {
        self.pcurves.push(MeshPCurve::new(face, orientation));
        self.pcurves.len() - 1
    }

    /// Pcurve for the given discrete face with the matching orientation.
    /// Mirrors `IMeshData_Edge::GetPCurve(face, orientation)`: when several
    /// pcurves exist for the face, the one with the requested orientation is
    /// returned, otherwise the last one.
    pub fn pcurve_for(&self, face: usize, orientation: Orientation) -> Option<&MeshPCurve> {
        let candidates: Vec<&MeshPCurve> =
            self.pcurves.iter().filter(|p| p.face() == face).collect();
        candidates
            .iter()
            .find(|p| p.orientation() == orientation)
            .or_else(|| candidates.last())
            .copied()
    }

    /// Indices of all pcurves assigned to the given discrete face.
    pub fn pcurves_for(&self, face: usize) -> Vec<usize> {
        self.pcurves
            .iter()
            .enumerate()
            .filter(|(_, p)| p.face() == face)
            .map(|(i, _)| i)
            .collect()
    }

    /// True when the edge is free, i.e. it has no pcurves.
    pub fn is_free(&self) -> bool {
        self.pcurves.is_empty()
    }

    /// Clears the curve discretization and all pcurves.
    pub fn clear(&mut self, keep_end_points: bool) {
        self.discretization.clear(keep_end_points);
        for pc in &mut self.pcurves {
            pc.clear(keep_end_points);
        }
    }

    /// Same-parameter flag (from `BRep_TEdge`).
    pub fn same_param(&self) -> bool {
        self.same_param
    }

    /// Sets the same-parameter flag.
    pub fn set_same_param(&mut self, value: bool) {
        self.same_param = value;
    }

    /// Same-range flag (from `BRep_TEdge`).
    pub fn same_range(&self) -> bool {
        self.same_range
    }

    /// Sets the same-range flag.
    pub fn set_same_range(&mut self, value: bool) {
        self.same_range = value;
    }

    /// Degenerated flag (from `BRep_TEdge`).
    pub fn degenerated(&self) -> bool {
        self.degenerated
    }

    /// Sets the degenerated flag.
    pub fn set_degenerated(&mut self, value: bool) {
        self.degenerated = value;
    }

    /// Angular deflection of the edge.
    pub fn angular_deflection(&self) -> f64 {
        self.angular_deflection
    }

    /// Sets the angular deflection of the edge.
    pub fn set_angular_deflection(&mut self, value: f64) {
        self.angular_deflection = value;
    }

    /// Deflection of the edge (`IMeshData_TessellatedShape`).
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Sets the deflection of the edge.
    pub fn set_deflection(&mut self, value: f64) {
        self.deflection = value;
    }

    /// Current status flags.
    pub fn status(&self) -> MeshStatus {
        self.status
    }

    /// Adds status flags.
    pub fn set_status(&mut self, value: MeshStatus) {
        self.status.set_status(value);
    }

    /// Removes status flags.
    pub fn unset_status(&mut self, value: MeshStatus) {
        self.status.unset_status(value);
    }

    /// True if the given flag is set.
    pub fn is_status(&self, value: MeshStatus) -> bool {
        self.status.is_set(value)
    }
}

/// Discrete model of a wire — an ordered chain of edges.
/// Source: `BRepMeshData_Wire.hxx` + `IMeshData_Wire.hxx`.
#[derive(Debug, Clone)]
pub struct MeshWire {
    wire: Wire,
    edges: Vec<usize>,
    orientations: Vec<Orientation>,
    deflection: f64,
    status: MeshStatus,
}

impl MeshWire {
    /// Builds an empty discrete wire from the topological wire.
    pub fn new(wire: Wire) -> Self {
        Self {
            wire,
            edges: Vec::new(),
            orientations: Vec::new(),
            deflection: REAL_LAST,
            status: MeshStatus::NO_ERROR,
        }
    }

    /// Topological wire attached to this discrete wire.
    pub fn wire(&self) -> &Wire {
        &self.wire
    }

    /// Number of edges in the wire chain.
    pub fn edges_nb(&self) -> usize {
        self.edges.len()
    }

    /// Appends a discrete edge (model index) with the given orientation;
    /// returns the index of the added edge in the chain.
    pub fn add_edge(&mut self, edge_index: usize, orientation: Orientation) -> usize {
        self.edges.push(edge_index);
        self.orientations.push(orientation);
        self.edges.len() - 1
    }

    /// Model index of the edge at the given chain position.
    pub fn edge(&self, index: usize) -> Result<usize, String> {
        self.edges
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshWire::edge", index, self.edges.len()))
    }

    /// Orientation of the edge at the given chain position.
    pub fn edge_orientation(&self, index: usize) -> Result<Orientation, String> {
        self.orientations
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshWire::edge_orientation", index, self.orientations.len()))
    }

    /// All edge indices of the chain.
    pub fn edges(&self) -> &[usize] {
        &self.edges
    }

    /// All orientations of the chain.
    pub fn orientations(&self) -> &[Orientation] {
        &self.orientations
    }

    /// Deflection of the wire.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Sets the deflection of the wire.
    pub fn set_deflection(&mut self, value: f64) {
        self.deflection = value;
    }

    /// Current status flags.
    pub fn status(&self) -> MeshStatus {
        self.status
    }

    /// Adds status flags.
    pub fn set_status(&mut self, value: MeshStatus) {
        self.status.set_status(value);
    }

    /// Removes status flags.
    pub fn unset_status(&mut self, value: MeshStatus) {
        self.status.unset_status(value);
    }

    /// True if the given flag is set.
    pub fn is_status(&self, value: MeshStatus) -> bool {
        self.status.is_set(value)
    }
}

/// Discrete model of a face. Source: `BRepMeshData_Face.hxx` + `IMeshData_Face.hxx`.
#[derive(Clone)]
pub struct MeshFace {
    face: Face,
    surface: Option<Arc<dyn Surface>>,
    wires: Vec<usize>,
    points: Vec<GpPnt>,
    deflection: f64,
    status: MeshStatus,
}

impl MeshFace {
    /// Builds a discrete face from the topological face, reading the surface
    /// geometry from the registry.
    pub fn new(face: Face) -> Self {
        let surface = BRepTool::face_surface(&face);
        Self {
            face,
            surface,
            wires: Vec::new(),
            points: Vec::new(),
            deflection: REAL_LAST,
            status: MeshStatus::NO_ERROR,
        }
    }

    /// Topological face attached to this discrete face.
    pub fn face(&self) -> &Face {
        &self.face
    }

    /// Surface geometry of the face.
    pub fn surface(&self) -> Option<Arc<dyn Surface>> {
        self.surface.clone()
    }

    /// Replaces the surface geometry.
    pub fn set_surface(&mut self, surface: Option<Arc<dyn Surface>>) {
        self.surface = surface;
    }

    /// Number of wires of the face. The first wire is always the outer one.
    pub fn wires_nb(&self) -> usize {
        self.wires.len()
    }

    /// Appends a wire (model index); returns the index of the added wire.
    pub fn add_wire(&mut self, wire_index: usize) -> usize {
        self.wires.push(wire_index);
        self.wires.len() - 1
    }

    /// Model index of the wire at the given position.
    pub fn wire(&self, index: usize) -> Result<usize, String> {
        self.wires
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshFace::wire", index, self.wires.len()))
    }

    /// All wire indices of the face.
    pub fn wires(&self) -> &[usize] {
        &self.wires
    }

    /// Number of internal (in-face) points.
    pub fn points_nb(&self) -> usize {
        self.points.len()
    }

    /// Appends an internal face point.
    pub fn add_point(&mut self, point: GpPnt) {
        self.points.push(point);
    }

    /// Internal face point at the given index.
    pub fn point(&self, index: usize) -> Result<GpPnt, String> {
        self.points
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("MeshFace::point", index, self.points.len()))
    }

    /// All internal face points.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }

    /// Deflection of the face.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Sets the deflection of the face.
    pub fn set_deflection(&mut self, value: f64) {
        self.deflection = value;
    }

    /// Current status flags.
    pub fn status(&self) -> MeshStatus {
        self.status
    }

    /// Adds status flags.
    pub fn set_status(&mut self, value: MeshStatus) {
        self.status.set_status(value);
    }

    /// Removes status flags.
    pub fn unset_status(&mut self, value: MeshStatus) {
        self.status.unset_status(value);
    }

    /// True if the given flag is set.
    pub fn is_status(&self, value: MeshStatus) -> bool {
        self.status.is_set(value)
    }

    /// True when the face discrete model is valid — status is `NoError`,
    /// `ReMesh` or `UnorientedWire` (mirrors `IMeshData_Face::IsValid`).
    pub fn is_valid(&self) -> bool {
        self.status.is_equal(MeshStatus::NO_ERROR)
            || self.status.is_equal(MeshStatus::REMESH)
            || self.status.is_equal(MeshStatus::UNORIENTED_WIRE)
    }
}

/// Discrete model of a shape. Source: `BRepMeshData_Model.hxx` + `IMeshData_Model.hxx`.
#[derive(Clone, Default)]
pub struct MeshModel {
    shape: Option<TopoShape>,
    max_size: f64,
    faces: Vec<MeshFace>,
    edges: Vec<MeshEdge>,
    wires: Vec<MeshWire>,
}

impl MeshModel {
    /// Builds an empty discrete model for the given shape.
    pub fn new(shape: TopoShape) -> Self {
        Self {
            shape: Some(shape),
            max_size: 0.0,
            faces: Vec::new(),
            edges: Vec::new(),
            wires: Vec::new(),
        }
    }

    /// Shape the model was built for.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.shape.as_ref()
    }

    /// Replaces the reference shape.
    pub fn set_shape(&mut self, shape: TopoShape) {
        self.shape = Some(shape);
    }

    /// Maximum size of the shape's bounding box.
    pub fn max_size(&self) -> f64 {
        self.max_size
    }

    /// Sets the maximum size of the shape's bounding box.
    pub fn set_max_size(&mut self, value: f64) {
        self.max_size = value;
    }

    /// Number of discrete faces.
    pub fn faces_nb(&self) -> usize {
        self.faces.len()
    }

    /// Adds a new discrete face; returns its model index.
    pub fn add_face(&mut self, face: Face) -> usize {
        self.faces.push(MeshFace::new(face));
        self.faces.len() - 1
    }

    /// Discrete face with the given index.
    pub fn face(&self, index: usize) -> Result<&MeshFace, String> {
        self.faces
            .get(index)
            .ok_or_else(|| out_of_range("MeshModel::face", index, self.faces.len()))
    }

    /// Mutable access to the discrete face with the given index.
    pub fn face_mut(&mut self, index: usize) -> Result<&mut MeshFace, String> {
        if index >= self.faces.len() {
            return Err(out_of_range("MeshModel::face_mut", index, self.faces.len()));
        }
        Ok(&mut self.faces[index])
    }

    /// Removes and returns the face at the given index (shifts following indices).
    pub fn remove_face(&mut self, index: usize) -> Result<MeshFace, String> {
        if index >= self.faces.len() {
            return Err(out_of_range("MeshModel::remove_face", index, self.faces.len()));
        }
        Ok(self.faces.remove(index))
    }

    /// All discrete faces.
    pub fn faces(&self) -> &[MeshFace] {
        &self.faces
    }

    /// Number of discrete edges.
    pub fn edges_nb(&self) -> usize {
        self.edges.len()
    }

    /// Adds a new discrete edge; returns its model index.
    pub fn add_edge(&mut self, edge: Edge) -> usize {
        self.edges.push(MeshEdge::new(edge));
        self.edges.len() - 1
    }

    /// Discrete edge with the given index.
    pub fn edge(&self, index: usize) -> Result<&MeshEdge, String> {
        self.edges
            .get(index)
            .ok_or_else(|| out_of_range("MeshModel::edge", index, self.edges.len()))
    }

    /// Mutable access to the discrete edge with the given index.
    pub fn edge_mut(&mut self, index: usize) -> Result<&mut MeshEdge, String> {
        if index >= self.edges.len() {
            return Err(out_of_range("MeshModel::edge_mut", index, self.edges.len()));
        }
        Ok(&mut self.edges[index])
    }

    /// Removes and returns the edge at the given index (shifts following indices).
    pub fn remove_edge(&mut self, index: usize) -> Result<MeshEdge, String> {
        if index >= self.edges.len() {
            return Err(out_of_range("MeshModel::remove_edge", index, self.edges.len()));
        }
        Ok(self.edges.remove(index))
    }

    /// All discrete edges.
    pub fn edges(&self) -> &[MeshEdge] {
        &self.edges
    }

    /// Number of discrete wires.
    pub fn wires_nb(&self) -> usize {
        self.wires.len()
    }

    /// Adds a new discrete wire; returns its model index.
    pub fn add_wire(&mut self, wire: Wire) -> usize {
        self.wires.push(MeshWire::new(wire));
        self.wires.len() - 1
    }

    /// Discrete wire with the given index.
    pub fn wire(&self, index: usize) -> Result<&MeshWire, String> {
        self.wires
            .get(index)
            .ok_or_else(|| out_of_range("MeshModel::wire", index, self.wires.len()))
    }

    /// Mutable access to the discrete wire with the given index.
    pub fn wire_mut(&mut self, index: usize) -> Result<&mut MeshWire, String> {
        if index >= self.wires.len() {
            return Err(out_of_range("MeshModel::wire_mut", index, self.wires.len()));
        }
        Ok(&mut self.wires[index])
    }

    /// Removes and returns the wire at the given index (shifts following indices).
    pub fn remove_wire(&mut self, index: usize) -> Result<MeshWire, String> {
        if index >= self.wires.len() {
            return Err(out_of_range("MeshModel::remove_wire", index, self.wires.len()));
        }
        Ok(self.wires.remove(index))
    }

    /// All discrete wires.
    pub fn wires(&self) -> &[MeshWire] {
        &self.wires
    }

    /// Aggregated status mask over all faces, edges and wires.
    pub fn status_mask(&self) -> u32 {
        let mut mask = 0u32;
        for f in &self.faces {
            mask |= f.status.bits();
        }
        for e in &self.edges {
            mask |= e.status.bits();
        }
        for w in &self.wires {
            mask |= w.status.bits();
        }
        mask
    }

    /// True if any contained entity has the given status flag set.
    pub fn has_status(&self, value: MeshStatus) -> bool {
        self.status_mask() & value.bits() != 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpPln, GpPnt, GpPnt2d};

    /// Release registry entries for a shape tree so tests don't leave stale
    /// geometry keyed by a freed Arc address in the process-wide side-table.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn sample_edge(b: &TopoBuilder) -> Edge {
        b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0))
    }

    fn sample_wire(b: &TopoBuilder) -> Wire {
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        b.make_wire(&[e1, e2])
    }

    fn sample_face(b: &TopoBuilder) -> Face {
        b.make_face_plane(&GpPln::new(GpAx3::standard()))
    }

    #[test]
    fn mesh_status_bit_flags() {
        let mut s = MeshStatus::NO_ERROR;
        assert!(s.is_empty());
        assert_eq!(s.status_mask(), 0x0);
        s.set_status(MeshStatus::OPEN_WIRE);
        assert!(s.is_set(MeshStatus::OPEN_WIRE));
        assert!(!s.is_set(MeshStatus::FAILURE));
        s.set_status(MeshStatus::FAILURE);
        assert_eq!(s.status_mask(), 0x1 | 0x4);
        assert_eq!(s.union(MeshStatus::REMESH).status_mask(), 0x1 | 0x4 | 0x8);
        s.unset_status(MeshStatus::OPEN_WIRE);
        assert!(!s.is_set(MeshStatus::OPEN_WIRE));
        assert_eq!(s.status_mask(), 0x4);
    }

    #[test]
    fn model_add_remove_query() {
        let b = TopoBuilder::new();
        let mut model = MeshModel::new(TopoShape::new(crate::abs::ShapeType::Compound));
        assert_eq!(model.faces_nb(), 0);
        assert_eq!(model.edges_nb(), 0);
        assert_eq!(model.wires_nb(), 0);

        let e1 = sample_edge(&b);
        let e2 = sample_edge(&b);
        let ei0 = model.add_edge(e1);
        let ei1 = model.add_edge(e2);
        assert_eq!(model.edges_nb(), 2);

        let fi0 = model.add_face(sample_face(&b));
        assert_eq!(model.faces_nb(), 1);

        let wi0 = model.add_wire(sample_wire(&b));
        assert_eq!(model.wires_nb(), 1);

        // Wire references both edges; face references the wire.
        model.wire_mut(wi0).unwrap().add_edge(ei0, Orientation::Forward);
        model.wire_mut(wi0).unwrap().add_edge(ei1, Orientation::Reversed);
        model.face_mut(fi0).unwrap().add_wire(wi0);
        assert_eq!(model.wire(wi0).unwrap().edges_nb(), 2);
        assert_eq!(model.face(fi0).unwrap().wires_nb(), 1);
        assert_eq!(model.wire(wi0).unwrap().edge(0).unwrap(), ei0);
        assert_eq!(
            model.wire(wi0).unwrap().edge_orientation(1).unwrap(),
            Orientation::Reversed
        );

        // Query with out-of-range index fails.
        assert!(model.edge(99).is_err());
        assert!(model.face(99).is_err());

        // Removal returns the removed item and shrinks the collection.
        let removed = model.remove_edge(ei0).unwrap();
        assert_eq!(removed.pcurves_nb(), 0);
        assert_eq!(model.edges_nb(), 1);

        let removed_face = model.remove_face(fi0).unwrap();
        assert_eq!(removed_face.wires_nb(), 1);
        assert_eq!(model.faces_nb(), 0);

        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn model_status_aggregation() {
        let b = TopoBuilder::new();
        let mut model = MeshModel::new(TopoShape::new(crate::abs::ShapeType::Compound));
        assert_eq!(model.status_mask(), 0);

        let fi = model.add_face(sample_face(&b));
        let ei = model.add_edge(sample_edge(&b));
        model.face_mut(fi).unwrap().set_status(MeshStatus::FAILURE);
        model.edge_mut(ei).unwrap().set_status(MeshStatus::OPEN_WIRE);
        model.edge_mut(ei).unwrap().set_status(MeshStatus::UNORIENTED_WIRE);

        assert_eq!(model.status_mask(), 0x4 | 0x1 | 0x10);
        assert!(model.has_status(MeshStatus::FAILURE));
        assert!(model.has_status(MeshStatus::OPEN_WIRE));
        assert!(!model.has_status(MeshStatus::REUSED));

        // Unsetting a flag on one entity is reflected in the aggregation.
        model.edge_mut(ei).unwrap().unset_status(MeshStatus::OPEN_WIRE);
        assert!(!model.has_status(MeshStatus::OPEN_WIRE));

        // A face that is only ReMesh / UnorientedWire stays valid.
        model.face_mut(fi).unwrap().unset_status(MeshStatus::FAILURE);
        model.face_mut(fi).unwrap().set_status(MeshStatus::REMESH);
        assert!(model.face(fi).unwrap().is_valid());

        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn curve_and_pcurve_creation() {
        let mut c = MeshCurve::new();
        assert_eq!(c.parameters_nb(), 0);
        c.add_point(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        c.add_point(GpPnt::new(1.0, 0.0, 0.0), 1.0);
        c.add_point(GpPnt::new(2.0, 0.0, 0.0), 2.0);
        c.insert_point(1, GpPnt::new(0.5, 0.0, 0.0), 0.5).unwrap();
        assert_eq!(c.parameters_nb(), 4);
        assert_eq!(c.get_parameter(1).unwrap(), 0.5);
        assert!(c.get_point(1).unwrap().is_equal(&GpPnt::new(0.5, 0.0, 0.0)));
        assert_eq!(c.points().len(), c.parameters().len());

        // Out-of-range access is an error.
        assert!(c.get_point(99).is_err());

        // Clear keeps only the end points when requested.
        c.clear(true);
        assert_eq!(c.parameters_nb(), 2);
        assert_eq!(c.get_parameter(0).unwrap(), 0.0);
        assert_eq!(c.get_parameter(1).unwrap(), 2.0);
        c.clear(false);
        assert_eq!(c.parameters_nb(), 0);

        // Deflection attribute.
        c.set_deflection(0.01);
        assert_eq!(c.deflection(), 0.01);

        let mut pc = MeshPCurve::new(3, Orientation::Reversed);
        pc.add_point(GpPnt2d::new(0.0, 0.0), 0.0);
        pc.add_point(GpPnt2d::new(1.0, 1.0), 1.0);
        assert_eq!(pc.parameters_nb(), 2);
        assert_eq!(pc.get_point(1).unwrap(), GpPnt2d::new(1.0, 1.0));
        assert_eq!(pc.get_index(0).unwrap(), 0);
        *pc.get_index_mut(1).unwrap() = 7;
        assert_eq!(pc.get_index(1).unwrap(), 7);
        assert!(!pc.is_forward());
        assert!(!pc.is_internal());
        assert_eq!(pc.face(), 3);
        assert_eq!(pc.orientation(), Orientation::Reversed);
        assert_eq!(pc.indices().len(), pc.points().len());

        pc.set_deflection(0.02);
        assert_eq!(pc.deflection(), 0.02);

        let mut fwd = MeshPCurve::new(0, Orientation::Forward);
        assert!(fwd.is_forward());
    }

    #[test]
    fn mesh_edge_reads_topological_flags_and_curve() {
        let b = TopoBuilder::new();
        let e = sample_edge(&b);
        let mut me = MeshEdge::new(e);
        assert!(me.curve().is_some());
        assert_eq!(me.first_parameter(), 0.0);
        assert_eq!(me.last_parameter(), 1.0);
        assert!(me.same_param());
        assert!(me.same_range());
        assert!(!me.degenerated());
        assert!(me.is_free());
        assert_eq!(me.angular_deflection(), f64::MAX);

        // Add pcurves and look them up per face / orientation.
        let fi0 = me.add_pcurve(0, Orientation::Forward);
        me.add_pcurve(0, Orientation::Reversed);
        me.add_pcurve(1, Orientation::Forward);
        assert_eq!(me.pcurves_nb(), 3);
        assert_eq!(me.pcurve(fi0).unwrap().orientation(), Orientation::Forward);
        assert_eq!(me.pcurve_for(0, Orientation::Reversed).unwrap().orientation(), Orientation::Reversed);
        assert_eq!(me.pcurve_for(1, Orientation::Forward).unwrap().face(), 1);
        assert_eq!(me.pcurves_for(0), vec![0, 1]);
        assert!(!me.is_free());

        clear_tree(&me.edge().0);
    }
}
