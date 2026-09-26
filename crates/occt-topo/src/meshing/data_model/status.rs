use super::prelude::*;


/// OCCT `RealLast()` — largest finite double, used as the "unset" deflection
/// sentinel (matches `IMeshData_TessellatedShape` default).

pub(super) const REAL_LAST: f64 = f64::MAX;

/// Out-of-range error message, mirroring OCCT's `Standard_OutOfRange`.
pub(super) fn out_of_range(op: &str, index: usize, len: usize) -> String {
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
    pub(super) points: Vec<GpPnt>,
    pub(super) parameters: Vec<f64>,
    pub(super) deflection: f64,
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
    pub(super) face: usize,
    pub(super) orientation: Orientation,
    pub(super) points: Vec<GpPnt2d>,
    pub(super) parameters: Vec<f64>,
    pub(super) indices: Vec<i32>,
    pub(super) deflection: f64,
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

    /// `IMeshData_PCurve::IsForward` (`hxx:52`): anything but `REVERSED`.
    pub fn is_forward(&self) -> bool {
        self.orientation != Orientation::Reversed
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
    pub(super) edge: Edge,
    pub(super) curve: Option<Arc<dyn Curve>>,
    pub(super) first: f64,
    pub(super) last: f64,
    pub(super) discretization: MeshCurve,
    pub(super) pcurves: Vec<MeshPCurve>,
    pub(super) same_param: bool,
    pub(super) same_range: bool,
    pub(super) degenerated: bool,
    pub(super) angular_deflection: f64,
    pub(super) deflection: f64,
    pub(super) status: MeshStatus,
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

    /// Mutable pcurve with the given index (populated by the edge-discretization
    /// step).
    pub fn pcurve_mut(&mut self, index: usize) -> Result<&mut MeshPCurve, String> {
        let len = self.pcurves.len();
        self.pcurves
            .get_mut(index)
            .ok_or_else(|| out_of_range("MeshEdge::pcurve_mut", index, len))
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
    pub(super) wire: Wire,
    pub(super) edges: Vec<usize>,
    pub(super) orientations: Vec<Orientation>,
    /// `ShapeAnalysis_WireOrder::Ordered < 0`: walk this slot's CheckOrder
    /// pcurve backwards. Not `Edge.Reverse()` — on a seam that selects PCurve2
    /// instead of reversing the chord (`ShapeExtend_WireData::Edge(signed)`).
    pub(super) reverse_walk: Vec<bool>,
    pub(super) deflection: f64,
    pub(super) status: MeshStatus,
}

impl MeshWire {
    /// Builds an empty discrete wire from the topological wire.
    pub fn new(wire: Wire) -> Self {
        Self {
            wire,
            edges: Vec::new(),
            orientations: Vec::new(),
            reverse_walk: Vec::new(),
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
        self.add_edge_ordered(edge_index, orientation, false)
    }

    /// `AddEdge` plus the `Ordered` sign: `reverse` walks the same pcurve
    /// backwards (`WireOrder` negative index).
    pub fn add_edge_ordered(
        &mut self,
        edge_index: usize,
        orientation: Orientation,
        reverse: bool,
    ) -> usize {
        self.edges.push(edge_index);
        self.orientations.push(orientation);
        self.reverse_walk.push(reverse);
        self.edges.len() - 1
    }

    /// True when `Ordered` walks this slot backwards.
    pub fn edge_reverse_walk(&self, index: usize) -> bool {
        self.reverse_walk.get(index).copied().unwrap_or(false)
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
