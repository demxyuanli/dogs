use super::prelude::*;
use super::*;

/// Discrete model of a face. Source: `BRepMeshData_Face.hxx` + `IMeshData_Face.hxx`.
#[derive(Clone)]
pub struct MeshFace {
    pub(super) face: Face,
    pub(super) surface: Option<Arc<dyn Surface>>,
    pub(super) wires: Vec<usize>,
    pub(super) points: Vec<GpPnt>,
    pub(super) deflection: f64,
    pub(super) status: MeshStatus,
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
    pub(super) shape: Option<TopoShape>,
    pub(super) max_size: f64,
    pub(super) faces: Vec<MeshFace>,
    pub(super) edges: Vec<MeshEdge>,
    pub(super) wires: Vec<MeshWire>,
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
