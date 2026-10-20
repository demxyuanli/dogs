//! `Default` and the `DiscretRoot` accessor impl, plus the free-function
//! entry point `incremental_mesh_to_shape_mesh`.

use super::*;

impl Default for IncrementalMesh {
    fn default() -> Self {
        Self::new()
    }
}

// DiscretRoot accessors re-exposed on IncrementalMesh for callers that do not
// want to reach through `.root`.
impl IncrementalMesh {
    /// Set the shape to triangulate (before [`IncrementalMesh::perform`]).
    pub fn set_shape(&mut self, shape: &TopoShape) {
        self.root.set_shape(shape);
        self.modified = false;
    }

    /// The shape being triangulated, if set.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.root.shape()
    }

    /// `true` if triangulation was performed and succeeded.
    pub fn is_done(&self) -> bool {
        self.root.is_done()
    }
}

/// Mesh a shape into a [`ShapeMesh`] through the OCCT-style BRepMesh pipeline
/// with the given linear deflection.
///
/// Convenience entry wrapping [`IncrementalMesh::perform`]; the returned
/// [`ShapeMesh`] holds the assembled triangle soup (vertices + triangles).
pub fn incremental_mesh_to_shape_mesh(shape: &TopoShape, deflection: f64) -> Result<ShapeMesh, String> {
    let mut inc = IncrementalMesh::new();
    inc.set_shape(shape);
    inc.change_parameters().deflection = deflection;
    inc.perform()
}
