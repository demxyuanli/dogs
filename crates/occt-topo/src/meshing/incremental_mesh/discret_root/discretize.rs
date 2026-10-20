//! Public edge/face discretization helpers.

use super::*;

impl IncrementalMesh {
    /// Discretize a 3D edge into a deflection-bounded polyline.
    ///
    /// Port of `BRepMesh_EdgeDiscret`/`GCPnts_UniformDeflection`: the edge's
    /// curve is sampled adaptively so every chord deviates from the curve by at
    /// most `deflection`. Both endpoints are always included. Returns an empty
    /// polyline when the edge has no registered curve or an unbounded range.
    pub fn discretize_edge(&self, edge: &Edge, deflection: f64) -> Vec<GpPnt> {
        let Some(curve) = BRepTool::edge_curve_adaptor_world(edge) else {
            return Vec::new();
        };
        let (first, last) = BRepTool::edge_parameters(edge);
        if !(first.is_finite() && last.is_finite() && last >= first) {
            return Vec::new();
        }
        EdgeDiscret::discretize_edge(&curve, first, last, deflection.max(1e-9))
    }

    /// Discretize a single face into a 3D triangle soup through the OCCT-style
    /// pipeline (a temporary single-face model is triangulated).
    ///
    /// The OCCT equivalent is `BRepMesh_IncrementalMesh(aFace, theDeflection,
    /// false, theAngle)`: that constructor stores the deflection in
    /// `myParameters` and calls `Perform`
    /// (`BRepMesh_IncrementalMesh.cxx:42-57`), which runs `initParameters()`
    /// once (`:86`) and then hands the same `myParameters` to the whole
    /// `IMeshTools_Context` chain (`:90-91`). So one initialised parameter set
    /// drives *every* stage below — `ModelBuilder`, the boundary/heal/pre-process
    /// stages inside `build_shape_mesh`, and the face triangulation. Reading the
    /// driver's own `self.parameters` (the caller's, still carrying the
    /// `-1` sentinels that `initParameters` is meant to resolve) for part of the
    /// chain is not something OCCT does.
    ///
    /// Returns an empty soup for invalid parameters or an untriangulatable face:
    /// OCCT has no fallback tessellator, a face whose model cannot be built or
    /// whose Delaunay fails keeps `IMeshData_Failure`
    /// (`BRepMesh_BaseMeshAlgo.cxx:40-62`) and yields no triangles.
    pub fn discretize_face(&self, face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
        let mut params = self.parameters.clone();
        params.deflection = deflection;
        // `BRepMesh_IncrementalMesh::initParameters` (`BRepMesh_IncrementalMesh.hxx:79-108`).
        if init_parameters(&mut params).is_err() {
            return (Vec::new(), Vec::new());
        }
        let Ok(mut model) = ModelBuilder::build_model(&face.0, &params) else {
            return (Vec::new(), Vec::new());
        };
        // `IMeshTools_Context::PreProcessModel` (`hxx:85`) runs inside `build_shape_mesh`
        // through `prepare_model_boundaries`, which is the only place
        // `IncrementalMesh::perform` runs it too. Running it here as well would
        // execute the seam-edge amplifier twice on the same model.
        // Same driver shape as [`Self::build_boundary_model`]: the whole chain
        // must see the initialised set, so run it on a driver carrying that set
        // instead of on `self`.
        let mut driver = IncrementalMesh { parameters: params, ..IncrementalMesh::new() };
        match driver.build_shape_mesh(&mut model) {
            Ok(m) => (m.vertices, m.triangles),
            Err(_) => (Vec::new(), Vec::new()),
        }
    }

}
