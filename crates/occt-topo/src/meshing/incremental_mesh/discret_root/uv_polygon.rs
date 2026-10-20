//! Face boundary UV polygon assembly (`EdgeDiscret` step).

use super::*;

impl IncrementalMesh {
    /// Build the UV polygon of a model face from its boundary edges.
    ///
    /// Each wire edge is discretized in 3D ([`Self::discretize_edge`],
    /// `EdgeDiscret`) and its polyline is projected onto the face surface to
    /// give the boundary UV points; edges are stitched so the polygon is a
    /// continuous closed chain regardless of the wire's edge orientations.
    ///
    /// The first edge's direction is ambiguous (the wire stores every edge in
    /// its natural curve direction, not necessarily the traversal direction), so
    /// the chain is built twice — natural and flipped — and the closed one is
    /// kept.
    pub(super) fn build_face_uv_polygon(&self, model: &MeshModel, face_index: usize) -> Result<UvFace, String> {
        let face = model
            .face(face_index)
            .map_err(|e| format!("IncrementalMesh::build_face_uv_polygon: {e}"))?;
        let surface = face.surface().ok_or("IncrementalMesh::build_face_uv_polygon: face has no surface")?;

        let mut outer_wire: Vec<GpPnt2d> = Vec::new();
        let mut inner_wires: Vec<Vec<GpPnt2d>> = Vec::new();
        for (w_pos, &wire_index) in face.wires().iter().enumerate() {
            let mut chain = self.wire_uv_chain(model, wire_index, surface.as_ref(), false);
            if !chain_closed(&chain) {
                chain = self.wire_uv_chain(model, wire_index, surface.as_ref(), true);
            }
            if chain.len() >= 2 {
                if w_pos == 0 {
                    outer_wire = chain;
                } else {
                    inner_wires.push(chain);
                }
            }
        }

        if outer_wire.is_empty() {
            return Err(format!(
                "IncrementalMesh::build_face_uv_polygon: face {face_index} has no boundary UV points"
            ));
        }
        Ok(UvFace {
            outer_wire,
            inner_wires,
            deflection: face.deflection(),
            id: face_index,
        })
    }

    /// UV boundary chain of one wire on a face surface.
    ///
    /// Each edge is discretized in 3D and projected to UV, then stitched into a
    /// continuous chain. `reverse_first` flips the first edge's direction,
    /// resolving the wire-start ambiguity.
    pub(super) fn wire_uv_chain(
        &self,
        model: &MeshModel,
        wire_index: usize,
        surface: &dyn Surface,
        reverse_first: bool,
    ) -> Vec<GpPnt2d> {
        let wire = match model.wire(wire_index) {
            Ok(w) => w,
            Err(_) => return Vec::new(),
        };
        let mut chain: Vec<GpPnt2d> = Vec::new();
        for j in 0..wire.edges_nb() {
            let edge_index = match wire.edge(j) {
                Ok(e) => e,
                Err(_) => return Vec::new(),
            };
            let edge = match model.edge(edge_index) {
                Ok(e) => e,
                Err(_) => return Vec::new(),
            };
            let mut uv = self.edge_uv_points(edge, surface);
            if reverse_first && j == 0 {
                uv.reverse();
            }
            stitch_chain(&mut chain, uv);
        }
        chain
    }

    /// UV points of one edge on a face.
    ///
    /// The 3D edge is discretized into a deflection-bounded polyline
    /// ([`Self::discretize_edge`]) and every polyline point is projected onto
    /// the face surface.
    pub(super) fn edge_uv_points(&self, edge: &MeshEdge, surface: &dyn Surface) -> Vec<GpPnt2d> {
        let deflection = edge.deflection();
        let deflection = if deflection.is_finite() && deflection > 0.0 && deflection < f64::MAX {
            deflection
        } else {
            self.parameters.deflection
        };
        self.discretize_edge(edge.edge(), deflection)
            .into_iter()
            .map(|p| {
                // Analytic ProjLib inverse first (exact even on unbounded-v
                // surfaces); Newton fallback for non-analytic (B-spline) faces.
                crate::pcurve_full::surface_project_point(surface, &p).unwrap_or_else(|| {
                    let (u, v) = project_uv(surface, &p);
                    GpPnt2d::new(u, v)
                })
            })
            .collect()
    }

}
