//! `BRepMesh_EdgeDiscret` — edge discretization driver.

// ---------------------------------------------------------------------------
// EdgeDiscret
// ---------------------------------------------------------------------------

use super::prelude::*;
use super::*;

/// Discretizes 3D edges into deflection-bounded point sequences.
///
/// Port of `BRepMesh_EdgeDiscret`: computes a per-edge deflection, creates a
/// [`CurveTessellator`] and emits the 3D polyline (optionally mirrored into 2D
/// pcurve points by [`EdgeParameterProvider`]).
#[derive(Debug, Clone)]
pub struct EdgeDiscret {
    params: MeshParameters,
}

impl EdgeDiscret {
    pub fn new(params: MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &MeshParameters {
        &self.params
    }

    /// Set the parameters.
    pub fn set_parameters(&mut self, params: MeshParameters) {
        self.params = params;
    }

    /// EdgeDeflection used for the given edge: the parameters' deflection, adjusted
    /// for a relative deflection against `max_shape_size`.
    pub fn edge_deflection(&self, max_shape_size: f64) -> f64 {
        if self.params.relative {
            EdgeDeflection::compute_absolute_deflection(self.params.deflection, max_shape_size)
        } else {
            self.params.deflection
        }
    }

    /// Create a tessellator for `curve` over `[first, last]`, seeding it with
    /// the minimum point count in effect. Port of `CreateEdgeTessellator`.
    pub fn create_edge_tessellator(
        &self,
        curve: &Arc<dyn Curve>,
        first: f64,
        last: f64,
    ) -> CurveTessellator {
        CurveTessellator::from_range(
            curve.clone(),
            first,
            last,
            self.edge_deflection(0.0),
            2,
        )
    }

    /// Discretize a 3D curve into a deflection-bounded polyline.
    ///
    /// Arc-length/parameter adaptive: the curve is seeded uniformly and refined
    /// while any chord deviates from the curve by more than `deflection`.
    /// Endpoints are always present.
    pub fn discretize_edge(
        curve: &Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
    ) -> Vec<GpPnt> {
        CurveTessellator::from_range(curve.clone(), first, last, deflection, 2).points
    }

    /// Discretize an edge of the mesh model into a 3D polyline using the edge's
    /// own deflection. Degenerated edges collapse to their endpoints.
    pub fn discretize_edge_mesh(&self, edge: &DiscreteMeshEdge) -> Vec<GpPnt> {
        let def = if edge.deflection > 0.0 {
            edge.deflection
        } else {
            self.edge_deflection(0.0)
        };
        let tess = CurveTessellator::from_range(
            edge.curve.clone(),
            edge.first,
            edge.last,
            def,
            2,
        );
        if edge.degenerated {
            // Degenerated edges contribute their start point only.
            return tess.points.first().copied().map(|p| vec![p]).unwrap_or_default();
        }
        tess.points
    }

    /// Port of `Tessellate3d(theUpdateEnds=true)`: the first and last points are
    /// the edge's vertex points (`BRep_Tool::Pnt(firstVertex)` /
    /// `Pnt(lastVertex)`), the interior points come from the tessellator
    /// (OCCT's indices 2..PointsNb−1). Degenerated edges keep only the start
    /// vertex.
    pub fn tessellate_3d(
        edge: &DiscreteMeshEdge,
        tessellator: &CurveTessellator,
        first_vertex: GpPnt,
        last_vertex: GpPnt,
    ) -> Vec<GpPnt> {
        let n = tessellator.points_nb();
        if n == 0 {
            return Vec::new();
        }
        let mut pts = Vec::with_capacity(n);
        pts.push(first_vertex);
        if !edge.degenerated {
            // OCCT: `for (i = 2; i < PointsNb(); ++i)` — interior points only.
            for i in 1..n.saturating_sub(1) {
                if let Some(p) = tessellator.point(i) {
                    if pts.last().is_none_or(|q: &GpPnt| q.distance(&p) > 1e-12) {
                        pts.push(p);
                    }
                }
            }
        }
        pts.push(last_vertex);
        pts
    }

    /// Port of `Tessellate2d(theUpdateEnds=true)`: evaluate the 2D pcurve at the
    /// 3D tessellation parameters (mapped through the provider). The provider
    /// passes a SameParameter edge's parameters through verbatim and linearly
    /// rescales a non-SameParameter edge's stored parameters onto the pcurve
    /// range (`BRepMesh_EdgeParameterProvider::Parameter` without the 3D-point
    /// projection — see `EdgeParameterProvider`).
    pub fn tessellate_2d(
        provider: &EdgeParameterProvider,
        tessellator: &CurveTessellator,
        pcurve: &dyn Fn(f64) -> GpPnt2d,
    ) -> Vec<GpPnt2d> {
        let mut out = Vec::with_capacity(tessellator.points_nb());
        for &u in tessellator.params() {
            let actual = provider.remap(u);
            let p2 = pcurve(actual);
            if out.last().is_none_or(|q: &GpPnt2d| q.distance(&p2) > 1e-12) {
                out.push(p2);
            }
        }
        out
    }
}
