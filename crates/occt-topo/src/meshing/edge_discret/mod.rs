//! Port of OCCT edge_discret — Wave 1 BRepMesh.
//!
//! `EdgeDiscret` turns a 3D edge into a deflection-bounded polyline
//! (`BRepMesh_EdgeDiscret` + `BRepMesh_CurveTessellator` + `GCPnts_*EdgeDeflection`),
//! `EdgeParameterProvider` maps stored edge parameters to the actual pcurve
//! parameters (`BRepMesh_EdgeParameterProvider`), and `CurveTessellator`
//! flattens an `Arc<dyn Curve>` with adaptive chord-deviation control
//! (`BRepMesh_CurveTessellator`).
//!
//! `MeshParameters` is imported from the sibling `parameters.rs`; `DiscreteMeshEdge`,
//! `MeshFace` and `DiscreteMeshModel` are local UV-polygon stand-ins (see below) until
//! the `data_model.rs` topology-indexed types are stable.

mod prelude {
    pub(crate) use std::sync::Arc;
    pub(crate) use occt_core::elib::{clib, slib};
    pub(crate) use occt_core::gcpnts::{perform_linear, perform_tangential_curve, CurveSecondDeriv};
    pub(crate) use occt_core::gp::{
        GpAx1, GpAx2, GpCirc, GpDir, GpDir2d, GpElips, GpHypr, GpLin, GpParab, GpPln, GpPnt, GpPnt2d,
        GpVec,
    };
    pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
    pub(crate) use occt_geom::{geom_api, Curve, Surface};
    pub(crate) use occt_geom2d::curve::Curve2d;
    pub(crate) use crate::abs::{Orientation, ShapeType};
    pub(crate) use crate::brep_tool::BRepTool;
    pub(crate) use crate::shape::{Edge, Vertex};
}

use prelude::*;

// Meshing parameters — owned by the sibling `parameters.rs` stub.
pub use super::parameters::MeshParameters;

// ---------------------------------------------------------------------------
// Sibling-contract stand-ins
// ---------------------------------------------------------------------------
// ponytail: temporary UV-polygon stand-ins mirroring the sibling-stub contract.
// The real `data_model::{DiscreteMeshEdge,MeshFace}` are topology-index based and not
// stable yet; swap to `use super::data_model::*;` once the discretizers are
// rewritten against the model's wire/edge/pcurve collections.

/// Discrete edge carrying a 3D curve + parameter range + discretization data.
#[derive(Clone)]
pub struct DiscreteMeshEdge {
    pub curve: Arc<dyn Curve>,
    pub first: f64,
    pub last: f64,
    pub deflection: f64,
    pub angular_deflection: f64,
    pub same_param: bool,
    pub same_range: bool,
    pub degenerated: bool,
    /// Discretized 2D pcurves (one entry per face the edge bounds).
    pub pcurves: Vec<Vec<GpPnt2d>>,
}

impl DiscreteMeshEdge {
    pub fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve,
            first,
            last,
            deflection: 0.1,
            angular_deflection: 0.5,
            same_param: true,
            same_range: true,
            degenerated: false,
            pcurves: Vec::new(),
        }
    }
}

/// Discrete face — boundary wires as UV polygons + mesh control values.
#[derive(Debug, Clone)]
pub struct UvFace {
    /// Outer wire as a closed UV polygon (first == last point).
    pub outer_wire: Vec<GpPnt2d>,
    /// Inner (hole) wires as closed UV polygons.
    pub inner_wires: Vec<Vec<GpPnt2d>>,
    pub deflection: f64,
    pub id: usize,
}

impl UvFace {
    pub fn new(outer_wire: Vec<GpPnt2d>) -> Self {
        Self { outer_wire, inner_wires: Vec::new(), deflection: 0.1, id: 0 }
    }
}

/// Discrete model — the mesh data shared between the two discretizers.
#[derive(Default, Clone)]
pub struct DiscreteMeshModel {
    pub edges: Vec<DiscreteMeshEdge>,
    pub faces: Vec<UvFace>,
    pub max_size: f64,
}
// ---------------------------------------------------------------------------
// EdgeParameterProvider
// ---------------------------------------------------------------------------

/// Maps stored polygon parameters to the actual curve parameters.
///
/// Port of `BRepMesh_EdgeParameterProvider`. When an edge is *not* SameParameter
/// the stored parameters are rescaled linearly onto `[first, last]`; for
/// SameParameter edges stored parameters are used verbatim. Also provides the
/// uniform parameter grids used to seed edge discretization.
#[derive(Debug, Clone)]
pub struct EdgeParameterProvider {
    first: f64,
    last: f64,
    is_same_param: bool,
    old_first: f64,
    old_last: f64,
    scale: f64,
    /// `myCurParam` / `myFoundParam` (`hxx:86, 119-136`).
    cur_param: f64,
    found_param: f64,
}

impl EdgeParameterProvider {
    /// Provider over `[first, last]` assuming the edge is SameParameter.
    pub fn new(first: f64, last: f64) -> Self {
        Self {
            first,
            last,
            is_same_param: true,
            old_first: first,
            old_last: last,
            scale: 1.0,
            cur_param: first,
            found_param: first,
        }
    }

    /// Provider aware of the stored parameter range `(old_first, old_last)`.
    /// The linear scale factor `(last - first) / (old_last - old_first)` maps a
    /// stored parameter onto the actual range, mirroring the OCCT constructor.
    pub fn with_stored(first: f64, last: f64, old_first: f64, old_last: f64) -> Self {
        let denom = old_last - old_first;
        let scale = if denom.abs() > f64::EPSILON && first.is_finite() && last.is_finite() {
            (last - first) / denom
        } else {
            1.0
        };
        Self {
            first,
            last,
            is_same_param: false,
            old_first,
            old_last,
            scale,
            cur_param: first,
            found_param: first,
        }
    }

    /// First parameter of the actual range.
    pub fn first(&self) -> f64 {
        self.first
    }

    /// Last parameter of the actual range.
    pub fn last(&self) -> f64 {
        self.last
    }

    /// `(first, last)` of the actual range.
    pub fn range(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// Whether the edge is SameParameter.
    pub fn is_same_param(&self) -> bool {
        self.is_same_param
    }

    /// Linear scale factor between the stored and the actual parameter ranges.
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Length of the actual parameter range.
    pub fn length(&self) -> f64 {
        self.last - self.first
    }

    /// Clamp `u` into `[first, last]`.
    pub fn clamp(&self, u: f64) -> f64 {
        if self.last >= self.first {
            u.clamp(self.first, self.last)
        } else {
            u
        }
    }

    /// Map a stored parameter onto the actual range. SameParameter edges pass
    /// the value through unchanged.
    ///
    /// Projection-free linear rescale: the `Parameter` result before the
    /// `Extrema_LocateExtPC` refinement. Kept for callers that carry only 2D
    /// data (`tessellate_2d`, `parameter`). Callers that do have the 3D point
    /// and the curve-on-surface must use `parameter_of` below, which implements
    /// the full `BRepMesh_EdgeParameterProvider::Parameter` refinement and IS
    /// the path taken by the live pipeline (`incremental_mesh/discret_root.rs:921,1229`).
    pub fn remap(&self, stored: f64) -> f64 {
        if self.is_same_param {
            return stored;
        }
        self.first + self.scale * (stored - self.old_first)
    }

    /// `BRepMesh_EdgeParameterProvider::Parameter` (`hxx:109-139`).
    /// SameParameter returns `stored`. Otherwise scale, then
    /// `Extrema_LocateExtPC` on the face `BRepAdaptor_Curve` (curve-on-surface)
    /// with the period/regression guard. The projector is initialized with
    /// `Precision::PConfusion()` (`hxx:99-103`), which drives the iteration stop
    /// criterion `abs(Un - Un-1) < TolF` (`Extrema_GLocateExtPC.hxx:71-72`).
    pub fn parameter_of(&mut self, stored: f64, point3d: &GpPnt, cos: &dyn Curve) -> f64 {
        if self.is_same_param {
            return stored;
        }
        let prev_param = self.cur_param;
        self.cur_param = self.first + self.scale * (stored - self.old_first);
        let prev_found = self.found_param;
        self.found_param += self.cur_param - prev_param;
        if let Some((found, _)) =
            crate::int_tools_vertex_line::extrema_locate_ext_pc_with_tol(
                cos,
                point3d,
                self.found_param,
                self.first,
                self.last,
                PCONFUSION,
            )
        {
            if (prev_found < self.found_param && prev_found < found)
                || (prev_found > self.found_param && prev_found > found)
            {
                self.found_param = found;
            }
        }
        self.found_param
    }

    /// `n` parameters uniformly spaced across the actual range, inclusive of
    /// both endpoints. `n` is clamped to at least 2.
    pub fn uniform_parameters(&self, n: usize) -> Vec<f64> {
        let n = n.max(2);
        let span = self.length();
        (0..n)
            .map(|i| self.first + span * i as f64 / (n - 1) as f64)
            .collect()
    }

    /// `n - 1` parameter segments spanning `[first, last]`.
    pub fn segments(&self, n: usize) -> Vec<(f64, f64)> {
        self.uniform_parameters(n).windows(2).map(|w| (w[0], w[1])).collect()
    }

    /// Parameter of the `index`-th point of an `nb`-point grid. When the edge is
    /// not SameParameter the stored grid is remapped onto the actual range.
    pub fn parameter(&self, index: usize, nb: usize) -> f64 {
        let n = nb.max(2);
        let i = index.min(n - 1);
        let stored = self.old_first + (self.old_last - self.old_first) * i as f64 / (n - 1) as f64;
        self.remap(stored)
    }
}

mod curve_on_surface;
mod curve_tessellator;
mod discret;
mod edge_deflection;
mod split_2d;
#[cfg(test)]
mod tests;

pub use self::curve_on_surface::*;
pub use self::curve_tessellator::*;
pub use self::discret::*;
pub use self::edge_deflection::*;
pub use self::split_2d::*;

fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}
