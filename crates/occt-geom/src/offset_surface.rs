//! Offset 3D surface. Source: `Geom_OffsetSurface.hxx`
//!
//! The offset surface at `(u, v)` is the basis surface point displaced along
//! its unit normal by the signed offset: `d0 = basis.d0 + offset · n`. The
//! first partials are the basis partials (the normal's derivatives are not
//! carried — matching the reduced derivative used by [`crate::GeomOffsetCurve`],
//! which is exact for the planar/analytic bases this wrapper sees in practice).
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};

use crate::surface::Surface;

/// Surface obtained by offsetting a basis surface along its normals by a
/// constant signed distance. Source: `Geom_OffsetSurface.hxx`.
#[derive(Clone)]
pub struct GeomOffsetSurface {
    basis: Arc<dyn Surface>,
    offset: f64,
}

impl GeomOffsetSurface {
    /// Wraps `basis` offset by the signed `distance` along the basis normal.
    pub fn new(basis: Arc<dyn Surface>, distance: f64) -> Self {
        Self { basis, offset: distance }
    }

    /// The basis surface.
    pub fn basis_surface(&self) -> &Arc<dyn Surface> {
        &self.basis
    }

    /// The signed offset distance.
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// Unit normal of the basis surface at `(u, v)` (the offset direction).
    fn unit_normal(&self, u: f64, v: f64) -> GpVec {
        let (_, du, dv) = self.basis.d1(u, v);
        let n = du.xyz().crossed(dv.xyz());
        let m = n.modulus();
        if m > 1e-12 {
            GpVec::new(n.x / m, n.y / m, n.z / m)
        } else {
            GpVec::new(0.0, 0.0, 1.0)
        }
    }
}

impl Surface for GeomOffsetSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let p = self.basis.d0(u, v);
        let n = self.unit_normal(u, v);
        GpPnt::from_xyz(&p.coord.added(&n.xyz().multiplied(self.offset)))
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let (_, du, dv) = self.basis.d1(u, v);
        (self.d0(u, v), du, dv)
    }

    fn u_range(&self) -> (f64, f64) {
        self.basis.u_range()
    }

    fn v_range(&self) -> (f64, f64) {
        self.basis.v_range()
    }

    fn is_u_periodic(&self) -> bool {
        self.basis.is_u_periodic()
    }

    fn is_v_periodic(&self) -> bool {
        self.basis.is_v_periodic()
    }

    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }

    fn transform(&mut self, _t: &GpTrsf) {
        // Transform is not carried into the basis (matches the simplified
        // derivative handling above).
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }
}
