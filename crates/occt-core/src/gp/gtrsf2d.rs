//! 2D general (affine) transform. Source: `gp_GTrsf2d.hxx/.cxx`.
//!
//! Unlike `gp_Trsf2d` this may be non-orthogonal (an affinity); it is the type
//! `ShapeBuild_Edge::TransformPCurve` uses to scale the U parameter of a pcurve
//! (`ShapeBuild_Edge.cxx:620-698`).

use crate::gp::{GpAx2d, GpMat2d, GpPnt2d, GpXY};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpGTrsf2d {
    pub matrix: GpMat2d,
    pub loc: GpXY,
}

impl GpGTrsf2d {
    pub fn identity() -> Self {
        Self { matrix: GpMat2d::identity(), loc: GpXY::zero() }
    }

    /// `gp_GTrsf2d::SetAffinity(A, Ratio)` (`gp_GTrsf2d.cxx:24-38`).
    pub fn set_affinity(&mut self, a: &GpAx2d, ratio: f64) {
        let (ax, ay) = (a.direction().x, a.direction().y);
        self.matrix.data[0][0] = (1.0 - ratio) * ax * ax + ratio;
        self.matrix.data[1][1] = (1.0 - ratio) * ay * ay + ratio;
        self.matrix.data[0][1] = (1.0 - ratio) * ax * ay;
        self.matrix.data[1][0] = self.matrix.data[0][1];
        let base = *a.location().xy();
        let mut l = base;
        l.reverse();
        l.multiply_mat2d(&self.matrix);
        l.add(&base);
        self.loc = l;
    }

    /// `gp_GTrsf2d::Transforms(gp_XY&)`: `coord = matrix * coord + loc`.
    pub fn transforms(&self, p: &GpPnt2d) -> GpPnt2d {
        let mut xy = *p.xy();
        xy.multiply_mat2d(&self.matrix);
        xy.add(&self.loc);
        GpPnt2d::from_xy(xy)
    }
}
