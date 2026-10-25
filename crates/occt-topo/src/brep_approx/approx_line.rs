//! `BRepApprox_ApproxLine` — the line the fit reads: either three B-spline
//! curves of poles (3D, UV on surface 1, UV on surface 2) or a point list.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_ApproxLine.hxx:22-56`,
//! `.cxx:24-94` (OCCT 8.0.0).
//!
//! `Point(Index)` (`cxx:68-94`) prefers `myLineOn2S` whenever it holds points,
//! and otherwise reads the `Index`-th pole of the three curves. `NbPnts()`
//! (`cxx:46-66`) prefers the curves: the 3D one when present, then UV1, then
//! UV2, and only falls back to the point list. That asymmetry is why the
//! `MakeMLBetween` result line (`ApproxInt_MultiLine.gxx:510`) is built with the
//! point-list constructor: it must be read as points even when the original
//! line was curve-backed.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;

use crate::approx_int::LineSource;
use crate::int_tools_wline::PntOn2S;

/// `BRepApprox_ApproxLine` (`hxx:36-53`).
#[derive(Clone, Default)]
pub struct ApproxLine {
    /// `myCurveXYZ` (`cxx:30`).
    curve_xyz: Option<Arc<dyn Curve>>,
    /// `myCurveUV1`.
    curve_uv1: Option<Arc<dyn Curve2d>>,
    /// `myCurveUV2`.
    curve_uv2: Option<Arc<dyn Curve2d>>,
    /// `myLineOn2S` (`cxx:38`): `None` is OCCT's null handle.
    line_on2s: Option<Vec<PntOn2S>>,
}

impl ApproxLine {
    /// `BRepApprox_ApproxLine(CurveXYZ, CurveUV1, CurveUV2)` (`cxx:24-31`).
    pub fn new(
        curve_xyz: Option<Arc<dyn Curve>>,
        curve_uv1: Option<Arc<dyn Curve2d>>,
        curve_uv2: Option<Arc<dyn Curve2d>>,
    ) -> Self {
        Self {
            curve_xyz,
            curve_uv1,
            curve_uv2,
            line_on2s: None,
        }
    }

    /// `BRepApprox_ApproxLine(lin, false)` (`cxx:35-39`).
    pub fn from_line_on2s(points: Vec<PntOn2S>) -> Self {
        Self {
            curve_xyz: None,
            curve_uv1: None,
            curve_uv2: None,
            line_on2s: Some(points),
        }
    }

    /// `NbPnts()` (`cxx:46-66`).
    pub fn nb_pnts(&self) -> i32 {
        if let Some(c) = &self.curve_xyz {
            return c.bspline_poles().map_or(0, |p| p.len() as i32);
        }
        if let Some(c) = &self.curve_uv1 {
            return c.poles2d().map_or(0, |p| p.len() as i32);
        }
        if let Some(c) = &self.curve_uv2 {
            return c.poles2d().map_or(0, |p| p.len() as i32);
        }
        self.line_on2s.as_ref().map_or(0, |l| l.len() as i32)
    }

    /// `Point(Index)` (`cxx:68-94`): 1-based.
    pub fn point(&self, index: i32) -> PntOn2S {
        if let Some(line) = &self.line_on2s {
            if !line.is_empty() {
                return line[(index - 1) as usize];
            }
        }

        let p = self
            .curve_xyz
            .as_ref()
            .and_then(|c| c.bspline_poles())
            .and_then(|poles| poles.get((index - 1) as usize).copied())
            .unwrap_or_else(GpPnt::zero);
        let p1 = self
            .curve_uv1
            .as_ref()
            .and_then(|c| c.poles2d())
            .and_then(|poles| poles.get((index - 1) as usize).copied())
            .unwrap_or_else(|| GpPnt2d::new(0.0, 0.0));
        let p2 = self
            .curve_uv2
            .as_ref()
            .and_then(|c| c.poles2d())
            .and_then(|poles| poles.get((index - 1) as usize).copied())
            .unwrap_or_else(|| GpPnt2d::new(0.0, 0.0));

        // `IntSurf_PntOn2S aPntOn2S; aPntOn2S.SetValue(P, P1.X(), P1.Y(), P2.X(), P2.Y())`
        let mut point = PntOn2S {
            p: GpPnt::zero(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };
        point.set_value(p, p1.x(), p1.y(), p2.x(), p2.y());
        point
    }
}

impl LineSource for ApproxLine {
    fn nb_pnts(&self) -> i32 {
        ApproxLine::nb_pnts(self)
    }

    fn point(&self, index: i32) -> PntOn2S {
        ApproxLine::point(self, index)
    }

    fn from_points(points: Vec<PntOn2S>) -> Self {
        ApproxLine::from_line_on2s(points)
    }

    fn empty() -> Self {
        ApproxLine::from_line_on2s(Vec::new())
    }
}
