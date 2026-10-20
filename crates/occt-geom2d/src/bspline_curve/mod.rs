//! 2D B-spline curve (rational and non-rational). Source: `Geom2d_BSplineCurve.hxx`

mod curve;
mod eval;
mod knot_ops;
mod periodic;
#[cfg(test)]
mod tests;

pub use self::curve::Geom2dBSplineCurve;

/// `Geom2d_BSplineCurve::MaxDegree()` = `BSplCLib::MaxDegree()`
/// (`BSplCLib.lxx:24-27`), same value as the 3D port's `MAX_DEGREE`.
pub(super) const MAX_DEGREE: usize = 25;
