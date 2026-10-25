//! `ChFi3d_FilletShape` — the three fillet surface families.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilletShape.hxx:18-33`.

/// `ChFi3d_FilletShape` (`ChFi3d_FilletShape.hxx:28-33`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletShape {
    /// `ChFi3d_Rational` — the default NURBS representation of circles.
    #[default]
    Rational,
    /// `ChFi3d_QuasiAngular` — NURBS whose parameter matches the circle's.
    QuasiAngular,
    /// `ChFi3d_Polynomial` — polynomial approximation of circles.
    Polynomial,
}
