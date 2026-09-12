//! Analytic intersection of primitives. Port of the `IntAna` package
//! (TKGeomBase): `IntAna_Int3Pln`, `IntAna_QuadQuadGeo` (the analytic
//! quadric-quadric cases: plane-plane, plane-sphere, sphere-sphere,
//! plane-cylinder, plane-cone, cylinder-cylinder, cylinder-sphere,
//! sphere-cone, cone-cone, cylinder-cone, and torus pairs) and
//! `IntAna_IntLinTorus`.
//!
//! Shape-level intersection (`IntCurvesFace`, TKTopAlgo) is deliberately out
//! of scope; this module is the geometric kernel only.
mod prelude {

pub(crate) use std::cmp::Ordering;
pub(crate) use std::f64::consts::PI;

pub(crate) use occt_core::elib::clib;
pub(crate) use occt_core::gp::dir::DirAxis;
pub(crate) use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpDir2d, GpElips, GpHypr, GpLin,
    GpParab, GpPln, GpPnt, GpPnt2d, GpSphere, GpTorus, GpVec, GpVec2d,
};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION};

}

#[path = "../intana_torus.rs"]

mod torus;
pub use torus::{
    quadric_quadric_cone_torus, quadric_quadric_cylinder_torus, quadric_quadric_plane_torus,
    quadric_quadric_sphere_torus, quadric_quadric_torus_torus, TorusIntersection,
};

#[path = "../intana_quadric.rs"]
mod quadric;
#[path = "../intana_trig.rs"]
mod trig;
#[path = "../intana_curve.rs"]
mod curve_ana;
#[path = "../intana_intquadquad.rs"]
mod intquadquad;

pub use curve_ana::IntAnaCurve;
pub use intquadquad::IntQuadQuad;
pub use quadric::IntAnaQuadric;

mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
