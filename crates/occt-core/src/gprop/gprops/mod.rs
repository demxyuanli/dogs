//! GProp — global properties framework. Source: `GProp/`.
//!
//! Ports the OCCT `GProp` package: the composite global-properties
//! accumulator [`GProps`] (`GProp_GProps`), the principal-properties
//! presentation [`PrincipalProps`] (`GProp_PrincipalProps`), the element
//! frameworks `PGProps` / `SelGProps` / `VelGProps` / `CelGProps`, the plane
//! equation [`PEquation`], the Huygens operator [`h_operator`] and the
//! [`ValueType`] enumeration.
//!
//! `GProps` composes the global properties (mass / length / area / volume,
//! centre of mass, quadratic inertia matrix) of a *compound geometric system*.
//! Elementary pieces are added either directly ([`GProps::add_point_mass`]) or
//! through the element frameworks, each of which accumulates one family of
//! geometric cells into a `GProps`:
//!
//! * [`PGProps`] — point set (`add_point`);
//! * [`SelGProps`] — surface (triangles);
//! * [`VelGProps`] — solid (tetrahedra / boxes / octahedra);
//! * [`CelGProps`] — curve (segments).
//!
//! All inertia is accumulated about the framework's reference point `loc`
//! (origin by default). [`GProps::matrix_of_inertia`] shifts it to the centre
//! of mass by the parallel-axis (Huygens) theorem, and [`h_operator`] provides
//! the parallel-axis term `m·(d²δ − d·dᵀ)` used everywhere.
//!
//! Each element contribution is accumulated exactly: a point mass uses
//! `r²δ − r·rᵀ`, a segment / triangle / tetrahedron uses the closed-form
//! integral of `r⊗r` over the cell (barycentric formulas), so a box
//! decomposed into tetrahedra reproduces the analytic box inertia tensor
//! exactly (up to floating-point rounding).
mod prelude {

pub(crate) use crate::gp::{GpAx1, GpAx3, GpDir, GpMat, GpPln, GpPnt, GpVec, GpXyz};
pub(crate) use crate::precision::CONFUSION;
pub(crate) use std::ops::{Deref, DerefMut};

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
