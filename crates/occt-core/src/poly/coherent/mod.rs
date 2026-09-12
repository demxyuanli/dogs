//! Coherent triangulation: nodes / links / triangles with adjacency and
//! boundary-loop enumeration. Pure data structures.
//! Source: `Poly_CoherentNode.hxx`, `Poly_CoherentLink.hxx`,
//! `Poly_CoherentTriangle.hxx`, `Poly_CoherentTriPtr.hxx`,
//! `Poly_CoherentTriangulation.hxx`.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::slice::Iter;

pub(crate) use crate::gp::{GpPnt, GpPnt2d, GpXyz};
pub(crate) use crate::precision::{CONFUSION, INFINITE};

pub(crate) use super::super::triangulation_full::PolyTriangulation;

}


mod p01;
pub use p01::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
