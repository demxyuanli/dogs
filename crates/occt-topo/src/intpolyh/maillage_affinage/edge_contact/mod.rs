//! `CalculPtsInterTriEdgeCoplanaires` (`IntPolyh_MaillageAffinage.cxx:2051-2546`)
//! and `TriangleEdgeContact` (`cxx:2550-3138`): the two routines that turn a
//! triangle edge crossing the other triangle into start points.
//!
//! The two halves live in `coplanar` and `triangle_edge`; the second one is an
//! `impl MaillageAffinage` block, so nothing has to be re-exported.

mod coplanar;
mod triangle_edge;

pub(crate) use coplanar::calcul_pts_inter_tri_edge_coplanaires;
