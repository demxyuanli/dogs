//! Variable-radius and chain rolling-ball fillets — a simplified port of
//! `BRepFilletAPI_MakeFillet` / `ChFi3d` with radius laws.
//!
//! Where `fillet_edge` replaces a sharp convex edge with a *constant-radius*
//! quarter-cylinder blend, this module lets the ball radius vary along the
//! edge (`VarFilletSpec` + `RadiusLaw`, mirroring OCCT's `Law_Function`), and
//! applies the fillet to a *chain* of edges (mirroring `ChFi3d`'s edge-chain
//! handling).
//!
//! The edge is sampled at `DEFAULT_VAR_SAMPLES` points. At each sample the
//! local tangency points, the blend-circle center and the radius are computed
//! with the same in-wedge geometry as `fillet_edge`; the blend face is a
//! degree-1 tensor-product B-spline (`occt_geom::bspline_surface::fit_surface_grid`)
//! interpolating the sampled cross-section arcs, bounded by the two tangency
//! polylines (shared with the trimmed adjacent faces) and the two end arcs
//! (shared with the rebuilt end faces). The adjacent planar faces are trimmed
//! with a *curved* clip boundary `R(t)·cot(θ/2)` along the tangency curve.
//!
//! The chain functions apply a fillet per edge sequentially (the result of one
//! step feeds the next). Consecutive edges that share a vertex meet when
//! `specs[k].r_start` equals `specs[k-1].r_end` at the shared vertex — the
//! caller encodes the shared-radius continuity; the code does not silently
//! rewrite a mismatch.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::f64::consts::SQRT_2;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpAx3, GpDir, GpPnt, GpSphere, GpVec};
pub(crate) use occt_geom::bspline_surface::fit_surface_grid;
pub(crate) use occt_geom::{GeomSphere, Surface};

pub(crate) use crate::brep_surface::face_plane;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::fillet_edge::faces_touching_edge;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex};
pub(crate) use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, wires_of_face,
};

}


mod spec;
mod variable_fillet;
mod end_faces;
pub use spec::*;
pub use variable_fillet::*;
pub use end_faces::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
