//! `ElSLib::*D3` for the surface adaptor.
//!
//! Source: `ModelingData/TKGeomBase/ElSLib/ElSLib.cxx:1147-1543` (OCCT 8.0.0),
//! ported in [`occt_core::elib::surface_eval_d3`]. This module only re-exports
//! the five elementary jets plus [`SurfD3`] so [`crate::Surface::d3`] can build
//! its dispatch without naming `occt-core` internals at the call site.
//!
//! `Geom_Surface::EvalD3` for the non-elementary surfaces is
//! [`crate::Surface::try_eval_dn`]; there is no `ElSLib` entry for them.

pub use occt_core::elib::surface_eval_d3::{
    cone_d3, cylinder_d3, plane_d3, sphere_d3, torus_d3, SurfD3,
};
