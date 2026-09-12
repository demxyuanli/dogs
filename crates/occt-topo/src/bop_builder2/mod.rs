//! General Fuse builder + Boolean wrapper — the building phase of the exact
//! NURBS boolean (Phase 20).
//!
//! Port of `BOPAlgo_Builder` (`BOPAlgo_Builder.cxx` + `BOPAlgo_Builder_1.cxx`)
//! and `BOPAlgo_BOP` (`BOPAlgo_BOP.hxx/.cxx`) from TKBO/BOPAlgo.
//!
//! [`BopBuilder`] is the *General Fuse* algorithm — the base algorithm of the
//! Boolean Component. It consumes a [`PaveFiller`] whose intersection phase
//! already filled the [`crate::bopds::BopdsDS`] with every participating
//! shape, its pave blocks and the section vertices/edges, then:
//!
//! 1. [`BopBuilder::perform`] runs the intersection phase (`PaveFiller::perform`)
//!    and the GF building phase;
//! 2. FillImages maps each source shape to its split pieces
//!    ([`crate::bop_hist::BopHistory`] + origins + same-domain map);
//! 3. [`BopBuilder::build_result`] assembles every argument image
//!    (`BOPAlgo_Builder::BuildResult` / `BOPAlgo_BOP::BuildResult` — no
//!    obj/tool IN/OUT filter);
//! 4. Boolean runs then call [`crate::bop_bop::build_shape`]
//!    (`BOPAlgo_BOP::BuildShape`: BuildRC / BuildSolid / open-solid BuildBOP);
//! 5. [`BopBuilder::post_treat`] fixes the vertex tolerances of the result.
//!
//! Face/solid reconstruction is `crate::bop_build_faces` and
//! `crate::bop_images_solids` (GF only). [`BoolOp2`] + [`builder_bop`] convert
//! the operation into object/tool states consumed by BuildShape, not by
//! FillImagesSolids.
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::sync::Arc;

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::algo_tools::AlgoTools;
pub(crate) use crate::bop_build_common::{is_split_to_reverse, BopBuildOps};
pub(crate) use crate::bop_build_faces::BopBuilderLike;
pub(crate) use crate::bop_hist::{is_supported_type, BopHistory};
pub(crate) use crate::bop_occt_util::BopSolidHost;
pub(crate) use crate::bopds::{BopdsDS, BopdsInterf, BopdsPaveBlock};
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::fclass2d::FaceState;
pub(crate) use crate::pave_filler::PaveFiller;
pub(crate) use crate::shape::TopoShape;
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::all_subshapes;

}


mod p01;
mod p02;
mod p03;
pub use p01::*;
pub use p02::*;
pub use p03::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
