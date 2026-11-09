//! `BRepExtrema_DistanceSS` — the distance between one pair of sub-shapes.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_DistanceSS.cxx`
//! (1206 lines) and `.hxx:32-133` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(S1, S2, B1, B2, DstRef, Deflection, Flag, Algo)` | `hxx:50-64` | [`DistanceSS::new`] |
//! | `Perform(S1, S2, B1, B2)` type switch | `cxx:656-775` | [`DistanceSS::perform`] |
//! | `Perform(V, V)` | `cxx:780-802` | [`DistanceSS::perform_vertex_vertex`] |
//! | `Perform(V, E)` | `cxx:806-864` | [`DistanceSS::perform_vertex_edge`] |
//! | `Perform(V, F)` | `cxx:868-925` | [`DistanceSS::perform_vertex_face`] |
//! | `Perform(E, E)` | `cxx:929-1023` | [`DistanceSS::perform_edge_edge`] |
//! | `Perform(E, F)` | `cxx:1026-1187` | [`DistanceSS::perform_edge_face`] |
//! | `Perform(F, F)` | `cxx:1191-1206` | [`DistanceSS::perform_face_face`] |
//! | `TRI_SOLUTION` | `cxx:60-75` | [`tri_solution`] |
//! | `MIN_SOLUTION` | `cxx:79-97` | [`min_solution`] |
//! | `TRIM_INFINIT_EDGE` | `cxx:102-245` | [`trim_infinit_edge`] |
//! | `TRIM_INFINIT_FACE` | `cxx:249-484` | [`trim_infinit_face`] |
//! | `PERFORM_C0` | `cxx:488-630` | [`perform_c0`] |
//! | `isOnBoundary` | `cxx:637-653` | [`is_on_boundary`] |
//!
//! # Translation notes
//!
//! * `BRepClass_FaceClassifier` is [`crate::brep_class::FaceClassifier`]; the
//!   accepted states are `IN` only in the vertex-face (`cxx:907`),
//!   edge-face (`cxx:1080`) and face-face arms, exactly as the `.cxx` reads.
//! * `Extrema_ExtPC` is [`occt_geom::extrema_pc::extrema_ext_pc_range`] plus
//!   its `IsDone()`-carrying twin; the wrapper class lives in
//!   [`super::brep_extrema_ext_pc`].
//! * `GeomAdaptor_Curve::NbIntervals(GeomAbs_C1)` / `Intervals(...)` are
//!   [`occt_geom::extrema_pc::extrema_nb_intervals`] /
//!   [`occt_geom::extrema_pc::extrema_intervals`] (`s = 1` is `GeomAbs_C1`).
//! * `BRep_Tool::Curve(E, F, L)` is [`BRepTool::edge_curve`] +
//!   [`BRepTool::edge_parameters`].
//! * `BRepBuilderAPI_MakeVertex(P)` is
//!   [`crate::builder::TopoBuilder::make_vertex`] with
//!   `Precision::Confusion()` (the constructor's default tolerance).
//! * `TRIM_INFINIT_EDGE` rebuilds an edge with `BRepBuilderAPI_MakeEdge(new
//!   Geom_TrimmedCurve(pCurv, Umin, Umax))` (`cxx:240-241`); the port uses
//!   [`crate::builder::TopoBuilder::make_edge`] over a
//!   [`occt_geom::trimmed::GeomTrimmedCurve`].
//! * `TRIM_INFINIT_FACE` rebuilds a face with `BRepBuilderAPI_MakeFace(new
//!   Geom_RectangularTrimmedSurface(...), Precision::Confusion())`
//!   (`cxx:471-472`); the port uses
//!   [`crate::builder::TopoBuilder::make_face_with_tolerance`].
//! * The `NCollection_Sequence` pairs stay two parallel `Vec`s, like OCCT's
//!   two sequences.

/// Crate imports shared by `mod.rs` and every submodule.
mod prelude {
    pub(super) use std::sync::Arc;
    pub(super) use occt_core::bnd::BndBox;
    pub(super) use occt_core::gp::{GpPnt, GpPnt2d};
    pub(super) use occt_core::precision::Precision;
    pub(super) use occt_geom::curve::Curve;
    pub(super) use occt_geom::extrema_pc;
    pub(super) use crate::brep_class::FaceClassifier;
    pub(super) use crate::fclass2d::FaceState;
    pub(super) use crate::brep_bnd_lib;
    pub(super) use crate::brep_tool::BRepTool;
    pub(super) use crate::brep_extrema_ext_cc::ExtCC;
    pub(super) use crate::brep_extrema_ext_cf::ExtCF;
    pub(super) use crate::brep_extrema_ext_ff::ExtFF;
    pub(super) use crate::brep_extrema_ext_pc::ExtPC;
    pub(super) use crate::brep_extrema_ext_pf::ExtPF;
    pub(super) use crate::brep_extrema_solution_elem::{SolutionElem, SupportType};
    pub(super) use crate::builder::TopoBuilder;
    pub(super) use crate::abs::ShapeType;
    pub(super) use crate::shape::{Edge, Face, TopoShape, Vertex};
    pub(super) use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};
}

use self::prelude::*;

/// `BRepExtrema_DistanceSS` (`BRepExtrema_DistanceSS.hxx:32-133`).
pub struct DistanceSS {
    /// `myDstRef` (`hxx:126`).
    my_dst_ref: f64,
    /// `myModif` (`hxx:127`) — `IsDone()`.
    my_modif: bool,
    /// `myEps` (`hxx:128`).
    my_eps: f64,
    /// `mySeqSolShape1` (`hxx:130`).
    my_seq_sol_shape1: Vec<SolutionElem>,
    /// `mySeqSolShape2` (`hxx:131`).
    my_seq_sol_shape2: Vec<SolutionElem>,
    /// `myFlag` / `myAlgo` (`hxx:132-133`).
    my_flag: ExtremaExtFlag,
    my_algo: ExtremaExtAlgo,
}

impl DistanceSS {
    /// `BRepExtrema_DistanceSS(S1, S2, Box1, Box2, DstRef, Deflection, Flag,
    /// Algo)` (`hxx:50-64`): stores the reference distance and runs `Perform`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        the_s1: &TopoShape,
        the_s2: &TopoShape,
        the_box1: &BndBox,
        the_box2: &BndBox,
        the_dst_ref: f64,
        the_deflection: f64,
        the_ext_flag: ExtremaExtFlag,
        the_ext_algo: ExtremaExtAlgo,
    ) -> Self {
        let mut this = Self {
            my_dst_ref: the_dst_ref,
            my_modif: false,
            my_eps: the_deflection,
            my_seq_sol_shape1: Vec::new(),
            my_seq_sol_shape2: Vec::new(),
            my_flag: the_ext_flag,
            my_algo: the_ext_algo,
        };
        this.perform(the_s1, the_s2, the_box1, the_box2);
        this
    }
}

// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod perform;
mod perform_c0;
mod perform_pairs;
mod statics;
mod trim_infinite;
