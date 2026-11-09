//! Bi-point assembly of the outline update. Source: `HLRBRep_PolyAlgo.cxx:878-2978`
//! (`InitBiPointsWithConnexity`, the three `Interpolation` overloads and the
//! three `MoveOrInsertPoint` overloads).
//!
//! The C++ overloads of `Interpolation` / `MoveOrInsertPoint` differ only in how
//! many faces bound the edge run being split: one face (an edge with a single
//! adjacent face), two faces (an edge between two faces, one interpolated point
//! on each), and two faces with *two* interpolated points. The port keeps the
//! same three levels: [`HlrBRepPolyAlgo::interpolation_one`] /
//! [`HlrBRepPolyAlgo::move_or_insert_one`], `interpolation_two` /
//! `move_or_insert_two` and `move_or_insert_two_two`.
//!
//! The long in/out parameter lists (`theX1 .. theU2`) become [`Segment`], the
//! `(theNod11*, theNod12*, theI1p1, theI1p2, theI1, thePid1)` bundle of one side
//! becomes [`FaceSide`], and the `(X3..ZTI3, Coef3, U3, insP3, mP3P1)` out bundle
//! of one interpolation becomes [`Interp`]. OCCT hands out `NodeIndices&` /
//! `NodeData&` references into a reallocating array; the port passes the node
//! index instead and re-reads the slot, which is exactly what OCCT's post-
//! `AddNode` re-fetch (`cxx:2039-2041`, `:2349-2354`, `:2819-2830`) does.
//!
//! UNPORTED / deviations (each with its OCCT line):
//! * `BRep_Tool::PolygonOnTriangulation` needs a `BRep_PolygonOnTriangulation`
//!   representation on the edge (`cxx:900`, `:1070`, `:1079`). The producer is
//!   `BRepMesh_ModelPostProcessor::PolygonCommitter`
//!   (`BRepMesh_ModelPostProcessor.cxx:42-157`), ported as
//!   [`occt_topo::meshing::model_post_processor::PolygonCommitter`] and driven by
//!   `IncrementalMesh::perform`. An edge that still has no representation (for
//!   instance when the mesh skipped its wire) takes the same null-handle path
//!   OCCT takes (the `#ifdef OCCT_DEBUG` branch at `cxx:1055-1063`,
//!   `:1301-1309`).
//! * A polygon without parameters (`!HasParameters()`): OCCT dereferences the
//!   null `par` handle at `cxx:908` / `:928` (undefined behaviour); the port
//!   skips the pairing for that edge instead of reproducing a crash.
//! * `BRep_Tool::Continuity(theEdge, aF1, aF2)` (`cxx:1080`) has no direct
//!   counterpart on a `BRep_TEdge`: the port's continuity lives in
//!   [`occt_topo::brep_lib_continuity::stored_continuity`] (what
//!   `BRepLib::EncodeRegularity` records, `BRepLib.cxx:279-330`), which answers
//!   `GeomAbs_C0` when no record matches - the same fallback as
//!   `BRep_TEdge::Continuity`. The decoded `GeomAbs_Shape` reaches this module
//!   as `the_rg` (`GeomAbs_G1` -> 1, `GeomAbs_G2` -> 3, `GeomAbs_Shape.hxx:25-29`)
//!   and turns into the flag of `cxx:1512-1519`.
//! * `myBCurv.Value(theU3)` (`cxx:1874`) needs the edge's 3D curve. When the
//!   edge has none, OCCT's `BRepAdaptor_Curve` has no curve to evaluate (the
//!   `Initialize` at `cxx:891` already left it null); the port answers "no
//!   interpolated point", i.e. the same `false` the enclosing `NewNode` failure
//!   arm produces.
//! * `TopoDS_Edge::Reverse()` (`cxx:1076`) only flips the local copy's
//!   orientation: the port threads the flag into the second polygon lookup
//!   instead of mutating a shape.
//! * `HLRAlgo_PolyInternalData::IncPINod` / `IncTData` / `IncPISeg` reallocate
//!   *two* arrays at once in OCCT (`PolyInternalData.cxx:923-1038`); the port's
//!   arrays belong to one face object, so the second out-parameter disappears
//!   (same deviation T-408 recorded).

use std::sync::Arc;

use occt_core::gp::GpXyz;
use occt_core::poly::Triangulation;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use super::{nmask, HlrBRepPolyAlgo, PolyInternalDataArray};
use crate::hlr_algo::{
    HlrAlgoBiPoint, HlrAlgoPolyInternalData, PolyInternalNodeData, PolyInternalNodeIndices,
};
use occt_topo::shape::{Edge, TopoShape};

/// `GeomAbs_G1` (`GeomAbs_Shape.hxx:26`).
pub(super) const GEOMABS_G1: u8 = 1;
/// `GeomAbs_G2` (`GeomAbs_Shape.hxx:28`).
pub(super) const GEOMABS_G2: u8 = 3;

/// The `theX1 .. theU2` block of `Interpolation` / `MoveOrInsertPoint`.
///
/// `x1/y1/z1` and `x2/y2/z2` are the **projected** endpoints (node space), the
/// `xt*/yt*/zt*` fields the 3D ones; `u1`/`u2` are the edge parameters. Every
/// field is in/out: `MoveOrInsertPoint` overwrites an endpoint when the
/// interpolated point snaps onto it (`cxx:1951-1959`, `:2001-2009`).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Segment {
    pub x1: f64,
    pub y1: f64,
    pub z1: f64,
    pub x2: f64,
    pub y2: f64,
    pub z2: f64,
    pub xt1: f64,
    pub yt1: f64,
    pub zt1: f64,
    pub xt2: f64,
    pub yt2: f64,
    pub zt2: f64,
    pub u1: f64,
    pub u2: f64,
}

/// One side of a bi-point pair: `(theNod1*, theNod2*, theI1p1, theI1p2, theI1,
/// thePid1)`. `face` is `theI1` (1-based) and doubles as the
/// [`PolyInternalDataArray`] slot, because `thePid1` is always `thePID(theI1)`
/// (`cxx:902`, `:1081`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct FaceSide {
    pub face: usize,
    pub p1: i32,
    pub p2: i32,
}

/// The `(X3, Y3, Z3, XTI3, YTI3, ZTI3, Coef3, U3, insP3, mP3P1)` out bundle of
/// `Interpolation(theU1, theU2, ...)` (`cxx:1857-1884`).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Interp {
    /// Projected point (`TMultiply` of the 3D point, `cxx:1878`).
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// 3D point (`myBCurv.Value(theU3)`, `cxx:1874-1877`).
    pub xt: f64,
    pub yt: f64,
    pub zt: f64,
    pub coef: f64,
    pub u: f64,
    pub ins: bool,
    /// `mP3P1`: the crossing sits nearer `P1` than `P2`.
    pub m_p3_p1: bool,
}

/// `BRep_Tool::Triangulation(face, loc)` for a bare [`TopoShape`].
fn face_triangulation(face: &TopoShape) -> Option<Arc<Triangulation>> {
    face.tshape.read().unwrap().face_triangulation().cloned()
}

/// `BRepAdaptor_Curve2d(theEdge, theFace)` (`BRepAdaptor_Curve2d.cxx:39-53`):
/// the pcurve `BRep_Tool::CurveOnSurface` answers for that face, with the
/// `IsPeriodic()` / `Period()` the `U2` wrap at `cxx:961-964` needs.
///
/// `reversed` selects the second (reversed) pcurve of a seam edge, the way
/// `BRep_Tool::CurveOnSurface` does when the edge orientation and the face
/// orientation disagree.
fn pcurve_of(edge: &TopoShape, face: &TopoShape, reversed: bool) -> Option<Arc<dyn Curve2d>> {
    let face_key = occt_topo::tgeometry::GeometryRegistry::shape_key(face);
    let pcs = occt_topo::tgeometry::GeometryRegistry::global().edge_pcurves(edge, face_key);
    if pcs.is_empty() {
        return None;
    }
    let pos = if reversed && pcs.len() >= 2 { 1 } else { 0 };
    Some(pcs[pos].clone())
}

fn node_data(data: &HlrAlgoPolyInternalData, i: i32) -> PolyInternalNodeData {
    *data.pi_nod()[i as usize].data()
}

fn node_indices(data: &HlrAlgoPolyInternalData, i: i32) -> PolyInternalNodeIndices {
    *data.pi_nod()[i as usize].indices()
}

/// The `theNod1.RValues.Point = gp_XYZ(theX3, ...)`, `PCu1`/`PCu2`, `Scal = 0`
/// and `Flag |= NMsk_OutL` block that `MoveOrInsertPoint` repeats for every face
/// side (`cxx:1958-1977`, `:2008-2027`, and their two-face twins).
fn move_node_onto_point(
    data: &mut HlrAlgoPolyInternalData,
    node: i32,
    the_i_edge: i32,
    x3: f64,
    y3: f64,
    z3: f64,
    u3: f64,
) {
    let (edg1, edg2) = {
        let ind = node_indices(data, node);
        (ind.edg1, ind.edg2)
    };
    let nd = data.change_pi_nod()[node as usize].change_data();
    nd.point = GpXyz::new(x3, y3, z3);
    // `cxx:1959-1973`: the parameter lands on whichever slot holds this edge.
    if edg1 == the_i_edge {
        nd.pcu1 = u3;
    } else if edg2 == the_i_edge {
        nd.pcu2 = u3;
    }
    nd.scal = 0.0;
    data.change_pi_nod()[node as usize].change_indices().flag |= nmask::OUT_L;
}

mod connexity;
mod one_face;
mod two_faces;

