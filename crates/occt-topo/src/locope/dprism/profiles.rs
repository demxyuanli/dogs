//! The profile wire of `LocOpe_DPrism` (`cxx:70-110`, `cxx:366-390`).
//!
//! Both constructors build three segments in the frame `gp::XOY()` the sweep
//! runs in: `myProfile2` is the swept wall, `myProfile3` the "first" extension
//! and `myProfile1` the "last" one. The end-cap heights are turned into the
//! tilted points of the reference (`y = H*sin(Angle)`, `z = H*cos(Angle)`), and
//! each extension is long enough to leave the spine's UV box, twice the larger
//! side plus its own offset (`cxx:90-92`, `cxx:102-104`).
//!
//! `myProfile` is `BRepLib_MakeWire(myProfile1, myProfile2, myProfile3)`
//! (`cxx:108`, `cxx:388`), i.e. three `Add` calls of `BRepLib_MakeWire`
//! (`cxx:110-113`) - not the plain `BRep_Builder` append, because the class
//! reverses and re-connects the edges by proximity.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, REAL_SMALL};

use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Vertex, Wire};

/// `myProfile1` / `myProfile2` / `myProfile3` and the wire built from them.
pub(crate) struct Profile {
    /// `myProfile1` (`hxx:75`).
    pub edge1: Edge,
    /// `myProfile2` (`hxx:76`).
    pub edge2: Edge,
    /// `myProfile3` (`hxx:77`).
    pub edge3: Edge,
    /// `myProfile` (`hxx:74`).
    pub wire: Wire,
}

/// `LocOpe_DPrism(Spine, Height1, Height2, Angle)` profile (`cxx:70-109`).
///
/// `None` is the `!IsDone()` arm of one of the `BRepLib_MakeEdge` /
/// `BRepLib_MakeWire` calls (identical points, `cxx:189-193`, or a wire that
/// cannot be connected).
pub(crate) fn dual_profile(
    spine: &Face,
    height1: f64,
    height2: f64,
    angle: f64,
) -> Option<Profile> {
    let b = TopoBuilder::new();

    // `cxx:72-73`.
    let y = height1 * angle.sin();
    let z = height1 * angle.cos();
    // `cxx:76`.
    let vert2 = make_vertex(&b, GpPnt::new(0.0, y, z));

    // `cxx:78-79`.
    let y1 = -height2 * angle.sin();
    let z1 = -height2 * angle.cos();
    // `cxx:81`.
    let vert1 = make_vertex(&b, GpPnt::new(0.0, y1, z1));

    // `cxx:82`.
    let edge2 = segment(&b, &vert1, &vert2)?;

    // `cxx:84-92`.
    let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(spine);
    let mut deltay = (umax - umin).max(vmax - vmin) + y.abs();
    deltay *= 2.0;

    // `cxx:94`.
    let vert3 = make_vertex(&b, GpPnt::new(0.0, y + deltay, z));
    // `cxx:95`.
    let edge3 = segment(&b, &vert2, &vert3)?;

    // `cxx:97-104`: the reference asks for the same UV box a second time.
    let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(spine);
    let mut deltay1 = (umax - umin).max(vmax - vmin) + y1.abs();
    deltay1 *= 2.0;

    // `cxx:106`.
    let vert4 = make_vertex(&b, GpPnt::new(0.0, y1 + deltay1, z1));
    // `cxx:107`.
    let edge1 = segment(&b, &vert4, &vert1)?;

    // `cxx:109`.
    let wire = make_wire(&[&edge1, &edge2, &edge3])?;
    Some(Profile {
        edge1,
        edge2,
        edge3,
        wire,
    })
}

/// `LocOpe_DPrism(Spine, Height, Angle)` profile (`cxx:366-388`).
pub(crate) fn single_profile(spine: &Face, height: f64, angle: f64) -> Option<Profile> {
    let b = TopoBuilder::new();

    // `cxx:368-369`.
    let y = height * angle.sin();
    let z = height * angle.cos();

    // `cxx:371-373`.
    let vert1 = make_vertex(&b, GpPnt::new(0.0, 0.0, 0.0));
    let vert2 = make_vertex(&b, GpPnt::new(0.0, y, z));
    // `cxx:374`.
    let edge2 = segment(&b, &vert1, &vert2)?;

    // `cxx:376-380`.
    let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(spine);
    let mut deltay = (umax - umin).max(vmax - vmin) + y.abs();
    deltay *= 2.0;

    // `cxx:382-384`.
    let vert3 = make_vertex(&b, GpPnt::new(0.0, y + deltay, z));
    let edge3 = segment(&b, &vert2, &vert3)?;

    // `cxx:385-387`.
    let vert4 = make_vertex(&b, GpPnt::new(0.0, deltay, 0.0));
    let edge1 = segment(&b, &vert4, &vert1)?;

    // `cxx:388`.
    let wire = make_wire(&[&edge1, &edge2, &edge3])?;
    Some(Profile {
        edge1,
        edge2,
        edge3,
        wire,
    })
}

/// `BRepLib_MakeVertex(const gp_Pnt&)` (`BRepLib_MakeVertex.cxx:30-40`): a
/// vertex carrying `BRepLib::Precision()` (`cxx:76` `:81` `:93` `:106`).
fn make_vertex(b: &TopoBuilder, p: GpPnt) -> Vertex {
    b.make_vertex(p, CONFUSION)
}

/// `BRepLib_MakeEdge(const TopoDS_Vertex&, const TopoDS_Vertex&)`
/// (`BRepLib_MakeEdge.cxx:185-198`): the segment `V1 -> V2` over `[0, l]`.
/// `None` is `BRepLib_LineThroughIdenticPoints` (`cxx:189-193`).
fn segment(b: &TopoBuilder, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    if p1.distance(&p2) <= REAL_SMALL {
        return None;
    }
    Some(b.make_edge_segment_with_vertices(&p1, &p2, v1, v2))
}

/// `BRepLib_MakeWire(E1, E2, E3)` (`BRepLib_MakeWire.cxx:66-78`): three `Add`
/// calls; any `myError` other than `WireDone` leaves `IsDone()` false.
fn make_wire(edges: &[&Edge]) -> Option<Wire> {
    let mut mw = MakeWire::new();
    for e in edges {
        mw.add(e).ok()?;
    }
    if !mw.is_done() {
        return None;
    }
    Some(mw.wire())
}
