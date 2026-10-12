//! `BRepLib_MakeEdge` (`BRepLib_MakeEdge.cxx`).
//!
//! Only the two `Init` overloads the IGES read reaches **with both vertices
//! null** are ported - the form every `BRepBuilderAPI_MakeEdge` call in
//! `IGESToBRep_TopoCurve` (`cxx:1046`, `cxx:1263`),
//! `IGESToBRep_BasicCurve`, `IGESToBRep_IGESBoundary` and
//! `ShapeBuild_Edge::MakeEdge` (`ShapeBuild_Edge.cxx:793`, `:858`) makes:
//!
//! * `Init(C, p1, p2)` (`cxx:512-517`) -> `Init(C, V1, V2, p1, p2)`
//!   (`cxx:606-796`).
//! * `Init(C2d, S, p1, p2)` (`cxx:807-816`) -> `Init(C2d, S, V1, V2, p1, p2)`
//!   (`cxx:905-1094`).
//!
//! UNPORTED, and unreachable from the IGES read: the `gp_Pnt` / `TopoDS_Vertex`
//! overloads (`cxx:412-440`, `cxx:519-539`, `cxx:819-841`) and the
//! vertex-taking `Init`s (`cxx:542-572` `Init(C, V1, V2)`, whose `Project`
//! (`cxx:1133-1206`) is also unported). `degenerated` stays false in the
//! ported path: `cxx:719-727` only sets it in the closed arm with **both**
//! vertices given (`V1.IsSame(V2)` plus a mid-parameter check), so with null
//! vertices `B.Degenerated(E, false)` (`cxx:791`) is what `Init` writes.
//!
//! The one thing `Init` writes that this module does not: the edge tolerance
//! `preci` of `B.MakeEdge(E, C, preci)` (`cxx:782`) /
//! `B.UpdateEdge(E, C, S, Location(), preci)` (`cxx:1077`). `preci` is
//! `BRepLib::Precision()` = `Precision::Confusion()` (`BRepLib.cxx:84`,
//! `:124-127`) and is the same value the endpoint vertices get here, so a
//! caller that wants the full `Init` writes it with
//! `GeometryRegistry::set_edge_tolerance`.

use std::sync::Arc;

use occt_core::elib::clib2d;
use occt_core::gp::GpPnt;
use occt_core::precision::{Precision, CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::curve::Curve;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::builder::TopoBuilder;
use crate::shape::Vertex;

/// The decisions `BRepLib_MakeEdge::Init` makes before it fills the edge:
/// the checked parameter range, the two endpoint vertices (unified when the
/// curve is closed) and the two flags the caller writes back.
pub struct EdgeBuild {
    /// `p1` / `p2` after the `ElCLib::AdjustPeriodic` adjustment
    /// (`cxx:630`) or the `p1 > p2` swap (`cxx:645-652` / `cxx:948-956`).
    /// These are the values `B.Range(E, p1, p2)` stores (`cxx:790`,
    /// `cxx:1087`), **not** the arguments.
    pub first: f64,
    pub last: f64,
    /// `myVertex1` / `myVertex2` (`cxx:778-779`, `cxx:1072-1073`). `None` is
    /// the infinite-parameter arm, where OCCT adds no vertex at all
    /// (`cxx:734-775`, `cxx:1028-1065`). A closed curve gives two `TopoDS_Vertex`
    /// handles onto **one** `TShape` (`V2 = V1`, `cxx:695` / `cxx:999`),
    /// which is what `TopExp::FirstVertex` / `LastVertex` then report as the
    /// same vertex.
    pub v1: Option<Vertex>,
    pub v2: Option<Vertex>,
    /// `closed` (`cxx:687` / `cxx:987`): `P1.Distance(P2) <= preci`. Written
    /// by `E.Closed(closed)` (`cxx:792`); the 2D `Init` has no such write.
    pub closed: bool,
    /// The 2D `Init`'s `reverse` (`cxx:932`, `:957`, `:1089-1092`): a
    /// non-periodic input with `p1 >= p2` is swapped and the resulting edge is
    /// stored `REVERSED`. Always false for the 3D `Init`.
    pub reverse: bool,
}

/// `cxx:610-616`: `while (CT = down_cast<Geom_TrimmedCurve>(C)) C = CT->BasisCurve();`
/// so the range check, `IsPeriodic()` and `Value()` all see the basis curve.
fn basis_3d(mut c: Arc<dyn Curve>) -> Arc<dyn Curve> {
    while let Some(b) = c.trim_basis() {
        c = b;
    }
    c
}

/// `cxx:917-922`: the same trim-kill for `Geom2d_TrimmedCurve`.
fn basis_2d<'a>(mut c: &'a dyn Curve2d) -> &'a dyn Curve2d {
    while let Some(b) = c.trimmed_basis() {
        c = b;
    }
    c
}

/// `BRepLib_MakeEdge::Init(C, V1=null, V2=null, p1, p2)` (`cxx:606-796`).
/// `None` is the `!IsDone()` of `ParameterOutOfRange` (`cxx:655-658`) and
/// `LineThroughIdenticPoints` (`cxx:660-663`).
pub fn init_3d(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<EdgeBuild> {
    let c = basis_3d(curve.clone());
    let cf = c.first_parameter();
    let cl = c.last_parameter();
    let mut p1 = first;
    let mut p2 = last;
    if c.is_periodic() {
        // `cxx:630`: adjust in period.
        clib2d::adjust_periodic(cf, cl, PCONFUSION, &mut p1, &mut p2);
    } else {
        // `cxx:639-652`: reorder so `p1 < p2`. With `V1` and `V2` null the
        // swap of the vertex handles has no effect.
        if p1 >= p2 {
            std::mem::swap(&mut p1, &mut p2);
        }
        // `cxx:655-658`: check range.
        if cf - p1 > PCONFUSION || p2 - cl > PCONFUSION {
            return None;
        }
        // `cxx:660-663`: check ponctuallity, `gp::Resolution()`.
        if p2 - p1 <= REAL_SMALL {
            return None;
        }
    }
    // `cxx:667-677`: the points come from `GeomAdaptor_Curve::Value`, i.e. the
    // basis curve's `Value` for a `Geom_Curve` adaptor.
    let p1_inf = Precision::is_negative_infinite(p1);
    let p2_inf = Precision::is_positive_infinite(p2);
    let q1 = if p1_inf { None } else { Some(c.d0(p1)) };
    let q2 = if p2_inf { None } else { Some(c.d0(p2)) };
    // `cxx:680`: `preci = BRepLib::Precision()`.
    let preci = CONFUSION;
    // `cxx:687`: `closed = (P1.Distance(P2) <= preci)`.
    let closed = match (q1, q2) {
        (Some(a), Some(b)) => a.distance(&b) <= preci,
        _ => false,
    };
    let b = TopoBuilder::new();
    let (v1, v2) = if closed {
        // `cxx:692-697`: both vertices null -> one vertex, `V2 = V1`.
        let v = b.make_vertex(q1.unwrap_or_else(|| c.d0(p1)), preci);
        (Some(v.clone()), Some(v))
    } else {
        // `cxx:734-775`: the not-closed arm, null-vertex branches only.
        (
            q1.map(|p| b.make_vertex(p, preci)),
            q2.map(|p| b.make_vertex(p, preci)),
        )
    };
    Some(EdgeBuild { first: p1, last: p2, v1, v2, closed, reverse: false })
}

/// `BRepLib_MakeEdge::Init(C2d, S, V1=null, V2=null, p1, p2)` (`cxx:905-1094`).
/// `None` is the `!IsDone()` of `ParameterOutOfRange` (`cxx:960-964`). Note
/// there is no punctuality check in the 2D overload.
pub fn init_2d(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    first: f64,
    last: f64,
) -> Option<EdgeBuild> {
    let c = basis_2d(pc.as_ref());
    let cf = c.first_parameter();
    let cl = c.last_parameter();
    let mut p1 = first;
    let mut p2 = last;
    let mut reverse = false;
    if c.is_periodic() {
        // `cxx:934-938`: adjust in period.
        clib2d::adjust_periodic(cf, cl, PCONFUSION, &mut p1, &mut p2);
    } else {
        // `cxx:942-956`: reorder; the swap is recorded because the edge ends
        // up stored REVERSED (`cxx:1089-1092`).
        if p1 >= p2 {
            std::mem::swap(&mut p1, &mut p2);
            reverse = true;
        }
        // `cxx:960-964`: check range.
        if cf - p1 > PCONFUSION || p2 - cl > PCONFUSION {
            return None;
        }
    }
    // `cxx:968-981`: `P = S->Value(C->Value(p))`.
    let p1_inf = Precision::is_negative_infinite(p1);
    let p2_inf = Precision::is_positive_infinite(p2);
    let on_surface = |u: f64| {
        let uv = c.d0(u);
        surf.d0(uv.x(), uv.y())
    };
    let q1 = if p1_inf { None } else { Some(on_surface(p1)) };
    let q2 = if p2_inf { None } else { Some(on_surface(p2)) };
    // `cxx:983`: `preci = BRepLib::Precision()`.
    let preci = CONFUSION;
    // `cxx:987`: `closed = (P1.Distance(P2) <= preci)`.
    let closed = match (q1, q2) {
        (Some(a), Some(b)) => a.distance(&b) <= preci,
        _ => false,
    };
    let b = TopoBuilder::new();
    let (v1, v2) = if closed {
        // `cxx:994-1000`: both vertices null -> one vertex, `V2 = V1`.
        let v = b.make_vertex(q1.unwrap_or_else(|| on_surface(p1)), preci);
        (Some(v.clone()), Some(v))
    } else {
        // `cxx:1028-1065`: the not-closed arm, null-vertex branches only.
        (
            q1.map(|p| b.make_vertex(p, preci)),
            q2.map(|p| b.make_vertex(p, preci)),
        )
    };
    Some(EdgeBuild { first: p1, last: p2, v1, v2, closed, reverse })
}

/// Attach the endpoint vertices to `edge` the way `cxx:783-791` / `cxx:1078-1087`
/// do: `B.Add(E, V1)` with `V1` FORWARD (`cxx:777` / `cxx:1070`) and then
/// `B.Add(E, V2)` with `V2` REVERSED (`cxx:778` / `cxx:1071`), each only when
/// the vertex is not null. A closed curve's two handles share one `TShape`, so
/// the edge ends up with the same vertex stored FORWARD and REVERSED - which is
/// exactly what makes `TopExp::FirstVertex` and `TopExp::LastVertex` agree
/// (`TopExp.cxx:182-207`).
pub fn add_endpoint_vertices(b: &TopoBuilder, edge: &mut crate::shape::Edge, build: &EdgeBuild) {
    use crate::abs::Orientation;
    if let Some(v1) = &build.v1 {
        b.add(&mut edge.0, &v1.0.oriented(Orientation::Forward));
    }
    if let Some(v2) = &build.v2 {
        b.add(&mut edge.0, &v2.0.oriented(Orientation::Reversed));
    }
}
