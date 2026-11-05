//! `BRepLib_MakeEdge` (`BRepLib_MakeEdge.cxx`).
//!
//! Ported:
//!
//! * `Init(C, p1, p2)` (`cxx:505-517`) -> `Init(C, V1, V2, p1, p2)`
//!   (`cxx:603-796`) with **both vertices null** - the form every
//!   `BRepBuilderAPI_MakeEdge` call in `IGESToBRep_TopoCurve` (`cxx:1046`,
//!   `cxx:1263`), `IGESToBRep_BasicCurve`, `IGESToBRep_IGESBoundary` and
//!   `ShapeBuild_Edge::MakeEdge` (`ShapeBuild_Edge.cxx:793`, `:858`) makes.
//! * `Init(C2d, S, p1, p2)` (`cxx:808-819`) -> `Init(C2d, S, V1, V2, p1, p2)`
//!   (`cxx:907-1094`), including the vertex-taking arms of `cxx:992-1067`
//!   (`BRepAlgo::ConvertWire` is the caller that passes both vertices,
//!   `BRepAlgo.cxx:190-205`).
//! * `Init(C, V1, V2)` (`cxx:543-572`), reached through
//!   [`init_3d_projecting`], and the vertex-taking arm of
//!   `Init(C, V1, V2, p1, p2)` (`cxx:603-796`) it feeds; the
//!   `Project(C, V, p)` (`cxx:55-109`) it needs is [`project`]. This is the
//!   form `TopOpeBRepTool_FuseEdges::BuildListResultEdges`
//!   (`TopOpeBRepTool_FuseEdges.cxx:334`) and
//!   `TopOpeBRepBuild_FuseFace::MakeEdge` (`FuseFace.cxx:1163`) use.
//!
//! `degenerated` is only set in the closed arm with **both** vertices given
//! (`cxx:719-727`: `V1.IsSame(V2)` plus a mid-parameter check), so with null
//! vertices `B.Degenerated(E, false)` (`cxx:791`) is what `Init` writes.
//!
//! * `BRepLib_MakeEdge(V1, V2)` (`cxx:185-198`), the vertex-to-vertex straight
//!   segment -> [`init_line_vertices`]. It is the `BRepLib_MakeEdge` used by
//!   `BRepFill_TrimShellCorner::MakeFacesNonSec` (`BRepFill_TrimShellCorner.cxx:528`,
//!   `:542`) to close the corner between two edges.
//!
//! UNPORTED: the `Init(C, P1, P2, p1, p2)` overload (`cxx:575-595`) and the 2D
//! vertex-taking pair (`cxx:821-843` `Init(C2d, S, V1, V2)` with its
//! `Project(C2d, S, V, p)` at `cxx:115-176`, `cxx:845-876`
//! `Init(C2d, S, P1, P2)`), plus the `gp_Lin`/`gp_Circ`/... `BRepLib_MakeEdge`
//! constructors (`cxx:217-386`). The 3D `Init(C, P1, P2)` (`cxx:522-541`) is
//! [`init_curve_points`].
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
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::{Precision, CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Vertex;
use crate::topo_tools::is_same;

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
    /// `degenerated` (`cxx:722`): only the closed arm with **two** given
    /// vertices can set it, when the mid-parameter point sits on `P1`.
    pub degenerated: bool,
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

/// `GeomAdaptor_Curve` over `C` (`GeomAdaptor_Curve.cxx:252-254`): the untrimmed
/// basis plus the interval the trim leaves, expressed in that basis's
/// parameters. OCCT's `Project` builds the adaptor from the handle it is given,
/// so both the end parameters and the `Extrema_ExtPC` search run on the basis
/// restricted to the trim interval — the same space the full `Init`
/// (`cxx:610-616`, which kills the trimmed curves first) expects `p1`/`p2` in.
///
/// The port's [`GeomTrimmedCurve`](occt_geom::trimmed::GeomTrimmedCurve)
/// reparameterises its trim onto `[0, 1]`, so `first_parameter()` /
/// `last_parameter()` of a trimmed handle are **not** the adaptor's range;
/// `untrimmed_basis()` is the view that reports it in basis parameters.
pub(crate) fn adaptor_curve(curve: &Arc<dyn Curve>) -> (Arc<dyn Curve>, f64, f64) {
    curve
        .untrimmed_basis()
        .unwrap_or_else(|| (curve.clone(), curve.first_parameter(), curve.last_parameter()))
}

/// `Project(C, V, p)` (`cxx:55-108`): the parameter of `V` on `C`, accepted
/// only when the squared distance stays within `V`'s squared tolerance.
///
/// The two interval ends are tested first (`cxx:70-79`), then `Extrema_ExtPC`
/// picks the closest solution (`cxx:83-105`). `None` is the C++ `false`. The
/// returned parameter is in the basis curve of `C`, as OCCT's adaptor is.
pub fn project(curve: &Arc<dyn Curve>, v: &Vertex) -> Option<f64> {
    let eps2 = {
        let t = crate::brep_tool::BRepTool::vertex_tolerance(v);
        t * t
    };
    let p = crate::brep_tool::BRepTool::vertex_point_world(v);
    // `cxx:62`: `GeomAdaptor_Curve GAC(C)`.
    let (c, first, last) = adaptor_curve(curve);
    // `cxx:66-79`: the ends win when they are the closer one and inside `Eps2`.
    let d1 = if Precision::is_infinite(first) {
        f64::MAX
    } else {
        c.d0(first).square_distance(&p)
    };
    let d2 = if Precision::is_infinite(last) {
        f64::MAX
    } else {
        c.d0(last).square_distance(&p)
    };
    if d1 < d2 && d1 <= eps2 {
        return Some(first);
    } else if d2 < d1 && d2 <= eps2 {
        return Some(last);
    }
    // `cxx:83-105`
    let (sols, done) = occt_geom::extrema_pc::extrema_ext_pc_range_done(c.as_ref(), &p, first, last);
    if !done {
        return None;
    }
    let mut best: Option<&occt_geom::extrema_pc::ExtPcSolution> = None;
    for s in &sols {
        if best.map_or(true, |b| s.sq_dist < b.sq_dist) {
            best = Some(s);
        }
    }
    let best = best?;
    if best.sq_dist <= eps2 {
        Some(best.u)
    } else {
        None
    }
}

/// `BRepLib_MakeEdge::Init(C, V1, V2)` (`cxx:542-572`): `Project` each given
/// vertex, then the full `Init`. `None` is `BRepLib_PointProjectionFailed`.
pub fn init_3d_projecting(
    curve: &Arc<dyn Curve>,
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
) -> Option<EdgeBuild> {
    // `cxx:552` / `:562`: `C->FirstParameter()` / `C->LastParameter()` of the
    // handle as given, i.e. the trim ends in basis parameters.
    let (_, first, last) = adaptor_curve(curve);
    let p1 = match v1 {
        Some(v) => project(curve, v)?,
        None => first,
    };
    let p2 = match v2 {
        Some(v) => project(curve, v)?,
        None => last,
    };
    init_3d_vertices(curve, v1, v2, p1, p2)
}

/// `BRepLib_MakeEdge::Init(C, V1, V2, p1, p2)` (`cxx:606-796`): the whole
/// decision table, including the vertex arms the null-vertex
/// [`init_3d`] skips.
///
/// `None` is any `myError` arm (`cxx:653-657` `ParameterOutOfRange`,
/// `:660-663` `LineThroughIdenticPoints`, `:709-718` `DifferentPointsOnClosedCurve`,
/// `:735-739` / `:754-758` `PointWithInfiniteParameter`,
/// `:747-751` / `:767-771` `DifferentsPointAndParameter`): OCCT leaves
/// `IsDone()` false and the caller stops.
pub fn init_3d_vertices(
    curve: &Arc<dyn Curve>,
    vv1: Option<&Vertex>,
    vv2: Option<&Vertex>,
    first: f64,
    last: f64,
) -> Option<EdgeBuild> {
    // `cxx:610-616`
    let c = basis_3d(curve.clone());
    let cf = c.first_parameter();
    let cl = c.last_parameter();
    let epsilon = PCONFUSION;
    let mut p1 = first;
    let mut p2 = last;
    let (mut v1, mut v2) = (vv1.cloned(), vv2.cloned());
    if c.is_periodic() {
        // `cxx:630`: adjust in period; the vertex handles keep their roles.
        clib2d::adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
    } else {
        // `cxx:638-650`: reorder, swapping the vertex handles too.
        if !(p1 < p2) {
            std::mem::swap(&mut v2, &mut v1);
            std::mem::swap(&mut p1, &mut p2);
        }
        // `cxx:653-657`: check range.
        if cf - p1 > epsilon || p2 - cl > epsilon {
            return None;
        }
        // `cxx:660-663`: check ponctuallity, `gp::Resolution()`.
        if p2 - p1 <= REAL_SMALL {
            return None;
        }
    }

    // `cxx:667-677`
    let p1_inf = Precision::is_negative_infinite(p1);
    let p2_inf = Precision::is_positive_infinite(p2);
    let q1 = if p1_inf { None } else { Some(c.d0(p1)) };
    let q2 = if p2_inf { None } else { Some(c.d0(p2)) };

    // `cxx:680-733`
    let preci = CONFUSION;
    let mut degenerated = false;
    let closed = match (&q1, &q2) {
        (Some(a), Some(b)) => a.distance(b) <= preci,
        _ => false,
    };
    if closed {
        let (qa, qb) = (q1.unwrap(), q2.unwrap());
        if v1.is_none() && v2.is_none() {
            let v = TopoBuilder::new().make_vertex(qa, preci);
            v1 = Some(v.clone());
            v2 = Some(v);
        } else if v1.is_none() {
            v1 = v2.clone();
        } else if v2.is_none() {
            v2 = v1.clone();
        } else {
            let (va, vb) = (v1.clone().unwrap(), v2.clone().unwrap());
            if !crate::topo_tools::is_same(&va.0, &vb.0) {
                return None;
            }
            let pv = crate::brep_tool::BRepTool::vertex_point_world(&va);
            if qa.distance(&pv)
                > preci.max(crate::brep_tool::BRepTool::vertex_tolerance(&va))
            {
                return None;
            }
            // `cxx:722-725`: a mid-parameter point on `P1` means the edge is
            // degenerated, not a full turn.
            if let Some(pm) = q1.map(|_| c.d0((p1 + p2) / 2.0)) {
                if qa.distance(&pm) < preci {
                    degenerated = true;
                }
            }
        }
    } else {
        // `cxx:732-775`: not closed.
        if p1_inf {
            if v1.is_some() {
                return None;
            }
        } else if v1.is_none() {
            v1 = q1.map(|p| TopoBuilder::new().make_vertex(p, preci));
        } else {
            let va = v1.clone().unwrap();
            let pv = crate::brep_tool::BRepTool::vertex_point_world(&va);
            if q1.map_or(false, |p| {
                p.distance(&pv) > preci.max(crate::brep_tool::BRepTool::vertex_tolerance(&va))
            }) {
                return None;
            }
        }
        if p2_inf {
            if v2.is_some() {
                return None;
            }
        } else if v2.is_none() {
            v2 = q2.map(|p| TopoBuilder::new().make_vertex(p, preci));
        } else {
            let vb = v2.clone().unwrap();
            let pv = crate::brep_tool::BRepTool::vertex_point_world(&vb);
            if q2.map_or(false, |p| {
                p.distance(&pv) > preci.max(crate::brep_tool::BRepTool::vertex_tolerance(&vb))
            }) {
                return None;
            }
        }
    }

    // `cxx:777-778`: `V1` FORWARD, `V2` REVERSED. The orientation is applied
    // by `add_endpoint_vertices`, which is where OCCT's `B.Add` reads it.
    Some(EdgeBuild { first: p1, last: p2, v1, v2, closed, reverse: false, degenerated })
}

/// `BRepLib_MakeEdge::Init(C, V1=null, V2=null, p1, p2)` (`cxx:606-796`).
/// `None` is the `!IsDone()` of `ParameterOutOfRange` (`cxx:655-658`) and
/// `LineThroughIdenticPoints` (`cxx:660-663`).
pub fn init_3d(curve: &Arc<dyn Curve>, first: f64, last: f64) -> Option<EdgeBuild> {
    init_3d_vertices(curve, None, None, first, last)
}

/// `BRepLib_MakeEdge::BRepLib_MakeEdge(V1, V2)` (`cxx:185-198`): the straight
/// segment between the two vertex points. `P1 == P2` within `gp::Resolution()`
/// gives `BRepLib_LineThroughIdenticPoints` (`cxx:189-193`) -> `None`; otherwise
/// a `Geom_Line` through `P1` with direction `P1 -> P2` is run through
/// `Init(GL, V1, V2, 0, l)` (`cxx:196`).
///
/// Both vertices keep their orientation role of the `Init` contract: `V1` is the
/// FORWARD end at parameter `0`, `V2` the REVERSED end at `l`.
pub fn init_line_vertices(p1: &GpPnt, p2: &GpPnt, v1: &Vertex, v2: &Vertex) -> Option<EdgeBuild> {
    let l = p1.distance(p2);
    // `cxx:189`: `gp::Resolution()`.
    if l <= REAL_SMALL {
        return None;
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).ok()?;
    let line = GeomLine::new(GpLin::from_pnt_dir(*p1, dir));
    init_3d_vertices(&(Arc::new(line) as Arc<dyn Curve>), Some(v1), Some(v2), 0.0, l)
}

/// `BRepLib_MakeEdge::Init(C, P1, P2)` (`cxx:522-541`): the two vertices are
/// built at `BRepLib::Precision()` (`cxx:526`), `V2` aliases `V1` when the
/// points are closer than that tolerance (`cxx:529-533`), and the pair then
/// goes through `Init(C, V1, V2)` (`cxx:538`), i.e. [`init_3d_projecting`].
///
/// `None` is `BRepLib_PointProjectionFailed` (`cxx:555` / `:566`), raised when
/// [`project`] cannot place a vertex on the curve.
pub fn init_curve_points(
    curve: &Arc<dyn Curve>,
    p1: &GpPnt,
    p2: &GpPnt,
) -> Option<EdgeBuild> {
    let b = TopoBuilder::new();
    let preci = CONFUSION;
    let v1 = b.make_vertex(*p1, preci);
    let v2 = if p1.distance(p2) < preci {
        v1.clone()
    } else {
        b.make_vertex(*p2, preci)
    };
    init_3d_projecting(curve, Some(&v1), Some(&v2))
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
    init_2d_vertices(pc, surf, None, None, first, last)
}

/// `BRepLib_MakeEdge::Init(C2d, S, VV1, VV2, p1, p2)` (`cxx:905-1094`), the
/// vertex-taking overload `BRepBuilderAPI_MakeEdge(C, S, V1, V2, p1, p2)`
/// forwards to (`BRepBuilderAPI_MakeEdge.cxx:483-496`).
///
/// `None` is `!IsDone()`; besides `ParameterOutOfRange` (`cxx:960-964`) the
/// vertex arms raise `DifferentPointsOnClosedCurve` (`cxx:1013`, `:1018`),
/// `PointWithInfiniteParameter` (`cxx:1032`, `:1052`) and
/// `DifferentsPointAndParameter` (`cxx:1042`, `:1063`).
pub fn init_2d_vertices(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    vv1: Option<&Vertex>,
    vv2: Option<&Vertex>,
    first: f64,
    last: f64,
) -> Option<EdgeBuild> {
    let c = basis_2d(pc.as_ref());
    let cf = c.first_parameter();
    let cl = c.last_parameter();
    let mut p1 = first;
    let mut p2 = last;
    let mut reverse = false;
    // `TopoDS_Vertex V1, V2` (`cxx:931`), filled by the two arms below.
    let mut v1 = vv1.cloned();
    let mut v2 = vv2.cloned();
    if c.is_periodic() {
        // `cxx:934-938`: adjust in period; the vertices keep their order.
        clib2d::adjust_periodic(cf, cl, PCONFUSION, &mut p1, &mut p2);
    } else {
        // `cxx:942-956`: reorder; the swap is recorded because the edge ends
        // up stored REVERSED (`cxx:1089-1092`).
        if p1 >= p2 {
            // `V2 = VV1; V1 = VV2` (`cxx:951-952`).
            std::mem::swap(&mut v1, &mut v2);
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
    let closed = match (&q1, &q2) {
        (Some(a), Some(b)) => a.distance(b) <= preci,
        _ => false,
    };
    let b = TopoBuilder::new();
    if closed {
        // `cxx:994-1021`.
        match (v1.is_none(), v2.is_none()) {
            (true, true) => {
                // `cxx:996-999`: both null -> one vertex, `V2 = V1`.
                let v = b.make_vertex(q1.unwrap_or_else(|| on_surface(p1)), preci);
                v1 = Some(v.clone());
                v2 = Some(v);
            }
            (true, false) => {
                // `cxx:1001-1004`: `V1 = V2`.
                v1 = v2.clone();
            }
            (false, true) => {
                // `cxx:1005-1008`: `V2 = V1`.
                v2 = v1.clone();
            }
            (false, false) => {
                // `cxx:1010-1021`.
                let (a, bb) = (v1.as_ref()?, v2.as_ref()?);
                if !is_same(&a.0, &bb.0) {
                    // `BRepLib_DifferentPointsOnClosedCurve`
                    return None;
                }
                let dist = q1.as_ref()?.distance(&BRepTool::vertex_point_world(a));
                if dist > preci.max(BRepTool::vertex_tolerance(a)) {
                    return None;
                }
            }
        }
    } else {
        // `cxx:1024-1067`: the not-closed arm.
        if p1_inf {
            // `cxx:1028-1034`.
            if v1.is_some() {
                return None;
            }
        } else {
            match v1.as_ref() {
                // `cxx:1037-1039`.
                None => v1 = Some(b.make_vertex(q1?, preci)),
                // `cxx:1040-1046`.
                Some(a) => {
                    let dist = q1.as_ref()?.distance(&BRepTool::vertex_point_world(a));
                    if dist > preci.max(BRepTool::vertex_tolerance(a)) {
                        return None;
                    }
                }
            }
        }
        if p2_inf {
            // `cxx:1048-1054`.
            if v2.is_some() {
                return None;
            }
        } else {
            match v2.as_ref() {
                // `cxx:1057-1059`.
                None => v2 = Some(b.make_vertex(q2?, preci)),
                // `cxx:1060-1066`.
                Some(a) => {
                    let dist = q2.as_ref()?.distance(&BRepTool::vertex_point_world(a));
                    if dist > preci.max(BRepTool::vertex_tolerance(a)) {
                        return None;
                    }
                }
            }
        }
    }
    Some(EdgeBuild { first: p1, last: p2, v1, v2, closed, reverse, degenerated: false })
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
