//! The file statics of `BRepFill_CompatibleWires.cxx`.
//!
//! Ported: `AddNewEdge` (`cxx:131-150`), `SeqOfVertices` (`cxx:152-172`),
//! `PlaneOfWire` (`cxx:174-298`), `TrimEdge` (`cxx:367-447`), `SearchRoot`
//! (`cxx:449-487`), `SearchVertex` (`cxx:489-525`), `Transform` (`cxx:682-734`)
//! and `BuildConnectedEdges` (`cxx:736-770`).
//!
//! Not ported: `EdgesFromVertex` (`cxx:63-129`). It sits inside
//! `#ifdef OCCT_DEBUG_EFV` together with its only caller (`ComputeOrigin`'s
//! debug block, `cxx:2128-2446`), so a default build never defines it.
//!
//! [`trim_edge`] also carries the body of `BRepFill.cxx`'s file static of the
//! same name (`BRepFill.cxx:143-207`): the two differ only by the
//! `Precision::Confusion` guards on a zero-length piece, which the
//! `reject_empty_piece` flag switches.

use std::collections::HashMap;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_gprop_full::{linear_properties, GaussProperties};
use crate::brep_lib_find_surface::FindSurface;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::{
    edge_vertices_cum_ori, first_vertex_cum_ori, last_vertex_cum_ori,
};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, is_same_occurrence, map_shapes_and_ancestors,
    shape_map_key, ShapeMapKey, wire_vertices,
};

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>` over
/// vertices, keeping both the key shape and the insertion order: `SearchRoot`
/// has to answer the *key* vertex, so a bare `HashMap<ShapeMapKey, _>` is not
/// enough.
pub(crate) type VertexMap = Vec<(Vertex, Vec<TopoShape>)>;

/// `Bind` of `NCollection_DataMap::Bind` — a bound key is overwritten.
pub(crate) fn vertex_map_bind(map: &mut VertexMap, key: &Vertex, value: Vec<TopoShape>) {
    let k = shape_map_key(&key.0);
    match map.iter_mut().find(|(existing, _)| shape_map_key(&existing.0) == k) {
        Some(slot) => slot.1 = value,
        None => map.push((key.clone(), value)),
    }
}

/// `Map(Key)` — the mutable value list of a bound key.
pub(crate) fn vertex_map_get_mut<'a>(
    map: &'a mut VertexMap,
    key: &Vertex,
) -> Option<&'a mut Vec<TopoShape>> {
    let k = shape_map_key(&key.0);
    map.iter_mut()
        .find(|(existing, _)| shape_map_key(&existing.0) == k)
        .map(|(_, value)| value)
}

/// `Map(Key)` on a shared map.
pub(crate) fn vertex_map_get<'a>(map: &'a VertexMap, key: &Vertex) -> Option<&'a [TopoShape]> {
    let k = shape_map_key(&key.0);
    map.iter()
        .find(|(existing, _)| shape_map_key(&existing.0) == k)
        .map(|(_, value)| value.as_slice())
}

/// `AddNewEdge(theEdge, theEdgeNewEdges, ListNewEdges)` (`cxx:131-150`): the
/// leaves of the split tree rooted at `edge`, in `NCollection_Sequence` order.
pub(crate) fn add_new_edge(
    edge: &TopoShape,
    edge_new_edges: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    list: &mut Vec<TopoShape>,
) {
    match edge_new_edges.get(&shape_map_key(edge)) {
        Some(new_edges) => {
            for e in new_edges {
                add_new_edge(e, edge_new_edges, list);
            }
        }
        None => list.push(edge.clone()),
    }
}

/// `SeqOfVertices(W, S)` (`cxx:152-172`): the wire's vertices in
/// `TopExp_Explorer` order, first occurrence of each kept.
pub(crate) fn seq_of_vertices(w: &Wire) -> Vec<Vertex> {
    let mut out: Vec<Vertex> = Vec::new();
    for s in explorer_shapes(&w.0, ShapeType::Vertex) {
        if !out.iter().any(|v| is_same(&v.0, &s)) {
            out.push(Vertex(s));
        }
    }
    out
}

/// `PlaneOfWire(W, P)` (`cxx:174-298`) answering both the `bool` result and the
/// out-parameter `P`.
///
/// `BRepLib_FindSurface::Init(W, -1, true)` (`cxx:237`) is the plane source for
/// a planar wire, moved onto the barycenter. A non-planar wire takes the
/// `cxx:250-293` arm: OCCT fills `P` from an inertia-axis frame but still
/// answers `false`. The symmetric sub-arm (`cxx:258-262`) leaves `P` untouched,
/// i.e. at the caller's `gp_Pln()`, which is returned here.
pub(crate) fn plane_of_wire_full(w: &Wire) -> (GpPln, bool) {
    // `cxx:188-198`: `W.Closed()` or the two boundary vertices equal.
    let mut w_closed = w.0.is_closed();
    if !w_closed {
        let (v1, v2) = wire_vertices(w);
        if let (Some(v1), Some(v2)) = (v1, v2) {
            w_closed = is_same(&v1.0, &v2.0);
        }
    }

    let mut bary: Option<GpPnt> = None;
    if w_closed {
        // `TopoDS_Iterator anIter(W)` (`cxx:205-214`): the direct children.
        let children = crate::bop_occt_util::iter_children(&w.0);
        if children.len() == 1 {
            // `GeomAdaptor_Curve AdC; AdC.Load(BRep_Tool::Curve(E, first, last))`
            // (`cxx:212-216`): a `Geom_TrimmedCurve` forwards the type of its
            // basis, so `gp_circ` / `gp_ellipse` answer on the trim too.
            let e = Edge(children[0].clone());
            if let Some(c) = BRepTool::edge_curve_world(&e) {
                if let Some(circ) = c.gp_circ() {
                    bary = Some(circ.location());
                } else if let Some(el) = c.gp_ellipse() {
                    bary = Some(el.location());
                }
            }
        }
    }

    // `cxx:231-235`: `BRepGProp::LinearProperties(W, GP); Bary = GP.CentreOfMass()`.
    let props = linear_properties(&w.0).unwrap_or_else(|_| {
        // An empty accumulator is `GProp_GProps(P)` before any `Add`
        // (`BRepGProp.cxx:126-129`), i.e. the origin.
        GaussProperties::new(GpPnt::zero())
    });
    let bary = bary.unwrap_or_else(|| props.center());

    // `findPlanarSurf.Init(W, -1, true)` (`cxx:237`); `OnlyClosed` keeps its
    // `false` default (`BRepLib_FindSurface.hxx:76-79`).
    let find = FindSurface::new(&w.0, -1.0, true, false);
    let Some(surface) = find.surface() else {
        // `cxx:250-293`: not planar.
        return (non_planar_frame(&props).unwrap_or_default(), false);
    };
    // `cxx:240-247`: `L` is the identity location here, so the surface is used
    // as found; `P.SetLocation(Bary)` then moves the plane onto the barycenter.
    match surface.gp_pln() {
        Some(mut p) => {
            p.set_location(&bary);
            (p, true)
        }
        // `down_cast<Geom_Plane>` on a non-plane surface: `gp_Pln` stays as the
        // caller's default.
        None => (GpPln::default(), false),
    }
}

/// `PlaneOfWire(W, P)` with only the `bool` result.
pub(crate) fn plane_of_wire(w: &Wire) -> Option<GpPln> {
    let (p, is_plane) = plane_of_wire_full(w);
    if is_plane {
        Some(p)
    } else {
        None
    }
}

/// The `cxx:250-293` inertia frame of a non-planar wire, or `None` when the
/// section is symmetric (all three radii of gyration equal, `cxx:258-262`).
pub(crate) fn non_planar_frame(props: &GaussProperties) -> Option<GpPln> {
    let pp = props.principal_properties().ok()?;
    let (r1, r2, r3) = pp.radius_of_gyration();
    let tol = CONFUSION;
    let r_max = r1.max(r2).max(r3);
    if (r_max - r1).abs() < tol && (r_max - r2).abs() < tol
        || (r_max - r1).abs() < tol && (r_max - r3).abs() < tol
        || (r_max - r2).abs() < tol && (r_max - r3).abs() < tol
    {
        return None;
    }
    // `cxx:265-289`: normal from the largest-moment axis, X from the smallest.
    let n_dir = if r1 >= r2 && r1 >= r3 {
        pp.first_axis_of_inertia()
    } else if r2 >= r1 && r2 >= r3 {
        pp.second_axis_of_inertia()
    } else {
        pp.third_axis_of_inertia()
    };
    let x_dir = if r3 <= r2 && r3 <= r1 {
        pp.third_axis_of_inertia()
    } else if r2 <= r1 && r2 <= r3 {
        pp.second_axis_of_inertia()
    } else {
        pp.first_axis_of_inertia()
    };
    let n_dir = GpDir::from_vec(&n_dir).ok()?;
    let x_dir = GpDir::from_vec(&x_dir).ok()?;
    // `gp_Ax3 repere(Bary, NDir, XDir); Geom_Plane GPlan(repere); P = GPlan.Pln();`
    let bary = props.center();
    let ax3 = GpAx3::new(bary, n_dir, &x_dir).ok()?;
    Some(GpPln::new(ax3))
}

/// `TrimEdge(CurrentEdge, CutValues, t0, t1, SeqOrder, S)`
/// (`BRepFill_CompatibleWires.cxx:367-447`, `BRepFill.cxx:143-207`): cut
/// `current_edge` at every value of `cut_values`, which are abscissae on the
/// wire's `[t0, t1]` total length and are mapped onto the edge's own
/// `[first, last]` range.
///
/// `reject_empty_piece` is `true` for the `BRepFill_CompatibleWires` copy,
/// whose `abs(m0 - m1) < Precision::Confusion()` guards return early
/// (`cxx:392-395`, `:404-407`, `:424-427`, `:436-439`); the `BRepFill.cxx` copy
/// has no such guard.
///
/// Deviation: `BRepLib_MakeEdge` can leave `IsDone()` false (`ParameterOutOfRange`,
/// `BRepLib_MakeEdge.cxx:655-658`), in which case OCCT appends a null edge to
/// `S` and reads its null vertices on the next iteration. A null edge cannot be
/// represented here, so the cut stops with the pieces produced so far.
pub(crate) fn trim_edge(
    current_edge: &Edge,
    cut_values: &[f64],
    t0: f64,
    t1: f64,
    seq_order: bool,
    reject_empty_piece: bool,
) -> Vec<Edge> {
    let mut out: Vec<Edge> = Vec::new();
    let Some(curve) = BRepTool::edge_curve_world(current_edge) else {
        return out;
    };
    let (first, last) = BRepTool::edge_parameters(current_edge);
    let (vf, vl) = edge_vertices(current_edge);
    let (Some(vf), Some(vl)) = (vf, vl) else {
        return out;
    };
    let current_orient = current_edge.0.orientation();
    let n_dec = cut_values.len();
    let builder = TopoBuilder::new();

    // `BRepLib_MakeEdge(C, V1, V2, p1, p2)`; `None` is `!IsDone()`.
    let make = |v1: Option<&Vertex>, v2: Option<&Vertex>, p1: f64, p2: f64| -> Option<Edge> {
        let build = init_3d_vertices(&curve, v1, v2, p1, p2)?;
        let mut e = builder.make_edge(curve.clone(), build.first, build.last);
        add_endpoint_vertices(&builder, &mut e, &build);
        Some(e)
    };
    // `CutE.Orientation(CurrentOrient)`.
    let oriented = |mut e: Edge| -> Edge {
        e.0.set_orientation(current_orient);
        e
    };

    if seq_order {
        // from first to last (`cxx:387-416`)
        let mut m0 = first;
        let mut v0 = vf;
        for j in 0..n_dec {
            let m1 = (cut_values[j] - t0) * (last - first) / (t1 - t0) + first;
            if reject_empty_piece && (m0 - m1).abs() < CONFUSION {
                return out;
            }
            let Some(cut_e) = make(Some(&v0), None, m0, m1) else {
                return out;
            };
            let cut_e = oriented(cut_e);
            out.push(cut_e.clone());
            m0 = m1;
            if let Some(v) = last_vertex_cum_ori(&cut_e) {
                v0 = v;
            }
            if j == n_dec - 1 {
                if reject_empty_piece && (m0 - last).abs() < CONFUSION {
                    return out;
                }
                let Some(last_e) = make(Some(&v0), Some(&vl), m0, last) else {
                    return out;
                };
                out.push(oriented(last_e));
            }
        }
    } else {
        // from last to first (`cxx:420-447`)
        let mut m1 = last;
        let mut v1 = vl;
        for j in (0..n_dec).rev() {
            let m0 = (cut_values[j] - t0) * (last - first) / (t1 - t0) + first;
            if reject_empty_piece && (m0 - m1).abs() < CONFUSION {
                return out;
            }
            let Some(cut_e) = make(None, Some(&v1), m0, m1) else {
                return out;
            };
            let cut_e = oriented(cut_e);
            out.push(cut_e.clone());
            m1 = m0;
            if let Some(v) = first_vertex_cum_ori(&cut_e) {
                v1 = v;
            }
            if j == 0 {
                if reject_empty_piece && (first - m1).abs() < CONFUSION {
                    return out;
                }
                let Some(last_e) = make(Some(&vf), Some(&v1), first, m1) else {
                    return out;
                };
                out.push(oriented(last_e));
            }
        }
    }
    out
}

/// `SearchRoot(V, Map, VRoot)` (`cxx:449-487`): the key of `map` whose value
/// list contains `v`.
pub(crate) fn search_root(map: &VertexMap, v: &Vertex) -> Option<Vertex> {
    for (key, list) in map {
        if list.iter().any(|s| is_same(s, &v.0)) {
            return Some(key.clone());
        }
    }
    None
}

/// `SearchVertex(List, W, VonW)` (`cxx:489-525`): the first vertex of `w` (in
/// [`seq_of_vertices`] order) that also belongs to `list`.
pub(crate) fn search_vertex(list: &[TopoShape], w: &Wire) -> Option<Vertex> {
    for vi in seq_of_vertices(w) {
        if list.iter().any(|s| is_same(s, &vi.0)) {
            return Some(vi);
        }
    }
    None
}

/// `Transform(WithRotation, P, Pos1, Ax1, Pos2, Ax2, Pnew)` (`cxx:682-734`):
/// carry `p` from the plane `(pos1, ax1)` to the plane `(pos2, ax2)`.
pub(crate) fn transform(
    with_rotation: bool,
    p: &GpPnt,
    pos1: &GpPnt,
    ax1: &GpVec,
    pos2: &GpPnt,
    ax2: &GpVec,
) -> GpPnt {
    let mut pnew = p.translated_pnts(pos1, pos2);
    let mut axe1 = *ax1;
    let mut axe2 = *ax2;
    if !axe1.is_parallel_ang(&axe2, 1.0e-4) {
        let vtrans = GpVec::from_pnts(pos1, pos2);
        // `cxx:697-706`: orient both axes along the translation. OCCT's
        // `alpha` / `beta` are recomputed from `Vsign` at `cxx:711-713`, so the
        // first pair only feeds the two sign tests above them.
        if vtrans.dot(&axe1) < -1.0e-7 {
            axe1 = axe1.reversed();
        }
        if vtrans.dot(&axe2) < 1.0e-7 {
            axe2 = axe2.reversed();
        }
        let norm2 = axe1.crossed(&axe2);
        // `Vsign.SetLinearForm(Vtrans.Dot(axe1), axe2, -Vtrans.Dot(axe2), axe1)`
        let vsign =
            GpVec::zero().set_linear_form_2(vtrans.dot(&axe1), &axe2, -vtrans.dot(&axe2), &axe1);
        let alpha = vsign.dot(&axe1);
        let beta = vsign.dot(&axe2);
        let pasnul = alpha.abs() > 1.0e-4 && beta.abs() > 1.0e-4;
        let mut sign = 1.0;
        if alpha * beta > 0.0 && pasnul {
            sign = -1.0;
        }
        let dir = GpDir::from_vec(&norm2).ok();
        let Some(dir) = dir else {
            return pnew;
        };
        let norm = GpAx1::new(*pos2, dir);
        let mut ang = axe1.angle_with_ref(&axe2, &norm2);
        if !with_rotation {
            if ang > std::f64::consts::FRAC_PI_2 {
                ang -= std::f64::consts::PI;
            }
            if ang < -std::f64::consts::FRAC_PI_2 {
                ang += std::f64::consts::PI;
            }
        }
        ang *= sign;
        pnew = pnew.rotated(&norm, ang);
    }
    pnew
}

/// `BuildConnectedEdges(aWire, StartEdge, StartVertex, ConnectedEdges)`
/// (`cxx:736-770`): walk the wire from `start_edge` through shared vertices,
/// away from `start_vertex`, until the edge that closes the loop is found.
///
/// Deviation: when `MapVE` has no other edge at `CurVertex`, OCCT's `for(;;)`
/// never terminates (`FindFromKey` on an unbound key throws, an empty list
/// spins); this port stops the walk there.
pub(crate) fn build_connected_edges(
    wire: &Wire,
    start_edge: &Edge,
    start_vertex: &Vertex,
) -> Vec<TopoShape> {
    let map = map_shapes_and_ancestors(&wire.0, ShapeType::Vertex, ShapeType::Edge);
    let mut cur_edge = start_edge.clone();
    let mut cur_vertex = start_vertex.clone();
    let (v1, v2) = edge_vertices(start_edge);
    let origin = match (&v1, &v2) {
        (Some(a), Some(b)) => {
            if is_same(&a.0, &start_vertex.0) {
                Some(b.clone())
            } else {
                Some(a.clone())
            }
        }
        _ => None,
    };
    let mut out: Vec<TopoShape> = Vec::new();
    loop {
        let mut advanced = false;
        if let Some(list) = map.get(&shape_map_key(&cur_vertex.0)) {
            for e in list {
                if !is_same_occurrence(e, &cur_edge.0) {
                    out.push(e.clone());
                    let (v1, v2) = edge_vertices(&Edge(e.clone()));
                    cur_vertex = match (&v1, &v2) {
                        (Some(a), Some(b)) => {
                            if is_same(&a.0, &cur_vertex.0) {
                                b.clone()
                            } else {
                                a.clone()
                            }
                        }
                        _ => cur_vertex.clone(),
                    };
                    cur_edge = Edge(e.clone());
                    advanced = true;
                    break;
                }
            }
        }
        if let Some(o) = &origin {
            if is_same(&cur_vertex.0, &o.0) {
                break;
            }
        }
        if !advanced {
            break;
        }
    }
    out
}

/// One `BRepExtrema_DistShapeShape(ME, W)` solution: the pair of a point on the
/// line edge `ME` and a point on the wire.
struct WireSolution {
    distance: f64,
    /// `PointOnShape2(isol)`.
    point2: GpPnt,
    /// `SupportTypeShape2(isol) == BRepExtrema_IsVertex`.
    support_is_vertex: bool,
    /// `SupportOnShape2(isol)` when the support is a wire edge.
    support_edge: Option<Edge>,
    /// `ParOnEdgeS2(isol)` when the support is a wire edge.
    param: f64,
}

/// `BRepExtrema_DistShapeShape(ME, W)` (`BRepExtrema_DistShapeShape.cxx`), for
/// the one edge / wire pair [`edge_intersect_on_wire`] needs, keeping the
/// solutions inside `Precision::Confusion()` of the minimum
/// (`cxx:966-973`).
///
/// Re-derivation, matching [`crate::brep_fill::search_origin`]'s treatment of
/// the vertex / wire pair: the class explodes both shapes into their vertices
/// and edges (`BRepExtrema_ShapeList`) and computes vertex-vertex,
/// vertex-edge and edge-edge pairs. Here the three families are
///
/// * every wire vertex against `ME` (`Extrema_ExtPC` over `[0, max_param]`),
///   whose support on the wire is that vertex;
/// * both endpoints of `ME` against every wire edge (`Extrema_ExtPC` over the
///   edge's own range), whose support on the wire is that edge;
/// * `ME` against every wire edge (`Extrema_ExtCC` over both ranges), whose
///   support on the wire is that edge.
///
/// The vertex-vertex family of OCCT is the zero-length case of the first one
/// and adds no solution of its own. `NbSolution()` order is the order these
/// families append, which only matters for an exact tie in the caller's
/// "closest to P2" loop.
fn dist_edge_wire(
    line: &occt_geom::line::GeomLine,
    max_param: f64,
    wire: &Wire,
) -> Vec<WireSolution> {
    use occt_geom::extrema_cc::curve_curve_extrema_all_range;
    use occt_geom::extrema_pc::point_curve_extrema_range;

    let mut out: Vec<WireSolution> = Vec::new();
    let line_curve: &dyn occt_geom::curve::Curve = line;

    // Vertex family: `PointOnShape2` is the vertex itself.
    for s in explorer_shapes(&wire.0, ShapeType::Vertex) {
        let p = BRepTool::vertex_point_world(&Vertex(s));
        let Some(pair) = point_curve_extrema_range(line_curve, &p, 0.0, max_param) else {
            continue;
        };
        out.push(WireSolution {
            distance: pair.distance,
            point2: p,
            support_is_vertex: true,
            support_edge: None,
            param: pair.u1,
        });
    }

    // Edge families.
    let me_start = line.d0(0.0);
    let me_end = line.d0(max_param);
    for s in explorer_shapes(&wire.0, ShapeType::Edge) {
        let e = Edge(s);
        let (first, last) = BRepTool::edge_parameters(&e);
        let Some(c) = BRepTool::edge_curve_world(&e) else {
            continue;
        };
        for p in [&me_start, &me_end] {
            if let Some(pair) = point_curve_extrema_range(c.as_ref(), p, first, last) {
                out.push(WireSolution {
                    distance: pair.distance,
                    point2: pair.p2,
                    support_is_vertex: false,
                    support_edge: Some(e.clone()),
                    param: pair.u1,
                });
            }
        }
        for pair in curve_curve_extrema_all_range(line_curve, c.as_ref(), 0.0, max_param, first, last)
        {
            out.push(WireSolution {
                distance: pair.distance,
                point2: pair.p2,
                support_is_vertex: false,
                support_edge: Some(e.clone()),
                param: pair.u2,
            });
        }
    }
    out
}

/// The outcome of [`edge_intersect_on_wire`]: OCCT's `bool` return plus the two
/// in / out parameters.
pub(crate) struct EdgeIntersect {
    /// `NewVertex`.
    pub new_vertex: bool,
    /// `Vsol`, null in OCCT when the solution has no vertex support.
    pub vsol: Option<Vertex>,
    /// `newW`, only set when `new_vertex` is true.
    pub new_wire: Option<Wire>,
}

/// `EdgeIntersectOnWire(P1, P2, percent, Map, W, Vsol, newW, theEdgeNewEdges)`
/// (`BRepFill_CompatibleWires.cxx:527-680`): where the line from `P1` to `P2`
/// meets the wire `W`, cutting the met edge when that point is not already a
/// vertex of `W` (or a vertex already present in `Map`).
pub(crate) fn edge_intersect_on_wire(
    p1: &GpPnt,
    p2: &GpPnt,
    percent: f64,
    map: &VertexMap,
    w: &Wire,
    edge_new_edges: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> EdgeIntersect {
    let mut result = EdgeIntersect { new_vertex: false, vsol: None, new_wire: None };

    // `gp_Lin droite(P1, gp_Dir(gp_Vec(P1, P2)))` (`cxx:544`); a null direction
    // is `gp_Dir`'s `Standard_ConstructionError`, and every caller guards the
    // distance before calling.
    let Some(dir) = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).ok() else {
        return result;
    };
    let line = occt_geom::line::GeomLine::new(occt_core::gp::GpLin::from_pnt_dir(*p1, dir));

    // `Bnd_Box B; BRepBndLib::Add(W, B); ... diag = BP1.Distance(BP2)`
    // (`cxx:552-560`).
    let mut b = occt_core::bnd::box3d::BndBox::new();
    crate::brep_bnd_lib::add(&w.0, &mut b, false);
    let (bp1, bp2) = (b.corner_min(), b.corner_max());
    let dernier_param = bp1.distance(&bp2);

    // `BRepExtrema_DistShapeShape DSS(ME.Edge(), W)` (`cxx:565`).
    let solutions = dist_edge_wire(&line, dernier_param, w);
    let Some(min_d) = solutions.iter().map(|s| s.distance).reduce(f64::min) else {
        // `!DSS.IsDone()`: no shape pair, no solution.
        return result;
    };

    // `cxx:568-585`: keep the solution closest to `P2`.
    let mut best: Option<&WireSolution> = None;
    let mut dss = f64::INFINITY;
    for s in &solutions {
        if s.distance > min_d + CONFUSION {
            continue;
        }
        let d = p2.distance(&s.point2);
        if d < dss {
            dss = d;
            best = Some(s);
        }
    }
    let Some(best) = best else {
        return result;
    };

    // `NewVertex = (DSS.SupportTypeShape2(isol) != BRepExtrema_IsVertex)`
    // (`cxx:587`).
    let mut new_vertex = !best.support_is_vertex;
    let mut vsol: Option<Vertex> = None;
    let param = best.param;
    if new_vertex {
        let Some(e) = best.support_edge.clone() else {
            return result;
        };
        // `cxx:590-608`: an endpoint hit is not a new vertex.
        let (first, last) = BRepTool::edge_parameters(&e);
        let tol = PCONFUSION.max(percent * (last - first).abs());
        let (vf, vl) = edge_vertices_cum_ori(&e);
        if (first - param).abs() < tol {
            new_vertex = false;
            vsol = vf;
        } else if (last - param).abs() < tol {
            new_vertex = false;
            vsol = vl;
        }
        if !new_vertex {
            // `if (SearchRoot(Vsol, Map, VRoot)) NewVertex = true;` (`cxx:609-612`).
            if let Some(v) = &vsol {
                if search_root(map, v).is_some() {
                    new_vertex = true;
                }
            }
        }
    } else {
        // `Vsol = TopoDS::Vertex(DSS.SupportOnShape2(isol))` (`cxx:618-619`).
        // The support is the wire vertex the solution was built from, which is
        // the point itself.
        vsol = explorer_shapes(&w.0, ShapeType::Vertex)
            .into_iter()
            .map(Vertex)
            .find(|v| BRepTool::vertex_point_world(v).distance(&best.point2) <= CONFUSION);
    }

    result.vsol = vsol.clone();
    if !new_vertex {
        result.new_wire = Some(w.clone());
        return result;
    }

    // `cxx:625-671`: cut the supporting edge and rebuild the wire.
    let Some(e) = best.support_edge.clone() else {
        return result;
    };
    edge_new_edges.insert(shape_map_key(&e.0), Vec::new());
    let (first, last) = BRepTool::edge_parameters(&e);
    let (ef, el) = edge_vertices_cum_ori(&e);
    let mut mw = crate::brep_lib_make_wire::MakeWire::new();
    let mut exp = crate::brep_tools_wire_explorer::WireExplorer::new();
    exp.init(w, None);
    while exp.more() {
        let Some(cur) = exp.current_edge() else {
            exp.next();
            continue;
        };
        if is_same(&cur.0, &e.0) {
            // `SO = anExp.CurrentVertex().IsSame(TopExp::FirstVertex(E))`
            let so = match (exp.current_vertex(), &ef) {
                (Some(cv), Some(v0)) => is_same(&cv.0, &v0.0),
                _ => false,
            };
            let se = trim_edge(&e, &[param], first, last, so, true);
            if se.is_empty() {
                // `return false` (`cxx:647-650`) — the pieces produced so far
                // are dropped with the wire, as in OCCT.
                return EdgeIntersect { new_vertex: false, vsol, new_wire: None };
            }
            edge_new_edges.insert(shape_map_key(&e.0), se.iter().map(|x| x.0.clone()).collect());
            // `cxx:651-660`: `Vsol` becomes the endpoint the cut did not keep.
            let (vv1, vv2) = edge_vertices(&se[0]);
            if let (Some(vv1), Some(vv2)) = (vv1, vv2) {
                let first_hit = |v: &Vertex| {
                    ef.as_ref().map_or(false, |f| is_same(&f.0, &v.0))
                        || el.as_ref().map_or(false, |l| is_same(&l.0, &v.0))
                };
                if first_hit(&vv1) {
                    vsol = Some(vv2.clone());
                }
                if first_hit(&vv2) {
                    vsol = Some(vv1.clone());
                }
            }
            for piece in &se {
                let _ = mw.add(piece);
            }
        } else {
            let _ = mw.add(&cur);
        }
        exp.next();
    }

    EdgeIntersect { new_vertex: true, vsol, new_wire: Some(mw.wire()) }
}

