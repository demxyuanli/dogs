//! `ShapeAnalysis_Wire::CheckLoop` (`ShapeAnalysis_Wire.cxx:2205-2340`) and
//! `ShapeFix_Face::FindNext` / `isClosed2D` / `FixLoopWire`
//! (`ShapeFix_Face.cxx:2398-2641`).
//!
//! A wire whose path returns to a vertex it already passed and leaves again is
//! not a single boundary loop: the face it bounds has more than one boundary
//! wire. `ShapeFix_Face::Perform`'s second wire round (`ShapeFix_Face.cxx:583`)
//! calls `FixLoopWire` for every boundary wire and replaces the wire by the
//! split result.
//!
//! UNPORTED: `SendWarning` (`cxx:587`, needs `Message_Msg`) and the
//! `ShapeFix_Face::myStatus` bits `FixLoopWire` sets (`cxx:2578` / `cxx:2616`).

use std::collections::{HashMap, HashSet};

use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::ShapeType;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_surface::SurfaceKind;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::pcurve_full::classify_surface_kind;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, is_same};

use super::wire_fix::{check_small, first_vertex, fix_reorder_wire_3d, is_seam_use, last_vertex};

fn key(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// `ShapeAnalysis_Wire::CheckLoop` outputs (`ShapeAnalysis_Wire.cxx:2228-2311`):
/// the vertices reached by more than two "significant" edges, the edges stored
/// per vertex, and the small / seam edge sets the loop test ignores.
pub struct LoopCheck {
    /// `aMapVertices` (`NCollection_IndexedMap<TopoDS_Shape>`), in `Add` order.
    pub loop_vertices: Vec<TopoShape>,
    /// `aMapVertexEdges`: every edge of the wire per endpoint vertex, in the
    /// order `CheckLoop` appended them — a closed edge twice (`cxx:2270-2273`).
    pub vertex_edges: HashMap<usize, Vec<TopoShape>>,
    /// `aMapSmallEdges`: degenerated edges and closed edges `CheckSmall` accepts.
    pub small_edges: HashSet<usize>,
    /// `aMapSeemEdges`: edges the wire visits twice (`IsSeam`).
    pub seam_edges: HashSet<usize>,
}

/// `ShapeAnalysis_Wire::isMultiVertex` (`ShapeAnalysis_Wire.cxx:2205-2227`):
/// true when the edges at the vertex, ignoring the small and seam ones, are
/// more than two.
fn is_multi_vertex(list: &[TopoShape], small: &HashSet<usize>, seam: &HashSet<usize>) -> bool {
    let nb_not_account = list
        .iter()
        .filter(|e| small.contains(&key(e)) || seam.contains(&key(e)))
        .count();
    list.len() - nb_not_account > 2
}

/// `ShapeAnalysis_Wire::CheckLoop(aMapVertices, aMapVertexEdges, aMapSmallEdges,
/// aMapSeemEdges)` (`ShapeAnalysis_Wire.cxx:2228-2340`). `None` is the
/// function's `false` (`cxx:2236-2238` empty wire, `cxx:2249-2252` `FAIL2`).
pub fn check_loop(wire: &Wire, face: &Face) -> Option<LoopCheck> {
    let edges = edges_of_wire(wire);
    if edges.len() < 2 {
        return None;
    }
    // `cxx:2255-2256` `myWire->IsSeam(i)` (`ShapeExtend_WireData::IsSeam`).
    let seam_flags: Vec<bool> = (0..edges.len()).map(|i| is_seam_use(&edges, i)).collect();

    let mut small_edges: HashSet<usize> = HashSet::new();
    let mut seam_edges: HashSet<usize> = HashSet::new();
    let mut vertex_edges: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let mut loop_vertices: Vec<TopoShape> = Vec::new();
    let mut loop_keys: HashSet<usize> = HashSet::new();

    for (i, e) in edges.iter().enumerate() {
        // `cxx:2246-2247`: `TopExp::Vertices(aedge, aV1, aV2)`, CumOri on.
        let (Some(v1), Some(v2)) = (first_vertex(e), last_vertex(e)) else {
            return None; // cxx:2249-2252 ShapeExtend_FAIL2
        };
        let same = is_same(&v1.0, &v2.0);
        let ek = key(e);
        // `cxx:2255-2265`.
        if seam_flags[i] {
            seam_edges.insert(ek);
        } else if BRepTool::is_degenerated(e) {
            small_edges.insert(ek);
        } else if same && check_small(wire, face, BRepTool::vertex_tolerance(&v1), i + 1).ok {
            small_edges.insert(ek);
        }
        // `cxx:2264-2304`: register the edge under both endpoint vertices and
        // record a vertex that reaches more than two significant edges.
        let mut push_at = |vkey: usize, twice: bool, vtx: &TopoShape| {
            let list = vertex_edges.entry(vkey).or_default();
            list.push(e.0.clone());
            if twice {
                list.push(e.0.clone());
            }
            if list.len() > 2
                && is_multi_vertex(list.as_slice(), &small_edges, &seam_edges)
                && loop_keys.insert(vkey)
            {
                loop_vertices.push(vtx.clone());
            }
        };
        if same {
            // `cxx:2266-2275`: a closed edge reaches `aV1` twice.
            push_at(key(&v1.0), true, &v1.0);
        } else {
            push_at(key(&v1.0), false, &v1.0);
            push_at(key(&v2.0), false, &v2.0);
        }
    }

    if loop_vertices.is_empty() {
        return None;
    }
    Some(LoopCheck { loop_vertices, vertex_edges, small_edges, seam_edges })
}

/// `ShapeFix_Face::FindNext` (`ShapeFix_Face.cxx:2398-2456`): walk on from
/// `init_edge` through `vert` towards the far endpoint, appending every edge
/// (and its reversed view for a seam) to the wire being collected.
fn find_next(
    vert: &TopoShape,
    init_edge: &TopoShape,
    loop_keys: &HashSet<usize>,
    check: &LoopCheck,
    used: &mut HashSet<usize>,
    wd: &mut Vec<Edge>,
) {
    // `cxx:2409-2420`: the other vertex of `init_edge` (`TopoDS_Iterator` with
    // the default CumOri/CumLoc), staying on `vert` for a closed edge.
    let mut is_find = false;
    let mut next_vert = vert.clone();
    for c in cumulated_children(init_edge) {
        if c.shape_type() != ShapeType::Vertex {
            continue;
        }
        if !is_same(&c, vert) {
            is_find = true;
            next_vert = c;
            break;
        }
    }
    // `cxx:2421-2424`.
    if !is_find && !check.small_edges.contains(&key(init_edge)) {
        return;
    }
    // `cxx:2425-2428`.
    if is_find && loop_keys.contains(&key(&next_vert)) {
        return;
    }
    // `cxx:2430-2456`.
    let Some(aledges) = check.vertex_edges.get(&key(&next_vert)) else {
        return;
    };
    for e in aledges {
        let ek = key(e);
        if used.contains(&ek) || is_same(e, init_edge) {
            continue;
        }
        wd.push(Edge(e.clone()));
        if check.seam_edges.contains(&ek) {
            // `cxx:2440-2443`.
            let mut rev = e.clone();
            rev.reverse();
            wd.push(Edge(rev));
        }
        used.insert(ek);
        find_next(&next_vert, e, loop_keys, check, used, wd);
        break;
    }
}

/// `TopExp::Vertices(wire, Vfirst, Vlast)` (`TopExp.cxx:255-318`): the wire's
/// end vertices, where a vertex reached an even number of times cancels out.
/// A closed wire yields the same vertex for both, and more than two leftovers
/// yield none (the `vmap.Extent() > 2` case leaves both null).
fn wire_end_vertices(wire: &Wire) -> (Option<Vertex>, Option<Vertex>) {
    let mut odd: Vec<usize> = Vec::new();
    let mut last_last: Option<Vertex> = None;
    for e in edges_of_wire(wire) {
        let (Some(v1), Some(v2)) = (first_vertex(&e), last_vertex(&e)) else {
            return (None, None);
        };
        // `TopExp.cxx:272-286`: V1 enters FORWARD, V2 enters REVERSED, and
        // `TopTools_ShapeMapHasher` keys with `IsSame`, so the orientation only
        // decides the surviving view, never the parity.
        toggle(&mut odd, key(&v1.0));
        toggle(&mut odd, key(&v2.0));
        last_last = Some(v2);
    }
    if odd.is_empty() {
        // `TopExp.cxx:287-294` closed: both ends are the last edge's last vertex.
        let v = last_last;
        return (v.clone(), v);
    }
    if odd.len() == 2 {
        // `TopExp.cxx:296-316` open. Only `IsSame` is asked of the pair
        // downstream (`cxx:2534-2570`), so the views are reported as they are.
        let mut out: Vec<Vertex> = Vec::new();
        for k in &odd {
            for e in edges_of_wire(wire) {
                let (Some(v1), Some(v2)) = (first_vertex(&e), last_vertex(&e)) else {
                    return (None, None);
                };
                let v = if key(&v1.0) == *k {
                    Some(v1)
                } else if key(&v2.0) == *k {
                    Some(v2)
                } else {
                    None
                };
                if let Some(v) = v {
                    out.push(v);
                    break;
                }
            }
        }
        if out.len() == 2 {
            return (Some(out[0].clone()), Some(out[1].clone()));
        }
    }
    (None, None)
}

fn toggle(list: &mut Vec<usize>, k: usize) {
    if let Some(p) = list.iter().position(|&x| x == k) {
        list.remove(p);
    } else {
        list.push(k);
    }
}

/// `TopoDS_Shape::IsSame` on two possibly null views (`TopoDS_Shape.cxx`):
/// two nulls compare same.
fn same_view(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        _ => false,
    }
}

/// `ShapeFix_Wire`-level merge of `asfw->Load(asewd); asfw->FixReorder();
/// asfw->Wire()` (`cxx:2569-2577` / `cxx:2608-2615`): concatenate the two
/// wires' traversal edges and reorder the result in 3D.
fn merge_wires(a: &Wire, b: &Wire, builder: &TopoBuilder) -> Wire {
    let mut edges: Vec<Edge> = edges_of_wire(a);
    edges.extend(edges_of_wire(b));
    let w = builder.make_wire(&edges);
    fix_reorder_wire_3d(&w);
    w
}

/// `ShapeAnalysis_Wire::CheckGap2d(num)` status (`ShapeAnalysis_Wire.cxx:1156-1192`).
/// `Ok` is the `ShapeExtend_OK` (no flag) outcome — the 2D gap at that vertex
/// is within the surface resolution; `Done1` is the real gap; `Fail1` is a
/// missing face / surface / pcurve.
#[derive(Clone, Copy, PartialEq, Eq)]
enum GapStatus {
    Ok,
    Done1,
    Fail1,
}

/// `ShapeAnalysis_Wire::CheckGap2d(num)` (`ShapeAnalysis_Wire.cxx:1156-1192`).
/// `my_precision` is the analyzer precision (`SetPrecision`).
fn check_gap_2d(face: &Face, wire: &Wire, my_precision: f64, num: usize) -> GapStatus {
    let Some(surf) = BRepTool::face_surface(face) else {
        return GapStatus::Fail1; // cxx:1158-1162
    };
    let edges = edges_of_wire(wire);
    if edges.is_empty() {
        return GapStatus::Ok; // cxx:1166-1169
    }
    let nb = edges.len();
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    // `cxx:1177-1181`: `ShapeAnalysis_Edge::PCurve(E, myFace, C, f, l)`.
    let (Some((c1, _f1, l1)), Some((c2, f2, _l2))) = (
        curve_on_surface_oriented(&edges[n1 - 1], face, true),
        curve_on_surface_oriented(&edges[n2 - 1], face, true),
    ) else {
        return GapStatus::Fail1;
    };
    let p1 = c1.d0(l1);
    let p2 = c2.d0(f2);
    // `cxx:1184-1190`: `SA.UResolution(myPrecision)` / `VResolution`.
    let res = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), my_precision)
        .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), my_precision));
    if p1.distance(&p2) > res + PCONFUSION {
        GapStatus::Done1
    } else {
        GapStatus::Ok
    }
}

/// `ShapeFix_Face::isClosed2D(face, wire)` (`ShapeFix_Face.cxx:2458-2476`):
/// the wire is closed in 2D with the tolerance of each edge's first vertex,
/// checked by `ShapeAnalysis_Wire::CheckGap2d` on every edge, stopping at the
/// first gap.
fn is_closed_2d(face: &Face, wire: &Wire) -> bool {
    let edges = edges_of_wire(wire);
    for (i, e) in edges.iter().enumerate() {
        // `cxx:2464-2468`: `sae.FirstVertex(edge1)` then `SetPrecision`.
        let prec = match first_vertex(e) {
            Some(v) => BRepTool::vertex_tolerance(&v),
            None => CONFUSION,
        };
        // `cxx:2470-2471`: `LastCheckStatus(ShapeExtend_OK)` is `myStatus == 0`.
        if check_gap_2d(face, wire, prec, i + 1) != GapStatus::Ok {
            return false;
        }
    }
    true
}

/// `ShapeFix_Face::FixLoopWire(aResWires)` (`ShapeFix_Face.cxx:2478-2641`).
/// `None` is the function's `false`: no loop vertex, or the split wires are not
/// closed in 2D on a non-planar face.
pub fn fix_loop_wire(wire: &Wire, face: &Face) -> Option<Vec<Wire>> {
    let check = check_loop(wire, face)?;
    let loop_keys: HashSet<usize> = check.loop_vertices.iter().map(key).collect();
    let builder = TopoBuilder::new();
    let mut used: HashSet<usize> = HashSet::new();
    let mut res: Vec<Wire> = Vec::new();
    let mut seq: Vec<Wire> = Vec::new();

    // `cxx:2496-2548`: every loop vertex seeds one wire, extended by `FindNext`.
    for v in &check.loop_vertices {
        let Some(aledges) = check.vertex_edges.get(&key(v)).cloned() else {
            continue;
        };
        for e in &aledges {
            let ek = key(e);
            if used.contains(&ek) {
                continue;
            }
            let mut wd: Vec<Edge> = vec![Edge(e.clone())];
            if check.seam_edges.contains(&ek) {
                // `cxx:2511-2514`.
                let mut rev = e.clone();
                rev.reverse();
                wd.push(Edge(rev));
            }
            used.insert(ek);
            find_next(v, e, &loop_keys, &check, &mut used, &mut wd);
            // `cxx:2526-2529`: a lone small edge is no loop.
            if wd.len() == 1 && check.small_edges.contains(&key(&wd[0].0)) {
                continue;
            }
            let w = builder.make_wire(&wd);
            let (v1, v2) = wire_end_vertices(&w);
            if same_view(&v1, &v2) {
                // `cxx:2536-2544`: a closed loop keeps only its reordered form.
                let w2 = w.clone();
                fix_reorder_wire_3d(&w2);
                res.push(w2);
            } else {
                seq.push(w);
            }
        }
    }

    // `cxx:2550-2623`.
    if seq.len() == 1 {
        res.push(seq[0].clone());
    } else {
        // `cxx:2553-2586`: glue open wires sharing both end vertices.
        let mut i = 0usize;
        while i < seq.len() {
            let (v1, v2) = wire_end_vertices(&seq[i]);
            let mut matched = None;
            for j in (i + 1)..seq.len() {
                let (w1, w2) = wire_end_vertices(&seq[j]);
                if (same_view(&v1, &w1) || same_view(&v1, &w2))
                    && (same_view(&v2, &w1) || same_view(&v2, &w2))
                {
                    matched = Some(j);
                    break;
                }
            }
            match matched {
                // `cxx:2574-2581`: `aSeqWires.Remove(j--)` then `Remove(i--)`;
                // the outer `i++` leaves `i` on the element shifted into its slot.
                Some(j) => {
                    let merged = merge_wires(&seq[i], &seq[j], &builder);
                    res.push(merged);
                    seq.remove(j);
                    seq.remove(i);
                }
                None => i += 1,
            }
        }
        if seq.len() < 3 {
            // `cxx:2587-2592`.
            for w in &seq {
                res.push(w.clone());
            }
        } else {
            // `cxx:2593-2621`: with three or more wires left, glue by one
            // common end vertex, accumulating into the current wire.
            let mut i = 0usize;
            while i < seq.len() {
                let mut a_wire = seq[i].clone();
                let (mut v1, mut v2) = wire_end_vertices(&a_wire);
                let mut j = i + 1;
                while j < seq.len() {
                    let (w1, w2) = wire_end_vertices(&seq[j]);
                    if same_view(&v1, &w1)
                        || same_view(&v1, &w2)
                        || same_view(&v2, &w1)
                        || same_view(&v2, &w2)
                    {
                        a_wire = merge_wires(&a_wire, &seq[j], &builder);
                        let (n1, n2) = wire_end_vertices(&a_wire);
                        v1 = n1;
                        v2 = n2;
                        seq.remove(j);
                    } else {
                        j += 1;
                    }
                }
                res.push(a_wire);
                i += 1;
            }
        }
    }

    // `cxx:2626-2639`: on a non-planar face every split wire must be closed in
    // 2D, otherwise `FixLoopWire` reports false and the caller keeps the wire.
    let mut is_closed = true;
    if let Some(surf) = BRepTool::face_surface(face) {
        if classify_surface_kind(surf.as_ref()) != SurfaceKind::Plane {
            for w in &res {
                if !is_closed_2d(face, w) {
                    is_closed = false;
                    break;
                }
            }
        }
    }
    if res.is_empty() || !is_closed {
        return None;
    }
    Some(res)
}

/// `ShapeFix_Face::Perform`'s loop-wire step (`ShapeFix_Face.cxx:583-598`): run
/// [`fix_loop_wire`] on every boundary wire of `face` and, when a wire splits,
/// store the resulting wires in its place. OCCT rebuilds `tmpFace` with
/// `B.Add(tmpFace, aLoopWires.Value(k))` (`cxx:596`), i.e. the split wires are
/// added FORWARD; the port rewrites the face's child list in place so the
/// face's own `BRep_TFace` (surface, tolerance) and its edge pcurves stay
/// valid. Returns true when at least one wire was split.
pub fn split_loop_wires(face: &Face) -> bool {
    let children: Vec<TopoShape> = {
        let ts = face.0.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    let mut changed = false;
    let mut out: Vec<TopoShape> = Vec::with_capacity(children.len());
    for child in children {
        if child.shape_type() != ShapeType::Wire {
            out.push(child); // `cxx:533-537`: non-wire children pass through.
            continue;
        }
        let wire = Wire(child);
        match fix_loop_wire(&wire, face) {
            Some(wires) => {
                changed = true;
                for w in wires {
                    out.push(w.0);
                }
            }
            // `cxx:599-601`: keep the wire as it stands.
            None => out.push(wire.0),
        }
    }
    if changed {
        let mut ts = face.0.tshape.write().expect("poisoned TShape lock");
        ts.children = out;
    }
    changed
}
