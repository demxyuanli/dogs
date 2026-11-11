use super::next_param::wire_area_2d;
use super::pcurve::pcurve_and_range_on_face;
use super::wire_walk::wire_next_by_uv;

use super::super::prelude::*;

/// `ShapeAnalysis::OuterWire` over an already-oriented wire list.
pub fn outer_of_wires(wires: &[Wire], face: &Face) -> Option<Wire> {
    if wires.is_empty() {
        return None;
    }
    for (i, w) in wires.iter().enumerate() {
        if i == wires.len() - 1 {
            return Some(w.clone());
        }
        if wire_area_2d(w, face) >= 0.0 {
            return Some(w.clone());
        }
    }
    None
}

/// `ShapeExtend_WireData(wire, chained=true)` edge list: FORWARD/REVERSED
/// only; a REVERSED wire prepends so iterator order is reversed
/// (`ShapeExtend_WireData.cxx:115-122`).
pub(super) fn wire_edges_sewd(wire: &Wire) -> Vec<Edge> {
    let mut stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if wire.0.orientation().is_reversed() {
        stored.reverse();
    }
    stored
}

pub(super) fn edge_oriented_ends(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (a, b) = edge_vertices(e);
    if e.0.orientation().is_reversed() {
        (b, a)
    } else {
        (a, b)
    }
}

/// `BRepTools_WireExplorer` — the wire's edges in connection order on `face`.
///
/// Port of `BRepTools_WireExplorer::Init` (`BRepTools_WireExplorer.cxx:124-382`)
/// and `Next` (`cxx:394-705`). Init keys every edge by its cumulated first
/// vertex (`TopExp::Vertices(E, V1, V2, true)`, `cxx:244`) and starts from the
/// wire's first vertex (`cxx:303-373`). Next advances through
/// `TopExp::LastVertex(myEdge, true)` and takes the next edge from that vertex's
/// list. Storage order only decides the *first* edge: when several edges start
/// at the current vertex (`cxx:410-705`) OCCT runs the list through the
/// closed-edge-only pass, then the open-edge pass, and inside a pass picks the
/// edge whose pcurve start lies within `(myTolU, myTolV)` of the current 2D
/// point and whose direction turns least from the current edge's backward
/// tangent (`cxx:536-695`). Walking by vertex identity alone picks the wrong
/// occurrence of a seam edge on a closed surface (the cylinder lateral wire
/// would visit `seam(u = 2*pi)` instead of the bottom circle), which is what
/// makes the sampled UV ring self-intersect.
pub(crate) fn wire_edges_explorer(wire: &Wire, face: &Face) -> Vec<Edge> {
    // `BRepTools_WireExplorer::Init` leaves `myEdge` null in every path that
    // bails out here (`BRepTools_WireExplorer.cxx:160-164`, `:368-372`), so the
    // caller sees no edge at all.
    macro_rules! bail {
        ($why:expr) => {{
            let _ = $why;
            return Vec::new();
        }};
    }
    let edges: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o != Orientation::Internal && o != Orientation::External
        })
        .collect();
    if edges.is_empty() {
        bail!("no edges");
    }

    // `Init` (`cxx:139-217`): UV tolerances from the largest vertex tolerance of
    // the wire (falling back to the edge tolerances below Confusion).
    let mut df_vert_toler = 0.0f64;
    for e in &edges {
        let (a, b) = edge_oriented_ends(e);
        for v in [a, b].into_iter().flatten() {
            df_vert_toler = df_vert_toler.max(BRepTool::vertex_tolerance(&v));
        }
    }
    if df_vert_toler < CONFUSION {
        for e in &edges {
            df_vert_toler = df_vert_toler.max(BRepTool::edge_tolerance(e));
        }
        if df_vert_toler < CONFUSION {
            // empty wire (`cxx:160-164`)
            bail!("vert toler under confusion");
        }
    }
    let Some(surf) = BRepTool::face_surface_world(&face) else {
        bail!("no face surface");
    };
    let mut tol_u =
        2.0 * occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), df_vert_toler);
    let mut tol_v =
        2.0 * occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), df_vert_toler);
    {
        // The two-argument `BRepTools_WireExplorer(W, F)` entry point
        // (`cxx:100-116`) fills the UV bounds only for Cone / BSpline / Bezier
        // faces, from `BRepTools::UVBounds`; every other surface type keeps them
        // at 0.
        let (umin, umax, vmin, vmax) = match mesh_surface_type(surf.as_ref()) {
            MeshSurfaceType::Cone
            | MeshSurfaceType::BSplineSurface
            | MeshSurfaceType::BezierSurface => crate::brep_tools::uv_bounds(face),
            _ => (0.0, 0.0, 0.0, 0.0),
        };
        match mesh_surface_type(surf.as_ref()) {
            MeshSurfaceType::Cone => {
                // `cxx:170-200`.
                let maxtol = 0.0005 * (umax - umin);
                let a1 = surf.d1(umin, vmin).1.magnitude();
                let tol1 = if a1 <= CONFUSION {
                    maxtol
                } else {
                    maxtol.min(df_vert_toler / a1)
                };
                let a2 = surf.d1(umin, vmax).1.magnitude();
                let tol2 = if a2 <= CONFUSION {
                    maxtol
                } else {
                    maxtol.min(df_vert_toler / a2)
                };
                tol_u = 2.0 * tol1.max(tol2);
            }
            MeshSurfaceType::BSplineSurface | MeshSurfaceType::BezierSurface => {
                // `cxx:202-217`.
                let mut max_tol = tol_u.max(tol_v);
                let (_, du, dv) = surf.d1((umax - umin) / 2.0, (vmax - vmin) / 2.0);
                let modulus = (du.dot(&du) + dv.dot(&dv)).sqrt();
                if modulus > REAL_SMALL {
                    if modulus * max_tol / df_vert_toler < 1.5 {
                        max_tol = 1.5 * df_vert_toler / modulus;
                    }
                    tol_u = max_tol;
                    tol_v = max_tol;
                }
            }
            _ => {}
        }
    }

    // `myMap` / `vmap` (`cxx:226-302`): `myMap(V)` is the list of edges whose
    // cumulated first vertex is `V`, appended in wire order; `vmap` keeps the
    // vertices met an odd number of times (the open ends of the wire).
    let mut map: Vec<(usize, Vec<Edge>)> = Vec::new();
    let mut vmap: Vec<(usize, Orientation)> = Vec::new();
    for e in &edges {
        let (v1, v2) = edge_oriented_ends(e);
        if let Some(v1) = v1 {
            let key = vertex_key(&v1);
            match map.iter_mut().find(|(k, _)| *k == key) {
                Some((_, l)) => l.push(e.clone()),
                None => map.push((key, vec![e.clone()])),
            }
            vmap_toggle(&mut vmap, key, Orientation::Forward);
        }
        if let Some(v2) = v2 {
            vmap_toggle(&mut vmap, vertex_key(&v2), Orientation::Reversed);
        }
    }

    // First vertex (`cxx:303-373`): the wire's open end when it has one, else the
    // first vertex of the first edge.
    let first_key = if !vmap.is_empty() {
        vmap.iter()
            .find(|(_, o)| *o == Orientation::Forward)
            .map(|(k, _)| *k)
    } else {
        if let Some(e) = edges.iter().find(|e| {
            let (v1, v2) = edge_oriented_ends(e);
            v1.is_none() || v2.is_none()
        }) {
            // Infinite edge (`cxx:326-348`): the walk is that single edge.
            {
                let (f, l) = BRepTool::edge_parameters(e);
                let rev = e.0.orientation().is_reversed();
                let infinite = if rev {
                    l == f64::NEG_INFINITY || l == f64::INFINITY
                } else {
                    f == f64::NEG_INFINITY || f == f64::INFINITY
                };
                if infinite {
                    return vec![e.clone()];
                }
            }
        }
        edge_oriented_ends(&edges[0]).0.map(|v| vertex_key(&v))
    };
    let Some(first_key) = first_key else {
        bail!("no first vertex");
    };
    let Some(entry) = map.iter_mut().find(|(k, _)| *k == first_key) else {
        bail!("first vertex not in map");
    };
    let mut out = vec![entry.1.remove(0)];

    while let Some(lv) = edge_oriented_ends(out.last().unwrap()).1 {
        let Some(entry) = map.iter_mut().find(|(k, _)| *k == vertex_key(&lv)) else {
            break;
        };
        let list = &mut entry.1;
        if list.is_empty() {
            break;
        }
        if list.len() == 1 {
            // `cxx:414-482`: a single candidate, with the pcurve-end check for a
            // closed candidate edge.
            let candidate = list[0].clone();
            let (v1, v2) = edge_oriented_ends(&candidate);
            let Some(v1) = v1 else { break };
            if !same_vertex(&v1, &lv) {
                break;
            }
            // `cxx:429-475`: a closed candidate (`aV1.IsSame(aV2)`) is only
            // accepted when both pcurves exist and the current edge's 2D
            // endpoint is nearer the candidate's 2D start than its 2D end.
            if let Some(v2) = v2 {
                if same_vertex(&v1, &v2) {
                    let (Some((prev_pc, p11, p12)), Some((next_pc, p21, p22))) = (
                        pcurve_and_range_on_face(out.last().unwrap(), face),
                        pcurve_and_range_on_face(&candidate, face),
                    ) else {
                        break;
                    };
                    let prev_par = if out.last().unwrap().0.orientation().is_reversed() {
                        p11
                    } else {
                        p12
                    };
                    let (next_fpar, next_lpar) = if candidate.0.orientation().is_reversed() {
                        (p22, p21)
                    } else {
                        (p21, p22)
                    };
                    let p_prev = prev_pc.d0(prev_par);
                    let p_next_f = next_pc.d0(next_fpar);
                    let p_next_l = next_pc.d0(next_lpar);
                    if p_prev.square_distance(&p_next_f) > p_prev.square_distance(&p_next_l) {
                        break;
                    }
                }
            }
            list.clear();
            out.push(candidate);
        } else {
            // `cxx:521-705`: several candidates start at the current vertex.
            let current = out.last().unwrap().clone();
            let Some(idx) = wire_next_by_uv(&current, face, list, tol_u, tol_v) else {
                break;
            };
            out.push(list.remove(idx));
        }
    }
    out
}

/// `TopTools_ShapeMapHasher::IsSame` key for a vertex (`TopoDS_Shape`
/// identity = `TShape` + location).
pub(super) fn vertex_key(v: &Vertex) -> usize {
    GeometryRegistry::shape_key(&v.0)
}
pub(super) fn same_vertex(a: &Vertex, b: &Vertex) -> bool {
    vertex_key(a) == vertex_key(b)
}

/// `vmap.Add` / `RemoveKey` (`BRepTools_WireExplorer.cxx:252-268`): the map
/// compares shapes by `IsSame`, so a key added twice is removed again and only
/// the last orientation written survives.
pub(super) fn vmap_toggle(vmap: &mut Vec<(usize, Orientation)>, key: usize, ori: Orientation) {
    if let Some(pos) = vmap.iter().position(|(k, _)| *k == key) {
        vmap.remove(pos);
    } else {
        vmap.push((key, ori));
    }
}
