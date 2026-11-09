use super::prelude::*;
use super::acr::{cut_edge_at, make_wire_from};


/// One `BRepExtrema_DistShapeShape(V, W)` solution: the distance, whether the
/// support on the wire is a vertex, and the supporting sub-shape.
struct DistSolution {
    /// `P.Distance(DSS.PointOnShape2(isol))` (`cxx:908`): the solution's
    /// `PointOnShape2` is only ever read through this distance.
    distance: f64,
    is_vertex: bool,
    shape: TopoShape,
}


/// Re-derivation of `BRepExtrema_DistShapeShape(V, W)` for a vertex `V` and a
/// wire `W` (`BRepFill.cxx:904`). See the module notes: the two solution
/// families are the wire's vertices (vertex-vertex) and the wire's edges
/// (point-curve extrema). The minimum wins; a vertex solution wins a tie, which
/// is the order `BRepExtrema_ShapeList` builds for a wire (vertices before
/// edges).
fn dist_vertex_wire(p: &GpPnt, wire: &Wire) -> Option<DistSolution> {
    let mut best: Option<DistSolution> = None;
    for shape in crate::topo_tools::explorer_shapes(&wire.0, crate::abs::ShapeType::Vertex) {
        let v = Vertex(shape.clone());
        let point = BRepTool::vertex_point_world(&v);
        let distance = point.distance(p);
        let better = match &best {
            None => true,
            Some(b) => distance < b.distance,
        };
        if better {
            best = Some(DistSolution {
                distance,
                is_vertex: true,
                shape,
            });
        }
    }
    for shape in crate::topo_tools::explorer_shapes(&wire.0, crate::abs::ShapeType::Edge) {
        let e = Edge(shape.clone());
        let (u, point) = crate::brep_extrema::closest_point_on_edge(&e, p, 0);
        if !u.is_finite() {
            continue;
        }
        let distance = point.distance(p);
        let better = match &best {
            None => true,
            Some(b) => distance < b.distance,
        };
        if better {
            best = Some(DistSolution {
                distance,
                is_vertex: false,
                shape,
            });
        }
    }
    best
}


/// `BRepFill::SearchOrigin(W, P, Dir, Tol)` (`BRepFill.cxx:885-1023`): move the
/// seam of the closed wire `w` onto `p`, so the wire starts where the frame
/// search should begin, and orient the parse direction towards `dir`.
///
/// `w` must be closed; OCCT throws `Standard_NoSuchObject` otherwise. `p` is
/// projected onto the wire. When the projection lands inside an edge, that edge
/// is cut at the projected parameter (`CutEdge`), otherwise the vertex that
/// already sits at `p` becomes the seam.
pub fn search_origin(w: &mut Wire, p: &GpPnt, dir: &GpVec, tol: f64) -> Result<(), String> {
    if !w.0.is_closed() {
        // `throw Standard_NoSuchObject("BRepFill::SearchOrigin : the wire must
        //  be closed")` (`cxx:887-890`).
        return Err("Standard_NoSuchObject: BRepFill::SearchOrigin : the wire must be closed".into());
    }

    let mut new_vertex = false;
    let mut theparam = 1.0e101;
    let mut v: Option<Vertex> = None;
    let mut e_ref: Option<Edge> = None;

    // `W.Orientation(TopAbs_FORWARD)` (`cxx:900`).
    w.0.set_orientation(Orientation::Forward);

    // `B.MakeVertex(V, P, Tol); BRepExtrema_DistShapeShape DSS(V, W)`
    // (`cxx:903-904`).
    if let Some(sol) = dist_vertex_wire(p, w) {
        // `cxx:907-916`: the loop only keeps the minimum solution.
        let _dss = sol.distance;
        if sol.is_vertex {
            // `V = TopoDS::Vertex(supp)` (`cxx:918-921`).
            v = Some(Vertex(sol.shape));
        } else {
            // `cxx:923-946`.
            let e = Edge(sol.shape);
            let (vf, vl) = edge_vertices(&e);
            let (Some(vf), Some(vl)) = (vf, vl) else {
                return Err("Standard_NullObject: TopExp::Vertices".to_string());
            };
            let dist = p.distance(&BRepTool::vertex_point_world(&vf));
            if dist < tol {
                v = Some(vl.clone());
            }
            let d = p.distance(&BRepTool::vertex_point_world(&vl));
            if d < tol && d < dist {
                v = Some(vf.clone());
                // `dist = d;` -- OCCT keeps this local, it feeds `NewVertex`.
            }
            let dist = if d < tol && d < dist { d } else { dist };
            new_vertex = dist > tol;
            if new_vertex {
                // `DSS.ParOnEdgeS2(isol, theparam)` (`cxx:944`).
                let (u, _) = crate::brep_extrema::closest_point_on_edge(&e, p, 0);
                theparam = u;
            }
            e_ref = Some(e);
        }
    }

    // `cxx:955-968`: count the edges; a new vertex adds one (the cut).
    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(w, None);
    let mut nb_edges = 0usize;
    while wexp.more() {
        nb_edges += 1;
        wexp.next();
    }
    if new_vertex {
        nb_edges += 1;
    }

    // `cxx:971-992`: the edge table and `rangdeb`.
    let mut rangdeb = 0usize;
    let mut edges: Vec<TopoShape> = Vec::with_capacity(nb_edges);
    wexp.init(w, None);
    let mut ii = 1usize;
    while wexp.more() {
        let Some(e) = wexp.current_edge() else {
            wexp.next();
            continue;
        };
        let split = match &e_ref {
            Some(r) => new_vertex && is_same(&e.0, &r.0),
            None => false,
        };
        if split {
            let current_vertex = wexp.current_vertex().cloned();
            let Some(current_vertex) = current_vertex else {
                return Err("Standard_NullObject: BRepTools_WireExplorer::CurrentVertex".into());
            };
            let (e1, e2) = cut_edge_at(&e, theparam, &current_vertex)?;
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e1.0.clone());
            ii += 1;
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e2.0.clone());
            rangdeb = ii;
        } else {
            if ii > nb_edges {
                return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges".to_string());
            }
            edges.push(e.0.clone());
        }
        if !new_vertex {
            if let (Some(v), Some(cv)) = (&v, wexp.current_vertex()) {
                if is_same(&v.0, &cv.0) {
                    rangdeb = ii;
                }
            }
        }
        ii += 1;
        wexp.next();
    }
    if rangdeb == 0 {
        rangdeb = nb_edges;
    }

    // `cxx:998-1018`: the direction of parsing.
    if rangdeb == 0 || rangdeb > edges.len() {
        return Err("Standard_OutOfRange: BRepFill::SearchOrigin Edges(rangdeb)".to_string());
    }
    let e = Edge(edges[rangdeb - 1].clone());
    if !new_vertex {
        // `theparam = BRep_Tool::Parameter(V, E)` (`cxx:1003`).
        let Some(v) = &v else {
            return Err("Standard_NullObject: BRep_Tool::Parameter".to_string());
        };
        theparam = BRepTool::parameter_of_vertex(&e.0, &v.0);
    }
    // `BRepAdaptor_Curve AC(E); AC.D1(theparam, Pe, Ve)` (`cxx:1005-1008`).
    let curve = BRepTool::edge_curve_world(&e)
        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
    let (_pe, mut ve) = curve.d1(theparam);
    if e.0.orientation() == Orientation::Reversed {
        ve = ve.multiplied_scalar(-1.0);
    }
    let mut angle = ve.angle(dir);
    if angle > std::f64::consts::PI {
        angle = 2.0 * std::f64::consts::PI - angle;
    }
    let forward = angle <= std::f64::consts::PI / 2.0;

    // `MakeWire(Edges, rangdeb, forward, W); W.Closed(true)` (`cxx:1020-1022`).
    make_wire_from(&edges, rangdeb, forward, w)?;
    Ok(())
}
