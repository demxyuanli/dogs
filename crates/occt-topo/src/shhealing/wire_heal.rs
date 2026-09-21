use super::prelude::*;
use super::*;


/// Whether every wire of the shape is geometrically closed after
/// [`heal_shape`] runs with the given tolerance (used for both the welding and
/// the closing tolerance, and as the minimum-edge length).
pub fn wire_is_closed_after_heal(shape: &TopoShape, tol: f64) -> bool {
    let (healed, _) = heal_shape(shape, tol, tol);
    wires_of(&healed)
        .iter()
        .all(|w| wire_geometrically_closed_with_tol(w, tol))
}
/// Relocate a vertex to `to` and update its incident edges' geometry
/// (`ShapeFix_Vertex`).
///
/// The moved vertex is replaced by a fresh `Vertex` instance registered at the
/// target point (the input shape is untouched). Every edge incident to it is
/// rebuilt as a straight segment between the moved endpoint and its other
/// endpoint, so the edge curves pass through the new location. All other
/// vertices and edges are preserved. Errors when `vertex` is not found in
/// `shape`.
pub fn move_vertex(shape: &TopoShape, vertex: &Vertex, to: &GpPnt) -> Result<TopoShape, String> {
    let target_ptr = ptr(&vertex.0);
    let verts = vertices_of(shape);
    if !verts.iter().any(|v| ptr(&v.0) == target_ptr) {
        return Err("move_vertex: vertex not found in shape".into());
    }

    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();

    // The moved vertex becomes a new instance at the target point.
    let mut vertex_map: HashMap<usize, Vertex> = HashMap::new();
    for v in &verts {
        let vptr = ptr(&v.0);
        if vptr == target_ptr {
            let nv = builder.make_vertex(*to, reg.vertex_tolerance(v));
            vertex_map.insert(vptr, nv);
        } else {
            vertex_map.insert(vptr, v.clone());
        }
    }

    // Rebuild every incident edge so its curve runs through the new point.
    let edges = edges_of(shape);
    let mut edge_map: HashMap<usize, Edge> = HashMap::new();
    for e in &edges {
        let eptr = ptr(&e.0);
        let (a, b) = edge_vertices(e);
        match (a, b) {
            (Some(va), Some(vb)) => {
                let ia_ptr = ptr(&va.0);
                let ib_ptr = ptr(&vb.0);
                let incident = ia_ptr == target_ptr || ib_ptr == target_ptr;
                if !incident {
                    edge_map.insert(eptr, e.clone());
                    continue;
                }
                let ca = vertex_map.get(&ia_ptr).unwrap().clone();
                let cb = vertex_map.get(&ib_ptr).unwrap().clone();
                let pa = reg.vertex_point(&ca.0);
                let pb = reg.vertex_point(&cb.0);
                let ne = make_edge_with_vertices(&builder, &pa, &pb, &ca, &cb);
                edge_map.insert(eptr, ne);
            }
            _ => {
                edge_map.insert(eptr, e.clone());
            }
        }
    }

    let healed = rebuild_shape(
        shape,
        &HealCtx::new(builder, vertex_map, edge_map, HashSet::new(), HashMap::new()),
    );
    Ok(healed)
}
