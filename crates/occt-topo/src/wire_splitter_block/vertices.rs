use super::prelude::*;
use super::SmartMap;
use super::angles::*;

pub(super) fn dump_map(face: &Face, map: &SmartMap, vert_closed: &HashMap<usize, bool>) {
    for slot in 0..map.verts.len() {
        let v = &map.verts[slot];
        let p = BRepTool::vertex_point(&Vertex(v.clone()));
        let vk = GeometryRegistry::shape_key(v);
        eprintln!(
            "[sb] V key={} p=({},{},{}) closed={:?} n={}",
            sk(v),
            p.x(),
            p.y(),
            p.z(),
            vert_closed.get(&vk).copied().unwrap_or(false) as u8,
            map.infos[slot].len()
        );
        for ei in &map.infos[slot] {
            let mut vv = v.clone();
            vv.set_orientation(if ei.is_in {
                Orientation::Reversed
            } else {
                Orientation::Forward
            });
            let tv = vertex_parameter_on_face(&vv, &ei.edge);
            let mut pv = vv.clone();
            pv.set_orientation(Orientation::Forward);
            let uv = coord2d(&pv, &ei.edge, face);
            eprintln!(
                "[sb]   E key={} in={} ins={} tv={} uv=({},{}) ang={}",
                sk(&ei.edge),
                ei.is_in as u8,
                ei.is_inside as u8,
                tv,
                uv.x(),
                uv.y(),
                ei.angle
            );
        }
    }
}

pub(super) fn sk(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

pub(super) fn trace() -> bool {
    std::env::var("ZZ_SB").is_ok()
}

pub(super) fn edge_closed_on_face(e: &Edge, _face: &Face, face_key: usize) -> bool {
    if BRepTool::is_degenerated(e) {
        return true;
    }
    GeometryRegistry::global().edge_pcurves(&e.0, face_key).len() >= 2
}

pub(super) fn stored_vertices(e: &Edge) -> Vec<TopoShape> {
    e.0.tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
}

/// The edge's stored vertices with the `(FORWARD, REVERSED)` orientation pair
/// OCCT guarantees (`BRep_Builder::Add`); the port also builds edges whose two
/// children are stored Forward, so the last child is normalised to Reversed
/// there. No composition with the edge's own orientation.
pub(super) fn stored_vertices_oriented(e: &Edge) -> Vec<TopoShape> {
    let stored = stored_vertices(e);
    let mut oris: Vec<Orientation> = stored.iter().map(|v| v.orientation()).collect();
    if oris.len() >= 2
        && oris[0] == Orientation::Forward
        && oris[1] == Orientation::Forward
    {
        oris[1] = Orientation::Reversed;
    }
    stored
        .into_iter()
        .zip(oris)
        .map(|(mut v, o)| {
            v.set_orientation(o);
            v
        })
        .collect()
}

/// Edge vertices with `TopoDS_Iterator(cumOri=true)` orientation.
///
/// When both children are stored Forward (`make_edge_segment` in this port),
/// the last child is treated as Reversed, matching `BRepLib_MakeEdge`.
pub(super) fn oriented_vertices(e: &Edge) -> Vec<TopoShape> {
    stored_vertices_oriented(e)
        .into_iter()
        .map(|mut v| {
            v.set_orientation(Orientation::compose(e.0.orientation(), v.orientation()));
            v
        })
        .collect()
}
