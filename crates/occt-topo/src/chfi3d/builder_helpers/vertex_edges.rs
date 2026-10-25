//! Vertex/edge neighbourhood searches: `ChFi3d_cherche_face1`,
//! `ChFi3d_cherche_element`, `ChFi3d_cherche_edge`, `ChFi3d_cherche_vertex`,
//! `ChFi3d_nbface`, `ChFi3d_edge_common_faces` and `ChFi3d_AngleEdge`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:5330-5549`.

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, edges_of, is_same};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// `ChFi3d_cherche_face1(map, F1, F)` (`cxx:5334-5353`): the first face of the
/// list that is not `F1`.
///
/// The list is the edge→faces entry, so it holds the two faces of an edge
/// occurrence (the same face twice for a seam).
pub fn cherche_face1(map: &[TopoShape], f1: &Face) -> Result<Face, String> {
    for value in map {
        let Some(cur) = Face::wrap(value.clone()) else {
            continue;
        };
        if !is_same(&cur.0, &f1.0) {
            return Ok(cur);
        }
    }
    Err("Standard_ConstructionError: ChFi3d_cherche_face1".to_string())
}

/// `ChFi3d_cherche_element(V, E1, F1, E, Vtx)` (`cxx:5361-5405`): the edge of
/// `F1` other than `E1` that carries `V`, together with its other vertex.
///
/// Only edges with exactly two vertices are considered (`cxx:5380`).
pub fn cherche_element(
    v: &Vertex,
    e1: &Edge,
    f1: &Face,
) -> Result<(Edge, Vertex), String> {
    for cur in edges_of(&f1.0) {
        if is_same(&cur.0, &e1.0) {
            continue;
        }
        let (v1, v2) = edge_vertices(&cur);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            continue;
        };
        // `TopExp::MapShapes(Ecur, TopAbs_VERTEX, MapV)` counts *distinct*
        // vertices; `edge_vertices` returns the two stored occurrences, which
        // differ for a closed edge. `MapV.Extent() == 2` therefore needs the
        // two vertices to differ as shapes.
        if is_same(&v1.0, &v2.0) {
            continue;
        }
        if is_same(&v1.0, &v.0) {
            return Ok((cur, v2));
        } else if is_same(&v2.0, &v.0) {
            return Ok((cur, v1));
        }
    }
    Err("Standard_ConstructionError: ChFi3d_cherche_element".to_string())
}

/// `ChFi3d_cherche_edge(V, E1, F1, E, Vtx)` (`cxx:5410-5462`): like
/// [`cherche_element`] but skipping every edge of the `e1` array.
pub fn cherche_edge(
    v: &Vertex,
    e1: &[TopoShape],
    f1: &Face,
) -> Result<(Edge, Vertex), String> {
    for cur in edges_of(&f1.0) {
        let same = e1.iter().any(|excluded| is_same(&cur.0, excluded));
        if same {
            continue;
        }
        let (v1, v2) = edge_vertices(&cur);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            continue;
        };
        if is_same(&v1.0, &v2.0) {
            continue;
        }
        if is_same(&v1.0, &v.0) {
            return Ok((cur, v2));
        } else if is_same(&v2.0, &v.0) {
            return Ok((cur, v1));
        }
    }
    Err("Standard_ConstructionError: ChFi3d_cherche_edge".to_string())
}

/// `ChFi3d_cherche_vertex(E1, E2, vertex, trouve)` (`cxx:5716-5744`): the
/// vertex the two edges share, `None` when they do not share one.
pub fn cherche_vertex(e1: &Edge, e2: &Edge) -> Option<Vertex> {
    let verts1 = crate::topo_tools::vertices_of(&e1.0);
    let verts2 = crate::topo_tools::vertices_of(&e2.0);
    for v1 in &verts1 {
        for v2 in &verts2 {
            if is_same(&v2.0, &v1.0) {
                return Some(v1.clone());
            }
        }
    }
    None
}

/// `ChFi3d_nbface(mapVF)` (`cxx:5468-5490`): how many *distinct* faces the
/// vertex map lists.
///
/// The OCCT loop counts an entry when no earlier entry is the same shape, i.e.
/// first occurrences only — this is not `TopTools_MapOfShape::Extent` because
/// the order and the `IsSame` test are the list's.
pub fn nbface(map_vf: &[TopoShape]) -> usize {
    let mut nface = 0usize;
    for (j, cur) in map_vf.iter().enumerate() {
        let mut kf = 0usize;
        for earlier in map_vf.iter().take(j) {
            if is_same(cur, earlier) {
                break;
            }
            kf += 1;
        }
        if kf == j {
            nface += 1;
        }
    }
    nface
}

/// `ChFi3d_edge_common_faces(mapEF, F1, F2)` (`cxx:5498-5521`): the two faces
/// of an edge, `F2 == F1` when the edge has only one.
///
/// `None` when the list is empty — OCCT dereferences the iterator's first value
/// unconditionally (`cxx:5505`).
pub fn edge_common_faces(map_ef: &[TopoShape]) -> Option<(Face, Face)> {
    let f1 = Face::wrap(map_ef.first()?.clone())?;
    for value in map_ef {
        let Some(f) = Face::wrap(value.clone()) else {
            continue;
        };
        if !is_same(&f.0, &f1.0) {
            return Some((f1, f));
        }
    }
    Some((f1.clone(), f1))
}

/// `ChFi3d_AngleEdge(Vtx, E1, E2)` (`cxx:5527-5549`): the unsigned angle
/// between the two edge tangents, both oriented away from `Vtx`.
pub fn angle_edge(vtx: &Vertex, e1: &Edge, e2: &Edge) -> Result<f64, String> {
    let par1 = BRepTool::parameter_of_vertex(&e1.0, &vtx.0);
    let par2 = BRepTool::parameter_of_vertex(&e2.0, &vtx.0);
    let c1 = CurveAdaptor::of_edge(e1);
    let c2 = CurveAdaptor::of_edge(e2);
    let (_, mut dir1) = c1
        .d1(par1)
        .ok_or_else(|| "ChFi3d_AngleEdge: no curve on E1".to_string())?;
    let (_, mut dir2) = c2
        .d1(par2)
        .ok_or_else(|| "ChFi3d_AngleEdge: no curve on E2".to_string())?;
    let (first1, _) = edge_vertices(e1);
    if !first1.is_some_and(|f| is_same(&f.0, &vtx.0)) {
        dir1.reverse();
    }
    let (first2, _) = edge_vertices(e2);
    if !first2.is_some_and(|f| is_same(&f.0, &vtx.0)) {
        dir2.reverse();
    }
    Ok(dir1.angle(&dir2).abs())
}
