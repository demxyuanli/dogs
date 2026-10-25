//! Edge counting around a vertex: `ChFi3d_ChercheBordsLibres`,
//! `ChFi3d_NbNotDegeneratedEdges`, `ChFi3d_NbSharpEdges`,
//! `ChFi3d_NumberOfEdges` and `ChFi3d_NumberOfSharpEdges`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:5557-5708`.

use occt_core::kernel::geomabs::Shape;

use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_helpers::faces::conex_faces;
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::map::ChFiDSMap;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::is_same;

/// `ChFi3d_ChercheBordsLibres(myVEMap, V1, bordlibre, edgelibre1,
/// edgelibre2)` (`cxx:5557-5611`).
///
/// A free border edge is a non-degenerated entry of the vertex→edges map that
/// occurs only once there (`cxx:5570-5585`); the map lists a junction edge
/// twice, so a `nboccur` of 1 means the edge really has one face side. The
/// second pass looks for a second such edge different from the first
/// (`cxx:5588-5611`).
///
/// Returns `(bordlibre, edgelibre1, edgelibre2)`; the two edges are `None`
/// whenever `bordlibre` is `false`.
pub fn cherche_bords_libres(
    ve_map: &ChFiDSMap,
    v1: &Vertex,
) -> (bool, Option<Edge>, Option<Edge>) {
    let list: &[TopoShape] = ve_map.get_from_key(&v1.0);
    let mut edge_libre1: Option<Edge> = None;
    let mut bord_libre = false;
    for value in list {
        let Some(cur) = Edge::wrap(value.clone()) else {
            continue;
        };
        let mut nboccur = 0usize;
        if !BRepTool::is_degenerated(&cur) {
            for other in list {
                if is_same(other, &cur.0) {
                    nboccur += 1;
                }
            }
        }
        if nboccur == 1 {
            edge_libre1 = Some(cur);
            bord_libre = true;
        }
        if bord_libre {
            break;
        }
    }
    if !bord_libre {
        return (false, None, None);
    }
    let first = edge_libre1.expect("bord_libre is set with it");
    let mut edge_libre2: Option<Edge> = None;
    bord_libre = false;
    for value in list {
        let Some(cur) = Edge::wrap(value.clone()) else {
            continue;
        };
        let mut nboccur = 0usize;
        if !BRepTool::is_degenerated(&cur) && !is_same(&cur.0, &first.0) {
            for other in list {
                if is_same(other, &cur.0) {
                    nboccur += 1;
                }
            }
        }
        if nboccur == 1 {
            edge_libre2 = Some(cur);
            bord_libre = true;
        }
        if bord_libre {
            break;
        }
    }
    (bord_libre, Some(first), edge_libre2)
}

/// `ChFi3d_NbNotDegeneratedEdges(Vtx, VEMap)` (`cxx:5619-5631`): the number of
/// non-degenerated entries of `VEMap(Vtx)`, junction edges counted twice.
pub fn nb_not_degenerated_edges(vtx: &Vertex, ve_map: &ChFiDSMap) -> i32 {
    let list = ve_map.get_from_key(&vtx.0);
    let mut nba = list.len() as i32;
    for value in list {
        if let Some(cur) = Edge::wrap(value.clone()) {
            if BRepTool::is_degenerated(&cur) {
                nba -= 1;
            }
        }
    }
    nba
}

/// `ChFi3d_NbSharpEdges(Vtx, VEMap, EFMap)` (`cxx:5639-5661`): the number of
/// entries of `VEMap(Vtx)` that are neither degenerated nor `G2`-tangent to
/// their opposite face.
pub fn nb_sharp_edges(
    vtx: &Vertex,
    ve_map: &ChFiDSMap,
    ef_map: &ChFiDSMap,
) -> Result<i32, String> {
    let list = ve_map.get_from_key(&vtx.0);
    let mut nba = list.len() as i32;
    for value in list {
        let Some(cur) = Edge::wrap(value.clone()) else {
            continue;
        };
        if BRepTool::is_degenerated(&cur) {
            nba -= 1;
        } else {
            let (f1, f2) = conex_faces(&cur, ef_map);
            if let Some(f2) = f2 {
                let f1 = f1.expect("conex_faces returns F1 with F2");
                if is_tangent_faces(&cur, &f1, &f2, Shape::G2)? {
                    nba -= 1;
                }
            }
        }
    }
    Ok(nba)
}

/// `ChFi3d_NumberOfEdges(Vtx, VEMap)` (`cxx:5668-5685`): the number of edges
/// arriving at the vertex, degenerated edges excluded and junction edges
/// collapsed.
///
/// The integer arithmetic is OCCT's: `(nba - 2) / 2 + 2` with a free border
/// (the two free ends count as one edge each) and `nba / 2` otherwise, both
/// truncating toward zero.
pub fn number_of_edges(vtx: &Vertex, ve_map: &ChFiDSMap) -> i32 {
    let nba = nb_not_degenerated_edges(vtx, ve_map);
    let (bord_libre, _, _) = cherche_bords_libres(ve_map, vtx);
    if bord_libre {
        (nba - 2) / 2 + 2
    } else {
        nba / 2
    }
}

/// `ChFi3d_NumberOfSharpEdges(Vtx, VEMap, EFmap)` (`cxx:5691-5708`): the same
/// fold as [`number_of_edges`] over the sharp-edge count.
pub fn number_of_sharp_edges(
    vtx: &Vertex,
    ve_map: &ChFiDSMap,
    ef_map: &ChFiDSMap,
) -> Result<i32, String> {
    let nba = nb_sharp_edges(vtx, ve_map, ef_map)?;
    let (bord_libre, _, _) = cherche_bords_libres(ve_map, vtx);
    Ok(if bord_libre {
        (nba - 2) / 2 + 2
    } else {
        nba / 2
    })
}
