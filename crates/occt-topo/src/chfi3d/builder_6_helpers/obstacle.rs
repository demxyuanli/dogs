//! `IsVois` / `IsObst` / `nbedconnex` — the obstacle test of a fillet
//! extremity against the vertex/edge map of the shape.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:107-205`
//! (OCCT 8.0.0, file-static helpers of `ChFi3d_Builder`).
//!
//! `IsObst(CP, Vref, VEMap)` answers `true` when the arc of `CP` is **not**
//! connected to `Vref` by a chain of edges that all meet vertex `Vref`, walking
//! the vertex/edge map up to `prof = 4` "wide" vertices
//! (`nbedconnex(L) > 2`). `IsVois` is that walk:
//!
//! ```text
//! if (prof > profmax) return false;          // depth budget
//! if (DONE.Contains(E)) return false;        // each edge once
//! TopExp::Vertices(E, V1, V2);
//! if (Vref.IsSame(V1) || Vref.IsSame(V2)) return true;
//! DONE.Add(E);
//! i1 = nbedconnex(VEMap(V1));
//! for (curE : VEMap(V1)) if (i1 <= 2 ? IsVois(curE, ..., prof)   : IsVois(curE, ..., prof + 1)) return true;
//! for (curE : VEMap(V2)) if (i1 <= 2 ? IsVois(curE, ..., prof)   : IsVois(curE, ..., prof + 1)) return true;
//! return false;
//! ```
//!
//! Note that OCCT reuses `i1` — computed from `V1`'s list — for the `V2` loop
//! too (`cxx:180`); the port keeps that, it is not a transcription slip.
//!
//! Divergences that are not heuristics:
//! - `ChFiDS_Map::operator()(S)` is [`ChFiDSMap::get_from_key`]; OCCT raises
//!   `Standard_NoSuchObject` for an absent key, the port answers an empty
//!   slice, so an isolated vertex contributes no neighbour instead of throwing,
//! - `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` is a
//!   `Vec<TopoShape>` scanned with [`topo_tools::is_same`], the same identity
//!   `IsSame` uses,
//! - `TopExp::Vertices(E, V1, V2)` on an edge with fewer than two vertices
//!   leaves null handles in OCCT, which `Vref.IsSame(...)` answers `false` for;
//!   the port's [`topo_tools::edge_vertices`] answers `None`, handled the same
//!   way,
//! - `int prof` / `profmax` are `i32`.

use crate::chfi_ds::{ChFiDSMap, CommonPoint};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools;

/// `nbedconnex(L)` (`cxx:107-130`): the number of distinct shapes of `l`, with
/// `IsSame` as the identity.
pub fn nbedconnex(l: &[TopoShape]) -> usize {
    let mut nb = 0usize;
    for (i, curs) in l.iter().enumerate() {
        let dejavu = l[..i]
            .iter()
            .any(|earlier| topo_tools::is_same(curs, earlier));
        if !dejavu {
            nb += 1;
        }
    }
    nb
}

/// `IsVois(E, Vref, VEMap, DONE, prof, profmax)` (`cxx:132-193`).
///
/// `done` is OCCT's `DONE` set, grown in place; seed it with an empty vector.
pub fn is_vois(
    e: &TopoShape,
    vref: &TopoShape,
    ve_map: &ChFiDSMap,
    done: &mut Vec<TopoShape>,
    prof: i32,
    profmax: i32,
) -> bool {
    if prof > profmax {
        return false;
    }
    if done.iter().any(|d| topo_tools::is_same(d, e)) {
        return false;
    }
    // `TopoDS::Edge(E)` (`cxx:159`, `:179`); the caller only ever passes edges.
    let Some(edge) = Edge::wrap(e.clone()) else {
        return false;
    };
    let (v1, v2) = topo_tools::edge_vertices(&edge);
    // `if (Vref.IsSame(V1) || Vref.IsSame(V2)) return true;` (`cxx:149-152`).
    let same_ref = |v: &Option<crate::shape::Vertex>| match v {
        Some(v) => topo_tools::is_same(&v.0, vref),
        None => false,
    };
    if same_ref(&v1) || same_ref(&v2) {
        return true;
    }
    // `DONE.Add(E)` (`cxx:153`).
    done.push(e.clone());

    // `const NCollection_List<TopoDS_Shape>& L1 = VEMap(V1);`
    // `int i1 = nbedconnex(L1);` (`cxx:154-155`).
    let l1: &[TopoShape] = match &v1 {
        Some(v) => ve_map.get_from_key(&v.0),
        None => &[],
    };
    let i1 = nbedconnex(l1);
    for cur_e in l1 {
        if i1 <= 2 {
            if is_vois(cur_e, vref, ve_map, done, prof, profmax) {
                return true;
            }
        } else if is_vois(cur_e, vref, ve_map, done, prof + 1, profmax) {
            return true;
        }
    }

    // `const NCollection_List<TopoDS_Shape>& L2 = VEMap(V2);` (`cxx:172`).
    let l2: &[TopoShape] = match &v2 {
        Some(v) => ve_map.get_from_key(&v.0),
        None => &[],
    };
    for cur_e in l2 {
        // `cxx:180` reuses `i1`, as OCCT does.
        if i1 <= 2 {
            if is_vois(cur_e, vref, ve_map, done, prof, profmax) {
                return true;
            }
        } else if is_vois(cur_e, vref, ve_map, done, prof + 1, profmax) {
            return true;
        }
    }
    false
}

/// `IsObst(CP, Vref, VEMap)` (`cxx:195-205`): `false` unless `CP` lies on an
/// arc, in which case it is `!IsVois(arc, Vref, VEMap, {}, 0, 4)`.
pub fn is_obst(cp: &CommonPoint, vref: &TopoShape, ve_map: &ChFiDSMap) -> bool {
    if !cp.is_on_arc() {
        return false;
    }
    // `const TopoDS_Edge& E = CP.Arc();` (`cxx:201`).
    let Some(arc) = cp.arc() else {
        return false;
    };
    let mut done: Vec<TopoShape> = Vec::new();
    // `int prof = 4;` (`cxx:203`).
    let prof = 4;
    !is_vois(arc, vref, ve_map, &mut done, 0, prof)
}
