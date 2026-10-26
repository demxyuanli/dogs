//! `PerformIntersectionAtEnd` — the two end edges of the corner,
//! `ChFi3d_Builder_C1.cxx:2136-2280`.
//!
//! Ported blocks, in OCCT order:
//! - `cxx:2146-2160` — `F1` / `F2` re-read from the DS, the
//!   `couture || bordlibre` bump of `nface` and the `nbface` split;
//! - `cxx:2161-2165` — the `!CV1.IsOnArc() || !CV2.IsOnArc()` dispatch to
//!   `PerformMoreThreeCorner(Index, 1)`;
//! - `cxx:2167-2170` — the two end edges and their `tabedg` entries;
//! - `cxx:2175-2227` / `cxx:2229-2280` — the `CV1.IsVertex()` /
//!   `CV2.IsVertex()` re-selection of the end edge by the distance from the
//!   corner vertex, and the `SetArc` that goes with it.
//!
//! Divergences that are not heuristics:
//! - `Edge[nn]` is a C array of default-constructed null shapes in OCCT. The
//!   port keeps `Option<Edge>`; reading an entry the algorithm never wrote is
//!   the null-shape dereference OCCT would take, answered as `Err` with the
//!   source line. The literal `Edge[2]` of the `CV2.IsVertex()` block is kept
//!   as `2` on purpose: with `nface > 3` OCCT reads an entry that only the
//!   `nface == 3` path fills (`nbface == 2`), and that quirk is part of the
//!   translated control flow.
//! - `CV1` / `CV2` are `ChFiDS_CommonPoint&` into the `SurfData` in OCCT. The
//!   port takes them out as copies, mutates them here and writes them back
//!   with [`super::store_end_points`].

use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_c1::contain_e;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

use super::super::ChFi3dBuilder;
use super::{end_points, fd_of, EndIntersection, EndIntersectionExit, NN};

/// `Edge[i]` of the C array: `None` where OCCT holds a null shape.
pub(super) fn edge_at(st: &EndIntersection, i: usize, cxx: u32) -> Result<Edge, String> {
    st.edges
        .get(i)
        .and_then(|e| e.clone())
        .ok_or_else(|| format!("PerformIntersectionAtEnd: Edge[{i}] is null (cxx:{cxx})"))
}

/// `Face[i]` of the C array.
pub(super) fn face_at(st: &EndIntersection, i: usize, cxx: u32) -> Result<Face, String> {
    st.faces
        .get(i)
        .and_then(|f| f.clone())
        .ok_or_else(|| format!("PerformIntersectionAtEnd: Face[{i}] is null (cxx:{cxx})"))
}

/// `TopExp::MapShapes(E, TopAbs_VERTEX, MapV); MapV.Extent() == 2` followed by
/// "neither is `Vtx`" (`cxx:2316-2337`).
///
/// `MapV` is an `IndexedMap`, so a closed edge whose two occurrences are the
/// same vertex counts as extent 1.
pub(super) fn edge_has_two_other_vertices(e: &Edge, vtx: &Vertex) -> bool {
    let (v1, v2) = edge_vertices(e);
    match (v1, v2) {
        (Some(a), Some(b)) => {
            !is_same(&a.0, &b.0) && !is_same(&a.0, &vtx.0) && !is_same(&b.0, &vtx.0)
        }
        _ => false,
    }
}

/// `for (ex.Init(E, TopAbs_VERTEX); ex.More(); ex.Next()) if
/// (Vtx.IsSame(ex.Current()))` (`cxx:2457-2463`).
pub(super) fn edge_carries_vertex(e: &Edge, vtx: &Vertex) -> bool {
    let (v1, v2) = edge_vertices(e);
    v1.is_some_and(|v| is_same(&v.0, &vtx.0)) || v2.is_some_and(|v| is_same(&v.0, &vtx.0))
}

/// `cxx:2136-2170`.
pub(super) fn start(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
) -> Result<Option<EndIntersectionExit>, String> {
    // `F1 = TopoDS::Face(DStr.Shape(Fd->IndexOfS1())); F2 =
    //  TopoDS::Face(DStr.Shape(Fd->IndexOfS2())); F3 = F1;` (`cxx:2146-2148`).
    let fd = fd_of(st)?;
    let f1 = ds_face(b, fd.index_of_s1(), 1)?;
    let f2 = ds_face(b, fd.index_of_s2(), 2)?;
    st.f1 = Some(f1);
    st.f2 = Some(f2);

    // `if (couture || bordlibre) nface = nface + 1;` (`cxx:2149-2152`).
    if st.couture || st.bord_libre {
        st.nface += 1;
    }
    // `if (nface == 3) nbface = 2; else nbface = nface - 2;` (`cxx:2153-2160`).
    // `nface < 2` makes OCCT index `Edge[nbface]` with a negative value; the
    // port reports that instead of taking the out-of-bounds write.
    let nbface = if st.nface == 3 {
        2
    } else {
        i32::try_from(st.nface)
            .ok()
            .and_then(|n| n.checked_sub(2))
            .filter(|n| *n >= 0)
            .ok_or_else(|| {
                "PerformIntersectionAtEnd: nface < 2 makes Edge[nbface] out of range".to_string()
            })?
    };
    if nbface as usize >= NN {
        return Err("IntersectionAtEnd : the max number of faces reached".to_string());
    }
    st.nbface = nbface;

    // `if (!CV1.IsOnArc() || !CV2.IsOnArc()) { PerformMoreThreeCorner(Index, 1);
    //  return; }` (`cxx:2161-2165`).
    let (cv1, cv2) = end_points(st)?;
    if !cv1.is_on_arc() || !cv2.is_on_arc() {
        return Ok(Some(EndIntersectionExit::MoreThreeCorner(1)));
    }
    let arc1 = cv1
        .arc()
        .ok_or_else(|| "PerformIntersectionAtEnd: CV1.Arc() is null".to_string())?
        .clone();
    let arc2 = cv2
        .arc()
        .ok_or_else(|| "PerformIntersectionAtEnd: CV2.Arc() is null".to_string())?
        .clone();
    // `Edge[0] = CV1.Arc(); Edge[nbface] = CV2.Arc(); tabedg.SetValue(0, ...);
    //  tabedg.SetValue(nbface, ...);` (`cxx:2167-2170`).
    st.edges[0] = Some(Edge(arc1.clone()));
    st.edges[nbface as usize] = Some(Edge(arc2.clone()));
    st.tabedg = vec![arc1, arc2];
    Ok(None)
}

/// `DStr.Shape(Fd->IndexOfS1()/IndexOfS2())` cast to a face (`cxx:2146-2147`).
fn ds_face(b: &mut ChFi3dBuilder, index: i32, on_s: i32) -> Result<Face, String> {
    let shape = b.my_ds.change_ds().shape(index, false).ok_or_else(|| {
        format!("PerformIntersectionAtEnd: DStr.Shape(IndexOfS{on_s}) is null")
    })?;
    Face::wrap(shape)
        .ok_or_else(|| format!("PerformIntersectionAtEnd: IndexOfS{on_s} is not a face"))
}

/// `cxx:2175-2280`.
///
/// OCCT binds `ChFiDS_CommonPoint& CV1 = Fd->ChangeVertex(...)`, so the
/// `SetArc` calls mutate the `SurfData`.
pub(super) fn resolve_end_edges(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    cv1: &mut crate::chfi_ds::CommonPoint,
    cv2: &mut crate::chfi_ds::CommonPoint,
) -> Result<(), String> {
    // `if (CV1.IsVertex())` (`cxx:2176-2227`).
    if let Some(v) = cv1.vertex().and_then(|s| Vertex::wrap(s.clone())) {
        // `trouve = false; for (It3.Initialize(myVEMap(V)); It3.More() &&
        //  !trouve; It3.Next()) { E = ...; if (!E.IsSame(Edge[0]) &&
        //  containE(F1, E)) trouve = true; }` (`cxx:2178-2187`).
        let e0 = edge_at(st, 0, 2179)?;
        let f1 = st
            .f1
            .clone()
            .ok_or_else(|| "PerformIntersectionAtEnd: F1 is null".to_string())?;
        let e = first_other_edge(b, &v, &e0, &f1)
            .ok_or_else(|| "PerformIntersectionAtEnd: no edge of myVEMap(V) in F1".to_string())?;

        // `V1 = FirstVertex(Edge[0]); V2 = LastVertex(Edge[0]); if
        //  (V.IsSame(V1)) Vt = V2; else Vt = V1;
        //  dist1 = Pnt(Vt).Distance(Pnt(Vtx));` (`cxx:2188-2199`).
        let (v1, v2) = two_vertices(&e0, "Edge[0]")?;
        let dist1 = far_vertex_distance(&v, &v1, &v2, st);

        // `V3 = FirstVertex(E); V4 = LastVertex(E); if (V.IsSame(V3)) Vt = V4;
        //  else Vt = V3; dist2 = ...;` (`cxx:2200-2210`).
        let (v3, v4) = two_vertices(&e, "E")?;
        let dist2 = far_vertex_distance(&v, &v3, &v4, st);

        // `if (dist2 < dist1) { Edge[0] = E; ori = (V2.IsSame(V3) ||
        //  V1.IsSame(V4)) ? CV1.TransitionOnArc() :
        //  TopAbs::Reverse(CV1.TransitionOnArc()); par =
        //  BRep_Tool::Parameter(V, Edge[0]); CV1.SetArc(tol, Edge[0], par,
        //  ori); }` (`cxx:2211-2226`).
        if dist2 < dist1 {
            st.edges[0] = Some(e.clone());
            let trans = cv1
                .transition_on_arc()
                .ok_or_else(|| "PerformIntersectionAtEnd: CV1 is not on an arc".to_string())?;
            let ori = if is_same(&v2.0, &v3.0) || is_same(&v1.0, &v4.0) {
                trans
            } else {
                trans.reversed()
            };
            let par = BRepTool::parameter_of_vertex(&e.0, &v.0);
            let tol = cv1.tolerance();
            cv1.set_arc(tol, &e.0, par, ori);
        }
    }

    // `if (CV2.IsVertex())` (`cxx:2229-2280`). OCCT reaches `Edge[2]` here, not
    // `Edge[nbface]` (see the module note).
    if let Some(v) = cv2.vertex().and_then(|s| Vertex::wrap(s.clone())) {
        let e2 = edge_at(st, 2, 2236)?;
        let f2 = st
            .f2
            .clone()
            .ok_or_else(|| "PerformIntersectionAtEnd: F2 is null".to_string())?;
        let e = first_other_edge(b, &v, &e2, &f2)
            .ok_or_else(|| "PerformIntersectionAtEnd: no edge of myVEMap(V) in F2".to_string())?;

        let (v1, v2) = two_vertices(&e2, "Edge[2]")?;
        let dist1 = far_vertex_distance(&v, &v1, &v2, st);

        let (v3, v4) = two_vertices(&e, "E")?;
        let dist2 = far_vertex_distance(&v, &v3, &v4, st);

        if dist2 < dist1 {
            st.edges[2] = Some(e.clone());
            let trans = cv2
                .transition_on_arc()
                .ok_or_else(|| "PerformIntersectionAtEnd: CV2 is not on an arc".to_string())?;
            let ori = if is_same(&v2.0, &v3.0) || is_same(&v1.0, &v4.0) {
                trans
            } else {
                trans.reversed()
            };
            let par = BRepTool::parameter_of_vertex(&e.0, &v.0);
            let tol = cv2.tolerance();
            cv2.set_arc(tol, &e.0, par, ori);
        }
    }
    Ok(())
}

/// The `for (It3.Initialize(myVEMap(V)); ...) if (!E.IsSame(Edge) &&
/// containE(F1, E))` search of `cxx:2180-2187` / `cxx:2233-2240`.
fn first_other_edge(
    b: &ChFi3dBuilder,
    v: &Vertex,
    edge: &Edge,
    f: &Face,
) -> Option<Edge> {
    for value in b.my_ve_map.get_from_key(&v.0) {
        let Some(e) = Edge::wrap(value.clone()) else {
            continue;
        };
        if !is_same(&e.0, &edge.0) && contain_e(f, &e) {
            return Some(e);
        }
    }
    None
}

/// `TopExp::FirstVertex` / `LastVertex`, both required.
fn two_vertices(e: &Edge, what: &str) -> Result<(Vertex, Vertex), String> {
    let (v1, v2) = edge_vertices(e);
    let v1 = v1.ok_or_else(|| format!("PerformIntersectionAtEnd: {what} has no first vertex"))?;
    let v2 = v2.ok_or_else(|| format!("PerformIntersectionAtEnd: {what} has no last vertex"))?;
    Ok((v1, v2))
}

/// `Vt = (V.IsSame(A) ? B : A); dist = Pnt(Vt).Distance(Pnt(Vtx));`
/// (`cxx:2191-2199` and its `cxx:2202-2210` twin).
fn far_vertex_distance(v: &Vertex, a: &Vertex, b: &Vertex, st: &EndIntersection) -> f64 {
    let vt = if is_same(&v.0, &a.0) { b } else { a };
    BRepTool::vertex_point(vt).distance(&BRepTool::vertex_point(&Vertex(st.vtx.clone())))
}
