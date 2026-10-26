//! `PerformIntersectionAtEnd` — the faces intersected at the end,
//! `ChFi3d_Builder_C1.cxx:2281-2580`.
//!
//! Ported blocks, in OCCT order:
//! - `cxx:2281-2446` — the `nface == 3` walk: `findonf1` / `findonf2`, the
//!   `ChFiDS_OnSame` sub-case (`ChFi3d_cherche_element` x2 +
//!   `ChFi3d_EdgeState`), the two face chains with their `cherche_edge1`
//!   limiting edges, and the `isOnSame2` close;
//! - `cxx:2448-2573` — the "more than three edges" search over
//!   `ChFi3d_ChercheBordsLibres`' free edges, `cherche_face` /
//!   `ChFi3d_cherche_edge`, the `IntersectMoreCorner(Index)` dispatch and the
//!   `nbarete == 4` / `nbarete == 5` collapses;
//! - `cxx:2575-2580` — the `onecorner` answer.
//!
//! Divergences that are not heuristics:
//! - `Edge[nn]` / `Face[nn]` are C arrays of default-constructed null shapes in
//!   OCCT; the port keeps `Option<...>` and reports a null read as `Err` with
//!   the source line (see [`super::corner_faces`]).
//! - `ChFi3d_cherche_vertex(E1, E2, Vcom, trouve)` sets `trouve = false` on
//!   entry (`ChFi3d_Builder_0.cxx:5723`) and writes `Vcom` only when both edges
//!   share a vertex, so `Option<Vertex>` carries the same information. `ang1`
//!   is **not** reset between the two `nbarete == 4` probes (`cxx:2514`,
//!   `:2538`), so the second probe may test the angle the first one computed;
//!   that carry is kept.
//! - `ChFi3d_cherche_edge(Vtx, tabedg, Face[nb - 1], Edge[nb], V)` receives
//!   `tabedg` whose unset entries are null shapes; `IsSame` never matches a
//!   null shape, so passing the set entries only is equivalent.
//! - `ChFi3d_AngleEdge` raising (null curve on either edge) and
//!   `ChFi3d::IsTangentFaces` raising are `Err` here, like the port's other
//!   `ChFi3d_C1` callers.

use occt_core::kernel::geomabs::Shape;

use crate::chfi_ds::State;
use crate::chfi3d::builder_c1::{cherche_edge1, cherche_face};
use crate::chfi3d::builder_helpers::{
    angle_edge, cherche_edge, cherche_element, cherche_face1, cherche_vertex, edge_state,
};
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::is_same;

use super::super::ChFi3dBuilder;
use super::corner_faces::{edge_at, edge_carries_vertex, edge_has_two_other_vertices, face_at};
use super::{end_points, EndIntersection, EndIntersectionExit, NN};

/// `cxx:2281-2580`: the faces intersected at the end and the edges limiting
/// those intersections.
pub(super) fn collect_faces(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
) -> Result<Option<EndIntersectionExit>, String> {
    if st.one_corner {
        // `nbface = 1; Face[0] = face; Edge[1] = Edge[2];` (`cxx:2575-2580`).
        st.nbface = 1;
        st.faces[0] = st.face.clone();
        st.edges[1] = st.edges[2].clone();
        return Ok(None);
    }
    // `if (nface == 3)` (`cxx:2304`).
    if st.nface == 3 {
        three_face_walk(b, st)?;
        Ok(None)
    } else {
        many_edge_walk(b, st)
    }
}

/// `cxx:2304-2446`: the three-face case.
fn three_face_walk(b: &mut ChFi3dBuilder, st: &mut EndIntersection) -> Result<(), String> {
    let vtx = Vertex(st.vtx.clone());
    // `if (CV1.IsVertex()) findonf1 = true; if (CV2.IsVertex()) findonf2 =
    //  true;` (`cxx:2306-2313`).
    let (cv1, cv2) = end_points(st)?;
    if cv1.is_vertex() {
        st.find_on_f1 = true;
    }
    if cv2.is_vertex() {
        st.find_on_f2 = true;
    }
    // `cxx:2314-2337`: an end edge whose two vertices both differ from `Vtx`.
    if !st.find_on_f1 {
        let e0 = edge_at(st, 0, 2317)?;
        if edge_has_two_other_vertices(&e0, &vtx) {
            st.find_on_f1 = true;
        }
    }
    if !st.find_on_f2 {
        let e2 = edge_at(st, 2, 2329)?;
        if edge_has_two_other_vertices(&e2, &vtx) {
            st.find_on_f2 = true;
        }
    }

    // `if (state == ChFiDS_OnSame)` (`cxx:2340-2365`).
    if st.state == State::OnSame {
        let f1 = need(st.f1.clone(), "F1")?;
        let f2 = need(st.f2.clone(), "F2")?;
        let edge_spine = spine_edge(st)?;
        // `ChFi3d_cherche_element(Vtx, EdgeSpine, F1, threeE[0], V2);
        //  ChFi3d_cherche_element(Vtx, EdgeSpine, F2, threeE[1], V2);
        //  threeE[2] = EdgeSpine;` (`cxx:2343-2345`).
        let three0 = cherche_element(&vtx, &edge_spine, &f1)?.0;
        let three1 = cherche_element(&vtx, &edge_spine, &f2)?.0;
        // `if (ChFi3d_EdgeState(threeE, myEFMap) == ChFiDS_OnSame)`
        // (`cxx:2346-2360`).
        if edge_state([&three0, &three1, &edge_spine], &b.my_ef_map)? == State::OnSame {
            st.is_on_same1 = true;
            st.nb = 1;
            st.edges[0] = Some(three0.clone());
            st.faces[0] = Some(cherche_face1(b.my_ef_map.get_from_key(&three0.0), &f1)?);
            if st.find_on_f2 {
                // `findonf1 = true; // not to look for Face[0] again`
                st.find_on_f1 = true;
            } else {
                st.edges[1] = Some(Edge(arc_of(&cv2, "CV2")?));
            }
        } else {
            st.is_on_same2 = true;
        }
    }

    // `if (findonf1 && !isOnSame1)` (`cxx:2371-2401`).
    if st.find_on_f1 && !st.is_on_same1 {
        let f1 = need(st.f1.clone(), "F1")?;
        let arc1 = Edge(arc_of(&cv1, "CV1")?);
        // `if (CV1.TransitionOnArc() == TopAbs_FORWARD) V1 = FirstVertex else
        //  V1 = LastVertex` (`cxx:2373-2380`).
        let (fwd, rev) = crate::topo_tools::edge_vertices(&arc1);
        let trans1 = cv1
            .transition_on_arc()
            .ok_or_else(|| "PerformIntersectionAtEnd: CV1 is not on an arc".to_string())?;
        let mut v1 = if trans1.is_forward() { fwd } else { rev }
            .ok_or_else(|| "PerformIntersectionAtEnd: CV1.Arc() has no such vertex".to_string())?;

        // `ChFi3d_cherche_face1(myEFMap(CV1.Arc()), F1, Face[0]); nb = 1;
        //  Ei = Edge[0];` (`cxx:2381-2383`).
        st.faces[0] = Some(cherche_face1(b.my_ef_map.get_from_key(&arc1.0), &f1)?);
        st.nb = 1;
        let mut ei = edge_at(st, 0, 2383)?;
        // `while (!V1.IsSame(Vtx))` (`cxx:2384-2396`).
        while !is_same(&v1.0, &vtx.0) {
            // `ChFi3d_cherche_element(V1, Ei, F1, E, V2); V1 = V2; Ei = E;`
            let (e, v2) = cherche_element(&v1, &ei, &f1)?;
            v1 = v2;
            ei = e;
            // `ChFi3d_cherche_face1(myEFMap(E), F1, Face[nb]);`
            let f = cherche_face1(b.my_ef_map.get_from_key(&ei.0), &f1)?;
            st.faces[st.nb] = Some(f.clone());
            // `cherche_edge1(Face[nb - 1], Face[nb], Edge[nb]);`
            let prev = face_at(st, st.nb - 1, 2390)?;
            st.edges[st.nb] = Some(cherche_edge1(&prev, &f)?);
            st.nb += 1;
            if st.nb >= NN {
                return Err("IntersectionAtEnd : the max number of faces reached".to_string());
            }
        }
        if !st.find_on_f2 {
            st.edges[st.nb] = Some(Edge(arc_of(&cv2, "CV2")?));
        }
    }

    // `if (findonf2 && !isOnSame2)` (`cxx:2402-2438`).
    if st.find_on_f2 && !st.is_on_same2 {
        let f1 = need(st.f1.clone(), "F1")?;
        let f2 = need(st.f2.clone(), "F2")?;
        // `if (!findonf1) nb = 1;` (`cxx:2404-2407`).
        if !st.find_on_f1 {
            st.nb = 1;
        }
        let mut v1 = vtx.clone();
        let arc2 = Edge(arc_of(&cv2, "CV2")?);
        // `if (CV2.TransitionOnArc() == TopAbs_FORWARD) Vfin =
        //  LastVertex(CV2.Arc()); else Vfin = FirstVertex(CV2.Arc());`
        // (`cxx:2409-2416`).
        let (fwd2, rev2) = crate::topo_tools::edge_vertices(&arc2);
        let trans2 = cv2
            .transition_on_arc()
            .ok_or_else(|| "PerformIntersectionAtEnd: CV2 is not on an arc".to_string())?;
        let vfin = if trans2.is_forward() { rev2 } else { fwd2 }
            .ok_or_else(|| "PerformIntersectionAtEnd: CV2.Arc() has no such vertex".to_string())?;
        // `if (!findonf1) ChFi3d_cherche_face1(myEFMap(CV1.Arc()), F1,
        //  Face[nb - 1]);` (`cxx:2417-2420`).
        if !st.find_on_f1 {
            let arc1 = arc_of(&cv1, "CV1")?;
            st.faces[st.nb - 1] = Some(cherche_face1(b.my_ef_map.get_from_key(&arc1), &f1)?);
        }
        // `ChFi3d_cherche_element(V1, EdgeSpine, F2, E, V2); Ei = E; V1 = V2;`
        // (`cxx:2421-2423`).
        let edge_spine = spine_edge(st)?;
        let (mut e, v2) = cherche_element(&v1, &edge_spine, &f2)?;
        let mut ei = e.clone();
        v1 = v2;
        // `while (!V1.IsSame(Vfin))` (`cxx:2424-2436`).
        while !is_same(&v1.0, &vfin.0) {
            let (e_next, v2) = cherche_element(&v1, &ei, &f2)?;
            e = e_next;
            ei = e.clone();
            v1 = v2;
            let f = cherche_face1(b.my_ef_map.get_from_key(&ei.0), &f2)?;
            st.faces[st.nb] = Some(f.clone());
            let prev = face_at(st, st.nb - 1, 2430)?;
            st.edges[st.nb] = Some(cherche_edge1(&prev, &f)?);
            st.nb += 1;
            if st.nb >= NN {
                return Err("IntersectionAtEnd : the max number of faces reached".to_string());
            }
        }
        st.edges[st.nb] = Some(arc2);
    }

    // `if (isOnSame2) { cherche_edge1(Face[nb - 1], F2, Edge[nb]);
    //  Face[nb] = F2; }` (`cxx:2439-2443`).
    if st.is_on_same2 {
        let f2 = need(st.f2.clone(), "F2")?;
        let prev = face_at(st, st.nb - 1, 2441)?;
        st.edges[st.nb] = Some(cherche_edge1(&prev, &f2)?);
        st.faces[st.nb] = Some(f2);
    }

    // `nbface = nb;` (`cxx:2445`).
    st.nbface = i32::try_from(st.nb)
        .map_err(|_| "PerformIntersectionAtEnd: nb does not fit in an int".to_string())?;
    Ok(())
}

/// `cxx:2448-2573`: "a top has more than three edges".
fn many_edge_walk(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
) -> Result<Option<EndIntersectionExit>, String> {
    let vtx = Vertex(st.vtx.clone());
    let (cv1, cv2) = end_points(st)?;
    let f1 = need(st.f1.clone(), "F1")?;
    let f2 = need(st.f2.clone(), "F2")?;
    let arc1 = Edge(arc_of(&cv1, "CV1")?);
    let arc2 = Edge(arc_of(&cv2, "CV2")?);

    // `for (ex.Init(CV1.Arc(), TopAbs_VERTEX); ...) if (Vtx.IsSame(...))
    //  possible1 = true;` and the same for `CV2` (`cxx:2457-2470`).
    let possible1 = edge_carries_vertex(&arc1, &vtx);
    let possible2 = edge_carries_vertex(&arc2, &vtx);

    // `if ((possible1 && possible2) || (!possible1 && !possible2) ||
    //  (nbarete > 4))` (`cxx:2471`).
    if !((possible1 && possible2) || (!possible1 && !possible2) || (st.nbarete > 4)) {
        // `IntersectMoreCorner(Index); return;` (`cxx:2508-2509`).
        return Ok(Some(EndIntersectionExit::IntersectMoreCorner));
    }

    // `trouve = false; while (!trouve) { nb++; if (nb >= nn) throw; if (nb !=
    //  1) F3 = Face[nb - 2]; Face[nb - 1] = F3; ... }` (`cxx:2454`,
    //  `:2473-2503`). `F3` starts at `F1` (`cxx:2148`).
    let mut trouve = false;
    let mut f3 = f1.clone();
    while !trouve {
        st.nb += 1;
        if st.nb >= NN {
            return Err("IntersectionAtEnd : the max number of faces reached".to_string());
        }
        if st.nb != 1 {
            f3 = face_at(st, st.nb - 2, 2482)?;
        }
        st.faces[st.nb - 1] = Some(f3.clone());
        // `if (CV1.Arc().IsSame(edgelibre1)) cherche_face(myVFMap(Vtx),
        //  edgelibre2, ...); else if (CV1.Arc().IsSame(edgelibre2))
        //  cherche_face(myVFMap(Vtx), edgelibre1, ...); else
        //  cherche_face(myVFMap(Vtx), Edge[nb - 1], ...)`
        //  (`cxx:2485-2496`).
        let probe = if st
            .edge_libre1
            .as_ref()
            .is_some_and(|e| is_same(&arc1.0, &e.0))
        {
            need(st.edge_libre2.clone(), "edgelibre2")?
        } else if st
            .edge_libre2
            .as_ref()
            .is_some_and(|e| is_same(&arc1.0, &e.0))
        {
            need(st.edge_libre1.clone(), "edgelibre1")?
        } else {
            edge_at(st, st.nb - 1, 2495)?
        };
        let f = cherche_face(b.my_vf_map.get_from_key(&vtx.0), &probe, &f1, &f2, &f3)?;
        st.faces[st.nb - 1] = Some(f.clone());
        // `ChFi3d_cherche_edge(Vtx, tabedg, Face[nb - 1], Edge[nb], V);
        //  tabedg.SetValue(nb, Edge[nb]);` (`cxx:2497-2498`).
        let e_nb = cherche_edge(&vtx, &st.tabedg, &f)?.0;
        st.edges[st.nb] = Some(e_nb.clone());
        st.tabedg.push(e_nb.0.clone());
        // `if (Edge[nb].IsSame(CV2.Arc())) trouve = true;` (`cxx:2499-2502`).
        if is_same(&e_nb.0, &arc2.0) {
            trouve = true;
        }
    }
    // `nbface = nb;` (`cxx:2504`).
    st.nbface = i32::try_from(st.nb)
        .map_err(|_| "PerformIntersectionAtEnd: nb does not fit in an int".to_string())?;

    if st.nbarete == 4 {
        // `cxx:2511-2549`.
        collapse_four(st, &vtx)?;
    } else if st.nbarete == 5 {
        // `cxx:2550-2572`.
        collapse_five(st, &f1)?;
    }
    Ok(None)
}

/// `if (nbarete == 4)` (`cxx:2511-2549`): two consecutive G1 edges collapse to
/// one intersection face.
fn collapse_four(st: &mut EndIntersection, vtx: &Vertex) -> Result<(), String> {
    // `double ang1 = 0.0; trouve = false; ChFi3d_cherche_vertex(Edge[0],
    //  Edge[1], Vcom, trouve); if (Vcom.IsSame(Vtx)) ang1 =
    //  ChFi3d_AngleEdge(Vtx, Edge[0], Edge[1]);` (`cxx:2514-2521`).
    let mut ang1 = 0.0_f64;
    let e0 = edge_at(st, 0, 2517)?;
    let e1 = edge_at(st, 1, 2517)?;
    if cherche_vertex(&e0, &e1).is_some_and(|vcom| is_same(&vcom.0, &vtx.0)) {
        ang1 = angle_edge(vtx, &e0, &e1)?;
    }
    // `if (abs(ang1 - M_PI) < 0.01)` (`cxx:2522-2530`).
    if (ang1 - std::f64::consts::PI).abs() < 0.01 {
        st.one_intersection1 = true;
        st.face_sau = face_at(st, 0, 2525)?.into();
        st.edge_sau = Some(e1);
        st.faces[0] = face_at(st, 1, 2527)?.into();
        st.edges[1] = edge_at(st, 2, 2528)?.into();
        st.nbface = 1;
    }

    // `if (!oneintersection1)` (`cxx:2532-2548`). `ang1` is **not** reset
    // before the second probe, matching the source.
    if !st.one_intersection1 {
        let e1 = edge_at(st, 1, 2535)?;
        let e2 = edge_at(st, 2, 2535)?;
        if cherche_vertex(&e1, &e2).is_some_and(|vcom| is_same(&vcom.0, &vtx.0)) {
            ang1 = angle_edge(vtx, &e1, &e2)?;
        }
        if (ang1 - std::f64::consts::PI).abs() < 0.01 {
            st.one_intersection2 = true;
            st.face_sau = face_at(st, 1, 2543)?.into();
            st.edge_sau = Some(e1);
            st.edges[1] = e2.into();
            st.nbface = 1;
        }
    }
    Ok(())
}

/// `else if (nbarete == 5)` (`cxx:2550-2572`, `pro15368`).
fn collapse_five(st: &mut EndIntersection, f1: &Face) -> Result<(), String> {
    let e0 = edge_at(st, 0, 2554)?;
    let e1 = edge_at(st, 1, 2555)?;
    let e2 = edge_at(st, 2, 2556)?;
    let face0 = face_at(st, 0, 2555)?;
    let face1 = face_at(st, 1, 2556)?;
    let face2 = face_at(st, 2, 2556)?;
    // `isTangent0 = IsTangentFaces(Edge[0], F1, Face[0]); isTangent1 =
    //  IsTangentFaces(Edge[1], Face[0], Face[1]); isTangent2 =
    //  IsTangentFaces(Edge[2], Face[1], Face[2]);` (`cxx:2554-2556`).
    let is_tangent0 = is_tangent_faces(&e0, f1, &face0, Shape::G1)?;
    let is_tangent1 = is_tangent_faces(&e1, &face0, &face1, Shape::G1)?;
    let is_tangent2 = is_tangent_faces(&e2, &face1, &face2, Shape::G1)?;
    // `if ((isTangent0 || isTangent2) && isTangent1)` (`cxx:2557-2571`).
    if (is_tangent0 || is_tangent2) && is_tangent1 {
        st.face_sau = Some(face0);
        st.edge_sau = Some(e0);
        st.nbface = 1;
        st.edges[1] = edge_at(st, 3, 2568)?.into();
        st.faces[0] = Some(face2);
        st.one_intersection1 = true;
    }
    Ok(())
}

/// `EdgeSpine` (`cxx:1901` / `:1907`), re-read from the spine.
fn spine_edge(st: &EndIntersection) -> Result<Edge, String> {
    let guard = st
        .spine
        .read()
        .map_err(|_| "PerformIntersectionAtEnd: poisoned spine".to_string())?;
    let e = if st.is_first {
        guard.edges(1)
    } else {
        guard.edges(st.nbedge)
    };
    let e = e
        .cloned()
        .ok_or_else(|| "PerformIntersectionAtEnd: spine has no such edge".to_string())?;
    Ok(Edge(e))
}

/// `CV.Arc()` or the null dereference OCCT would take.
fn arc_of(cv: &crate::chfi_ds::CommonPoint, what: &str) -> Result<crate::shape::TopoShape, String> {
    cv.arc()
        .cloned()
        .ok_or_else(|| format!("PerformIntersectionAtEnd: {what}.Arc() is null"))
}

/// A local OCCT reads unconditionally that the port keeps as `Option`.
fn need<T>(value: Option<T>, what: &str) -> Result<T, String> {
    value.ok_or_else(|| format!("PerformIntersectionAtEnd: {what} is null"))
}
