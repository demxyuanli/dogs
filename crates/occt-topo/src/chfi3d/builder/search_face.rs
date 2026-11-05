//! `ChFi3d_Builder::SearchFace` — the face the band leaves the guide through
//! at an extremity.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:1523-1627`
//! (OCCT 8.0.0).
//!
//! The three branches of the member function are all kept verbatim, including
//! the two that never report success:
//!
//! - a common point that is not on an arc answers `false` *before* the
//!   out-parameter is nullified (`cxx:1529-1532`), so a caller that reuses its
//!   `FVoi` keeps the previous value,
//! - a vertex carrying its own vector asks [`search_face_on_v`] for one or two
//!   faces (`cxx:1538-1544`),
//! - a vertex without a vector walks the guideline tangent and the `G1`
//!   neighbours (`cxx:1545-1620`) but falls through to `return false` at
//!   `cxx:1626` — `FVoi` may be set and the answer is still `false`. That is
//!   OCCT's control flow, not a transcription slip: the caller only uses the
//!   face to test `FVoi.IsNull()` (`_6.cxx:2192`, `:2324`),
//! - an arc answers [`is_g1`] on `Pc.Arc()` (`cxx:1624`).
//!
//! Divergences that are not heuristics:
//! - `Spine->D1` needs a mutable guideline in the port
//!   (`SpineHandle` is a lock), so the write guard is taken for the call only,
//! - a null `TopoDS_Face` / `TopoDS_Edge` has no counterpart in the port: the
//!   `FRef` parameter is an `Option` and a missing arc answers `false`.

use super::ctx::ChFi3dCtx;
use super::ChFi3dBuilder;

use crate::chfi_ds::{CommonPoint, SpineHandle};
use crate::shape::{Edge, Face};
use crate::topo_tools;

use crate::chfi3d::builder_2_helpers::{is_g1, is_input, search_face_on_v};

impl ChFi3dBuilder {
    /// `SearchFace(Spine, Pc, FRef, FVoi)` (`ChFi3d_Builder_2.cxx:1523-1627`).
    ///
    /// The body reads only `myVEMap` / `myEFMap`, so it lives on
    /// [`ChFi3dCtx`] and this is the plain-base forwarding the callers use.
    pub fn search_face(
        &self,
        spine: Option<&SpineHandle>,
        pc: &CommonPoint,
        f_ref: Option<&Face>,
        f_voi: &mut Option<Face>,
    ) -> bool {
        self.ctx().search_face(spine, pc, f_ref, f_voi)
    }
}

impl<'a> ChFi3dCtx<'a> {
    /// `SearchFace(Spine, Pc, FRef, FVoi)` (`ChFi3d_Builder_2.cxx:1523-1627`),
    /// as the subclass virtuals call it (they hold a [`ChFi3dCtx`], not the
    /// base object).
    ///
    /// `f_voi` is the out-parameter `FVoi`; it is set on every path that reached
    /// `FVoi.Nullify()` (`cxx:1533`), even when the answer is `false`.
    #[allow(unused_assignments)]
    pub fn search_face(
        &self,
        spine: Option<&SpineHandle>,
        pc: &CommonPoint,
        f_ref: Option<&Face>,
        f_voi: &mut Option<Face>,
    ) -> bool {
        // `if (!Pc.IsOnArc()) return false;` (`cxx:1529-1532`) — before the
        // nullification of `FVoi`.
        if !pc.is_on_arc() {
            return false;
        }
        // `FVoi.Nullify()` (`cxx:1533`).
        *f_voi = None;

        if pc.is_vertex() {
            if pc.has_vector() {
                // `nb_faces = SearchFaceOnV(...); return (nb_faces > 0);`
                // (`cxx:1540-1543`).
                let mut fbis = None;
                let nb_faces = search_face_on_v(
                    pc,
                    f_ref,
                    &self.ve_map,
                    &self.ef_map,
                    f_voi,
                    &mut fbis,
                );
                return nb_faces > 0;
            }

            // Processing using the spine (`cxx:1546-1620`).
            let Some(spine) = spine else {
                // `Spine->D1` below needs a guideline; OCCT dereferences the
                // handle.
                return false;
            };
            let Some(vertex) = pc.vertex() else {
                return false;
            };
            // `Spine->D1(Pc.Parameter(), Point, VecSpine)` (`cxx:1550`).
            let d1 = {
                let Ok(mut guard) = spine.write() else {
                    return false;
                };
                guard.d1(pc.parameter())
            };
            let Some((_point, vec_spine)) = d1 else {
                return false;
            };

            // `FindFace = IsInput(VecSpine, Pc.Vertex(), FRef); if (FindFace) {
            // VecSpine.Reverse(); FindFace = IsInput(...); }` (`cxx:1553-1558`).
            let mut find_face = is_input(&vec_spine, vertex, f_ref);
            if find_face {
                find_face = is_input(&vec_spine.reversed(), vertex, f_ref);
            }
            // `if (FindFace) { FVoi = FRef; return true; }` (`cxx:1560-1564`).
            if find_face {
                *f_voi = f_ref.cloned();
                return true;
            }

            // `for (ItE.Initialize(myVEMap(Pc.Vertex())); ItE.More() &&
            // (!FindFace); ItE.Next())` (`cxx:1569-1619`).
            for value in self.ve_map.get_from_key(vertex) {
                if find_face {
                    break;
                }
                let Some(e) = Edge::wrap(value.clone()) else {
                    continue;
                };
                let mut trouve = f_ref.is_some_and(|f_ref| {
                    self.ef_map
                        .get_from_key(&e.0)
                        .iter()
                        .any(|f| topo_tools::is_same(f, &f_ref.0))
                });
                if trouve {
                    find_face = is_g1(&self.ef_map, &e, f_ref, f_voi);
                }
                if !find_face {
                    continue;
                }
                find_face = false;
                // `if (Spine.IsNull()) return false;` (`cxx:1589-1593`): the null
                // guideline already answered `false` at the `Spine->D1` above.
                // It is checked if the selected face actually possesses edges of
                // the spine containing the vertex on its front (`cxx:1595-1617`).
                let Ok(guard) = spine.read() else {
                    return false;
                };
                trouve = false;
                for ie in 1..=guard.nb_edges() {
                    let Some(edge) = guard.edges(ie).cloned().map(Edge) else {
                        continue;
                    };
                    let (vf, vl) = topo_tools::edge_vertices(&edge);
                    let touches = vf
                        .as_ref()
                        .is_some_and(|v| topo_tools::is_same(&v.0, vertex))
                        || vl
                            .as_ref()
                            .is_some_and(|v| topo_tools::is_same(&v.0, vertex));
                    if !touches {
                        continue;
                    }
                    if self.ef_map.get_from_key(&edge.0).iter().any(|f| {
                        f_voi
                            .as_ref()
                            .is_some_and(|f_voi| topo_tools::is_same(f, &f_voi.0))
                    }) {
                        trouve = true;
                        break;
                    }
                }
                find_face = trouve;
            }
        } else {
            // `return IsG1(myEFMap, Pc.Arc(), FRef, FVoi)` (`cxx:1624`).
            let Some(arc) = pc.arc() else {
                return false;
            };
            let Some(arc) = Edge::wrap(arc.clone()) else {
                return false;
            };
            return is_g1(&self.ef_map, &arc, f_ref, f_voi);
        }
        // `return false;` (`cxx:1626`).
        false
    }
}
