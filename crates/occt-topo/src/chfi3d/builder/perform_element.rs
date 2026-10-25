//! `ChFi3d_Builder::PerformElement`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_1.cxx:887-1179`
//! (OCCT 8.0.0): starting from the first edge of a spine, walk the downstream
//! extremity and then the upstream one, appending every edge whose opposing
//! faces are tangent, and finally record the two `ChFiDS_State`s.
//!
//! The file statics the walk uses live in
//! [`super::perform_helpers`](crate::chfi3d::builder::perform_helpers):
//! `ReorderFaces`, `MakeOffsetEdge`, `TangentOnVertex`.
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - `theFirstFace` is `Option<&Face>` because a null `TopoDS_Face` has no
//!   counterpart in the port,
//! - `TopExp::Vertices` / `TopExp::CommonVertex` raise `Standard_NoSuchObject`
//!   on an edge without vertices and answer `false` when two edges share none;
//!   the port answers `Err` for both, which is the point where OCCT would keep
//!   walking with null vertex handles,
//! - `ChFi3d_conexfaces` may leave one of its two faces null; the downstream /
//!   upstream walks hand both to `ReorderFaces` and `MakeOffsetEdge`, so the
//!   port answers `Err` there,
//! - `Spine->SetEdges` / `PutInFirst` / `SetOffsetEdges` / `PutInFirstOffset`
//!   go through the shared-mutable [`SpineHandle`],
//! - the `Wl = Wf` store of both walks (`:1035`, `:1141`) is dead in OCCT too:
//!   `Wl` is recomputed at the top of the next iteration, so it is not ported.

use std::f64::consts::PI;

use occt_core::kernel::geomabs::Shape;

use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_helpers::conex_faces;
use crate::chfi3d::connect_type::define_connect_type;
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::{SpineHandle, State};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

use super::perform_helpers::{make_offset_edge, reorder_faces, tangent_on_vertex};
use super::ChFi3dBuilder;

/// `TopExp::CommonVertex(E1, E2, V)` (`TopExp.cxx:323-340`): the first endpoint
/// of `E1` that is an endpoint of `E2`, `None` for OCCT's `return false`.
fn common_vertex(e1: &Edge, e2: &Edge) -> Option<TopoShape> {
    let (f1, l1) = edge_vertices(e1);
    let (f2, l2) = edge_vertices(e2);
    let (f1, l1) = (f1?, l1?);
    let (f2, l2) = (f2?, l2?);
    if is_same(&f1.0, &f2.0) || is_same(&f1.0, &l2.0) {
        return Some(f1.0);
    }
    if is_same(&l1.0, &f2.0) || is_same(&l1.0, &l2.0) {
        return Some(l1.0);
    }
    None
}

impl ChFi3dBuilder {
    /// `PerformElement(Spine, Offset, theFirstFace)` (`cxx:887-1179`).
    ///
    /// `Ok(false)` covers OCCT's three early `return false` arms: a degenerated
    /// first edge, a first edge without two opposing faces, and a first edge
    /// whose faces are already tangent.
    ///
    /// `Err` carries the failures OCCT would raise (`Standard_OutOfRange` on an
    /// empty spine, `Standard_NoSuchObject` from `TopExp::Vertices`, the
    /// `ChFi3d::IsTangentFaces` / `DefineConnectType` exceptions).
    pub fn perform_element(
        &mut self,
        spine: &SpineHandle,
        offset: f64,
        the_first_face: Option<&Face>,
    ) -> Result<bool, String> {
        let ta = self.angular;
        let first_edge = {
            let guard = spine
                .read()
                .map_err(|_| "PerformElement: poisoned spine".to_string())?;
            guard
                .edges(1)
                .cloned()
                .ok_or_else(|| "PerformElement: ChFiDS_Spine::Edges(1) out of range".to_string())?
        };
        let mut ec = Edge::wrap(first_edge)
            .ok_or_else(|| "PerformElement: first spine shape is not an edge".to_string())?;
        if BRepTool::is_degenerated(&ec) {
            return Ok(false);
        }
        // it is checked if the edge is a cut edge
        let (ff1, ff2) = conex_faces(&ec, &self.my_ef_map);
        let (Some(ff1), Some(ff2)) = (ff1, ff2) else {
            return Ok(false);
        };
        if is_tangent_faces(&ec, &ff1, &ff2, Shape::G1)? {
            return Ok(false);
        }

        let mut ff1 = ff1;
        let mut ff2 = ff2;
        let mut first_face = ff1.clone();
        if let Some(the_first_face) = the_first_face {
            if is_same(&ff2.0, &the_first_face.0) {
                std::mem::swap(&mut ff1, &mut ff2);
                first_face = ff1.clone();
            }
        }
        self.bind_edge_first_face(&ec.0, &first_face.0);

        // Define concavity
        let type_of_concavity = define_connect_type(&ec, &ff1, &ff2, 1.0e-5, true)?;
        {
            let mut guard = spine
                .write()
                .map_err(|_| "PerformElement: poisoned spine".to_string())?;
            guard.set_type_of_concavity(type_of_concavity);
        }

        if offset > 0.0 {
            // `BRepAdaptor_Surface Sb1(ff1, ToRestrict)` with `ToRestrict =
            // Offset > 0` (`cxx:934-935`).
            let sb1 = SurfaceAdaptor::of_face_restricted(&ff1);
            let sb2 = SurfaceAdaptor::of_face_restricted(&ff2);
            if let Some(mut offset_edge) = make_offset_edge(&ec, offset, &sb1, &sb2) {
                offset_edge.0.set_orientation(ec.0.orientation());
                let mut guard = spine
                    .write()
                    .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                guard.set_offset_edges(&offset_edge.0);
            }
        }

        let cur_orientation = ec.0.orientation();
        let (v_start, lv_ec) = edge_vertices(&ec);
        let (Some(v_start), Some(mut lv_ec)) = (v_start, lv_ec) else {
            return Err("PerformElement: TopExp::Vertices on the first edge".to_string());
        };
        let mut fini = false;
        let mut cur_st = State::Closed;
        if is_same(&v_start.0, &lv_ec.0) {
            // case if only one edge is closed
            let cec = CurveAdaptor::of_edge(&ec);
            let wl = BRepTool::parameter_of_vertex(&ec.0, &v_start.0);
            let (_, v1) = cec
                .d1(wl)
                .ok_or_else(|| "PerformElement: CEc.D1".to_string())?;
            let wl = BRepTool::parameter_of_vertex(&ec.0, &lv_ec.0);
            let (_, v2) = cec
                .d1(wl)
                .ok_or_else(|| "PerformElement: CEc.D1".to_string())?;
            let is_face_tangency = self.face_tangency(&ec, &ec, &v_start)?;
            if v1.is_parallel_ang(&v2, ta) || is_face_tangency {
                if is_face_tangency {
                    cur_st = State::Closed;
                } else {
                    cur_st = State::BreakPoint;
                }
            } else {
                cur_st = State::BreakPoint;
            }
            let mut guard = spine
                .write()
                .map_err(|_| "PerformElement: poisoned spine".to_string())?;
            guard.set_last_status(cur_st);
            guard.set_first_status(cur_st);
        } else {
            // Downstream progression
            let mut cur_orientation = cur_orientation;
            while !fini {
                cur_st = State::FreeBoundary;
                let wl = BRepTool::parameter_of_vertex(&ec.0, &lv_ec.0);
                let degene_on_ec = tangent_on_vertex(&lv_ec.0, &ec, &self.my_ef_map, ta)?;
                let cec = CurveAdaptor::of_edge(&ec);
                let (_, v1) = cec
                    .d1(wl)
                    .ok_or_else(|| "PerformElement: CEc.D1".to_string())?;
                let nb = {
                    let guard = spine
                        .read()
                        .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                    guard.nb_edges()
                };

                for value in self.my_ve_map.get_from_key(&lv_ec.0).to_vec() {
                    // `TopoDS::Edge(It.Value())` would raise `TypeMismatch`; the
                    // VE map is filled from `TopAbs_EDGE` ancestors only.
                    let Some(ev) = Edge::wrap(value) else {
                        continue;
                    };
                    if is_same(&ev.0, &ec.0) || BRepTool::is_degenerated(&ev) {
                        continue;
                    }
                    let (fv_ev, lv_ev) = edge_vertices(&ev);
                    let (Some(mut fv_ev), Some(mut lv_ev)) = (fv_ev, lv_ev) else {
                        return Err(
                            "PerformElement: TopExp::Vertices on a vertex edge".to_string()
                        );
                    };
                    let or1;
                    if is_same(&lv_ec.0, &lv_ev.0) {
                        std::mem::swap(&mut fv_ev, &mut lv_ev);
                        or1 = crate::abs::Orientation::Reversed;
                    } else {
                        or1 = crate::abs::Orientation::Forward;
                    }

                    let wf = BRepTool::parameter_of_vertex(&ev.0, &fv_ev.0);
                    let cev = CurveAdaptor::of_edge(&ev);
                    let (_, v2) = cev
                        .d1(wf)
                        .ok_or_else(|| "PerformElement: CEv.D1".to_string())?;
                    let av1v2 = v1.angle(&v2);
                    let rev = or1 != cur_orientation;
                    let mut on_ajoute = false;
                    if self.face_tangency(&ec, &ev, &Vertex(fv_ev.0.clone()))? {
                        // there is no need of tolerance to make a decision
                        // (PRO9486) the regularity is enough. However, the
                        // absence of turn-back is checked (PRO9810)
                        on_ajoute = (!rev && av1v2 < PI / 2.0) || (rev && av1v2 > PI / 2.0);
                        // mate attention to the single case (cf CTS21610_1)
                        if on_ajoute
                            && (degene_on_ec
                                || tangent_on_vertex(&lv_ec.0, &ev, &self.my_ef_map, ta)?)
                        {
                            on_ajoute =
                                (!rev && av1v2 < ta) || (rev && (PI - av1v2) < ta);
                        }
                    }
                    if on_ajoute {
                        fini = false; // If this can be useful (Cf PRO14713)
                        let Some(common) = common_vertex(&ec, &ev) else {
                            return Err(
                                "PerformElement: TopExp::CommonVertex failed".to_string()
                            );
                        };
                        let prev_edge = ec.clone();
                        ec = ev;
                        ec.0.set_orientation(or1);
                        lv_ec = lv_ev.clone();
                        {
                            let mut guard = spine.write().map_err(|_| {
                                "PerformElement: poisoned spine".to_string()
                            })?;
                            guard.set_edges(&ec.0);
                        }
                        let (cur_f1, cur_f2) = conex_faces(&ec, &self.my_ef_map);
                        let (Some(mut cur_f1), Some(mut cur_f2)) = (cur_f1, cur_f2) else {
                            return Err(
                                "PerformElement: ChFi3d_conexfaces on a walked edge".to_string()
                            );
                        };
                        reorder_faces(
                            &mut cur_f1,
                            &mut cur_f2,
                            &first_face,
                            &prev_edge,
                            &common,
                            &self.my_ef_map,
                        );
                        self.bind_edge_first_face(&ec.0, &cur_f1.0);
                        if offset > 0.0 {
                            // `BRepAdaptor_Surface CurSb1(CurF1)` — no restriction
                            // here, unlike the two at `cxx:934-935`.
                            let cur_sb1 = SurfaceAdaptor::of_face(&cur_f1);
                            let cur_sb2 = SurfaceAdaptor::of_face(&cur_f2);
                            if let Some(mut an_offset_edge) =
                                make_offset_edge(&ec, offset, &cur_sb1, &cur_sb2)
                            {
                                an_offset_edge.0.set_orientation(or1);
                                let mut guard = spine.write().map_err(|_| {
                                    "PerformElement: poisoned spine".to_string()
                                })?;
                                guard.set_offset_edges(&an_offset_edge.0);
                            }
                        }
                        first_face = cur_f1;
                        cur_orientation = or1;
                        if is_same(&v_start.0, &lv_ev.0) {
                            let first = {
                                let guard = spine.read().map_err(|_| {
                                    "PerformElement: poisoned spine".to_string()
                                })?;
                                guard.edges(1).cloned()
                            };
                            let Some(first) = first else {
                                return Err(
                                    "PerformElement: ChFiDS_Spine::Edges(1) out of range"
                                        .to_string(),
                                );
                            };
                            if self.face_tangency(&ec, &Edge(first), &Vertex(lv_ev.0.clone()))? {
                                cur_st = State::Closed;
                                fini = true;
                            } else {
                                cur_st = State::BreakPoint;
                                fini = true;
                            }
                        }
                        break;
                    } else {
                        let nbface = self.my_ef_map.get_from_key(&ev.0).len();
                        if nbface > 1 {
                            cur_st = State::BreakPoint;
                        }
                        fini = (!rev && av1v2 < ta) || (rev && (PI - av1v2) < ta);
                    }
                }
                let nb_after = {
                    let guard = spine
                        .read()
                        .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                    guard.nb_edges()
                };
                fini = fini || (nb == nb_after);
            }
            {
                let mut guard = spine
                    .write()
                    .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                guard.set_last_status(cur_st);
            }
            if cur_st == State::Closed {
                let mut guard = spine
                    .write()
                    .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                guard.set_first_status(cur_st);
            } else {
                // Upstream progression
                fini = false;
                let first_spine_edge = {
                    let guard = spine
                        .read()
                        .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                    guard.edges(1).cloned()
                }
                .ok_or_else(|| {
                    "PerformElement: ChFiDS_Spine::Edges(1) out of range".to_string()
                })?;
                let mut ec = Edge::wrap(first_spine_edge).ok_or_else(|| {
                    "PerformElement: first spine shape is not an edge".to_string()
                })?;
                let mut first_face = {
                    let bound = self
                        .edge_first_face(&ec.0)
                        .cloned()
                        .ok_or_else(|| {
                            "PerformElement: myEdgeFirstFace has no entry for the first edge"
                                .to_string()
                        })?;
                    Face::wrap(bound).ok_or_else(|| {
                        "PerformElement: myEdgeFirstFace is not a face".to_string()
                    })?
                };
                let mut cur_orientation = ec.0.orientation();
                let mut fv_ec = v_start.clone();
                while !fini {
                    cur_st = State::FreeBoundary;
                    let wl = BRepTool::parameter_of_vertex(&ec.0, &fv_ec.0);
                    let degene_on_ec = tangent_on_vertex(&fv_ec.0, &ec, &self.my_ef_map, ta)?;
                    let cec = CurveAdaptor::of_edge(&ec);
                    let (_, v1) = cec
                        .d1(wl)
                        .ok_or_else(|| "PerformElement: CEc.D1".to_string())?;
                    let nb = {
                        let guard = spine.read().map_err(|_| {
                            "PerformElement: poisoned spine".to_string()
                        })?;
                        guard.nb_edges()
                    };

                    for value in self.my_ve_map.get_from_key(&fv_ec.0).to_vec() {
                        let Some(ev) = Edge::wrap(value) else {
                            continue;
                        };
                        if is_same(&ev.0, &ec.0) || BRepTool::is_degenerated(&ev) {
                            continue;
                        }
                        let (fv_ev, lv_ev) = edge_vertices(&ev);
                        let (Some(mut fv_ev), Some(mut lv_ev)) = (fv_ev, lv_ev) else {
                            return Err(
                                "PerformElement: TopExp::Vertices on a vertex edge".to_string()
                            );
                        };
                        let or1;
                        if is_same(&fv_ec.0, &fv_ev.0) {
                            std::mem::swap(&mut fv_ev, &mut lv_ev);
                            or1 = crate::abs::Orientation::Reversed;
                        } else {
                            or1 = crate::abs::Orientation::Forward;
                        }
                        let wf = BRepTool::parameter_of_vertex(&ev.0, &lv_ev.0);
                        let cev = CurveAdaptor::of_edge(&ev);
                        let (_, v2) = cev
                            .d1(wf)
                            .ok_or_else(|| "PerformElement: CEv.D1".to_string())?;
                        let av1v2 = v1.angle(&v2);
                        let rev = or1 != cur_orientation;
                        let mut on_ajoute = false;
                        if self.face_tangency(&ec, &ev, &Vertex(lv_ev.0.clone()))? {
                            on_ajoute =
                                (!rev && av1v2 < PI / 2.0) || (rev && av1v2 > PI / 2.0);
                            if on_ajoute
                                && (degene_on_ec
                                    || tangent_on_vertex(&fv_ec.0, &ev, &self.my_ef_map, ta)?)
                            {
                                on_ajoute =
                                    (!rev && av1v2 < ta) || (rev && (PI - av1v2) < ta);
                            }
                        }
                        if on_ajoute {
                            let Some(common) = common_vertex(&ec, &ev) else {
                                return Err(
                                    "PerformElement: TopExp::CommonVertex failed".to_string()
                                );
                            };
                            let prev_edge = ec.clone();
                            ec = ev;
                            ec.0.set_orientation(or1);
                            fv_ec = fv_ev;
                            {
                                let mut guard = spine.write().map_err(|_| {
                                    "PerformElement: poisoned spine".to_string()
                                })?;
                                guard.put_in_first(&ec.0);
                            }
                            let (cur_f1, cur_f2) = conex_faces(&ec, &self.my_ef_map);
                            let (Some(mut cur_f1), Some(mut cur_f2)) = (cur_f1, cur_f2)
                            else {
                                return Err(
                                    "PerformElement: ChFi3d_conexfaces on a walked edge"
                                        .to_string(),
                                );
                            };
                            reorder_faces(
                                &mut cur_f1,
                                &mut cur_f2,
                                &first_face,
                                &prev_edge,
                                &common,
                                &self.my_ef_map,
                            );
                            self.bind_edge_first_face(&ec.0, &cur_f1.0);
                            if offset > 0.0 {
                                let cur_sb1 = SurfaceAdaptor::of_face(&cur_f1);
                                let cur_sb2 = SurfaceAdaptor::of_face(&cur_f2);
                                if let Some(mut an_offset_edge) =
                                    make_offset_edge(&ec, offset, &cur_sb1, &cur_sb2)
                                {
                                    an_offset_edge.0.set_orientation(or1);
                                    let mut guard = spine.write().map_err(|_| {
                                        "PerformElement: poisoned spine".to_string()
                                    })?;
                                    guard.put_in_first_offset(&an_offset_edge.0);
                                }
                            }
                            first_face = cur_f1;
                            cur_orientation = or1;
                            break;
                        } else {
                            let nbface = self.my_ef_map.get_from_key(&ev.0).len();
                            if nbface > 1 {
                                cur_st = State::BreakPoint;
                            }
                            fini = (!rev && av1v2 < ta) || (rev && (PI - av1v2) < ta);
                        }
                    }
                    let nb_after = {
                        let guard = spine.read().map_err(|_| {
                            "PerformElement: poisoned spine".to_string()
                        })?;
                        guard.nb_edges()
                    };
                    fini = fini || (nb == nb_after);
                }
                let mut guard = spine
                    .write()
                    .map_err(|_| "PerformElement: poisoned spine".to_string())?;
                guard.set_first_status(cur_st);
            }
        }
        Ok(true)
    }
}
