//! `ChFi3d_Builder::PerformExtremity` and the `ConexFaces` member it drives.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_1.cxx:730-885`
//! and `ChFi3d_Builder_2.cxx:859-880` (OCCT 8.0.0).
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - `Spine->SetFirstStatus` / `SetLastStatus` become
//!   [`Spine::set_status`](crate::chfi_ds::Spine::set_status) on the
//!   shared-mutable handle,
//! - `ConexFaces` hands out two `BRepAdaptor_Surface` through null handles; the
//!   port answers `Option` for the case where `ChFi3d_conexfaces` finds no face
//!   pair, which OCCT would follow into a null-shape dereference,
//! - `E[3]` is an array of `Option<Edge>` because a `TopoDS_Edge` default
//!   constructs to the null shape, which the port has no value for.

use occt_core::kernel::geomabs::Shape;

use crate::abs::Orientation;
use crate::chfi3d::builder_helpers::{conex_faces, edge_state};
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::{SpineHandle, State};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::perform_helpers::tangent_extremity;
use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `ConexFaces(Spine, IEdge, HS1, HS2)` (`ChFi3d_Builder_2.cxx:859-880`):
    /// the pair of faces of `Spine->Edges(IEdge)`, ordered so that the first one
    /// is `myEdgeFirstFace(anEdge)`.
    pub fn cone_x_faces(
        &self,
        spine: &SpineHandle,
        i_edge: usize,
    ) -> Option<(SurfaceAdaptor, SurfaceAdaptor)> {
        let an_edge = {
            let guard = spine.read().ok()?;
            guard.edges(i_edge)?.clone()
        };
        let edge = Edge::wrap(an_edge.clone())?;
        let (ff1, ff2) = conex_faces(&edge, &self.my_ef_map);
        let mut ff1 = ff1?;
        let mut ff2 = ff2?;
        if let Some(first_face) = self.edge_first_face(&an_edge) {
            let first_face = Face::wrap(first_face.clone())?;
            if is_same(&ff2.0, &first_face.0) {
                std::mem::swap(&mut ff1, &mut ff2);
            }
        }
        Some((
            SurfaceAdaptor::of_face(&ff1),
            SurfaceAdaptor::of_face(&ff2),
        ))
    }

    /// `PerformExtremity(Spine)` (`ChFi3d_Builder_1.cxx:730-885`): after
    /// `PerformElement` may have left an extremity at `ChFiDS_BreakPoint`, the
    /// state is refined from the concavities of the faces around the vertex.
    ///
    /// `Err` carries the `ChFi3d::IsTangentFaces` / `ChFi3d_EdgeState` failure
    /// OCCT would throw.
    pub fn perform_extremity(&mut self, spine: &SpineHandle) -> Result<(), String> {
        let mut nb_g1_connections = 0usize;

        for ii in 1..=2 {
            // `ChFiDS_State sst = ChFiDS_OnSame; int iedge = 0; TopoDS_Vertex v;`
            // followed by the two-branch fill of `cxx:757-770`.
            let (mut sst, iedge, v) = {
                let guard = spine
                    .read()
                    .map_err(|_| "PerformExtremity: poisoned spine".to_string())?;
                if ii == 1 {
                    (guard.first_status(), 1, guard.first_vertex())
                } else {
                    (guard.last_status(), guard.nb_edges(), guard.last_vertex())
                }
            };
            // `Spine->Edges(iedge)` raises `Standard_OutOfRange` on a guideline
            // without edges; the port skips the extremity.
            let e0 = {
                let guard = spine
                    .read()
                    .map_err(|_| "PerformExtremity: poisoned spine".to_string())?;
                match guard.edges(iedge) {
                    Some(e) => e.clone(),
                    None => continue,
                }
            };
            let Some(v) = v else {
                continue;
            };

            // Before all it is checked if the tangency is not dead.
            if let Some((hs1, hs2)) = self.cone_x_faces(spine, iedge) {
                if tangent_extremity(&v, &Edge(e0.clone()), &hs1, &hs2, self.angular)? {
                    if let Ok(mut guard) = spine.write() {
                        guard.set_tangency_extremity(true, ii == 1);
                    }
                }
            }

            if sst == State::BreakPoint {
                let mut a_loc_nb_g1_connections = 0usize;
                let mut sommet_pourri = false;
                let mut edges_of_v: Vec<TopoShape> = Vec::new();
                let mut edges: Vec<TopoShape> = Vec::new();
                let mut e: [Option<Edge>; 3] = [Some(Edge(e0.clone())), None, None];
                edges.push(e0.clone());
                edges_of_v.push(e0.clone());
                let mut ind_of_e = 0usize;
                let value_edges: Vec<TopoShape> = self.my_ve_map.get_from_key(&v).to_vec();
                for value in value_edges {
                    let Some(an_edge) = Edge::wrap(value) else {
                        continue;
                    };
                    if crate::brep_tool::BRepTool::is_degenerated(&an_edge) {
                        continue;
                    }
                    let (f1, f2) = conex_faces(&an_edge, &self.my_ef_map);
                    if let (Some(f1), Some(f2)) = (&f1, &f2) {
                        // smooth edge
                        if is_tangent_faces(&an_edge, f1, f2, Shape::G2)? {
                            if !is_same(&f1.0, &f2.0) {
                                nb_g1_connections += 1;
                                a_loc_nb_g1_connections += 1;
                            }
                            continue;
                        }
                    }

                    if !edges.iter().any(|x| is_same(x, &an_edge.0)) {
                        edges.push(an_edge.0.clone());
                        edges_of_v.push(an_edge.0.clone());
                        if ind_of_e < 2 {
                            ind_of_e += 1;
                            e[ind_of_e] = Some(an_edge);
                        }
                    } else {
                        // edge is closed - two ends of the edge in the vertex
                        let (v1, v2) = edge_vertices(&an_edge);
                        let closed = match (&v1, &v2) {
                            (Some(v1), Some(v2)) => is_same(&v1.0, &v2.0),
                            _ => false,
                        };
                        if closed {
                            let mut an_ind = edges_of_v
                                .iter()
                                .position(|x| is_same(x, &an_edge.0));
                            if an_ind.is_none() {
                                an_ind = edges_of_v.iter().position(|x| {
                                    is_same(x, &an_edge.0.oriented(Orientation::Reversed))
                                });
                            }
                            if let Some(an_ind) = an_ind {
                                let reversed = Edge(
                                    edges_of_v[an_ind]
                                        .clone()
                                        .oriented(Orientation::Reversed),
                                );
                                // `EdgesOfV.Add(anEdge)` after `anEdge.Reverse()`:
                                // `TopTools_ShapeMapHasher` compares with
                                // `IsSame`, so the reversed handle is already in
                                // the map and `Add` answers false — in OCCT too,
                                // which keeps this arm from ever filling `E`.
                                if !edges_of_v.iter().any(|x| is_same(x, &reversed.0)) {
                                    edges_of_v.push(reversed.0.clone());
                                    if ind_of_e < 2 {
                                        ind_of_e += 1;
                                        e[ind_of_e] = Some(reversed);
                                    }
                                }
                            }
                        }
                    }
                }

                if edges_of_v.len() != 3 {
                    sommet_pourri = true;
                }

                if !sommet_pourri && a_loc_nb_g1_connections < 4 {
                    // `E[0..2]` are all set when three shapes were collected.
                    let (Some(e0), Some(e1), Some(e2)) =
                        (e[0].as_ref(), e[1].as_ref(), e[2].as_ref())
                    else {
                        return Err(
                            "PerformExtremity: ChFi3d_EdgeState with a null edge".to_string()
                        );
                    };
                    sst = edge_state([e0, e1, e2], &self.my_ef_map)?;
                }
                if let Ok(mut guard) = spine.write() {
                    guard.set_status(sst, ii == 1);
                }
            }
        }

        let is_periodic = {
            let guard = spine
                .read()
                .map_err(|_| "PerformExtremity: poisoned spine".to_string())?;
            guard.is_periodic()
        };
        if !is_periodic {
            let distinct_faces = |list: &[TopoShape]| {
                let mut nbf: i64 = 0;
                for (jf, cur) in list.iter().enumerate() {
                    let mut kf = 0usize;
                    for earlier in list.iter().take(jf) {
                        kf += 1;
                        if is_same(cur, earlier) {
                            break;
                        }
                    }
                    if kf == jf {
                        nbf += 1;
                    }
                }
                nbf
            };
            let first_vertex = {
                let guard = spine
                    .read()
                    .map_err(|_| "PerformExtremity: poisoned spine".to_string())?;
                guard.first_vertex()
            };
            if let Some(first_vertex) = first_vertex {
                let nbf = distinct_faces(self.my_vf_map.get_from_key(&first_vertex))
                    - nb_g1_connections as i64;
                if nbf > 3 {
                    if let Ok(mut guard) = spine.write() {
                        guard.set_first_status(State::BreakPoint);
                    }
                }
            }
            let last_vertex = {
                let guard = spine
                    .read()
                    .map_err(|_| "PerformExtremity: poisoned spine".to_string())?;
                guard.last_vertex()
            };
            if let Some(last_vertex) = last_vertex {
                let nbf = distinct_faces(self.my_vf_map.get_from_key(&last_vertex))
                    - nb_g1_connections as i64;
                if nbf > 3 {
                    if let Ok(mut guard) = spine.write() {
                        guard.set_last_status(State::BreakPoint);
                    }
                }
            }
        }

        Ok(())
    }
}
