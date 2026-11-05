//! The topology stage of `PerformMoreThreeCorner` (`cxx:1237-1581`): the vertex
//! counts, the sewing edge, the edge/face ring parse and the `regul` flags.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1237-1581`
//! (OCCT 8.0.0).
//!
//! The ring is walked one `Enext` at a time: `ChFi3d_cherche_edge` finds the
//! edge of the current face that carries `V1`, the face on the other side of
//! that edge becomes the next current face, and so on until `nedge` entries are
//! filled. A free border (an edge of `ChFi3d_ChercheBordsLibres`) contributes
//! two entries, both carrying the current face.
//!
//! Translation notes:
//! - `ChFi3d_cherche_edge` / `ChFi3d_cherche_face1` write through an
//!   out-parameter that OCCT does **not** clear on failure, and only raise when
//!   the answer is still null afterwards (`ChFi3d_Builder_0.cxx:5350-5353`,
//!   `:5458-5461`). [`cherche_edge_into`] restores exactly that behaviour; the
//!   `Fnext` call site keeps its pre-set `Fcur` value when nothing is found,
//!   which is why the port ignores the `Result` there,
//! - the OCCT body keeps `Enext` across iterations, so a failed search in the
//!   second iteration and onwards silently reuses the previous edge. The port
//!   keeps the same variable across the loop for the same reason.

use std::f64::consts::PI;

use occt_core::kernel::geomabs::Shape;

use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi3d::builder_helpers::{
    angle_edge, cherche_bords_libres, cherche_edge, cherche_face1, couture_on_vertex,
    edge_common_faces, nb_not_degenerated_edges,
};
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::StripeHandle;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

use super::remove_surf_data::remove_surf_data;
use super::surface_normal::summarize_normal;
use super::surf_index::{indices, surf_index, SurfType};
use super::{spine_edge_at, ChFi3dBuilder, MoreThreeCorner};

/// `TopExp::FirstVertex(E)`.
fn first_vertex(e: &Edge) -> Option<Vertex> {
    edge_vertices(e).0
}

/// `ChFi3d_cherche_edge(V, E1, F1, E, Vtx)` with the OCCT out-parameter
/// semantics: on failure `E` keeps the value it already held, and only a null
/// `E` after the call raises `Standard_ConstructionError`
/// (`ChFi3d_Builder_0.cxx:5458-5461`).
pub(super) fn cherche_edge_into(
    v: &Vertex,
    excluded: &[TopoShape],
    f1: &Face,
    e: &mut Option<Edge>,
) -> Result<(), String> {
    match cherche_edge(v, excluded, f1) {
        Ok((found, _vtx)) => {
            *e = Some(found);
            Ok(())
        }
        Err(error) => {
            if e.is_some() {
                Ok(())
            } else {
                Err(error)
            }
        }
    }
}

/// `spine->Edges(1)` / `spine->Edges(NbEdges())` of `cxx:1363-1370`.
fn spine_edge(stripe: &StripeHandle, first: bool) -> Result<Edge, String> {
    spine_edge_at(stripe, first)
}

impl ChFi3dBuilder {
    /// `PerformMoreThreeCorner(Jndex, nconges)` up to `cxx:1581`.
    pub(super) fn more_three_corner_topology(
        &mut self,
        index: usize,
        st: &mut MoreThreeCorner,
    ) -> Result<(), String> {
        // `const TopoDS_Vertex& V1 = myVDataMap.FindKey(Jndex);` (`cxx:1246`).
        let vtx_shape = self
            .my_v_data_map
            .find_key(index)
            .cloned()
            .ok_or_else(|| "PerformMoreThreeCorner: myVDataMap index out of range".to_string())?;
        let vertex = Vertex::wrap(vtx_shape.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: the map key is not a vertex".to_string())?;
        st.vtx = vtx_shape;

        // `nedge = ChFi3d_NbNotDegeneratedEdges(V1, myVEMap);` (`cxx:1253`).
        let mut nedge = nb_not_degenerated_edges(&vertex, &self.my_ve_map);

        // `ChFi3d_ChercheBordsLibres(myVEMap, V1, bordlibre, edgelibre1,
        // edgelibre2);` (`cxx:1254`).
        let (bord_libre, edge_libre1, edge_libre2) =
            cherche_bords_libres(&self.my_ve_map, &vertex);

        // `cxx:1255-1263`.
        let mut droit = false;
        if bord_libre {
            nedge = (nedge - 2) / 2 + 2;
            let e1 = edge_libre1
                .as_ref()
                .ok_or_else(|| "PerformMoreThreeCorner: edgelibre1 is null".to_string())?;
            let e2 = edge_libre2
                .as_ref()
                .ok_or_else(|| "PerformMoreThreeCorner: edgelibre2 is null".to_string())?;
            let angedg = angle_edge(&vertex, e1, e2)?.abs();
            droit = (angedg - PI).abs() < 0.01;
        } else {
            nedge /= 2;
        }
        if nedge < 1 {
            // `NCollection_Array1(0, size)` with a negative size raises in OCCT.
            return Err(
                "PerformMoreThreeCorner: NCollection_Array1 invalid range (nedge < 1)".to_string(),
            );
        }
        st.nedge = nedge;
        st.bord_libre = bord_libre;
        st.edge_libre1 = edge_libre1;
        st.edge_libre2 = edge_libre2;
        st.droit = droit;

        // `size = nedge * 2;` and the array allocations (`cxx:1264-1304`).
        let size = (nedge * 2) as usize;
        st.allocate(size);

        // `cxx:1328-1340`: the sewing edge of the vertex.
        let mut couture = false;
        let mut face_couture: Option<Face> = None;
        let mut edge_couture: Option<Edge> = None;
        for value in self.my_vf_map.get_from_key(&st.vtx) {
            let Some(fcur) = Face::wrap(value.clone()) else {
                continue;
            };
            if let Some(e) = couture_on_vertex(&fcur, &vertex) {
                couture = true;
                edge_couture = Some(e);
                face_couture = Some(fcur);
                break;
            }
        }
        st.couture = couture;
        st.face_couture = face_couture.clone();
        st.edge_couture = edge_couture.clone();

        // `RemoveSurfData(myVDataMap, myEFMap, edgecouture, facecouture, V1);`
        // (`cxx:1342-1343`).
        remove_surf_data(
            &self.my_v_data_map,
            &self.my_ef_map,
            edge_couture.as_ref(),
            face_couture.as_ref(),
            &vertex,
        )?;

        // `cxx:1346-1373`: the first stripe of the vertex gives the starting
        // edge and the starting face.
        let stripes: Vec<StripeHandle> = self
            .my_v_data_map
            .find_from_index(index)
            .map(|l| l.to_vec())
            .unwrap_or_default();
        let mut fcur: Option<Face> = None;
        let mut ecur: Option<Edge> = None;
        for cnext in &stripes {
            st.cd[0] = Some(cnext.clone());
            let mut sense = 0;
            let index0 = ch_fi3d_index_of_surf_data(&st.vtx, cnext, &mut sense)?;
            st.index[0] = index0;
            st.sens[0] = sense;
            st.numfa[0][1] = surf_index(&st.cd, 0, index0 as usize, SurfType::Face2)?;
            st.numfa[1][0] = st.numfa[0][1];
            let f = self.ds_face(st.numfa[0][1])?;
            st.fvive[0][1] = Some(f.0.clone());
            st.fvive[1][0] = Some(f.0.clone());
            st.jf[0] = 2;
            let e = spine_edge(cnext, sense == 1)?;
            st.evive[0] = Some(e.0.clone());
            let excluded = st.excluded_edges();
            cherche_edge_into(&vertex, &excluded, &f, &mut ecur)?;
            fcur = Some(f);
            break;
        }
        let (mut ecur, mut fcur) = match (ecur, fcur) {
            (Some(e), Some(f)) => (e, f),
            _ => {
                return Err(
                    "PerformMoreThreeCorner: pb in the parsing of edges and faces".to_string(),
                )
            }
        };

        // `SummarizeNormal(V1, Fcur, Ecur, SumFaceNormalAtV1);` (`cxx:1375`).
        summarize_normal(&mut st.sum_face_normal_at_v1, &st.vtx, &fcur, &ecur)?;

        // `cxx:1377-1540`: the ring walk.
        let nedge_us = nedge as usize;
        let mut nbcouture = 0;
        let mut enext: Option<Edge> = None;
        let mut ii = 1usize;
        while ii < nedge_us {
            if st
                .face_couture
                .as_ref()
                .is_some_and(|fc| is_same(&fcur.0, &fc.0))
                && nbcouture == 0
            {
                enext = st.edge_couture.clone();
                nbcouture += 1;
            } else {
                let excluded = st.excluded_edges();
                cherche_edge_into(&vertex, &excluded, &fcur, &mut enext)?;
            }
            let edge_next = enext.clone().ok_or_else(|| {
                "PerformMoreThreeCorner: pb in the parsing of edges and faces".to_string()
            })?;

            let is_free = st
                .edge_libre1
                .as_ref()
                .is_some_and(|e| is_same(&edge_next.0, &e.0))
                || st
                    .edge_libre2
                    .as_ref()
                    .is_some_and(|e| is_same(&edge_next.0, &e.0));
            if is_free {
                // `cxx:1394-1416`: the free border itself.
                st.cd[ii] = None;
                st.index[ii] = 0;
                st.sens[ii] = -1;
                if first_vertex(&edge_next).is_some_and(|v| is_same(&v.0, &st.vtx)) {
                    st.sens[ii] = 1;
                }
                st.sharp[ii] = true;
                st.evive[ii] = Some(edge_next.0.clone());
                st.jf[ii] = 0;
                let (icplus, _icmoins) = indices(nedge_us, ii);
                st.fvive[ii][icplus] = Some(fcur.0.clone());
                st.fvive[icplus][ii] = Some(fcur.0.clone());
                let num = self.my_ds.change_ds().add_shape(&fcur.0);
                st.numfa[ii][icplus] = num;
                st.numfa[icplus][ii] = num;
                ii += 1;

                // `cxx:1417-1438`: the second free border, on the same face.
                let ecur2 = if st
                    .edge_libre1
                    .as_ref()
                    .is_some_and(|e| is_same(&edge_next.0, &e.0))
                {
                    st.edge_libre2.clone()
                } else {
                    st.edge_libre1.clone()
                };
                let ecur2 = ecur2.ok_or_else(|| {
                    "PerformMoreThreeCorner: the second free border is null".to_string()
                })?;
                let _ = edge_common_faces(self.my_ef_map.get_from_key(&ecur2.0));
                let (icplus, _icmoins) = indices(nedge_us, ii);
                st.fvive[ii][icplus] = Some(fcur.0.clone());
                st.fvive[icplus][ii] = Some(fcur.0.clone());
                let num = self.my_ds.change_ds().add_shape(&fcur.0);
                st.numfa[ii][icplus] = num;
                st.numfa[icplus][ii] = num;
                st.cd[ii] = None;
                st.index[ii] = 0;
                st.sens[ii] = -1;
                if first_vertex(&ecur2).is_some_and(|v| is_same(&v.0, &st.vtx)) {
                    st.sens[ii] = 1;
                }
                st.sharp[ii] = true;
                st.evive[ii] = Some(ecur2.0.clone());
                st.jf[ii] = 0;
            } else {
                // `cxx:1442-1469`: is `Enext` the end edge of one of the
                // stripes of the vertex?
                let mut trouve = false;
                let mut found: Option<StripeHandle> = None;
                let mut found_index = 0i32;
                let mut found_sense = 0i32;
                for it in &stripes {
                    let mut sense = 0;
                    let index_it = ch_fi3d_index_of_surf_data(&st.vtx, it, &mut sense)?;
                    let ee = spine_edge(it, sense == 1)?;
                    if is_same(&edge_next.0, &ee.0) {
                        found = Some(it.clone());
                        found_index = index_it;
                        found_sense = sense;
                        trouve = true;
                        break;
                    }
                }
                // `cxx:1462-1485`.
                if trouve {
                    st.cd[ii] = found;
                    st.index[ii] = found_index;
                    st.sens[ii] = found_sense;
                    st.sharp[ii] = false;
                    st.evive[ii] = Some(edge_next.0.clone());
                } else {
                    st.cd[ii] = None;
                    st.index[ii] = 0;
                    st.sens[ii] = -1;
                    if first_vertex(&edge_next).is_some_and(|v| is_same(&v.0, &st.vtx)) {
                        st.sens[ii] = 1;
                    }
                    st.sharp[ii] = true;
                    st.evive[ii] = Some(edge_next.0.clone());
                    st.jf[ii] = 0;
                }

                // `cxx:1486-1491`: the face of `Enext` other than `Fcur`.
                let mut fnext = fcur.clone();
                if let Ok(f) = cherche_face1(self.my_ef_map.get_from_key(&edge_next.0), &fcur) {
                    fnext = f;
                }
                let (icplus, icmoins) = indices(nedge_us, ii);
                st.fvive[ii][icplus] = Some(fnext.0.clone());
                st.fvive[icplus][ii] = Some(fnext.0.clone());
                let num = self.my_ds.change_ds().add_shape(&fnext.0);
                st.numfa[ii][icplus] = num;
                st.numfa[icplus][ii] = num;

                // `cxx:1494-1535`: `jf` and the face cells of the two sides.
                if trouve {
                    let surf_data = st.index[ii] as usize;
                    let numface2 = surf_index(&st.cd, ii, surf_data, SurfType::Face2)?;
                    if numface2 == st.numfa[ii][icplus] {
                        st.jf[ii] = 2;
                    } else {
                        let numface1 = surf_index(&st.cd, ii, surf_data, SurfType::Face1)?;
                        if numface1 == st.numfa[ii][icplus] {
                            st.jf[ii] = 1;
                        } else {
                            if numface1 == st.numfa[icmoins][ii] {
                                st.jf[ii] = 2;
                                let f2 = self.ds_face(numface2)?;
                                st.fvive[ii][icplus] = Some(f2.0.clone());
                                st.fvive[icplus][ii] = Some(f2.0.clone());
                                let n = self.my_ds.change_ds().add_shape(&f2.0);
                                st.numfa[ii][icplus] = n;
                                st.numfa[icplus][ii] = n;
                            }
                            if numface2 == st.numfa[icmoins][ii] {
                                st.jf[ii] = 1;
                                let f1 = self.ds_face(numface1)?;
                                st.fvive[ii][icplus] = Some(f1.0.clone());
                                st.fvive[icplus][ii] = Some(f1.0.clone());
                                let n = self.my_ds.change_ds().add_shape(&f1.0);
                                st.numfa[ii][icplus] = n;
                                st.numfa[icplus][ii] = n;
                            }
                        }
                    }
                }

                ecur = edge_next;
                fcur = fnext;
                // `SummarizeNormal(V1, Fcur, Ecur, SumFaceNormalAtV1);`
                // (`cxx:1538`).
                summarize_normal(&mut st.sum_face_normal_at_v1, &st.vtx, &fcur, &ecur)?;
            }
            ii += 1;
        }

        // `cxx:1542-1557`: `regul` is the `G1` tangency of the two faces of a
        // sharp edge, unless that edge is the sewing edge.
        for ic in 0..nedge_us {
            if !st.sharp[ic] {
                continue;
            }
            let edge = st
                .evive[ic]
                .clone()
                .ok_or_else(|| "PerformMoreThreeCorner: Evive is null".to_string())?;
            if st
                .edge_couture
                .as_ref()
                .is_some_and(|c| is_same(&edge, &c.0))
            {
                continue;
            }
            let (f1, f2) = edge_common_faces(self.my_ef_map.get_from_key(&edge))
                .ok_or_else(|| "PerformMoreThreeCorner: the edge has no face".to_string())?;
            let e = Edge::wrap(edge)
                .ok_or_else(|| "PerformMoreThreeCorner: Evive is not an edge".to_string())?;
            st.regul[ic] = is_tangent_faces(&e, &f1, &f2, Shape::G1)?;
        }

        // `cxx:1559-1581`: a regular edge tangent to another edge of the vertex
        // is not regular (cts60072).
        for ic in 0..nedge_us {
            if !st.regul[ic] {
                continue;
            }
            let ereg = st
                .evive[ic]
                .clone()
                .ok_or_else(|| "PerformMoreThreeCorner: Evive is null".to_string())?;
            let ereg = Edge::wrap(ereg)
                .ok_or_else(|| "PerformMoreThreeCorner: Evive is not an edge".to_string())?;
            for ind in 0..nedge_us {
                if ind == ic {
                    continue;
                }
                let ecur_ind = st
                    .evive[ind]
                    .clone()
                    .ok_or_else(|| "PerformMoreThreeCorner: Evive is null".to_string())?;
                let ecur_ind = Edge::wrap(ecur_ind)
                    .ok_or_else(|| "PerformMoreThreeCorner: Evive is not an edge".to_string())?;
                let ang = angle_edge(&vertex, &ecur_ind, &ereg)?.abs();
                if ang < 0.01 || (ang - PI).abs() < 0.01 {
                    st.regul[ic] = false;
                    st.tangentregul[ic] = true;
                    break;
                }
            }
        }

        // `Fcur` / `Ecur` stay live for the stages from `cxx:1583` on.
        Ok(())
    }

    /// `DStr.Shape(index)` (`TopOpeBRepDS_DataStructure::Shape`), read as the
    /// face the index refers to.
    fn ds_face(&mut self, index: i32) -> Result<Face, String> {
        let shape = self
            .my_ds
            .change_ds()
            .shape(index, false)
            .ok_or_else(|| format!("PerformMoreThreeCorner: DStr.Shape({index}) is null"))?;
        Face::wrap(shape)
            .ok_or_else(|| format!("PerformMoreThreeCorner: DStr.Shape({index}) is not a face"))
    }
}
