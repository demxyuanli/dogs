//! `GetInvalidEdgesByBounds` (`cxx:8672-8990`): the split filter that
//! intersects the new splits with the bounds of the offset-face splits and
//! decides which splits are inside the solid.
//!
//! Six passes, in the reference's order:
//!
//! 1. `cxx:8694-8700`: a `BOPAlgo_Section` between the splits and the bounds.
//! 2. `cxx:8705-8775`: the `InterfEE` records. A record without a new shape
//!    marks the checked-edge `aE1` invalid when its common part is an edge
//!    (`cxx:8719-8725`); a record with one marks it invalid when `aE2` is an
//!    invalid edge (`cxx:8728-8731`), and - unless both of `aE2`'s vertices are
//!    reachable in the split map (`cxx:8734-8743`) - adds every vertex of the
//!    images of `aE1`'s origins that is not old (`cxx:8746-8773`).
//! 3. `cxx:8784-8834`: an image that the section did not reproduce *whole*
//!    (`!aMSSec.Contains(aEIm)`) with a new bounding vertex gets a fresh
//!    midpoint vertex, mapped back to its source edge (`cxx:8811-8831`).
//! 4. `cxx:8838-8852`: every edge of `theMECheckExt` gets the same midpoint
//!    vertex treatment.
//! 5. `cxx:8854-8896`: a vertex that touches only new or only old edges is
//!    invalid outright.
//! 6. `cxx:8898-8999`: every remaining vertex is checked against the splits of
//!    its faces (a vertex that *is* a split vertex is fine, otherwise its
//!    projection must fall inside one of the face's splits), and - for the
//!    vertices that pass the vertex/edge classification - against the solids of
//!    `mySolids`. Survivors extend `theVertsToAvoid` and mark all their edges
//!    invalid.
//!
//! The `ComputeVF` probe uses [`crate::int_tools_full::IntToolsContext::compute_vf_uv`]
//! because the reference reuses the returned `(U, V)` for the per-split
//! `IsPointInOnFace` test (`cxx:8946-8958`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_builder2::BopBuilder;
use crate::brep_class3d::SolidClassifier;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Face, TopoShape, Vertex};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, explorer_shapes_avoiding, is_same,
    shape_map_key, ShapeMapKey,
};

use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void GetInvalidEdgesByBounds(const Shape& theSplits, const Shape&
    /// theBounds, const Map& theMVOld, const Map& theMENew, const DataMap&
    /// theDMEOr, const DataMap& theMELF, const DataMap& theEImages, Map&
    /// theMECheckExt, const Map& theMEInvOnArt, Map& theVertsToAvoid, Map&
    /// theMEInv)` (`cxx:8672-8990`).
    #[allow(clippy::too_many_arguments)]
    pub fn get_invalid_edges_by_bounds(
        &mut self,
        the_splits: &TopoShape,
        the_bounds: &TopoShape,
        the_mv_old: &ShapeIndexedMap,
        the_menew: &ShapeIndexedMap,
        the_dmeor: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_melf: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_e_images: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_me_check_ext: &ShapeIndexedMap,
        the_me_inv_on_art: &ShapeIndexedMap,
        the_verts_to_avoid: &mut ShapeIndexedMap,
        the_me_inv: &mut ShapeIndexedMap,
    ) {
        // `cxx:8690-8692`: vertex -> its edges, in `FindKey` order.
        let mut a_dmve = Dmvv::new();
        for a_e in explorer_shapes(the_splits, ShapeType::Edge) {
            for a_v in explorer_shapes(&a_e, ShapeType::Vertex) {
                a_dmve.add(&a_v, &a_e);
            }
        }
        // `TopExp::MapShapesAndAncestors` also binds the sub-shapes that have no
        // ancestor; `explorer_shapes_avoiding` is that pass.
        for a_v in explorer_shapes_avoiding(the_splits, ShapeType::Vertex, ShapeType::Edge) {
            if !a_dmve.contains(&a_v) {
                a_dmve.bind_list(&a_v, &[]);
            }
        }

        // `cxx:8695-8700`.
        let mut a_sec = BopBuilder::new();
        a_sec.add_argument(the_splits);
        a_sec.add_argument(the_bounds);
        a_sec.perform_section(&[the_splits.clone()], &[the_bounds.clone()])
            .ok();

        // `cxx:8702-8705`.
        let mut a_mv_inv = ShapeIndexedMap::new();
        let mut a_mv_check_add = ShapeIndexedMap::new();

        // `cxx:8710-8775`: the edge/edge interferences.
        let a_ees: Vec<(usize, usize, i64, bool)> = a_sec
            .ds()
            .interf_ee()
            .iter()
            .map(|a_ee| {
                let b_edge_part = a_ee
                    .common_part
                    .as_ref()
                    .map(|a_cp| a_cp.is_edge_part())
                    .unwrap_or(false);
                (a_ee.index1, a_ee.index2, a_ee.index_new, b_edge_part)
            })
            .collect();
        for (an_i1, an_i2, an_index_new, b_edge_part) in a_ees {
            let Some(a_e1) = a_sec.ds().shape(an_i1).cloned() else {
                continue;
            };
            let Some(a_e2) = a_sec.ds().shape(an_i2).cloned() else {
                continue;
            };
            // `cxx:8719-8726`.
            if an_index_new < 0 {
                if the_me_check_ext.contains(&a_e1) && b_edge_part {
                    the_me_inv.add(&a_e1);
                }
                continue;
            }
            // `cxx:8728-8731`.
            if self.my_invalid_edges.contains(&a_e2) {
                the_me_inv.add(&a_e1);
            }
            // `cxx:8734-8743`: the artificial-invalid edge's own vertices do not
            // need checking when both of them are known split vertices.
            if the_me_inv_on_art.contains(&a_e2) {
                let a_verts = explorer_shapes(&a_e2, ShapeType::Vertex);
                let b_all_known = !a_verts.is_empty()
                    && a_verts.iter().all(|a_v| a_dmve.contains(a_v));
                if b_all_known {
                    continue;
                }
            }
            // `cxx:8746-8773`: every vertex of the images of `aE1`'s origins.
            let Some(a_le_or) = the_dmeor.get(&shape_map_key(&a_e1)).cloned() else {
                continue;
            };
            for a_e_or in &a_le_or {
                let Some(a_le_im) = the_e_images.get(&shape_map_key(a_e_or)) else {
                    continue;
                };
                for a_e_im in a_le_im {
                    for a_v in children_keep_orientation(a_e_im) {
                        if !the_mv_old.contains(&a_v) {
                            a_mv_inv.add(&a_v);
                            a_mv_check_add.add(&a_v);
                        }
                    }
                }
            }
        }

        // `cxx:8784-8790`: the section result and its images.
        let a_sec_r = a_sec.result().clone();
        let mut a_mssec = ShapeIndexedMap::new();
        for a_s in explorer_shapes(&a_sec_r, ShapeType::Vertex) {
            a_mssec.add(&a_s);
        }
        for a_s in explorer_shapes(&a_sec_r, ShapeType::Edge) {
            a_mssec.add(&a_s);
        }
        for a_s in explorer_shapes(&a_sec_r, ShapeType::Wire) {
            a_mssec.add(&a_s);
        }
        for a_s in explorer_shapes(&a_sec_r, ShapeType::Face) {
            a_mssec.add(&a_s);
        }

        // `cxx:8791-8834`: an image the section only reproduced in part.
        for a_e in explorer_shapes(the_splits, ShapeType::Edge) {
            // `cxx:8793-8797`.
            if a_sec.history().is_deleted(&a_e) {
                continue;
            }
            // `cxx:8799-8803`.
            let Some(a_le_im) = a_sec
                .history()
                .images()
                .get(&shape_map_key(&a_e))
                .cloned()
            else {
                continue;
            };
            for a_e_im in &a_le_im {
                if a_mssec.contains(a_e_im) {
                    continue;
                }
                // `cxx:8816-8820`: only an image with a new bounding vertex is
                // suspicious.
                let a_verts = explorer_shapes(a_e_im, ShapeType::Vertex);
                let b_any_new = a_verts.iter().any(|a_v| !the_mv_old.contains(a_v));
                if !b_any_new {
                    continue;
                }
                // `cxx:8821-8830`: a fresh midpoint vertex of the image.
                let Some(a_curve) = BRepTool::edge_curve(&convert_edge(a_e_im)) else {
                    break;
                };
                let (a_first, a_last) = BRepTool::edge_parameters(&convert_edge(a_e_im));
                let a_p = a_curve.d0((a_first + a_last) * 0.5);
                let a_v = TopoBuilder::new().make_vertex(a_p, occt_core::precision::CONFUSION);
                a_dmve.add(&a_v.0, &a_e);
                a_mv_inv.add(&a_v.0);
                break;
            }
        }

        // `cxx:8838-8852`: the edges under extended check get the same probe.
        let a_me_check_ext_list: Vec<TopoShape> = the_me_check_ext.iter().cloned().collect();
        for a_e in &a_me_check_ext_list {
            let Some(a_curve) = BRepTool::edge_curve(&convert_edge(a_e)) else {
                continue;
            };
            let (a_first, a_last) = BRepTool::edge_parameters(&convert_edge(a_e));
            let a_p = a_curve.d0((a_first + a_last) * 0.5);
            let a_v = TopoBuilder::new().make_vertex(a_p, occt_core::precision::CONFUSION);
            a_dmve.add(&a_v.0, a_e);
            a_mv_inv.add(&a_v.0);
        }

        // `cxx:8854-8896`: a vertex that touches only new or only old edges.
        let a_nb = a_dmve.extent();
        for i in 1..=a_nb {
            let Some(a_v) = a_dmve.keys().get(i - 1).cloned() else {
                continue;
            };
            if the_mv_old.contains(&a_v) {
                continue;
            }
            let mut b_new = false;
            let mut b_old = false;
            let a_lex: Vec<TopoShape> = a_dmve.seek(&a_v).map(|l| l.to_vec()).unwrap_or_default();
            for a_e in &a_lex {
                // `cxx:8869-8872`.
                if the_me_check_ext.contains(a_e) {
                    continue;
                }
                if the_menew.contains(a_e) {
                    b_new = true;
                } else {
                    b_old = true;
                }
                if b_new && b_old {
                    break;
                }
            }
            // `cxx:8889-8893`.
            if !b_new || !b_old {
                a_mv_inv.add(&a_v);
                a_mv_check_add.remove(&a_v);
            }
        }

        // `cxx:8898-8989`: the per-vertex classification.
        let a_nb_iv = a_mv_inv.extent();
        for iv in 1..=a_nb_iv {
            let Some(a_v) = a_mv_inv.value(iv).cloned() else {
                continue;
            };
            // `cxx:8900-8904`.
            if the_mv_old.contains(&a_v) {
                continue;
            }
            // `cxx:8906-8910`.
            let Some(a_le_inv) = a_dmve.seek(&a_v).map(|l| l.to_vec()) else {
                continue;
            };

            // `cxx:8912-8922`: the faces behind the vertex's edges.
            let mut a_mf = ShapeIndexedMap::new();
            for a_e in &a_le_inv {
                // `theMELF.Find(aE)` raises on an unbound edge in the
                // reference; the port skips it.
                let Some(a_lf) = the_melf.get(&shape_map_key(a_e)) else {
                    continue;
                };
                for a_f in a_lf {
                    a_mf.add(a_f);
                }
            }

            // `cxx:8924-8963`.
            let a_v_vertex = convert_vertex(&a_v);
            let mut b_invalid = true;
            let a_nb_f = a_mf.extent();
            for i in 1..=a_nb_f {
                if !b_invalid {
                    break;
                }
                let Some(a_f) = a_mf.value(i).cloned() else {
                    continue;
                };
                // `cxx:8933`: `myOFImages.FindFromKey(aF)` raises on an unbound
                // face in the reference; the port skips it.
                let a_lf_im: Vec<TopoShape> = match self.my_of_images.seek(&a_f) {
                    Some(a_l) => a_l.to_vec(),
                    None => continue,
                };
                // `cxx:8936-8943`: a vertex that already bounds one of the
                // splits is valid.
                for a_f_im in &a_lf_im {
                    for a_vf in explorer_shapes(a_f_im, ShapeType::Vertex) {
                        b_invalid = !is_same(&a_vf, &a_v);
                        if !b_invalid {
                            break;
                        }
                    }
                    if !b_invalid {
                        break;
                    }
                }
                // `cxx:8946-8962`: otherwise the vertex must project inside the
                // face and land in one of its splits.
                if b_invalid {
                    let (a_status, a_u, a_v_param, _a_tol) =
                        self.my_context.compute_vf_uv(&a_v_vertex, &convert_face(&a_f), 0.0);
                    if a_status == 0 {
                        for a_f_im in &a_lf_im {
                            b_invalid = !self
                                .my_context
                                .is_point_in_on_face_2d(&convert_face(a_f_im), (a_u, a_v_param));
                            if !b_invalid {
                                break;
                            }
                        }
                    }
                }
            }

            // `cxx:8965-8980`: vertices that passed the vertex/edge
            // classification are also checked against the solids.
            if b_invalid && a_mv_check_add.contains(&a_v) {
                let a_p = BRepTool::vertex_point(&a_v_vertex);
                let a_tol_v = BRepTool::vertex_tolerance(&a_v_vertex);
                if let Some(a_solids) = self.my_solids.clone() {
                    for a_sol in explorer_shapes(&a_solids, ShapeType::Solid) {
                        if !b_invalid {
                            break;
                        }
                        let a_state = SolidClassifier::classify(&a_sol, &a_p, a_tol_v);
                        b_invalid = a_state == FaceState::Out;
                    }
                }
            }

            // `cxx:8982-8989`.
            if b_invalid {
                the_verts_to_avoid.add(&a_v);
                for a_e in &a_le_inv {
                    the_me_inv.add(a_e);
                }
            }
        }
    }
}

/// `TopoDS::Edge(S)` - the reference's cast of a shape known to be an edge.
fn convert_edge(s: &TopoShape) -> crate::shape::Edge {
    crate::shape::Edge(s.clone())
}

/// `TopoDS::Face(S)`.
fn convert_face(s: &TopoShape) -> Face {
    Face(s.clone())
}

/// `TopoDS::Vertex(S)`.
fn convert_vertex(s: &TopoShape) -> Vertex {
    Vertex(s.clone())
}
