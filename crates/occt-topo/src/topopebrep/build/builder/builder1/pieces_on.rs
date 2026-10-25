//! `TopOpeBRepBuild_Builder1` -- the ON-piece comparison of two same-domain
//! edges.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `PerformPieceOn2D(aPieceObj, aFaceObj, anEdgeObj, pieces, faces, out2d)` | `1531-1653` |
//! | `TwoPiecesON(aSeq, aListOfPieces, aListOfFaces, aListOfPiecesOut2d)` | `1656-2110` |
//!
//! Both run inside `GFillEdgeSameDomWES` (`cxx:886-948`), which is why the
//! four out-lists exist: the accepted pieces carry the faces they belong to and
//! the pieces that need a 2d translation of their p-curve.
//!
//! `myProcessedPartsON2d` is written only here (`cxx:1985`, `:1993`, `:2062`),
//! exactly as in OCCT.
//!
//! UNPORTED: the `// mine` / `// ???` commented-out code OCCT keeps beside the
//! live statements (`cxx:1803-1813`, `:1938-1941`, `:1620-1621`, `:1629-1631`).

use occt_core::gp::GpVec;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_same, map_shapes_and_ancestors, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::tools::{
    get_adjacent_face, get_normal_in_nearest_point, get_normal_to_face_on_edge,
    get_tangent_to_edge_edge,
};

use super::super::Builder;

/// `gp_Dir(V)`: the direction of `V`. OCCT raises on a null vector; the port's
/// [`GpVec::normalized`] leaves it untouched instead.
fn gp_dir(v: GpVec) -> GpVec {
    v.normalized()
}

/// `GetNormalToFaceOnEdge(F, E, out)`: OCCT answers `false` without touching
/// `out`, so a failed query keeps the previous normal (the `anyN` chain of
/// `TwoPiecesON` relies on that).
fn normal_on_edge(f: &Face, e: &Edge, out: &mut GpVec) {
    if let Some(n) = get_normal_to_face_on_edge(f, e) {
        *out = n;
    }
}

/// `GetNormalInNearestPoint(F, E, out)`, same keep-on-failure rule.
fn normal_in_nearest_point(f: &Face, e: &Edge, out: &mut GpVec) {
    if let Some(n) = get_normal_in_nearest_point(f, e) {
        *out = n;
    }
}

impl Builder {
    /// `PerformPieceOn2D(aPieceObj, aFaceObj, anEdgeObj, aListOfPieces,
    /// aListOfFaces, aListOfPiecesOut2d)` (`cxx:1531-1653`): runs
    /// [`Builder1::two_pieces_on`] for every same-domain partner of
    /// `anEdgeObj` on the same-domain faces of `aFaceObj`.
    ///
    /// `aMapOfShapeWithStateTool` (`cxx:1550-1553`) is the map of the *other*
    /// operand: the tool map when `iRef == 1`, the object map otherwise. The
    /// port selects exactly that map instead of probing both.
    pub(crate) fn perform_piece_on_2d(
        &mut self,
        a_piece_obj: &TopoShape,
        a_face_obj: &TopoShape,
        an_edge_obj: &TopoShape,
        a_list_of_pieces: &mut Vec<TopoShape>,
        a_list_of_faces: &mut Vec<TopoShape>,
        a_list_of_pieces_out_2d: &mut Vec<TopoShape>,
    ) -> i32 {
        let i_ref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(a_face_obj))
            .unwrap_or(0);

        let (same_domain, edge_same_domain) = {
            let Some(hds) = self.data_structure() else {
                return -1;
            };
            if !hds.has_same_domain(a_face_obj, true) {
                return -1;
            }
            let faces = hds.same_domain(a_face_obj).to_vec();
            let edges = if hds.has_same_domain(an_edge_obj, true) {
                hds.same_domain(an_edge_obj).to_vec()
            } else {
                Vec::new()
            };
            (faces, edges)
        };

        let mut flag = 0;
        let mut cases: Vec<i32> = Vec::new();

        // `cxx:1556-1637`: every same-domain face of `aFaceObj` with its own
        // edge map, cross-checked against the same-domain edges of `anEdgeObj`.
        for a_face_tool in same_domain.iter() {
            let an_edges_tool: Vec<TopoShape> =
                crate::topo_tools::explorer_shapes(a_face_tool, ShapeType::Edge);

            for an_edge_tool in edge_same_domain.iter() {
                if !an_edges_tool
                    .iter()
                    .any(|e| shape_map_key(e) == shape_map_key(an_edge_tool))
                {
                    continue;
                }

                for an_exp_edge_tool in
                    crate::topo_tools::explorer_shapes(a_face_tool, ShapeType::Edge)
                {
                    if !is_same(&an_exp_edge_tool, an_edge_tool) {
                        continue;
                    }

                    let mut an_edge_tool = an_edge_tool.clone();
                    an_edge_tool.set_orientation(an_exp_edge_tool.orientation());

                    // `aMapOfShapeWithStateTool.FindFromKey(anEdgeTool).Part(ON)`
                    let a_part_on_tool: Vec<TopoShape> = {
                        let Some(hds) = self.data_structure() else {
                            break;
                        };
                        let ds = hds.ds();
                        let map = if i_ref == 1 {
                            ds.map_of_shape_with_state_tool()
                        } else {
                            ds.map_of_shape_with_state_obj()
                        };
                        map.find_from_key(&an_edge_tool)
                            .map(|sws| sws.part(State::On).to_vec())
                            .unwrap_or_default()
                    };

                    for mut a_piece_tool in a_part_on_tool.into_iter() {
                        a_piece_tool.set_orientation(an_edge_tool.orientation());

                        let is_deg_flag = BRepTool::is_degenerated(&Edge(a_piece_obj.clone()))
                            && BRepTool::is_degenerated(&Edge(a_piece_tool.clone()));

                        let a_is_same_cnd = if is_deg_flag {
                            crate::topopebrep::build::tools::is_deg_edges_the_same(
                                a_piece_obj,
                                &a_piece_tool,
                            )
                        } else {
                            is_same(a_piece_obj, &a_piece_tool)
                        };

                        if !a_is_same_cnd {
                            continue;
                        }

                        // `cxx:1609-1622`: the six-shape sequence, then the pair.
                        let a_seq: Vec<TopoShape> = vec![
                            a_face_obj.clone(),
                            an_edge_obj.clone(),
                            a_piece_obj.clone(),
                            a_face_tool.clone(),
                            an_edge_tool.clone(),
                            a_piece_tool.clone(),
                        ];

                        flag += 1;
                        let priz = self.two_pieces_on(
                            &a_seq,
                            a_list_of_pieces,
                            a_list_of_faces,
                            a_list_of_pieces_out_2d,
                        );
                        if !cases.contains(&priz) {
                            cases.push(priz);
                        }
                        break;
                    }
                }
            }
        }

        // `cxx:1638-1650`: a piece found twice with both 14 and 12 (or 11 and
        // 13) is not kept.
        if flag > 1 {
            let has = |v: i32| cases.contains(&v);
            if has(14) && has(12) && self.opefus() {
                a_list_of_pieces.clear();
            }
            if has(11) && has(13) && (self.opec12() || self.opec21()) {
                a_list_of_pieces.clear();
            }
        }

        flag
    }

    /// `TwoPiecesON(aSeq, aListOfPieces, aListOfFaces,
    /// aListOfPiecesOut2d)` (`cxx:1656-2110`): decides, for one pair of
    /// same-domain ON pieces, whether each piece belongs to the result.
    ///
    /// The return value is the OCCT "case" number (`1`/`2`/`3`/`4` with usable
    /// adjacent faces, `11`..`14` without, `10` for the "nothing meaningful"
    /// arm, `0` for the unreachable tail, `-2` for a short sequence).
    pub(crate) fn two_pieces_on(
        &mut self,
        a_seq: &[TopoShape],
        a_list_of_pieces: &mut Vec<TopoShape>,
        a_list_of_faces: &mut Vec<TopoShape>,
        a_list_of_pieces_out_2d: &mut Vec<TopoShape>,
    ) -> i32 {
        if a_seq.len() < 6 {
            return -2;
        }
        let mut a_face_obj = a_seq[0].clone();
        let mut an_e_obj = a_seq[1].clone();
        let mut a_piece_obj = a_seq[2].clone();
        let mut a_face_tool = a_seq[3].clone();
        let mut an_e_tool = a_seq[4].clone();
        let mut a_piece_tool = a_seq[5].clone();

        let i_ref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(&a_face_obj))
            .unwrap_or(0);

        // `cxx:1675-1694`: both arms build the maps identically, then the
        // `iRef != 1` arm swaps object and tool.
        let mut an_edge_face_map_obj: std::collections::HashMap<ShapeMapKey, Vec<TopoShape>> =
            std::collections::HashMap::new();
        let mut an_edge_face_map_tool: std::collections::HashMap<ShapeMapKey, Vec<TopoShape>> =
            std::collections::HashMap::new();
        if let Some(s1) = self.my_shape1.clone() {
            an_edge_face_map_obj = map_shapes_and_ancestors(&s1, ShapeType::Edge, ShapeType::Face);
        }
        if let Some(s2) = self.my_shape2.clone() {
            an_edge_face_map_tool = map_shapes_and_ancestors(&s2, ShapeType::Edge, ShapeType::Face);
        }
        if i_ref != 1 {
            std::mem::swap(&mut a_face_obj, &mut a_face_tool);
            std::mem::swap(&mut a_piece_obj, &mut a_piece_tool);
            std::mem::swap(&mut an_e_obj, &mut an_e_tool);
        }

        let is_faces_dif_oriented;
        let is_edges_rev_sense;
        let an_ad1;
        let an_ad2;
        let mut a_sc_pr_flag1 = false;
        let mut a_sc_pr_flag2 = false;
        let mut rejected1 = true;
        let mut rejected2 = true;

        let a_tol = 1.0e-5;

        // ```cxx:1705```: `gp_Vec anyN` starts null and every query either
        // writes it or leaves the previous value.
        let mut any_n = GpVec::new(0.0, 0.0, 0.0);

        let a_f_obj = Face::wrap(a_face_obj.clone());
        let a_f_tool = Face::wrap(a_face_tool.clone());
        let (Some(a_f_obj), Some(a_f_tool)) = (a_f_obj, a_f_tool) else {
            return -2;
        };

        let an_edge_obj = Edge(a_piece_obj.clone());
        let an_edge_tool = Edge(a_piece_tool.clone());
        let a_ori_e_obj = Edge(an_e_obj.clone());
        let a_ori_e_tool = Edge(an_e_tool.clone());

        // Normals to the two faces (`cxx:1720-1733`).
        normal_on_edge(&a_f_obj, &an_edge_obj, &mut any_n);
        if a_f_obj.0.orientation() == Orientation::Reversed {
            any_n.reverse();
        }
        let a_dn_obj = gp_dir(any_n);

        normal_on_edge(&a_f_tool, &an_edge_tool, &mut any_n);
        if a_f_tool.0.orientation() == Orientation::Reversed {
            any_n.reverse();
        }
        let a_dn_tool = gp_dir(any_n);

        is_faces_dif_oriented = a_dn_obj.dot(&a_dn_tool) < 0.0;

        // Sense of the pieces (`cxx:1736-1755`).
        let rev_obj = match get_tangent_to_edge_edge(&an_edge_obj, &a_ori_e_obj) {
            Some((t, rev)) => {
                any_n = t;
                rev
            }
            None => false,
        };
        if rev_obj {
            a_piece_obj.reverse();
            an_e_obj.reverse();
        }
        let a_dt_obj = gp_dir(any_n);

        let rev_tool = match get_tangent_to_edge_edge(&an_edge_tool, &a_ori_e_tool) {
            Some((t, rev)) => {
                any_n = t;
                rev
            }
            None => false,
        };
        if rev_tool {
            a_piece_tool.reverse();
            an_e_tool.reverse();
        }
        let a_dt_tool = gp_dir(any_n);

        is_edges_rev_sense = a_dt_obj.dot(&a_dt_tool) < 0.0;

        // `cxx:1757-1760`: the adjacent faces, when present.
        let an_adj_face_obj = get_adjacent_face(&a_face_obj, &an_e_obj, &an_edge_face_map_obj);
        let an_adj_face_tool = get_adjacent_face(&a_face_tool, &an_e_tool, &an_edge_face_map_tool);
        an_ad1 = an_adj_face_obj.is_some();
        an_ad2 = an_adj_face_tool.is_some();

        // `cxx:1702`: both products start at `0.` and are only computed when
        // the two adjacent faces exist.
        let (a_sc_product_obj, a_sc_product_tool) = if an_ad1 && an_ad2 {
            let an_adj_face_obj = an_adj_face_obj.as_ref().unwrap();
            let an_adj_face_tool = an_adj_face_tool.as_ref().unwrap();

            normal_on_edge(&Face(an_adj_face_obj.clone()), &an_edge_obj, &mut any_n);
            if an_adj_face_obj.orientation() == Orientation::Reversed {
                any_n.reverse();
            }
            let mut a_dn_a_obj = gp_dir(any_n);

            normal_on_edge(
                &Face(an_adj_face_tool.clone()),
                &an_edge_tool,
                &mut any_n,
            );
            if an_adj_face_tool.orientation() == Orientation::Reversed {
                any_n.reverse();
            }
            let mut a_dn_a_tool = gp_dir(any_n);

            let mut a_sc_pr_obj = a_dt_obj.crossed(&a_dn_tool).dot(&a_dn_a_obj);
            let mut a_sc_pr_tool = a_dt_tool.crossed(&a_dn_obj).dot(&a_dn_a_tool);

            if a_sc_pr_obj.abs() <= a_tol {
                normal_in_nearest_point(
                    &Face(an_adj_face_obj.clone()),
                    &an_edge_obj,
                    &mut any_n,
                );
                if an_adj_face_obj.orientation() == Orientation::Reversed {
                    any_n.reverse();
                }
                a_dn_a_obj = gp_dir(any_n);
                a_sc_pr_obj = a_dt_obj.crossed(&a_dn_tool).dot(&a_dn_a_obj);
            }

            if a_sc_pr_tool.abs() <= a_tol {
                normal_in_nearest_point(
                    &Face(an_adj_face_tool.clone()),
                    &an_edge_tool,
                    &mut any_n,
                );
                if an_adj_face_tool.orientation() == Orientation::Reversed {
                    any_n.reverse();
                }
                a_dn_a_tool = gp_dir(any_n);
                a_sc_pr_tool = a_dt_tool.crossed(&a_dn_obj).dot(&a_dn_a_tool);
            }

            let sc_product_obj = a_sc_pr_obj;
            let sc_product_tool = a_sc_pr_tool;

            if sc_product_tool.abs() > a_tol {
                a_sc_pr_flag1 = true;
            }
            if sc_product_obj.abs() > a_tol {
                a_sc_pr_flag2 = true;
            }
            (sc_product_obj, sc_product_tool)
        } else {
            (0.0, 0.0)
        };

        // Management (`cxx:1824-2109`).
        if !an_ad1 || !an_ad2 || !a_sc_pr_flag1 || !a_sc_pr_flag2 {
            // case a.  No==Nt , To!=Tt
            if is_edges_rev_sense && !is_faces_dif_oriented {
                if self.opec12() {
                    let poisc = BRepTool::is_closed_edge_face(&an_edge_obj, &a_f_obj);
                    if !poisc {
                        a_list_of_pieces.push(a_piece_obj.clone());
                        a_list_of_faces.push(a_face_obj.clone());
                    }
                }
                return 11;
            }

            // case b.  No!=Nt , To!=Tt
            if is_edges_rev_sense && is_faces_dif_oriented {
                if self.opec12() {
                    a_list_of_pieces.push(a_piece_obj.clone());
                    a_list_of_faces.push(a_face_obj.clone());
                }
                if !an_ad1 || !an_ad2 {
                    return 12;
                }
                return 10;
            }

            // case c.  No==Nt , To==Tt
            if !is_edges_rev_sense && !is_faces_dif_oriented {
                if an_ad1 && an_ad2 {
                    if !self.opecom() {
                        if !a_sc_pr_flag2 {
                            a_list_of_pieces.push(a_piece_obj.clone());
                            a_list_of_faces.push(a_face_obj.clone());
                        }
                        if !a_sc_pr_flag1 {
                            a_list_of_pieces.push(a_piece_tool.clone());
                            a_list_of_faces.push(a_face_tool.clone());
                        }
                    }
                } else if self.opefus() {
                    a_list_of_pieces.push(a_piece_obj.clone());
                    a_list_of_faces.push(a_face_obj.clone());
                }
                return 13;
            }

            // case d.  No!=Nt , To==Tt
            if !is_edges_rev_sense && is_faces_dif_oriented {
                if an_ad1 && an_ad2 {
                    if !self.opecom() {
                        if !a_sc_pr_flag2 {
                            a_list_of_pieces.push(a_piece_obj.clone());
                            a_list_of_faces.push(a_face_obj.clone());
                        }
                        if !a_sc_pr_flag1 {
                            a_list_of_pieces.push(a_piece_tool.clone());
                            a_list_of_faces.push(a_face_tool.clone());
                        }
                    }
                } else if self.opefus() {
                    a_list_of_pieces.push(a_piece_obj.clone());
                    a_list_of_faces.push(a_face_obj.clone());
                }
                if !an_ad1 || !an_ad2 {
                    return 14;
                }
                return 10;
            }
            return 10;
        }

        // We can use the adjacents (`cxx:1944-2109`).
        let a_state_obj = if a_sc_product_obj > 0.0 {
            State::In
        } else {
            State::Out
        };
        let a_state_tool = if a_sc_product_tool > 0.0 {
            State::In
        } else {
            State::Out
        };

        // case I  RevSense && DifOriented
        if is_edges_rev_sense && is_faces_dif_oriented {
            if self.opec12() {
                a_list_of_pieces.push(a_piece_obj.clone());
                a_list_of_faces.push(a_face_obj.clone());
            }
            return 1;
        }

        // case III SameSense && !DifOriented
        if !is_edges_rev_sense && !is_faces_dif_oriented {
            if !self.opec12() {
                a_list_of_pieces.push(a_piece_obj.clone());
                a_list_of_faces.push(a_face_obj.clone());
            }
            return 3;
        }

        let st_out = State::Out;
        let st_in = State::In;

        // case II  RevSense && !DifOriented
        if is_edges_rev_sense && !is_faces_dif_oriented {
            if self.opefus() {
                if a_state_obj == st_out && a_state_tool == st_in {
                    rejected1 = false;
                } else if a_state_obj == st_in && a_state_tool == st_out {
                    rejected2 = false;
                } else if a_state_obj == st_in && a_state_tool == st_in {
                    if !self.b1.my_processed_parts_on_2d.contains(&a_piece_obj) {
                        self.b1.my_processed_parts_on_2d.add(&a_piece_obj);
                        self.is_same_2d(a_seq, a_list_of_pieces_out_2d);
                    }
                } else if a_state_obj == st_out && a_state_tool == st_out {
                    if !self.b1.my_processed_parts_on_2d.contains(&a_piece_obj) {
                        self.b1.my_processed_parts_on_2d.add(&a_piece_obj);
                        self.is_same_2d(a_seq, a_list_of_pieces_out_2d);
                    }
                }
            }

            if self.opecom() {
                if a_state_obj == st_out && a_state_tool == st_in {
                    rejected2 = false;
                } else if a_state_obj == st_in && a_state_tool == st_out {
                    rejected1 = false;
                }
            }

            if self.opec12() {
                if a_state_obj == st_out && a_state_tool == st_out {
                    rejected2 = false;
                } else if a_state_obj == st_in && a_state_tool == st_in {
                    rejected1 = false;
                } else if a_state_obj == st_out && a_state_tool == st_in {
                    rejected1 = false;
                    rejected2 = false;
                }
            }
            if !rejected1 {
                a_list_of_pieces.push(a_piece_obj.clone());
                a_list_of_faces.push(a_face_obj.clone());
            }
            if !rejected2 {
                a_list_of_pieces.push(a_piece_tool.clone());
                a_list_of_faces.push(a_face_tool.clone());
            }
            return 2;
        }

        // case IV   !RevSense && DifOriented
        if !is_edges_rev_sense && is_faces_dif_oriented {
            if self.opefus() {
                if a_state_obj == st_out && a_state_tool == st_out {
                    rejected1 = false;
                    rejected2 = false;
                } else if a_state_obj == st_out && a_state_tool == st_in {
                    rejected2 = false;
                } else if a_state_obj == st_in && a_state_tool == st_out {
                    rejected1 = false;
                } else if a_state_obj == st_in && a_state_tool == st_in {
                    if !self.b1.my_processed_parts_on_2d.contains(&a_piece_obj) {
                        self.b1.my_processed_parts_on_2d.add(&a_piece_obj);
                        self.is_same_2d(a_seq, a_list_of_pieces_out_2d);
                    }
                }
            }

            if self.opecom() {
                if a_state_obj == st_in && a_state_tool == st_in {
                    rejected1 = false;
                    rejected2 = false;
                } else if a_state_obj == st_out && a_state_tool == st_in {
                    rejected1 = false;
                } else if a_state_obj == st_in && a_state_tool == st_out {
                    rejected2 = false;
                }
            }

            if self.opec12() {
                if a_state_obj == st_out && a_state_tool == st_out {
                    rejected1 = false;
                } else if a_state_obj == st_in && a_state_tool == st_in {
                    rejected2 = false;
                }
            }
            if !rejected1 {
                a_list_of_pieces.push(a_piece_obj.clone());
                a_list_of_faces.push(a_face_obj.clone());
            }
            if !rejected2 {
                a_list_of_pieces.push(a_piece_tool.clone());
                a_list_of_faces.push(a_face_tool.clone());
            }
            return 4;
        }

        // Unknown case for existing adjacents (`cxx:2107-2108`).
        0
    }
}
