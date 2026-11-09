//! `ProcessCommonEdges` (`cxx:6975-7106`), the file-local `FindOrigins`
//! (`cxx:7113-7152`) and `UpdateIntersectedFaces` (`cxx:7159-7227`).
//!
//! `ProcessCommonEdges` keeps the common edges that really belong to the
//! invalidity: an edge that is invalid but not valid is dropped
//! (`cxx:6998-7001`), an edge touching the invalidity either directly
//! (`cxx:7005-7009`) or through a vertex (`cxx:7010-7019`) is kept. When no
//! edge survives, either the whole list is reused (`theForceUse`,
//! `cxx:7027-7029`), or a single edge of the list is kept - provided none of the
//! trimmed parts sharing its infinite edge reaches an invalidity
//! (`cxx:7032-7063`). The kept edges contribute their vertices to
//! `theMECV` (`cxx:7075-7087`), bind infinite edge -> trimmed part in
//! `theDMEETrim` (`cxx:7089-7097`), and give `theLFEi` / `theLFEj` / `theMEToInt`
//! their entries (`cxx:7102-7105`).
//!
//! `FindOrigins` collects the origins of the edges shared by two split lists
//! (`cxx:7122-7150`), and `UpdateIntersectedFaces` adds the common edges to
//! `theMEToInt` and merges the found origins into `myEdgesOrigins`
//! (`cxx:7202-7224`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, is_same, shape_map_key, ShapeMapKey,
};

use super::statics::{append_to_list, find_common_parts};
use super::BRepOffsetBuildOffsetFaces;

/// `static void FindOrigins(const List& theLFIm1, const List& theLFIm2, const
/// IndexedMap& theME, const DataMap& theOrigins, List& theLEOr)`
/// (`cxx:7113-7152`).
pub(super) fn find_origins(
    the_lf_im1: &[TopoShape],
    the_lf_im2: &[TopoShape],
    the_me: &ShapeIndexedMap,
    the_origins: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_le_or: &mut Vec<TopoShape>,
) {
    // `cxx:7121`.
    let mut a_mfence = ShapeIndexedMap::new();
    for i in 0..2 {
        // `cxx:7125-7126`.
        let a_lf: &[TopoShape] = if i == 0 { the_lf_im1 } else { the_lf_im2 };
        for a_f in a_lf {
            for a_e in explorer_shapes(a_f, ShapeType::Edge) {
                // `cxx:7136-7149`.
                if the_me.contains(&a_e) {
                    if let Some(a_le_or_src) = the_origins.get(&shape_map_key(&a_e)) {
                        for a_e_or in a_le_or_src.clone() {
                            // `cxx:7142-7145`.
                            let a_added = !a_mfence.contains(&a_e_or);
                            a_mfence.add(&a_e_or);
                            if a_added && a_e_or.shape_type() == ShapeType::Edge {
                                append_to_list(the_le_or, &a_e_or);
                            }
                        }
                    }
                }
            }
        }
    }
}

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::ProcessCommonEdges(const List&
    /// theLEC, const IndexedMap& theME, const DataMap& theMEInfETrim, const
    /// Map& theAllInvs, const bool theForceUse, IndexedMap& theMECV, Map&
    /// theMECheckExt, DataMap& theDMEETrim, List& theLFEi, List& theLFEj,
    /// IndexedMap& theMEToInt)` (`cxx:6975-7106`).
    #[allow(clippy::too_many_arguments)]
    pub fn process_common_edges(
        &self,
        the_lec: &[TopoShape],
        the_me: &ShapeIndexedMap,
        the_me_inf_e_trim: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_all_invs: &ShapeIndexedMap,
        the_force_use: bool,
        the_mecv: &mut ShapeIndexedMap,
        the_me_check_ext: &mut ShapeIndexedMap,
        the_dmeet_trim: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_lfei: &mut Vec<TopoShape>,
        the_lfej: &mut Vec<TopoShape>,
        the_me_to_int: &mut ShapeIndexedMap,
    ) {
        // `cxx:6990-7020`: the common edges connected to the invalidity.
        let mut a_lec: Vec<TopoShape> = Vec::new();
        for a_ec in the_lec {
            // `cxx:6998-7001`.
            if self.my_invalid_edges.contains(a_ec) && !self.my_valid_edges.contains(a_ec) {
                continue;
            }
            // `cxx:7004-7009`.
            if the_me.contains(a_ec) {
                append_to_list(&mut a_lec, a_ec);
                continue;
            }
            // `cxx:7011-7019`: `TopoDS_Iterator` over the edge.
            for a_ve in children_keep_orientation(a_ec) {
                if the_me.contains(&a_ve) {
                    append_to_list(&mut a_lec, a_ec);
                    break;
                }
            }
        }

        // `cxx:7022-7065`.
        let b_use_only_inf = a_lec.is_empty();
        if b_use_only_inf {
            if the_force_use {
                // `cxx:7027-7029`.
                a_lec = the_lec.to_vec();
            } else {
                for a_ec in the_lec {
                    // `myETrimEInf->Find` raises on an unbound edge in the
                    // reference; the port skips it.
                    let Some(a_e_int) = self.my_e_trim_e_inf.get(&shape_map_key(a_ec)).cloned()
                    else {
                        continue;
                    };
                    // `theMEInfETrim.Find` raises on an unbound edge in the
                    // reference; the port reads it as an empty list.
                    let a_lve = the_me_inf_e_trim
                        .get(&shape_map_key(&a_e_int))
                        .cloned()
                        .unwrap_or_default();
                    for a_ecx in &a_lve {
                        // `cxx:7044-7048`.
                        if the_all_invs.contains(a_ecx) || self.my_invalid_edges.contains(a_ecx) {
                            return;
                        }
                        // `cxx:7050-7055`.
                        for a_v in children_keep_orientation(a_ecx) {
                            if the_all_invs.contains(&a_v) {
                                return;
                            }
                        }
                    }
                    // `cxx:7058-7061`: use only one element.
                    if a_lec.is_empty() {
                        append_to_list(&mut a_lec, a_ec);
                    }
                }
            }
        }

        // `cxx:7067-7105`.
        for a_ec in &a_lec {
            // `myETrimEInf->Find` raises on an unbound edge in the reference;
            // the port skips it.
            let Some(a_e_int) = self.my_e_trim_e_inf.get(&shape_map_key(a_ec)).cloned() else {
                continue;
            };
            if !b_use_only_inf {
                // `cxx:7073-7087`: the vertices of the parts of the same
                // original edge.
                let a_lve = the_me_inf_e_trim
                    .get(&shape_map_key(&a_e_int))
                    .cloned()
                    .unwrap_or_default();
                for a_ecx in &a_lve {
                    let p_le_or = self.my_oe_origins.get(&shape_map_key(a_ecx));
                    // `cxx:7084-7087`.
                    match p_le_or {
                        Some(l) if l.len() != 1 => {}
                        _ => {
                            for a_v in explorer_shapes(a_ecx, ShapeType::Vertex) {
                                the_mecv.add(&a_v);
                            }
                        }
                    }
                }
                // `cxx:7089-7097`: bind the unlimited edge to its trimmed
                // part.
                append_to_list(
                    the_dmeet_trim.entry(shape_map_key(&a_e_int)).or_default(),
                    a_ec,
                );
            } else if !the_force_use {
                // `cxx:7099-7102`.
                the_me_check_ext.add(&a_e_int);
            }
            // `cxx:7104-7106`.
            append_to_list(the_lfei, &a_e_int);
            append_to_list(the_lfej, &a_e_int);
            the_me_to_int.add(&a_e_int);
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::UpdateIntersectedFaces(const
    /// TopoDS_Shape& theFInv, const TopoDS_Shape& theFi, const TopoDS_Shape&
    /// theFj, const List& theLFInv, const List& theLFImi, const List&
    /// theLFImj, const List& theLFEi, const List& theLFEj, IndexedMap&
    /// theMEToInt)` (`cxx:7159-7227`).
    #[allow(clippy::too_many_arguments)]
    pub fn update_intersected_faces(
        &mut self,
        the_f_inv: &TopoShape,
        the_fi: &TopoShape,
        the_fj: &TopoShape,
        the_lf_inv: &[TopoShape],
        the_lf_imi: &[TopoShape],
        the_lf_imj: &[TopoShape],
        the_lfei: &[TopoShape],
        the_lfej: &[TopoShape],
        the_me_to_int: &mut ShapeIndexedMap,
    ) {
        // `cxx:7171-7178`: the common edges of the two lists.
        let mut a_mei = ShapeIndexedMap::new();
        for a_e in the_lfei {
            a_mei.add(a_e);
        }

        // `cxx:7180-7188`: the parts of the invalid splits shared with each
        // face.
        let mut a_le_to_find_origins: Vec<TopoShape> = Vec::new();
        if !is_same(the_fi, the_f_inv) {
            find_common_parts(
                the_lf_imi,
                the_lf_inv,
                &mut a_le_to_find_origins,
                ShapeType::Edge,
            );
        }
        if !is_same(the_fj, the_f_inv) {
            find_common_parts(
                the_lf_imj,
                the_lf_inv,
                &mut a_le_to_find_origins,
                ShapeType::Edge,
            );
        }

        // `cxx:7191-7198`.
        let mut a_me_to_find_origins = ShapeIndexedMap::new();
        for a_ec in &a_le_to_find_origins {
            a_me_to_find_origins.add(a_ec);
        }
        let mut a_le_or_init: Vec<TopoShape> = Vec::new();
        find_origins(
            the_lf_imi,
            the_lf_imj,
            &a_me_to_find_origins,
            &self.my_edges_origins,
            &mut a_le_or_init,
        );

        // `cxx:7201-7225`.
        for a_e in the_lfej {
            if !a_mei.contains(a_e) {
                continue;
            }
            the_me_to_int.add(a_e);
            if a_le_or_init.is_empty() {
                continue;
            }
            // `cxx:7210-7222`: merge the origins of the shared edges.
            let a_le_or = self
                .my_edges_origins
                .entry(shape_map_key(a_e))
                .or_default();
            for a_e_or in &a_le_or_init {
                append_to_list(a_le_or, a_e_or);
            }
        }
    }
}
