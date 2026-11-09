//! `LocOpe_Gluer::Perform` (`LocOpe_Gluer.cxx:154-333`).
//!
//! Stage order, as the reference writes it:
//!
//! 1. bind every entry of `myMapEF`: an edge key becomes a
//!    `LocOpe_WiresOnShape::Bind(edge, face)`, a face key becomes a
//!    `LocOpe_GluedShape::GlueOnFace` (`cxx:171-188`);
//! 2. replay `myMapEE` as edge-to-edge binds (`cxx:190-194`);
//! 3. `BindAll` the projection, and split the basis with `LocOpe_Spliter`
//!    (`cxx:196-207`);
//! 4. seed `myDescF` with the split's descendants for every basis face, and
//!    with the face itself for every oriented face of the glued shape
//!    (`cxx:209-225`);
//! 5. rebuild with `LocOpe_Generator` and map every recorded descendant through
//!    `DescendantFace` (`cxx:227-258`);
//! 6. re-encode regularity on the edges that moved to another face pair
//!    (`cxx:260-298`);
//! 7. record the section edges that meet two target faces, forcing `G1` on
//!    them (`cxx:300-332`).

use crate::abs::ShapeType;
use crate::brep_lib_continuity::{set_stored_continuity, stored_continuity};
use crate::locope::generated_shape::GeneratedShape;
use crate::locope::glued_shape::GluedShape;
use crate::locope::generator::Generator;
use crate::locope::spliter::Spliter;
use crate::locope::tgt_faces::tgt_faces;
use crate::locope::wires_on_shape::WiresOnShape;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, AncestorMap};

use super::statics::{add_edges, contains};
use super::{Gluer, CONSTRUCTION_ERROR};

/// `GeomAbs_G1` (`GeomAbs_Shape.hxx:26`): the code
/// `BRep_Builder::Continuity` / `BRep_Tool::Continuity` store in the
/// `BRep_TEdge` continuity slot (`brep_lib_continuity`).
const GEOM_ABS_G1: u8 = 1;

impl Gluer {
    /// `Perform()` (`cxx:154-333`).
    ///
    /// `Err` is the `Standard_ConstructionError` of `cxx:165` (no basis or glued
    /// shape, nothing bound, or an operation the binds invalidated) and the
    /// `LocOpe_Spliter` / `LocOpe_WiresOnShape` failures it propagates.
    pub fn perform(&mut self) -> Result<(), String> {
        // `cxx:158-162`.
        if self.my_done {
            return Ok(());
        }
        // `cxx:163-166`.
        if self.my_sb.is_none()
            || self.my_sn.is_none()
            || self.my_map_ef.is_empty()
            || self.my_ope.is_invalid()
        {
            return Err(CONSTRUCTION_ERROR.to_string());
        }

        let Some(my_sb) = self.my_sb.clone() else {
            return Err(CONSTRUCTION_ERROR.to_string());
        };
        let Some(my_sn) = self.my_sn.clone() else {
            return Err(CONSTRUCTION_ERROR.to_string());
        };

        // `cxx:168-169`.
        let mut the_wons = WiresOnShape::new(my_sb.clone());
        let mut the_gs = GluedShape::with_shape(&my_sn);

        // `cxx:171-188`.
        let lmap = self.my_map_ef.extent();
        for ind in 1..=lmap {
            let s = self.my_map_ef.find_key(ind).clone();
            if s.shape_type() == ShapeType::Edge {
                let s2 = self.my_map_ef.find_from_index(ind).clone();
                // `cxx:178-184`: a nullified entry is an edge that borders two
                // glued faces bound to different basis faces; skip it.
                let Some(s2) = s2 else { continue };
                let Some(edge) = Edge::wrap(s) else {
                    return Err(CONSTRUCTION_ERROR.to_string());
                };
                let Some(face) = Face::wrap(s2) else {
                    return Err(CONSTRUCTION_ERROR.to_string());
                };
                the_wons.bind_edge(&edge, &face)?;
            } else {
                // `cxx:186-188`: TopAbs_FACE.
                let Some(face) = Face::wrap(s) else {
                    return Err(CONSTRUCTION_ERROR.to_string());
                };
                the_gs.glue_on_face(&face)?;
            }
        }

        // `cxx:190-194`.
        for (key, value) in self.my_map_ee.iter() {
            let (Some(e_wir), Some(e_fac)) = (Edge::wrap(key.clone()), Edge::wrap(value.clone()))
            else {
                return Err(CONSTRUCTION_ERROR.to_string());
            };
            the_wons.bind_edge_to_edge(&e_wir, &e_fac);
        }

        // `cxx:196-201`.
        the_wons.bind_all()?;
        if !the_wons.is_done() {
            return Ok(());
        }

        // `cxx:203-207`.
        let mut the_split = Spliter::with_shape(&my_sb);
        the_split.perform(&mut the_wons)?;
        if !the_split.is_done() {
            return Ok(());
        }

        // `cxx:209-217`: the split's descendants for every basis face.
        for fac in explorer_shapes(&my_sb, ShapeType::Face) {
            let desc = the_split.descendant_shapes(&fac);
            self.my_desc_f.bind(&fac, desc);
        }
        // `cxx:219-225`: the glued faces that come back correctly oriented are
        // their own descendant.
        for fac in explorer_shapes(&my_sn, ShapeType::Face) {
            self.my_desc_f.bind(&fac, Vec::new());
            if contains(the_gs.oriented_faces(), &fac) {
                if let Some(list) = self.my_desc_f.change(&fac) {
                    list.push(fac.clone());
                }
            }
        }

        // `cxx:227-231`.
        let Some(res_split) = the_split.resulting_shape().cloned() else {
            return Ok(());
        };
        let mut the_gen = Generator::with_shape(&res_split);
        the_gen.perform(&the_gs)?;
        self.my_done = the_gen.is_done();
        if self.my_done {
            // `cxx:233-235`.
            self.my_res = the_gen.resulting_shape().ok().cloned();
            add_edges();

            // `cxx:237-258`: every recorded descendant is pushed through the
            // generator's own descendant table.
            let keys: Vec<TopoShape> = self.my_desc_f.iter().map(|(k, _)| k.clone()).collect();
            for key in keys {
                let children = self.my_desc_f.find(&key).cloned().unwrap_or_default();
                let mut new_desc = Vec::new();
                for child in children {
                    let Some(child_face) = Face::wrap(child) else {
                        continue;
                    };
                    // `cxx:249-253`: an unbound face raises in the reference;
                    // `Generator::descendant_face` answers `Ok(None)` for it and
                    // a null descendant never reaches this port.
                    if let Ok(Some(list)) = the_gen.descendant_face(&child_face) {
                        new_desc.extend(list.iter().cloned());
                    }
                }
                if let Some(slot) = self.my_desc_f.change(&key) {
                    *slot = new_desc;
                }
            }
        }

        // `cxx:260-332` run unconditionally in the reference, but with `myRes`
        // still null when the generator failed: `TopExp::MapShapesAndAncestors`
        // then yields an empty map, no edge is matched, and neither list grows.
        let Some(my_res) = self.my_res.clone() else {
            return Ok(());
        };

        // `cxx:260-298`: re-encode regularity across the glued boundaries.
        let the_map_ef1 = AncestorMap::map_shapes_and_ancestors(&my_sn, ShapeType::Edge, ShapeType::Face);
        let the_map_ef2 = AncestorMap::map_shapes_and_ancestors(&my_res, ShapeType::Edge, ShapeType::Face);
        for ind in 0..the_map_ef1.length() {
            let edg = the_map_ef1.find_key(ind).clone();
            let ll = the_map_ef1.ancestors(ind);
            if ll.len() != 2 {
                continue;
            }
            let fac1 = ll[0].clone();
            let fac2 = ll[ll.len() - 1].clone();
            let thecont = stored_continuity(
                &Edge(edg.clone()),
                &Face(fac1.clone()),
                &Face(fac2.clone()),
            );
            if thecont < GEOM_ABS_G1 {
                continue;
            }
            // `cxx:279-295`: only recode when the pair actually changed.
            let Some(ind2) = the_map_ef2.find_index(&edg) else {
                continue;
            };
            let ll2 = the_map_ef2.ancestors(ind2);
            if ll2.len() != 2 {
                continue;
            }
            let ff1 = ll2[0].clone();
            let ff2 = ll2[ll2.len() - 1].clone();
            let same_pair = (is_same(&ff1, &fac1)
                && is_same(&ff2, &fac2))
                || (is_same(&ff1, &fac2)
                    && is_same(&ff2, &fac1));
            if same_pair {
                continue;
            }
            set_stored_continuity(&Edge(edg), &Face(ff1), &Face(ff2), thecont);
        }

        // `cxx:300-332`: the section edges that meet two target faces.
        the_wons.init_edge_iterator()?;
        while the_wons.more_edge() {
            let edg = the_wons.edge();
            for ind in 0..the_map_ef2.length() {
                if !is_same(the_map_ef2.find_key(ind), &edg.0) {
                    continue;
                }
                self.my_edges.push(edg.0.clone());
                let l = the_map_ef2.ancestors(ind);
                if l.len() != 2 {
                    continue;
                }
                let fac1 = Face(l[0].clone());
                let fac2 = Face(l[l.len() - 1].clone());
                if tgt_faces(&edg, &fac1, &fac2) {
                    self.my_tgt_edges.push(edg.0.clone());
                    let thecont = stored_continuity(&edg, &fac1, &fac2);
                    if thecont < GEOM_ABS_G1 {
                        set_stored_continuity(&edg, &fac1, &fac2, GEOM_ABS_G1);
                    }
                }
            }
            the_wons.next_edge();
        }
        Ok(())
    }
}
