//! `TopOpeBRepTool_REGUW::SplitEds` and `MapS`.
//!
//! Source: `TopOpeBRepTool_REGUW.cxx:187-276` (`SplitEds`) and `:322-480`
//! (`MapS`).

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::super::closed::{
    closed_e, is_closing_e_in_wire, ori_in_sor, ori_in_sor_closed, CLOSING, INTERNAL,
};
use super::super::split_e::split_e;
use super::super::connexity::Connexity;
use super::Reguw;

impl Reguw {
    /// `SplitEds()` (`cxx:187-276`): split every edge that carries an INTERNAL
    /// vertex, then hand each split edge to its vertices' connexity.
    pub fn split_eds(&mut self) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");

        // `cxx:199-216`: `mehasIv` = the edges that `IsInternal` reports for at
        // least one vertex.
        let mut me_has_iv: Vec<TopoShape> = Vec::new();
        let mut seen: HashSet<ShapeMapKey> = HashSet::new();
        let v_keys = self.my_v_order.clone();
        for key in v_keys {
            let Some(co) = self.my_map_v_eds.get(&key) else {
                continue;
            };
            let mut loe: Vec<TopoShape> = Vec::new();
            let ni = co.is_internal(&mut loe);
            if ni == 0 {
                continue;
            }
            for e in loe {
                if seen.insert(shape_map_key(&e)) {
                    me_has_iv.push(e);
                }
            }
        }

        // `cxx:218-274`
        for e in me_has_iv {
            let e_key = shape_map_key(&e);
            let mut splits: Vec<TopoShape> = Vec::new();
            let mut issp = false;
            let isdone = self.my_e_splits.contains_key(&e_key);
            if isdone {
                // `cxx:225-228`
                splits = self.my_e_splits.get(&e_key).cloned().unwrap_or_default();
            } else {
                // `cxx:229-232`
                issp = split_e(&e, &mut splits);
            }
            if issp {
                self.hasnew_splits = true;
            }
            if !(issp || isdone) {
                // `cxx:239-242`: `nnnnnnnnnnnnnnnnnnnnnnnnnnnnyi`
                continue;
            }

            // `cxx:245-267`: each split edge `esp` takes the place of `e` in
            // the connexity of its vertices.
            for esp in splits {
                for vv in explorer_shapes(&esp, ShapeType::Vertex) {
                    let vv_key = shape_map_key(&vv);
                    let Some(co) = self.my_map_v_eds.get_mut(&vv_key) else {
                        continue;
                    };
                    if !co.remove_item_ori(INTERNAL, &e) {
                        continue;
                    }
                    let ivv = ori_in_sor_closed(&vv, &esp);
                    co.add_item(ivv, &esp);
                }
            }
        }
        true
    }

    /// `MapS()` (`cxx:322-480`): bind every vertex of the wire to the edges it
    /// belongs to, by orientation. A vertex of a singular edge (degenerated,
    /// closing, or the closure vertex of a closed edge) is bound through the
    /// `CLOSING` slot instead of its own orientation.
    pub fn map_s(&mut self) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");

        // `cxx:362`: `const TopoDS_Shape& CS = myCORRISO.S();`
        let Some(cs) = self.my_corr_iso.s().cloned() else {
            return false;
        };
        let f_ref = self.my_corr_iso.f_ref().clone();

        // `cxx:363-412`
        for ed in explorer_shapes(&cs, ShapeType::Edge) {
            let ed_key = shape_map_key(&ed);
            let _ = ed_key;

            let isdg_e = BRepTool::is_degenerated(&Edge(ed.clone()));
            // `cxx:368-369`
            let mut isc_e = is_closing_e_in_wire(&ed, &cs, &f_ref);
            isc_e = isc_e && !isdg_e;

            // `cxx:370-372`
            let vcl = closed_e(&ed);
            let isvcl = vcl.is_some() && !isdg_e;

            for v in explorer_shapes(&ed, ShapeType::Vertex) {
                let key = shape_map_key(&v);
                if !self.my_map_v_eds.contains_key(&key) {
                    self.my_map_v_eds.insert(key.clone(), Connexity::with_key(&v));
                    self.my_v_order.push(key.clone());
                }
                let Some(co) = self.my_map_v_eds.get_mut(&key) else {
                    continue;
                };

                if isdg_e {
                    // `cxx:387-392`: avoid adding it twice.
                    co.remove_item_ori(CLOSING, &ed);
                    co.add_item(CLOSING, &ed);
                    break;
                } else if isc_e {
                    // `cxx:393-397`: the FORWARD and the REVERSED edge are added.
                    co.add_item(CLOSING, &ed);
                } else if isvcl && matches!(&vcl, Some(vcl) if is_same(&v, vcl)) {
                    // `cxx:398-403`: avoid adding it twice.
                    co.remove_item_ori(CLOSING, &ed);
                    co.add_item(CLOSING, &ed);
                } else {
                    // `cxx:404-408`: `iov != 0`
                    let iov = ori_in_sor(&v, &ed, false);
                    co.add_item(iov, &ed);
                }
            }
        }

        // `cxx:454-458`
        let n_v = self.my_map_v_eds.len();
        if n_v < 1 {
            return false;
        }

        // `cxx:460-479`
        let v_keys = self.my_v_order.clone();
        for key in v_keys {
            let Some(co) = self.my_map_v_eds.get(&key) else {
                continue;
            };
            if co.is_faulty() {
                return false;
            }
            if co.is_multiple() && self.my_map_v_multiple.insert(key) {
                if let Some(v) = co.key().cloned() {
                    self.my_list_v_multiple.push(v);
                }
            }
        }
        true
    }
}
