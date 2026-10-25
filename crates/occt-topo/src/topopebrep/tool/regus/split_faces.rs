//! `TopOpeBRepTool_REGUS` face splitting.
//!
//! Source (OCCT 8.0.0, `TopOpeBRepTool_REGUS.cxx`):
//! - `WireToFace(Fanc, nWs, nFs)` (`cxx:230-258`)
//! - `SplitF(Fanc, FSplits)` (`cxx:262-397`)
//! - `SplitFaces()` (`cxx:400-471`)

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key};
use crate::topopebrep::tool::classi::Classi;
use crate::topopebrep::tool::closed::{ori_in_sor, INTERNAL};
use crate::topopebrep::tool::mk_shell::remove_from_list;
use crate::topopebrep::tool::reguw::Reguw;
use crate::topopebrep::tool::wire_to_faces::wire_to_face;

use super::Regus;

impl Regus {
    /// `WireToFace(Fanc, nWs, nFs)` (`cxx:230-258`): classify the wires, then
    /// rebuild one face per boundary wire through `TOOL::WireToFace`.
    ///
    /// `cxx:247-256` returns `false` when the 2d classification fails;
    /// `TOOL::WireToFace` itself always answers `true` (`TOOL.cxx:1972`).
    fn wire_to_faces(&self, f_anc: &TopoShape, n_ws: &[TopoShape]) -> Option<Vec<TopoShape>> {
        let mut classi = Classi::new();
        let fwd = f_anc.clone().oriented(Orientation::Forward);
        classi.init_2d(&fwd);

        let mut map_wlow = std::collections::HashMap::new();
        if !classi.classilist(n_ws, &mut map_wlow) {
            return None;
        }

        // `WireToFace` walks the `NCollection_DataMap`; the port keeps the
        // input wire order so the faces come out in a stable order.
        let pairs: Vec<(TopoShape, Vec<TopoShape>)> = n_ws
            .iter()
            .filter_map(|w| {
                map_wlow
                    .get(&shape_map_key(w))
                    .map(|holes| (w.clone(), holes.clone()))
            })
            .collect();

        let mut n_fs = Vec::new();
        wire_to_face(f_anc, &pairs, &mut n_fs);
        Some(n_fs)
    }

    /// `SplitF(Fanc, FSplits)` (`cxx:262-397`): regularize the INTERNAL edges
    /// of each wire of `f_anc` and rebuild the face on the resulting wires.
    ///
    /// `None` covers the four `return false` arms: a `REGUW` connexity update
    /// (`cxx:330-355`), `REGUW::REGU` (`cxx:365-369`) and the wire rebuild
    /// (`cxx:383-387`).
    pub fn split_f(&self, f_anc: &TopoShape) -> Option<Vec<TopoShape>> {
        let o_anc = f_anc.orientation();
        let a_face = f_anc.clone().oriented(Orientation::Forward);

        let mut reguw = Reguw::new(&a_face);
        let mut n_ws: Vec<TopoShape> = Vec::new();
        let mut has_sp = false;

        // `cxx:277-374`: one wire at a time.
        for w in explorer_shapes(&a_face, ShapeType::Wire) {
            reguw.init(&w);
            let _ = reguw.map_s();

            // `cxx:284-293`: the INTERNAL edges of `w`.
            let mut e_is: Vec<TopoShape> = Vec::new();
            for e in explorer_shapes(&w, ShapeType::Edge) {
                if e.orientation() == Orientation::Internal {
                    e_is.push(e);
                }
            }

            // `cxx:298-358`: every INTERNAL edge with a multiple connexity at
            // both bounds is replaced by its two oriented occurrences.
            let mut i = 0usize;
            while i < e_is.len() {
                let e_i = e_is[i].clone();
                let e_fwd = e_i.clone().oriented(Orientation::Forward);
                let (vf, vl) = edge_vertices(&crate::shape::Edge(e_fwd));
                let nef = match vf.as_ref().and_then(|v| reguw.connexity(&v.0)) {
                    Some(co) => {
                        let mut lef = Vec::new();
                        co.all_items(&mut lef)
                    }
                    None => 0,
                };
                let nel = match vl.as_ref().and_then(|v| reguw.connexity(&v.0)) {
                    Some(co) => {
                        let mut lel = Vec::new();
                        co.all_items(&mut lel)
                    }
                    None => 0,
                };
                if nef <= 1 || nel <= 1 {
                    // `cxx:315-318`: `eIs.Remove(ite)`, the iterator stays on
                    // the following element.
                    e_is.remove(i);
                    continue;
                }

                // `cxx:321-327`: the reversed and forward occurrences.
                let e_r = e_i.clone().oriented(Orientation::Reversed);
                let e_f = e_i.clone().oriented(Orientation::Forward);

                // `cxx:329-355`.
                for v in explorer_shapes(&e_i, ShapeType::Vertex) {
                    if !reguw.remove_old_connexity(&v, INTERNAL, &e_i) {
                        return None;
                    }
                    let iv_f = ori_in_sor(&v, &e_f, false);
                    if !reguw.add_new_connexity(&v, iv_f, &e_f) {
                        return None;
                    }
                    let iv_r = ori_in_sor(&v, &e_r, false);
                    if !reguw.add_new_connexity(&v, iv_r, &e_r) {
                        return None;
                    }
                    if !reguw.update_multiple(&v) {
                        return None;
                    }
                }
                i += 1;
            }

            // `cxx:363-373`: `REGU()` then `GetSplits`.
            if !reguw.regu_all() {
                return None;
            }
            let mut sp_w: Vec<TopoShape> = Vec::new();
            if reguw.get_splits(&mut sp_w) && !sp_w.is_empty() {
                n_ws.extend(sp_w);
                has_sp = true;
            }
        }

        // `cxx:377-381`.
        if !has_sp {
            return None;
        }
        let n_fs = self.wire_to_faces(&a_face, &n_ws)?;
        // `cxx:391-395`: every split is restored to the ancestor orientation.
        Some(n_fs.into_iter().map(|f| f.oriented(o_anc)).collect())
    }

    /// `SplitFaces()` (`cxx:400-471`): split every face of `myS` and keep
    /// `mymapeFs` / `mymapemult` in step with the splits.
    pub fn split_faces(&mut self) -> bool {
        let Some(s) = self.my_s.clone() else {
            return true;
        };
        for f in explorer_shapes(&s, ShapeType::Face) {
            let Some(lfsp) = self.split_f(&f) else {
                continue;
            };
            self.my_f_splits.insert(shape_map_key(&f), lfsp.clone());

            // `cxx:426-457`: one ancestor face becomes `lfsp.len()`.
            self.my_n_f -= 1;
            for fsp in &lfsp {
                self.my_n_f += 1;
                for e in explorer_shapes(fsp, ShapeType::Edge) {
                    if !self.is_bound(&e) {
                        return false;
                    }
                    if let Some(lof) = self.change_find(&e) {
                        remove_from_list(lof, &f);
                        lof.push(fsp.clone());
                        if lof.len() > 2 {
                            self.my_mape_mult.insert(shape_map_key(&e));
                        }
                    }
                }
            }
        }
        true
    }
}
