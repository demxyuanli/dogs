//! `TopOpeBRepTool_REGUS::MapS` (`cxx:164-227`).
//!
//! Builds `mymapeFs` (edge -> untouched ancestor faces), `mymapeFsstatic`
//! (edge -> every ancestor face of `myS`) and `mymapemult` (edges with more
//! than two ancestor faces).

use crate::abs::ShapeType;
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::Regus;

impl Regus {
    /// `MapS()` (`cxx:164-227`): always answers `true`.
    pub fn map_s(&mut self) -> bool {
        // `cxx:174-199`: `mymapeFs`, `mymapeFsstatic`, `myoldnF`.
        self.my_old_n_f = 0;
        let Some(s) = self.my_s.clone() else {
            return true;
        };
        for f in explorer_shapes(&s, ShapeType::Face) {
            self.my_old_n_f += 1;
            for e in explorer_shapes(&f, ShapeType::Edge) {
                self.bind_face_to_edge(&e, &f);
                let key = shape_map_key(&e);
                match self.my_mape_fs_static.get_mut(&key) {
                    Some(lof) => lof.push(f.clone()),
                    None => {
                        self.my_mape_fs_static.insert(key, vec![f.clone()]);
                    }
                }
            }
        }
        self.my_n_f = self.my_old_n_f;

        // `cxx:202-226`: `mymapemult`.
        for key in self.edge_keys() {
            let nf = self.my_mape_fs.get(&key).map(|l| l.len()).unwrap_or(0);
            if nf > 2 {
                self.my_mape_mult.insert(key);
            }
        }
        true
    }
}
