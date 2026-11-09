//! `Store` (`BRepOffset_Inter3d.cxx:1536-1554`), `SetDone` (`:1500-1516`) and
//! `IsDone` (`:1519-1533`) of `BRepOffset_Inter3d`.

use crate::shape::{Face, TopoShape};
use crate::topo_tools::{is_same, shape_map_key};

use super::BRepOffsetInter3d;

impl BRepOffsetInter3d {
    /// `Store(F1, F2, LInt1, LInt2)` (`cxx:1536-1554`): when `LInt1` is not
    /// empty, both faces become touched, `LInt1`/`LInt2` are added to their
    /// descendant lists and every `LInt1` entry becomes a new edge; the pair is
    /// marked done either way.
    pub(crate) fn store(
        &mut self,
        f1: &Face,
        f2: &Face,
        l_int1: &[TopoShape],
        l_int2: &[TopoShape],
    ) {
        if !l_int1.is_empty() {
            self.touched.add(&f1.0);
            self.touched.add(&f2.0);
            self.as_des.add_list(&f1.0, l_int1);
            self.as_des.add_list(&f2.0, l_int2);
            for e in l_int1 {
                self.new_edges.add(e);
            }
        }
        self.set_done(f1, f2);
    }

    /// `SetDone(F1, F2)` (`cxx:1500-1516`): binds the pair both ways.
    pub(crate) fn set_done(&mut self, f1: &Face, f2: &Face) {
        self.done
            .entry(shape_map_key(&f1.0))
            .or_default()
            .push(f2.0.clone());
        self.done
            .entry(shape_map_key(&f2.0))
            .or_default()
            .push(f1.0.clone());
    }

    /// `IsDone(F1, F2)` (`cxx:1519-1533`): `F2` occurs in `F1`'s list (the
    /// binding is symmetric, so the one direction is enough).
    pub(crate) fn is_done(&self, f1: &Face, f2: &Face) -> bool {
        match self.done.get(&shape_map_key(&f1.0)) {
            Some(list) => list.iter().any(|s| is_same(s, &f2.0)),
            None => false,
        }
    }
}
