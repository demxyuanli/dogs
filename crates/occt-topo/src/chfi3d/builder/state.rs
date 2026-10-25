//! `ChFi3d_Builder` state transitions: `Reset`, `Generated`, `UpdateTolesp`
//! and the `myEVIMap` / `myEdgeFirstFace` maintenance.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.cxx:924-977`,
//! `ChFi3d_Builder_C2.cxx:814-824`, `ChFi3d_Builder_2.cxx:3801-3806`,
//! `ChFi3d_Builder_1.cxx:927` / `:1040` / `:1089` / `:1146` (OCCT 8.0.0).
//! Split out of `builder/mod.rs` together with the accessors.

use super::{spine_of, ChFi3dBuilder};
use crate::abs::ShapeType;
use crate::chfi_ds::StripeHandle;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;

impl ChFi3dBuilder {
    /// `Reset()` (`ChFi3d_Builder.cxx:924-948`): the state reset of the
    /// acquisition, dropping the stripes that lost their spine.
    pub fn reset(&mut self) {
        self.done = false;
        self.my_v_data_map.clear();
        self.my_regul.clear();
        self.my_evi_map.clear();
        self.bad_stripes.clear();
        self.bad_vertices.clear();

        let mut i = 0;
        while i < self.my_list_stripe.len() {
            if spine_of(&self.my_list_stripe[i]).is_some() {
                if let Ok(mut guard) = self.my_list_stripe[i].write() {
                    guard.reset();
                }
                i += 1;
            } else {
                self.my_list_stripe.remove(i);
            }
        }
    }

    /// `Generated(EouV)` (`ChFi3d_Builder.cxx:950-977`): the faces built by the
    /// topological operation for an edge or a vertex, through `myEVIMap`.
    ///
    /// `e_or_v` is an `Option` because a null `TopoDS_Shape` handle has no
    /// counterpart in the port.
    pub fn generated(&mut self, e_or_v: Option<&TopoShape>) -> &[TopoShape] {
        self.my_generated.clear();
        let Some(eou_v) = e_or_v else {
            return &self.my_generated;
        };
        if eou_v.shape_type() != ShapeType::Edge && eou_v.shape_type() != ShapeType::Vertex {
            return &self.my_generated;
        }
        let Some(list) = self.evi_indices(eou_v).map(|l| l.to_vec()) else {
            return &self.my_generated;
        };
        for i in list {
            if let Some(faces) = self.my_coup.new_faces(i) {
                self.my_generated.extend_from_slice(faces);
            }
        }
        &self.my_generated
    }

    /// `UpdateTolesp()` (`ChFi3d_Builder_C2.cxx:814-824`): lowers `tolesp` to
    /// the `tolesp` of every stripe's guideline.
    pub fn update_tolesp(&mut self) {
        for stripe in &self.my_list_stripe {
            let Some(spine) = spine_of(stripe) else {
                continue;
            };
            let Ok(guard) = spine.read() else {
                continue;
            };
            let current_stripe_tolesp = guard.tolesp();
            if self.tolesp > current_stripe_tolesp {
                self.tolesp = current_stripe_tolesp;
            }
        }
    }

    /// `IsBound(E)` / `Bind(E, L)` / `ChangeFind(E)` of `myEVIMap`
    /// (`ChFi3d_Builder_2.cxx:3801-3806`, `ChFi3d_Builder_C1.cxx:2817-2834`):
    /// the index list of `s`, created empty on first sight of `s`.
    pub fn change_evi(&mut self, s: &TopoShape) -> &mut Vec<i32> {
        if let Some(k) = self
            .my_evi_map
            .iter()
            .position(|(key, _)| is_same(key, s))
        {
            return &mut self.my_evi_map[k].1;
        }
        self.my_evi_map.push((s.clone(), Vec::new()));
        let last = self.my_evi_map.len() - 1;
        &mut self.my_evi_map[last].1
    }

    /// `myEdgeFirstFace.Bind(Ec, Face)` (`ChFi3d_Builder_1.cxx:927`, `:1040`,
    /// `:1146`).
    pub fn bind_edge_first_face(&mut self, e: &TopoShape, f: &TopoShape) {
        if let Some(k) = self
            .my_edge_first_face
            .iter()
            .position(|(key, _)| is_same(key, e))
        {
            self.my_edge_first_face[k].1 = f.clone();
            return;
        }
        self.my_edge_first_face.push((e.clone(), f.clone()));
    }

    /// `myEdgeFirstFace(Ec)` (`ChFi3d_Builder_1.cxx:1089`,
    /// `ChFi3d_Builder_2.cxx:836`, `:879`).
    pub fn edge_first_face(&self, e: &TopoShape) -> Option<&TopoShape> {
        self.my_edge_first_face
            .iter()
            .find(|(key, _)| is_same(key, e))
            .map(|(_, f)| f)
    }

    /// The `I`-th `myListStripe` entry, `I` `1`-based.
    pub(super) fn stripe_at(&self, ic: usize) -> Option<StripeHandle> {
        self.my_list_stripe.get(ic.checked_sub(1)?).cloned()
    }

    /// `myEVIMap.Find(EouV)`, without creating the entry.
    fn evi_indices(&self, s: &TopoShape) -> Option<&[i32]> {
        self.my_evi_map
            .iter()
            .find(|(key, _)| is_same(key, s))
            .map(|(_, list)| list.as_slice())
    }
}
