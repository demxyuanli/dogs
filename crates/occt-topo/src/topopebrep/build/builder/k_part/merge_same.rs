//! `TopOpeBRepBuild_Builder` — the `isfafa` and `issoso` KPart cases.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! (OCCT 8.0.0):
//! - `MergeKPartisfafa` (`:1034-1101`),
//! - `MergeKPartissoso` (`:1105-1192`).
//!
//! Both find the same-domain partners of the first face / solid of each
//! operand, pick the `GTopo` of the operation, map the two operands and hand
//! the lists to `GMergeFaces` / `GMergeSolids`. Those two are UNPORTED
//! (`TopOpeBRepBuild_GridFF.cxx:335-463`, `GridSS.cxx:282-349`): the port stops
//! right after `GMapShapes`, which is where OCCT's merge would leave the merged
//! lists empty, and the compound bookkeeping below still runs.
//!
//! `L1 = ChangeMerged(F1, myState1)` (`cxx:1091`, `:1179`) assigns the
//! *reference* in OCCT, so the compound and `F1` share one `NCollection_List`;
//! the port's map holds a `Vec` per shape, so the two lists are made equal by
//! copying the content, which is what every later reader of the map sees.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::g_tool::GTool;

use super::super::k_part::{KP_ISAFA, KP_ISSOSO};
use super::Builder;

/// `GFindSamDom(S, L1, L2)` (`Grid.cxx:182-190`): `L1 = {S}` then
/// `GFindSamDom(L1, L2)` (`Grid.cxx:197-258`).
fn g_find_sam_dom(b: &Builder, s: &TopoShape, l1: &mut Vec<TopoShape>, l2: &mut Vec<TopoShape>) {
    l1.clear();
    l2.clear();
    l1.push(s.clone());
    b.find_same_domain(l1, l2);
}

impl Builder {
    /// `MergeKPartisfafa()` (`KPart.cxx:1034-1101`).
    pub fn merge_k_part_isfafa(&mut self) {
        // cxx:1042
        if self.my_is_k_part != KP_ISAFA {
            return;
        }
        let (Some(shape1), Some(shape2)) = (self.my_shape1.clone(), self.my_shape2.clone()) else {
            return;
        };
        let st1 = self.my_state1;
        let st2 = self.my_state2;

        // cxx:1044-1057
        let Some(f1) = explorer_shapes(&shape1, ShapeType::Face).first().cloned() else {
            return;
        };
        let Some(f2) = explorer_shapes(&shape2, ShapeType::Face).first().cloned() else {
            return;
        };

        // cxx:1059-1060
        let mut lf1: Vec<TopoShape> = Vec::new();
        let mut lf2: Vec<TopoShape> = Vec::new();
        g_find_sam_dom(self, &f1, &mut lf1, &mut lf2);

        // cxx:1062-1083
        let tf = ShapeType::Face;
        let g = if self.opec12() {
            GTool::g_cut_same(tf, tf)
        } else if self.opec21() {
            GTool::g_cut_same(tf, tf).copy_permuted()
        } else if self.opecom() {
            GTool::g_com_same(tf, tf)
        } else if self.opefus() {
            GTool::g_fus_same(tf, tf)
        } else {
            return;
        };

        // cxx:1085-1086
        self.map_shapes(Some(&shape1), Some(&shape2));
        // `GMergeFaces(lf1, lf2, g)` (`GridFF.cxx:335-463`) is UNPORTED.
        let _ = (&g, &lf2);

        // cxx:1088-1095
        if shape1.shape_type() == ShapeType::Compound {
            let l1 = self.change_merged(&f1, st1).clone();
            *self.change_merged(&shape1, st1) = l1;
        }

        // cxx:1097-1101
        if shape2.shape_type() == ShapeType::Compound {
            let l2 = self.change_merged(&f2, st2).clone();
            *self.change_merged(&shape2, st2) = l2;
        }
    }

    /// `MergeKPartissoso()` (`KPart.cxx:1105-1192`).
    pub fn merge_k_part_issoso(&mut self) {
        // cxx:1112
        if self.my_is_k_part != KP_ISSOSO {
            return;
        }
        let shape1 = self.my_shape1.clone();
        let shape2 = self.my_shape2.clone();
        let st1 = self.my_state1;
        let st2 = self.my_state2;

        // cxx:1115-1127
        let mut so1: Option<TopoShape> = None;
        if let Some(shape1) = &shape1 {
            match explorer_shapes(shape1, ShapeType::Solid).first() {
                Some(s) => so1 = Some(s.clone()),
                None => return,
            }
        }

        // cxx:1129-1138
        let mut so2: Option<TopoShape> = None;
        if let Some(shape2) = &shape2 {
            match explorer_shapes(shape2, ShapeType::Solid).first() {
                Some(s) => so2 = Some(s.clone()),
                None => return,
            }
        }

        // cxx:1140-1143
        let Some(so1) = so1 else {
            return;
        };

        // cxx:1145-1146
        let mut lso1: Vec<TopoShape> = Vec::new();
        let mut lso2: Vec<TopoShape> = Vec::new();
        g_find_sam_dom(self, &so1, &mut lso1, &mut lso2);

        // cxx:1148-1169
        let tf = ShapeType::Face; // NYI TopAbs_SOLID
        let g = if self.opec12() {
            GTool::g_cut_same(tf, tf)
        } else if self.opec21() {
            GTool::g_cut_same(tf, tf).copy_permuted()
        } else if self.opecom() {
            GTool::g_com_same(tf, tf)
        } else if self.opefus() {
            GTool::g_fus_same(tf, tf)
        } else {
            return;
        };

        // cxx:1171-1172
        self.map_shapes(shape1.as_ref(), shape2.as_ref());
        // `GMergeSolids(lso1, lso2, g)` (`GridSS.cxx:282-349`) is UNPORTED.
        let _ = (&g, &lso2);

        // cxx:1174-1181
        if let Some(shape1) = &shape1 {
            if shape1.shape_type() == ShapeType::Compound {
                let l1 = self.change_merged(&so1, st1).clone();
                *self.change_merged(shape1, st1) = l1;
            }
        }

        // cxx:1183-1191
        if let (Some(shape2), Some(so2)) = (&shape2, &so2) {
            if shape2.shape_type() == ShapeType::Compound {
                let l2 = self.change_merged(so2, st2).clone();
                *self.change_merged(shape2, st2) = l2;
            }
        }
    }
}
