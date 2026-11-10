//! `TopOpeBRepBuild_Builder` — the `MergeShapes` drivers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Merge.cxx:174-336`
//! (`MergeShapes`), `:356-361` (`MergeSolids`), `:368-371` (`MergeSolid`)
//! (OCCT 8.0.0).
//!
//! `MergeShapes` is the top of the merge side of the class: it records the two
//! operands and their states, maps their sub-shapes (which is what
//! `ShapeRank` answers from), splits every shape of each operand into a single
//! `ShellFaceSet`, and — when an operand is a `COMPOUND` — makes the compound
//! carry the merged list of the first sub-shape it explored.
//!
//! UNPORTED, with its OCCT location:
//! - the `IsKPart()` branch (`cxx:203-208`): `FindIsKPart`
//!   (`KPart.cxx:111-166`) and `MergeKPart` (`KPart.cxx:176-220`) need the five
//!   `KP*` predicates, none of which is translated. `myIsKPart` therefore stays
//!   `0` and the branch is not taken, which is the OCCT answer when no KPart
//!   case is recognised.

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::is_equal;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

use super::super::shell_face_set::ShellFaceSet;
use super::Builder;

impl Builder {
    /// `MergeShapes(S1, ToBuild1, S2, ToBuild2)` (`Merge.cxx:174-336`).
    ///
    /// A null operand is `Option::None`; OCCT keeps it in `myShape1` /
    /// `myShape2` and skips its exploration (`cxx:227`, `:267`).
    pub fn merge_shapes(
        &mut self,
        s1: Option<&TopoShape>,
        to_build1: State,
        s2: Option<&TopoShape>,
        to_build2: State,
    ) {
        // cxx:179-187: `S1.IsEqual(S2)` returns immediately; two null shapes
        // are equal in OCCT (`TopoDS_Shape.lxx:43-52`).
        let les_memes = match (s1, s2) {
            (None, None) => true,
            (Some(a), Some(b)) => is_equal(a, b),
            _ => false,
        };
        if les_memes {
            return;
        }

        // cxx:193-196
        self.my_state1 = to_build1;
        self.my_state2 = to_build2;
        self.my_shape1 = s1.cloned();
        self.my_shape2 = s2.cloned();
        let _s1null = s1.is_none();
        let _s2null = s2.is_none();

        // cxx:199-200
        self.map_shapes(s1, s2);
        self.split_section_edges();

        // cxx:203-208
        if self.is_k_part() != 0 {
            self.merge_k_part();
            self.clear_maps();
            return;
        }

        // cxx:211-212
        let rev_ori1 = Self::reverse(to_build1, to_build2);
        let rev_ori2 = Self::reverse(to_build2, to_build1);

        // cxx:216
        let mut sfs = ShellFaceSet::new();

        // ----- S1 (`cxx:225-261`) -----
        let mut t1 = ShapeType::Compound;
        let mut tex1 = ShapeType::Compound;
        if let Some(shape1) = s1 {
            t1 = Self::top_type(shape1);
            tex1 = t1;
            if t1 == ShapeType::Compound {
                // cxx:232-250: the deepest of SOLID / SHELL / FACE / EDGE that
                // the compound actually holds.
                tex1 = ShapeType::Solid;
                let mut ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                if !ex1.more() {
                    tex1 = ShapeType::Shell;
                    ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                    if !ex1.more() {
                        tex1 = ShapeType::Face;
                        ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                        if !ex1.more() {
                            tex1 = ShapeType::Edge;
                            ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                        }
                    }
                }
                self.split_shapes(&mut ex1, to_build1, to_build2, &mut sfs, rev_ori1);
            } else if t1 == ShapeType::Wire {
                // cxx:251-255
                tex1 = ShapeType::Edge;
                let mut ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                self.split_shapes(&mut ex1, to_build1, to_build2, &mut sfs, rev_ori1);
            } else {
                // cxx:256-259
                let mut ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                self.split_shapes(&mut ex1, to_build1, to_build2, &mut sfs, rev_ori1);
            }
        }

        // ----- S2 (`cxx:263-300`) -----
        let mut t2 = ShapeType::Compound;
        let mut tex2 = ShapeType::Compound;
        if let Some(shape2) = s2 {
            t2 = Self::top_type(shape2);
            tex2 = t2;
            if t2 == ShapeType::Compound {
                tex2 = ShapeType::Solid;
                let mut ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                if !ex2.more() {
                    tex2 = ShapeType::Shell;
                    ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                    if !ex2.more() {
                        tex2 = ShapeType::Face;
                        ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                        if !ex2.more() {
                            tex2 = ShapeType::Edge;
                            ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                        }
                    }
                }
                self.split_shapes(&mut ex2, to_build2, to_build1, &mut sfs, rev_ori2);
            } else if t2 == ShapeType::Wire {
                tex2 = ShapeType::Edge;
                let mut ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                self.split_shapes(&mut ex2, to_build2, to_build1, &mut sfs, rev_ori2);
            } else {
                let mut ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                self.split_shapes(&mut ex2, to_build2, to_build1, &mut sfs, rev_ori2);
            }
        }

        // ----- `cxx:302-333`: a COMPOUND takes the merged list of the first
        // sub-shape of the type that was explored. -----
        if let Some(shape1) = s1 {
            if t1 == ShapeType::Compound {
                let ex1 = ShapeExplorer::with_shape(shape1, tex1, ShapeType::Shape);
                if ex1.more() {
                    if let Some(ss1) = ex1.current().cloned() {
                        let lss1 = self.merged(&ss1, to_build1).to_vec();
                        *self.change_merged(shape1, to_build1) = lss1;
                    }
                }
            }
        }
        if let Some(shape2) = s2 {
            if t2 == ShapeType::Compound {
                let ex2 = ShapeExplorer::with_shape(shape2, tex2, ShapeType::Shape);
                if ex2.more() {
                    if let Some(ss2) = ex2.current().cloned() {
                        let lss2 = self.merged(&ss2, to_build2).to_vec();
                        *self.change_merged(shape2, to_build2) = lss2;
                    }
                }
            }
        }

        // cxx:335
        self.clear_maps();
    }

    /// `MergeSolids(S1, ToBuild1, S2, ToBuild2)` (`Merge.cxx:356-361`).
    pub fn merge_solids(
        &mut self,
        s1: &TopoShape,
        to_build1: State,
        s2: &TopoShape,
        to_build2: State,
    ) {
        self.merge_shapes(Some(s1), to_build1, Some(s2), to_build2);
    }

    /// `MergeSolid(S, ToBuild)` (`Merge.cxx:368-371`): a null `S2`.
    pub fn merge_solid(&mut self, s: &TopoShape, to_build: State) {
        self.merge_shapes(Some(s), to_build, None, to_build);
    }
}
