//! `TopOpeBRepBuild_Builder` — the "keep" filters of `TopOpeBRepBuild_keep.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_keep.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GKeepShape(S, LSclass, TB)` | `41-47` |
//! | `GKeepShape1(S, LSclass, TB, pos)` | `49-114` |
//! | `GKeepShapes(S, LSclass, TB, Lin, Lou)` | `125-202` |
//!
//! A shape is classified against `LSclass` only when the list is not empty
//! (`cxx:61`); an empty list keeps everything. `GKeepShape1` additionally
//! lets `ChangeClassify` disable the classification of a face that the data
//! structure does not carry (`cxx:57-60`).
//!
//! UNPORTED: the `OCCT_DEBUG` trace block (`cxx:71-112`, `:137-197`) — dump
//! only, `S` is unused outside it.

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;

use super::Builder;

impl Builder {
    /// `GKeepShape(S, LSclass, TB)` (`cxx:41-47`).
    pub fn g_keep_shape(&mut self, s: &TopoShape, ls_class: &[TopoShape], tb: State) -> bool {
        let mut pos = State::Unknown;
        self.g_keep_shape1(s, ls_class, tb, &mut pos)
    }

    /// `GKeepShape1(S, LSclass, TB, pos)` (`cxx:49-114`).
    pub fn g_keep_shape1(
        &mut self,
        s: &TopoShape,
        ls_class: &[TopoShape],
        tb: State,
        pos: &mut State,
    ) -> bool {
        let mut keep = true;
        *pos = State::Unknown;

        // `cxx:58-61`: `myClassifyDef` false leaves `toclassify` true, which is
        // exactly what `Classify()` answers.
        let mut to_classify = true;
        if s.shape_type() == ShapeType::Face
            && !self
                .data_structure()
                .map(|h| h.has_shape(s, true))
                .unwrap_or(false)
        {
            to_classify = self.classify();
        }

        to_classify = to_classify && !ls_class.is_empty();
        if to_classify {
            *pos = self.shape_position(s, ls_class);
            if *pos != tb {
                keep = false;
            }
        }

        keep
    }

    /// `GKeepShapes(S, LSclass, TB, Lin, Lou)` (`cxx:125-202`): `Lou` is
    /// appended to, never cleared.
    pub fn g_keep_shapes(
        &mut self,
        _s: &TopoShape,
        ls_class: &[TopoShape],
        tb: State,
        l_in: &[TopoShape],
        l_out: &mut Vec<TopoShape>,
    ) {
        for sl in l_in {
            let mut keep = true;
            if !ls_class.is_empty() {
                let pos = self.shape_position(sl, ls_class);
                if pos != tb {
                    keep = false;
                }
            }
            if keep {
                l_out.push(sl.clone());
            }
        }
    }
}
