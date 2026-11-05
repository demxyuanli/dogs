//! `ChFi3d_Builder::PerformFilletOnVertex` - dispatch one vertex corner.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.cxx:759-920`
//! (OCCT 8.0.0).
//!
//! The stripe list of the vertex is walked first to answer two questions:
//! whether every stripe's end `SurfData` is degenerated (`nondegenere`) and
//! whether none of them is (`toujoursdegenere`). `i` is left at the number of
//! stripes and `sp` / `isfirst` at the values of the **last** stripe - the
//! `case 1:` arm then reads `sp->Status(isfirst)` of that last stripe, which
//! is the OCCT behaviour and is kept.
//!
//! Dispatch:
//!
//! | `i` | `nba` | call |
//! |---|---|---|
//! | 1 | - | `PerformIntersectionAtEnd` when `nba > 3`, else `MoreSurfdata` / `PerformMoreSurfdata` / `PerformOneCorner` |
//! | 2 | `> 3` | `PerformMoreThreeCorner` |
//! | 2 | `<= 3` | `PerformTwoCorner` |
//! | 3 | `> 3` | `PerformMoreThreeCorner` |
//! | 3 | `<= 3` | `PerformThreeCorner` |
//! | other | - | `PerformMoreThreeCorner` |
//!
//! `nba` is `ChFi3d_NumberOfSharpEdges(Vtx, myVEMap, myEFMap)`; the commented
//! face-count block of `cxx:793-813` is not translated. `PerformIntersectionAtEnd`
//! is the base body itself: its three internal exits are dispatched here from the
//! returned [`EndIntersectionExit`] (`PerformMoreSurfdata` / `IntersectMoreCorner`
//! / `PerformMoreThreeCorner` are all base members, the first two fully ported
//! and the third up to `cxx:1581` — see `perform_more_three_corner`). The
//! subclass bodies of the other corners (`PerformTwoCorner` /
//! `PerformThreeCorner` / `Extent*Corner`) live on [`ChFi3dCornerOps`], so the
//! caller supplies them; `Extent*Corner` and `SetRegul` are ported for the
//! fillet builder (`chfi3d/builder/fil_builder/extent_corner.rs` /
//! `set_regul.rs`), the two `Perform*Corner` bodies are not.

use crate::chfi3d::builder_helpers::number_of_sharp_edges;
use crate::chfi_ds::State;
use crate::shape::Vertex;

use super::{ChFi3dBuilder, ChFi3dCornerOps, EndIntersectionExit};

impl ChFi3dBuilder {
    /// `PerformFilletOnVertex(Index)` (`cxx:759-920`).
    pub fn perform_fillet_on_vertex(
        &mut self,
        ops: &mut dyn ChFi3dCornerOps,
        index: usize,
    ) -> Result<(), String> {
        // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);`
        let vtx = self
            .my_v_data_map
            .find_key(index)
            .cloned()
            .ok_or_else(|| "PerformFilletOnVertex: myVDataMap index out of range".to_string())?;
        let stripes = self
            .my_v_data_map
            .find_from_index(index)
            .map(|l| l.to_vec())
            .unwrap_or_default();

        let mut nondegenere = true;
        let mut toujoursdegenere = true;
        // `bool isfirst = false;` (`cxx:771`): OCCT reassigns it at `cxx:779`
        // before any read, so the initializer is never used.
        #[allow(unused_assignments)]
        let mut isfirst = false;
        let mut sp = None;
        let mut count = 0usize;

        for stripe in &stripes {
            // `sp = stripe->Spine();`
            sp = stripe.read().ok().and_then(|g| g.spine().cloned());

            // `num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
            // isfirst = (sens == 1); Fd = ...Sequence().Value(num);`
            let mut sens = 0;
            let num = crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data(
                &vtx, stripe, &mut sens,
            )?;
            isfirst = sens == 1;
            let fd = {
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformFilletOnVertex: poisoned stripe".to_string())?;
                let idx = usize::try_from(num)
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or_else(|| "PerformFilletOnVertex: bad SurfData index".to_string())?;
                guard.set_of_surf_data().get(idx).cloned().ok_or_else(|| {
                    "PerformFilletOnVertex: SurfData index out of range".to_string()
                })?
            };

            // `const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1);`
            // `if (CV1.Point().IsEqual(CV2.Point(), 0)) nondegenere = false;
            //  else toujoursdegenere = false;` (`cxx:783-791`).
            let cv1 = fd.vertex(isfirst, 1);
            let cv2 = fd.vertex(isfirst, 2);
            if cv1.point().is_equal_with_tolerance(cv2.point(), 0.0) {
                nondegenere = false;
            } else {
                toujoursdegenere = false;
            }
            count += 1;
        }

        // `int nba = ChFi3d_NumberOfSharpEdges(Vtx, myVEMap, myEFMap);`
        let vertex = Vertex::wrap(vtx.clone())
            .ok_or_else(|| "PerformFilletOnVertex: myVDataMap key is not a vertex".to_string())?;
        let nba = number_of_sharp_edges(&vertex, &self.my_ve_map, &self.my_ef_map)?;

        if nondegenere {
            match count {
                1 => {
                    // `if (sp->Status(isfirst) == ChFiDS_FreeBoundary) return;`
                    let status = match &sp {
                        Some(sp) => Some(
                            sp.read()
                                .map_err(|_| {
                                    "PerformFilletOnVertex: poisoned spine".to_string()
                                })?
                                .status(isfirst),
                        ),
                        None => None,
                    };
                    if status == Some(State::FreeBoundary) {
                        return Ok(());
                    }
                    if nba > 3 {
                        // `PerformIntersectionAtEnd(Index);` (`cxx:826`).
                        self.dispatch_end_intersection_exit(index)?;
                    } else if self.more_surfdata(index) {
                        // `PerformMoreSurfdata(Index);` (`cxx:839-841`).
                        ops.perform_more_surfdata(self, index)?;
                    } else {
                        // `PerformOneCorner(Index);` (`cxx:846`).
                        ops.perform_one_corner(self, index, false)?;
                    }
                }
                2 => {
                    if nba > 3 {
                        self.perform_more_three_corner(index, count)?;
                    } else {
                        ops.perform_two_corner(self, index)?;
                    }
                }
                3 => {
                    if nba > 3 {
                        self.perform_more_three_corner(index, count)?;
                    } else {
                        ops.perform_three_corner(self, index)?;
                    }
                }
                n => {
                    self.perform_more_three_corner(index, n)?;
                }
            }
        } else if toujoursdegenere {
            // `PerformSingularCorner(Index);` (`cxx:912`).
            self.perform_singular_corner(index)?;
        } else {
            // `PerformMoreThreeCorner(Index, i); // Last chance...` (`cxx:917`).
            self.perform_more_three_corner(index, count)?;
        }
        Ok(())
    }

    /// `PerformIntersectionAtEnd(Index)` followed by the three exits that body
    /// takes (`cxx:2006-2010` → `PerformMoreSurfdata`, `cxx:2161-2165` →
    /// `PerformMoreThreeCorner(1)`, `cxx:2508-2509` → `IntersectMoreCorner`).
    ///
    /// Both OCCT callers of `PerformIntersectionAtEnd` dispatch those exits
    /// identically: `PerformFilletOnVertex` at `cxx:826` and `PerformOneCorner`
    /// at `cxx:821-825` / `cxx:910-914`. All three targets are non-virtual base
    /// members, so the inherent methods are called directly.
    pub(crate) fn dispatch_end_intersection_exit(&mut self, index: usize) -> Result<(), String> {
        match self.perform_intersection_at_end(index)? {
            EndIntersectionExit::Done => {}
            EndIntersectionExit::MoreSurfdata => {
                self.perform_more_surfdata(index)?;
            }
            EndIntersectionExit::MoreThreeCorner(n) => {
                self.perform_more_three_corner(index, n)?;
            }
            EndIntersectionExit::IntersectMoreCorner => {
                self.intersect_more_corner(index)?;
            }
        }
        Ok(())
    }
}
