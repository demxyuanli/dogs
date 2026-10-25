//! `TopOpeBRepBuild_Builder` — the `iskole` KPart predicate and its result
//! table.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_kpkole.cxx`
//! (OCCT 8.0.0):
//! - `KPiskolesh` (`:1085-1157`), the per-operand shape test,
//! - `KPiskole` (`:576-732`), which fills `myKPMAPf1f2` and the DS section
//!   edges,
//! - `KPiskoleanalyse` (`:734-843`), the `(ires, icla1, icla2)` table.
//!
//! `FUNKP_KPiskolesh` (`:1033-1075`) already lives in `predicates.rs` as
//! `fun_kp_iskolesh`; `KPiskolesh` is the wrapper that adds the planar /
//! single-wire / no-edge-geometry constraints.

use crate::abs::{ShapeType, State};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::tool::geometry::fun_tool_plane;

use super::super::k_part::{
    RESFACE1, RESFACE2, RESNEWSHE, RESNULL, RESSHAPE1, RESSHAPE2, RESUNDEF, SHEAUCU,
    SHECLASAUTR, SHEGARDAUTR, SHEUNDEF,
};
use super::Builder;

impl Builder {
    /// `KPiskolesh(Sarg, lShsd, lfhsd)` (`kpkole.cxx:1085-1157`).
    ///
    /// `Sarg` is a shape the glue case handles: one same-domain solid whose
    /// faces are all planar, each with a single wire and no edge geometry.
    pub fn kp_iskolesh(
        &self,
        sarg: &TopoShape,
        lshsd: &mut Vec<TopoShape>,
        lfhsd: &mut Vec<TopoShape>,
    ) -> bool {
        // kpkole.cxx:1094-1097
        if !self.fun_kp_iskolesh(sarg, lshsd) {
            return false;
        }

        // kpkole.cxx:1101-1102
        self.kp_lhsd_into(sarg, ShapeType::Face, lfhsd);
        for fac in lfhsd.clone() {
            let Some(face) = Face::wrap(fac.clone()) else {
                return false;
            };
            // kpkole.cxx:1107-1111
            if !fun_tool_plane(&face) {
                return false;
            }
            // kpkole.cxx:1113-1116
            if self.kp_ls(&fac, ShapeType::Wire) > 1 {
                return false;
            }
            // kpkole.cxx:1122-1126
            if self.kp_lhg(&fac, ShapeType::Edge) != 0 {
                return false;
            }
        }
        // kpkole.cxx:1156
        true
    }

    /// `KPiskole()` (`kpkole.cxx:576-732`).
    ///
    /// Fills `myKPMAPf1f2` with the `(OUT face, {IN faces})` pairs of the glue
    /// configuration and registers the IN faces' edges as DS section edges.
    pub fn kp_iskole(&mut self) -> i32 {
        // kpkole.cxx:583-591
        let Some(shape1) = self.my_shape1.clone() else {
            return 0;
        };
        let mut lshsd1: Vec<TopoShape> = Vec::new();
        let mut lfhsd1: Vec<TopoShape> = Vec::new();
        if !self.kp_iskolesh(&shape1, &mut lshsd1, &mut lfhsd1) {
            return 0;
        }
        if lfhsd1.is_empty() {
            return 0;
        }

        // kpkole.cxx:595-602
        let Some(shape2) = self.my_shape2.clone() else {
            return 0;
        };
        let mut lshsd2: Vec<TopoShape> = Vec::new();
        let mut lfhsd2: Vec<TopoShape> = Vec::new();
        if !self.kp_iskolesh(&shape2, &mut lshsd2, &mut lfhsd2) {
            return 0;
        }
        if lfhsd2.is_empty() {
            return 0;
        }

        // kpkole.cxx:611-615: several solids on either side fall back to the
        // general case.
        if lshsd1.len() > 1 || lshsd2.len() > 1 {
            return 0;
        }

        // kpkole.cxx:617-720
        let mut les: Vec<TopoShape> = Vec::new();
        for f1 in lfhsd1.clone() {
            let mut lf1 = vec![f1.clone()];
            let mut lf2: Vec<TopoShape> = Vec::new();
            self.kp_same_domain(&mut lf1, &mut lf2);

            for f2 in lf2 {
                // kpkole.cxx:656-660
                let Some((classok, state1, state2)) = self.kp_iskole_ff_states(&f1, &f2) else {
                    return 0;
                };
                if !classok {
                    return 0;
                }

                // kpkole.cxx:663-697: rebuild the OUT face keyed on the other
                // side's IN faces.
                if state1 == State::Out && state2 == State::In {
                    self.kp_map_f1f2_bind(&f1);
                    self.kp_map_f1f2_append(&f1, &f2);
                    self.kp_map_f1f2_bind(&f2);
                } else if state2 == State::Out && state1 == State::In {
                    self.kp_map_f1f2_bind(&f2);
                    self.kp_map_f1f2_append(&f2, &f1);
                    self.kp_map_f1f2_bind(&f1);
                }

                // kpkole.cxx:699-718: the IN face's edges are section edges.
                let fw = if state1 == State::In {
                    Some(f1.clone())
                } else if state2 == State::In {
                    Some(f2.clone())
                } else {
                    None
                };
                let Some(fw) = fw else {
                    continue;
                };
                les.extend(explorer_shapes(&fw, ShapeType::Edge));
            }
        }

        // kpkole.cxx:722-729
        let section_edges: Vec<TopoShape> = les
            .into_iter()
            .filter(|e| e.shape_type() == ShapeType::Edge)
            .collect();
        self.kp_set_section_edges(&section_edges);

        // kpkole.cxx:731
        1
    }

    /// `KPiskoleanalyse(Stfac1, Stfac2, Stsol1, Stsol2, ires, icla1, icla2)`
    /// (`kpkole.cxx:734-843`).
    pub fn kp_iskole_analyse(
        &self,
        stfac1: State,
        stfac2: State,
        stsol1: State,
        stsol2: State,
    ) -> (i32, i32, i32) {
        let mut ires = RESUNDEF;
        let mut icla1 = SHEUNDEF;
        let mut icla2 = SHEUNDEF;

        if self.opefus() {
            // kpkole.cxx:748-764
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESNEWSHE;
                icla1 = SHEGARDAUTR;
                icla2 = SHEGARDAUTR;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESSHAPE1;
                icla1 = SHECLASAUTR;
                icla2 = SHEAUCU;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESSHAPE2;
                icla1 = SHEAUCU;
                icla2 = SHECLASAUTR;
            }
        } else if self.opec12() {
            // kpkole.cxx:769-786
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE1;
                icla1 = SHEGARDAUTR;
                icla2 = SHEAUCU;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNEWSHE;
                icla1 = SHECLASAUTR;
                icla2 = SHEAUCU;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNULL;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            }
        } else if self.opec21() {
            // kpkole.cxx:789-806
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE2;
                icla1 = SHEAUCU;
                icla2 = SHEGARDAUTR;
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNULL;
                icla1 = SHEAUCU;
                icla2 = SHEAUCU;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESNEWSHE;
                icla1 = SHEAUCU;
                icla2 = SHECLASAUTR;
            }
        } else if self.opecom() {
            // kpkole.cxx:809-834
            if stsol1 == State::Out && stsol2 == State::Out {
                if stfac1 == State::In {
                    ires = RESFACE1;
                    icla1 = SHEAUCU;
                    icla2 = SHEAUCU;
                }
                if stfac2 == State::In {
                    ires = RESFACE2;
                    icla1 = SHEAUCU;
                    icla2 = SHEAUCU;
                }
            } else if stsol1 == State::Out && stsol2 == State::In {
                ires = RESSHAPE2;
                icla1 = SHECLASAUTR;
                icla2 = SHEGARDAUTR;
            } else if stsol1 == State::In && stsol2 == State::Out {
                ires = RESSHAPE1;
                icla1 = SHEGARDAUTR;
                icla2 = SHECLASAUTR;
            }
        }
        (ires, icla1, icla2)
    }
}
