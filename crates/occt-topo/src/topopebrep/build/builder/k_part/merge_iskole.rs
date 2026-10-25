//! `TopOpeBRepBuild_Builder` — `MergeKPartiskole`, the `iskole` (glued planar
//! faces) KPart case.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_kpkole.cxx`
//! (OCCT 8.0.0):
//! - `MergeKPartiskole` (`:56-572`),
//! - `FUNKP_KPmakefaces` (`:845-941`),
//! - `KPmakeface` (`:945-1031`).
//!
//! The `ires` / `icla1` / `icla2` values come from `KPiskoleanalyse`
//! (`:734-843`) in `iskole.rs`.

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::build::wire_to_face::WireToFace;
use crate::topopebrep::ds::build_tool::BuildTool;

use super::super::k_part::{
    RESFACE1, RESFACE2, RESNEWSHE, RESNULL, RESSHAPE1, RESSHAPE2, RESUNDEF, SHECLASAUTR,
    SHEGARDAUTR, SHEUNDEF,
};
use super::Builder;

/// `TopoDS_Shape::Complement()` (`TopoDS_Shape.hxx`): the orientation is
/// complemented — `FORWARD` <-> `REVERSED`, `INTERNAL` <-> `EXTERNAL`
/// (`TopAbs::Complement`, `TopAbs.hxx:110-119`).
fn complement(s: &mut TopoShape) {
    s.set_orientation(s.orientation().complement());
}

impl Builder {
    /// `MergeKPartiskole()` (`kpkole.cxx:56-572`).
    pub fn merge_k_part_iskole(&mut self) {
        // kpkole.cxx:65-70
        if self.my_is_k_part != 1 {
            return;
        }

        // kpkole.cxx:72
        let (shape1, shape2) = match (self.my_shape1.clone(), self.my_shape2.clone()) {
            (Some(s1), Some(s2)) => (s1, s2),
            _ => return,
        };
        self.map_shapes(Some(&shape1), Some(&shape2));

        // kpkole.cxx:78-79
        let sta1 = self.my_state1;
        let sta2 = self.my_state2;
        self.change_merged(&shape1, sta1);
        self.change_merged(&shape2, sta2);

        // kpkole.cxx:81-85
        let mut lshsd1: Vec<TopoShape> = Vec::new();
        let mut lfhsd1: Vec<TopoShape> = Vec::new();
        let mut lshsd2: Vec<TopoShape> = Vec::new();
        let mut lfhsd2: Vec<TopoShape> = Vec::new();
        self.kp_iskolesh(&shape1, &mut lshsd1, &mut lfhsd1);
        self.kp_iskolesh(&shape2, &mut lshsd2, &mut lfhsd2);

        // kpkole.cxx:87-89. OCCT reads `lShsd1.First()` unchecked; both lists
        // are non-empty here because `KPiskole` returned 1, which implies
        // `KPiskolesh` answered true on both operands.
        let (Some(sol1), Some(sol2)) = (lshsd1.first().cloned(), lshsd2.first().cloned()) else {
            return;
        };
        self.change_merged(&sol1, sta1);
        self.change_merged(&sol2, sta2);

        // kpkole.cxx:92-97
        if self.kp_map_f1f2().is_empty() {
            return;
        }

        // kpkole.cxx:123-152: `pfOU` is the first key that carries a non-empty
        // list, `plfIN` that list, `pfIN` its first element.
        let mut pfin: Option<TopoShape> = None;
        let mut plfin: Option<Vec<TopoShape>> = None;
        let mut lfin: Vec<TopoShape> = Vec::new();
        let mut pfou: Option<TopoShape> = None;
        for (f, los) in self.kp_map_f1f2() {
            if los.is_empty() {
                continue;
            }
            if plfin.is_none() {
                plfin = Some(los.clone());
            }
            if pfou.is_none() {
                pfou = Some(f.clone());
            }
            if pfin.is_none() {
                pfin = los.first().cloned();
            }
            lfin.extend(los.iter().cloned());
        }

        // kpkole.cxx:154-165
        let (Some(pfou), Some(plfin), Some(pfin)) = (pfou, plfin, pfin) else {
            return;
        };

        // kpkole.cxx:182-190
        let rankpfou = self.shape_rank(&pfou);
        let rankpfin = self.shape_rank(&pfin);
        if rankpfou != 1 && rankpfou != 2 {
            return;
        }
        if rankpfin != 1 && rankpfin != 2 {
            return;
        }

        // kpkole.cxx:196-211
        let solfou = if rankpfou == 1 {
            sol1.clone()
        } else {
            sol2.clone()
        };
        let solfin = if rankpfin == 1 {
            sol1.clone()
        } else {
            sol2.clone()
        };

        // kpkole.cxx:214-215
        let stsolfou = self.kp_clas_ss_av(&solfou, Some(&pfou), &solfin);
        let stsolfin = self.kp_clas_ss_ls(&solfin, &lfin, &solfou);

        // kpkole.cxx:216-238
        let (stsol1, stfac1, fac1, stsol2, stfac2, fac2) = if rankpfou == 1 {
            (
                stsolfou,
                State::Out,
                pfou.clone(),
                stsolfin,
                State::In,
                pfin.clone(),
            )
        } else {
            (
                stsolfin,
                State::In,
                pfin.clone(),
                stsolfou,
                State::Out,
                pfou.clone(),
            )
        };

        // kpkole.cxx:239-248
        let (ires, icla1, icla2) = self.kp_iskole_analyse(stfac1, stfac2, stsol1, stsol2);
        if ires == RESUNDEF {
            return;
        }
        if icla1 == SHEUNDEF || icla2 == SHEUNDEF {
            return;
        }

        // kpkole.cxx:250-266: the shells accessing `fac1` in `sol1` and `fac2`
        // in `sol2`.
        let map1 = map_shapes_and_ancestors(&sol1, ShapeType::Face, ShapeType::Shell);
        let Some(she1) = map1
            .get(&shape_map_key(&fac1))
            .and_then(|l| l.first().cloned())
        else {
            return;
        };
        let map2 = map_shapes_and_ancestors(&sol2, ShapeType::Face, ShapeType::Shell);
        let Some(she2) = map2
            .get(&shape_map_key(&fac2))
            .and_then(|l| l.first().cloned())
        else {
            return;
        };

        // kpkole.cxx:268-269
        self.change_merged(&she1, sta1);
        self.change_merged(&she2, sta2);

        // kpkole.cxx:286-301
        let mut newshe: Option<TopoShape> = None;
        if ires == RESNULL {
            return;
        } else if ires == RESSHAPE1 {
            let mut sh = TopoShape::new(ShapeType::Shell);
            self.my_build_tool.make_shell(&mut sh);
            newshe = Some(she1.clone());
        } else if ires == RESSHAPE2 {
            let mut sh = TopoShape::new(ShapeType::Shell);
            self.my_build_tool.make_shell(&mut sh);
            newshe = Some(she2.clone());
        } else if ires == RESFACE1 {
            // kpkole.cxx:304-317
            if rankpfou == 1 {
                self.change_merged(&shape1, sta1).push(pfou.clone());
                self.change_merged(&fac2, sta2).push(pfou.clone());
            } else if rankpfou == 2 {
                let lmergesha1 = self.change_merged(&shape1, sta1);
                Builder::copy_list(&plfin, lmergesha1);
                let lfac2 = self.change_merged(&fac2, sta2);
                Builder::copy_list(&plfin, lfac2);
            }
            return;
        } else if ires == RESFACE2 {
            // kpkole.cxx:320-335
            if rankpfou == 2 {
                self.change_merged(&shape1, sta1).push(pfou.clone());
                self.change_merged(&fac1, sta1).push(pfou.clone());
            } else if rankpfou == 1 {
                let lmergesha1 = self.change_merged(&shape1, sta1);
                Builder::copy_list(&plfin, lmergesha1);
                let lfac1 = self.change_merged(&fac1, sta1);
                Builder::copy_list(&plfin, lfac1);
            }
            return;
        } else if ires == RESNEWSHE {
            // kpkole.cxx:339-477
            let pairs = self.kp_map_f1f2().to_vec();
            // `NCollection_DataMap<Shape, Shape>` used as an ordered set of
            // the faces already placed in `newshe`.
            let mut addedfaces: Vec<TopoShape> = Vec::new();

            for (f1, lf2) in pairs {
                if lf2.is_empty() {
                    continue;
                }
                let f2 = lf2[0].clone();

                // kpkole.cxx:366-374
                let rankf1 = self.shape_rank(&f1);
                let rankf2 = self.shape_rank(&f2);
                if rankf1 == 0 || rankf2 == 0 {
                    continue;
                }

                // kpkole.cxx:377-390
                let (stf1, stf2) = self.kp_class_ff(&f1, &f2);
                let ires = if rankf1 == 1 {
                    self.kp_iskole_analyse(stf1, stf2, stsol1, stsol2).0
                } else if rankf1 == 2 {
                    self.kp_iskole_analyse(stf2, stf1, stsol2, stsol1).0
                } else {
                    RESUNDEF
                };
                if ires == RESUNDEF {
                    continue;
                }

                // kpkole.cxx:392-403
                let r1 = stsol1 == State::In;
                let r2 = stsol2 == State::In;
                let fac = if rankf1 == 1 {
                    self.kp_makeface(&f1, &lf2, stf1, r1, r2)
                } else {
                    self.kp_makeface(&f1, &lf2, stf1, r2, r1)
                };
                let Some(fac) = fac else {
                    continue;
                };
                // kpkole.cxx:409-411
                if !addedfaces.iter().any(|f| crate::topo_tools::is_same(f, &fac)) {
                    addedfaces.push(fac.clone());
                }

                // kpkole.cxx:414-419
                let statemergef1 = if rankf1 == 1 { sta1 } else { sta2 };
                let statemergef2 = if rankf2 == 2 { sta2 } else { sta1 };
                self.change_merged(&f1, statemergef1).push(fac.clone());
                for f in &lf2 {
                    self.change_merged(f, statemergef2).push(fac.clone());
                }

                // kpkole.cxx:421-439: `she1`'s faces minus the tangents and
                // those already added.
                for facur in explorer_shapes(&she1, ShapeType::Face) {
                    let isfsd = self.kp_map_f1f2_contains(&facur);
                    let isadded = addedfaces
                        .iter()
                        .any(|f| crate::topo_tools::is_same(f, &facur));
                    if isfsd || isadded {
                        continue;
                    }
                    let mut fori = facur.clone();
                    if stsol1 == State::In {
                        complement(&mut fori);
                    }
                    addedfaces.push(fori);
                }

                // kpkole.cxx:441-461: the same for `she2`.
                for facur in explorer_shapes(&she2, ShapeType::Face) {
                    let isfsd = self.kp_map_f1f2_contains(&facur);
                    let isadded = addedfaces
                        .iter()
                        .any(|f| crate::topo_tools::is_same(f, &facur));
                    if isfsd || isadded {
                        continue;
                    }
                    let mut fori = facur.clone();
                    if stsol2 == State::In {
                        complement(&mut fori);
                    }
                    addedfaces.push(fori);
                }
            }

            // kpkole.cxx:463-477
            if !addedfaces.is_empty() {
                let mut sh = TopoShape::new(ShapeType::Shell);
                self.my_build_tool.make_shell(&mut sh);
                // `cxx:470`: `MyBuildTool.Closed(newshe, true)` — the exact
                // closedness of the shell is not checked ("NYI").
                self.my_build_tool.closed(&sh, true);
                for f in &addedfaces {
                    self.my_build_tool.add_shell_face(&mut sh, f);
                }
                newshe = Some(sh);
            }
        }

        // kpkole.cxx:488-492
        let mut newsol: Option<TopoShape> = None;
        if let Some(she) = newshe.as_ref() {
            let mut sol = TopoShape::new(ShapeType::Solid);
            self.my_build_tool.make_solid(&mut sol);
            self.my_build_tool.add_solid_shell(&mut sol, she);
            newsol = Some(sol);
        }

        // kpkole.cxx:494-528
        if icla1 == SHECLASAUTR || icla1 == SHEGARDAUTR {
            let mut loshe1: Vec<TopoShape> = Vec::new();
            for shecur in explorer_shapes(&sol1, ShapeType::Shell) {
                if crate::topo_tools::is_equal(&she1, &shecur) {
                    continue;
                }
                if icla1 == SHECLASAUTR {
                    let state1 = self.kp_clas_ss_av(&shecur, Some(&fac1), &sol2);
                    if state1 == sta1 {
                        loshe1.push(shecur);
                    }
                } else if icla1 == SHEGARDAUTR {
                    loshe1.push(shecur);
                }
            }
            // kpkole.cxx:523-528
            if let Some(sol) = newsol.as_mut() {
                for shecur in &loshe1 {
                    self.my_build_tool.add_solid_shell(sol, shecur);
                }
            }
        }

        // kpkole.cxx:530-564
        if icla2 == SHECLASAUTR || icla2 == SHEGARDAUTR {
            let mut loshe2: Vec<TopoShape> = Vec::new();
            for shecur in explorer_shapes(&sol2, ShapeType::Shell) {
                if crate::topo_tools::is_equal(&she2, &shecur) {
                    continue;
                }
                if icla2 == SHECLASAUTR {
                    let state2 = self.kp_clas_ss_av(&shecur, Some(&fac2), &sol1);
                    if state2 == sta2 {
                        loshe2.push(shecur);
                    }
                } else if icla2 == SHEGARDAUTR {
                    loshe2.push(shecur);
                }
            }
            // kpkole.cxx:558-563
            if let Some(sol) = newsol.as_mut() {
                for shecur in &loshe2 {
                    self.my_build_tool.add_solid_shell(sol, shecur);
                }
            }
        }

        // kpkole.cxx:566-570
        if let Some(sol) = newsol {
            self.change_merged(&shape1, sta1).push(sol);
        }
    }

    /// `KPmakeface(Fac1, LF2, Stfac1, Stfac2, R1, R2)` (`kpkole.cxx:945-1031`):
    /// rebuilds the face holding the tangent one, by copying `Fac1` and adding
    /// `Fac1`'s wires then `LF2`'s.
    ///
    /// OCCT never reads `Stfac2` (`cxx:947`) and throws `Standard_ProgramError`
    /// when `Stfac1 != OUT` or the IN rank resolves to 0 (`cxx:980-987`); the
    /// port drops the unused argument and answers `None`, which is what the
    /// `cxx:403` null test observes.
    pub fn kp_makeface(
        &self,
        fac1: &TopoShape,
        lf2: &[TopoShape],
        stfac1: State,
        r1: bool,
        r2: bool,
    ) -> Option<TopoShape> {
        // kpkole.cxx:954-958
        if stfac1 != State::Out {
            return None;
        }
        // kpkole.cxx:971-979
        let mut fac = fac1.empty_copied();
        let rankf = self.shape_rank(fac1);
        let rankin = if rankf != 0 {
            if rankf == 1 {
                2
            } else {
                1
            }
        } else {
            0
        };
        // kpkole.cxx:985-987
        if rankin == 0 {
            return None;
        }
        let mut lfso: Vec<TopoShape> = Vec::new();
        let mut lfdo: Vec<TopoShape> = Vec::new();
        self.find_same_domain_same_orientation_from(fac1, &mut lfso, &mut lfdo);

        // kpkole.cxx:989-991: `LFIN` is filled and never read.
        let mut lfin: Vec<TopoShape> = Vec::new();
        self.find_same_rank(&lfso, rankin, &mut lfin);
        self.find_same_rank(&lfdo, rankin, &mut lfin);

        // kpkole.cxx:1004-1013
        for wicur in explorer_shapes(fac1, ShapeType::Wire) {
            let mut wori = wicur.clone();
            if r1 {
                complement(&mut wori);
            }
            self.my_build_tool.add_face_wire(&mut fac, &wori);
        }

        // kpkole.cxx:1015-1027
        for fac2 in lf2 {
            for wicur in explorer_shapes(fac2, ShapeType::Wire) {
                let mut wori = wicur.clone();
                if r2 {
                    complement(&mut wori);
                }
                self.my_build_tool.add_face_wire(&mut fac, &wori);
            }
        }

        // kpkole.cxx:1030
        Some(fac)
    }

    /// `FUNKP_KPmakefaces(BU, Fac1, LF2, Stfac1, Stfac2, R1, R2, Lres)`
    /// (`kpkole.cxx:845-941`): the `WireToFace` variant, which rebuilds the
    /// wires' faces through `FaceBuilder` rather than copying `Fac1`.
    ///
    /// UNPORTED: no OCCT caller reaches this `Standard_EXPORT` helper from the
    /// build pipeline (it is a debugging/extension entry point), so the port
    /// keeps it out of the merge path. `WireToFace::MakeFaces` below is what it
    /// needs and is translated in `build/wire_to_face.rs`.
    pub fn fun_kp_makefaces(
        &self,
        fac1: &TopoShape,
        lf2: &[TopoShape],
        r1: bool,
        r2: bool,
        lres: &mut Vec<TopoShape>,
    ) {
        let mut wtof = WireToFace::new();
        for wicur in explorer_shapes(fac1, ShapeType::Wire) {
            let mut wori = wicur.clone();
            if r1 {
                complement(&mut wori);
            }
            wtof.add_wire(&wori);
        }
        for fac2 in lf2 {
            for wicur in explorer_shapes(fac2, ShapeType::Wire) {
                let mut wori = wicur.clone();
                if r2 {
                    complement(&mut wori);
                }
                wtof.add_wire(&wori);
            }
        }
        let mut b = Builder::new(BuildTool::new());
        wtof.make_faces(fac1, &mut b, lres);
    }
}
