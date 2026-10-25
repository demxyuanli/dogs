//! `ChFi3d_Builder::ComputeData` — the face / face walking driver.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:1676-2507`
//! (OCCT 8.0.0). This is the overload that walks with
//! `BRepBlend_Walking` (face / face) instead of `BRepBlend_SurfRstLineBuilder`
//! / `BRepBlend_RstRstLineBuilder`: the line it produces is the band between
//! two faces, and the trailing half extends the two ends with the `Gd1` / `Gd2`
//! / `Gf1` / `Gf2` flags and truncates the line when only one side ran into a
//! face (`cxx:2336-2505`).
//!
//! The control flow is kept statement by statement:
//!
//! | OCCT | here |
//! |---|---|
//! | offset guide lookup, `ConstThroatWithPenetrationChamfer` (`:1705-1720`) | `offset_index` |
//! | extension values reset (`:1722-1725`) | [`SurfData::set_first_extension_value`] |
//! | `F1` / `F2` through `occ::down_cast<BRepAdaptor_Surface>` (`:1727-1738`) | the `f1` / `f2` parameters |
//! | `BRepBlend_Walking` + `Check2d(false)` (`:1744-1749`) | [`Walking`] |
//! | `Target` from `intf` / `intl` (`:1756-1775`) | inlined |
//! | the two singular extremities (`:1777-1817`) | [`comp_blend_point`] |
//! | `PerformFirstSection` (`:1819-1846`) | [`Walking::perform_first_section_recadre`] |
//! | the `again < 3` retry loop (`:1848-2004`) | kept verbatim |
//! | twist flags (`:2006-2019`) | [`Walking::twist_on_s1`] |
//! | the beginning extension (`:2038-2420`) | [`ChFi3dBuilder::search_face`], [`is_obst`] |
//! | the end extension (`:2422-2505`) | same, mirrored |
//! | the two `Gd` / `Gf` ladders and the truncations (`:2336-2505`) | [`ChFi3dBuilder::apply_extensions`] in `extension.rs` |
//!
//! Divergences that are not heuristics:
//! - `F1` / `F2` come from `occ::down_cast<BRepAdaptor_Surface>(S1)->Face()`.
//!   The port's `Adaptor3d_Surface` is the bare `Arc<dyn Surface>` and carries
//!   no face, so the caller passes the two faces next to the adaptors;
//!   `Option<&Face>` is the null `TopoDS_Face` OCCT keeps when the adaptor is
//!   not a `BRepAdaptor_Surface`, and every consumer of the face answers the
//!   null-face behaviour (`is_input` sees no wire, `SearchFace` answers
//!   `false`),
//! - `occ::handle<ChFiDS_ElSpine>` is `Arc<ElSpine>` in the port, i.e. shared
//!   and immutable once handed out. The two mutations OCCT performs on the
//!   guideline (`FirstParameter` / `LastParameter` / `SetOrigin`,
//!   `:1889-1897`) go through `Arc::make_mut` on the caller's slot, so the
//!   caller must pass the slot the `ChFiDS_Spine` owns
//!   (`Spine::el_spines_mut()[k]`) for the new range to be what later
//!   `Spine::el_spine(w)` lookups return. The `BRepBlend_Walking` built a few
//!   lines above holds its own `Arc`, so it keeps the pre-`SetOrigin` curve;
//!   `Complete` / `Continu` only read the *saved* parameters and `D1(param)`,
//!   which `SetOrigin` does not move, so the walk is unaffected there,
//! - `int& intf` / `int& intl` are `&mut i32`; `intf` is written with `0` / `1`
//!   and only ever read as `>= 2`, `>= 3`, `>= 5` afterwards (`:2028-2031`,
//!   `:2402-2419`), so the truncation to `0` / `1` is the C++ one,
//! - `Lin->TransitionOnS1()` raises `StdFail_NotDone` when the line carries no
//!   transition; the port answers `false` for the whole call there (the walk
//!   always stores both transitions, `BRepBlend_Walking.cxx:2258-2270`),
//! - `Spine->FirstVertex()` / `LastVertex()` raise on a guideline that was
//!   never loaded; the port leaves the obstacle flag `false` there,
//! - `SearchIndex` / `UpdateLine` / `CompBlendPoint` / `IsObst` are the
//!   file statics of the same `.cxx`, ported at
//!   [`crate::chfi3d::builder_6_helpers`].
//! - the extension half lives in `extension.rs` (same `.cxx` range, see the
//!   table above); the split is for file length only.

mod extension;

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::Surface;
use occt_math::MathVector;

use crate::blend::{BlendFuncInv, BlendFunction};
use crate::brep_blend::line::Line;
use crate::brep_blend::walking::Walking;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{ChamfMode, ElSpine, SpineHandle, SurfData};
use crate::chfi3d::builder_6_helpers::{comp_blend_point, is_obst};
use crate::chfi3d::builder_helpers::fil_common_point;
use crate::shape::{Edge, Face, Vertex};

use extension::ExtensionState;

use super::compute_data::point_parameter;
use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `ComputeData(Data, HGuide, Spine, Lin, S1, I1, S2, I2, Func, FInv,
    /// PFirst, MaxStep, Fleche, tolguide, First, Last, Inside, Appro, Forward,
    /// Soldep, intf, intl, Gd1, Gd2, Gf1, Gf2, RecOnS1, RecOnS2)`
    /// (`ChFi3d_Builder_6.cxx:1676-2507`).
    ///
    /// Representation of elementary fillet by path between two faces.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn compute_data_face_face(
        &mut self,
        data: &mut SurfData,
        hguide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        lin: &mut Option<Rc<RefCell<Line>>>,
        s1: Arc<dyn Surface>,
        i1: Rc<RefCell<TopolTool>>,
        f1: Option<&Face>,
        s2: Arc<dyn Surface>,
        i2: Rc<RefCell<TopolTool>>,
        f2: Option<&Face>,
        func: &mut dyn BlendFunction,
        finv: &mut dyn BlendFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
        gd1: &mut bool,
        gd2: &mut bool,
        gf1: &mut bool,
        gf2: &mut bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
    ) -> bool {
        // Get offset guide if exists (`cxx:1705-1720`). OCCT walks the two lists
        // with two iterators and keeps the offset guide of the `HGuide` entry;
        // the port keeps its index.
        let mut offset_index: Option<usize> = None;
        if let Some(spine) = spine {
            if let Ok(guard) = spine.read() {
                if guard.mode() == ChamfMode::ConstThroatWithPenetrationChamfer {
                    for (k, hel) in guard.el_spines().iter().enumerate() {
                        if Arc::ptr_eq(hel, &*hguide) {
                            offset_index = Some(k);
                        }
                    }
                }
            }
        }

        // The extensions are created in case of output of two domains directly
        // and not by path (`cxx:1722-1725`).
        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        // Path framing variables (`cxx:1740-1742`).
        let mut tol_guide = tol_guide;
        let nbptmin = 4;

        // `BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide)` (`cxx:1744`).
        let mut the_walk = Walking::new(
            Arc::clone(&s1),
            Arc::clone(&s2),
            i1,
            i2,
            Arc::clone(hguide),
        );
        // `TheWalk.Check2d(false)` (`cxx:1749`).
        the_walk.check2d(false);

        let mut ms = max_step;
        let mut sp_first = hguide.first_parameter();
        let mut sp_last = hguide.last_parameter();
        let reverse = !forward || inside;
        let mut target;
        if reverse {
            // `Target = SpFirst; if (!intf) Target = Last;` (`cxx:1760-1767`).
            target = sp_first;
            if *intf == 0 {
                target = *last;
            }
        } else {
            // `Target = SpLast + std::abs(SpLast); if (!intl) Target = Last;`
            // (`cxx:1768-1775`).
            target = sp_last + sp_last.abs();
            if *intl == 0 {
                target = *last;
            }
        }

        // In case the singularity is pre-determined, the path is indicated
        // (`cxx:1777-1817`).
        if let Some(spine) = spine {
            let (tgt_first, tgt_last) = match spine.read() {
                Ok(guard) => (guard.is_tangency_extremity(true), guard.is_tangency_extremity(false)),
                Err(_) => (false, false),
            };
            for is_first in [true, false] {
                let tangency = if is_first { tgt_first } else { tgt_last };
                if !tangency {
                    continue;
                }
                let (vertex, edge, param) = {
                    let Ok(guard) = spine.read() else {
                        continue;
                    };
                    if is_first {
                        (
                            guard.first_vertex(),
                            guard.edges(1).cloned(),
                            Some(guard.first_parameter()),
                        )
                    } else {
                        (
                            guard.last_vertex(),
                            guard.edges(guard.nb_edges()).cloned(),
                            guard.last_parameter(),
                        )
                    }
                };
                let (Some(vertex), Some(edge), Some(param)) = (vertex, edge, param) else {
                    continue;
                };
                let (Some(vertex), Some(edge)) = (Vertex::wrap(vertex), Edge::wrap(edge)) else {
                    continue;
                };
                // `if (CompBlendPoint(V, E, param, F1, F2, BP))`
                // (`cxx:1787`): a null face carries no p-curve, which is the
                // `false` of `CompBlendPoint`.
                let (Some(f1), Some(f2)) = (f1, f2) else {
                    continue;
                };
                let Some(bp) = comp_blend_point(&vertex, &edge, param, f1, f2) else {
                    continue;
                };
                // `math_Vector vec(1, 4); BP.ParametersOnS1(vec(1), vec(2));
                // BP.ParametersOnS2(vec(3), vec(4));` (`cxx:1789-1791`).
                let (Ok((u1, v1)), Ok((u2, v2))) =
                    (bp.parameters_on_s1(), bp.parameters_on_s2())
                else {
                    continue;
                };
                let mut vec = MathVector::with_init(1, 4, 0.0);
                vec.set_value(1, u1);
                vec.set_value(2, v1);
                vec.set_value(3, u2);
                vec.set_value(4, v2);
                // `Func.Set(param); if (Func.IsSolution(vec, tolapp3d))`
                // (`cxx:1792-1796`).
                func.set(param);
                if func.is_solution(&vec, self.tolapp3d) {
                    the_walk.add_singular_point(bp);
                }
            }
        }

        // The starting solution is reframed if necessary (`cxx:1819-1846`).
        let mut par_sol = MathVector::new(1, 4);
        let mut new_first = p_first;
        if rec_on_s1 || rec_on_s2 {
            if !the_walk.perform_first_section_recadre(
                func,
                finv,
                p_first,
                target,
                soldep,
                self.tolapp3d,
                tol_guide,
                rec_on_s1,
                rec_on_s2,
                &mut new_first,
                &mut par_sol,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        // First the valid part is calculated, without caring for the extensions
        // (`cxx:1848-2004`).
        let mut again = 0;
        let mut tchernobyl = false;
        let mut u1sov = 0.0;
        let mut u2sov = 0.0;
        // `TopoDS_Face bif;` (`cxx:1853`): the last face `SearchFace` answered.
        let mut bif: Option<Face> = None;
        // Max step is relevant, but too great; the vector is required to detect
        // the twists (`cxx:1854-1860`).
        if (*last - *first).abs() <= ms * 5.0
            && (*last - *first).abs() >= 0.01 * (new_first - target).abs()
        {
            ms = (*last - *first).abs() * 0.2;
        }

        let mut line = Line::new();
        let mut nbpnt;
        while again < 3 {
            // Path (`cxx:1864-1875`).
            if again == 0 && ms < 5.0 * tol_guide {
                ms = 5.0 * tol_guide;
            } else if 5.0 * tol_guide > ms {
                tol_guide = ms / 5.0;
            }
            the_walk.perform(
                func,
                finv,
                new_first,
                target,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol,
                fleche,
                approx,
            );
            if !the_walk.is_done() {
                return false;
            }
            // `Lin = TheWalk.Line()` (`cxx:1884`).
            let Ok(walk_line) = the_walk.line() else {
                return false;
            };
            line = walk_line.clone();

            // `if (HGuide->IsPeriodic() && Inside)` (`cxx:1885-1898`).
            if hguide.is_periodic() && inside {
                sp_first = point_parameter(&line, 1);
                sp_last = sp_first + hguide.period().unwrap_or(0.0);
                {
                    let hg = Arc::make_mut(hguide);
                    hg.set_first_parameter(sp_first);
                    hg.set_last_parameter(sp_last);
                    let _ = hg.set_origin(sp_first);
                }
                // `if (!OffsetHGuide.IsNull())` (`cxx:1892-1897`).
                if let (Some(spine), Some(k)) = (spine, offset_index) {
                    if let Ok(mut guard) = spine.write() {
                        if let Some(slot) = guard.offset_el_spines_mut().get_mut(k) {
                            let off = Arc::make_mut(slot);
                            off.set_first_parameter(sp_first);
                            off.set_last_parameter(sp_last);
                            let _ = off.set_origin(sp_first);
                        }
                    }
                }
            }

            // `bool complmnt = true; if (Inside) complmnt =
            // TheWalk.Complete(Func, FInv, SpLast); if (!complmnt) return false;`
            // (`cxx:1899-1910`).
            let mut complmnt = true;
            if inside {
                complmnt = the_walk.complete(func, finv, sp_last).unwrap_or(false);
            }
            if !complmnt {
                return false;
            }

            // The result is controlled using two criteria (`cxx:1912-2003`):
            // enough points, and one went far enough. `Complete` may have
            // appended points, so the line is read again.
            let Ok(walk_line) = the_walk.line() else {
                return false;
            };
            line = walk_line.clone();
            nbpnt = line.nb_points();
            if nbpnt == 0 {
                return false;
            }
            let fpointpar = point_parameter(&line, 1);
            let lpointpar = point_parameter(&line, nbpnt);

            let factor = 1.0 / (nbptmin as f64 + 1.0);
            let mut okdeb = forward && !inside;
            let mut okfin = !forward && !inside;
            if !okdeb {
                let narc1 = line.start_point_on_first().nb_point_on_rst();
                let narc2 = line.start_point_on_second().nb_point_on_rst();
                okdeb = narc1 > 0 || narc2 > 0 || (fpointpar - *first) < 10.0 * tol_guide;
            }
            if !okfin {
                let narc1 = line.end_point_on_first().nb_point_on_rst();
                let narc2 = line.end_point_on_second().nb_point_on_rst();
                okfin = narc1 > 0 || narc2 > 0 || (*last - lpointpar) < 10.0 * tol_guide;
            }
            if !okdeb || !okfin || nbpnt == 1 {
                // It drags: the controls are extended, and a satisfactory
                // maximum step is expected. When that was already done, quit.
                if tchernobyl {
                    return false;
                }
                tchernobyl = true;
                the_walk.check(false);
                if nbpnt == 1 {
                    ms *= 0.01;
                } else {
                    ms = (lpointpar - fpointpar) / nbpnt as f64;
                }
            } else if nbpnt < nbptmin as usize {
                if again == 0 {
                    u1sov = fpointpar;
                    u2sov = lpointpar;
                    ms = (lpointpar - fpointpar) * factor;
                } else if again == 1 {
                    if (fpointpar - u1sov).abs() >= tol_guide
                        || (lpointpar - u2sov).abs() >= tol_guide
                    {
                        ms = (lpointpar - fpointpar) * factor;
                    } else {
                        return false;
                    }
                }
                again += 1;
            } else {
                again = 3;
            }
        }

        // `if (TheWalk.TwistOnS1()) Data->TwistOnS1(true);` and S2
        // (`cxx:2006-2019`).
        if the_walk.twist_on_s1() {
            data.set_twist_on_s1(true);
        }
        if the_walk.twist_on_s2() {
            data.set_twist_on_s2(true);
        }

        // Here there is a more or less presentable result (`cxx:2021-2024`).
        *gd1 = false;
        *gd2 = false;
        *gf1 = false;
        *gf2 = false;

        let unseulsuffitdeb = *intf >= 2;
        let unseulsuffitfin = *intl >= 2;
        let noproldeb = *intf >= 3;
        let noprolfin = *intl >= 3;

        let rab = 0.03 * (sp_last - sp_first);

        let mut debarc1 = false;
        let mut debarc2 = false;
        let mut debcas1 = false;
        let mut debcas2 = false;
        let mut debobst1 = false;
        let mut debobst2 = false;

        let mut finarc1 = false;
        let mut finarc2 = false;
        let mut fincas1 = false;
        let mut fincas2 = false;
        let mut finobst1 = false;
        let mut finobst2 = false;

        let mut backw_continue_failed = false;
        if reverse && *intf != 0 {
            let mut narc1 = line.start_point_on_first().nb_point_on_rst() as i32;
            let mut narc2 = line.start_point_on_second().nb_point_on_rst() as i32;
            if narc1 != 0 {
                let Ok(tran) = line.transition_on_s1() else {
                    return false;
                };
                let sp = line.start_point_on_first().clone();
                fil_common_point(&sp, tran, true, data.vertex_mut(true, 1), self.tolapp3d);
                debarc1 = true;
                if !self.search_face(spine, data.vertex_first_on_s1(), f1, &mut bif) {
                    // It is checked if there is not an obstacle.
                    debcas1 = true;
                    if let Some(sp) = spine {
                        if let Ok(guard) = sp.read() {
                            if guard.is_periodic() {
                                debobst1 = true;
                            } else if let Some(v) = guard.first_vertex() {
                                debobst1 = is_obst(data.vertex_first_on_s1(), &v, &self.my_ve_map);
                            }
                        }
                    }
                }
            }
            if narc2 != 0 {
                let Ok(tran) = line.transition_on_s2() else {
                    return false;
                };
                let sp = line.start_point_on_second().clone();
                fil_common_point(&sp, tran, true, data.vertex_mut(true, 2), self.tolapp3d);
                debarc2 = true;
                if !self.search_face(spine, data.vertex_first_on_s2(), f2, &mut bif) {
                    debcas2 = true;
                    if let Some(sp) = spine {
                        if let Ok(guard) = sp.read() {
                            if guard.is_periodic() {
                                debobst2 = true;
                            } else if let Some(v) = guard.first_vertex() {
                                debobst2 = is_obst(data.vertex_first_on_s2(), &v, &self.my_ve_map);
                            }
                        }
                    }
                }
            }
            let mut oncontinue = !noproldeb && (narc1 != 0 || narc2 != 0);
            if debobst1 || debobst2 {
                oncontinue = false;
            } else if debcas1 && debcas2 {
                oncontinue = false;
            } else if (!debcas1 && debarc1) || (!debcas2 && debarc2) {
                oncontinue = false;
            }

            if oncontinue {
                the_walk.classification_on_s1(!debarc1);
                the_walk.classification_on_s2(!debarc2);
                // It should be strict (PMN) (`cxx:2118`).
                the_walk.check2d(true);
                let _ = the_walk.continu(func, finv, target);
                the_walk.classification_on_s1(true);
                the_walk.classification_on_s2(true);
                the_walk.check2d(false);
                if let Ok(walk_line) = the_walk.line() {
                    line = walk_line.clone();
                }
                narc1 = line.start_point_on_first().nb_point_on_rst() as i32;
                narc2 = line.start_point_on_second().nb_point_on_rst() as i32;
                // modified by eap Fri Feb 8 11:43:48 2002 (`cxx:2126-2185`).
                if !debarc1 {
                    if narc1 == 0 {
                        backw_continue_failed =
                            line.start_point_on_first().parameter_on_guide() > target;
                    } else {
                        let Ok(tran) = line.transition_on_s1() else {
                            return false;
                        };
                        let sp = line.start_point_on_first().clone();
                        fil_common_point(&sp, tran, true, data.vertex_mut(true, 1), self.tolapp3d);
                        debarc1 = true;
                        if !self.search_face(spine, data.vertex_first_on_s1(), f1, &mut bif) {
                            // The obstacle test of `cxx:2144-2152` is commented
                            // out in OCCT, so only `debcas1` is raised.
                            debcas1 = true;
                        }
                    }
                }
                if !debarc2 {
                    if narc2 == 0 {
                        backw_continue_failed =
                            line.start_point_on_second().parameter_on_guide() > target;
                    } else {
                        let Ok(tran) = line.transition_on_s2() else {
                            return false;
                        };
                        let sp = line.start_point_on_second().clone();
                        fil_common_point(&sp, tran, true, data.vertex_mut(true, 2), self.tolapp3d);
                        debarc2 = true;
                        if !self.search_face(spine, data.vertex_first_on_s2(), f2, &mut bif) {
                            debcas2 = true;
                        }
                    }
                }
                if backw_continue_failed {
                    // If `backwContinueFailed` is left as is, the walk stops in
                    // this direction; it is continued when there are no more
                    // faces on the side with the arc (`cxx:2186-2196`).
                    let a_cp = if debarc1 {
                        data.vertex_first_on_s1()
                    } else {
                        data.vertex_first_on_s2()
                    };
                    if a_cp.is_on_arc() && bif.is_none() {
                        backw_continue_failed = false;
                    }
                }
            }
        }

        let mut forw_continue_failed = false;
        // modified by eap Fri Feb 8 11:44:11 2002 (`cxx:2199-2331`).
        if forward && *intl != 0 {
            target = sp_last;
            let mut narc1 = line.end_point_on_first().nb_point_on_rst() as i32;
            let mut narc2 = line.end_point_on_second().nb_point_on_rst() as i32;
            if narc1 != 0 {
                let Ok(tran) = line.transition_on_s1() else {
                    return false;
                };
                let sp = line.end_point_on_first().clone();
                fil_common_point(&sp, tran, false, data.vertex_mut(false, 1), self.tolapp3d);
                finarc1 = true;
                if !self.search_face(spine, data.vertex_last_on_s1(), f1, &mut bif) {
                    fincas1 = true;
                    if let Some(sp) = spine {
                        if let Ok(guard) = sp.read() {
                            if let Some(v) = guard.last_vertex() {
                                finobst1 = is_obst(data.vertex_last_on_s1(), &v, &self.my_ve_map);
                            }
                        }
                    }
                }
            }
            if narc2 != 0 {
                let Ok(tran) = line.transition_on_s2() else {
                    return false;
                };
                let sp = line.end_point_on_second().clone();
                fil_common_point(&sp, tran, false, data.vertex_mut(false, 2), self.tolapp3d);
                finarc2 = true;
                if !self.search_face(spine, data.vertex_last_on_s2(), f2, &mut bif) {
                    fincas2 = true;
                    if let Some(sp) = spine {
                        if let Ok(guard) = sp.read() {
                            if let Some(v) = guard.last_vertex() {
                                finobst2 = is_obst(data.vertex_last_on_s2(), &v, &self.my_ve_map);
                            }
                        }
                    }
                }
            }
            let mut oncontinue = !noprolfin && (narc1 != 0 || narc2 != 0);
            if finobst1 || finobst2 {
                oncontinue = false;
            } else if fincas1 && fincas2 {
                oncontinue = false;
            } else if (!fincas1 && finarc1) || (!fincas2 && finarc2) {
                oncontinue = false;
            }

            if oncontinue {
                the_walk.classification_on_s1(!finarc1);
                the_walk.classification_on_s2(!finarc2);
                the_walk.check2d(true);
                let _ = the_walk.continu(func, finv, target);
                the_walk.classification_on_s1(true);
                the_walk.classification_on_s2(true);
                the_walk.check2d(false);
                if let Ok(walk_line) = the_walk.line() {
                    line = walk_line.clone();
                }
                narc1 = line.end_point_on_first().nb_point_on_rst() as i32;
                narc2 = line.end_point_on_second().nb_point_on_rst() as i32;
                // modified by eap Fri Feb 8 11:44:57 2002 (`cxx:2267-2329`).
                if !finarc1 {
                    if narc1 == 0 {
                        forw_continue_failed =
                            line.end_point_on_first().parameter_on_guide() < target;
                    } else {
                        let Ok(tran) = line.transition_on_s1() else {
                            return false;
                        };
                        let sp = line.end_point_on_first().clone();
                        fil_common_point(&sp, tran, false, data.vertex_mut(false, 1), self.tolapp3d);
                        finarc1 = true;
                        if !self.search_face(spine, data.vertex_last_on_s1(), f1, &mut bif) {
                            fincas1 = true;
                        }
                    }
                }
                if !finarc2 {
                    if narc2 == 0 {
                        forw_continue_failed =
                            line.end_point_on_second().parameter_on_guide() < target;
                    } else {
                        let Ok(tran) = line.transition_on_s2() else {
                            return false;
                        };
                        let sp = line.end_point_on_second().clone();
                        fil_common_point(&sp, tran, false, data.vertex_mut(false, 2), self.tolapp3d);
                        finarc2 = true;
                        if !self.search_face(spine, data.vertex_last_on_s2(), f2, &mut bif) {
                            fincas2 = true;
                        }
                    }
                }
                if forw_continue_failed {
                    let a_cp = if finarc1 {
                        data.vertex_last_on_s1()
                    } else {
                        data.vertex_last_on_s2()
                    };
                    if a_cp.is_on_arc() && bif.is_none() {
                        forw_continue_failed = false;
                    }
                }
            }
        }

        // `Nbpnt = Lin->NbPoints(); First = Lin->Point(1).Parameter(); Last =
        // Lin->Point(Nbpnt).Parameter();` (`cxx:2332-2334`).
        nbpnt = line.nb_points();
        *first = point_parameter(&line, 1);
        *last = point_parameter(&line, nbpnt);

        // ============= POSSIBLE INVALIDATION =============
        // Preparation of the extensions by tangent plane (`cxx:2336-2505`).
        self.apply_extensions(&mut ExtensionState {
            line: &mut line,
            data,
            first,
            last,
            sp_first,
            sp_last,
            rab,
            reverse,
            forward,
            intf,
            intl,
            gd1,
            gd2,
            gf1,
            gf2,
            debarc1,
            debarc2,
            debcas1,
            debcas2,
            finarc1,
            finarc2,
            fincas1,
            fincas2,
            backw_continue_failed,
            forw_continue_failed,
            unseulsuffitdeb,
            unseulsuffitfin,
        });

        // `Lin` is the line the driver owns, which is dropped at the end of the
        // call: the port hands the caller the state it ended with.
        *lin = Some(Rc::new(RefCell::new(line)));
        true
    }
}
