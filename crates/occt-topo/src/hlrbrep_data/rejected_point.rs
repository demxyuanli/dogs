//! Interference classification of `HLRBRep_Data` - `RejectedInterference`
//! (`HLRBRep_Data.cxx:1486-1523`), `AboveInterference` (`:1525-1530`) and
//! `RejectedPoint` (`:2330-2595`), plus the local `AdjustParameter`
//! (`cxx:498-516`).
//!
//! `RejectedInterference` maps the flat cursor `iInterf` of the current
//! edge/face pair back to either a single point (`iInterf <= NbPoints`) or one
//! end of an `On` segment (`On` segments occupy **two** cursor slots, hence
//! `n & 1` and `n >> 1`), and asks `RejectedPoint` whether that point can be
//! turned into an `HLRAlgo_Interference` of `myIntf`.
//!
//! `RejectedPoint` is the depth test of the hiding pass:
//!
//! * the 3D parameters are recovered from the 2D ones (`Parameter3d`), and the
//!   depth difference `dz = Z(p1) - Z(p2)` decides everything: `dz >= TolZ`
//!   means the point is **above** the face and is rejected
//!   (`myAboveIntf = true`); `dz <= -TolZ` means `IN` (behind the face),
//!   anything closer is `ON`;
//! * for an auto-intersection (`myLE == myFE`) the two parameters are swapped
//!   when OCCT's own sorting put the wrong end first, so the same depth sign
//!   convention holds;
//! * the transition of `myIntf` is derived from the first transition and the
//!   orientation of `myFE` (complemented when `iFaceBack`), the "on-boundary"
//!   orientation from the position on the curve (`Head`/`Middle`/`End`), and
//!   both ends snapped inward by `CutBig = 10%` of the edge range when the
//!   corresponding `VerAtSta`/`VerAtEnd` flag is set - a vertical end has an
//!   undefined tangent, so `D2` supplies the normal as a stand-in;
//! * `ON` points additionally record the two 3D states through `EdgeState`;
//! * `Level` (the `decal` value) is `2` on an internal face edge, `0` for a
//!   forward crossing of forward-facing edge and face, and `1` otherwise.
//!
//! Deviations:
//! * `myLEData` / `myFEData` / `myLEGeom` / `myFEGeom` are raw pointers in
//!   OCCT; the port re-reads the entity from its index, so no two borrows of
//!   `myEData` are live at once and the scalar caches (`myLETol`) are used
//!   where OCCT would dereference the pointer.
//! * `t1` / `t2` are written by `AdjustParameter` and never read again in
//!   OCCT; the port keeps them only because the signature is part of the
//!   control flow.
//! * `myClassifier` / `mySLProps` failures answer `OUT` where OCCT would
//!   dereference a null handle.

use occt_core::gp::GpVec2d;
use occt_core::intres2d::{
    IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dSituation, IntRes2dTypeTrans,
};

use crate::abs::{Orientation, State};
use crate::hlr_algo::HlrAlgoInterference;
use crate::hlrbrep_edge_data::HlrBRepEdgeData;

use super::HlrBRepData;

/// `static const double CutBig = 1.e-1` (`cxx:57`): the fraction of an edge
/// range cut away at a vertical end.
const CUT_BIG: f64 = 1.0e-1;

/// `static const double DERIVEE_PREMIERE_NULLE = 0.000000000001` (`cxx:60`).
const DERIVEE_PREMIERE_NULLE: f64 = 0.000_000_000_001;

/// `static void AdjustParameter(HLRBRep_EdgeData* E, const bool h, double& p,
/// float& t)` (`cxx:498-516`): replace `p` / `t` by the corresponding bound of
/// `E` and, when that end is vertical, pull the parameter `CutBig` of the
/// range inward.
pub(crate) fn adjust_parameter(edge: &HlrBRepEdgeData, h: bool, p: &mut f64, t: &mut f32) {
    let (p1, t1, p2, t2) = edge.status().bounds();
    if h {
        // `E->Status().Bounds(p, t, p2, t2);`
        *p = p1;
        *t = t1;
        if edge.ver_at_sta() {
            *p += (p2 - *p) * CUT_BIG;
        }
    } else {
        // `E->Status().Bounds(p1, t1, p, t);`
        *p = p2;
        *t = t2;
        if edge.ver_at_end() {
            *p -= (*p - p1) * CUT_BIG;
        }
    }
}

impl HlrBRepData {
    /// `RejectedInterference()` (`cxx:1486-1523`).
    pub fn rejected_interference(&mut self) -> bool {
        if self.i_interf <= self.my_nb_points {
            // `myIntersector.Point(iInterf)`, OCCT 1-based.
            let p = self.my_intersector.point((self.i_interf - 1) as usize);
            return self.rejected_point(&p, Orientation::External, 0);
        }

        // `n` counts within the segment part: two cursor slots per segment.
        let n = self.i_interf - self.my_nb_points;
        let mut first_point = (n & 1) != 0;
        let mut nseg = n >> 1;
        if first_point {
            nseg += 1;
        }

        // `HLRBRep_Curve::Parameter3d` of both `ParamOnFirst` values decides
        // which end comes first in 3D (`cxx:1500-1508`).
        let seg = self.my_intersector.segment((nseg - 1) as usize);
        let le = self.my_le as usize;
        let pf = self.my_e_data[le]
            .geometry()
            .parameter3d(seg.first_point().param_on_first());
        let pl = self.my_e_data[le]
            .geometry()
            .parameter3d(seg.last_point().param_on_first());
        if pf > pl {
            first_point = !first_point;
        }

        if first_point {
            let p = *seg.first_point();
            self.rejected_point(&p, Orientation::Forward, nseg)
        } else {
            let p = *seg.last_point();
            self.rejected_point(&p, Orientation::Reversed, -nseg)
        }
    }

    /// `AboveInterference()` (`cxx:1525-1530`): the last `RejectedPoint`
    /// answered "the point is above the face".
    pub fn above_interference(&self) -> bool {
        self.my_above_intf
    }

    /// `Interference()` (`hxx:141`): the interference `RejectedPoint` fills.
    pub fn interference(&self) -> &HlrAlgoInterference {
        &self.my_intf
    }

    /// The non-const `Interference()` (`hxx:141`).
    pub fn interference_mut(&mut self) -> &mut HlrAlgoInterference {
        &mut self.my_intf
    }

    /// `RejectedPoint(PInter, BoundOri, NumSeg)` (`cxx:2330-2595`): `true`
    /// when the point must not become an interference.
    pub fn rejected_point(
        &mut self,
        p_inter: &IntRes2dIntersectionPoint,
        bound_ori: Orientation,
        num_seg: i32,
    ) -> bool {
        let le = self.my_le as usize;
        let fe = self.my_fe as usize;

        let mut ind = 0;
        let mut t1 = 0.0f32;
        let mut t2 = 0.0f32;
        let mut orie = Orientation::Forward;
        let mut or2 = Orientation::Internal;
        let mut inverted = false;
        // `TolZ = myBigSize * 0.00001` (`cxx:2347`).
        let tol_z = self.my_big_size * 0.000_01;

        let mut p1 = self.my_e_data[le]
            .geometry()
            .parameter3d(p_inter.param_on_first());
        let mut p2 = self.my_e_data[fe]
            .geometry()
            .parameter3d(p_inter.param_on_second());
        let mut dz = self.my_e_data[le].geometry().z(p1) - self.my_e_data[fe].geometry().z(p2);

        if self.my_le == self.my_fe {
            // `cxx:2352-2360`: an auto-intersection can be inverted.
            if dz >= tol_z {
                inverted = true;
                std::mem::swap(&mut p1, &mut p2);
                dz = -dz;
            }
        }

        if dz >= tol_z {
            self.my_above_intf = true;
            return true;
        }
        self.my_above_intf = false;
        let st = if dz <= -tol_z { State::In } else { State::On };

        // `cxx:2371-2380`: the transitions follow the inversion.
        let (tr1, tr2) = if inverted {
            (p_inter.transition_of_second(), p_inter.transition_of_first())
        } else {
            (p_inter.transition_of_first(), p_inter.transition_of_second())
        };

        // `cxx:2382-2405`.
        if self.i_face_test {
            if self.my_le == self.my_fe {
                if st == State::In {
                    self.my_e_data[le].set_simple(false);
                }
            } else if self.my_same_vertex
                && (st == State::On
                    || tr1.position_on_curve() != IntRes2dPosition::Middle
                    || tr2.position_on_curve() != IntRes2dPosition::Middle)
            {
                return true;
            }
            if st == State::In {
                self.i_face_smpl = false;
            }
        }

        // `cxx:2407-2430`: the interference transition from the first one.
        match tr1.transition_type() {
            IntRes2dTypeTrans::In => {
                orie = if self.my_fe_ori == Orientation::Reversed {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                };
            }
            IntRes2dTypeTrans::Out => {
                orie = if self.my_fe_ori == Orientation::Reversed {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
            }
            IntRes2dTypeTrans::Touch => match tr1.situation() {
                IntRes2dSituation::Inside => {
                    orie = if self.my_fe_ori == Orientation::Reversed {
                        Orientation::External
                    } else {
                        Orientation::Internal
                    };
                }
                IntRes2dSituation::Outside => {
                    orie = if self.my_fe_ori == Orientation::Reversed {
                        Orientation::Internal
                    } else {
                        Orientation::External
                    };
                }
                IntRes2dSituation::Unknown => return true,
            },
            IntRes2dTypeTrans::Undecided => return true,
        }

        // `cxx:2432-2435`: a back face complements the transition.
        if self.i_face_back {
            orie = orie.complement();
        }

        // `cxx:2436-2448`: the "on-boundary" orientation from the position.
        let ori = match tr1.position_on_curve() {
            IntRes2dPosition::Head => Orientation::Forward,
            IntRes2dPosition::Middle => Orientation::Internal,
            IntRes2dPosition::End => Orientation::Reversed,
        };

        if st != State::Out {
            // `cxx:2451-2510`: correct the transition of `myFE` at a
            // non-middle position.
            if tr2.position_on_curve() != IntRes2dPosition::Middle {
                if self.my_same_vertex {
                    // "si intersection a une extremite verticale !"
                    return true;
                }

                let mut douteux = false;
                let psav = p2;
                let mut tgsav = GpVec2d::zero();
                if tr2.position_on_curve() == IntRes2dPosition::Head {
                    ind = self.my_e_data[fe].v_sta();
                    or2 = Orientation::Forward;
                    adjust_parameter(&self.my_e_data[fe], true, &mut p2, &mut t2);
                    if self.my_e_data[fe].ver_at_sta() {
                        douteux = true;
                        let (_, tg, nm) = self.my_e_data[fe].geometry().d2(psav);
                        tgsav = if tg.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
                            nm
                        } else {
                            tg
                        };
                    }
                } else {
                    ind = self.my_e_data[fe].v_end();
                    or2 = Orientation::Reversed;
                    adjust_parameter(&self.my_e_data[fe], false, &mut p2, &mut t2);
                    if self.my_e_data[fe].ver_at_end() {
                        douteux = true;
                        let (_, tg, nm) = self.my_e_data[fe].geometry().d2(psav);
                        tgsav = if tg.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
                            nm
                        } else {
                            tg
                        };
                    }
                }

                let (_, tgfe) = self.my_e_data[fe].geometry().d1(p2);
                if douteux && tgfe.dot(&tgsav) < 0.0 {
                    if orie == Orientation::Forward {
                        orie = Orientation::Reversed;
                    } else if orie == Orientation::Reversed {
                        orie = Orientation::Forward;
                    }
                }
                self.my_intf.change_boundary().set_2d(self.my_fe, p2);
            }

            // `cxx:2511-2559`: the same correction on `myLE`, skipped when the
            // point is in the middle of the edge.
            if ori != Orientation::Internal {
                let mut douteux = false;
                let psav = p1;
                let mut tgsav = GpVec2d::zero();
                if ori == Orientation::Forward {
                    adjust_parameter(&self.my_e_data[le], true, &mut p1, &mut t1);
                    if self.my_e_data[le].ver_at_sta() {
                        douteux = true;
                        let (_, tg, nm) = self.my_e_data[le].geometry().d2(psav);
                        tgsav = if tg.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
                            nm
                        } else {
                            tg
                        };
                    }
                } else {
                    adjust_parameter(&self.my_e_data[le], false, &mut p1, &mut t1);
                    if self.my_e_data[le].ver_at_end() {
                        douteux = true;
                        let (_, tg, nm) = self.my_e_data[le].geometry().d2(psav);
                        tgsav = if tg.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
                            nm
                        } else {
                            tg
                        };
                    }
                }

                if douteux {
                    let (_, tgle) = self.my_e_data[le].geometry().d1(p1);
                    if tgle.dot(&tgsav) < 0.0 {
                        if orie == Orientation::Forward {
                            orie = Orientation::Reversed;
                        } else if orie == Orientation::Reversed {
                            orie = Orientation::Forward;
                        }
                    }
                }
            }

            // `cxx:2561-2566`: an `ON` point also records the two 3D states.
            if st == State::On {
                let (stbef, staft) = self.edge_state(p1, p2);
                self.my_intf.change_boundary().set_state_3d(stbef, staft);
            }
        }

        // `cxx:2569-2580`.
        let decal = if self.my_fe_internal {
            2
        } else if st == State::In && ori == Orientation::Forward && orie == Orientation::Forward {
            0
        } else {
            1
        };

        // `cxx:2581-2593`.
        let inter = self.my_intf.change_intersection();
        inter.set_orientation(ori);
        inter.set_level(decal);
        inter.set_seg_index(num_seg);
        inter.set_index(ind);
        inter.set_parameter(p1);
        inter.set_tolerance(self.my_le_tol);
        inter.set_state(st);
        self.my_intf.set_orientation(or2);
        self.my_intf.set_transition(orie);
        self.my_intf.set_boundary_transition(bound_ori);
        self.my_intf.change_boundary().set_2d(self.my_fe, p2);
        false
    }
}
