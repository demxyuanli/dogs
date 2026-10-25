//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap6` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:2442-2808`
//! (OCCT 8.0.0).
//!
//! `yap6` keeps `EspON` when `FOR` is same-domain with `FCX` and `EG` is not
//! same-domain with `FOR`, deciding on `FUN_ds_shareG` and — for a fuse — on the
//! matter angles of `(FS, FCX)` and `(FOR, FF)` around `EG` and `Esd`
//! (`cxx:2494-2588`). `SO` re-evaluates the transition on the 2-D matter of
//! `FCX` (`cxx:2598-2628`); `DO` never keeps the split ON and only keeps the
//! IN/IN or the OUT part of a cut (`cxx:2633-2708`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:2444-2458`, `:2472-2474`, `:2677-2679`,
//! `:2757-2764`, `:2799-2806`) are not translated, neither are the commented
//! `if (shareG) ...` block (`cxx:2572-2588`) nor the `// nyiFUNRAISE` /
//! `// NYIRaise` markers, which carry no control flow.

use occt_core::gp::GpPnt2d;
use occt_core::precision::ANGULAR;

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::complete_ds::fun_ds_share_g;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::tool::matter::matter_faces;
use crate::topopebrep::tool::project_curve::{fun_tool_par_e, fun_tool_par_e_tol};
use crate::topopebrep::tool::regus::nearest_f::fun_tool_getxx;
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_par_on_e_f};

use super::super::Builder;
use super::helpers::{M_FORWARD, M_REVERSED};
use super::wes2::Wes2Ctx;
use super::wes2_fcx::{orient_ein_f_forward_or_external, FcxCtx};

impl Wes2Ctx {
    /// `yap6` (`cxx:2442-2808`).
    pub(super) fn branch_yap6(&self, fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = pb.data_structure() else {
            return;
        };
        let bds = hds.ds();

        // `cxx:2477-2482`: `FUN_ds_shareG(BDS, iFOR, iFCX, GI, Edge(EspON), shareG)`.
        let Some(share_g) = fun_ds_share_g(&hds, self.i_for, fcx.i_fcx, self.gi, &self.esp_on)
        else {
            return;
        };

        // `cxx:2467-2468`
        let so = fcx.ffin_sd_so;
        let do_ = fcx.ffin_sd_do;

        // `cxx:2463-2465`
        let sta_fcx = fcx.sta_fs;
        let mut b3d = false;
        let mut b2d = false;
        let mut b = false;

        // `cxx:2484-2632`: `FOR` and `FCX` share their geometric domain.
        if so {
            // `cxx:2488`
            let mut force_keep = self.ope_fus && share_g && self.rank_eg == 1;

            if self.ope_fus {
                // `cxx:2492-2503`: the faces connected to `FOR` across `Esd`.
                let Some(esd) = bds.shape(self.iesd, true) else {
                    return;
                };
                let mut lfor: Vec<TopoShape> = Vec::new();
                {
                    let mut connex = Connexity::new();
                    connex.prepare(&hds);
                    connex.face_edge_connex_faces(&self.for_, &esd, &hds, &mut lfor);
                }
                if lfor.is_empty() {
                    // `cxx:2498-2501`: `Esd` is the closing edge of `FOR`.
                    b = share_g;
                } else if lfor.len() > 1 {
                    // `cxx:2503-2506`
                    return;
                } else {
                    // `cxx:2508-2568`
                    let ff = lfor[0].clone();
                    let tola1 = ANGULAR * 1.0e2; // nyitolxpu

                    let Some(par_eg) =
                        fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
                    else {
                        return;
                    };
                    let Some(mat_fs) = matter_faces(&self.fs, &fcx.fcx, &self.eg, par_eg, tola1)
                    else {
                        return;
                    };

                    let tol_esd = BRepTool::edge_tolerance(&Edge(esd.clone())) * self.factor;
                    let Some(par_esd) =
                        fun_tool_par_e_tol(&self.esp_on, self.par_on, &esd, tol_esd)
                    else {
                        return;
                    };
                    let Some(mat_for) = matter_faces(&self.for_, &ff, &esd, par_esd, tola1)
                    else {
                        return;
                    };

                    let sum = mat_fs + mat_for;
                    let sum_is_pi = (sum - std::f64::consts::PI).abs() < tola1;
                    let fs_inf_pi = mat_fs < std::f64::consts::PI;
                    let for_inf_pi = mat_for < std::f64::consts::PI;
                    if sum_is_pi {
                        b = false;
                    } else if sum < std::f64::consts::PI {
                        b = true;
                    } else if fs_inf_pi && for_inf_pi {
                        b = false;
                    } else if (mat_fs - mat_for).abs() < tola1 {
                        b = false;
                    } else {
                        b = mat_fs < mat_for;
                    }
                    if b {
                        force_keep = true;
                    }
                }
            } else if self.ope_cut {
                // `cxx:2590-2593`
                b = !share_g;
            } else if self.ope_com {
                // `cxx:2594-2597`
                b = share_g && self.rank_eg == 1;
            }

            if b {
                // `cxx:2598-2628`: the transition of `FOR` against the 2-D
                // matter of `FCX` instead of the 3-D matter of `FS`.
                let mut oeg_fcxf = fcx.oeg_fcxf;
                if self.o_for != fcx.o_fcx {
                    oeg_fcxf = oeg_fcxf.complement();
                }
                let t1 = crate::topopebrep::ds::transition::Transition::from_orientation(oeg_fcxf);
                let t_aft = t1.after_state();
                let t_bef = t1.before_state();
                if oeg_fcxf == M_FORWARD {
                    b3d = self.tfe_aft == sta_fcx;
                } else if oeg_fcxf == M_REVERSED {
                    b3d = self.tfe_bef == sta_fcx;
                }
                if oeg_fcxf == M_FORWARD {
                    b2d = t_aft == State::In;
                } else if oeg_fcxf == M_REVERSED {
                    b2d = t_bef == State::In;
                }
                b = b3d && b2d;
            }
            if force_keep {
                // `cxx:2629-2632` (newnew)
                b = true;
            }
        }

        // `cxx:2633-2708`: `FOR` and `FCX` do not share their geometric
        // domain; the split ON is never kept.
        if do_ {
            if self.ope_fus {
                // `cxx:2635-2638` (cto001D2 spOU(f18))
                b = !share_g;
            }
            if self.ope_cut {
                let rank_ref = if self.opec12 { 1 } else { 2 };
                let mut binin = false;
                let mut binou = false;
                if self.tb1 == self.tb2 {
                    // `cxx:2643-2646`: the IN/IN part.
                    binin = self.rank_eg == rank_ref;
                } else {
                    // `cxx:2647-2649` (xpu291098 cto001D3 (f17ou))
                    binou = self.rank_for == 1 && self.tb_for == State::Out && !share_g;
                    if binou {
                        // `cxx:2650-2700` (xpu181198 CTS21302(FOR22,EG9,FCX5))
                        let Some(esd) = bds.shape(self.iesd, true) else {
                            return;
                        };
                        let mut lfcfk: Vec<TopoShape> = Vec::new();
                        {
                            let mut connex = Connexity::new();
                            connex.prepare(&hds);
                            connex.face_edge_connex_faces(&self.for_, &esd, &hds, &mut lfcfk);
                        }
                        if lfcfk.len() != 1 {
                            // `cxx:2657-2660`
                            return;
                        }
                        let Some(par_esd) =
                            fun_tool_par_e_tol(&self.esp_on, self.par_on, &esd, self.tol_eg)
                        else {
                            return;
                        };
                        let oofor = lfcfk[0].clone();
                        // `cxx:2666-2670`: OCCT ignores the `FUN_tool_paronEF`
                        // result, so a failed projection leaves the default
                        // `(0, 0)` point.
                        let uv = fun_tool_par_on_e_f(&esd, par_esd, &oofor)
                            .unwrap_or(GpPnt2d::new(0.0, 0.0));
                        let Some(mut nt_oofor) = fun_tool_ng_geom_f(&uv, &Face(oofor.clone()))
                        else {
                            return;
                        };
                        if oofor.orientation() == M_REVERSED {
                            nt_oofor.reverse();
                        }
                        // `cxx:2679-2684`
                        let Some(par_eg) = fun_tool_par_e(&esd, par_esd, &self.eg) else {
                            return;
                        };
                        let Some(xx_fcx) = fun_tool_getxx(&fcx.fcx, &self.eg, par_eg) else {
                            return;
                        };
                        // `cxx:2691-2700`
                        let dot = nt_oofor.dot(&xx_fcx);
                        if dot.abs() < self.tol_a {
                            binou = true; // nyixpu181198
                        } else if dot < 0.0 {
                            binou = false;
                        }
                    }
                }
                b = binin || binou;
            }
        }

        if !b {
            // `cxx:2701-2704`
            return;
        }

        // `cxx:2706-2755`
        let mut newori_e = self.tfe.orientation(self.tb1, ShapeType::Face);
        if fcx.fcx_is_ref && !self.eg_bound_for {
            newori_e = orient_ein_f_forward_or_external(&self.eg, &fcx.fcx);
            let rev = Builder::reverse(sta_fcx, fcx.sta_for);
            newori_e = Builder::orient(newori_e, rev);
        } else if self.for_is_ref && !self.eg_bound_for {
            // `cxx:2719-2731` (xpu280798 cto902C1)
            let Some(esd) = bds.shape(self.iesd, true) else {
                return;
            };
            newori_e = orient_ein_f_forward_or_external(&esd, &self.for_);
            let reverse = (newori_e == M_FORWARD || newori_e == M_REVERSED) && !self.e_on_so_esd;
            if reverse {
                newori_e = newori_e.complement();
            }
            let rev = Builder::reverse(fcx.sta_for, sta_fcx);
            newori_e = Builder::orient(newori_e, rev);
        } else if !self.eg_bound_for {
            // `cxx:2732-2754` (xpu210898 CTS20875)
            let rankref_for = bds.ancestor_rank(self.iref_for);
            if rankref_for == fcx.rank_fcx {
                newori_e = orient_ein_f_forward_or_external(&self.eg, &fcx.fcx);
            } else if rankref_for == self.rank_for {
                let Some(esd) = bds.shape(self.iesd, true) else {
                    return;
                };
                newori_e = orient_ein_f_forward_or_external(&esd, &self.for_);
                let reverse =
                    (newori_e == M_FORWARD || newori_e == M_REVERSED) && !self.e_on_so_esd;
                if reverse {
                    newori_e = newori_e.complement();
                }
            }
            let rev = Builder::reverse(fcx.sta_for, sta_fcx);
            newori_e = Builder::orient(newori_e, rev);
        }

        // `cxx:2756-2766`
        let mut new_e = self.esp_on.clone();
        pb.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);

        // `cxx:2768-2808` (xpu290798 cto 902 C1): the split IN and OUT parts
        // generated on the reference surface must be differently oriented, and
        // they share `eON`, which is `FORWARD` in one of them.
        let mut add_twice = so && !share_g && self.ope_cut;
        add_twice = add_twice && (newori_e == M_FORWARD || newori_e == M_REVERSED);
        if add_twice {
            let f_to_sp_in = if self.tb_for == State::In {
                self.for_.clone()
            } else {
                fcx.fcx.clone()
            };
            add_twice = false;
            let mut ex = Explorer::new(&f_to_sp_in, ShapeType::Edge);
            while ex.more() {
                let ee = ex.current().clone();
                ex.next();
                if !bds.has_shape(&ee, false) {
                    continue;
                }
                let issp = pb.is_split(&ee, State::In);
                if issp {
                    add_twice = !pb.splits(&ee, State::In).is_empty();
                }
                if add_twice {
                    break;
                }
            }
        }
        if add_twice {
            newori_e = newori_e.complement();
            pb.my_build_tool.set_orientation(&mut new_e, newori_e);
            wes.add_start_element(&new_e);
        }
    }
}
