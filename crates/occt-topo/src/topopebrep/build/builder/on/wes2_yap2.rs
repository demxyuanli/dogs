//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap2` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:2074-2439`
//! (OCCT 8.0.0).
//!
//! `yap2` keeps `EspON` when `FOR` is same-domain with `FCX` but `EG` is not
//! same-domain with `FOR` (`cxx:1701-1705`) and the 2-D matter of `FCX` agrees
//! with the transition of `FOR`. When `FOR` and `FS` are not same-domain the
//! decision follows the sign of the in-face direction `XX(FS)` against the
//! normal of `FOR` (`cxx:2182-2329`); otherwise it follows `b3d && b2d`, which
//! is always `false` here because OCCT evaluates `b` before writing either
//! (`cxx:2098` vs `:2100-2180`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:2076-2088`, `:2170-2175`,
//! `:2321-2326`, `:2430-2435`) are not translated, and neither are the
//! `FUN_coutmess("DEBUG ...")` lines (`cxx:2120`, `:2130`, `:2154`, `:2165`).

use occt_core::gp::GpDir;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::Face;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::export::fun_ds_sdm;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::curvature::get_xx;
use crate::topopebrep::tool::regus::nearest_f::fun_tool_getxx;
use crate::topopebrep::tool::sc::{fsc_get_psc, fsc_state_eon_face};
use crate::topopebrep::tool::{
    fun_tool_getdxx, fun_tool_ng_geom_f, fun_tool_par_e_tol, fun_tool_par_on_e_f_tol,
};

use super::super::Builder;
use super::helpers::{fun_build_tb, M_REVERSED};
use super::wes2::Wes2Ctx;
use super::wes2_fcx::{orient_ein_f_forward_or_external, FcxCtx};

impl Wes2Ctx {
    /// `yap2` (`cxx:2074-2439`).
    /// `b3d`, `b2d` and the first `neworiE` are assigned before they are read.
    /// That is what the C++ does (`cxx:2098`, `:2349-2351`); see the note on
    /// `b`.
    #[allow(unused_assignments)]
    pub(super) fn branch_yap2(&self, fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        let hds = pb.data_structure();
        let bds = hds.as_ref().map(|h| h.ds());

        // `cxx:2095-2098`
        let sta_fcx = fcx.sta_fs;
        let mut b3d = false;
        let mut b2d = false;
        // `cxx:2098`: `b` is evaluated before the two blocks that write `b3d`
        // and `b2d`, i.e. from two `false`s, and OCCT never writes it again.
        let b = b3d && b2d;

        // `cxx:2100-2134`: `FOR` same-oriented with `FCX`.
        let mut oeg_fcxf = fcx.oeg_fcxf;
        if fcx.ffin_sd_so {
            if self.o_for != fcx.o_fcx {
                oeg_fcxf = oeg_fcxf.complement();
            }
            let t1 = Transition::from_orientation(oeg_fcxf);
            let t_aft = t1.after_state();
            let t_bef = t1.before_state();
            if oeg_fcxf == Orientation::Forward {
                b3d = self.tfe_aft == sta_fcx;
            } else if oeg_fcxf == Orientation::Reversed {
                b3d = self.tfe_bef == sta_fcx;
            }
            if oeg_fcxf == Orientation::Forward {
                b2d = t_aft == State::In;
            } else if oeg_fcxf == Orientation::Reversed {
                b2d = t_bef == State::In;
            }
        }

        // `cxx:2135-2169`: `FOR` differently oriented with `FCX`.
        if fcx.ffin_sd_do {
            if self.o_for == fcx.o_fcx {
                oeg_fcxf = oeg_fcxf.complement();
            }
            let t2 = Transition::from_orientation(oeg_fcxf);
            let t_aft = t2.after_state();
            let t_bef = t2.before_state();
            if oeg_fcxf == Orientation::Forward {
                b3d = self.tfe_aft == sta_fcx;
            } else if oeg_fcxf == Orientation::Reversed {
                b3d = self.tfe_bef == sta_fcx;
            }
            if oeg_fcxf == Orientation::Forward {
                b2d = t_aft == State::In;
            } else if oeg_fcxf == Orientation::Reversed {
                b2d = t_bef == State::In;
            }
        }

        // `cxx:2177-2180` (xpu200598)
        if self.com_of_cut {
            b3d = true;
            b2d = true;
        }

        // `cxx:2182-2329`
        let mut ng_fs: Option<GpDir> = None;
        let mut par_eg = 0.0f64;
        let for_shape = bds.as_ref().and_then(|ds| ds.shape(self.i_for, true));
        let fs_shape = bds.as_ref().and_then(|ds| ds.shape(self.i_fs, true));
        let (Some(for_shape), Some(fs_shape), Some(ds)) = (for_shape, fs_shape, bds.as_ref()) else {
            return;
        };
        let sdm = fun_ds_sdm(ds, &for_shape, &fs_shape);
        // `cxx:2186`: `isfafa` — such interferences are computed in the fafa case.
        let isfafa = sdm;
        if !isfafa {
            // `cxx:2190-2196`: `parEG`.
            let Some(par) = fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
            else {
                return;
            };
            par_eg = par;

            // `cxx:2197-2203`: `uvFS` and `ngFS`.
            let Some(uv_fs) = fun_tool_par_on_e_f_tol(&self.eg, par_eg, &self.fs, self.tol_fs)
            else {
                return;
            };
            let Some(ng) = fun_tool_ng_geom_f(&uv_fs, &Face(self.fs.clone())) else {
                return;
            };
            ng_fs = Some(ng);

            // `cxx:2205-2209`: `xxFS`.
            let Some(mut xx_fs) = get_xx(&self.eg, &self.fs, par_eg, &ng) else {
                return;
            };

            // `cxx:2210-2221`: `ntFOR`.
            let Some(uv_for) = fun_tool_par_on_e_f_tol(&self.eg, par_eg, &self.for_, self.tol_for)
            else {
                return;
            };
            let Some(mut nt_for) = fun_tool_ng_geom_f(&uv_for, &Face(self.for_.clone())) else {
                return;
            };
            if self.o_for == M_REVERSED {
                nt_for.reverse();
            }

            let mut dot = xx_fs.dot(&nt_for);

            // `cxx:2223-2246` (xpu040698, CTS20578): `FS` tangent to `FCX` on
            // `EG` — retry with a slightly translated `uvFS`.
            let nulldot = dot.abs() < self.tol_a;
            if nulldot {
                let Some(dxx) = fun_tool_getdxx(&self.fs, &self.eg) else {
                    return;
                };
                let new_uv_fs = uv_fs.translated_vec(&dxx.multiplied_scalar(0.25));
                let Some(new_ng_fs) = fun_tool_ng_geom_f(&new_uv_fs, &Face(self.fs.clone()))
                else {
                    return;
                };
                let Some(new_xx_fs) = get_xx(&self.eg, &self.fs, par_eg, &new_ng_fs) else {
                    return;
                };
                xx_fs = new_xx_fs;
                dot = xx_fs.dot(&nt_for);
            }

            // `cxx:2248-2329`
            let tbfs = fun_build_tb(pb, self.rank_fs);
            let positive = dot > 0.0;
            let mut keep = true;

            // `cxx:2250-2272`
            if fcx.ffin_sd_so {
                if self.ope_fus {
                    keep = positive;
                }
                if self.ope_com {
                    keep = true;
                }
                if self.ope_cut {
                    keep = if positive {
                        tbfs == State::Out
                    } else {
                        tbfs == State::In
                    };
                }
            }

            // `cxx:2273-2320`
            if fcx.ffin_sd_do {
                if self.ope_fus {
                    keep = positive;
                    if keep {
                        let rk_ou = if self.tb1 == State::Out { 1 } else { 2 };
                        let rkfsp_ou = if self.rank_for == rk_ou { 1 } else { 2 };
                        keep = self.rank_for == rkfsp_ou;
                    }
                }
                if self.ope_com {
                    keep = !positive;
                }
                if self.ope_cut {
                    let for_moins_fs = tbfs == State::In;
                    if positive {
                        keep = if for_moins_fs { true } else { self.com_of_cut };
                    } else {
                        keep = if !for_moins_fs { true } else { self.com_of_cut };
                    }
                }
            }

            if !keep {
                return;
            }
        } else if !b {
            // `cxx:2330-2333`: always taken, see the note on `b` above.
            return;
        }

        // `cxx:2340-2348`
        let mut psc = fsc_get_psc(&self.for_);
        if fsc_state_eon_face(&self.esp_on, 0.345, &self.for_, &mut psc) != State::In {
            return;
        }

        // `cxx:2349-2352`
        let mut newori_e = self.tfe.orientation(self.tb1, ShapeType::Face);
        newori_e = orient_ein_f_forward_or_external(&self.eg, &fcx.fcx);
        let mut ok = newori_e == Orientation::Forward || newori_e == Orientation::Reversed;

        // `cxx:2355-2360`
        if fcx.ffin_sd_so && ok && self.ope_cut {
            newori_e = newori_e.complement();
        }

        // `cxx:2363-2408`
        if fcx.ffin_sd_do && ok {
            if self.ope_cut {
                if !self.com_of_cut {
                    newori_e = newori_e.complement();
                }
            }
            if self.ope_fus {
                newori_e = newori_e.complement();
            }
            if self.ope_com {
                // `cxx:2386-2406`: the direction `XX(FCX)` against the normal of
                // `FS` decides.
                let Some(xx_fcx) = fun_tool_getxx(&fcx.fcx, &self.eg, par_eg) else {
                    return;
                };
                let Some(ng_fs) = ng_fs else {
                    return;
                };
                let mut nt_fs = ng_fs;
                if self.o_fs == M_REVERSED {
                    nt_fs.reverse();
                }
                if xx_fcx.dot(&nt_fs) > 0.0 {
                    newori_e = newori_e.complement();
                }
                ok = true;
            }
        }

        // `cxx:2410-2424` (xpu270598): `reverse = false` means `FCX` has the
        // same geometric orientation.
        let mut reverse = false;
        if !fcx.fcx_is_ref {
            if fcx.c_fcx == Config::DiffOriented {
                reverse = true;
            }
            if fcx.o_fcx != fcx.oref_fcx {
                reverse = !reverse;
            }
        }
        if reverse {
            newori_e = newori_e.complement();
        }

        // `cxx:2426-2439`
        let mut new_e = self.esp_on.clone();
        pb.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
        let _ = ok;
    }
}
