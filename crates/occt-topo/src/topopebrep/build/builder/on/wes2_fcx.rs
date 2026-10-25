//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `FCX` block and the
//! `yap1` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:1605-2071`
//! (OCCT 8.0.0).
//!
//! Once the `yapc3` arm has returned, `GFillONPartsWES2` switches to the
//! connex face `FCX` of `FS` across `EG` (`cxx:1610`) and to the same-domain
//! membership of `FOR` in it. [`FcxCtx`] is that block (`cxx:1610-1787`) and
//! [`Wes2Ctx::branch_yap1`] is the first consumer (`cxx:1790-2071`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:1661-1663`, `:1783-1787`, `:1792-1806`,
//! `:2063-2067`) are not translated, and neither are the
//! `FUN_coutmess("DEBUG ...")` lines (`cxx:1880`, `:1890`, `:1910`, `:1921`,
//! `:1949`, `:1959`, `:1978`, `:1989`) — they only report an orientation that
//! is neither `FORWARD` nor `REVERSED`.

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::{Face, TopoShape};
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::export::{fds_data, fun_ds_sdm};
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::project::fun_tool_proj_pon_bounded_f;
use crate::topopebrep::tool::sc::{fsc_get_psc, fsc_state_eon_face};
use crate::topopebrep::tool::{fun_tool_orient_ein_f_forward, fun_tool_value};
use crate::topo_tools::is_same;

use super::super::Builder;
use super::helpers::{M_FORWARD, M_REVERSED};
use super::wes2::Wes2Ctx;

/// The `FCX` block of `GFillONPartsWES2` (`cxx:1610-1787`).
///
/// The fields are the C++ locals of that block. `c_fcx`, `oref_fcx` and
/// `ffin_sd` are only read by the branches that are still UNPORTED
/// (`cxx:2074-3123`), hence the `dead_code` allowance.
#[allow(dead_code)]
pub(super) struct FcxCtx {
    pub(super) fcx: TopoShape,
    pub(super) i_fcx: i32,
    pub(super) rank_fcx: i32,
    pub(super) c_fcx: crate::topopebrep::ds::kind::Config,
    pub(super) iref_fcx: i32,
    pub(super) oref_fcx: Orientation,
    pub(super) fcx_is_ref: bool,
    pub(super) o_fcx: Orientation,
    pub(super) ffin_sd_so: bool,
    pub(super) ffin_sd_do: bool,
    pub(super) ffin_sd: bool,
    pub(super) eg_bound_fcx: bool,
    /// `oegFCXF` (`cxx:1665-1666`): the orientation of `EG` in `FCX` oriented
    /// `FORWARD`, i.e. `EXTERNAL` when `EG` is not a bound of `FCX`
    /// (`FUN_tool_orientEinFFORWARD` writes its out-parameter in both cases).
    pub(super) oeg_fcxf: Orientation,
    pub(super) yap1: bool,
    pub(super) yap2: bool,
    pub(super) yap3: bool,
    pub(super) yap4: bool,
    pub(super) yap5: bool,
    pub(super) yap6: bool,
    pub(super) yap7: bool,
    pub(super) sta_for: State,
    pub(super) sta_fs: State,
}

/// `FUN_tool_orientEinFFORWARD(E, F, ori)` (`TOOL.cxx:421-439`): the
/// out-parameter is written `EXTERNAL` when `E` is not a bound of `F`, so the
/// port's `None` maps to that value rather than to "leave it alone".
pub(super) fn orient_ein_f_forward_or_external(e: &TopoShape, f: &TopoShape) -> Orientation {
    fun_tool_orient_ein_f_forward(e, f).unwrap_or(Orientation::External)
}

impl Wes2Ctx {
    /// `cxx:1605-1787`: the `FCX` block and the `yap1..yap7` predicates.
    ///
    /// `None` covers every early `return` of the block: `nlfcx == 0`
    /// (`cxx:1606-1609`) and the two `FUN_tool_value` / `FUN_tool_projPonboundedF`
    /// guards of `cxx:1670-1692`. The C++ leaves the function there, which adds
    /// no start element — the same as a `false` branch.
    pub(super) fn fcx_block(&self, b: &Builder) -> Option<FcxCtx> {
        // `cxx:1606-1609`
        if self.lfcx.is_empty() {
            return None;
        }
        let hds = b.data_structure()?;
        let bds = hds.ds();

        // `cxx:1610-1611`
        let fcx = self.lfcx[0].clone();
        let i_fcx = bds.shape_index(&fcx, true);

        // `cxx:1613-1618`: the same-domain faces of `FCX`, split by rank.
        let mut lfso: Vec<TopoShape> = Vec::new();
        let mut lfdo: Vec<TopoShape> = Vec::new();
        b.find_same_domain_same_orientation_from(&fcx, &mut lfso, &mut lfdo);
        let rank_fcx = b.shape_rank(&fcx);
        let rank_x = if rank_fcx != 0 {
            if rank_fcx == 1 {
                2
            } else {
                1
            }
        } else {
            0
        };

        // `cxx:1622-1625`
        let mut lfso1: Vec<TopoShape> = Vec::new();
        let mut lfdo1: Vec<TopoShape> = Vec::new();
        let mut lfso2: Vec<TopoShape> = Vec::new();
        let mut lfdo2: Vec<TopoShape> = Vec::new();
        b.find_same_rank(&lfso, rank_fcx, &mut lfso1);
        b.find_same_rank(&lfdo, rank_fcx, &mut lfdo1);
        b.find_same_rank(&lfso, rank_x, &mut lfso2);
        b.find_same_rank(&lfdo, rank_x, &mut lfdo2);

        // `cxx:1627-1653`
        let mut ffin_sd_so = false;
        for f in &lfso2 {
            if is_same(f, &self.for_) {
                ffin_sd_so = true;
                break;
            }
        }
        let mut ffin_sd_do = false;
        for f in &lfdo2 {
            if is_same(f, &self.for_) {
                ffin_sd_do = true;
                break;
            }
        }
        let ffin_sd = ffin_sd_so || ffin_sd_do;

        // `cxx:1655-1659`
        let o_fcx = bds.shape(i_fcx, true)?.orientation();
        let c_fcx = bds.same_domain_ori(i_fcx);
        let iref_fcx = bds.same_domain_ref_of(&fcx);
        let oref_fcx = bds.shape(iref_fcx, true)?.orientation();
        let fcx_is_ref = iref_fcx == i_fcx;

        // `cxx:1665-1668`: `oegFCXF` and its binding flag. The C++ also reads
        // `FUN_tool_orientEinF(EG, FCX, oegFCX)` here and never uses it.
        let oeg_fcxf = orient_ein_f_forward_or_external(&self.eg, &fcx);
        let eg_bound_fcx = fun_tool_orient_ein_f_forward(&self.eg, &fcx).is_some();

        // `cxx:1670-1692`: with `EG` bound to neither `FOR` nor its same-domain
        // partner, `EspON` must project onto `FOR`.
        if !self.eg_bound_for && !self.esp_on_esd {
            let pt_on = fun_tool_value(self.par_on, &self.esp_on)?;
            let (_u, _v, d) = fun_tool_proj_pon_bounded_f(&pt_on, &Face(self.for_.clone()))?;
            // xpu291098 / xpu051198 / xpu040998: `* 1.e2` on the max tolerance.
            let tol_on = self.tol_eg.max(self.tol_for) * 1.0e2;
            if !(d < tol_on) {
                return None;
            }
        }

        // `cxx:1694-1698`
        let mut yap1 = true;
        yap1 = yap1 && ffin_sd;
        yap1 = yap1 && !self.eg_bound_for;
        yap1 = yap1 && eg_bound_fcx;
        yap1 = yap1 && self.eghassd3d_on;

        // `cxx:1701-1705`
        let mut yap2 = true;
        yap2 = yap2 && ffin_sd;
        yap2 = yap2 && !self.eg_bound_for;
        yap2 = yap2 && eg_bound_fcx;
        yap2 = yap2 && !self.eghassd;

        // `cxx:1707-1718`
        let mut yap1b = true;
        yap1b = yap1b && ffin_sd;
        yap1b = yap1b && !self.eg_bound_for;
        yap1b = yap1b && eg_bound_fcx;
        yap1b = yap1b && self.eghassd3d;
        yap1b = yap1b && !self.eghassd3d_on;
        yap1b = yap1b && self.eon_fs;
        if yap1b {
            yap2 = true;
        }

        // `cxx:1719-1725`
        let mut yap6 = true;
        yap6 = yap6 && ffin_sd;
        yap6 = yap6 && !self.eg_bound_for;
        yap6 = yap6 && eg_bound_fcx;
        yap6 = yap6 && self.eghassd_on;
        yap6 = yap6 && !self.eghassd3d;

        // `cxx:1727-1740`
        let mut yap6b = true;
        yap6b = yap6b && ffin_sd;
        yap6b = yap6b && !self.eg_bound_for;
        yap6b = yap6b && eg_bound_fcx;
        yap6b = yap6b && self.eghassd;
        yap6b = yap6b && !self.eghassd_on;
        yap6b = yap6b && self.eon_fs;
        yap6b = yap6b && !self.eghassd3d;
        if yap6b {
            yap2 = true;
        }

        // `cxx:1742-1749` (CTS20205: `spOUT(f30)`, `e4 = eghassd3dON`, yet `f30`
        // is not same-domain with `lfcx(e4)`).
        let mut eghassd3fcx = self.eghassd_on;
        eghassd3fcx = eghassd3fcx && !fun_ds_sdm(bds, &self.for_, &self.fs);
        eghassd3fcx = eghassd3fcx && !ffin_sd;
        let mut yap5 = true;
        yap5 = yap5 && !self.eg_bound_for;
        yap5 = yap5 && eghassd3fcx;

        // `cxx:1751-1755`
        let mut yap3 = true;
        yap3 = yap3 && !ffin_sd;
        yap3 = yap3 && !self.eg_bound_for;
        yap3 = yap3 && self.eghassd3d_on;

        // `cxx:1757-1760`
        let mut yap4 = true;
        yap4 = yap4 && !ffin_sd;
        yap4 = yap4 && self.is_rest;
        yap4 = yap4 && !yap5;
        yap4 = yap4 && !self.hsd_for;

        // `cxx:1762-1771`
        let mut yap3b = true;
        yap3b = yap3b && !ffin_sd;
        yap3b = yap3b && !self.eg_bound_for;
        yap3b = yap3b && self.eghassd3d;
        yap3b = yap3b && !self.eghassd3d_on;
        yap3b = yap3b && self.eon_fs;
        if yap3b {
            yap4 = true;
        }

        // `cxx:1773-1778`
        let mut yap7 = true;
        yap7 = yap7 && !ffin_sd;
        yap7 = yap7 && self.is_rest;
        yap7 = yap7 && !yap5;
        yap7 = yap7 && self.hsd_for;

        // `cxx:1780-1781`
        let sta_for = self.tb1;
        let sta_fs = if self.rank_fs == self.rank_for {
            self.tb1
        } else {
            self.tb2
        };

        Some(FcxCtx {
            fcx,
            i_fcx,
            rank_fcx,
            c_fcx,
            iref_fcx,
            oref_fcx,
            fcx_is_ref,
            o_fcx,
            ffin_sd_so,
            ffin_sd_do,
            ffin_sd,
            eg_bound_fcx,
            oeg_fcxf,
            yap1,
            yap2,
            yap3,
            yap4,
            yap5,
            yap6,
            yap7,
            sta_for,
            sta_fs,
        })
    }

    /// `yap1` (`cxx:1790-2071`): `FOR` is same-domain with `FCX` across `EG`;
    /// the transition of `FOR` against the 2-D matter of `FCX` decides whether
    /// `EspON` is kept, and the start element carries the orientation of `EG`
    /// in `FOR` (or in `FCX` / `Esd` depending on the reference face).
    pub(super) fn branch_yap1(&self, fcx: &FcxCtx, b: &Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = b.data_structure() else {
            return;
        };
        let bds = hds.ds();

        // `cxx:1813-1818`
        let sta_fcx = fcx.sta_fs;
        let mut b3d = false;
        let mut b2d = false;
        let mut b3de3 = false;
        let mut b2de3 = false;

        // `cxx:1820-1821`: `e3 = le3.First()`. An empty `le3` cannot happen
        // here: `eghassd3d` is true, so `FDS_HasSameDomain3d` filled the list.
        let Some(e3) = self.le3.first().cloned() else {
            return;
        };
        let i_e3 = bds.shape_index(&e3, true);

        // `cxx:1822-1847`: the shape-shape interference `EDGE(ie3) / FACE` of
        // `FCX`, whose transition is the one of `e3`.
        let mut ssif = false;
        let mut tfe_e3 = Transition::default();
        for handle in bds.shape_interferences(&fcx.fcx, true) {
            let it = handle.read().expect("poisoned Interference lock");
            // `cxx:1829-1833`: `down_cast<TopOpeBRepDS_ShapeShapeInterference>`
            // followed by `ssi.IsNull()` -> `continue`.
            if !it.is_shape_shape() {
                continue;
            }
            let (gt1, g1, st1, _s1) = fds_data(&it);
            let cond = gt1 == Kind::Edge && st1 == Kind::Face && g1 == i_e3;
            // NYI in OCCT too (`cxx:1840`): `S1` must be a face connected to
            // `FOR` through `ie3`.
            if cond {
                ssif = true;
                tfe_e3 = it.transition();
                break;
            }
        }

        // `cxx:1850-1858`
        let mut tfe_e3_bef = tfe_e3.before_state();
        let mut tfe_e3_aft = tfe_e3.after_state();
        let mut oe3_forf = Orientation::Forward;
        if ssif {
            tfe_e3_bef = tfe_e3.before_state();
            tfe_e3_aft = tfe_e3.after_state();
            oe3_forf = orient_ein_f_forward_or_external(&e3, &self.for_);
        }

        // `cxx:1860-1926`: `FOR` same-oriented with `FCX`.
        let mut oeg_fcxf = fcx.oeg_fcxf;
        if fcx.ffin_sd_so {
            if self.o_for != fcx.o_fcx {
                oeg_fcxf = oeg_fcxf.complement();
                oe3_forf = oe3_forf.complement();
            }
            let t1 = Transition::from_orientation(oeg_fcxf);
            let t1_aft = t1.after_state();
            let t1_bef = t1.before_state();
            if oeg_fcxf == M_FORWARD {
                b3d = self.tfe_aft == sta_fcx;
            } else if oeg_fcxf == M_REVERSED {
                b3d = self.tfe_bef == sta_fcx;
            }
            if oeg_fcxf == M_FORWARD {
                b2d = t1_aft == State::In;
            } else if oeg_fcxf == M_REVERSED {
                b2d = t1_bef == State::In;
            }
            if ssif {
                let t11 = Transition::from_orientation(oe3_forf);
                let t11_aft = t11.after_state();
                let t11_bef = t11.before_state();
                if oe3_forf == M_FORWARD {
                    b3de3 = tfe_e3_aft == fcx.sta_for;
                } else if oe3_forf == M_REVERSED {
                    b3de3 = tfe_e3_bef == fcx.sta_for;
                }
                if oe3_forf == M_FORWARD {
                    b2de3 = t11_aft == State::In;
                } else if oe3_forf == M_REVERSED {
                    b2de3 = t11_bef == State::In;
                }
            }
        }

        // `cxx:1928-1992`: `FOR` differently oriented with `FCX`. The
        // complement condition is `oFOR == oFCX` here, the opposite of the
        // same-oriented arm above.
        if fcx.ffin_sd_do {
            if self.o_for == fcx.o_fcx {
                oeg_fcxf = oeg_fcxf.complement();
                oe3_forf = oe3_forf.complement();
            }
            let t2 = Transition::from_orientation(oeg_fcxf);
            let t_aft = t2.after_state();
            let t_bef = t2.before_state();
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
            if ssif {
                let t22 = Transition::from_orientation(oe3_forf);
                let t22_aft = t22.after_state();
                let t22_bef = t22.before_state();
                if oe3_forf == M_FORWARD {
                    b3de3 = tfe_e3_aft == fcx.sta_for;
                } else if oe3_forf == M_REVERSED {
                    b3de3 = tfe_e3_bef == fcx.sta_for;
                }
                if oe3_forf == M_FORWARD {
                    b2de3 = t22_aft == State::In;
                } else if oe3_forf == M_REVERSED {
                    b2de3 = t22_bef == State::In;
                }
            }
        }

        // `cxx:1994-1997`
        if !(b3d && b2d && b3de3 && b2de3) {
            return;
        }

        // `cxx:1999-2008`: `EspON` must be IN the 2-D matter of `FOR`.
        let mut psc = fsc_get_psc(&self.for_);
        let state2d = fsc_state_eon_face(&self.esp_on, 0.345, &self.for_, &mut psc);
        if state2d != State::In {
            return;
        }

        // `cxx:2010-2011`
        let tbeg = self.tb1;
        let mut newori_e = self.tfe.orientation(tbeg, ShapeType::Face);

        // `cxx:2013-2056`
        if fcx.fcx_is_ref && !self.eg_bound_for {
            newori_e = orient_ein_f_forward_or_external(&self.eg, &fcx.fcx);
            let rev = Builder::reverse(sta_fcx, fcx.sta_for);
            newori_e = Builder::orient(newori_e, rev);
        } else if self.for_is_ref && !self.eg_bound_for {
            // `cxx:2026-2037`
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
            // `cxx:2038-2055` (xpu210898)
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

        // `cxx:2058-2071`
        let mut new_e = self.esp_on.clone();
        b.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
    }
}
