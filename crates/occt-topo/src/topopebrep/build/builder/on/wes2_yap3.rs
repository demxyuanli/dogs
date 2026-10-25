//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap3` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:2811-2914`
//! (OCCT 8.0.0).
//!
//! `yap3` is the 3-D same-domain arm of `!FFinSD(so|do)` and `eghassd3dON`: the
//! 3-D partner `e3d` of `EG` gives the start orientation, unless `EG` and `e3d`
//! are not same-oriented geometrically, in which case `e3d`'s orientation is
//! complemented (`cxx:2849-2863`). The split is kept only when the transition
//! of `FOR` against `TB1` agrees with that orientation (`cxx:2864-2884`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:2813-2825`, `:2858-2862`,
//! `:2886-2893`, `:2903-2908`) are not translated; `TopOpeBRepBuild_GetcontextNOSG`
//! is a debug switch whose `#ifdef` guard is always closed in a release build,
//! so the `samegeom` complement is unconditional.

use crate::abs::Orientation;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::export::fds_config3d;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::project_curve::fun_tool_par_e_tol;
use crate::topopebrep::tool::tangent::curves_same_oriented_point;

use super::super::Builder;
use super::helpers::fun_find_e_sd;
use super::wes2::Wes2Ctx;
use super::wes2_fcx::FcxCtx;

impl Wes2Ctx {
    /// `yap3` (`cxx:2811-2914`).
    pub(super) fn branch_yap3(&self, _fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = pb.data_structure() else {
            return;
        };
        let bds = hds.ds();

        // `cxx:2826-2832`
        let mut oe3d = Orientation::External;
        let ie3d = fun_find_e_sd(bds, &self.esp_on, &self.eg, &self.for_, &mut oe3d, 3);
        if ie3d == 0 {
            return;
        }
        let Some(e3d) = bds.shape(ie3d, true) else {
            return;
        };

        // `cxx:2834-2840`
        let Some(cf) = fds_config3d(&self.esp_on, &e3d) else {
            return;
        };

        // `cxx:2842-2847`
        let Some(par_eg) =
            fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
        else {
            return;
        };

        // `cxx:2849-2863`
        let Some(same_geom) = curves_same_oriented_point(&self.eg, par_eg, &e3d) else {
            return;
        };
        let mut oe3dk = oe3d;
        if !same_geom {
            oe3dk = oe3dk.complement();
        }

        // `cxx:2864-2884`: the state to build on `FOR` is `TB1`.
        let keep3d = match oe3dk {
            Orientation::Forward => self.tfe_aft == self.tb1,
            Orientation::Reversed => self.tfe_bef == self.tb1,
            Orientation::Internal => true,
            Orientation::External => false,
        };

        // `cxx:2895-2911`
        if keep3d {
            let mut newori_e = oe3d;
            if cf == Config::DiffOriented {
                newori_e = newori_e.complement();
            }
            let mut new_e = self.esp_on.clone();
            pb.my_build_tool.set_orientation(&mut new_e, newori_e);
            wes.add_start_element(&new_e);
        }
    }
}
