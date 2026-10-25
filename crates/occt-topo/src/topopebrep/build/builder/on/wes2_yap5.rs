//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap5` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:2917-3010`
//! (OCCT 8.0.0).
//!
//! `yap5` is the 2-D counterpart of `yap3`: `EG` is not same-domain with `FOR`
//! in 3-D, yet it has a 2-D partner `e2d` in `FOR`. Note that OCCT calls
//! `FDS_Config3d` on the 2-D partner too (`cxx:2939`) — the name of the
//! function does not follow the dimension of its argument list.
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:2919-2931`, `:2958-2962`,
//! `:2995-3006`) are not translated; the `GetcontextNOSG` guard is a debug
//! switch, as in `yap3`.

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
    /// `yap5` (`cxx:2917-3010`).
    pub(super) fn branch_yap5(&self, _fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = pb.data_structure() else {
            return;
        };
        let bds = hds.ds();

        // `cxx:2932-2938`
        let mut oe2d = Orientation::External;
        let ie2d = fun_find_e_sd(bds, &self.esp_on, &self.eg, &self.for_, &mut oe2d, 2);
        if ie2d == 0 {
            return;
        }
        let Some(e2d) = bds.shape(ie2d, true) else {
            return;
        };

        // `cxx:2939-2945`
        let Some(cf) = fds_config3d(&self.esp_on, &e2d) else {
            return;
        };

        // `cxx:2947-2955`
        let Some(par_eg) =
            fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
        else {
            return;
        };

        // `cxx:2956-2968`
        let Some(same_geom) = curves_same_oriented_point(&self.eg, par_eg, &e2d) else {
            return;
        };
        let mut oe2dk = oe2d;
        if !same_geom {
            oe2dk = oe2dk.complement();
        }

        // `cxx:2971-2988`
        let keep2d = match oe2dk {
            Orientation::Forward => self.tfe_aft == self.tb1,
            Orientation::Reversed => self.tfe_bef == self.tb1,
            Orientation::Internal => true,
            Orientation::External => false,
        };

        // `cxx:2990-3007`
        if keep2d {
            let mut newori_e = oe2d;
            if cf == Config::DiffOriented {
                newori_e = newori_e.complement();
            }
            let mut new_e = self.esp_on.clone();
            pb.my_build_tool.set_orientation(&mut new_e, newori_e);
            wes.add_start_element(&new_e);
        }
    }
}
