//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yapc1`, `yapc2` and
//! `yapc3` branches.
//!
//! Sources (`TopOpeBRepBuild_BuilderON.cxx`, OCCT 8.0.0):
//! - `yapc1` (`:1201-1297`): closed edge, 3-D same-domain partner found; the
//!   keep decision reads `FDS_Config3d(EspON, e3d)`.
//! - `yapc2` (`:1299-1424`): closed edge; `FUN_keepEON` decides, and an
//!   `INTERNAL` orientation of a closed edge adds the ON part twice (forward
//!   and reversed).
//! - `yapc3` (`:1426-1603`): closed edge whose same-domain partner is closed on
//!   `FOR`; when the transition is `INTERNAL`/`EXTERNAL` the side is taken from
//!   the normals of `FOR` and `FS`, otherwise from `FUN_tool_curvesSO(EspON,
//!   EG)`.
//!
//! The `#ifdef OCCT_DEBUG` traces (`:1206-1216`, `:1364-1412`, `:1431-1445`,
//! `:1564-1590`) are not translated; `TopOpeBRepBuild_GetcontextNOSG()` only
//! guards a debug switch (`cxx:1242-1252`, `:1548-1557`), so the release-mode
//! `if (!samegeom)` body is kept.

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::Face;
use crate::topopebrep::ds::export::fds_config3d;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::sc::{fsc_get_psc, fsc_state_eon_face};
use crate::topopebrep::tool::tangent::curves_same_oriented;
use crate::topopebrep::tool::{
    curves_same_oriented_point, fun_tool_ng_geom_f, fun_tool_orient_ein_f_forward,
    fun_tool_par_e_tol, fun_tool_par_f_tol, fun_tool_par_on_e_f_tol,
};

use super::helpers::{fun_keep_eon, M_EXTERNAL, M_FORWARD, M_INTERNAL, M_REVERSED};
use super::wes2::Wes2Ctx;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;

impl Wes2Ctx {
    /// `yapc1` (`cxx:1201-1297`).
    pub(super) fn branch_yapc1(&self, b: &super::super::Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = b.data_structure() else {
            return;
        };
        // `cxx:1218`: `e3d = BDS.Shape(ie3d)`.
        let Some(e3d) = hds.ds().shape(self.ie3d, true) else {
            return;
        };
        // `cxx:1220-1225`
        let Some(cf) = fds_config3d(&self.esp_on, &e3d) else {
            return;
        };

        // `cxx:1227-1246`
        let mut oe3dk = self.oe3d;
        let Some(par_eg) = fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
        else {
            return;
        };
        let Some(samegeom) = curves_same_oriented_point(&self.eg, par_eg, &e3d) else {
            return;
        };
        if !samegeom {
            oe3dk = oe3dk.complement();
        }

        // `cxx:1255-1268`
        let keep3d = match oe3dk {
            Orientation::Forward => self.tfe_aft == self.tb1,
            Orientation::Reversed => self.tfe_bef == self.tb1,
            Orientation::Internal => true,
            Orientation::External => false,
        };

        // `cxx:1275-1294`
        if keep3d {
            let mut newori_e = self.oe3d;
            if cf == Config::DiffOriented {
                newori_e = newori_e.complement();
            }
            let mut new_e = self.esp_on.clone();
            b.my_build_tool.set_orientation(&mut new_e, newori_e);
            wes.add_start_element(&new_e);
        }
    }

    /// `yapc2` (`cxx:1299-1424`).
    pub(super) fn branch_yapc2(&self, b: &super::super::Builder, wes: &mut WireEdgeSet) {
        // `cxx:1315-1325`
        let keep = if self.eg_bound_for {
            fun_keep_eon(&self.eg, &self.for_, &self.tfe, self.tb1)
        } else {
            true
        };
        if !keep {
            return;
        }

        // `cxx:1329-1336`
        let mut newori_e = self.tfe.orientation(self.tb1, ShapeType::Face);
        let giveo_ein_for = self.eg_bound_for
            && !self.is_closed_ff
            && (newori_e == M_EXTERNAL || newori_e == M_INTERNAL);
        if giveo_ein_for {
            newori_e = fun_tool_orient_ein_f_forward(&self.eg, &self.for_)
                .unwrap_or(Orientation::External);
        }

        // `cxx:1338-1347`: the 2d state of the ON part on `FOR`.
        if self.yapc2b && newori_e != M_EXTERNAL {
            let mut psc = fsc_get_psc(&self.for_);
            let state2d = fsc_state_eon_face(&self.esp_on, 0.345, &self.for_, &mut psc);
            if state2d != State::In {
                return;
            }
        }

        // `cxx:1349-1362`
        let mut add_for_rev = false;
        let mut newori_ek = Orientation::External;
        if !self.is_closed_ff {
            newori_ek = newori_e;
        } else if newori_e == M_FORWARD || newori_e == M_REVERSED {
            newori_ek = newori_e;
        } else if newori_e == M_INTERNAL {
            add_for_rev = true;
        }

        let mut new_e = self.esp_on.clone();
        b.my_build_tool.set_orientation(&mut new_e, newori_ek);

        // `cxx:1414-1423`
        if add_for_rev {
            b.my_build_tool
                .set_orientation(&mut new_e, Orientation::Forward);
            wes.add_start_element(&new_e);
            b.my_build_tool
                .set_orientation(&mut new_e, Orientation::Reversed);
            wes.add_start_element(&new_e);
        } else {
            wes.add_start_element(&new_e);
        }
    }

    /// `yapc3` (`cxx:1426-1603`).
    pub(super) fn branch_yapc3(&self, b: &super::super::Builder, wes: &mut WireEdgeSet) {
        // `cxx:1447-1455`
        let keep = if self.eg_bound_for {
            fun_keep_eon(&self.eg, &self.for_, &self.tfe, self.tb1)
        } else {
            true
        };
        if !keep {
            return;
        }

        // `cxx:1459-1465`
        let mut newori_e = self.tfe.orientation(self.tb1, ShapeType::Face);
        let giveo_ein_for = self.eg_bound_for
            && !self.is_closed_ff
            && (newori_e == M_EXTERNAL || newori_e == M_INTERNAL);
        if giveo_ein_for {
            newori_e = fun_tool_orient_ein_f_forward(&self.eg, &self.for_)
                .unwrap_or(Orientation::External);
        }

        let mut add_for_rev = false;

        // `cxx:1471-1538`
        if self.otfe == M_INTERNAL || self.otfe == M_EXTERNAL {
            let so;

            let Some(par_eg) =
                fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
            else {
                return;
            };
            let Some(uv_fs) = fun_tool_par_on_e_f_tol(&self.eg, par_eg, &self.fs, self.tol_fs)
            else {
                return;
            };
            let Some(ng_fs) = fun_tool_ng_geom_f(&uv_fs, &Face(self.fs.clone())) else {
                return;
            };
            let mut nt_fs = ng_fs.clone();
            if self.o_fs == Orientation::Reversed {
                nt_fs.reverse();
            }

            let Some(uv_for) = fun_tool_par_f_tol(&self.eg, par_eg, &self.for_, self.tol_for)
            else {
                return;
            };
            let Some(ng_for) = fun_tool_ng_geom_f(&uv_for, &Face(self.for_.clone())) else {
                return;
            };
            let mut nt_for = ng_for;
            if self.o_for == Orientation::Reversed {
                nt_for.reverse();
            }

            let dot = nt_for.dot(&nt_fs);
            let nulldot = dot.abs() < self.tol_a;
            if nulldot {
                let Some(xx_fs) = crate::topopebrep::tool::curvature::get_xx(
                    &self.eg,
                    &self.fs,
                    par_eg,
                    &ng_fs,
                ) else {
                    return;
                };
                let dot2 = nt_for.dot(&xx_fs);
                so = dot2 < 0.0;
            } else {
                so = dot > 0.0;
            }

            // `cxx:1526-1537`
            let mut unkeepclo = false;
            if so && self.otfe == M_EXTERNAL && self.ope_com {
                unkeepclo = true;
            }
            if so && self.otfe == M_INTERNAL && self.ope_fus {
                unkeepclo = true;
            }
            add_for_rev = !unkeepclo;
        }

        let mut new_e = self.esp_on.clone();
        if !add_for_rev {
            // `cxx:1539-1563`
            let mut newori_ek = newori_e;
            let Some(samegeom) = curves_same_oriented(&self.esp_on, &self.eg) else {
                return;
            };
            if !samegeom {
                newori_ek = newori_ek.complement();
            }
            b.my_build_tool.set_orientation(&mut new_e, newori_ek);
        }

        // `cxx:1592-1601`
        if add_for_rev {
            b.my_build_tool
                .set_orientation(&mut new_e, Orientation::Forward);
            wes.add_start_element(&new_e);
            b.my_build_tool
                .set_orientation(&mut new_e, Orientation::Reversed);
            wes.add_start_element(&new_e);
        } else {
            wes.add_start_element(&new_e);
        }
    }
}
