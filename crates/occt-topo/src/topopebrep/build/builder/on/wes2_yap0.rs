//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap00` and `yap0`
//! branches.
//!
//! Sources (`TopOpeBRepBuild_BuilderON.cxx`, OCCT 8.0.0):
//! - `yap00` (`:791-877`): the ON part replaces the edge outright; the edge is
//!   added to the `WireEdgeSet` with the orientation of `EG` in `FOR` (flipped
//!   when the curves are not same-oriented).
//! - `yap0` (`:879-1199`): the "KPart" test — when the edge is tangent to both
//!   faces (`ccx:960-1108`) the keep decision follows the curvature of the two
//!   faces in the section plane, otherwise it follows the normal of `FS`
//!   against the in-face direction `XX(FOR)`.
//!
//! The `#ifdef OCCT_DEBUG` traces (`:794-808`, `:881-892`, `:1190-1197`) are
//! not translated.

use occt_core::gp::{GpDir, GpPnt2d, GpVec};

use crate::abs::{Orientation, State};
use crate::shape::Face;
use crate::topopebrep::tool::{
    curves_same_oriented_point, fun_tool_ng_geom_f_par_tol, fun_tool_orient_ein_f_forward,
    fun_tool_par_e_tol, fun_tool_par_on_e_f_tol, fun_tool_proj_pon_f, fun_tool_value_uv,
};

use super::helpers::{
    compute_face_crvt_in_sec, fun_kpart, M_INTERNAL, M_REVERSED,
};
use super::wes2::Wes2Ctx;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::tool::SurfaceAdaptor;
use crate::topexp::Explorer;
use crate::abs::ShapeType;

impl Wes2Ctx {
    /// `add` of `cxx:809-828` and `cxx:895-912`: `TFEaft == TB` when `EG` runs
    /// forward in `FOR`, `TFEbef == TB` when it runs reversed, and `true` for
    /// an internal / external occurrence.
    fn add_for(oeff: Orientation, tfe_bef: State, tfe_aft: State, tb: State) -> bool {
        match oeff {
            Orientation::Forward => tfe_aft == tb,
            Orientation::Reversed => tfe_bef == tb,
            Orientation::Internal | Orientation::External => true,
        }
    }

    /// `yap00` (`cxx:791-877`).
    pub(super) fn branch_yap00(&self, b: &super::super::Builder, wes: &mut WireEdgeSet) {
        // `cxx:809`: `FUN_tool_orientEinFFORWARD(EG, FOR, oeff)`; the C++
        // ignores the result, so a missing edge leaves `EXTERNAL`.
        let oeff = fun_tool_orient_ein_f_forward(&self.eg, &self.for_)
            .unwrap_or(Orientation::External);
        // `cxx:810-828`
        if !Self::add_for(oeff, self.tfe_bef, self.tfe_aft, self.tb1) {
            return;
        }

        // `cxx:831-843`: `Esd` must be one of `FS`'s edges.
        let Some(hds) = b.data_structure() else {
            return;
        };
        let bds = hds.ds();
        let Some(esd) = bds.shape(self.iesd, true) else {
            return;
        };
        let mut is_on_fs = false;
        let mut ex = Explorer::new(&self.fs, ShapeType::Edge);
        while ex.more() {
            if crate::topo_tools::is_same(ex.current(), &esd) {
                is_on_fs = true;
                break;
            }
            ex.next();
        }
        if !is_on_fs {
            return;
        }

        // `cxx:848-859`
        let Some(par_eg) = fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
        else {
            return;
        };
        let Some(samegeom) = curves_same_oriented_point(&self.eg, par_eg, &self.esp_on) else {
            return;
        };
        let mut newori_e = oeff;
        if !samegeom
            && (newori_e == Orientation::Forward || newori_e == Orientation::Reversed)
        {
            newori_e = newori_e.complement();
        }

        // `cxx:863-876`
        let mut new_e = self.esp_on.clone();
        b.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
    }

    /// `yap0` (`cxx:879-1199`).
    pub(super) fn branch_yap0(&self, b: &super::super::Builder, wes: &mut WireEdgeSet) {
        let oeff = fun_tool_orient_ein_f_forward(&self.eg, &self.for_)
            .unwrap_or(Orientation::External);
        // `cxx:895-912`
        if !Self::add_for(oeff, self.tfe_bef, self.tfe_aft, self.tb1) {
            return;
        }

        // `cxx:918-926`
        let Some(par_eg) = fun_tool_par_e_tol(&self.esp_on, self.par_on, &self.eg, self.tol_eg)
        else {
            return;
        };

        // `cxx:927-947`: `ngFS`, `ngFOR`, `xxFOR`.
        let mut tol_on = self.tol_fs.max(self.tol_eg);
        tol_on *= 1.0e2;
        let Some(ng_fs) = fun_tool_ng_geom_f_par_tol(par_eg, &self.eg, &self.fs, tol_on) else {
            return;
        };
        tol_on = self.tol_for.max(self.tol_eg);
        tol_on *= 1.0e2;
        let Some(ng_for) = fun_tool_ng_geom_f_par_tol(par_eg, &self.eg, &self.for_, tol_on)
        else {
            return;
        };
        let Some(xx_for) = crate::topopebrep::tool::curvature::get_xx(
            &self.eg,
            &self.for_,
            par_eg,
            &ng_for,
        ) else {
            return;
        };

        // `cxx:948-959`
        let mut nt_fs = ng_fs;
        if self.o_fs == M_REVERSED {
            nt_fs.reverse();
        }
        let mut nt_for = ng_for;
        if self.o_for == M_REVERSED {
            nt_for.reverse();
        }
        let ntdot = nt_for.dot(&nt_fs);

        // `cxx:960-1108`: KPart, i.e. FOR tangent to FS at EG.
        let mut kpart =
            self.otfe == Orientation::External || self.otfe == M_INTERNAL;
        if kpart {
            let ktg = (1.0 - ntdot.abs()).abs() < self.tol_a * 1.0e2;
            kpart = ktg;
        }
        if kpart {
            if !self.kpart_keep(par_eg, ntdot, &xx_for, &nt_for) {
                return;
            }
        } else {
            // `cxx:1112-1139`
            let dot = nt_fs.dot(&xx_for);
            let positive = dot > 0.0;
            let mut keep = true;
            if self.ope_fus {
                keep = positive;
            }
            if self.ope_com {
                keep = !positive;
            }
            if self.ope_cut {
                if positive {
                    keep = self.tb_for == State::Out;
                } else {
                    keep = self.tb_for != State::Out;
                }
            }
            if !keep {
                return;
            }
        }

        // `cxx:1142-1197`
        let mut newori_e = fun_tool_orient_ein_f_forward(&self.eg, &self.for_)
            .unwrap_or(Orientation::External);

        if self.hsd_for {
            // `cxx:1149-1188`; the `Opefus` / `Opecom` / `Opec12` blocks are
            // empty in OCCT (they contain only `nyi` comments).
            let mut reverse = false;
            if !self.for_is_ref {
                if self.c_for == crate::topopebrep::ds::kind::Config::DiffOriented {
                    reverse = !reverse;
                }
                if self.o_for != self.oref_for {
                    reverse = !reverse;
                }
            }
            if reverse {
                newori_e = newori_e.complement();
            }
        }

        let mut new_e = self.esp_on.clone();
        b.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
    }

    /// The `if (Kpart)` block of `cxx:966-1108`: `false` is a `return` of the
    /// C++ function (no start element added).
    fn kpart_keep(
        &self,
        par_eg: f64,
        ntdot: f64,
        xx_for: &GpDir,
        nt_for: &GpDir,
    ) -> bool {
        // `cxx:971-978`
        let Some(uv_for) = fun_tool_par_on_e_f_tol(&self.eg, par_eg, &self.for_, self.tol_for)
        else {
            return false;
        };
        let Some(p_for) = fun_tool_value_uv(&uv_for, &Face(self.for_.clone())) else {
            return false;
        };
        let tol_on = self.tol_for.max(self.tol_fs).max(self.tol_eg) * 10.0;

        let st_for = SurfaceAdaptor::of_face(&Face(self.for_.clone())).get_type();
        let st_fs = SurfaceAdaptor::of_face(&Face(self.fs.clone())).get_type();
        let kpart = fun_kpart(st_for, st_fs);

        if kpart == 0 {
            // `cxx:988-1031`: different curvature in the section.
            let Some(crv_for) =
                compute_face_crvt_in_sec(&Face(self.for_.clone()), &uv_for, xx_for)
            else {
                return false;
            };
            let Some((u, v, _d)) = fun_tool_proj_pon_f(&p_for, &Face(self.fs.clone())) else {
                return false;
            };
            let uv_fs = GpPnt2d::new(u, v);
            let Some(mut crv_fs) =
                compute_face_crvt_in_sec(&Face(self.fs.clone()), &uv_fs, xx_for)
            else {
                return false;
            };
            let mut crv_for = crv_for;
            if self.o_for == M_REVERSED {
                crv_for = -crv_for;
            }
            if self.o_fs == M_REVERSED {
                crv_fs = -crv_fs;
            }
            if ntdot < 0.0 {
                crv_fs = -crv_fs;
            }
            let eps = occt_core::precision::CONFUSION;
            let abs_for = crv_for.abs();
            let abs_fs = crv_fs.abs();
            if abs_for <= eps && abs_fs <= eps {
                return false;
            }
            if abs_for > eps && abs_fs > eps {
                let tol_r = tol_on;
                let r_for = 1.0 / crv_for;
                let r_fs = 1.0 / crv_fs;
                if (r_for - r_fs).abs() <= tol_r {
                    return false;
                }
            }
        }

        // `cxx:1033-1051`: pick the side of the tangent plane from the two
        // projected points.
        let so = ntdot > 0.0;
        let mut p_for = p_for;
        p_for.translate_vec(&GpVec::new(
            xx_for.x() * tol_on * 10.0,
            xx_for.y() * tol_on * 10.0,
            xx_for.z() * tol_on * 10.0,
        ));
        let Some((u, v, _d)) = fun_tool_proj_pon_f(&p_for, &Face(self.for_.clone())) else {
            return false;
        };
        let uv_for = GpPnt2d::new(u, v);
        // `cxx:1043`: the C++ ignores the `FUN_tool_value` result, so a failed
        // evaluation keeps the previous `Pfor`.
        p_for = fun_tool_value_uv(&uv_for, &Face(self.for_.clone())).unwrap_or(p_for);
        let Some((u, v, _d)) = fun_tool_proj_pon_f(&p_for, &Face(self.fs.clone())) else {
            return false;
        };
        let uv_fs = GpPnt2d::new(u, v);
        let Some(p_fs) = fun_tool_value_uv(&uv_fs, &Face(self.fs.clone())) else {
            return false;
        };

        // `cxx:1052-1108`
        let Some(d_forfs) = GpDir::from_vec(&GpVec::from_pnts(&p_for, &p_fs)).ok() else {
            return false;
        };
        let mut keep = true;
        let dot = d_forfs.dot(nt_for);
        if so {
            let so2 = dot > 0.0;
            if self.ope_fus {
                keep = !so2;
            }
            if self.ope_com {
                keep = so2;
            }
            if self.ope_cut {
                if self.tb_for == State::Out {
                    keep = !so2;
                } else {
                    keep = so2;
                }
            }
        } else {
            let do1 = dot > 0.0;
            if self.ope_fus {
                keep = do1;
            }
            if self.ope_com {
                keep = !do1;
            }
            if self.ope_cut {
                if self.tb_for == State::Out {
                    keep = do1;
                } else {
                    keep = !do1;
                }
            }
        }
        keep
    }
}
