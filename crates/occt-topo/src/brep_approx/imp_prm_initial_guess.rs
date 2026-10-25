//! The root-search half of `ApproxInt_ImpPrmSvSurfaces`: where the 2D search
//! starts, and the raw `SeekPoint` that only reports the point it found.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_ImpPrmSvSurfaces.gxx`
//! - `SeekPoint(u1, v1, u2, v2, Point)` (`:776-861`)
//! - `FillInitialVectorOfSolution(…)` (`:864-1020`)
//!
//! Both are members of the same class as `Compute`; they live in this file so
//! that neither module exceeds the size budget, and share the
//! `ImpPrmSvSurfaces` namespace as ordinary `impl` blocks.
//!
//! `FillInitialVectorOfSolution` treats the parametric surface's range as a
//! box: a start point outside it is either shifted by whole periods (when the
//! direction is periodic) or rejected, and a start point on a boundary is
//! nudged inside by `max(Confusion, U/VResolution(Confusion))`. The nudge only
//! happens with `GetUseSolver()`, because without the solver the given
//! parameters are exact and must not move (`gxx:975-1019`).

use occt_core::precision::CONFUSION;
use occt_math::{MathFunctionSetRoot, MathVector};

use crate::approx_int::SvSurfaces;
use crate::brep_approx::ImpPrmSvSurfaces;
use crate::brep_approx::SurfaceTool;
use crate::int_tools_wline::PntOn2S;
use occt_core::kernel::geomabs::SurfaceType;

impl ImpPrmSvSurfaces {
    /// `SeekPoint(u1, v1, u2, v2, Point)` (`gxx:776-861`): one solver run from
    /// the same start point `Compute` uses, storing the found point and its
    /// four parameters. No tangent bookkeeping, no cache update of the tangent
    /// slots — only `MyPnt` and `MyHasBeenComputed`.
    pub(crate) fn seek_point(
        &mut self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        point: &mut PntOn2S,
    ) -> bool {
        let a_q_surf = self.my_zer_imp_func.isurface().clone();
        let a_p_surf = self.my_zer_imp_func.psurface().clone();

        let mut x = MathVector::new(1, 2);
        let mut born_inf = MathVector::new(1, 2);
        let mut born_sup = MathVector::new(1, 2);
        let mut tolerance = MathVector::new(1, 2);
        // --- ThePSurfaceTool::GetResolution(aPSurf,Tolerance(1),Tolerance(2));
        tolerance.set_value(1, 1.0e-8);
        tolerance.set_value(2, 1.0e-8);
        let binfu = SurfaceTool::first_u_parameter(&a_p_surf);
        let binfv = SurfaceTool::first_v_parameter(&a_p_surf);
        let bsupu = SurfaceTool::last_u_parameter(&a_p_surf);
        let bsupv = SurfaceTool::last_v_parameter(&a_p_surf);
        born_inf.set_value(1, binfu);
        born_sup.set_value(1, bsupu);
        born_inf.set_value(2, binfv);
        born_sup.set_value(2, bsupv);
        let mut translation_u = 0.0_f64;
        let mut translation_v = 0.0_f64;

        if !self.fill_initial_vector_of_solution(
            u1,
            v1,
            u2,
            v2,
            binfu,
            bsupu,
            binfv,
            bsupv,
            &mut x,
            &mut translation_u,
            &mut translation_v,
        ) {
            return false;
        }

        let mut rsnld = MathFunctionSetRoot::with_iterations(&self.my_zer_imp_func, 100);
        rsnld.set_tolerance(&tolerance);
        rsnld.perform_with_bounds(&mut self.my_zer_imp_func, &x, &born_inf, &born_sup, false);
        if !rsnld.is_done() {
            return false;
        }
        self.my_has_been_computed = true;
        x = rsnld.root();

        self.my_pnt = SurfaceTool::value(&a_p_surf, x.value(1), x.value(2));

        let (new_u1, new_v1, new_u2, new_v2) = if self.my_implicit_first {
            let new_u2 = x.value(1) - translation_u;
            let new_v2 = x.value(2) - translation_v;

            let (mut new_u1, new_v1) = a_q_surf.parameters(&self.my_pnt);
            // adjust U
            if a_q_surf.type_quadric() != SurfaceType::Plane {
                let sign = if new_u1 > u1 { -1.0 } else { 1.0 };
                while (u1 - new_u1).abs() > std::f64::consts::PI {
                    new_u1 += sign * (std::f64::consts::PI + std::f64::consts::PI);
                }
            }
            (new_u1, new_v1, new_u2, new_v2)
        } else {
            let new_u1 = x.value(1) - translation_u;
            let new_v1 = x.value(2) - translation_v;

            let (mut new_u2, new_v2) = a_q_surf.parameters(&self.my_pnt);
            // adjust U
            if a_q_surf.type_quadric() != SurfaceType::Plane {
                let sign = if new_u2 > u2 { -1.0 } else { 1.0 };
                while (u2 - new_u2).abs() > std::f64::consts::PI {
                    new_u2 += sign * (std::f64::consts::PI + std::f64::consts::PI);
                }
            }
            (new_u1, new_v1, new_u2, new_v2)
        };

        point.set_value(self.my_pnt, new_u1, new_v1, new_u2, new_v2);
        true
    }

    /// `FillInitialVectorOfSolution(u1, v1, u2, v2, binfu, bsupu, binfv, bsupv,
    /// X, TranslationU, TranslationV)` (`gxx:864-1020`).
    #[allow(clippy::too_many_arguments)]
    pub fn fill_initial_vector_of_solution(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        binfu: f64,
        bsupu: f64,
        binfv: f64,
        bsupv: f64,
        x: &mut MathVector,
        translation_u: &mut f64,
        translation_v: &mut f64,
    ) -> bool {
        let a_p_surf = self.my_zer_imp_func.psurface().clone();

        *translation_u = 0.0;
        *translation_v = 0.0;

        if self.my_implicit_first {
            if u2 < binfu - 0.0000000001 {
                if SurfaceTool::is_u_periodic(&a_p_surf) {
                    let d = SurfaceTool::u_period(&a_p_surf);
                    loop {
                        *translation_u += d;
                        if u2 + *translation_u >= binfu {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if u2 > bsupu + 0.0000000001 {
                if SurfaceTool::is_u_periodic(&a_p_surf) {
                    let d = SurfaceTool::u_period(&a_p_surf);
                    loop {
                        *translation_u -= d;
                        if u2 + *translation_u <= bsupu {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            if v2 < binfv - 0.0000000001 {
                if SurfaceTool::is_v_periodic(&a_p_surf) {
                    let d = SurfaceTool::v_period(&a_p_surf);
                    loop {
                        *translation_v += d;
                        if v2 + *translation_v >= binfv {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if v2 > bsupv + 0.0000000001 {
                if SurfaceTool::is_v_periodic(&a_p_surf) {
                    let d = SurfaceTool::v_period(&a_p_surf);
                    loop {
                        *translation_v -= d;
                        if v2 + *translation_v <= bsupv {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            x.set_value(1, u2 + *translation_u);
            x.set_value(2, v2 + *translation_v);
        } else {
            if u1 < binfu - 0.0000000001 {
                if SurfaceTool::is_u_periodic(&a_p_surf) {
                    let d = SurfaceTool::u_period(&a_p_surf);
                    loop {
                        *translation_u += d;
                        if u1 + *translation_u >= binfu {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if u1 > bsupu + 0.0000000001 {
                if SurfaceTool::is_u_periodic(&a_p_surf) {
                    let d = SurfaceTool::u_period(&a_p_surf);
                    loop {
                        *translation_u -= d;
                        if u1 + *translation_u <= bsupu {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            if v1 < binfv - 0.0000000001 {
                if SurfaceTool::is_v_periodic(&a_p_surf) {
                    let d = SurfaceTool::v_period(&a_p_surf);
                    loop {
                        *translation_v += d;
                        if v1 + *translation_v >= binfv {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            } else if v1 > bsupv + 0.0000000001 {
                if SurfaceTool::is_v_periodic(&a_p_surf) {
                    let d = SurfaceTool::v_period(&a_p_surf);
                    loop {
                        *translation_v -= d;
                        if v1 + *translation_v <= bsupv {
                            break;
                        }
                    }
                } else {
                    return false;
                }
            }
            x.set_value(1, u1 + *translation_u);
            x.set_value(2, v1 + *translation_v);
        }

        // Make a small step from the boundaries in order to avoid finding an
        // "outboundaried" solution (Rsnld -> NotDone).
        if self.get_use_solver() {
            let du = CONFUSION.max(SurfaceTool::u_resolution(&a_p_surf, CONFUSION));
            let dv = CONFUSION.max(SurfaceTool::v_resolution(&a_p_surf, CONFUSION));
            if x.value(1) - 0.0000000001 <= binfu {
                x.set_value(1, x.value(1) + du);
            }
            if x.value(1) + 0.0000000001 >= bsupu {
                x.set_value(1, x.value(1) - du);
            }
            if x.value(2) - 0.0000000001 <= binfv {
                x.set_value(2, x.value(2) + dv);
            }
            if x.value(2) + 0.0000000001 >= bsupv {
                x.set_value(2, x.value(2) - dv);
            }
        }

        true
    }
}
