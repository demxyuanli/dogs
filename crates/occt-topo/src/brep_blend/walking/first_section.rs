//! `BRepBlend_Walking::PerformFirstSection`, both overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:311-356`
//! (the overload without recadrage) and `:359-621` (the one that recadres the
//! first section on one or both restrictions), OCCT 8.0.0.
//!
//! Translation notes:
//! - The OCCT `throw Standard_Failure("... echec")` at `cxx:609` is the port's
//!   `false`: the caller only sees a failed first section.
//! - `Ext1.SetValue(P, U, V, Tol)` is the four-argument **curve** overload
//!   (`BRepBlend_Extremity.hxx:81-84`, `.cxx:117-127`), which stores
//!   `U` in `u` and `V` in `param`; those calls are written with
//!   [`crate::brep_blend::Extremity::set_on_curve`].

use occt_core::gp::GpPnt;
use occt_math::MathVector;

use crate::blend::{BlendFunction, BlendFuncInv, BlendStatus};
use crate::brep_blend::{Extremity, Line};
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::FaceState;

use super::{arc_curve, classify_domain, solve_root_set, Walking};

impl Walking {
    /// `PerformFirstSection(Func, Pdep, ParDep, Tol3d, TolGuide, Pos1, Pos2)`
    /// (`cxx:311-356`).
    ///
    /// `Pos1` / `Pos2` are the classification of the solution on the two
    /// domains, returned even when the method answers `false` (`cxx:325-347`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section(
        &mut self,
        func: &mut dyn BlendFunction,
        pdep: f64,
        par_dep: &mut MathVector,
        tol3d: f64,
        tol_guide: f64,
        pos1: &mut FaceState,
        pos2: &mut FaceState,
    ) -> bool {
        self.iscomplete = false;
        self.comptra = false;
        self.line = Some(Line::new());
        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();

        *pos1 = FaceState::Unknown;
        *pos2 = FaceState::Unknown;

        self.param = pdep;
        func.set(self.param);

        let mut tolerance = MathVector::with_init(1, 4, 0.0);
        let mut infbound = MathVector::with_init(1, 4, 0.0);
        let mut supbound = MathVector::with_init(1, 4, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let Some(root) = solve_root_set(func, par_dep, &tolerance, &infbound, &supbound, 30) else {
            return false;
        };
        self.sol = root;
        *par_dep = self.sol.clone();
        *pos1 = classify_domain(
            &self.domain1,
            self.sol.value(1),
            self.sol.value(2),
            tolerance.value(1).min(tolerance.value(2)),
            true,
        );
        *pos2 = classify_domain(
            &self.domain2,
            self.sol.value(3),
            self.sol.value(4),
            tolerance.value(3).min(tolerance.value(4)),
            true,
        );
        if *pos1 != FaceState::In || *pos2 != FaceState::In {
            return false;
        }

        self.test_arret(func, BlendStatus::Ok, false, true, false);
        true
    }

    /// `PerformFirstSection(Func, FuncInv, Pdep, Pmax, ParDep, Tol3d, TolGuide,
    /// Psol, ParSol, RecOnS1, RecOnS2)` (`cxx:359-621`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section_recadre(
        &mut self,
        func: &mut dyn BlendFunction,
        func_inv: &mut dyn BlendFuncInv,
        pdep: f64,
        pmax: f64,
        par_dep: &MathVector,
        tol3d: f64,
        tol_guide: f64,
        rec_on_s1: bool,
        rec_on_s2: bool,
        psol: &mut f64,
        par_sol: &mut MathVector,
    ) -> bool {
        self.iscomplete = false;
        self.comptra = false;
        self.line = Some(Line::new());

        self.tolpoint3d = tol3d;
        self.tolgui = tol_guide.abs();
        self.sens = if pmax - pdep >= 0.0 { 1.0 } else { -1.0 };
        // 2% of the parameter range (`cxx:390`).
        let extrapol = (pmax - pdep).abs() / 50.0;

        self.param = pdep;
        func.set(self.param);

        let mut tolerance = MathVector::with_init(1, 4, 0.0);
        let mut infbound = MathVector::with_init(1, 4, 0.0);
        let mut supbound = MathVector::with_init(1, 4, 0.0);
        let mut solrst1 = MathVector::with_init(1, 4, 0.0);
        let mut solrst2 = MathVector::with_init(1, 4, 0.0);
        let mut ext1 = Extremity::new();
        let mut ext2 = Extremity::new();
        let mut index1 = 0usize;
        let mut index2 = 0usize;
        let mut is_vtx1 = false;
        let mut is_vtx2 = false;
        let mut vtx1: Option<HVertex> = None;
        let mut vtx2: Option<HVertex> = None;
        let mut corrected_u = 0.0f64;
        let mut corrected_v = 0.0f64;
        let mut corrected_pnt = GpPnt::new(0.0, 0.0, 0.0);
        let mut corrected_param = 0.0f64;

        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        let Some(root) = solve_root_set(func, par_dep, &tolerance, &infbound, &supbound, 30) else {
            return false;
        };
        self.sol = root;

        let mut w1 = pmax;
        let mut w2 = pmax;

        let sol_snapshot = self.sol.clone();
        let mut recad1 = false;
        if rec_on_s1 {
            match self.recadre(func_inv, true, &sol_snapshot, extrapol) {
                Some(r) => {
                    recad1 = true;
                    w1 = r.solrst.value(2);
                    solrst1 = r.solrst;
                    index1 = r.index;
                    is_vtx1 = r.is_vtx;
                    vtx1 = r.vtx;
                }
                None => recad1 = false,
            }
        }
        let mut recad2 = false;
        if rec_on_s2 {
            match self.recadre(func_inv, false, &sol_snapshot, extrapol) {
                Some(r) => {
                    recad2 = true;
                    w2 = r.solrst.value(2);
                    solrst2 = r.solrst;
                    index2 = r.index;
                    is_vtx2 = r.is_vtx;
                    vtx2 = r.vtx;
                }
                None => recad2 = false,
            }
        }

        if !recad1 && !recad2 {
            return false;
        }

        let state: BlendStatus;
        if recad1 && recad2 {
            if (w1 - w2).abs() <= self.tolgui {
                // Solution on both surfaces at once (`cxx:441-449`).
                state = BlendStatus::OnRst12;
                self.param = w1;
                par_sol.set_value(1, solrst2.value(3));
                par_sol.set_value(2, solrst2.value(4));
                par_sol.set_value(3, solrst1.value(3));
                par_sol.set_value(4, solrst1.value(4));
            } else if self.sens * (w2 - w1) < 0.0 {
                // Keep the largest one: solution on the first surface.
                state = BlendStatus::OnRst1;
                self.param = w1;
                self.arc_point(true, index1, solrst1.value(1), par_sol, 1);
                par_sol.set_value(3, solrst1.value(3));
                par_sol.set_value(4, solrst1.value(4));
            } else {
                // Solution on the second surface (`cxx:471-486`).
                state = BlendStatus::OnRst2;
                self.param = w2;
                par_sol.set_value(1, solrst2.value(3));
                par_sol.set_value(2, solrst2.value(4));
                self.arc_point(false, index2, solrst2.value(1), par_sol, 3);
            }
        } else if recad1 {
            state = BlendStatus::OnRst1;
            self.param = w1;
            self.arc_point(true, index1, solrst1.value(1), par_sol, 1);
            par_sol.set_value(3, solrst1.value(3));
            par_sol.set_value(4, solrst1.value(4));
            let the_pnt_on_rst = self.surf1.d0(par_sol.value(1), par_sol.value(2));
            if let Some((u, v, pnt, param)) = self.correct_extremity_on_one_rst(
                1,
                par_sol.value(3),
                par_sol.value(4),
                self.param,
                &the_pnt_on_rst,
            ) {
                self.to_correct_on_rst1 = true;
                corrected_u = u;
                corrected_v = v;
                corrected_pnt = pnt;
                corrected_param = param;
            }
        } else {
            // `recad2` (`cxx:521-551`).
            state = BlendStatus::OnRst2;
            self.param = w2;
            par_sol.set_value(1, solrst2.value(3));
            par_sol.set_value(2, solrst2.value(4));
            self.arc_point(false, index2, solrst2.value(1), par_sol, 3);
            let the_pnt_on_rst = self.surf2.d0(par_sol.value(3), par_sol.value(4));
            if let Some((u, v, pnt, param)) = self.correct_extremity_on_one_rst(
                2,
                par_sol.value(1),
                par_sol.value(2),
                self.param,
                &the_pnt_on_rst,
            ) {
                self.to_correct_on_rst2 = true;
                corrected_u = u;
                corrected_v = v;
                corrected_pnt = pnt;
                corrected_param = param;
            }
        }
        let _ = corrected_param;

        *psol = self.param;
        self.sol = par_sol.clone();
        func.set(self.param);
        let state = self.test_arret(func, state, false, true, false);

        match state {
            BlendStatus::OnRst1 => {
                self.make_extremity(&mut ext1, true, index1, solrst1.value(1), is_vtx1, vtx1.as_ref());
                if self.to_correct_on_rst1 {
                    ext2.set_on_curve(corrected_pnt, corrected_u, corrected_v, self.tolpoint3d);
                } else {
                    ext2.set_on_curve(
                        *self.previous_p.point_on_s2(),
                        self.sol.value(3),
                        self.sol.value(4),
                        self.tolpoint3d,
                    );
                }
            }
            BlendStatus::OnRst2 => {
                if self.to_correct_on_rst2 {
                    ext1.set_on_curve(corrected_pnt, corrected_u, corrected_v, self.tolpoint3d);
                } else {
                    ext1.set_on_curve(
                        *self.previous_p.point_on_s1(),
                        self.sol.value(1),
                        self.sol.value(2),
                        self.tolpoint3d,
                    );
                }
                self.make_extremity(&mut ext2, false, index2, solrst2.value(1), is_vtx2, vtx2.as_ref());
            }
            BlendStatus::OnRst12 => {
                self.make_extremity(&mut ext1, true, index1, solrst1.value(1), is_vtx1, vtx1.as_ref());
                self.make_extremity(&mut ext2, false, index2, solrst2.value(1), is_vtx2, vtx2.as_ref());
            }
            // `default: throw Standard_Failure` (`cxx:609`).
            _ => return false,
        }

        if let Some(line) = self.line.as_mut() {
            if self.sens < 0.0 {
                line.set_end_points(ext1, ext2);
            } else {
                line.set_start_points(ext1, ext2);
            }
        }
        true
    }

    /// `p2d = BRepBlend_HCurve2dTool::Value(recdomain->Value(), solrst(1))`
    /// followed by the `ParSol` assignment of the two surface parameters
    /// (`cxx:462-466` on the first surface, `:534-538` on the second).
    ///
    /// `first_of_pair` is the 1-based index of `ParSol` where the arc's `(x, y)`
    /// are stored: 1 for the first surface, 3 for the second.
    fn arc_point(
        &self,
        on_first: bool,
        index: usize,
        arc_param: f64,
        par_sol: &mut MathVector,
        first_of_pair: usize,
    ) {
        let tool = if on_first {
            &self.recdomain1
        } else {
            &self.recdomain2
        };
        let Some(arc) = Self::domain_arc(tool, index) else {
            return;
        };
        let Some(curve) = arc_curve(&arc) else {
            return;
        };
        let p2d = curve.d0(arc_param);
        par_sol.set_value(first_of_pair, p2d.x());
        par_sol.set_value(first_of_pair + 1, p2d.y());
    }
}
