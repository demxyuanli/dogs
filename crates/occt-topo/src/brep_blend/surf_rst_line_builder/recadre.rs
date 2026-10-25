//! `BRepBlend_SurfRstLineBuilder::Recadre` — the three reframing overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstLineBuilder.cxx:1002-1218`
//! (OCCT 8.0.0).
//!
//! OCCT overloads one name over three argument lists; the port splits them:
//! [`SurfRstLineBuilder::recadre_on_curv`] (`Blend_SurfCurvFuncInv`, `X = (w,
//! wguide, wcurv)`), [`SurfRstLineBuilder::recadre_surf_rst`]
//! (`Blend_SurfRstFunction` + `Blend_FuncInv`, `X = (wrst, wguide, U, V)`) and
//! [`SurfRstLineBuilder::recadre_on_point`] (`Blend_SurfPointFuncInv`, `X =
//! (w, U, V)`).
//!
//! Each returns `false` where OCCT returns `false`; the `Vtx` output is only
//! assigned when a vertex of the arc matches the solved parameter, and the
//! trailing `if (!domain->MoreVertex()) IsVtx = false` of OCCT is kept, so a
//! vertex found on the last element of the iterator is discarded.

use occt_core::gp::GpPnt2d;
use occt_math::MathVector;

use crate::blend::{BlendFuncInv, BlendSurfCurvFuncInv, BlendSurfPointFuncInv, BlendSurfRstFunction};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::BlendTool;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::FaceState;

use super::SurfRstLineBuilder;

impl SurfRstLineBuilder {
    /// `Recadre(FinvC, Solinv, Arc, IsVtx, Vtx)` (`cxx:1002-1134`).
    pub(super) fn recadre_on_curv(
        &mut self,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        sol_inv: &mut MathVector,
        arc: &mut Option<Curve2dAdaptor>,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let sol = self.sol.clone();
        let mut pt2d = GpPnt2d::new(0.0, 0.0);
        let mut last_pt2d = GpPnt2d::new(0.0, 0.0);
        let mut pmin = 0.0;

        let mut index_sol = self.arc_to_recadre(&sol, 0, &mut last_pt2d, &mut pt2d, &mut pmin);

        *is_vtx = false;
        if index_sol == 0 {
            return false;
        }

        let Some(arc_val) = self.domain_arc(index_sol) else {
            return false;
        };
        *arc = Some(arc_val.clone());
        finv_c.set(super::super::walking::arc_curve(&arc_val));

        let mut toler = MathVector::with_init(1, 3, 0.0);
        let mut infb = MathVector::with_init(1, 3, 0.0);
        let mut supb = MathVector::with_init(1, 3, 0.0);
        // The reduced tolerance is deliberate: `cxx:1034-1037` explains it with
        // a test case whose model is scaled by 0.1.
        finv_c.get_tolerance(&mut toler, 0.1 * self.tolpoint3d);
        finv_c.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, sol.value(3));
        sol_inv.set_value(3, pmin);

        let mut recadre = false;
        let start = sol_inv.clone();
        if let Some(root) = super::super::walking::solve_root_set(finv_c, &start, &toler, &infb, &supb, 30)
        {
            *sol_inv = root;
            recadre = finv_c.is_solution(sol_inv, self.tolpoint3d);
        }

        // In case of failure another arc is tried, for the case of an exit
        // close to a vertex (`cxx:1059-1101`).
        if !recadre {
            index_sol = self.arc_to_recadre(&sol, index_sol, &mut last_pt2d, &mut pt2d, &mut pmin);
            if index_sol == 0 {
                return false;
            }
            let Some(arc_val) = self.domain_arc(index_sol) else {
                return false;
            };
            *arc = Some(arc_val.clone());
            finv_c.set(super::super::walking::arc_curve(&arc_val));

            finv_c.get_tolerance(&mut toler, self.tolpoint3d);
            finv_c.get_bounds(&mut infb, &mut supb);

            sol_inv.set_value(3, pmin);

            let start = sol_inv.clone();
            if let Some(root) =
                super::super::walking::solve_root_set(finv_c, &start, &toler, &infb, &supb, 30)
            {
                *sol_inv = root;
                recadre = finv_c.is_solution(sol_inv, self.tolpoint3d);
            }
        }

        if !recadre {
            return false;
        }

        let w = sol_inv.value(2);
        if w < self.rst.first_parameter() - toler.value(2)
            || w > self.rst.last_parameter() + toler.value(2)
        {
            return false;
        }

        let Some(arc_val) = arc.clone() else {
            return false;
        };
        let tool = self.domain1.clone();
        let mut tool = tool.borrow_mut();
        tool.initialize_curve(&arc_val);
        tool.init_vertex_iterator();
        *is_vtx = !tool.more_vertex();
        while !*is_vtx {
            match tool.vertex() {
                Some(v) => {
                    if (BlendTool::parameter(&v, &arc_val) - sol_inv.value(3)).abs()
                        <= BlendTool::tolerance(&v, &arc_val)
                    {
                        *vtx = Some(v);
                        *is_vtx = true;
                    } else {
                        tool.next_vertex();
                        *is_vtx = !tool.more_vertex();
                    }
                }
                None => {
                    tool.next_vertex();
                    *is_vtx = !tool.more_vertex();
                }
            }
        }
        if !tool.more_vertex() {
            *is_vtx = false;
        }
        true
    }

    /// `Recadre(Func, Finv, Solinv, IsVtx, Vtx)` (`cxx:1138-1218`).
    pub(super) fn recadre_surf_rst(
        &mut self,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let mut toler = MathVector::with_init(1, 4, 0.0);
        let mut infb = MathVector::with_init(1, 4, 0.0);
        let mut supb = MathVector::with_init(1, 4, 0.0);
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.sol.value(3));
        sol_inv.set_value(2, self.param);
        sol_inv.set_value(3, self.sol.value(1));
        sol_inv.set_value(4, self.sol.value(2));

        let start = sol_inv.clone();
        let Some(root) = super::super::walking::solve_root_set(finv, &start, &toler, &infb, &supb, 30)
        else {
            return false;
        };
        *sol_inv = root;

        if !finv.is_solution(sol_inv, self.tolpoint3d) {
            return false;
        }

        let p2d = GpPnt2d::new(sol_inv.value(3), sol_inv.value(4));
        let situ = self
            .domain1
            .borrow_mut()
            .classify(p2d, self.tolpoint2d, false);
        if situ != FaceState::In && situ != FaceState::On {
            return false;
        }

        let rst = self.rst.clone();
        let tool = self.domain2.clone();
        {
            let mut tool = tool.borrow_mut();
            tool.initialize_curve(&rst);
            tool.init_vertex_iterator();
            *is_vtx = !tool.more_vertex();
            while !*is_vtx {
                match tool.vertex() {
                    Some(v) => {
                        if (BlendTool::parameter(&v, &rst) - sol_inv.value(1)).abs()
                            <= BlendTool::tolerance(&v, &rst)
                        {
                            *vtx = Some(v);
                            *is_vtx = true;
                        } else {
                            tool.next_vertex();
                            *is_vtx = !tool.more_vertex();
                        }
                    }
                    None => {
                        tool.next_vertex();
                        *is_vtx = !tool.more_vertex();
                    }
                }
            }
            if !tool.more_vertex() {
                *is_vtx = false;
            }
        }

        // The section is recalculated by direct resolution, otherwise
        // incoherences between the parameter and the ground caused by the
        // departure are returned (`cxx:1191-1195`).
        let mut infbound = MathVector::with_init(1, 3, 0.0);
        let mut supbound = MathVector::with_init(1, 3, 0.0);
        let mut parinit = MathVector::with_init(1, 3, 0.0);
        let mut tolerance = MathVector::with_init(1, 3, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        parinit.set_value(1, sol_inv.value(3));
        parinit.set_value(2, sol_inv.value(4));
        parinit.set_value(3, sol_inv.value(1));
        func.set(sol_inv.value(2));
        let Some(root) =
            super::super::walking::solve_root_set(func, &parinit, &tolerance, &infbound, &supbound, 30)
        else {
            return false;
        };
        sol_inv.set_value(3, root.value(1));
        sol_inv.set_value(4, root.value(2));
        sol_inv.set_value(1, root.value(3));
        true
    }

    /// `Recadre(FinvP, Solinv, IsVtx, Vtx)` (`cxx:1221-1289`).
    pub(super) fn recadre_on_point(
        &mut self,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let firstrst = self.rst.first_parameter();
        let lastrst = self.rst.last_parameter();
        let mut wpoint = firstrst;
        if (self.sol.value(3) - firstrst) > (lastrst - self.sol.value(3)) {
            wpoint = lastrst;
        }
        let Some(p2drst) = self.rst.d0(wpoint) else {
            return false;
        };
        let thepoint = self.surf2.value(p2drst.x(), p2drst.y());

        finv_p.set(thepoint);
        let mut toler = MathVector::with_init(1, 3, 0.0);
        let mut infb = MathVector::with_init(1, 3, 0.0);
        let mut supb = MathVector::with_init(1, 3, 0.0);
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, self.sol.value(1));
        sol_inv.set_value(3, self.sol.value(2));

        let start = sol_inv.clone();
        let Some(root) =
            super::super::walking::solve_root_set(finv_p, &start, &toler, &infb, &supb, 30)
        else {
            return false;
        };
        *sol_inv = root;

        if !finv_p.is_solution(sol_inv, self.tolpoint3d) {
            return false;
        }

        let p2d = GpPnt2d::new(sol_inv.value(2), sol_inv.value(3));
        let situ = self
            .domain1
            .borrow_mut()
            .classify(p2d, self.tolpoint2d, false);
        if situ != FaceState::In && situ != FaceState::On {
            return false;
        }

        let rst = self.rst.clone();
        let tool = self.domain2.clone();
        let mut tool = tool.borrow_mut();
        tool.initialize_curve(&rst);
        tool.init_vertex_iterator();
        *is_vtx = !tool.more_vertex();
        while !*is_vtx {
            match tool.vertex() {
                Some(v) => {
                    if (BlendTool::parameter(&v, &rst) - wpoint).abs()
                        <= BlendTool::tolerance(&v, &rst)
                    {
                        *vtx = Some(v);
                        *is_vtx = true;
                    } else {
                        tool.next_vertex();
                        *is_vtx = !tool.more_vertex();
                    }
                }
                None => {
                    tool.next_vertex();
                    *is_vtx = !tool.more_vertex();
                }
            }
        }
        if !tool.more_vertex() {
            *is_vtx = false;
        }
        true
    }

    /// `domain1->Init(); while (nbarc < Index) { nbarc++; domain1->Next(); }
    /// domain1->Value()` (`cxx:1021-1029`).
    fn domain_arc(&mut self, index: i32) -> Option<Curve2dAdaptor> {
        let tool = self.domain1.clone();
        let mut tool = tool.borrow_mut();
        tool.init();
        let mut nbarc = 1;
        while nbarc < index {
            nbarc += 1;
            tool.next();
        }
        tool.value().cloned()
    }
}
