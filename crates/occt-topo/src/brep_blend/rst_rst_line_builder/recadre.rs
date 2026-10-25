//! `BRepBlend_RstRstLineBuilder::Recadre1` / `Recadre2` — the four reframing
//! overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstLineBuilder.cxx:1137-1448`
//! (OCCT 8.0.0).
//!
//! OCCT overloads two names over two argument lists; the port splits them:
//!
//! | OCCT | port | vector | reframes |
//! |---|---|---|---|
//! | `Recadre1(Func, Finv, ...)` `:1137-1220` | [`RstRstLineBuilder::recadre1_surf_curv`] | `(wguide, w2, w1)` | `rst1` on `surf1` |
//! | `Recadre2(Func, Finv, ...)` `:1224-1302` | [`RstRstLineBuilder::recadre2_surf_curv`] | `(wguide, w1, w2)` | `rst2` on `surf2` |
//! | `Recadre1(FinvP, ...)` `:1307-1374` | [`RstRstLineBuilder::recadre1_curv_point`] | `(wguide, w2)` | `rst1` against a point of `surf1` |
//! | `Recadre2(FinvP, ...)` `:1379-1447` | [`RstRstLineBuilder::recadre2_curv_point`] | `(wguide, w1)` | `rst2` against a point of `surf2` |
//!
//! The two `surf/curv` overloads check the solved parameter against the *other*
//! restriction's range (`rst2` range in `Recadre1` `:1169-1173`, `rst1` range in
//! `Recadre2` `:1254-1257`), scan the vertexes of the arc being reframed and
//! then re-solve the section by direct resolution (`:1196-1217`, `:1282-1299`)
//! so that the parameter and the ground stay consistent. The two `curv/point`
//! overloads instead pick the nearest end of their own restriction
//! (`:1314-1321`, `:1387-1394`) and classify the resulting parameter on the
//! opposite domain (`:1345-1353`, `:1417-1425`).
//!
//! The vertex scan is the same six-line loop in all four overloads
//! (`:1177-1195`, `:1258-1276`, `:1352-1370`, `:1425-1443`), including the
//! trailing `if (!MoreVertex()) IsVtx = false` that discards a vertex found on
//! the last element of the iterator; [`scan_vertex`] keeps that behaviour.

use occt_math::MathVector;

use crate::blend::{BlendCurvPointFuncInv, BlendRstRstFunction, BlendSurfCurvFuncInv};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::BlendTool;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;

use super::RstRstLineBuilder;

/// The vertex scan shared by the four `Recadre` overloads.
fn scan_vertex(
    tool: &mut TopolTool,
    arc: &Curve2dAdaptor,
    target: f64,
) -> (bool, Option<HVertex>) {
    let mut matched: Option<HVertex> = None;
    let mut is_vtx = !tool.more_vertex();
    while !is_vtx {
        match tool.vertex() {
            Some(v) => {
                if (BlendTool::parameter(&v, arc) - target).abs()
                    <= BlendTool::tolerance(&v, arc)
                {
                    matched = Some(v);
                    is_vtx = true;
                } else {
                    tool.next_vertex();
                    is_vtx = !tool.more_vertex();
                }
            }
            None => {
                tool.next_vertex();
                is_vtx = !tool.more_vertex();
            }
        }
    }
    if !tool.more_vertex() {
        is_vtx = false;
    }
    (is_vtx, matched)
}

impl RstRstLineBuilder {
    /// `Recadre1(Func, Finv, Solinv, IsVtx, Vtx)` (`cxx:1137-1220`).
    pub(super) fn recadre1_surf_curv(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv: &mut dyn BlendSurfCurvFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let mut toler = MathVector::with_init(1, 3, 0.0);
        let mut infb = MathVector::with_init(1, 3, 0.0);
        let mut supb = MathVector::with_init(1, 3, 0.0);
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, self.sol.value(2));
        sol_inv.set_value(3, self.sol.value(1));

        let start = sol_inv.clone();
        let Some(root) =
            super::super::walking::solve_root_set(finv, &start, &toler, &infb, &supb, 30)
        else {
            return false;
        };
        *sol_inv = root;

        // The function value has to satisfy the second restriction
        // (`cxx:1167-1173`).
        if !finv.is_solution(sol_inv, self.tolpoint3d) {
            return false;
        }
        let w = sol_inv.value(2);
        if w < self.rst2.first_parameter() - toler.value(2)
            || w > self.rst2.last_parameter() + toler.value(2)
        {
            return false;
        }

        {
            let rst1 = self.rst1.clone();
            let tool = self.domain1.clone();
            let mut tool = tool.borrow_mut();
            tool.initialize_curve(&rst1);
            tool.init_vertex_iterator();
            let (found, matched) = scan_vertex(&mut tool, &rst1, sol_inv.value(3));
            *is_vtx = found;
            if let Some(v) = matched {
                *vtx = Some(v);
            }
        }

        // The section is recalculated by direct solution, otherwise the
        // incoherences between the parameter and the ground caused by the
        // departure are returned (`cxx:1197-1199`).
        let mut infbound = MathVector::with_init(1, 2, 0.0);
        let mut supbound = MathVector::with_init(1, 2, 0.0);
        let mut parinit = MathVector::with_init(1, 2, 0.0);
        let mut tolerance = MathVector::with_init(1, 2, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        parinit.set_value(1, sol_inv.value(3));
        parinit.set_value(2, sol_inv.value(2));
        func.set(sol_inv.value(1));
        let Some(root) = super::super::walking::solve_root_set(
            func, &parinit, &tolerance, &infbound, &supbound, 30,
        ) else {
            return false;
        };
        sol_inv.set_value(2, root.value(2));
        sol_inv.set_value(3, root.value(1));
        true
    }

    /// `Recadre2(Func, Finv, Solinv, IsVtx, Vtx)` (`cxx:1224-1302`).
    pub(super) fn recadre2_surf_curv(
        &mut self,
        func: &mut dyn BlendRstRstFunction,
        finv: &mut dyn BlendSurfCurvFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let mut toler = MathVector::with_init(1, 3, 0.0);
        let mut infb = MathVector::with_init(1, 3, 0.0);
        let mut supb = MathVector::with_init(1, 3, 0.0);
        finv.get_tolerance(&mut toler, self.tolpoint3d);
        finv.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, self.sol.value(1));
        sol_inv.set_value(3, self.sol.value(2));

        let start = sol_inv.clone();
        let Some(root) =
            super::super::walking::solve_root_set(finv, &start, &toler, &infb, &supb, 30)
        else {
            return false;
        };
        *sol_inv = root;

        if !finv.is_solution(sol_inv, self.tolpoint3d) {
            return false;
        }
        let w = sol_inv.value(2);
        if w < self.rst1.first_parameter() - toler.value(2)
            || w > self.rst1.last_parameter() + toler.value(2)
        {
            return false;
        }

        {
            let rst2 = self.rst2.clone();
            let tool = self.domain2.clone();
            let mut tool = tool.borrow_mut();
            tool.initialize_curve(&rst2);
            tool.init_vertex_iterator();
            let (found, matched) = scan_vertex(&mut tool, &rst2, sol_inv.value(3));
            *is_vtx = found;
            if let Some(v) = matched {
                *vtx = Some(v);
            }
        }

        let mut infbound = MathVector::with_init(1, 2, 0.0);
        let mut supbound = MathVector::with_init(1, 2, 0.0);
        let mut parinit = MathVector::with_init(1, 2, 0.0);
        let mut tolerance = MathVector::with_init(1, 2, 0.0);
        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        parinit.set_value(1, sol_inv.value(2));
        parinit.set_value(2, sol_inv.value(3));
        func.set(sol_inv.value(1));
        let Some(root) = super::super::walking::solve_root_set(
            func, &parinit, &tolerance, &infbound, &supbound, 30,
        ) else {
            return false;
        };
        sol_inv.set_value(2, root.value(1));
        sol_inv.set_value(3, root.value(2));
        true
    }

    /// `Recadre1(FinvP, Solinv, IsVtx, Vtx)` (`cxx:1307-1374`).
    pub(super) fn recadre1_curv_point(
        &mut self,
        finv_p: &mut dyn BlendCurvPointFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        // The nearest end of the restriction is used, following the direction
        // of the processing (`cxx:1311-1321`).
        let firstrst1 = self.rst1.first_parameter();
        let lastrst1 = self.rst1.last_parameter();
        let mut upoint = firstrst1;
        if (self.sol.value(1) - firstrst1) > (lastrst1 - self.sol.value(1)) {
            upoint = lastrst1;
        }
        let Some(p2drst1) = self.rst1_point(upoint) else {
            return false;
        };
        let thepoint = self.surf1.value(p2drst1.x(), p2drst1.y());

        finv_p.set(thepoint);
        let mut toler = MathVector::with_init(1, 2, 0.0);
        let mut infb = MathVector::with_init(1, 2, 0.0);
        let mut supb = MathVector::with_init(1, 2, 0.0);
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, self.sol.value(2));

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

        let Some(p2drst2) = self.rst2_point(sol_inv.value(2)) else {
            return false;
        };
        let situ = self
            .domain2
            .borrow_mut()
            .classify(p2drst2, toler.value(2), false);
        if situ != FaceState::In && situ != FaceState::On {
            return false;
        }

        let rst1 = self.rst1.clone();
        let tool = self.domain1.clone();
        let mut tool = tool.borrow_mut();
        tool.initialize_curve(&rst1);
        tool.init_vertex_iterator();
        let (found, matched) = scan_vertex(&mut tool, &rst1, upoint);
        *is_vtx = found;
        if let Some(v) = matched {
            *vtx = Some(v);
        }
        true
    }

    /// `Recadre2(FinvP, Solinv, IsVtx, Vtx)` (`cxx:1379-1447`).
    pub(super) fn recadre2_curv_point(
        &mut self,
        finv_p: &mut dyn BlendCurvPointFuncInv,
        sol_inv: &mut MathVector,
        is_vtx: &mut bool,
        vtx: &mut Option<HVertex>,
    ) -> bool {
        let firstrst2 = self.rst2.first_parameter();
        let lastrst2 = self.rst2.last_parameter();
        let mut vpoint = firstrst2;
        if (self.sol.value(2) - firstrst2) > (lastrst2 - self.sol.value(2)) {
            vpoint = lastrst2;
        }
        let Some(p2drst2) = self.rst2_point(vpoint) else {
            return false;
        };
        let thepoint = self.surf2.value(p2drst2.x(), p2drst2.y());

        finv_p.set(thepoint);
        let mut toler = MathVector::with_init(1, 2, 0.0);
        let mut infb = MathVector::with_init(1, 2, 0.0);
        let mut supb = MathVector::with_init(1, 2, 0.0);
        finv_p.get_tolerance(&mut toler, self.tolpoint3d);
        finv_p.get_bounds(&mut infb, &mut supb);
        sol_inv.set_value(1, self.param);
        sol_inv.set_value(2, self.sol.value(1));

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

        let Some(p2drst1) = self.rst1_point(sol_inv.value(2)) else {
            return false;
        };
        let situ = self
            .domain1
            .borrow_mut()
            .classify(p2drst1, toler.value(2), false);
        if situ != FaceState::In && situ != FaceState::On {
            return false;
        }

        let rst2 = self.rst2.clone();
        let tool = self.domain2.clone();
        let mut tool = tool.borrow_mut();
        tool.initialize_curve(&rst2);
        tool.init_vertex_iterator();
        let (found, matched) = scan_vertex(&mut tool, &rst2, vpoint);
        *is_vtx = found;
        if let Some(v) = matched {
            *vtx = Some(v);
        }
        true
    }
}
