//! `ChFi3d_Builder::ComputeData` — the two walking drivers that build a
//! `BRepBlend_Line` between either a surface and a restriction or two
//! restrictions.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:1009-1342`
//! (OCCT 8.0.0). The two C++ overloads become
//! [`ChFi3dBuilder::compute_data_surf_rst`] (`cxx:1009-1171`) and
//! [`ChFi3dBuilder::compute_data_rst_rst`] (`cxx:1175-1342`).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<BRepBlend_Line>& Lin` is `&mut Option<Rc<RefCell<Line>>>`: the
//!   walking driver owns its line (`Option<Line>`, see
//!   [`SurfRstLineBuilder`] / [`RstRstLineBuilder`]), so the handle OCCT copies
//!   out of `TheWalk.Line()` is reproduced as a clone of that line, taken in
//!   the same place (`cxx:1123`, `:1294`). Reading `NbPoints()` / `Point(i)`
//!   from the copy is equivalent, because the driver is dropped at the end of
//!   the call,
//! - `double& First` / `double& Last` are `&mut f64`; `Last` is also read as
//!   the `Pmax` of `Perform`, exactly like the C++ out-parameter,
//! - `bool& Decroch` is `&mut bool`, and the two flags of the edge / edge
//!   overload are `&mut bool` as well,
//! - the `Adaptor3d_TopolTool` handles are the port's
//!   `Rc<RefCell<TopolTool>>`, and `Adaptor2d_Curve2d` is a
//!   [`Curve2dAdaptor`], as in the drivers it feeds.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::Surface;
use occt_math::MathVector;

use crate::blend::{
    BlendCurvPointFuncInv, BlendFuncInv, BlendRstRstFunction, BlendSurfCurvFuncInv,
    BlendSurfPointFuncInv, BlendSurfRstFunction,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::line::Line;
use crate::brep_blend::{RstRstLineBuilder, SurfRstLineBuilder};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{ElSpine, SurfData};

use super::ChFi3dBuilder;

/// `Lin->Point(1).Parameter()` / `Lin->Point(Nbpnt).Parameter()`
/// (`cxx:1140-1141`, `:1167-1168`, `:1310-1311`, `:1338-1339`).
pub(super) fn point_parameter(line: &Line, index: usize) -> f64 {
    line.point(index).map_or(0.0, |p| p.parameter())
}

impl ChFi3dBuilder {
    /// `ComputeData(Data, HGuide, Lin, S1, I1, S2, PC2, I2, Decroch, Func,
    /// FInv, FInvP, FInvC, PFirst, MaxStep, Fleche, TolGuide, First, Last,
    /// Soldep, Inside, Appro, Forward, RecP, RecS, RecRst)`
    /// (`ChFi3d_Builder_6.cxx:1009-1171`).
    ///
    /// Head of the path edge / face for the bypass of obstacles.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn compute_data_surf_rst(
        &mut self,
        data: &mut SurfData,
        hguide: &ElSpine,
        lin: &mut Option<Rc<RefCell<Line>>>,
        s1: Arc<dyn Surface>,
        i1: Rc<RefCell<TopolTool>>,
        s2: Arc<dyn Surface>,
        pc2: Curve2dAdaptor,
        i2: Rc<RefCell<TopolTool>>,
        decroch: &mut bool,
        func: &mut dyn BlendSurfRstFunction,
        finv: &mut dyn BlendFuncInv,
        finv_p: &mut dyn BlendSurfPointFuncInv,
        finv_c: &mut dyn BlendSurfCurvFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &MathVector,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
    ) -> bool {
        // `BRepBlend_SurfRstLineBuilder TheWalk(S1, I1, S2, PC2, I2)`
        // (`cxx:1038`).
        let mut the_walk = SurfRstLineBuilder::new(s1, i1, s2, pc2, i2);

        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        let reverse = !forward || inside;
        let sp_first = hguide.first_parameter();
        let sp_last = hguide.last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let nbptmin = 3;
        let mut nbpnt = 1;
        let mut last_line: Option<Line> = None;
        // The initial solution is reframed if necessary (`cxx:1058-1088`).
        let mut par_sol = MathVector::new(1, 3);
        let mut new_first = p_first;
        if rec_p || rec_s || rec_rst {
            if !the_walk.perform_first_section(
                func,
                finv,
                finv_p,
                finv_c,
                p_first,
                target,
                soldep,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                rec_rst,
                rec_p,
                rec_s,
                &mut new_first,
                &mut par_sol,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv,
                finv_p,
                finv_c,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                self.tolapp2d,
                tol_guide,
                &par_sol,
                fleche,
                approx,
            );
            if !the_walk.is_done() {
                return false;
            }
            if reverse {
                // `cxx:1114-1121` ignores the result of `Complete`.
                let _ = the_walk.complete(func, finv, finv_p, finv_c, sp_last);
            }

            // `Lin = TheWalk.Line()` (`cxx:1123`).
            let Ok(walk_line) = the_walk.line() else {
                return false;
            };
            let walk_line = walk_line.clone();
            nbpnt = walk_line.nb_points();
            *lin = Some(Rc::new(RefCell::new(walk_line.clone())));
            last_line = Some(walk_line.clone());

            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nbptmin && again == 0 {
                again += 1;
                let u1 = point_parameter(&walk_line, 1);
                let u2 = point_parameter(&walk_line, nbpnt);
                ms = (u2 - u1) / (nbptmin as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nbptmin {
                return false;
            } else {
                again = 2;
            }
        }

        *decroch = if forward {
            the_walk.decroch_end()
        } else {
            the_walk.decroch_start()
        };
        *last = point_parameter(last_line.as_ref().unwrap(), nbpnt);
        *first = point_parameter(last_line.as_ref().unwrap(), 1);
        true
    }

    /// `ComputeData(Data, HGuide, Lin, S1, PC1, I1, Decroch1, S2, PC2, I2,
    /// Decroch2, Func, FInv1, FInvP1, FInv2, FInvP2, PFirst, MaxStep, Fleche,
    /// TolGuide, First, Last, Soldep, Inside, Appro, Forward, RecP1, RecRst1,
    /// RecP2, RecRst2)` (`ChFi3d_Builder_6.cxx:1175-1342`).
    ///
    /// Head of the path edge / edge for the bypass of obstacles.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn compute_data_rst_rst(
        &mut self,
        data: &mut SurfData,
        hguide: &ElSpine,
        lin: &mut Option<Rc<RefCell<Line>>>,
        s1: Arc<dyn Surface>,
        pc1: Curve2dAdaptor,
        i1: Rc<RefCell<TopolTool>>,
        decroch1: &mut bool,
        s2: Arc<dyn Surface>,
        pc2: Curve2dAdaptor,
        i2: Rc<RefCell<TopolTool>>,
        decroch2: &mut bool,
        func: &mut dyn BlendRstRstFunction,
        finv1: &mut dyn BlendSurfCurvFuncInv,
        finv_p1: &mut dyn BlendCurvPointFuncInv,
        finv2: &mut dyn BlendSurfCurvFuncInv,
        finv_p2: &mut dyn BlendCurvPointFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        soldep: &MathVector,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
    ) -> bool {
        // `BRepBlend_RstRstLineBuilder TheWalk(S1, PC1, I1, S2, PC2, I2)`
        // (`cxx:1208`).
        let mut the_walk = RstRstLineBuilder::new(s1, pc1, i1, s2, pc2, i2);

        data.set_first_extension_value(0.0);
        data.set_last_extension_value(0.0);

        let reverse = !forward || inside;
        let sp_first = hguide.first_parameter();
        let sp_last = hguide.last_parameter();
        let mut target = sp_last;
        if reverse {
            target = sp_first;
        }
        let targetsov = target;

        let mut ms = max_step;
        let mut again = 0;
        let nbptmin = 3;
        let mut nbpnt = 0;
        let mut last_line: Option<Line> = None;
        // The initial solution is reframed if necessary (`cxx:1228-1258`).
        let mut par_sol = MathVector::new(1, 2);
        let mut new_first = p_first;
        if rec_p1 || rec_rst1 || rec_p2 || rec_rst2 {
            if !the_walk.perform_first_section(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                p_first,
                target,
                soldep,
                self.tolapp3d,
                tol_guide,
                rec_rst1,
                rec_p1,
                rec_rst2,
                rec_p2,
                &mut new_first,
                &mut par_sol,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        while again < 2 {
            the_walk.perform(
                func,
                finv1,
                finv_p1,
                finv2,
                finv_p2,
                new_first,
                *last,
                ms,
                self.tolapp3d,
                tol_guide,
                &par_sol,
                fleche,
                approx,
            );
            if !the_walk.is_done() {
                return false;
            }
            if reverse {
                // `cxx:1285-1292` ignores the result of `Complete`.
                let _ = the_walk.complete(func, finv1, finv_p1, finv2, finv_p2, sp_last);
            }

            // `Lin = TheWalk.Line()` (`cxx:1294`).
            let Ok(walk_line) = the_walk.line() else {
                return false;
            };
            let walk_line = walk_line.clone();
            nbpnt = walk_line.nb_points();
            *lin = Some(Rc::new(RefCell::new(walk_line.clone())));
            last_line = Some(walk_line.clone());

            if nbpnt <= 1 && again == 0 {
                again += 1;
                ms /= 50.0;
                target = targetsov;
            } else if nbpnt <= nbptmin && again == 0 {
                again += 1;
                let u1 = point_parameter(&walk_line, 1);
                let u2 = point_parameter(&walk_line, nbpnt);
                ms = (u2 - u1) / (nbptmin as f64 + 1.0);
                target = targetsov;
            } else if nbpnt <= nbptmin {
                return false;
            } else {
                again = 2;
            }
        }

        if forward {
            *decroch1 = the_walk.decroch1_end();
            *decroch2 = the_walk.decroch2_end();
        } else {
            *decroch1 = the_walk.decroch1_start();
            *decroch2 = the_walk.decroch2_start();
        }
        *last = point_parameter(last_line.as_ref().unwrap(), nbpnt);
        *first = point_parameter(last_line.as_ref().unwrap(), 1);
        true
    }
}
