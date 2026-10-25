//! `ChFi3d_Builder::SimulData` — the two simulation drivers that walk a
//! `BRepBlend_Line` without filling a `ChFiDS_SurfData`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:1346-1669`
//! (OCCT 8.0.0). They are the `ComputeData` loops with two differences that
//! are kept verbatim:
//!
//! - `Data` is unused (`/*Data*/` in both C++ signatures), so it is not
//!   forwarded — the extension values are *not* reset here,
//! - the minimum number of points is the caller's `NbSecMin` instead of the
//!   hard-coded `3` of `ComputeData` (`cxx:1465`, `:1305`).
//!
//! `SimulData` is declared in `ChFi3d_Builder.hxx:203-520` as a member the
//! subclasses call; it is not virtual, so it belongs to this class, not to
//! `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`.
//!
//! Divergences that are not heuristics: the same ones as [`super::compute_data`]
//! (`handle<BRepBlend_Line>&` -> `&mut Option<Rc<RefCell<Line>>>`, out
//! parameters -> `&mut`).

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
use crate::chfi_ds::ElSpine;

use super::compute_data::point_parameter;
use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `SimulData(Data, HGuide, Lin, S1, I1, S2, PC2, I2, Decroch, Func, FInv,
    /// FInvP, FInvC, PFirst, MaxStep, Fleche, TolGuide, First, Last, Soldep,
    /// NbSecMin, Inside, Appro, Forward, RecP, RecS, RecRst)`
    /// (`ChFi3d_Builder_6.cxx:1346-1501`).
    ///
    /// Simulation head of the path edge / face for the bypass of obstacles.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn simul_data_surf_rst(
        &mut self,
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
        nb_sec_min: i32,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
    ) -> bool {
        // `BRepBlend_SurfRstLineBuilder TheWalk(S1, I1, S2, PC2, I2)`
        // (`cxx:1377`).
        let mut the_walk = SurfRstLineBuilder::new(s1, i1, s2, pc2, i2);

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
        // `int Nbpnt = 0;` (`cxx:1390`); unlike `ComputeData` (`cxx:1055`, `= 1`)
        // the value is dead before the first `Lin->NbPoints()` write.
        let mut nbpnt = 0;
        let mut last_line: Option<Line> = None;
        // The starting solution is reframed if needed (`cxx:1391-1421`).
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
                // `cxx:1446-1453` ignores the result of `Complete`.
                let _ = the_walk.complete(func, finv, finv_p, finv_c, sp_last);
            }

            // `Lin = TheWalk.Line()` (`cxx:1455`).
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
            } else if (nbpnt as i32) <= nb_sec_min && again == 0 {
                again += 1;
                let u1 = point_parameter(&walk_line, 1);
                let u2 = point_parameter(&walk_line, nbpnt);
                ms = (u2 - u1) / (nb_sec_min as f64 + 1.0);
                target = targetsov;
            } else if (nbpnt as i32) <= nb_sec_min {
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

    /// `SimulData(Data, HGuide, Lin, S1, PC1, I1, Decroch1, S2, PC2, I2,
    /// Decroch2, Func, FInv1, FInvP1, FInv2, FInvP2, PFirst, MaxStep, Fleche,
    /// TolGuide, First, Last, Soldep, NbSecMin, Inside, Appro, Forward, RecP1,
    /// RecRst1, RecP2, RecRst2)` (`ChFi3d_Builder_6.cxx:1506-1669`).
    ///
    /// Simulation head of the path edge / edge for the bypass of obstacles.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn simul_data_rst_rst(
        &mut self,
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
        nb_sec_min: i32,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
    ) -> bool {
        // `BRepBlend_RstRstLineBuilder TheWalk(S1, PC1, I1, S2, PC2, I2)`
        // (`cxx:1540`).
        let mut the_walk = RstRstLineBuilder::new(s1, pc1, i1, s2, pc2, i2);

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
        let mut nbpnt = 0;
        let mut last_line: Option<Line> = None;
        // The starting solution is reframed if needed (`cxx:1544-1575`).
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
                // `cxx:1612-1618` ignores the result of `Complete`.
                let _ = the_walk.complete(func, finv1, finv_p1, finv2, finv_p2, sp_last);
            }

            // `Lin = TheWalk.Line()` (`cxx:1620`).
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
            } else if (nbpnt as i32) <= nb_sec_min && again == 0 {
                again += 1;
                let u1 = point_parameter(&walk_line, 1);
                let u2 = point_parameter(&walk_line, nbpnt);
                ms = (u2 - u1) / (nb_sec_min as f64 + 1.0);
                target = targetsov;
            } else if (nbpnt as i32) <= nb_sec_min {
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
