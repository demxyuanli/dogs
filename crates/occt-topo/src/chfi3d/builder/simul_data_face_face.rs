//! `ChFi3d_Builder::SimulData` — the face / face simulation driver.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:2511-2703`
//! (OCCT 8.0.0). It is the face / face `ComputeData` loop without the
//! extension half: `BRepBlend_Walking` is run until the line holds at least
//! the caller's `NbSecMin` points, and only the guideline range of the
//! periodic / inside case is reframed.
//!
//! Kept verbatim from the C++:
//!
//! - `Data` is unused (`/*Data*/` in the signature), so the extension values
//!   are not reset here,
//! - the periodic reframing only touches `HGuide` and `AdditionalHGuide`, and
//!   it does *not* call `SetOrigin` (`cxx:2612-2623`),
//! - the retry ladder of `cxx:2646-2698`, where the second retry divides the
//!   step by `100` instead of scaling it by `1 / (NbSecMin + 1)`.
//!
//! Divergences that are not heuristics: the same ones as
//! [`super::compute_data_face_face`] — `handle<ChFiDS_ElSpine>` is
//! `Arc<ElSpine>`, so the two range mutations go through `Arc::make_mut` on the
//! caller's slot and the walking driver keeps the handle it was built with.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::Surface;
use occt_math::MathVector;

use crate::blend::{BlendFuncInv, BlendFunction};
use crate::brep_blend::line::Line;
use crate::brep_blend::walking::Walking;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::ElSpine;

use super::compute_data::point_parameter;
use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `SimulData(Data, HGuide, AdditionalHGuide, Lin, S1, I1, S2, I2, Func,
    /// FInv, PFirst, MaxStep, Fleche, tolguide, First, Last, Inside, Appro,
    /// Forward, Soldep, NbSecMin, RecOnS1, RecOnS2)`
    /// (`ChFi3d_Builder_6.cxx:2511-2703`).
    ///
    /// Simulation head of the path face / face.
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn simul_data_face_face(
        &mut self,
        hguide: &mut Arc<ElSpine>,
        mut additional_hguide: Option<&mut Arc<ElSpine>>,
        lin: &mut Option<Rc<RefCell<Line>>>,
        s1: Arc<dyn Surface>,
        i1: Rc<RefCell<TopolTool>>,
        s2: Arc<dyn Surface>,
        i2: Rc<RefCell<TopolTool>>,
        func: &mut dyn BlendFunction,
        finv: &mut dyn BlendFuncInv,
        p_first: f64,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        soldep: &MathVector,
        nb_sec_min: i32,
        rec_on_s1: bool,
        rec_on_s2: bool,
    ) -> bool {
        // `BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide);
        // TheWalk.Check2d(false);` (`cxx:2535-2536`).
        let mut the_walk = Walking::new(
            Arc::clone(&s1),
            Arc::clone(&s2),
            i1,
            i2,
            Arc::clone(hguide),
        );
        the_walk.check2d(false);

        let mut ms = max_step;
        let mut tol_guide = tol_guide;
        let mut sp_first = hguide.first_parameter();
        let mut sp_last = hguide.last_parameter();
        let reverse = !forward || inside;
        // `double Target; if (reverse) Target = SpFirst; else Target = SpLast;`
        // (`cxx:2544-2552`).
        let mut target = if reverse { sp_first } else { sp_last };

        let targetsov = target;
        let mut u1sov = 0.0;
        let mut u2sov = 0.0;
        // The starting solution is reframed on demand (`cxx:2556-2582`).
        let mut par_sol = MathVector::new(1, 4);
        let mut new_first = p_first;
        if rec_on_s1 || rec_on_s2 {
            if !the_walk.perform_first_section_recadre(
                func,
                finv,
                p_first,
                target,
                soldep,
                self.tolapp3d,
                tol_guide,
                rec_on_s1,
                rec_on_s2,
                &mut new_first,
                &mut par_sol,
            ) {
                return false;
            }
        } else {
            par_sol = soldep.clone();
        }

        let mut again = 0;
        let mut line = Line::new();
        let mut nbpnt = 0usize;
        while again < 3 {
            // When the start point is inside, the path goes first to the left
            // to determine the `Last` for the periodicals (`cxx:2586-2598`).
            if again == 0 && ms < 5.0 * tol_guide {
                ms = 5.0 * tol_guide;
            } else if 5.0 * tol_guide > ms {
                tol_guide = ms / 5.0;
            }

            the_walk.perform(
                func,
                finv,
                new_first,
                target,
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
            // `Lin = TheWalk.Line()` (`cxx:2609`).
            let Ok(walk_line) = the_walk.line() else {
                return false;
            };
            line = walk_line.clone();
            *lin = Some(Rc::new(RefCell::new(line.clone())));

            if reverse {
                if hguide.is_periodic() {
                    sp_first = point_parameter(&line, 1);
                    sp_last = sp_first + hguide.period().unwrap_or(0.0);
                    {
                        let hg = Arc::make_mut(hguide);
                        hg.set_first_parameter(sp_first);
                        hg.set_last_parameter(sp_last);
                    }
                    // `if (!AdditionalHGuide.IsNull())` (`cxx:2618-2622`).
                    if let Some(additional) = additional_hguide.as_deref_mut() {
                        let ahg = Arc::make_mut(additional);
                        ahg.set_first_parameter(sp_first);
                        ahg.set_last_parameter(sp_last);
                    }
                }
                // `bool complmnt = true; if (Inside) complmnt =
                // TheWalk.Complete(Func, FInv, SpLast); ...` (`cxx:2624-2635`).
                let mut complmnt = true;
                if inside {
                    complmnt = the_walk.complete(func, finv, sp_last).unwrap_or(false);
                }
                if !complmnt {
                    return false;
                }
                // `Complete` may have appended points.
                if let Ok(walk_line) = the_walk.line() {
                    line = walk_line.clone();
                }
                *lin = Some(Rc::new(RefCell::new(line.clone())));
            }

            nbpnt = line.nb_points();
            let factor = 1.0 / (nb_sec_min as f64 + 1.0);
            if nbpnt == 0 {
                return false;
            } else if nbpnt == 1 && again == 0 {
                again += 1;
                ms *= 0.01;
                target = targetsov;
                u1sov = point_parameter(&line, 1);
                u2sov = u1sov;
            } else if (nbpnt as i32) < nb_sec_min && again == 0 {
                again += 1;
                let u1 = point_parameter(&line, 1);
                let u2 = point_parameter(&line, nbpnt);
                u1sov = u1;
                u2sov = u2;
                ms = (u2 - u1) * factor;
                target = targetsov;
            } else if (nbpnt as i32) < nb_sec_min && again == 1 {
                let u1 = point_parameter(&line, 1);
                let u2 = point_parameter(&line, nbpnt);
                if (u1 - u1sov).abs() >= tol_guide || (u2 - u2sov).abs() >= tol_guide {
                    again += 1;
                    ms /= 100.0;
                    target = targetsov;
                } else {
                    return false;
                }
            } else if (nbpnt as i32) < nb_sec_min {
                return false;
            } else {
                again = 3;
            }
        }

        // `First = Lin->Point(1).Parameter(); Last =
        // Lin->Point(Nbpnt).Parameter();` (`cxx:2700-2701`).
        *first = point_parameter(&line, 1);
        *last = point_parameter(&line, nbpnt);
        true
    }
}
