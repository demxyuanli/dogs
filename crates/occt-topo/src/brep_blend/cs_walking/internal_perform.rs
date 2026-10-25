//! `BRepBlend_CSWalking::InternalPerform`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_CSWalking.cxx:824-1118`
//! (OCCT 8.0.0).
//!
//! The marching loop: from `parprec` the guide parameter advances by `stepw`
//! along `sens`, `math_FunctionSetRoot` solves the section in the plane of the
//! new parameter, the domain classifies the surface point, `TestArret` accepts
//! or rejects the section, and the switch (`cxx:936-1106`) decides between
//! appending, halving the step, growing it by 1.5, or stopping.
//!
//! A step that leaves the domain does not recadre in OCCT 8.0.0: the `Recadre`
//! call is commented out (`cxx:894`) and its body sits inside a `/* */` block
//! (`cxx:589-738`), so `recad` is always `false` and `echecrecad` always `true`
//! for an `OUT`/`ON` point, which ends the line with `Blend_SamePoints`. The
//! `Blend_OnRst1` arm (`cxx:1061-1084`) is kept although unreachable.

use occt_math::MathVector;

use crate::blend::{BlendCsFunction, BlendStatus};
use crate::brep_blend::Extremity;
use crate::fclass2d::FaceState;

use super::CsWalking;

impl CsWalking {
    /// `InternalPerform(Func, Sol, Bound)` (`cxx:824-1118`).
    pub(super) fn internal_perform(
        &mut self,
        func: &mut dyn BlendCsFunction,
        mut sol: MathVector,
        bound: f64,
    ) {
        let mut stepw = self.pasmax;
        let mut parprec = self.param;

        let nbvar = func.nb_variables();
        let mut tolerance = MathVector::with_init(1, nbvar, 0.0);
        let mut infbound = MathVector::with_init(1, nbvar, 0.0);
        let mut supbound = MathVector::with_init(1, nbvar, 0.0);
        // `parinit = sol` (`cxx:854`).
        let mut parinit = sol.clone();
        // Only the unreachable recadrage path reads `solrst` (`cxx:897-924`).
        #[allow(unused_mut)]
        let mut solrst = MathVector::with_init(1, nbvar, 0.0);
        let mut exts = Extremity::new();
        let mut extc = Extremity::new();

        // `Index` / `Isvtx` / `Vtx` only feed the unreachable recadrage path
        // (`cxx:836-837`, `:846`); OCCT initialises them the same way.
        let index = 0usize;
        let mut isvtx = false;

        func.get_tolerance(&mut tolerance, self.tolpoint3d);
        func.get_bounds(&mut infbound, &mut supbound);

        self.param = parprec + self.sens * stepw;
        let mut arrive = self.sens * (self.param - bound) > 0.0;

        while !arrive {
            func.set(self.param);
            let root = super::super::walking::solve_root_set(
                func,
                &parinit,
                &tolerance,
                &infbound,
                &supbound,
                30,
            );

            let state;
            match root {
                None => state = BlendStatus::StepTooLarge,
                Some(root) => {
                    sol = root;

                    let situ = self.domain.borrow_mut().classify(
                        *func.pnt2d(),
                        tolerance.value(1).min(tolerance.value(2)),
                        true,
                    );

                    let mut w = bound;
                    // `recad = Recadre(FuncInv, sol, solrst, Index, Isvtx, Vtx)`
                    // is commented out at `cxx:894`; the body of `Recadre`
                    // (`cxx:589-738`) is inside a `/* */` block.
                    let recad = false;
                    let mut echecrecad = false;

                    if situ == FaceState::Out || situ == FaceState::On {
                        // Inverse problem on the surface: nothing to do while
                        // `Recadre` is disabled (`cxx:892-901`).
                        isvtx = false;
                        if recad {
                            w = solrst.value(2);
                        } else {
                            echecrecad = true;
                        }
                    }

                    if !echecrecad {
                        if recad {
                            // Solution on the restriction (`cxx:907-926`).
                            self.param = w;
                            {
                                let mut iter = self.domain.borrow_mut();
                                iter.init();
                                let mut nbarc = 1usize;
                                while nbarc < index {
                                    nbarc += 1;
                                    iter.next();
                                }
                                if let Some(arc) = iter.value() {
                                    if let Some(p2d) = arc.d0(solrst.value(1)) {
                                        sol.set_value(1, p2d.x());
                                        sol.set_value(2, p2d.y());
                                        if nbvar >= 3 {
                                            sol.set_value(3, solrst.value(3));
                                        }
                                    }
                                }
                            }
                            func.set(self.param);
                            state = self.test_arret(func, &sol, true, BlendStatus::OnRst1);
                        } else {
                            state = self.test_arret(func, &sol, true, BlendStatus::Ok);
                        }
                    } else {
                        // Recadrage failed: leave with coincident points
                        // (`cxx:931-935`).
                        state = BlendStatus::SamePoints;
                    }
                }
            }

            match state {
                BlendStatus::Ok => {
                    self.push_previous();
                    parinit = sol.clone();
                    parprec = self.param;

                    if self.param == bound {
                        arrive = true;
                        self.set_terminate_extremities(&mut exts, &mut extc);
                    } else {
                        self.param += self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }
                BlendStatus::StepTooLarge => {
                    stepw /= 2.0;
                    if stepw.abs() < self.tolgui {
                        self.set_terminate_extremities(&mut exts, &mut extc);
                        arrive = true;
                    } else {
                        // Cannot overrun the bound (`cxx:1010`).
                        self.param = parprec + self.sens * stepw;
                    }
                }
                BlendStatus::StepTooSmall => {
                    self.push_previous();
                    parinit = sol.clone();
                    parprec = self.param;

                    stepw = (1.5 * stepw).min(self.pasmax);
                    if self.param == bound {
                        arrive = true;
                        self.set_terminate_extremities(&mut exts, &mut extc);
                    } else {
                        self.param += self.sens * stepw;
                        if self.sens * (self.param - bound) > -self.tolgui {
                            self.param = bound;
                        }
                    }
                }
                BlendStatus::OnRst1 => {
                    self.push_previous();
                    let mut exts_here = std::mem::replace(&mut exts, Extremity::new());
                    self.make_extremity(&mut exts_here, index, solrst.value(1), isvtx, None);
                    exts = exts_here;
                    let (w, param) = (
                        self.previous_p.parameter_on_c().unwrap_or(0.0),
                        self.previous_p.parameter(),
                    );
                    extc.set_on_curve(
                        *self.previous_p.point_on_c(),
                        w,
                        param,
                        self.tolpoint3d,
                    );
                    arrive = true;
                }
                BlendStatus::SamePoints => {
                    // `std::cout << " Points confondus dans le cheminement"`
                    // (`cxx:1087`): OCCT's trace, no geometry.
                    self.set_terminate_extremities(&mut exts, &mut extc);
                    arrive = true;
                }
                _ => {}
            }

            if arrive {
                if let Some(line) = self.line.as_mut() {
                    if self.sens > 0.0 {
                        line.set_end_points(exts.clone(), extc.clone());
                    } else {
                        line.set_start_points(exts.clone(), extc.clone());
                    }
                }
            }
        }
    }

    /// `line->Append(previousP)` / `line->Prepend(previousP)` by the sign of
    /// `sens` (`cxx:948-956`, `:1025-1031`, `:1069-1075`).
    fn push_previous(&mut self) {
        let point = self.previous_p;
        if let Some(line) = self.line.as_mut() {
            if self.sens > 0.0 {
                line.append(point);
            } else {
                line.prepend(point);
            }
        }
    }

    /// `Exts.SetValue(previousP.PointOnS(), U, V, previousP.Parameter(),
    /// tolpoint3d)` and `Extc.SetValue(previousP.PointOnC(),
    /// previousP.ParameterOnC(), previousP.Parameter(), tolpoint3d)`, the pair
    /// repeated at `cxx:961-971`, `:995-998`, `:1036-1046` and `:1093-1098`.
    fn set_terminate_extremities(&self, exts: &mut Extremity, extc: &mut Extremity) {
        let (u, v) = self.previous_parameters_on_s();
        let w = self.previous_p.parameter_on_c().unwrap_or(0.0);
        let param = self.previous_p.parameter();
        exts.set_on_surface(*self.previous_p.point_on_s(), u, v, param, self.tolpoint3d);
        extc.set_on_curve(*self.previous_p.point_on_c(), w, param, self.tolpoint3d);
    }
}
