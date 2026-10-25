//! `BRepBlend_Walking::ArcToRecadre` and `Recadre` — moving a section back onto
//! a restriction curve of one of the supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Walking.cxx:1197-1600`
//! (OCCT 8.0.0).
//!
//! `ArcToRecadre` (`cxx:1197-1287`) picks the arc of the recadrage domain the
//! section has to be moved to: the arc closest to the previous point, an
//! intersection with the chord winning over a projection and breaking the
//! search when no previous index was given.
//!
//! `Recadre` (`cxx:1289-1600`) solves `Blend_FuncInv` on that arc. The initial
//! guess is the plain last solution (`cxx:1386-1397`) or, when `byinter`, an
//! interpolation along the chord (`cxx:1350-1385`) — unless the interval
//! between the two parameters contains a known singular point
//! (`cxx:1319-1353`), which is then recovered with its own periodic recadrage
//! (`cxx:1400-1458`); on failure near an end of the arc the search is retried on
//! the neighbouring arc (`cxx:1517-1598`), which is the "exit close to a vertex"
//! case.
//!
//! Translation notes:
//! - `BRepBlend_HCurve2dTool` sees an `Adaptor2d_Curve2d`, whose parameter range
//!   is the arc's, not the p-curve's; [`super::arc_curve`] rebuilds that range
//!   as a trimmed 2d curve for the [`BlendTool`] queries, which take a bare
//!   curve.
//! - The `handle(Adaptor3d_HVertex)& Vtx` out-parameter becomes
//!   `Option<HVertex>`: OCCT leaves it untouched when no vertex matches.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;
use occt_math::MathVector;

use crate::blend::{BlendFuncInv, BlendPoint};

use crate::brep_blend::BlendTool;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;

use super::Walking;

/// The out-parameters of `Recadre`: `Solrst`, `Indexsol`, `IsVtx` and `Vtx`
/// (`cxx:1289-1297`).
pub(super) struct RecadreSolution {
    pub solrst: MathVector,
    pub index: usize,
    pub is_vtx: bool,
    pub vtx: Option<HVertex>,
}

impl Walking {
    /// `ArcToRecadre(OnFirst, Sol, PrevIndex, lastpt2d, pt2d, ponarc)`
    /// (`cxx:1197-1287`). Returns `IndexSol`; `0` means "no arc".
    pub(super) fn arc_to_recadre(
        &self,
        on_first: bool,
        sol: &MathVector,
        prev_index: usize,
        last_pt2d: &mut GpPnt2d,
        pt2d: &mut GpPnt2d,
        ponarc: &mut f64,
    ) -> usize {
        let mut index_sol = 0usize;
        let mut nbarc = 0usize;
        // `okinter` is *not* reset per arc in OCCT (`cxx:1203`).
        let mut okinter = false;
        let by_inter = self.line.as_ref().map(|l| l.nb_points() != 0).unwrap_or(false);
        let mut distmin = f64::MAX;
        let mut prm = 0.0f64;
        let mut dist = 0.0f64;

        let tool = if on_first {
            Rc::clone(&self.recdomain1)
        } else {
            Rc::clone(&self.recdomain2)
        };

        let mut uprev = 0.0f64;
        let mut vprev = 0.0f64;
        if by_inter {
            let params = if on_first {
                self.previous_p.parameters_on_s1()
            } else {
                self.previous_p.parameters_on_s2()
            };
            if let Ok((u, v)) = params {
                uprev = u;
                vprev = v;
            }
        }
        if on_first {
            *pt2d = GpPnt2d::new(sol.value(1), sol.value(2));
        } else {
            *pt2d = GpPnt2d::new(sol.value(3), sol.value(4));
        }
        *last_pt2d = GpPnt2d::new(uprev, vprev);

        let mut iter = tool.borrow_mut();
        iter.init();
        while iter.more() {
            nbarc += 1;
            let mut ok = false;
            if let Some(cur) = iter.value().cloned() {
                if let Some(arc) = super::arc_curve(&cur) {
                    if by_inter {
                        ok = BlendTool::inters(pt2d, last_pt2d, arc.as_ref(), &mut prm, &mut dist);
                        okinter = ok;
                    }
                    if !ok {
                        ok = BlendTool::project(pt2d, arc.as_ref(), &mut prm, &mut dist);
                    }
                }
            }
            if ok && nbarc != prev_index && (dist < distmin || okinter) {
                distmin = dist;
                *ponarc = prm;
                index_sol = nbarc;
                if okinter && prev_index == 0 {
                    break;
                }
            }
            iter.next();
        }
        index_sol
    }

    /// `Recadre(FInv, OnFirst, Sol, Solrst, Indexsol, IsVtx, Vtx, Extrap)`
    /// (`cxx:1289-1600`): `None` when the recadrage fails.
    pub(super) fn recadre(
        &mut self,
        func_inv: &mut dyn BlendFuncInv,
        on_first: bool,
        sol: &MathVector,
        extrap: f64,
    ) -> Option<RecadreSolution> {
        let mut jalons_trouve = false;
        let mut recadre = true;
        let by_inter = self.line.as_ref().map(|l| l.nb_points() != 0).unwrap_or(false);
        let mut le_jalon = 0usize;

        let mut toler = MathVector::with_init(1, 4, 0.0);
        let mut infb = MathVector::with_init(1, 4, 0.0);
        let mut supb = MathVector::with_init(1, 4, 0.0);
        let mut solrst = MathVector::with_init(1, 4, 0.0);
        let mut prm = 0.0f64;
        let mut dist = 0.0f64;

        let tool = if on_first {
            Rc::clone(&self.recdomain1)
        } else {
            Rc::clone(&self.recdomain2)
        };

        let mut last_pt2d = GpPnt2d::zero();
        let mut pt2d = GpPnt2d::zero();
        let mut pmin = 0.0f64;
        let index_sol = self.arc_to_recadre(on_first, sol, 0, &mut last_pt2d, &mut pt2d, &mut pmin);
        if index_sol == 0 {
            return None;
        }

        let thearc = Self::domain_arc(&tool, index_sol)?;
        let thecur = super::arc_curve(&thearc)?;

        // The problem to solve.
        func_inv.set(on_first, Some(Arc::clone(&thecur)));
        func_inv.get_bounds(&mut infb, &mut supb);
        infb.set_value(2, infb.value(2) - extrap);
        supb.set_value(2, supb.value(2) + extrap);
        func_inv.get_tolerance(&mut toler, 0.1 * self.tolpoint3d);
        // The solver copies the tolerance vector; the vector is then relaxed by
        // ten for the tests (`cxx:1336-1339`).
        let solver_toler = toler.clone();
        toler.multiply_scalar(10.0);

        // First guess.
        let (ufirst, ulast) = BlendTool::bounds(thecur.as_ref());
        // Help the singular-corner search by snapping the parameter to an end.
        let span = (ulast - ufirst).abs();
        if (pmin - ufirst).abs() < span / 1000.0 {
            pmin = ufirst;
        }
        if (pmin - ulast).abs() < span / 1000.0 {
            pmin = ulast;
        }

        if by_inter {
            let last_param = self.previous_p.parameter();
            if !self.jalons.is_empty() {
                let (t1, t2) = if last_param < self.param {
                    (last_param, self.param)
                } else {
                    (self.param, last_param)
                };
                let mut cherche = true;
                let mut ii = 1usize;
                while ii <= self.jalons.len() && cherche {
                    let t = self.jalons[ii - 1].parameter();
                    if (t1 < t && t2 > t) || t == self.param {
                        jalons_trouve = true;
                        le_jalon = ii;
                        cherche = false; // Only works when both exit together.
                    } else {
                        cherche = t < t2; // Stop as soon as t >= t2.
                    }
                    ii += 1;
                }
            }
            if !jalons_trouve {
                // Initialisation by interpolation (`cxx:1435-1462`).
                let pnt = thearc.d0(pmin).unwrap_or_default();
                let (pnt1, pnt2) = if on_first {
                    let (u, v) = self.previous_p.parameters_on_s2().unwrap_or((0.0, 0.0));
                    (GpPnt2d::new(u, v), GpPnt2d::new(sol.value(3), sol.value(4)))
                } else {
                    let (u, v) = self.previous_p.parameters_on_s1().unwrap_or((0.0, 0.0));
                    (GpPnt2d::new(u, v), GpPnt2d::new(sol.value(1), sol.value(2)))
                };
                let mut lambda = pnt.distance(&last_pt2d);
                if lambda > 1.0e-12 {
                    lambda /= pnt.distance(&last_pt2d) + pnt.distance(&pt2d);
                } else {
                    lambda = 0.0;
                }
                solrst.set_value(1, pmin);
                solrst.set_value(2, (1.0 - lambda) * last_param + lambda * self.param);
                solrst.set_value(3, (1.0 - lambda) * pnt1.x() + lambda * pnt2.x());
                solrst.set_value(4, (1.0 - lambda) * pnt1.y() + lambda * pnt2.y());
            }
        } else {
            // Otherwise the last computed point is the guess (`cxx:1464-1478`).
            solrst.set_value(1, pmin);
            solrst.set_value(2, self.param);
            if on_first {
                solrst.set_value(3, sol.value(3));
                solrst.set_value(4, sol.value(4));
            } else {
                solrst.set_value(3, sol.value(1));
                solrst.set_value(4, sol.value(2));
            }
        }

        if jalons_trouve {
            // Take the singular point back (`cxx:1480-1553`).
            let mon_jalon: BlendPoint = self.jalons[le_jalon - 1];
            solrst.set_value(2, mon_jalon.parameter());
            let (p3, p4) = if on_first {
                mon_jalon.parameters_on_s2().unwrap_or((0.0, 0.0))
            } else {
                mon_jalon.parameters_on_s1().unwrap_or((0.0, 0.0))
            };
            solrst.set_value(3, p3);
            solrst.set_value(4, p4);
            let other = if on_first {
                Arc::clone(&self.surf2)
            } else {
                Arc::clone(&self.surf1)
            };

            // Eventual recadrage for the periodic case (`cxx:1500-1530`).
            if other.is_u_periodic() || other.is_v_periodic() {
                let lastpt2d = thearc.d0(pmin).unwrap_or_default();
                if other.is_u_periodic() {
                    let uperiod = other.u_period();
                    if solrst.value(3) - lastpt2d.x() > uperiod * 0.6 {
                        solrst.set_value(3, solrst.value(3) - uperiod);
                    }
                    if solrst.value(3) - lastpt2d.x() < -uperiod * 0.6 {
                        solrst.set_value(3, solrst.value(3) + uperiod);
                    }
                }
                if other.is_v_periodic() {
                    let vperiod = other.v_period();
                    if solrst.value(4) - lastpt2d.y() > vperiod * 0.6 {
                        solrst.set_value(4, solrst.value(4) - vperiod);
                    }
                    if solrst.value(4) - lastpt2d.y() < -vperiod * 0.6 {
                        solrst.set_value(4, solrst.value(4) + vperiod);
                    }
                }
            }

            // The parameter on the arc has to be projected (`cxx:1532-1551`).
            pt2d = GpPnt2d::new(solrst.value(3), solrst.value(4));
            let mut pnt = thearc.d0(ufirst).unwrap_or_default();
            dist = pt2d.distance(&pnt);
            solrst.set_value(1, ufirst);
            pnt = thearc.d0(ulast).unwrap_or_default();
            let distaux = pt2d.distance(&pnt);
            if distaux < dist {
                solrst.set_value(1, ulast);
                dist = distaux;
            }

            if dist > PCONFUSION {
                prm = pmin;
                let mut distaux = 0.0f64;
                let ok = BlendTool::project(&pt2d, thecur.as_ref(), &mut prm, &mut distaux);
                let on_curve = thearc.d0(prm).unwrap_or_default();
                if ok && pt2d.distance(&on_curve) < dist {
                    solrst.set_value(1, prm);
                } else {
                    solrst.set_value(1, pmin);
                }
            }
            // Verify the singular point.
            jalons_trouve = func_inv.is_solution(&solrst, self.tolpoint3d);
        }

        if !jalons_trouve {
            // Resolution.
            match super::solve_root_set(func_inv, &solrst, &solver_toler, &infb, &supb, 35) {
                None => recadre = false,
                Some(root) => {
                    solrst = root;
                    recadre = func_inv.is_solution(&solrst, self.tolpoint3d);
                }
            }
        }

        // On failure, look whether another arc fits the bill: the "exit close to
        // a vertex" case (`cxx:1557-1598`).
        dist = (ulast - ufirst) / 100.0;
        if !recadre && ((pmin - ulast).abs() < dist || (pmin - ufirst).abs() < dist) {
            let new_index =
                self.arc_to_recadre(on_first, sol, index_sol, &mut last_pt2d, &mut pt2d, &mut pmin);
            if new_index == 0 {
                return None;
            }
            let thearc = Self::domain_arc(&tool, new_index)?;
            let thecur = super::arc_curve(&thearc)?;
            solrst.set_value(1, pmin);

            // The problem to solve.
            func_inv.set(on_first, Some(Arc::clone(&thecur)));
            func_inv.get_bounds(&mut infb, &mut supb);
            func_inv.get_tolerance(&mut toler, 0.1 * self.tolpoint3d);
            let retry_toler = toler.clone();
            toler.multiply_scalar(10.0);

            match super::solve_root_set(func_inv, &solrst, &retry_toler, &infb, &supb, 35) {
                None => recadre = false,
                Some(root) => {
                    solrst = root;
                    recadre = func_inv.is_solution(&solrst, self.tolpoint3d);
                }
            }

            if !recadre {
                return None;
            }
            return self.classify_recadre(&tool, new_index, solrst, &toler);
        }

        if !recadre {
            return None;
        }
        self.classify_recadre(&tool, index_sol, solrst, &toler)
    }

    /// The topological classification tail of `Recadre` (`cxx:1560-1600`).
    fn classify_recadre(
        &self,
        tool: &Rc<RefCell<TopolTool>>,
        index_sol: usize,
        solrst: MathVector,
        toler: &MathVector,
    ) -> Option<RecadreSolution> {
        let thearc = Self::domain_arc(tool, index_sol)?;
        let ufirst = thearc.first_parameter();
        let ulast = thearc.last_parameter();

        let (is_vtx, vtx) = {
            let mut iter = tool.borrow_mut();
            iter.initialize_curve(&thearc);
            iter.init_vertex_iterator();
            let mut is_vtx = !iter.more_vertex();
            let mut vtx: Option<HVertex> = None;
            while !is_vtx {
                let Some(candidate) = iter.vertex() else {
                    break;
                };
                // An upper bound of the tolerance (`cxx:1583`).
                let mut vtol = 0.4 * (ulast - ufirst).abs();
                let surf_tol = BlendTool::tolerance(&candidate, &thearc).max(toler.value(1));
                if vtol > surf_tol {
                    vtol = surf_tol;
                }
                if (BlendTool::parameter(&candidate, &thearc) - solrst.value(1)).abs() <= vtol {
                    is_vtx = true; // Inside the vertex ball, or the vertex is
                                   // inside the recadrage ball.
                    vtx = Some(candidate);
                } else {
                    iter.next_vertex();
                    is_vtx = !iter.more_vertex();
                }
            }
            if !iter.more_vertex() {
                is_vtx = false;
            }
            (is_vtx, vtx)
        };

        Some(RecadreSolution {
            solrst,
            index: index_sol,
            is_vtx,
            vtx,
        })
    }
}
