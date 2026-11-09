//! `IntWalk_IWalking::ComputeOpenLine` (`IntWalk_IWalking.gxx:1404-1916`):
//! walk from every start point of `Pnts1` that is not passing and not
//! tangent and not yet processed, until a stop test fires.

use occt_core::gp::GpVec;

use crate::int_tools_wline::PntOn2S;

use super::line_ops::make_walking_point;
use super::types::{math_vec2, uv_on, PathPoint, StatusDeflection};
use super::{ContapIWLineOfTheIWalking, ContapTheIWalking};
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::MathFunctionSetRoot;

impl ContapTheIWalking {
    /// `ComputeOpenLine(Umult, Vmult, Pnts1, Func, Rajout)`
    /// (`gxx:1404-1916`).
    pub(super) fn compute_open_line(
        &mut self,
        umult: &[f64],
        vmult: &[f64],
        pnts1: &[PathPoint],
        func: &mut super::super::surf_function::ContapSurfFunction,
        rajout: &mut bool,
    ) {
        let mut born_inf = [self.um, self.vm];
        let mut born_sup = [self.umax, self.vmax];

        let mut pas_c = 0.0f64;
        let mut arrive;
        let mut cadre;
        let mut arret_ajout;
        let mut psol = PntOn2S {
            p: occt_core::gp::GpPnt::default(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };
        let mut tgtend;
        let mut a_status = StatusDeflection::Ok;
        let mut status_precedent = StatusDeflection::Ok;
        let mut nb_division;
        let mut step_sign;
        let mut n: i32;
        let mut save_n: i32;

        let nb_path = pnts1.len();

        let mut solver = MathFunctionSetRoot::new(func, &self.tolerance_vector(), -1);

        // `gxx:1466-1474`: the movement direction remembered for every point.
        let mut movementdirectioninfo = vec![0i32; nb_path + 1];

        let mut a_func_for_duplicate = func.clone();

        for i in 1..=nb_path {
            if !(self.wd1[i].etat > 11
                || (self.wd1[i].etat < -11 && movementdirectioninfo[i] != 0))
            {
                continue;
            }

            let path_pnt = pnts1[i - 1].clone();
            let mut uvap = [self.wd1[i].ustart, self.wd1[i].vstart];
            make_walking_point(self.reversed, 11, uvap[0], uvap[1], func, &mut self.previous_point);

            if self.is_point_on_line(
                &self.previous_point,
                &math_vec2(born_inf[0], born_inf[1]),
                &math_vec2(born_sup[0], born_sup[1]),
                &mut solver,
                &mut a_func_for_duplicate,
            ) {
                self.wd1[i].etat = -self.wd1[i].etat.abs();
                continue;
            }

            let mut current_line = ContapIWLineOfTheIWalking::new();
            current_line.set_tangency_at_beginning(false);
            tgtend = false;
            current_line.add_status_first_with_point(false, true, i as i32, path_pnt.clone());
            self.previous_d3d = match func.direction3d() {
                Ok(v) => *v,
                Err(_) => GpVec::zero(),
            };
            self.previous_d2d = match func.direction2d() {
                Ok(v) => v,
                Err(_) => self.previous_d2d,
            };
            current_line.add_point(self.previous_point);

            // `gxx:1504-1533`: the sign of the first step.
            if movementdirectioninfo[i] != 0 {
                if movementdirectioninfo[i] < 0 {
                    step_sign = -1;
                    current_line.set_tangent_vector(self.previous_d3d.reversed(), 1);
                } else {
                    step_sign = 1;
                    current_line.set_tangent_vector(self.previous_d3d, 1);
                }
            } else {
                let tyutuyt = match path_pnt.direction_3d() {
                    Some(d) => d.dot(&self.previous_d3d),
                    None => 0.0,
                };
                if tyutuyt < 0.0 {
                    step_sign = -1;
                    current_line.set_tangent_vector(self.previous_d3d.reversed(), 1);
                } else {
                    step_sign = 1;
                    current_line.set_tangent_vector(self.previous_d3d, 1);
                }
            }

            self.wd1[i].etat = -self.wd1[i].etat.abs();
            movementdirectioninfo[i] = if movementdirectioninfo[i] == 0 {
                step_sign
            } else {
                0
            };

            // `gxx:1538-1553`: the first step of advancement.
            let d2dx = self.previous_d2d.x().abs();
            let d2dy = self.previous_d2d.y().abs();
            if d2dx < self.vec_tolerance[0] {
                pas_c = self.pas * (self.vmax - self.vm) / d2dy;
            } else if d2dy < self.vec_tolerance[1] {
                pas_c = self.pas * (self.umax - self.um) / d2dx;
            } else {
                pas_c = self.pas
                    * ((self.umax - self.um) / d2dx).min((self.vmax - self.vm) / d2dy);
            }

            arrive = false;
            arret_ajout = false;
            nb_division = 0;
            status_precedent = StatusDeflection::Ok;
            n = 0;
            save_n = 0;
            let mut index_of_path_point_do_not_check = 0i32;
            let mut a_nb_iter = 10;

            while !arrive {
                cadre = self.cadrage(&mut born_inf, &mut born_sup, &mut uvap, &mut pas_c, step_sign);

                let start = math_vec2(uvap[0], uvap[1]);
                let inf = math_vec2(born_inf[0], born_inf[1]);
                let sup = math_vec2(born_sup[0], born_sup[1]);
                solver.perform_with_bounds(func, &start, &inf, &sup, false);

                if cadre {
                    born_inf = [self.um, self.vm];
                    born_sup = [self.umax, self.vmax];
                }

                if solver.is_done() {
                    if func.root().abs() > func.tolerance() {
                        // `gxx:1584-1601`: no solution for the tolerance.
                        pas_c /= 2.0;
                        let pas_cu = (pas_c * self.previous_d2d.x()).abs();
                        let pas_cv = (pas_c * self.previous_d2d.y()).abs();
                        if pas_cu <= self.vec_tolerance[0] && pas_cv <= self.vec_tolerance[1] {
                            if current_line.nb_points() == 1 {
                                break;
                            }
                            arrive = true;
                            current_line.add_status_last(false);
                            tgtend = true;
                            *rajout = true;
                            let idx = self.lines.len() as i32 + 1;
                            self.seq_alone.push(idx);
                            self.seq_ajout.push(idx);
                        }
                    } else {
                        let root = solver.root();
                        uvap[0] = root.value(1);
                        uvap[1] = root.value(2);

                        let mut passage =
                            self.test_arret_passage_open(umult, vmult, func, &mut uvap);
                        arrive = passage.is_some();
                        if arrive {
                            n = passage.take().unwrap_or(0);
                            cadre = false;
                        } else {
                            if *rajout {
                                let a = self.test_arret_ajout(func, &mut uvap, &mut psol);
                                arret_ajout = a.is_some();
                                save_n = a.unwrap_or(0);
                                n = save_n;
                                if arret_ajout {
                                    tgtend = self
                                        .lines
                                        .get((n.unsigned_abs() as usize).saturating_sub(1))
                                        .map(|l| l.is_tangent_at_end())
                                        .unwrap_or(false);
                                    n = -n;
                                }
                            }

                            if !(*rajout && arret_ajout) {
                                let (prev_up, prev_vp) = uv_on(&self.previous_point, self.reversed);
                                let found = super::stop_tests::test_passed_solution_with_negative_state(
                                    &self.wd1,
                                    umult,
                                    vmult,
                                    prev_up,
                                    prev_vp,
                                    &self.nb_multiplicities,
                                    self.vec_tolerance,
                                    func,
                                    &mut uvap,
                                );
                                if let Some(idx) = found {
                                    n = idx;
                                    arrive = true;
                                    cadre = false;
                                }
                            }

                            if !arret_ajout && cadre {
                                if current_line.nb_points() == 1 {
                                    break;
                                }
                                if let Some(idx) = self.test_arret_cadre(
                                    umult,
                                    vmult,
                                    &mut current_line,
                                    func,
                                    &mut uvap,
                                ) {
                                    n = idx;
                                } else {
                                    n = 0;
                                }
                                if n <= 0 {
                                    make_walking_point(self.reversed, 2, uvap[0], uvap[1], func, &mut psol);
                                    tgtend = func.is_tangent();
                                    n = -n;
                                }
                            }
                        }

                        a_status = self.test_deflection(
                            func,
                            arrive,
                            &uvap,
                            status_precedent,
                            &mut nb_division,
                            &mut pas_c,
                            step_sign,
                        );
                        status_precedent = a_status;

                        if a_status == StatusDeflection::PasTropGrand {
                            arrive = false;
                            arret_ajout = false;
                            tgtend = false;
                            let (pu, pv) = uv_on(&self.previous_point, self.reversed);
                            uvap[0] = pu;
                            uvap[1] = pv;
                        } else if arret_ajout || cadre {
                            arrive = true;
                            current_line.add_status_last(false);
                            current_line.add_point(psol);

                            // `gxx:1688-1697`: drop `SaveN` from `seqAlone`.
                            if let Some(pos) = self.seq_alone.iter().position(|&v| v == save_n) {
                                self.seq_alone.remove(pos);
                            }

                            if cadre && n == 0 {
                                *rajout = true;
                                let idx = self.lines.len() as i32 + 1;
                                self.seq_ajout.push(idx);
                            }
                        } else if a_status == StatusDeflection::ArretSurPointPrecedent {
                            if current_line.nb_points() == 1 {
                                arrive = false;
                                break;
                            }
                            arrive = true;
                            *rajout = true;
                            let idx = self.lines.len() as i32 + 1;
                            self.seq_alone.push(idx);
                            self.seq_ajout.push(idx);
                            current_line.add_status_last(false);
                            tgtend = true;
                        } else if arrive {
                            if current_line.nb_points() == 1
                                && (n == i as i32 || a_status == StatusDeflection::PointConfondu)
                            {
                                arrive = false;
                                break;
                            }

                            // `gxx:1729-1776`: the stop point is a start point.
                            let path_pnt_n = pnts1[(n.max(1) as usize) - 1].clone();
                            let etat1n = self.wd1[n.max(1) as usize].etat;

                            if etat1n.abs() < 11 {
                                if a_status == StatusDeflection::ArretSurPoint {
                                    current_line.add_status_last(false);
                                    tgtend = true;
                                } else {
                                    arrive = false;
                                }
                                current_line.add_index_passing(n);
                            } else {
                                if etat1n == 11 {
                                    tgtend = true;
                                }
                                current_line.add_status_last_with_point(
                                    true,
                                    n,
                                    path_pnt_n.clone(),
                                );
                            }

                            self.add_point_in_current_line(
                                n.max(1),
                                &path_pnt_n,
                                &mut current_line,
                            );

                            if etat1n != 1 && etat1n != 11 {
                                let idx1 = n.max(1) as usize;
                                self.wd1[idx1].etat = -etat1n.abs();
                                movementdirectioninfo[idx1] = if movementdirectioninfo[idx1] == 0 {
                                    step_sign
                                } else {
                                    0
                                };
                                if arrive && movementdirectioninfo[idx1] != 0 {
                                    index_of_path_point_do_not_check = n;
                                }
                                if arrive {
                                    *rajout = true;
                                    let idx = self.lines.len() as i32 + 1;
                                    self.seq_ajout.push(idx);
                                }
                            }
                        } else if a_status == StatusDeflection::ArretSurPoint {
                            arrive = true;
                            current_line.add_status_last(false);
                            tgtend = true;
                            make_walking_point(self.reversed, 1, uvap[0], uvap[1], func, &mut psol);
                            current_line.add_point(psol);
                            *rajout = true;
                            let idx = self.lines.len() as i32 + 1;
                            self.seq_alone.push(idx);
                            self.seq_ajout.push(idx);
                        } else if a_status == StatusDeflection::Ok {
                            make_walking_point(self.reversed, 2, uvap[0], uvap[1], func, &mut self.previous_point);
                            if let Ok(v) = func.direction3d() {
                                self.previous_d3d = *v;
                            }
                            if let Ok(v) = func.direction2d() {
                                self.previous_d2d = v;
                            }
                            current_line.add_point(self.previous_point);
                        } else if a_status == StatusDeflection::PointConfondu {
                            a_nb_iter -= 1;
                        }
                    }
                } else {
                    // `gxx:1801-1816`: no numerical solution.
                    pas_c /= 2.0;
                    let pas_cu = (pas_c * self.previous_d2d.x()).abs();
                    let pas_cv = (pas_c * self.previous_d2d.y()).abs();
                    if pas_cu <= self.vec_tolerance[0] && pas_cv <= self.vec_tolerance[1] {
                        if current_line.nb_points() == 1 {
                            break;
                        }
                        arrive = true;
                        current_line.add_status_last(false);
                        tgtend = true;
                        *rajout = true;
                        let idx = self.lines.len() as i32 + 1;
                        self.seq_alone.push(idx);
                        self.seq_ajout.push(idx);
                    }
                }

                if a_nb_iter < 0 {
                    break;
                }
            }

            if arrive {
                current_line.set_tangency_at_end(tgtend);
                self.lines.push(current_line);
                let last = self.lines.len() - 1;

                movementdirectioninfo[i] = 0;
                if self.wd1[i].etat > 0 {
                    self.wd1[i].etat = -self.wd1[i].etat;
                }

                // `gxx:1836-1910`: mark every start point that coincides with
                // the end or the beginning of the new line.
                for av in 1..=nb_path {
                    if !(self.wd1[av].etat > 11
                        || (av != i
                            && av != index_of_path_point_do_not_check as usize
                            && self.wd1[av].etat < -11
                            && movementdirectioninfo[av] != 0))
                    {
                        continue;
                    }

                    let mut uav = self.wd1[av].ustart;
                    let mut vav = self.wd1[av].vstart;

                    let last_pt = self.lines[last].nb_points();
                    let (uavp, vavp) = uv_on(self.lines[last].value(last_pt), self.reversed);
                    uav -= uavp;
                    vav -= vavp;
                    uav *= 0.001;
                    vav *= 0.001;
                    if uav.abs() < self.vec_tolerance[0] && vav.abs() < self.vec_tolerance[1] {
                        if self.wd1[av].etat < 0 {
                            movementdirectioninfo[av] = 0;
                        } else {
                            self.wd1[av].etat = -self.wd1[av].etat;
                            movementdirectioninfo[av] = step_sign;
                        }
                        self.lines[last].add_status_last_with_point(
                            true,
                            av as i32,
                            pnts1[av - 1].clone(),
                        );
                    }

                    let (uavp, vavp) = uv_on(self.lines[last].value(1), self.reversed);
                    let mut uav = self.wd1[av].ustart;
                    let mut vav = self.wd1[av].vstart;
                    uav -= uavp;
                    vav -= vavp;
                    uav *= 0.001;
                    vav *= 0.001;
                    if uav.abs() < self.vec_tolerance[0] && vav.abs() < self.vec_tolerance[1] {
                        if self.wd1[av].etat < 0 {
                            movementdirectioninfo[av] = 0;
                        } else {
                            self.wd1[av].etat = -self.wd1[av].etat;
                            movementdirectioninfo[av] = -step_sign;
                        }
                        self.lines[last].add_status_first_with_point(
                            false,
                            true,
                            av as i32,
                            pnts1[av - 1].clone(),
                        );
                    }
                }
            }
        }
    }
}
