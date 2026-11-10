//! `IntWalk_IWalking::ComputeCloseLine` (`IntWalk_IWalking.gxx:1993-2605`)
//! and `FillPntsInHoles` (`gxx:1023-1174`).

use occt_core::gp::{GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::{PCONFUSION, RESOLUTION};
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::MathFunctionSetRoot;

use occt_topo::int_tools_wline::PntOn2S;

use super::line_ops::make_walking_point;
use super::types::{math_vec2, uv_on, InteriorPoint, PathPoint, StatusDeflection};
use super::{ContapIWLineOfTheIWalking, ContapTheIWalking};

impl ContapTheIWalking {
    /// `ComputeCloseLine(Umult, Vmult, Pnts1, Pnts2, Func, Rajout)`
    /// (`gxx:1993-2605`): walk from every interior point of `Pnts2`, closing
    /// the line or opening it on a border.
    pub(super) fn compute_close_line(
        &mut self,
        umult: &[f64],
        vmult: &[f64],
        pnts1: &[PathPoint],
        pnts2: &[InteriorPoint],
        func: &mut super::super::surf_function::ContapSurfFunction,
        rajout: &mut bool,
    ) {
        let mut born_inf = [self.um, self.vm];
        let mut born_sup = [self.umax, self.vmax];

        let mut pas_c;
        let mut pas_sav;
        let mut arrive;
        let mut cadre;
        let mut arret_ajout;
        let mut mut_psol = PntOn2S {
            p: occt_core::gp::GpPnt::default(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };
        let mut tgtbeg;
        let mut tgtend;
        let mut step_sign;
        let mut a_status = StatusDeflection::Ok;
        let mut status_precedent;
        let mut nb_division;
        let mut n: i32;
        let mut save_n: i32;
        let mut ipass: i32;

        let nb_loop = pnts2.len();

        let mut solver = MathFunctionSetRoot::new(func, &self.tolerance_vector(), -1);

        // `gxx:2055-2094`: borders whose derivative along the other parameter
        // vanishes everywhere are degenerated.
        let mut is_left_degenerated_border = [true, true];
        let mut is_right_degenerated_border = [true, true];
        {
            let a_nb_sample_pnt = 10;
            let mut a_step = [0.0f64; 2];
            for k in 0..2 {
                a_step[k] = (born_sup[k] - born_inf[k]) / f64::from(a_nb_sample_pnt - 1);
            }
            let mut d = occt_math::MathMatrix::new(1, 1, 1, 2);
            for border_idx in 1..=2usize {
                let a_change_idx = if border_idx == 2 { 1 } else { 2 };
                let mut uv = [0.0f64; 2];

                uv[border_idx - 1] = born_inf[border_idx - 1];
                for param_idx in 0..a_nb_sample_pnt {
                    let a_param = born_inf[a_change_idx - 1] + f64::from(param_idx) * a_step[a_change_idx - 1];
                    uv[a_change_idx - 1] = a_param;
                    let _ = func.derivatives(&math_vec2(uv[0], uv[1]), &mut d);
                    if d.value(1, a_change_idx).abs() > occt_core::precision::CONFUSION {
                        is_left_degenerated_border[border_idx - 1] = false;
                        break;
                    }
                }

                uv[border_idx - 1] = born_sup[border_idx - 1];
                for param_idx in 0..a_nb_sample_pnt {
                    let a_param = born_inf[a_change_idx - 1] + f64::from(param_idx) * a_step[a_change_idx - 1];
                    uv[a_change_idx - 1] = a_param;
                    let _ = func.derivatives(&math_vec2(uv[0], uv[1]), &mut d);
                    if d.value(1, a_change_idx).abs() > occt_core::precision::CONFUSION {
                        is_right_degenerated_border[border_idx - 1] = false;
                        break;
                    }
                }
            }
        }

        let mut a_func_for_duplicate = func.clone();

        for i in 1..=nb_loop {
            if self.wd2[i].etat <= 12 {
                continue;
            }
            let loop_pnt = pnts2[i - 1];

            self.previous_point.set_point(*loop_pnt.value());
            self.previous_point
                .set_uv_on(!self.reversed, self.wd2[i].ustart, self.wd2[i].vstart);

            if self.is_point_on_line(
                &self.previous_point,
                &math_vec2(born_inf[0], born_inf[1]),
                &math_vec2(born_sup[0], born_sup[1]),
                &mut solver,
                &mut a_func_for_duplicate,
            ) {
                self.wd2[i].etat = -self.wd2[i].etat;
                continue;
            }

            self.previous_d3d = *loop_pnt.direction();
            self.previous_d2d = match occt_core::gp::GpDir2d::new(
                loop_pnt.direction_2d().x(),
                loop_pnt.direction_2d().y(),
            ) {
                Ok(d) => d,
                Err(_) => self.previous_d2d,
            };

            let mut current_line = ContapIWLineOfTheIWalking::new();
            current_line.add_point(self.previous_point);
            current_line.set_tangent_vector(self.previous_d3d, 1);
            tgtbeg = false;
            tgtend = false;
            let mut uvap = [self.wd2[i].ustart, self.wd2[i].vstart];

            step_sign = 1;

            // `gxx:2126-2141`: first step of advancement.
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

            pas_sav = pas_c;

            arrive = false;
            arret_ajout = false;
            nb_division = 0;
            status_precedent = StatusDeflection::Ok;
            n = 0;
            save_n = 0;
            // `gxx:2013` initializes `Ipass = 0`, but `gxx:2229`
            // (`TestArretPassage`) always writes it inside this loop before the
            // `gxx:2470` read, so the initial store is dead in OCCT as well.

            let mut a_nb_iter = 10;

            while !arrive {
                cadre = self.cadrage(&mut born_inf, &mut born_sup, &mut uvap, &mut pas_c, step_sign);

                let start = math_vec2(uvap[0], uvap[1]);
                let inf = math_vec2(born_inf[0], born_inf[1]);
                let sup = math_vec2(born_sup[0], born_sup[1]);
                solver.perform_with_bounds(func, &start, &inf, &sup, false);

                let mut is_on_degenerated_border = false;

                if cadre {
                    born_inf = [self.um, self.vm];
                    born_sup = [self.umax, self.vmax];
                }

                if solver.is_done() {
                    if func.root().abs() > func.tolerance() {
                        // `gxx:2174-2219`.
                        pas_c /= 2.0;
                        let pas_cu = (pas_c * self.previous_d2d.x()).abs();
                        let pas_cv = (pas_c * self.previous_d2d.y()).abs();
                        if pas_cu <= self.vec_tolerance[0] && pas_cv <= self.vec_tolerance[1] {
                            if current_line.nb_points() == 1 {
                                self.remove_two_end_points(i as i32);
                                break;
                            }
                            if self.wd2[i].etat > 12 {
                                self.wd2[i].etat = 12;
                                arret_ajout = false;
                                let idx = -((self.lines.len() + 1) as i32);
                                let psol = mut_psol;
                                self.open_line(idx, &psol, pnts1, func, &mut current_line);
                                step_sign = -1;
                                status_precedent = StatusDeflection::Ok;
                                arrive = false;
                                pas_c = pas_sav;
                                *rajout = true;
                                self.seq_alone.push(idx);
                                self.seq_ajout.push(idx);
                            } else {
                                arrive = true;
                                current_line.add_status_last(false);
                                *rajout = true;
                                let idx = self.lines.len() as i32 + 1;
                                self.seq_alone.push(idx);
                                self.seq_ajout.push(idx);
                                tgtend = true;
                            }
                        }
                    } else {
                        let root = solver.root();
                        uvap[0] = root.value(1);
                        uvap[1] = root.value(2);

                        // `gxx:2226-2264`: fix a point stuck on a degenerated
                        // border by linear extrapolation.
                        if current_line.nb_points() > 2 {
                            for a_coord_idx in 1..=2usize {
                                let on_left = is_left_degenerated_border[a_coord_idx - 1]
                                    && (uvap[a_coord_idx - 1] - born_inf[a_coord_idx - 1]).abs()
                                        < PCONFUSION;
                                let on_right = is_right_degenerated_border[a_coord_idx - 1]
                                    && (uvap[a_coord_idx - 1] - born_sup[a_coord_idx - 1]).abs()
                                        < PCONFUSION;
                                if !(on_left || on_right) {
                                    continue;
                                }

                                let nbp = current_line.nb_points();
                                let uvprev = uv_on(current_line.value(nbp - 1), self.reversed);
                                let uv = uv_on(current_line.value(nbp), self.reversed);
                                let uvprev = [uvprev.0, uvprev.1];
                                let uv = [uv.0, uv.1];

                                let mut a_scale_coeff = 0.0;
                                if a_status != StatusDeflection::PasTropGrand {
                                    if (uv[a_coord_idx - 1] - uvprev[a_coord_idx - 1]).abs()
                                        > RESOLUTION
                                    {
                                        a_scale_coeff = ((uvap[a_coord_idx - 1]
                                            - uv[a_coord_idx - 1])
                                            / (uv[a_coord_idx - 1] - uvprev[a_coord_idx - 1]))
                                            .abs();
                                    }
                                    let a_fix_idx = if a_coord_idx == 1 { 2 } else { 1 };
                                    uvap[a_fix_idx - 1] = uv[a_fix_idx - 1]
                                        + (uv[a_fix_idx - 1] - uvprev[a_fix_idx - 1]) * a_scale_coeff;
                                    is_on_degenerated_border = true;
                                }
                            }
                        }

                        let (passed, pass_idx) =
                            self.test_arret_passage_closed(umult, vmult, &uvap, i as i32);
                        arrive = passed;
                        ipass = pass_idx;
                        if arrive {
                            // `gxx:2270-2283`: reset the parameters to the
                            // first point of the line before the arrow test.
                            let first = *current_line.value(1);
                            let (u, v) = uv_on(&first, self.reversed);
                            uvap[0] = u;
                            uvap[1] = v;
                            cadre = false;
                        } else {
                            if *rajout {
                                let a = self.test_arret_ajout(func, &mut uvap, &mut mut_psol);
                                arret_ajout = a.is_some();
                                save_n = a.unwrap_or(0);
                                n = save_n;
                                if arret_ajout {
                                    if n > 0 {
                                        tgtend = self
                                            .lines
                                            .get((n.unsigned_abs() as usize).saturating_sub(1))
                                            .map(|l| l.is_tangent_at_end())
                                            .unwrap_or(false);
                                        n = -n;
                                    } else {
                                        tgtend = self
                                            .lines
                                            .get((n.unsigned_abs() as usize).saturating_sub(1))
                                            .map(|l| l.is_tangent_at_beginning())
                                            .unwrap_or(false);
                                    }
                                    arrive = self.wd2[i].etat == 12;
                                }
                            }

                            if !arret_ajout && cadre {
                                if current_line.nb_points() == 1 {
                                    self.remove_two_end_points(i as i32);
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
                                save_n = n;
                                if n <= 0 {
                                    make_walking_point(
                                        self.reversed,
                                        2,
                                        uvap[0],
                                        uvap[1],
                                        func,
                                        &mut mut_psol,
                                    );
                                    tgtend = func.is_tangent();
                                    if is_on_degenerated_border {
                                        tgtend = true;
                                    }
                                    n = -n;
                                }
                                arrive = self.wd2[i].etat == 12;
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

                        if is_on_degenerated_border && tgtend {
                            a_status = StatusDeflection::ArretSurPoint;
                        }
                        status_precedent = a_status;

                        if a_status == StatusDeflection::PasTropGrand {
                            arrive = false;
                            arret_ajout = false;
                            tgtend = false;
                            let (pu, pv) = uv_on(&self.previous_point, self.reversed);
                            uvap[0] = pu;
                            uvap[1] = pv;
                        } else if arret_ajout || cadre {
                            if arrive {
                                // `gxx:2350-2392`: the line becomes open.
                                current_line.add_status_last(false);
                                current_line.add_point(mut_psol);

                                if self.is_valid_end_point(i as i32, save_n) {
                                    if let Some(pos) =
                                        self.seq_alone.iter().position(|&v| v == save_n)
                                    {
                                        self.seq_alone.remove(pos);
                                    }
                                    if current_line.nb_points() <= 3 {
                                        if let Some(pos) =
                                            self.seq_ajout.iter().position(|&v| v == save_n)
                                        {
                                            self.seq_ajout.remove(pos);
                                        }
                                    }
                                } else {
                                    let target = -((self.lines.len() + 1) as i32);
                                    if self.seq_alone.last() == Some(&target) {
                                        self.seq_alone.pop();
                                        self.seq_ajout.pop();
                                    }
                                    self.remove_two_end_points(i as i32);
                                    arrive = false;
                                    break;
                                }

                                if cadre && n == 0 {
                                    *rajout = true;
                                    let idx = self.lines.len() as i32 + 1;
                                    self.seq_ajout.push(idx);
                                }
                            } else {
                                self.wd2[i].etat = 12;
                                tgtbeg = tgtend;
                                tgtend = false;
                                arret_ajout = false;
                                step_sign = -1;
                                status_precedent = StatusDeflection::Ok;
                                pas_c = pas_sav;

                                if *rajout || (!arret_ajout && cadre && save_n <= 0) {
                                    if a_status == StatusDeflection::ArretSurPointPrecedent {
                                        current_line.add_point(mut_psol);
                                        let psol = mut_psol;
                                        self.open_line(0, &psol, pnts1, func, &mut current_line);
                                    } else {
                                        let idx = -((self.lines.len() + 1) as i32);
                                        let psol = mut_psol;
                                        self.open_line(idx, &psol, pnts1, func, &mut current_line);
                                    }
                                }

                                if self.is_valid_end_point(i as i32, save_n) {
                                    if let Some(pos) =
                                        self.seq_alone.iter().position(|&v| v == save_n)
                                    {
                                        self.seq_alone.remove(pos);
                                    }
                                    if current_line.nb_points() <= 2 {
                                        if let Some(pos) =
                                            self.seq_ajout.iter().position(|&v| v == save_n)
                                        {
                                            self.seq_ajout.remove(pos);
                                        }
                                    }
                                } else {
                                    self.remove_two_end_points(i as i32);
                                    break;
                                }

                                if cadre && n == 0 {
                                    *rajout = true;
                                    let idx = -((self.lines.len() + 1) as i32);
                                    self.seq_ajout.push(idx);
                                }
                            }
                        } else if a_status == StatusDeflection::ArretSurPointPrecedent {
                            if current_line.nb_points() == 1 {
                                arrive = false;
                                self.remove_two_end_points(i as i32);
                                break;
                            }
                            if self.wd2[i].etat > 12 {
                                self.wd2[i].etat = 12;
                                arret_ajout = false;
                                let psol = mut_psol;
                                self.open_line(0, &psol, pnts1, func, &mut current_line);
                                step_sign = -1;
                                status_precedent = StatusDeflection::Ok;
                                arrive = false;
                                pas_c = pas_sav;
                                *rajout = true;
                                let idx = -((self.lines.len() + 1) as i32);
                                self.seq_alone.push(idx);
                                self.seq_ajout.push(idx);
                            } else {
                                arrive = true;
                                current_line.add_status_last(false);
                                *rajout = true;
                                let idx = self.lines.len() as i32 + 1;
                                self.seq_alone.push(idx);
                                self.seq_ajout.push(idx);
                            }
                        } else if arrive {
                            if self.wd2[i].etat > 12 {
                                current_line.add_status_first_last(true, false, false);
                                let first = *current_line.value(1);
                                current_line.add_point(first);
                            } else if n > 0 && pnts1.len() >= n as usize {
                                let path_pnt = pnts1[(n - 1) as usize].clone();
                                current_line.add_status_last_with_point(
                                    true,
                                    n,
                                    path_pnt.clone(),
                                );
                                self.add_point_in_current_line(
                                    n,
                                    &path_pnt,
                                    &mut current_line,
                                );
                            }
                        } else if a_status == StatusDeflection::ArretSurPoint {
                            if self.wd2[i].etat > 12 {
                                self.wd2[i].etat = 12;
                                tgtbeg = true;
                                tgtend = false;
                                n = -((self.lines.len() + 1) as i32);
                                let mut psol = mut_psol;
                                psol.set_point(*func.point());
                                psol.set_uv_on(!self.reversed, uvap[0], uvap[1]);
                                self.open_line(n, &psol, pnts1, func, &mut current_line);
                                step_sign = -1;
                                *rajout = true;
                                self.seq_alone.push(n);
                                self.seq_ajout.push(n);
                                status_precedent = StatusDeflection::Ok;
                                arrive = false;
                                pas_c = pas_sav;
                            } else {
                                arrive = true;
                                if ipass != 0 {
                                    let path_pnt = pnts1[(ipass - 1) as usize].clone();
                                    current_line.add_status_last_with_point(
                                        true,
                                        ipass,
                                        path_pnt.clone(),
                                    );
                                    self.add_point_in_current_line(
                                        ipass,
                                        &path_pnt,
                                        &mut current_line,
                                    );
                                } else {
                                    current_line.add_status_last(false);
                                    let mut new_p = mut_psol;
                                    new_p.set_point(*func.point());
                                    new_p.set_uv_on(!self.reversed, uvap[0], uvap[1]);
                                    current_line.add_point(new_p);
                                    *rajout = true;
                                    let idx = self.lines.len() as i32 + 1;
                                    self.seq_alone.push(idx);
                                    self.seq_ajout.push(idx);
                                }
                            }
                        } else if a_status == StatusDeflection::Ok {
                            if ipass != 0 {
                                current_line.add_index_passing(ipass);
                            }
                            self.previous_point.set_point(*func.point());
                            self.previous_point
                                .set_uv_on(!self.reversed, uvap[0], uvap[1]);
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
                    // `gxx:2544-2587`: no numerical solution.
                    pas_c /= 2.0;
                    let pas_cu = (pas_c * self.previous_d2d.x()).abs();
                    let pas_cv = (pas_c * self.previous_d2d.y()).abs();
                    if pas_cu <= self.vec_tolerance[0] && pas_cv <= self.vec_tolerance[1] {
                        if current_line.nb_points() == 1 {
                            self.remove_two_end_points(i as i32);
                            break;
                        }
                        if self.wd2[i].etat > 12 {
                            self.wd2[i].etat = 12;
                            arret_ajout = false;
                            let psol = mut_psol;
                            self.open_line(0, &psol, pnts1, func, &mut current_line);
                            step_sign = -1;
                            status_precedent = StatusDeflection::Ok;
                            arrive = false;
                            pas_c = pas_sav;
                            *rajout = true;
                            let idx = -((self.lines.len() + 1) as i32);
                            self.seq_alone.push(idx);
                            self.seq_ajout.push(idx);
                        } else {
                            arrive = true;
                            current_line.add_status_last(false);
                            tgtend = true;
                            *rajout = true;
                            let idx = self.lines.len() as i32 + 1;
                            self.seq_alone.push(idx);
                            self.seq_ajout.push(idx);
                        }
                    }
                }

                if a_nb_iter < 0 {
                    break;
                }
            }

            if arrive {
                current_line.set_tangency_at_beginning(tgtbeg);
                current_line.set_tangency_at_end(tgtend);
                self.lines.push(current_line);
                self.wd2[i].etat = -self.wd2[i].etat;
            }
        }
    }

    /// `FillPntsInHoles(sp, CopySeqAlone, PntsInHoles)` (`gxx:1023-1174`):
    /// try to close the free ends of `seqAlone` with an interior point.
    pub(super) fn fill_pnts_in_holes(
        &mut self,
        sp: &mut super::super::surf_function::ContapSurfFunction,
        copy_seq_alone: &mut [i32],
        pnts_in_holes: &mut Vec<InteriorPoint>,
    ) {
        let born_inf = math_vec2(self.um, self.vm);
        let born_sup = math_vec2(self.umax, self.vmax);
        self.point_line_line.clear();
        let mut seq_to_remove: Vec<i32> = Vec::new();
        let mut bad_solutions: Vec<i32> = Vec::new();

        let mut i = 1usize;
        while i < copy_seq_alone.len() {
            let irang1 = copy_seq_alone[i - 1];
            if irang1 == 0 {
                i += 1;
                continue;
            }
            let mut to_remove = false;
            let line1 = match self
                .lines
                .get((irang1.unsigned_abs() as usize).saturating_sub(1))
            {
                Some(l) => l,
                None => {
                    i += 1;
                    continue;
                }
            };
            let point_alone1 = if irang1 > 0 {
                *line1.value(line1.nb_points())
            } else {
                *line1.value(1)
            };
            let p2d1: GpPnt2d = {
                let (u, v) = uv_on(&point_alone1, self.reversed);
                GpPnt2d::new(u, v)
            };

            let mut min_sq_dist = f64::MAX;
            let mut min_rang = 0i32;
            let mut min_index = 0usize;

            for j in (i + 1)..=copy_seq_alone.len() {
                let irang2 = copy_seq_alone[j - 1];
                if irang2 == 0 || bad_solutions.contains(&irang2) {
                    continue;
                }
                let line2 = match self
                    .lines
                    .get((irang2.unsigned_abs() as usize).saturating_sub(1))
                {
                    Some(l) => l,
                    None => continue,
                };
                let point_alone2 = if irang2 > 0 {
                    *line2.value(line2.nb_points())
                } else {
                    *line2.value(1)
                };
                let (u, v) = uv_on(&point_alone2, self.reversed);
                let p2d2 = GpPnt2d::new(u, v);
                let a_sq_dist = p2d1.square_distance(&p2d2);
                if a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                    min_rang = irang2;
                    min_index = j;
                }
            }

            if min_rang == 0 {
                seq_to_remove.push(irang1);
                bad_solutions.clear();
                i += 1;
                continue;
            }

            if irang1.abs() == min_rang.abs()
                && self.lines[(irang1.unsigned_abs() as usize) - 1].nb_points() == 2
            {
                seq_to_remove.push(irang1);
                seq_to_remove.push(min_rang);
                copy_seq_alone[i - 1] = 0;
                copy_seq_alone[min_index - 1] = 0;
                bad_solutions.clear();
                i += 1;
                continue;
            }

            let line2 = &self.lines[(min_rang.unsigned_abs() as usize) - 1];
            let point_alone2 = if min_rang > 0 {
                *line2.value(line2.nb_points())
            } else {
                *line2.value(1)
            };
            let pnt1 = point_alone1.p;
            let pnt2 = point_alone2.p;
            let (u2, v2) = uv_on(&point_alone2, self.reversed);
            let p2d2 = GpPnt2d::new(u2, v2);
            let min_sq_dist3d = pnt1.square_distance(&pnt2);

            if min_sq_dist3d <= self.epsilon
                || ((p2d1.coord.x - p2d2.coord.x).abs() <= self.vec_tolerance[0]
                    && (p2d1.coord.y - p2d2.coord.y).abs() <= self.vec_tolerance[1])
            {
                to_remove = true;
            } else {
                let mut a_uvap = math_vec2(
                    (p2d1.coord.x + p2d2.coord.x) / 2.0,
                    (p2d1.coord.y + p2d2.coord.y) / 2.0,
                );
                let mut rsnld = MathFunctionSetRoot::new(sp, &self.tolerance_vector(), -1);
                rsnld.perform_with_bounds(sp, &a_uvap, &born_inf, &born_sup, false);

                if rsnld.is_done()
                    && sp.root().abs() <= sp.tolerance()
                    && !sp.is_tangent()
                {
                    let root = rsnld.root();
                    a_uvap = root;
                    let pmid = GpPnt2d::new(a_uvap.value(1), a_uvap.value(2));
                    let p1p2 = GpVec2d::new(
                        p2d2.coord.x - p2d1.coord.x,
                        p2d2.coord.y - p2d1.coord.y,
                    );
                    let p1_pmid = GpVec2d::new(
                        pmid.coord.x - p2d1.coord.x,
                        pmid.coord.y - p2d1.coord.y,
                    );
                    let p2_pmid = GpVec2d::new(
                        pmid.coord.x - p2d2.coord.x,
                        pmid.coord.y - p2d2.coord.y,
                    );
                    let scal_prod1 = p1p2.dot(&p1_pmid);
                    let scal_prod2 = p1p2.dot(&p2_pmid);
                    let mut is_pmid_valid = scal_prod1 > 0.0 && scal_prod2 < 0.0;
                    if is_pmid_valid {
                        for iline in 1..=self.lines.len() {
                            if self.is_point_on_line_2d(&pmid, iline as i32) {
                                is_pmid_valid = false;
                                break;
                            }
                        }
                    }

                    if is_pmid_valid {
                        let a_point = InteriorPoint::new(
                            *sp.point(),
                            a_uvap.value(1),
                            a_uvap.value(2),
                            sp.direction3d().map(|v| *v).unwrap_or_else(|_| GpVec::zero()),
                            GpVec2d::new(
                                sp.direction2d().map(|d| d.x()).unwrap_or(0.0),
                                sp.direction2d().map(|d| d.y()).unwrap_or(0.0),
                            ),
                        );
                        pnts_in_holes.push(a_point);
                        self.point_line_line
                            .insert(pnts_in_holes.len() as i32, vec![irang1, min_rang]);
                    } else {
                        bad_solutions.push(min_rang);
                        continue;
                    }
                } else {
                    bad_solutions.push(min_rang);
                    continue;
                }
            }

            copy_seq_alone[i - 1] = 0;
            copy_seq_alone[min_index - 1] = 0;
            if to_remove {
                seq_to_remove.push(irang1);
                seq_to_remove.push(min_rang);
            }
            bad_solutions.clear();
            i += 1;
        }

        for &to_rm in seq_to_remove.iter() {
            if let Some(pos) = self.seq_alone.iter().position(|&v| v == to_rm) {
                self.seq_alone.remove(pos);
            }
        }
    }
}
