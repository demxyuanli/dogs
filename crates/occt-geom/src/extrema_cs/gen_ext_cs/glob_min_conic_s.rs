//! `GlobMinConicS` — the 2-variable PSO seed for a conic against a
//! non-quadric surface (`Extrema_GenExtCS.cxx:522-708`). Split out of
//! `gen_ext_cs/mod.rs` to keep every file under the 800-line gate; the body is
//! unchanged.

use super::*;

impl<'a> GenExtCS<'a> {
    /// `GlobMinConicS` (`cxx:522-708`).
    pub(super) fn glob_min_conic_s(
        &mut self,
        the_c: &'a dyn Curve,
        the_nb_particles: usize,
        the_tuvinf: &MathVector,
        the_tuvsup: &MathVector,
        the_tuv: &mut MathVector,
    ) {
        let s = self.s.expect("Extrema_GenExtCS: no surface");
        let a_nb_var = 2;
        let mut an_uvinf = MathVector::new(1, a_nb_var);
        let mut an_uvsup = MathVector::new(1, a_nb_var);
        for i in 1..=a_nb_var {
            an_uvinf.set_value(i, the_tuvinf.value(i + 1));
            an_uvsup.set_value(i, the_tuvsup.value(i + 1));
        }

        let mut a_particles = PsoParticlesPool::new(the_nb_particles, a_nb_var);

        let mut a_min_uv = MathVector::new(1, a_nb_var);
        for i in 1..=a_nb_var {
            a_min_uv.set_value(
                i,
                an_uvinf.value(i) + (an_uvsup.value(i) - an_uvinf.value(i)) / A_BORDER_DIVISOR,
            );
        }
        let mut a_max_uv = MathVector::new(1, a_nb_var);
        for i in 1..=a_nb_var {
            a_max_uv.set_value(
                i,
                an_uvsup.value(i) - (an_uvsup.value(i) - an_uvinf.value(i)) / A_BORDER_DIVISOR,
            );
        }

        // Increase numbers of UV samples to improve searching global minimum —
        // cxx:547-550.
        let an_addsample = (self.tsample / 2).max(3);
        let an_usample = self.usample + an_addsample;
        let a_vsample = self.vsample + an_addsample;

        let a_step_s_u = (a_max_uv.value(1) - a_min_uv.value(1)) / an_usample as f64;
        let a_step_s_v = (a_max_uv.value(2) - a_min_uv.value(2)) / a_vsample as f64;

        let mut a_func = GlobOptFuncConicS::new(
            s,
            an_uvinf.value(1),
            an_uvsup.value(1),
            an_uvinf.value(2),
            an_uvsup.value(2),
        );
        a_func.load_conic(the_c, the_tuvinf.value(1), the_tuvsup.value(1));

        let mut an_uv = MathVector::new(1, a_nb_var);
        let mut a_particle = a_particles.worst_particle_index();
        let mut a_su = a_min_uv.value(1);
        for _ in 0..=an_usample {
            an_uv.set_value(1, a_su);
            let mut a_sv = a_min_uv.value(2);
            for _ in 0..=a_vsample {
                an_uv.set_value(2, a_sv);
                let mut a_sq_dist = 0.0;
                if !a_func.value(&an_uv, &mut a_sq_dist) {
                    a_sq_dist = occt_core::precision::INFINITE;
                }
                if a_sq_dist < a_particles.particle(a_particle).distance {
                    let p = a_particles.particle_mut(a_particle);
                    p.position[0] = a_su;
                    p.position[1] = a_sv;
                    p.best_position[0] = a_su;
                    p.best_position[1] = a_sv;
                    p.distance = a_sq_dist;
                    p.best_distance = a_sq_dist;
                    a_particle = a_particles.worst_particle_index();
                }
                a_sv += a_step_s_v;
            }
            a_su += a_step_s_u;
        }

        let mut a_step = MathVector::new(1, a_nb_var);
        a_step.set_value(1, a_step_s_u);
        a_step.set_value(2, a_step_s_v);

        let mut a_value = 0.0;
        let mut a_pso = Pso::new(an_uvinf.clone(), an_uvsup.clone(), a_step, 32, 100);
        a_pso.perform_with_pool(&a_func, &mut a_particles, the_nb_particles, &mut a_value, &mut an_uv);

        // cxx:599-608.
        let mut a_ct = a_func.conic_parameter(&an_uv);
        if the_c.is_periodic()
            && (a_ct < the_tuvinf.value(1) - PCONFUSION || a_ct > the_tuvsup.value(1) + PCONFUSION)
        {
            a_ct = clib::in_period(
                a_ct,
                the_tuvinf.value(1),
                the_tuvinf.value(1) + 2.0 * std::f64::consts::PI,
            );
        }

        the_tuv.set_value(1, a_ct);
        the_tuv.set_value(2, an_uv.value(1));
        the_tuv.set_value(3, an_uv.value(2));

        // cxx:612-707: reject a "bad" solution whose distance vector is
        // (nearly) perpendicular to the surface normal.
        let mut is_bad_sol = false;
        let (mut a_p_on_s, mut a_du, mut a_dv) = s.d1(an_uv.value(1), an_uv.value(2));
        let (a_p_on_c, a_dt) = the_c.d1(a_ct);
        let a_sq_dist = a_p_on_c.square_distance(&a_p_on_s);
        if a_sq_dist <= SQUARE_CONFUSION {
            return;
        }
        let mut a_n = a_du.crossed(&a_dv);
        if a_n.square_magnitude() < SQUARE_CONFUSION {
            return;
        }

        let mut pc_ps = GpVec::from_pnts(&a_p_on_c, &a_p_on_s);
        let an_ang_min = std::f64::consts::FRAC_PI_2 - std::f64::consts::FRAC_PI_2 / 10.0;
        let an_ang_max = std::f64::consts::FRAC_PI_2 + std::f64::consts::FRAC_PI_2 / 10.0;

        let mut an_ang_n = pc_ps.angle(&a_n);
        if an_ang_n >= an_ang_min && an_ang_n <= an_ang_max {
            is_bad_sol = true;
            for iu in -1..=1 {
                let mut u = an_uv.value(1) + iu as f64 * a_step_s_u;
                u = u.max(an_uvinf.value(1)).min(an_uvsup.value(1));
                for iv in -1..=1 {
                    let mut v = an_uv.value(2) + iv as f64 * a_step_s_v;
                    v = v.max(an_uvinf.value(2)).min(an_uvsup.value(2));
                    let (p, du, dv) = s.d1(u, v);
                    a_p_on_s = p;
                    a_du = du;
                    a_dv = dv;
                    if a_p_on_c.square_distance(&a_p_on_s) < SQUARE_CONFUSION {
                        is_bad_sol = false;
                        break;
                    }
                    a_n = a_du.crossed(&a_dv);
                    if a_n.square_magnitude() < SQUARE_CONFUSION {
                        is_bad_sol = false;
                        break;
                    }
                    pc_ps = GpVec::from_pnts(&a_p_on_c, &a_p_on_s);
                    an_ang_n = pc_ps.angle(&a_n);
                    if an_ang_n < an_ang_min || an_ang_n > an_ang_max {
                        is_bad_sol = false;
                        break;
                    }
                }
                if !is_bad_sol {
                    break;
                }
            }
        }

        if is_bad_sol {
            // Try to precise solution with help of Extrema PS — cxx:680-706.
            let mut a_f = MathVector::new(1, 3);
            a_f.set_value(1, pc_ps.dot(&a_dt));
            a_f.set_value(2, pc_ps.dot(&a_du));
            a_f.set_value(3, pc_ps.dot(&a_dv));
            let a_ff = a_f.norm2();

            let mut an_ext_ps = GenLocateExtPs::new(s, self.tol2, self.tol2);
            an_ext_ps.perform(&a_p_on_c, an_uv.value(1), an_uv.value(2), false);
            if an_ext_ps.is_done() {
                if let Some(a_pmin) = an_ext_ps.point() {
                    let (u, v) = a_pmin.parameter();
                    let mut a_tuv = the_tuv.clone();
                    a_tuv.set_value(2, u);
                    a_tuv.set_value(3, v);
                    self.f.value(&a_tuv, &mut a_f);
                    let a_ff1 = a_f.norm2();
                    let a_sd = an_ext_ps.square_distance().unwrap_or(f64::MAX);
                    if a_sd < a_sq_dist && a_ff1 <= 1.1 * a_ff {
                        the_tuv.set_value(2, a_tuv.value(2));
                        the_tuv.set_value(3, a_tuv.value(3));
                    }
                }
            }
        }
    }
}

