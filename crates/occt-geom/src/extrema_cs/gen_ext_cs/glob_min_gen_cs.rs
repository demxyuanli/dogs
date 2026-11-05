//! `GlobMinGenCS` — the 3-variable PSO seed for `Extrema_GenExtCS`
//! (`Extrema_GenExtCS.cxx:427-518`). Split out of `gen_ext_cs/mod.rs` to keep
//! every file under the 800-line gate; the body is unchanged.

use super::*;

impl<'a> GenExtCS<'a> {
    /// `GlobMinGenCS` (`cxx:427-518`).
    pub(super) fn glob_min_gen_cs(
        &self,
        the_c: &'a dyn Curve,
        the_nb_particles: usize,
        the_tuvinf: &MathVector,
        the_tuvsup: &MathVector,
        the_tuv: &mut MathVector,
    ) {
        let s = self.s.expect("Extrema_GenExtCS: no surface");
        let mut a_particles = PsoParticlesPool::new(the_nb_particles, 3);

        let mut a_min_tuv = MathVector::new(1, 3);
        for i in 1..=3 {
            a_min_tuv.set_value(
                i,
                the_tuvinf.value(i) + (the_tuvsup.value(i) - the_tuvinf.value(i)) / A_BORDER_DIVISOR,
            );
        }
        let mut a_max_tuv = MathVector::new(1, 3);
        for i in 1..=3 {
            a_max_tuv.set_value(
                i,
                the_tuvsup.value(i) - (the_tuvsup.value(i) - the_tuvinf.value(i)) / A_BORDER_DIVISOR,
            );
        }

        let mut a_step_c_u = (a_max_tuv.value(1) - a_min_tuv.value(1)) / self.tsample as f64;
        let a_step_s_u = (a_max_tuv.value(2) - a_min_tuv.value(2)) / self.usample as f64;
        let a_step_s_v = (a_max_tuv.value(3) - a_min_tuv.value(3)) / self.vsample as f64;

        // Correct number of curve samples in case of low resolution —
        // cxx:445-466.
        let mut a_new_csample = self.tsample;
        let a_scale_factor = 5.0;
        let a_resolution_c_u = a_step_c_u / the_c.resolution(1.0);
        let u_res = crate::approx_same_parameter::u_resolution(s, 1.0);
        let v_res = crate::approx_same_parameter::v_resolution(s, 1.0);
        let a_min_resolution = a_scale_factor
            * a_resolution_c_u.min((a_step_s_u / u_res).min(a_step_s_v / v_res));
        if a_min_resolution > occt_core::precision::epsilon(1.0) && a_resolution_c_u > a_min_resolution
        {
            const A_MAX_NB_NODES: usize = 50;
            a_new_csample = A_MAX_NB_NODES.min(real_to_int(
                self.tsample as f64 * a_resolution_c_u / a_min_resolution,
            ));
            a_step_c_u = (a_max_tuv.value(1) - a_min_tuv.value(1)) / a_new_csample as f64;
        }

        // Pre-compute curve sample points — cxx:468-475.
        let mut a_curv_pnts = Vec::with_capacity(a_new_csample + 1);
        let mut a_cu1 = a_min_tuv.value(1);
        for _ in 0..=a_new_csample {
            a_curv_pnts.push(the_c.d0(a_cu1));
            a_cu1 += a_step_c_u;
        }

        // Select specified number of particles from pre-computed set of
        // samples — cxx:477-506.
        let stride = self.vsample + 1;
        let mut a_particle = a_particles.worst_particle_index();
        let mut a_su = a_min_tuv.value(2);
        for a_sui in 0..=self.usample {
            let mut a_sv = a_min_tuv.value(3);
            for a_svi in 0..=self.vsample {
                let a_surf_pnt = self.surf_pnts[a_sui * stride + a_svi];
                let mut a_cu2 = a_min_tuv.value(1);
                for a_curv_pnt in a_curv_pnts.iter() {
                    let a_sq_dist = a_surf_pnt.square_distance(a_curv_pnt);
                    if a_sq_dist < a_particles.particle(a_particle).distance {
                        let p = a_particles.particle_mut(a_particle);
                        p.position[0] = a_cu2;
                        p.position[1] = a_su;
                        p.position[2] = a_sv;
                        p.best_position[0] = a_cu2;
                        p.best_position[1] = a_su;
                        p.best_position[2] = a_sv;
                        p.distance = a_sq_dist;
                        p.best_distance = a_sq_dist;
                        a_particle = a_particles.worst_particle_index();
                    }
                    a_cu2 += a_step_c_u;
                }
                a_sv += a_step_s_v;
            }
            a_su += a_step_s_u;
        }

        let mut a_step = MathVector::new(1, 3);
        a_step.set_value(1, a_step_c_u);
        a_step.set_value(2, a_step_s_u);
        a_step.set_value(3, a_step_s_v);

        // Find min approximation — cxx:512-518.
        let a_func = GlobOptFuncCS::new(the_c, s);
        let mut a_value = 0.0;
        let mut a_pso = Pso::new(the_tuvinf.clone(), the_tuvsup.clone(), a_step, 32, 100);
        a_pso.perform_with_pool(&a_func, &mut a_particles, the_nb_particles, &mut a_value, the_tuv);
    }
}

