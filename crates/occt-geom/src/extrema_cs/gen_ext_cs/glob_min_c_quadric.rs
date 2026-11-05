//! `GlobMinCQuadric` — the 1-variable PSO seed for any curve against a
//! quadric (`Extrema_GenExtCS.cxx:711-823`). Split out of `gen_ext_cs/mod.rs`
//! to keep every file under the 800-line gate; the body is unchanged.

use super::*;

impl<'a> GenExtCS<'a> {
    /// `GlobMinCQuadric` (`cxx:711-823`).
    pub(super) fn glob_min_c_quadric(
        &self,
        the_c: &'a dyn Curve,
        the_nb_particles: usize,
        the_tuvinf: &MathVector,
        the_tuvsup: &MathVector,
        the_tuv: &mut MathVector,
    ) {
        let s = self.s.expect("Extrema_GenExtCS: no surface");
        let a_nb_var = 1;
        let mut a_tinf = MathVector::new(1, a_nb_var);
        let mut a_tsup = MathVector::new(1, a_nb_var);
        a_tinf.set_value(1, the_tuvinf.value(1));
        a_tsup.set_value(1, the_tuvsup.value(1));

        let mut a_particles = PsoParticlesPool::new(the_nb_particles, a_nb_var);

        let mut a_min_t = MathVector::new(1, a_nb_var);
        a_min_t.set_value(
            1,
            a_tinf.value(1) + (a_tsup.value(1) - a_tinf.value(1)) / A_BORDER_DIVISOR,
        );
        let mut a_max_t = MathVector::new(1, a_nb_var);
        a_max_t.set_value(
            1,
            a_tsup.value(1) - (a_tsup.value(1) - a_tinf.value(1)) / A_BORDER_DIVISOR,
        );

        // Increase numbers of curve samples to improve searching global
        // minimum because dimension of optimisation task is reduced —
        // cxx:731-739.
        const A_MAX_NB_NODES: usize = 50;
        let mut a_new_csample = self.tsample;
        let an_addsample = (self.usample / 2).max(3);
        a_new_csample += an_addsample;
        a_new_csample = a_new_csample.min(A_MAX_NB_NODES);

        let mut a_step_c_t = (a_max_t.value(1) - a_min_t.value(1)) / a_new_csample as f64;
        let a_step_s_u = (the_tuvsup.value(2) - the_tuvinf.value(2)) / self.usample as f64;
        let a_step_s_v = (the_tuvsup.value(3) - the_tuvinf.value(3)) / self.vsample as f64;
        let a_scale_factor = 5.0;
        let a_resolution_c_u = a_step_c_t / the_c.resolution(1.0);
        let u_res = crate::approx_same_parameter::u_resolution(s, 1.0);
        let v_res = crate::approx_same_parameter::v_resolution(s, 1.0);
        let a_min_resolution = a_scale_factor
            * a_resolution_c_u.min((a_step_s_u / u_res).min(a_step_s_v / v_res));
        if a_min_resolution > occt_core::precision::epsilon(1.0) && a_resolution_c_u > a_min_resolution
        {
            a_new_csample = A_MAX_NB_NODES.min(real_to_int(
                a_new_csample as f64 * a_resolution_c_u / a_min_resolution,
            ));
            a_step_c_t = (a_max_t.value(1) - a_min_t.value(1)) / a_new_csample as f64;
        }

        let mut a_func = GlobOptFuncCQuadric::with_range(the_c, a_tinf.value(1), a_tsup.value(1));
        a_func.load_quad(
            s,
            the_tuvinf.value(2),
            the_tuvsup.value(2),
            the_tuvinf.value(3),
            the_tuvsup.value(3),
        );

        let mut a_t = MathVector::new(1, a_nb_var);
        let mut a_particle = a_particles.worst_particle_index();
        let mut a_ct = a_min_t.value(1);
        for _ in 0..=a_new_csample {
            a_t.set_value(1, a_ct);
            let mut a_sq_dist = 0.0;
            if !a_func.value(&a_t, &mut a_sq_dist) {
                a_sq_dist = occt_core::precision::INFINITE;
            }
            if a_sq_dist < a_particles.particle(a_particle).distance {
                let p = a_particles.particle_mut(a_particle);
                p.position[0] = a_ct;
                p.best_position[0] = a_ct;
                p.distance = a_sq_dist;
                p.best_distance = a_sq_dist;
                a_particle = a_particles.worst_particle_index();
            }
            a_ct += a_step_c_t;
        }

        let mut a_step = MathVector::new(1, a_nb_var);
        a_step.set_value(1, a_step_c_t);

        let mut a_value = 0.0;
        let mut a_pso = Pso::new(a_tinf.clone(), a_tsup.clone(), a_step, 32, 100);
        a_pso.perform_with_pool(&a_func, &mut a_particles, the_nb_particles, &mut a_value, &mut a_t);

        let (mut an_u, mut an_v) = a_func.quadric_parameters(a_t.value(1));
        if s.is_u_periodic()
            && (an_u < the_tuvinf.value(2) - PCONFUSION
                || an_u > the_tuvsup.value(2) + PCONFUSION)
        {
            an_u = clib::in_period(
                an_u,
                the_tuvinf.value(2),
                the_tuvinf.value(2) + 2.0 * std::f64::consts::PI,
            );
        }
        if s.is_v_periodic()
            && (an_v < the_tuvinf.value(3) - PCONFUSION
                || an_v > the_tuvsup.value(3) + PCONFUSION)
        {
            an_v = clib::in_period(
                an_v,
                the_tuvinf.value(3),
                the_tuvinf.value(3) + 2.0 * std::f64::consts::PI,
            );
        }

        the_tuv.set_value(1, a_t.value(1));
        the_tuv.set_value(2, an_u);
        the_tuv.set_value(3, an_v);
    }
}

