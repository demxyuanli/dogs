//! `IntWalk_PWalking` constructor, first point, and accessors.

use occt_core::gp::{GpDir, GpDir2d};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::int_tools_wline::{u_resolution, v_resolution, PntOn2S};

use super::int2s::TheInt2S;
use super::iso::ConstIso;

/// Marching intersection of two parametric surfaces.
pub struct PWalking<'a> {
    pub(crate) s1: &'a dyn Surface,
    pub(crate) s2: &'a dyn Surface,
    pub(crate) done: bool,
    pub(crate) line: Vec<PntOn2S>,
    pub(crate) close: bool,
    pub(crate) tgfirst: bool,
    pub(crate) tglast: bool,
    pub(crate) my_tangent_idx: i32,
    pub(crate) tgdir: GpDir,
    pub(crate) fleche: f64,
    pub(crate) pas_max: f64,
    pub(crate) tolconf: f64,
    pub(crate) my_tol_tang: f64,
    pub(crate) pasuv: [f64; 4],
    pub(crate) my_step_min: [f64; 4],
    pub(crate) pas_sav: [f64; 4],
    pub(crate) pas_init: [f64; 4],
    pub(crate) um1: f64,
    pub(crate) um_1: f64,
    pub(crate) vm1: f64,
    pub(crate) vm_1: f64,
    pub(crate) um2: f64,
    pub(crate) um_2: f64,
    pub(crate) vm2: f64,
    pub(crate) vm_2: f64,
    pub(crate) reso_u1: f64,
    pub(crate) reso_v1: f64,
    pub(crate) reso_u2: f64,
    pub(crate) reso_v2: f64,
    pub(crate) sens: i32,
    pub(crate) choix_iso_sav: ConstIso,
    pub(crate) previous_point: PntOn2S,
    pub(crate) previous_tg: bool,
    pub(crate) previous_d: GpDir,
    pub(crate) previous_d1: GpDir2d,
    pub(crate) previous_d2: GpDir2d,
    pub(crate) first_d1: GpDir2d,
    pub(crate) first_d2: GpDir2d,
    pub(crate) inter: TheInt2S<'a>,
    pub(crate) blocage: i32,
    pub(crate) precedent_inflexion: i32,
}

impl<'a> PWalking<'a> {
    /// `IntWalk_PWalking(S1, S2, TolTangency, Epsilon, Deflection, Increment)`.
    pub fn new(
        s1: &'a dyn Surface,
        s2: &'a dyn Surface,
        tol_tangency: f64,
        epsilon: f64,
        deflection: f64,
        increment: f64,
    ) -> Self {
        const KELARG: f64 = 20.0;
        let pas_max = increment * 0.2;
        let (mut um1, mut um_1) = s1.u_range();
        let (mut vm1, mut vm_1) = s1.v_range();
        let (mut um2, mut um_2) = s2.u_range();
        let (mut vm2, mut vm_2) = s2.v_range();
        let mut reso_u1 = u_resolution(s1, CONFUSION);
        let mut reso_v1 = v_resolution(s1, CONFUSION);
        let mut reso_u2 = u_resolution(s2, CONFUSION);
        let mut reso_v2 = v_resolution(s2, CONFUSION);
        bump_reso(&mut reso_u1, um1, um_1);
        bump_reso(&mut reso_u2, um2, um_2);
        bump_reso(&mut reso_v1, vm1, vm_1);
        bump_reso(&mut reso_v2, vm2, vm_2);

        let mut pasuv = [
            pas_max * (um_1 - um1).abs(),
            pas_max * (vm_1 - vm1).abs(),
            pas_max * (um_2 - um2).abs(),
            pas_max * (vm_2 - vm2).abs(),
        ];
        for p in &mut pasuv {
            if !p.is_finite() {
                *p = pas_max.max(0.01);
            }
        }
        if reso_u1 > 0.0001 * pasuv[0] {
            reso_u1 = 0.00001 * pasuv[0];
        }
        if reso_v1 > 0.0001 * pasuv[1] {
            reso_v1 = 0.00001 * pasuv[1];
        }
        if reso_u2 > 0.0001 * pasuv[2] {
            reso_u2 = 0.00001 * pasuv[2];
        }
        if reso_v2 > 0.0001 * pasuv[3] {
            reso_v2 = 0.00001 * pasuv[3];
        }
        expand_periodic(s1.is_u_periodic(), &mut um1, &mut um_1, pasuv[0], KELARG);
        expand_periodic(s1.is_v_periodic(), &mut vm1, &mut vm_1, pasuv[1], KELARG);
        expand_periodic(s2.is_u_periodic(), &mut um2, &mut um_2, pasuv[2], KELARG);
        expand_periodic(s2.is_v_periodic(), &mut vm2, &mut vm_2, pasuv[3], KELARG);

        let my_step_min = [
            100.0 * reso_u1,
            100.0 * reso_v1,
            100.0 * reso_u2,
            100.0 * reso_v2,
        ];
        for p in &mut pasuv {
            if *p > 10.0 {
                *p = 10.0;
            }
        }
        let pas_init = pasuv;
        Self {
            s1,
            s2,
            done: true,
            line: Vec::new(),
            close: false,
            tgfirst: false,
            tglast: false,
            my_tangent_idx: 0,
            tgdir: GpDir::default_dir(),
            fleche: deflection,
            pas_max,
            tolconf: epsilon,
            my_tol_tang: tol_tangency,
            pasuv,
            my_step_min,
            pas_sav: pas_init,
            pas_init,
            um1,
            um_1,
            vm1,
            vm_1,
            um2,
            um_2,
            vm2,
            vm_2,
            reso_u1,
            reso_v1,
            reso_u2,
            reso_v2,
            sens: 1,
            choix_iso_sav: ConstIso::UOnS1,
            previous_point: PntOn2S {
                p: occt_core::gp::GpPnt::zero(),
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            },
            previous_tg: false,
            previous_d: GpDir::default_dir(),
            previous_d1: GpDir2d::default(),
            previous_d2: GpDir2d::default(),
            first_d1: GpDir2d::default(),
            first_d2: GpDir2d::default(),
            inter: TheInt2S::new(s1, s2, tol_tangency),
            blocage: 0,
            precedent_inflexion: 0,
        }
    }

    pub(crate) fn compute_pas_init(
        &mut self,
        du1: f64,
        dv1: f64,
        du2: f64,
        dv2: f64,
    ) {
        const RANGE_PART: f64 = 0.01;
        let increment = 2.0 * self.pas_max;
        let deltas = [
            (self.um_1 - self.um1).abs(),
            (self.vm_1 - self.vm1).abs(),
            (self.um_2 - self.um2).abs(),
            (self.vm_2 - self.vm2).abs(),
        ];
        let the_delta = [du1.abs(), dv1.abs(), du2.abs(), dv2.abs()];
        for i in 0..4 {
            if deltas[i].is_finite() {
                self.pasuv[i] = self.pasuv[i]
                    .max(increment * the_delta[i].max(RANGE_PART * deltas[i]));
            } else {
                self.pasuv[i] = self.pasuv[i].max(increment * the_delta[i]);
            }
        }
        let reso = [
            u_resolution(self.s1, self.tolconf),
            v_resolution(self.s1, self.tolconf),
            u_resolution(self.s2, self.tolconf),
            v_resolution(self.s2, self.tolconf),
        ];
        for i in 0..4 {
            self.my_step_min[i] = self.my_step_min[i].max(2.0 * reso[i]);
            self.pasuv[i] = self.pasuv[i].max(self.my_step_min[i]);
        }
    }

    /// `PerformFirstPoint`.
    pub fn perform_first_point(&mut self, par_dep: [f64; 4]) -> Option<PntOn2S> {
        self.sens = 1;
        self.close = false;
        self.inter.perform(par_dep);
        if !self.inter.is_done() || self.inter.is_empty() {
            return None;
        }
        Some(self.inter.point())
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_closed(&self) -> bool {
        self.close
    }

    pub fn nb_points(&self) -> i32 {
        self.line.len() as i32
    }

    pub fn value(&self, index: i32) -> PntOn2S {
        self.line[(index as usize).saturating_sub(1)]
    }

    pub fn line(&self) -> &[PntOn2S] {
        &self.line
    }

    pub fn tangent_at_first(&self) -> bool {
        self.tgfirst
    }

    pub fn tangent_at_last(&self) -> bool {
        self.tglast
    }

    pub fn tangent_at_line(&self) -> (GpDir, i32) {
        (self.tgdir, self.my_tangent_idx.max(1))
    }

    pub fn add_a_point(&mut self, p: PntOn2S) {
        self.line.push(p);
    }

    pub fn reverse_line(&mut self) {
        self.line.reverse();
    }
}

fn bump_reso(reso: &mut f64, a: f64, b: f64) {
    let max_val = a.abs().max(b.abs());
    let new_reso = *reso * max_val;
    if new_reso > *reso && new_reso < 10.0 {
        *reso = new_reso;
    }
}

fn expand_periodic(periodic: bool, lo: &mut f64, hi: &mut f64, pas: f64, kelarg: f64) {
    if !periodic || !lo.is_finite() || !hi.is_finite() {
        return;
    }
    let tspan = *hi - *lo;
    let period = tspan.max(std::f64::consts::TAU);
    if tspan < period {
        let mut t = 0.5 * (period - tspan);
        t = t.min(kelarg * pas);
        *hi += t;
        *lo -= t;
    }
}
