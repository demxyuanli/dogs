//! `IntWalk_PWalking::Perform`.

use super::iso::{ConstIso, StatusDeflection};
use super::pwalking::PWalking;

impl<'a> PWalking<'a> {
    /// `Perform(ParDep)` using the natural UV box.
    pub fn perform(&mut self, par_dep: [f64; 4]) {
        self.perform_box(
            par_dep,
            self.um1,
            self.vm1,
            self.um2,
            self.vm2,
            self.um_1,
            self.vm_1,
            self.um_2,
            self.vm_2,
        );
    }

    /// `Perform(ParDep, u1min..v2max)`.
    pub fn perform_box(
        &mut self,
        par_dep: [f64; 4],
        u1min: f64,
        v1min: f64,
        u2min: f64,
        v2min: f64,
        u1max: f64,
        v1max: f64,
        u2max: f64,
        v2max: f64,
    ) {
        const SQ_DIST_MAX: f64 = 1.0e-14;
        const REJECT_MAX: i32 = 250_000;
        self.done = false;
        self.compute_pas_init(u1max - u1min, v1max - v1min, u2max - u2min, v2max - v2min);
        for i in 0..4 {
            if self.pasuv[i] > 10.0 {
                self.pasuv[i] = 10.0;
            }
            self.pas_init[i] = self.pasuv[i];
            self.pas_sav[i] = self.pasuv[i];
        }
        self.line.clear();
        let mut param = par_dep;
        let mut choix = self.inter.perform(param);
        if !self.inter.is_done() || self.inter.is_empty() || self.inter.is_tangent() {
            return;
        }
        let mut deja = false;
        let mut inc_key = 0i32;
        let mut reject_index = 0i32;
        self.previous_point = self.inter.point();
        self.previous_tg = false;
        self.previous_d = self.inter.direction();
        self.previous_d1 = self.inter.direction_on_s1();
        self.previous_d2 = self.inter.direction_on_s2();
        self.my_tangent_idx = 1;
        self.tgdir = self.previous_d;
        self.first_d1 = self.previous_d1;
        self.first_d2 = self.previous_d2;
        self.tgfirst = false;
        self.tglast = false;
        self.choix_iso_sav = choix;
        let pf = self.previous_point.p;
        let mut test_first = true;
        param = [
            self.previous_point.u1,
            self.previous_point.v1,
            self.previous_point.u2,
            self.previous_point.v2,
        ];
        if self.is_tangent_ext_check(param[0], param[1], param[2], param[3]) {
            return;
        }
        self.add_a_point(self.previous_point);

        let (u_first1, u_last1) = self.s1.u_range();
        let (v_first1, v_last1) = self.s1.v_range();
        let (u_first2, u_last2) = self.s2.u_range();
        let (v_first2, v_last2) = self.s2.v_range();
        let a_tol = [
            eps_of(self.um_1 - self.um1),
            eps_of(self.vm_1 - self.vm1),
            eps_of(self.um_2 - self.um2),
            eps_of(self.vm_2 - self.vm2),
        ];
        let mut a_status = StatusDeflection::Ok;
        let mut prev_status = StatusDeflection::Ok;
        let mut no_test = false;
        let mut empty_level = 0i32;
        let mut confondu_level = 0i32;
        let mut without_append = -1i32;
        let mut nb_ok = 0i32;
        let mut prev_nb = self.line.len();
        let mut arrive = false;

        while !arrive {
            prev_status = a_status;
            without_append += 1;
            if without_append > 20 {
                arrive = true;
                if deja {
                    break;
                }
                self.repartir_ou_diviser(&mut deja, &mut choix, &mut arrive);
                without_append = 0;
            }
            let mut f = match choix {
                ConstIso::UOnS1 => self.previous_d1.x().abs(),
                ConstIso::VOnS1 => self.previous_d1.y().abs(),
                ConstIso::UOnS2 => self.previous_d2.x().abs(),
                ConstIso::VOnS2 => self.previous_d2.y().abs(),
            };
            if f < 0.1 {
                f = 0.1;
            }
            param = [
                self.previous_point.u1,
                self.previous_point.v1,
                self.previous_point.u2,
                self.previous_point.v2,
            ];
            let sens = self.sens as f64;
            let mut dp = [
                sens * self.pasuv[0] * self.previous_d1.x() / f,
                sens * self.pasuv[1] * self.previous_d1.y() / f,
                sens * self.pasuv[2] * self.previous_d2.x() / f,
                sens * self.pasuv[3] * self.previous_d2.y() / f,
            ];
            let inc = 5.0 * inc_key as f64;
            const A_EPS: f64 = 1.0e-7;
            match choix {
                ConstIso::UOnS1 if dp[0].abs() < A_EPS => dp[0] *= inc,
                ConstIso::VOnS1 if dp[1].abs() < A_EPS => dp[1] *= inc,
                ConstIso::UOnS2 if dp[2].abs() < A_EPS => dp[2] *= inc,
                ConstIso::VOnS2 if dp[3].abs() < A_EPS => dp[3] *= inc,
                _ => {}
            }
            for i in 0..4 {
                param[i] += dp[i];
            }
            let sv = param;
            let mut try_n = 0;
            let mut best_iso = choix;
            loop {
                let mut bad = false;
                choix = self.inter.perform_iso(param, best_iso);
                if self.inter.is_done() && !self.inter.is_empty() {
                    let mut np = [
                        self.inter.point().u1,
                        self.inter.point().v1,
                        self.inter.point().u2,
                        self.inter.point().v2,
                    ];
                    let pmin = [self.um1, self.vm1, self.um2, self.vm2];
                    let pmax = [self.um_1, self.vm_1, self.um_2, self.vm_2];
                    for i in 0..4 {
                        if pmin[i].is_finite() && (np[i] - pmin[i]).abs() < a_tol[i] {
                            np[i] = pmin[i];
                        } else if pmax[i].is_finite() && (np[i] - pmax[i]).abs() < a_tol[i] {
                            np[i] = pmax[i];
                        }
                    }
                    if out_box(np, pmin, pmax) {
                        break;
                    }
                    {
                        let p = self.inter.point();
                        self.inter.change_point().set_value(p.p, np[0], np[1], np[2], np[3]);
                    }
                    let dist = [
                        (param[0] - dp[0] - np[0]).abs(),
                        (param[1] - dp[1] - np[1]).abs(),
                        (param[2] - dp[2] - np[2]).abs(),
                        (param[3] - dp[3] - np[3]).abs(),
                    ];
                    if dist[0] < self.reso_u1
                        && dist[1] < self.reso_v1
                        && dist[2] < self.reso_u2
                        && dist[3] < self.reso_v2
                        && a_status != StatusDeflection::PasTropGrand
                    {
                        bad = true;
                        best_iso = best_iso.next();
                    }
                }
                try_n += 1;
                if !bad || try_n > 4 {
                    break;
                }
            }

            if !self.inter.is_done() {
                arrive = false;
                param = sv;
                self.repartir_ou_diviser(&mut deja, &mut choix, &mut arrive);
                continue;
            }
            if self.inter.is_empty() {
                let (u1, v1, u2, v2) = self.previous_point.parameters();
                arrive = out_natural(u1, u_first1, u_last1)
                    || out_natural(u2, u_first2, u_last2)
                    || out_natural(v1, v_first1, v_last1)
                    || out_natural(v2, v_first2, v_last2);
                self.repartir_ou_diviser(&mut deja, &mut choix, &mut arrive);
                empty_level += 1;
                if empty_level > 10 {
                    self.pasuv = self.pas_sav;
                }
                continue;
            }
            if no_test {
                no_test = false;
            } else if {
                empty_level -= 1;
                empty_level <= 0
            } {
                empty_level = 0;
                if without_append < 10 {
                    a_status = self.test_deflection(choix, a_status);
                } else if a_status != StatusDeflection::StepTooSmall {
                    for p in &mut self.pasuv {
                        *p *= 0.5;
                    }
                }
            }
            if confondu_level > 5 {
                a_status = StatusDeflection::ArretSurPoint;
                confondu_level = 0;
            }
            if a_status == StatusDeflection::Ok {
                nb_ok += 1;
                if nb_ok >= 5 {
                    nb_ok = 0;
                    grow_steps(&mut self.pasuv, &mut self.pas_init, &mut self.pas_max);
                }
            } else {
                nb_ok = 0;
            }
            match a_status {
                StatusDeflection::ArretSurPointPrecedent => {
                    arrive = false;
                    self.repartir_ou_diviser(&mut deja, &mut choix, &mut arrive);
                }
                StatusDeflection::PasTropGrand => {
                    param = sv;
                    if without_append > 5 {
                        shrink_init(&mut self.pas_init, &self.pas_sav);
                        if prev_status != StatusDeflection::StepTooSmall
                            && self.line.len() != prev_nb
                        {
                            without_append = 0;
                        }
                        prev_nb = self.line.len();
                    }
                }
                StatusDeflection::PointConfondu => {
                    confondu_level += 1;
                    if confondu_level > 5 {
                        grow_steps(&mut self.pasuv, &mut self.pas_init, &mut self.pas_max);
                    }
                }
                StatusDeflection::StepTooSmall => {
                    let mut grew = false;
                    for i in 0..4 {
                        let ns = (1.5 * self.pasuv[i]).min(self.pas_init[i]);
                        if ns > self.pasuv[i] {
                            self.pasuv[i] = ns;
                            grew = true;
                        }
                    }
                    if grew {
                        param = sv;
                        if prev_status != StatusDeflection::PasTropGrand
                            && self.line.len() != prev_nb
                        {
                            without_append = 0;
                        }
                        prev_nb = self.line.len();
                    } else {
                        self.append_or_stop(
                            &mut arrive,
                            &mut deja,
                            &mut choix,
                            &mut param,
                            a_status,
                            prev_status,
                            &mut without_append,
                            &mut nb_ok,
                            &mut inc_key,
                            &mut reject_index,
                            &mut test_first,
                            pf,
                            SQ_DIST_MAX,
                            REJECT_MAX,
                        );
                    }
                }
                StatusDeflection::Ok | StatusDeflection::ArretSurPoint => {
                    self.append_or_stop(
                        &mut arrive,
                        &mut deja,
                        &mut choix,
                        &mut param,
                        a_status,
                        prev_status,
                        &mut without_append,
                        &mut nb_ok,
                        &mut inc_key,
                        &mut reject_index,
                        &mut test_first,
                        pf,
                        SQ_DIST_MAX,
                        REJECT_MAX,
                    );
                }
            }
        }
        self.done = self.line.len() >= 2;
    }

    fn append_or_stop(
        &mut self,
        arrive: &mut bool,
        deja: &mut bool,
        choix: &mut ConstIso,
        param: &mut [f64; 4],
        a_status: StatusDeflection,
        prev_status: StatusDeflection,
        without_append: &mut i32,
        nb_ok: &mut i32,
        inc_key: &mut i32,
        reject_index: &mut i32,
        test_first: &mut bool,
        pf: occt_core::gp::GpPnt,
        sq_max: f64,
        reject_max: i32,
    ) {
        *arrive = self.test_arret(*deja, param, choix);
        if !*arrive && a_status == StatusDeflection::ArretSurPoint {
            *arrive = true;
        }
        if *arrive {
            *nb_ok = -10;
        }
        if !*arrive {
            let (u1, v1, u2, v2) = self.inter.point().parameters();
            let valid = in_box(
                [u1, v1, u2, v2],
                [self.um1, self.vm1, self.um2, self.vm2],
                [self.um_1, self.vm_1, self.um_2, self.vm_2],
            );
            if valid {
                self.previous_point = self.inter.point();
                self.previous_tg = self.inter.is_tangent();
                if !self.previous_tg {
                    self.previous_d = self.inter.direction();
                    self.previous_d1 = self.inter.direction_on_s1();
                    self.previous_d2 = self.inter.direction_on_s2();
                }
                let (u1, v1, u2, v2) = self.previous_point.parameters();
                if in_box(
                    [u1, v1, u2, v2],
                    [self.um1, self.vm1, self.um2, self.vm2],
                    [self.um_1, self.vm_1, self.um_2, self.vm_2],
                ) {
                    let pl = self.previous_point.p;
                    if *test_first {
                        if pf.square_distance(&pl) < sq_max {
                            *inc_key += 1;
                            if *inc_key == 5000 {
                                *arrive = true;
                            }
                            return;
                        }
                        *test_first = false;
                    } else if pf.square_distance(&pl) < sq_max {
                        self.close = true;
                        *arrive = true;
                        if !self.line.is_empty() {
                            self.add_a_point(self.line[0]);
                        }
                        *without_append = 0;
                        return;
                    }
                    self.add_a_point(self.previous_point);
                    *reject_index += 1;
                    if *reject_index >= reject_max {
                        *arrive = true;
                        return;
                    }
                    *without_append = 0;
                }
            }
            if a_status == StatusDeflection::ArretSurPoint {
                self.repartir_ou_diviser(deja, choix, arrive);
            } else if self.line.len() == 2 {
                self.pas_sav = self.pasuv;
                if prev_status == StatusDeflection::PasTropGrand && *without_append > 0 {
                    self.pas_init = self.pasuv;
                }
            }
        } else if self.close {
            if !self.line.is_empty() {
                self.add_a_point(self.line[0]);
            }
            *without_append = 0;
        } else {
            self.repartir_ou_diviser(deja, choix, arrive);
        }
    }
}

fn eps_of(x: f64) -> f64 {
    if x.is_finite() {
        (x.abs() * f64::EPSILON).max(f64::EPSILON)
    } else {
        f64::EPSILON
    }
}

fn out_box(p: [f64; 4], lo: [f64; 4], hi: [f64; 4]) -> bool {
    for i in 0..4 {
        if lo[i].is_finite() && p[i] < lo[i] {
            return true;
        }
        if hi[i].is_finite() && p[i] > hi[i] {
            return true;
        }
    }
    false
}

fn in_box(p: [f64; 4], lo: [f64; 4], hi: [f64; 4]) -> bool {
    !out_box(p, lo, hi)
}

fn out_natural(x: f64, lo: f64, hi: f64) -> bool {
    (lo.is_finite() && x < lo) || (hi.is_finite() && x > hi)
}

fn grow_steps(pasuv: &mut [f64; 4], pas_init: &mut [f64; 4], pas_max: &mut f64) {
    loop {
        let mut too_small = true;
        for i in 0..4 {
            if pasuv[i] < pas_init[i] {
                let mut t = (pas_init[i] - pasuv[i]) * 0.25;
                if t > 0.1 * pas_init[i] {
                    t = 0.1 * pasuv[i];
                }
                pasuv[i] += t;
                too_small = false;
            }
        }
        if too_small {
            if *pas_max < 0.1 {
                *pas_max *= 1.1;
                for p in pas_init.iter_mut() {
                    *p *= 1.1;
                }
            } else {
                break;
            }
        } else {
            break;
        }
    }
}

fn shrink_init(pas_init: &mut [f64; 4], pas_sav: &[f64; 4]) {
    for i in 0..4 {
        if pas_sav[i] > pas_init[i] {
            continue;
        }
        let delta = (pas_init[i] - pas_sav[i]) * 0.25;
        if delta > f64::EPSILON * pas_init[i].abs() {
            pas_init[i] -= delta;
        }
    }
}
