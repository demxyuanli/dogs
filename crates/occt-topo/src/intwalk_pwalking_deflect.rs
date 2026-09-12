//! `IntWalk_PWalking::TestDeflection` / `TestArret` / `RepartirOuDiviser`.

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use super::iso::{ConstIso, StatusDeflection};
use super::pwalking::PWalking;

const COS_REF_2D: f64 = 0.939_692_620_785_908_4; // cos(pi/9)
const ANG_REF_2D: f64 = std::f64::consts::FRAC_PI_2;
const D_EXP: f64 = 7.0;

impl<'a> PWalking<'a> {
    /// Extra tangent check used at the first marching point.
    pub(crate) fn is_tangent_ext_check(&self, u1: f64, v1: f64, u2: f64, v2: f64) -> bool {
        let (p1, du1, dv1) = self.s1.d1(u1, v1);
        let (p2, du2, dv2) = self.s2.d1(u2, v2);
        let _ = (p1, p2);
        let n1 = du1.crossed(&dv1);
        let n2 = du2.crossed(&dv2);
        let sq1 = n1.square_magnitude();
        let sq2 = n2.square_magnitude();
        if sq1 < f64::EPSILON || sq2 < f64::EPSILON {
            return true;
        }
        let dp = n1.dot(&n2);
        if dp * dp < 0.9998 * sq1 * sq2 {
            return false;
        }
        let sq_tol = 4.0 * self.my_tol_tang * self.my_tol_tang;
        let par_u1 = [u1 + self.pasuv[0], u1 - self.pasuv[0], u1, u1];
        let par_v1 = [v1, v1, v1 + self.pasuv[1], v1 - self.pasuv[1]];
        let par_u2 = [u2 + self.pasuv[2], u2 - self.pasuv[2], u2, u2];
        let par_v2 = [v2, v2, v2 + self.pasuv[3], v2 - self.pasuv[3]];
        for i in 0..4 {
            let p = self.s1.d0(par_u1[i], par_v1[i]);
            if sq_dist_point_surface(&p, self.s2, u2, v2) > sq_tol {
                return false;
            }
        }
        for i in 0..4 {
            let p = self.s2.d0(par_u2[i], par_v2[i]);
            if sq_dist_point_surface(&p, self.s1, u1, v1) > sq_tol {
                return false;
            }
        }
        true
    }

    pub(crate) fn repartir_ou_diviser(
        &mut self,
        deja: &mut bool,
        choix: &mut ConstIso,
        arrive: &mut bool,
    ) {
        if *arrive {
            if !*deja {
                *arrive = false;
                *deja = true;
                self.reverse_restart(choix);
            }
            return;
        }
        if self.pasuv[0] * 0.5 < self.reso_u1
            && self.pasuv[1] * 0.5 < self.reso_v1
            && self.pasuv[2] * 0.5 < self.reso_u2
            && self.pasuv[3] * 0.5 < self.reso_v2
        {
            if !self.previous_tg {
                self.tglast = true;
            }
            if !*deja {
                *deja = true;
                self.reverse_restart(choix);
            } else {
                *arrive = true;
            }
        } else {
            for p in &mut self.pasuv {
                *p *= 0.5;
            }
        }
    }

    fn reverse_restart(&mut self, choix: &mut ConstIso) {
        if self.line.is_empty() {
            return;
        }
        self.previous_point = self.line[0];
        self.previous_tg = false;
        self.previous_d1 = self.first_d1;
        self.previous_d2 = self.first_d2;
        self.previous_d = self.tgdir;
        self.my_tangent_idx = self.line.len() as i32;
        self.tgdir.reverse();
        self.line.reverse();
        self.sens = -1;
        self.tgfirst = self.tglast;
        self.tglast = false;
        *choix = self.choix_iso_sav;
        let nn = self.line.len();
        if nn > 2 {
            let a = self.line[nn - 1].parameters();
            let b = self.line[nn - 2].parameters();
            self.pasuv = [
                (a.0 - b.0).abs(),
                (a.1 - b.1).abs(),
                (a.2 - b.2).abs(),
                (a.3 - b.3).abs(),
            ];
        }
    }

    pub(crate) fn test_arret(
        &mut self,
        deja: bool,
        param: &mut [f64; 4],
        choix: &mut ConstIso,
    ) -> bool {
        let eps = [self.reso_u1, self.reso_v1, self.reso_u2, self.reso_v2];
        let uvp = self.previous_point.parameters();
        let uvp = [uvp.0, uvp.1, uvp.2, uvp.3];
        let sol = self.inter.point().parameters();
        let sol = [sol.0, sol.1, sol.2, sol.3];
        let uvd = [self.um1, self.vm1, self.um2, self.vm2];
        let uvf = [self.um_1, self.vm_1, self.um_2, self.vm_2];
        let pair = [1, 0, 3, 2];
        let mut trouve = false;
        let mut duv = [-1.0; 4];
        let mut parc = *param;
        for i in 0..4 {
            let k = pair[i];
            if out_lo(param[i], sol[i], uvd[i], eps[i]) {
                trouve = true;
                let dpc = uvp[i] - param[i];
                let dpb = uvp[i] - uvd[i];
                parc[i] = uvd[i];
                let dv = param[k] - uvp[k];
                let dv2 = dv * dv;
                duv[i] = if dv2 > f64::EPSILON {
                    let t = dpc * dpb + dv2;
                    t * t / ((dpc * dpc + dv2) * (dpb * dpb + dv2))
                } else {
                    -1.0
                };
            } else if out_hi(param[i], sol[i], uvf[i], eps[i]) {
                trouve = true;
                let dpc = param[i] - uvp[i];
                let dpb = uvf[i] - uvp[i];
                parc[i] = uvf[i];
                let dv = param[k] - uvp[k];
                let dv2 = dv * dv;
                duv[i] = if dv2 > f64::EPSILON {
                    let t = dpc * dpb + dv2;
                    t * t / ((dpc * dpc + dv2) * (dpb * dpb + dv2))
                } else {
                    -1.0
                };
            } else {
                duv[i] = -1.0;
                parc[i] = param[i];
            }
        }
        if trouve {
            let mut ddv = -1.0;
            let mut k = -1i32;
            for i in 0..4 {
                param[i] = parc[i];
                if duv[i] > ddv {
                    ddv = duv[i];
                    k = i as i32;
                }
            }
            if k >= 0 {
                *choix = ConstIso::from_index(k);
            } else if on_bound(parc[0], uvd[0], uvf[0], eps[0]) {
                *choix = ConstIso::UOnS1;
            } else if on_bound(parc[1], uvd[1], uvf[1], eps[1]) {
                *choix = ConstIso::VOnS1;
            } else if on_bound(parc[2], uvd[2], uvf[2], eps[2]) {
                *choix = ConstIso::UOnS2;
            } else if on_bound(parc[3], uvd[3], uvf[3], eps[3]) {
                *choix = ConstIso::VOnS2;
            }
            self.close = false;
            return true;
        }
        if deja || self.line.is_empty() {
            self.close = false;
            return false;
        }
        let first = self.line[0];
        let prev = self.previous_point;
        let cur = self.inter.point();
        let close_s1 = close2d(
            first.u1, first.v1, prev.u1, prev.v1, cur.u1, cur.v1,
        );
        let close_s2 = close2d(
            first.u2, first.v2, prev.u2, prev.v2, cur.u2, cur.v2,
        );
        self.close = close_s1 && close_s2;
        self.close
    }

    pub(crate) fn test_deflection(
        &mut self,
        choix: ConstIso,
        _the_status: StatusDeflection,
    ) -> StatusDeflection {
        if self.line.len() == 1 {
            self.blocage = 0;
            self.precedent_inflexion = 0;
        }
        let mut a_status = StatusDeflection::Ok;
        if self.inter.is_tangent() {
            return StatusDeflection::ArretSurPoint;
        }
        let current = self.inter.point();
        let tg = self.inter.direction();
        let cos_tg = tg.xyz().x() * self.previous_d.xyz().x()
            + tg.xyz().y() * self.previous_d.xyz().y()
            + tg.xyz().z() * self.previous_d.xyz().z();
        if cos_tg < 0.0 {
            for p in &mut self.pasuv {
                *p *= 0.5;
            }
            self.precedent_inflexion += 3;
            if self.pasuv[0] < self.reso_u1
                && self.pasuv[1] < self.reso_v1
                && self.pasuv[2] < self.reso_u2
                && self.pasuv[3] < self.reso_v2
            {
                return StatusDeflection::ArretSurPointPrecedent;
            }
            return StatusDeflection::PasTropGrand;
        } else if self.precedent_inflexion > 0 {
            self.precedent_inflexion -= 1;
            return StatusDeflection::Ok;
        }

        let a_sq = self.previous_point.p.square_distance(&current.p);
        if a_sq < SQUARE_CONFUSION {
            for i in 0..4 {
                let reso = [self.reso_u1, self.reso_v1, self.reso_u2, self.reso_v2][i];
                self.pas_init[i] = self.pas_init[i].max(5.0 * reso);
                self.pasuv[i] = self.pasuv[i]
                    .max((1.5 * self.pasuv[i]).min(self.pas_init[i]));
            }
            a_status = StatusDeflection::PointConfondu;
        }

        let (up1, vp1, up2, vp2) = self.previous_point.parameters();
        let (uc1, vc1, uc2, vc2) = current.parameters();
        let du1 = uc1 - up1;
        let dv1 = vc1 - vp1;
        let du2 = uc2 - up2;
        let dv2 = vc2 - vp2;
        let abs_du1 = du1.abs();
        let abs_dv1 = dv1.abs();
        let abs_du2 = du2.abs();
        let abs_dv2 = dv2.abs();
        if abs_du1 < self.reso_u1
            && abs_dv1 < self.reso_v1
            && abs_du2 < self.reso_u2
            && abs_dv2 < self.reso_v2
        {
            self.pasuv = [self.reso_u1, self.reso_v1, self.reso_u2, self.reso_v2];
            return StatusDeflection::ArretSurPointPrecedent;
        }

        let mut tol_area = 100.0;
        if self.reso_u1 < PCONFUSION
            || self.reso_v1 < PCONFUSION
            || self.reso_u2 < PCONFUSION
            || self.reso_v2 < PCONFUSION
        {
            tol_area *= 2.0;
        }
        let cosi1 = du1 * self.previous_d1.x() + dv1 * self.previous_d1.y();
        let cosi2 = du2 * self.previous_d2.x() + dv2 * self.previous_d2.y();
        let duv1 = du1 * du1 + dv1 * dv1;
        let duv2 = du2 * du2 + dv2 * dv2;
        let reso_uv1 = self.reso_u1 * self.reso_u1 + self.reso_v1 * self.reso_v1;
        let reso_uv2 = self.reso_u2 * self.reso_u2 + self.reso_v2 * self.reso_v2;
        let min_div2 = CONFUSION * CONFUSION;
        let mut d1 = D_EXP;
        if duv1 > min_div2 {
            d1 = (reso_uv1 / duv1).abs().sqrt() * tol_area;
            d1 = d1.min(D_EXP);
        }
        let mut d2 = D_EXP;
        if duv2 > min_div2 {
            d2 = (reso_uv2 / duv2).abs().sqrt() * tol_area;
            d2 = d2.min(D_EXP);
        }
        let tol_c1 = d1.exp();
        let tol_c2 = d2.exp();
        let cos_ref1 = COS_REF_2D / tol_c1;
        let cos_ref2 = COS_REF_2D / tol_c2;
        if a_status != StatusDeflection::PointConfondu
            && (cosi1 * cosi1 < cos_ref1 * duv1 || cosi2 * cosi2 < cos_ref2 * duv2)
        {
            for p in &mut self.pasuv {
                *p *= 0.5;
            }
            if self.pasuv[0] < self.reso_u1
                && self.pasuv[1] < self.reso_v1
                && self.pasuv[2] < self.reso_u2
                && self.pasuv[3] < self.reso_v2
            {
                return StatusDeflection::ArretSurPointPrecedent;
            }
            for p in &mut self.pasuv {
                *p *= 0.5;
            }
            return StatusDeflection::PasTropGrand;
        }
        if a_status != StatusDeflection::PointConfondu {
            let tg1 = self.inter.direction_on_s1();
            let tg2 = self.inter.direction_on_s2();
            let c1 = du1 * tg1.x() + dv1 * tg1.y();
            let c2 = du2 * tg2.x() + dv2 * tg2.y();
            let ang1 = self.previous_d1.angle(&tg1).abs();
            let ang2 = self.previous_d2.angle(&tg2).abs();
            let ang_ref1 = ANG_REF_2D * tol_c1;
            let ang_ref2 = ANG_REF_2D * tol_c2;
            if c1 * c1 < cos_ref1 * duv1
                || c2 * c2 < cos_ref2 * duv2
                || ang1 > ang_ref1
                || ang2 > ang_ref2
            {
                for p in &mut self.pasuv {
                    *p *= 0.5;
                }
                if self.pasuv[0] < self.reso_u1
                    && self.pasuv[1] < self.reso_v1
                    && self.pasuv[2] < self.reso_u2
                    && self.pasuv[3] < self.reso_v2
                {
                    return StatusDeflection::ArretSurPoint;
                }
                return StatusDeflection::PasTropGrand;
            }
        }

        let prev_xyz = self.previous_d.xyz();
        let tg_xyz = tg.xyz();
        let dd = occt_core::gp::GpXyz::new(
            prev_xyz.x() - tg_xyz.x(),
            prev_xyz.y() - tg_xyz.y(),
            prev_xyz.z() - tg_xyz.z(),
        );
        let fleche_c = (dd.square_modulus() * a_sq).abs().sqrt() / 8.0;
        if fleche_c <= self.fleche * 0.5 {
            let mut ratio = if fleche_c > 1e-16 {
                0.5 * (self.fleche / fleche_c)
            } else {
                10.0
            };
            let pas_s = self.pasuv;
            self.pasuv[0] = self.pasuv[0].max(abs_du1).max(self.reso_u1);
            self.pasuv[1] = self.pasuv[1].max(abs_dv1).max(self.reso_v1);
            self.pasuv[2] = self.pasuv[2].max(abs_du2).max(self.reso_u2);
            self.pasuv[3] = self.pasuv[3].max(abs_dv2).max(self.reso_v2);
            let mut r = self.pas_init[0] / self.pasuv[0];
            for i in 1..4 {
                r = r.min(self.pas_init[i] / self.pasuv[i]);
            }
            if ratio > r {
                ratio = r;
            }
            for i in 0..4 {
                self.pasuv[i] = (ratio * self.pasuv[i]).min(self.pas_init[i]);
            }
            if self.pasuv != pas_s {
                self.blocage += 1;
                if self.blocage > 5 {
                    self.blocage = 0;
                    return StatusDeflection::PasTropGrand;
                }
            }
            if a_status == StatusDeflection::Ok {
                self.blocage = 0;
            }
            return a_status;
        }
        if fleche_c > self.fleche {
            let ratio = self.fleche / fleche_c;
            for p in &mut self.pasuv {
                *p *= ratio;
            }
            return StatusDeflection::PasTropGrand;
        }
        let _ = choix;
        a_status
    }
}

fn sq_dist_point_surface(p: &GpPnt, s: &dyn Surface, u0: f64, v0: f64) -> f64 {
    if let Some(q) = project_point_on_surface(s, p, 0.0) {
        let pq = s.d0(q.u, q.v);
        return p.square_distance(&pq);
    }
    let q = s.d0(u0, v0);
    p.square_distance(&q)
}

fn out_lo(param: f64, sol: f64, bound: f64, eps: f64) -> bool {
    bound.is_finite() && (param < bound - eps || sol < bound - eps)
}

fn out_hi(param: f64, sol: f64, bound: f64, eps: f64) -> bool {
    bound.is_finite() && (param > bound + eps || sol > bound + eps)
}

fn on_bound(p: f64, lo: f64, hi: f64, eps: f64) -> bool {
    (lo.is_finite() && p <= lo + eps) || (hi.is_finite() && p >= hi - eps)
}

fn close2d(u1: f64, v1: f64, up: f64, vp: f64, uc: f64, vc: f64) -> bool {
    let dx0 = u1 - up;
    let dy0 = v1 - vp;
    let dx1 = u1 - uc;
    let dy1 = v1 - vc;
    dx0 * dx1 + dy0 * dy1 < 0.0
}

#[allow(dead_code)]
fn _unused_dir(d: GpDir) -> GpVec {
    GpVec::from_xyz(d.xyz())
}
