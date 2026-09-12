//! `IntAna_IntQuadQuad` — cylinder/cone vs implicit quadric.
//! Source: `IntAna_IntQuadQuad.cxx` / `.hxx`.

use occt_core::elib::slib;
use occt_core::gp::{GpCone, GpCylinder, GpPnt};
use occt_core::precision::SQUARE_CONFUSION;

use super::curve_ana::IntAnaCurve;
use super::quadric::IntAnaQuadric;
use super::trig::{trig_function_roots, TrigRoots};

pub(crate) const MAX_CURVES: usize = 12;
const EPS: f64 = 1.0e-8;
const EPS_POLY: f64 = 1.0e-8;

#[path = "intana_intquadquad_cyl.rs"]
mod cyl;
#[path = "intana_intquadquad_cone.rs"]
mod cone;

/// Internal sorted trigonometric roots (`TrigonometricRoots` in IntQuadQuad.cxx).
#[derive(Clone)]
pub(crate) struct TrigPolyRoots {
    roots: [f64; 4],
    n: usize,
    done: bool,
    infinite: bool,
}

impl TrigPolyRoots {
    pub(crate) fn new(cc: f64, sc: f64, c: f64, s: f64, cte: f64, binf: f64, bsup: f64) -> Self {
        let two_pi = std::f64::consts::PI + std::f64::consts::PI;
        let mut out = Self {
            roots: [0.0; 4],
            n: 0,
            done: false,
            infinite: false,
        };
        match trig_function_roots(cc, sc, c, s, cte, binf, bsup) {
            TrigRoots::Fail => return out,
            TrigRoots::Infinite => {
                out.done = true;
                out.infinite = true;
                return out;
            }
            TrigRoots::Values(vals) => {
                out.done = true;
                out.n = vals.len().min(4);
                for (i, v) in vals.iter().take(4).enumerate() {
                    let mut r = *v;
                    if r < 0.0 {
                        r += two_pi;
                    }
                    if r > two_pi {
                        r -= two_pi;
                    }
                    out.roots[i] = r;
                }
            }
        }
        for i in 0..out.n {
            let co = out.roots[i].cos();
            let si = out.roots[i].sin();
            let y = co * (cc * co + (sc + sc) * si + c) + s * si + cte;
            if y.abs() > 1e-8 {
                out.done = false;
                return out;
            }
        }
        for i in 1..out.n {
            let mut j = i;
            while j > 0 && out.roots[j] < out.roots[j - 1] {
                out.roots.swap(j, j - 1);
                j -= 1;
            }
        }
        if out.n == 0 && (cc.abs() + sc.abs() + c.abs() + s.abs()) < 1e-10 && cte.abs() < 1e-10 {
            out.infinite = true;
        }
        out
    }

    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    pub(crate) fn infinite_roots(&self) -> bool {
        self.infinite
    }

    pub(crate) fn nb_solutions(&self) -> usize {
        self.n
    }

    pub(crate) fn value(&self, n1: usize) -> f64 {
        self.roots[n1.saturating_sub(1)]
    }

    pub(crate) fn is_a_root(&self, u: f64) -> bool {
        let two_pi = std::f64::consts::PI + std::f64::consts::PI;
        let eps = f64::EPSILON;
        self.roots[..self.n]
            .iter()
            .any(|&r| (u - r).abs() <= eps || (u - r - two_pi).abs() <= eps)
    }
}

/// `MyTrigonometricFunction` for the discriminant.
pub(crate) struct TrigFn {
    cc: f64,
    ss: f64,
    sc: f64,
    s: f64,
    c: f64,
    cte: f64,
}

impl TrigFn {
    pub(crate) fn new(cc: f64, ss: f64, sc: f64, c: f64, s: f64, cte: f64) -> Self {
        Self {
            cc,
            ss,
            sc,
            s,
            c,
            cte,
        }
    }

    pub(crate) fn value(&self, u: f64) -> f64 {
        let (si, co) = u.sin_cos();
        self.cc * co * co
            + self.ss * si * si
            + 2.0 * (si * (self.sc * co + self.s) + co * self.c)
            + self.cte
    }
}

pub(crate) fn add_special_points_cyl(quad: &IntAnaQuadric, cyl: &GpCylinder, t1: &mut f64, t2: &mut f64) {
    let period = std::f64::consts::PI + std::f64::consts::PI;
    let mut max_delta: f64 = 0.0;
    for pt in quad.special_points() {
        let (u, v) = slib::cylinder_parameters(&cyl.position(), pt);
        let proj = slib::cylinder_value(cyl, u, v);
        if pt.square_distance(&proj) > SQUARE_CONFUSION {
            continue;
        }
        let mut d1 = (*t1 - u).min(0.0);
        let mut d2 = (u - *t2).max(0.0);
        if d1 < -std::f64::consts::PI {
            d1 = 0.0;
        }
        if d2 > std::f64::consts::PI {
            d2 = 0.0;
        }
        max_delta = max_delta.max((-d1).max(d2));
    }
    if max_delta != 0.0 {
        *t1 -= max_delta;
        *t2 += max_delta;
        if *t2 - *t1 > period {
            *t2 = *t1 + period;
        }
    }
}

pub(crate) fn add_special_points_cone(quad: &IntAnaQuadric, cone: &GpCone, t1: &mut f64, t2: &mut f64) {
    let period = std::f64::consts::PI + std::f64::consts::PI;
    let mut max_delta: f64 = 0.0;
    for pt in quad.special_points() {
        let (u, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), pt);
        let proj = slib::cone_value(cone, u, v);
        if pt.square_distance(&proj) > SQUARE_CONFUSION {
            continue;
        }
        let mut d1 = (*t1 - u).min(0.0);
        let mut d2 = (u - *t2).max(0.0);
        if d1 < -std::f64::consts::PI {
            d1 = 0.0;
        }
        if d2 > std::f64::consts::PI {
            d2 = 0.0;
        }
        max_delta = max_delta.max((-d1).max(d2));
    }
    if max_delta != 0.0 {
        *t1 -= max_delta;
        *t2 += max_delta;
        if *t2 - *t1 > period {
            *t2 = *t1 + period;
        }
    }
}

/// `IntAna_IntQuadQuad`.
#[derive(Clone)]
pub struct IntQuadQuad {
    pub(crate) done: bool,
    pub(crate) identical: bool,
    pub(crate) curves: [IntAnaCurve; MAX_CURVES],
    pub(crate) previous: [i32; MAX_CURVES],
    pub(crate) next: [i32; MAX_CURVES],
    pub(crate) nb_curves: usize,
    pub(crate) nb_points: usize,
    pub(crate) points: [GpPnt; 2],
    pub(crate) epsilon: f64,
    pub(crate) epsilon_poly: f64,
}

impl Default for IntQuadQuad {
    fn default() -> Self {
        Self {
            done: false,
            identical: false,
            curves: std::array::from_fn(|_| IntAnaCurve::new()),
            previous: [0; MAX_CURVES],
            next: [0; MAX_CURVES],
            nb_curves: 0,
            nb_points: 0,
            points: [GpPnt::zero(), GpPnt::zero()],
            epsilon: EPS,
            epsilon_poly: EPS_POLY,
        }
    }
}

impl IntQuadQuad {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn perform_cylinder(&mut self, cyl: &GpCylinder, quad: &IntAnaQuadric, tol: f64) {
        cyl::perform(self, cyl, quad, tol);
    }

    pub fn perform_cone(&mut self, cone: &GpCone, quad: &IntAnaQuadric, tol: f64) {
        cone::perform(self, cone, quad, tol);
    }

    pub fn cylinder_quad(cyl: &GpCylinder, quad: &IntAnaQuadric, tol: f64) -> Self {
        let mut s = Self::new();
        s.perform_cylinder(cyl, quad, tol);
        s
    }

    pub fn cone_quad(cone: &GpCone, quad: &IntAnaQuadric, tol: f64) -> Self {
        let mut s = Self::new();
        s.perform_cone(cone, quad, tol);
        s
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn identical_elements(&self) -> bool {
        self.identical
    }

    pub fn nb_curve(&self) -> usize {
        if !self.done || self.identical {
            0
        } else {
            self.nb_curves
        }
    }

    pub fn nb_pnt(&self) -> usize {
        if !self.done || self.identical {
            0
        } else {
            self.nb_points
        }
    }

    pub fn curve(&self, n1: usize) -> &IntAnaCurve {
        &self.curves[n1.saturating_sub(1)]
    }

    pub fn point(&self, n1: usize) -> GpPnt {
        self.points[n1.saturating_sub(1)]
    }

    pub(crate) fn reset_links(&mut self) {
        self.previous = [0; MAX_CURVES];
        self.next = [0; MAX_CURVES];
    }

    pub(crate) fn internal_set_next_and_previous(&mut self) {
        let a_eps = 1.0e-7;
        let a_dist = 1.0e-10;
        self.reset_links();
        for c1 in 0..self.nb_curves {
            let (d_inf1, d_sup1) = self.curves[c1].domain();
            for c2 in 0..self.nb_curves {
                if c2 == c1 {
                    break;
                }
                let not_last2 = !self.curves[c2].is_last_open();
                let not_first2 = !self.curves[c2].is_first_open();
                let (d_inf2, d_sup2) = self.curves[c2].domain();
                if !self.curves[c1].is_first_open() {
                    if not_last2
                        && (d_inf1 - d_sup2).abs() <= a_eps
                        && self.curves[c1].value(d_inf1).distance(&self.curves[c2].value(d_sup2))
                            < a_dist
                    {
                        self.previous[c1] = (c2 + 1) as i32;
                        self.next[c2] = (c1 + 1) as i32;
                    }
                    if not_first2
                        && (d_inf1 - d_inf2).abs() <= a_eps
                        && self.curves[c1].value(d_inf1).distance(&self.curves[c2].value(d_inf2))
                            < a_dist
                    {
                        self.previous[c1] = -((c2 + 1) as i32);
                        self.previous[c2] = -((c1 + 1) as i32);
                    }
                }
                if !self.curves[c1].is_last_open() {
                    if not_last2
                        && (d_sup1 - d_sup2).abs() <= a_eps
                        && self.curves[c1].value(d_sup1).distance(&self.curves[c2].value(d_sup2))
                            < a_dist
                    {
                        self.next[c1] = -((c2 + 1) as i32);
                        self.next[c2] = -((c1 + 1) as i32);
                    }
                    if not_first2
                        && (d_sup1 - d_inf2).abs() <= a_eps
                        && self.curves[c1].value(d_sup1).distance(&self.curves[c2].value(d_inf2))
                            < a_dist
                    {
                        self.next[c1] = (c2 + 1) as i32;
                        self.previous[c2] = (c1 + 1) as i32;
                    }
                }
            }
        }
    }
}
