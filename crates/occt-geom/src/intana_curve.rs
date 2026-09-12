//! `IntAna_Curve` — parametric intersection of a cylinder/cone with a quadric.
//! Source: `IntAna_Curve.cxx` / `.hxx`.

use occt_core::elib::slib;
use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpPnt, GpSphere, GpTrsf, GpVec};
use occt_core::precision::{PCONFUSION, SQUARE_CONFUSION};
use std::f64::consts::PI;

use crate::curve::Curve;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QuadKind {
    Other,
    Cylinder,
    Cone,
    Sphere,
}

/// `IntAna_Curve`.
#[derive(Debug, Clone)]
pub struct IntAnaCurve {
    z0_cte: f64,
    z0_sin: f64,
    z0_cos: f64,
    z0_sin_sin: f64,
    z0_cos_cos: f64,
    z0_cos_sin: f64,
    z1_cte: f64,
    z1_sin: f64,
    z1_cos: f64,
    z1_sin_sin: f64,
    z1_cos_cos: f64,
    z1_cos_sin: f64,
    z2_cte: f64,
    z2_sin: f64,
    z2_cos: f64,
    z2_sin_sin: f64,
    z2_cos_cos: f64,
    z2_cos_sin: f64,
    two_curves: bool,
    take_z_positive: bool,
    tolerance: f64,
    domain_inf: f64,
    domain_sup: f64,
    restricted_inf: bool,
    restricted_sup: bool,
    first_bounded: bool,
    last_bounded: bool,
    kind: QuadKind,
    r_cyl: f64,
    angle: f64,
    ax3: GpAx3,
    my_first: f64,
    my_last: f64,
}

impl Default for IntAnaCurve {
    fn default() -> Self {
        Self {
            z0_cte: 0.0,
            z0_sin: 0.0,
            z0_cos: 0.0,
            z0_sin_sin: 0.0,
            z0_cos_cos: 0.0,
            z0_cos_sin: 0.0,
            z1_cte: 0.0,
            z1_sin: 0.0,
            z1_cos: 0.0,
            z1_sin_sin: 0.0,
            z1_cos_cos: 0.0,
            z1_cos_sin: 0.0,
            z2_cte: 0.0,
            z2_sin: 0.0,
            z2_cos: 0.0,
            z2_sin_sin: 0.0,
            z2_cos_cos: 0.0,
            z2_cos_sin: 0.0,
            two_curves: false,
            take_z_positive: false,
            tolerance: 0.0,
            domain_inf: 0.0,
            domain_sup: 0.0,
            restricted_inf: false,
            restricted_sup: false,
            first_bounded: false,
            last_bounded: false,
            kind: QuadKind::Other,
            r_cyl: 0.0,
            angle: 0.0,
            ax3: GpAx3::standard(),
            my_first: 0.0,
            my_last: 0.0,
        }
    }
}

impl IntAnaCurve {
    pub fn new() -> Self {
        Self::default()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn set_cone_quad_values(
        &mut self,
        cone: &GpCone,
        qxx: f64,
        qyy: f64,
        qzz: f64,
        qxy: f64,
        qxz: f64,
        qyz: f64,
        qx: f64,
        qy: f64,
        qz: f64,
        q1: f64,
        tol: f64,
        dom_inf: f64,
        dom_sup: f64,
        two_curves: bool,
        take_z_positive: bool,
    ) {
        self.ax3 = cone.position();
        self.r_cyl = cone.radius();
        self.angle = cone.semi_angle();
        let un_tg = 1.0 / cone.semi_angle().tan();
        self.kind = QuadKind::Cone;
        self.two_curves = two_curves;
        self.take_z_positive = take_z_positive;
        self.z0_cte = q1;
        self.z0_sin = 0.0;
        self.z0_cos = 0.0;
        self.z0_cos_cos = 0.0;
        self.z0_sin_sin = 0.0;
        self.z0_cos_sin = 0.0;
        self.z1_cte = 2.0 * un_tg * qz;
        self.z1_sin = qy + qy;
        self.z1_cos = qx + qx;
        self.z1_cos_cos = 0.0;
        self.z1_sin_sin = 0.0;
        self.z1_cos_sin = 0.0;
        self.z2_cte = qzz * un_tg * un_tg;
        self.z2_sin = (un_tg + un_tg) * qyz;
        self.z2_cos = (un_tg + un_tg) * qxz;
        self.z2_cos_cos = qxx;
        self.z2_sin_sin = qyy;
        self.z2_cos_sin = qxy;
        self.tolerance = tol;
        self.domain_inf = dom_inf;
        self.domain_sup = dom_sup;
        self.restricted_inf = true;
        self.restricted_sup = true;
        self.first_bounded = false;
        self.last_bounded = false;
        self.my_first = dom_inf;
        self.my_last = if two_curves {
            dom_sup + dom_sup - dom_inf
        } else {
            dom_sup
        };
    }

    #[allow(clippy::too_many_arguments)]
    pub fn set_cylinder_quad_values(
        &mut self,
        cyl: &GpCylinder,
        qxx: f64,
        qyy: f64,
        qzz: f64,
        qxy: f64,
        qxz: f64,
        qyz: f64,
        qx: f64,
        qy: f64,
        qz: f64,
        q1: f64,
        tol: f64,
        dom_inf: f64,
        dom_sup: f64,
        two_curves: bool,
        take_z_positive: bool,
    ) {
        self.ax3 = cyl.position();
        self.r_cyl = cyl.radius();
        self.kind = QuadKind::Cylinder;
        self.two_curves = two_curves;
        self.take_z_positive = take_z_positive;
        let r2 = self.r_cyl + self.r_cyl;
        self.z0_cte = q1;
        self.z0_sin = r2 * qy;
        self.z0_cos = r2 * qx;
        self.z0_cos_cos = qxx * self.r_cyl * self.r_cyl;
        self.z0_sin_sin = qyy * self.r_cyl * self.r_cyl;
        self.z0_cos_sin = self.r_cyl * self.r_cyl * qxy;
        self.z1_cte = qz + qz;
        self.z1_sin = r2 * qyz;
        self.z1_cos = r2 * qxz;
        self.z1_cos_cos = 0.0;
        self.z1_sin_sin = 0.0;
        self.z1_cos_sin = 0.0;
        self.z2_cte = qzz;
        self.z2_sin = 0.0;
        self.z2_cos = 0.0;
        self.z2_cos_cos = 0.0;
        self.z2_sin_sin = 0.0;
        self.z2_cos_sin = 0.0;
        self.tolerance = tol;
        self.domain_inf = dom_inf;
        self.domain_sup = dom_sup;
        self.restricted_inf = true;
        self.restricted_sup = true;
        self.first_bounded = false;
        self.last_bounded = false;
        self.my_first = dom_inf;
        self.my_last = if two_curves {
            dom_sup + dom_sup - dom_inf
        } else {
            dom_sup
        };
    }

    pub fn is_open(&self) -> bool {
        self.restricted_inf && self.restricted_sup
    }

    pub fn domain(&self) -> (f64, f64) {
        (self.my_first, self.my_last)
    }

    pub fn is_constant(&self) -> bool {
        false
    }

    pub fn is_first_open(&self) -> bool {
        self.first_bounded
    }

    pub fn is_last_open(&self) -> bool {
        self.last_bounded
    }

    pub fn set_is_first_open(&mut self, flag: bool) {
        self.first_bounded = flag;
    }

    pub fn set_is_last_open(&mut self, flag: bool) {
        self.last_bounded = flag;
    }

    pub fn set_domain(&mut self, first: f64, last: f64) {
        if last > first {
            self.my_first = first;
            self.my_last = last;
        }
    }

    pub fn value(&self, theta: f64) -> GpPnt {
        let uv = self.internal_uv(theta);
        self.internal_value(uv.0, uv.1)
    }

    pub fn d1u(&self, theta: f64) -> Option<(GpPnt, GpVec)> {
        let uv = self.internal_uv(theta);
        let pt = self.value(theta);
        if uv.2.abs() < 1.0e-7 || uv.3.abs() < 1.0e-10 {
            return None;
        }
        let mut dtheta = (self.domain_sup - self.domain_inf) * 1.0e-6;
        let mut theta2 = theta + dtheta;
        if theta2 < self.domain_inf
            || (theta2 > self.domain_sup && !self.two_curves)
            || theta2 > (self.domain_sup + self.domain_sup - self.domain_inf + 1.0e-14)
        {
            dtheta = -dtheta;
            theta2 = theta + dtheta;
        }
        let p2 = self.value(theta2);
        dtheta = 1.0 / dtheta;
        let vec = GpVec::new(
            (p2.x() - pt.x()) * dtheta,
            (p2.y() - pt.y()) * dtheta,
            (p2.z() - pt.z()) * dtheta,
        );
        let _ = uv;
        Some((pt, vec))
    }

    /// `IntAna_Curve::FindParameter`.
    pub fn find_parameter(&self, p: &GpPnt) -> Vec<f64> {
        let two_pi = 2.0 * PI;
        let eps_ang = 1.0e-8;
        let internal_prec = 1.0e-8;
        let mut theta = match self.kind {
            QuadKind::Cylinder => slib::cylinder_parameters(&self.ax3, p).0,
            QuadKind::Cone => slib::cone_parameters(&self.ax3, self.r_cyl, self.angle, p).0,
            _ => return Vec::new(),
        };
        if !self.first_bounded && self.domain_inf > theta && (self.domain_inf - theta) <= eps_ang {
            theta = self.domain_inf;
        } else if !self.last_bounded && theta > self.domain_sup && (theta - self.domain_sup) <= eps_ang
        {
            theta = self.domain_sup;
        }
        if theta < self.domain_inf {
            theta += two_pi;
        } else if theta > self.domain_sup {
            theta -= two_pi;
        }
        let mut params = [
            self.domain_inf,
            self.domain_sup,
            theta,
            if self.two_curves {
                self.domain_sup + self.domain_sup - theta
            } else {
                f64::MAX
            },
            if self.two_curves {
                self.domain_sup + self.domain_sup - self.domain_inf
            } else {
                f64::MAX
            },
        ];
        params[0..4].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut out = Vec::new();
        for i in 0..5 {
            if params[i] > self.my_last {
                break;
            }
            if params[i] < self.my_first {
                continue;
            }
            if i > 0 && (params[i] - params[i - 1]).abs() < PCONFUSION {
                continue;
            }
            let uv = self.internal_uv(params[i]);
            let q = self.internal_value(uv.0, uv.1);
            let sq_tol = if params[i] == theta
                || (self.two_curves && params[i] == self.domain_sup + self.domain_sup - theta)
            {
                internal_prec
            } else {
                SQUARE_CONFUSION
            };
            if q.square_distance(p) < sq_tol {
                out.push(params[i]);
            }
        }
        out
    }

    /// Returns (U, V, A, SigneSqrtDis).
    fn internal_uv(&self, theta_in: f64) -> (f64, f64, f64, f64) {
        let rel_p = 1.0 + f64::EPSILON;
        let rel_m = 1.0 - f64::EPSILON;
        let a_dt = 100.0 * f64::EPSILON * (self.domain_sup + self.domain_sup - self.domain_inf).abs().max(1.0);
        let mut theta = theta_in;
        let mut second = false;
        if theta < self.domain_inf * rel_m
            || (theta > self.domain_sup * rel_p && !self.two_curves)
            || theta > (self.domain_sup + self.domain_sup - self.domain_inf) * rel_p
        {
            theta = theta.clamp(self.my_first, self.my_last);
        }
        if (theta - self.domain_sup).abs() < a_dt {
            theta = self.domain_sup;
        } else if theta > self.domain_sup {
            theta = self.domain_sup + self.domain_sup - theta;
            second = true;
        }
        let param1 = theta;
        if !self.two_curves {
            second = self.take_z_positive;
        }
        let (sint, cost) = (theta.sin(), theta.cos());
        let sin2t = (theta + theta).sin();
        let cos2t = (theta + theta).cos();
        let a = self.z2_cte
            + sint * (self.z2_sin + sint * self.z2_sin_sin)
            + cost * (self.z2_cos + cost * self.z2_cos_cos)
            + self.z2_cos_sin * sin2t;
        let da = cost * self.z2_sin - sint * self.z2_cos
            + sin2t * (self.z2_sin_sin - self.z2_cos_cos)
            + cos2t * (self.z2_cos_sin * self.z2_cos_sin);
        let b = self.z1_cte
            + sint * (self.z1_sin + sint * self.z1_sin_sin)
            + cost * (self.z1_cos + cost * self.z1_cos_cos)
            + self.z1_cos_sin * sin2t;
        let db = self.z1_sin * cost - self.z1_cos * sint
            + sin2t * (self.z1_sin_sin - self.z1_cos_cos)
            + cos2t * (self.z1_cos_sin + self.z1_cos_sin);
        let c = self.z0_cte
            + sint * (self.z0_sin + sint * self.z0_sin_sin)
            + cost * (self.z0_cos + cost * self.z0_cos_cos)
            + self.z0_cos_sin * sin2t;
        let dc = self.z0_sin * cost - self.z0_cos * sint
            + sin2t * (self.z0_sin_sin - self.z0_cos_cos)
            + cos2t * (self.z0_cos_sin + self.z0_cos_sin);
        let mut disc = b * b - 4.0 * a * c;
        let tol_d = 2.0 * a_dt * (b * db - 2.0 * (a * dc + c * da)).abs();
        if disc < tol_d {
            disc = 0.0;
        }
        let mut signe = 0.0;
        let param2 = if a.abs() <= PCONFUSION {
            if b.abs() <= PCONFUSION {
                0.0
            } else {
                -c / b
            }
        } else {
            disc = disc.max(0.0);
            signe = if second { disc.sqrt() } else { -disc.sqrt() };
            (-b + signe) / (a + a)
        };
        (param1, param2, a, signe)
    }

    fn internal_value(&self, u: f64, v_in: f64) -> GpPnt {
        let v = v_in.clamp(-100000.0, 100000.0);
        match self.kind {
            QuadKind::Cone => {
                let co = GpCone {
                    pos: self.ax3,
                    radius: self.r_cyl,
                    semi_angle: self.angle,
                };
                slib::cone_value(&co, u, (v - self.r_cyl) / self.angle.sin())
            }
            QuadKind::Cylinder => {
                let cy = GpCylinder {
                    pos: self.ax3,
                    radius: self.r_cyl,
                };
                slib::cylinder_value(&cy, u, v)
            }
            QuadKind::Sphere => {
                let sp = GpSphere {
                    pos: self.ax3,
                    radius: self.r_cyl,
                };
                slib::sphere_value(&sp, u, v)
            }
            QuadKind::Other => GpPnt::new(0.0, 0.0, 0.0),
        }
    }
}

impl Curve for IntAnaCurve {
    fn d0(&self, u: f64) -> GpPnt {
        self.value(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.d1u(u).unwrap_or_else(|| (self.value(u), GpVec::zero()))
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1) = self.d1(u);
        (p, d1, GpVec::zero())
    }

    fn first_parameter(&self) -> f64 {
        self.my_first
    }

    fn last_parameter(&self) -> f64 {
        self.my_last
    }

    fn continuity(&self) -> u8 {
        1
    }

    fn transform(&mut self, _t: &GpTrsf) {}

    fn reverse(&mut self) {
        std::mem::swap(&mut self.my_first, &mut self.my_last);
    }

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
