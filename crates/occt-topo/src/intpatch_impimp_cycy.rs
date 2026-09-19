//! Cylinder-cylinder non-geometric walking (`CyCyNoGeometric`).
//! Source: `IntPatch_ImpImpIntersection.cxx` ComputationMethods,
//! WorkWithBoundaries, and the CyCyNoGeometric walker.

use std::f64::consts::PI;

use occt_core::gp::{GpAx1, GpCylinder, GpXyz};
use occt_core::precision::{ANGULAR, CONFUSION, Precision, REAL_SMALL};

use super::glines::PairOutcome;

#[path = "intpatch_impimp_cycy_bounds.rs"]
mod bounds;
#[path = "intpatch_impimp_cycy_step.rs"]
mod step;
#[path = "intpatch_impimp_cycy_wl.rs"]
mod wl;
#[path = "intpatch_impimp_cycy_walk.rs"]
mod walk;

pub(crate) const NUL_VALUE: f64 = 1.0e-11;
pub(crate) const PERIOD: f64 = 2.0 * PI;

/// UV box of one cylinder (`Bnd_Box2d` of First/Last U,V).
#[derive(Clone, Copy)]
pub(crate) struct CylUv {
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

/// Surface UV box used while walking U1.
#[derive(Clone, Copy)]
pub(crate) struct SurfDom {
    pub u1f: f64,
    pub u1l: f64,
    pub u2f: f64,
    pub u2l: f64,
    pub v1f: f64,
    pub v1l: f64,
    pub v2f: f64,
    pub v2l: f64,
    pub period: f64,
    pub tol2d: f64,
    pub tol3d: f64,
}

/// Equation coefficients of `ComputationMethods::stCoeffsValue`.
#[derive(Clone, Copy)]
pub(crate) struct Coeffs {
    pub vec_a1: GpXyz,
    pub vec_a2: GpXyz,
    pub vec_b1: GpXyz,
    pub vec_b2: GpXyz,
    pub vec_c1: GpXyz,
    pub vec_c2: GpXyz,
    pub vec_d: GpXyz,
    pub k21: f64,
    pub k11: f64,
    pub l21: f64,
    pub l11: f64,
    pub m1: f64,
    pub k22: f64,
    pub k12: f64,
    pub l22: f64,
    pub l12: f64,
    pub m2: f64,
    pub b: f64,
    pub c: f64,
    pub fi1: f64,
    pub fi2: f64,
}

fn xyz_comp(v: &GpXyz, i: usize) -> f64 {
    match i {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}

fn rotate_123(v: &mut GpXyz) {
    let (x, y, z) = (v.x, v.y, v.z);
    v.x = y;
    v.y = z;
    v.z = x;
}

fn swap_yz(v: &mut GpXyz) {
    let (y, z) = (v.y, v.z);
    v.y = z;
    v.z = y;
}

fn short_cos_form(cos_f: f64, sin_f: f64) -> (f64, f64) {
    let coeff = (cos_f * cos_f + sin_f * sin_f).sqrt();
    if coeff.abs() <= f64::EPSILON {
        return (0.0, 0.0);
    }
    let mut angle = (cos_f / coeff).abs().acos();
    if sin_f > 0.0 {
        if cos_f.abs() <= f64::EPSILON {
            angle = PI / 2.0;
        } else if cos_f < 0.0 {
            angle = PI - angle;
        }
    } else if sin_f.abs() <= f64::EPSILON {
        if cos_f < 0.0 {
            angle = PI;
        }
    }
    if sin_f < 0.0 {
        if cos_f > 0.0 {
            angle = 2.0 * PI - angle;
        } else if cos_f.abs() <= f64::EPSILON {
            angle = 3.0 * PI / 2.0;
        } else if cos_f < 0.0 {
            angle = PI + angle;
        }
    }
    (coeff, angle)
}

fn coeffs_from_cylinders(c1: &GpCylinder, c2: &GpCylinder) -> Option<Coeffs> {
    let mut a1 = c1.x_axis().direction().xyz().multiplied(-c1.radius());
    let mut a2 = c2.x_axis().direction().xyz().multiplied(c2.radius());
    let mut b1 = c1.y_axis().direction().xyz().multiplied(-c1.radius());
    let mut b2 = c2.y_axis().direction().xyz().multiplied(c2.radius());
    let mut c1v = *c1.axis().direction().xyz();
    let mut c2v = c2.axis().direction().xyz().reversed();
    let mut d = c2.location().coord.subtracted(&c1.location().coord);

    let delta1 = xyz_comp(&c1v, 0) * xyz_comp(&c2v, 1) - xyz_comp(&c1v, 1) * xyz_comp(&c2v, 0);
    let delta2 = xyz_comp(&c1v, 1) * xyz_comp(&c2v, 2) - xyz_comp(&c1v, 2) * xyz_comp(&c2v, 1);
    let delta3 = xyz_comp(&c1v, 0) * xyz_comp(&c2v, 2) - xyz_comp(&c1v, 2) * xyz_comp(&c2v, 0);
    let (abs1, abs2, abs3) = (delta1.abs(), delta2.abs(), delta3.abs());
    let (couple, det) = if abs1 >= abs2 {
        if abs3 > abs1 {
            (13, delta3)
        } else {
            (12, delta1)
        }
    } else if abs3 > abs2 {
        (13, delta3)
    } else {
        (23, delta2)
    };
    if det.abs() < ANGULAR {
        return None;
    }
    match couple {
        23 => {
            rotate_123(&mut a1);
            rotate_123(&mut a2);
            rotate_123(&mut b1);
            rotate_123(&mut b2);
            rotate_123(&mut c1v);
            rotate_123(&mut c2v);
            rotate_123(&mut d);
        }
        13 => {
            swap_yz(&mut a1);
            swap_yz(&mut a2);
            swap_yz(&mut b1);
            swap_yz(&mut b2);
            swap_yz(&mut c1v);
            swap_yz(&mut c2v);
            swap_yz(&mut d);
        }
        _ => {}
    }

    let k21 = (xyz_comp(&c2v, 1) * xyz_comp(&b2, 0) - xyz_comp(&c2v, 0) * xyz_comp(&b2, 1)) / det;
    let k11 = (xyz_comp(&c2v, 1) * xyz_comp(&b1, 0) - xyz_comp(&c2v, 0) * xyz_comp(&b1, 1)) / det;
    let l21 = (xyz_comp(&c2v, 1) * xyz_comp(&a2, 0) - xyz_comp(&c2v, 0) * xyz_comp(&a2, 1)) / det;
    let l11 = (xyz_comp(&c2v, 1) * xyz_comp(&a1, 0) - xyz_comp(&c2v, 0) * xyz_comp(&a1, 1)) / det;
    let m1 = (xyz_comp(&c2v, 1) * xyz_comp(&d, 0) - xyz_comp(&c2v, 0) * xyz_comp(&d, 1)) / det;

    let k22 = (xyz_comp(&c1v, 0) * xyz_comp(&b2, 1) - xyz_comp(&c1v, 1) * xyz_comp(&b2, 0)) / det;
    let k12 = (xyz_comp(&c1v, 0) * xyz_comp(&b1, 1) - xyz_comp(&c1v, 1) * xyz_comp(&b1, 0)) / det;
    let l22 = (xyz_comp(&c1v, 0) * xyz_comp(&a2, 1) - xyz_comp(&c1v, 1) * xyz_comp(&a2, 0)) / det;
    let l12 = (xyz_comp(&c1v, 0) * xyz_comp(&a1, 1) - xyz_comp(&c1v, 1) * xyz_comp(&a1, 0)) / det;
    let m2 = (xyz_comp(&c1v, 0) * xyz_comp(&d, 1) - xyz_comp(&c1v, 1) * xyz_comp(&d, 0)) / det;

    let a_a1 = xyz_comp(&c1v, 2) * k21 + xyz_comp(&c2v, 2) * k22 - xyz_comp(&b2, 2);
    let a_a2 = xyz_comp(&c1v, 2) * l21 + xyz_comp(&c2v, 2) * l22 - xyz_comp(&a2, 2);
    let a_b1 = xyz_comp(&b1, 2) - xyz_comp(&c1v, 2) * k11 - xyz_comp(&c2v, 2) * k12;
    let a_b2 = xyz_comp(&a1, 2) - xyz_comp(&c1v, 2) * l11 - xyz_comp(&c2v, 2) * l12;
    let mut m_c = xyz_comp(&d, 2) - xyz_comp(&c1v, 2) * m1 - xyz_comp(&c2v, 2) * m2;

    let (m_b, fi1) = short_cos_form(a_b2, a_b1);
    let (a_a, fi2) = short_cos_form(a_a2, a_a1);
    if a_a.abs() < NUL_VALUE {
        return None;
    }
    let m_b = m_b / a_a;
    m_c /= a_a;

    Some(Coeffs {
        vec_a1: a1,
        vec_a2: a2,
        vec_b1: b1,
        vec_b2: b2,
        vec_c1: c1v,
        vec_c2: c2v,
        vec_d: d,
        k21,
        k11,
        l21,
        l11,
        m1,
        k22,
        k12,
        l22,
        l12,
        m2,
        b: m_b,
        c: m_c,
        fi1,
        fi2,
    })
}

fn in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    let period = u_last - u_first;
    if period.abs() < f64::EPSILON {
        return u;
    }
    let mut x = u;
    while x < u_first {
        x += period;
    }
    while x >= u_last {
        x -= period;
    }
    x
}

/// `InscribePoint`. Writes the wrapped parameter; returns whether it lies in `[uf, ul]`.
pub(crate) fn inscribe_point(
    uf: f64,
    ul: f64,
    u: &mut f64,
    tol2d: f64,
    period: f64,
    force: bool,
) -> bool {
    if Precision::is_infinite(*u) {
        return false;
    }
    if uf - *u <= tol2d && *u - ul <= tol2d {
        if force {
            let mut tmp = *u + period;
            if uf - tmp <= tol2d && tmp - ul <= tol2d {
                *u = tmp;
                return true;
            }
            tmp = *u - period;
            if uf - tmp <= tol2d && tmp - ul <= tol2d {
                *u = tmp;
            }
        }
        return true;
    }
    let a_uf = uf - tol2d;
    let a_ul = a_uf + period;
    *u = in_period(*u, a_uf, a_ul);
    uf - *u <= tol2d && *u - ul <= tol2d
}

/// `CylCylComputeParameters` (U2). Optional `delta` is the acos error estimate.
pub(crate) fn compute_u2(u1: f64, wl: i32, c: &Coeffs, delta: Option<&mut f64>) -> Option<f64> {
    if wl < 0 || wl > 1 {
        return None;
    }
    let sign = if wl == 0 { 1.0 } else { -1.0 };
    let tol0 = (10.0 * f64::EPSILON * c.b).min(NUL_VALUE);
    let tol = 1.0 - tol0;
    let mut arg = c.b * (u1 - c.fi1).cos() + c.c;
    if arg >= tol {
        if let Some(d) = delta {
            *d = 0.0;
        }
        arg = 1.0;
    } else if arg <= -tol {
        if let Some(d) = delta {
            *d = 0.0;
        }
        arg = -1.0;
    } else if let Some(d) = delta {
        let dd = (1.0 - arg).min(1.0 + arg);
        // `IntPatch_ImpImpIntersection.cxx:5391`:
        // `Standard_DivideByZero_Raise_if((aDelta * aDelta < RealSmall()) || (aDelta >= 2.0), ...)`.
        if dd * dd < REAL_SMALL || dd >= 2.0 {
            return None;
        }
        *d = tol0 / (dd * (2.0 - dd)).sqrt();
    }
    Some(c.fi2 + sign * arg.acos())
}

pub(crate) fn compute_v(u1: f64, u2: f64, c: &Coeffs) -> (f64, f64) {
    let v1 = c.k21 * u2.sin() + c.k11 * u1.sin() + c.l21 * u2.cos() + c.l11 * u1.cos() + c.m1;
    let v2 = c.k22 * u2.sin() + c.k12 * u1.sin() + c.l22 * u2.cos() + c.l12 * u1.cos() + c.m2;
    (v1, v2)
}

pub(crate) fn compute_params(u1: f64, wl: i32, c: &Coeffs) -> Option<(f64, f64, f64)> {
    let u2 = compute_u2(u1, wl, c, None)?;
    let (v1, v2) = compute_v(u1, u2, c);
    Some((u2, v1, v2))
}

fn extrema_line_line(c1: &GpAx1, c2: &GpAx1, cos_a: f64, sq_sin_a: f64) -> (f64, f64) {
    let l1l2 = c2.location().coord.subtracted(&c1.location().coord);
    let d1l = c1.direction().xyz().dot(&l1l2);
    let d2l = c2.direction().xyz().dot(&l1l2);
    let par1 = (d1l - cos_a * d2l) / sq_sin_a;
    let par2 = (cos_a * d1l - d2l) / sq_sin_a;
    (par1, par2)
}

fn boundary_estimation(
    c1: &GpCylinder,
    c2: &GpCylinder,
    uv1: CylUv,
    uv2: CylUv,
) -> Option<(f64, f64)> {
    let ax1 = c1.axis();
    let ax2 = c2.axis();
    let d1 = *ax1.direction();
    let d2 = *ax2.direction();
    let r1 = c1.radius();
    let r2 = c2.radius();
    let cos_a = d1.dot(&d2);
    let sq_sin = d1.xyz().cross_square_magnitude(d2.xyz());
    if sq_sin < ANGULAR * ANGULAR {
        return None;
    }
    let sin_a = sq_sin.sqrt();
    let abs_cos = cos_a.abs();
    let hdv1 = (r1 * abs_cos + r2) / sin_a;
    let hdv2 = (r2 * abs_cos + r1) / sin_a;
    let (v01, v02) = extrema_line_line(&ax1, &ax2, cos_a, sq_sin);
    let mut v1a = v01 - hdv1 - CONFUSION;
    let mut v1b = v01 + hdv1 + CONFUSION;
    let mut v2a = v02 - hdv2 - CONFUSION;
    let mut v2b = v02 + hdv2 + CONFUSION;
    if uv1.v0.is_finite() {
        v1a = v1a.max(uv1.v0);
    }
    if uv1.v1.is_finite() {
        v1b = v1b.min(uv1.v1);
    }
    if uv2.v0.is_finite() {
        v2a = v2a.max(uv2.v0);
    }
    if uv2.v1.is_finite() {
        v2b = v2b.min(uv2.v1);
    }
    if v1a > v1b || v2a > v2b {
        return None;
    }
    Some((v1b - v1a, v2b - v2a))
}

fn line_line_distance(a: &GpAx1, b: &GpAx1) -> f64 {
    use occt_core::gp::GpVec;
    let d1 = GpVec::from_xyz(a.direction().xyz());
    let d2 = GpVec::from_xyz(b.direction().xyz());
    let n = d1.crossed(&d2);
    let nm = n.magnitude();
    let w = GpVec::from_pnts(a.location(), b.location());
    if nm < ANGULAR {
        w.subtracted(&d1.multiplied_scalar(w.dot(&d1))).magnitude()
    } else {
        w.dot(&n).abs() / nm
    }
}

/// OCCT `isGoodIntersection` block in `CyCyNoGeometric`.
fn good_intersection(c1: &GpCylinder, c2: &GpCylinder) -> Option<f64> {
    use std::f64::consts::{FRAC_PI_2, PI};
    let to_much = 3.0;
    let crit_ang = PI / 18.0;
    let r1 = c1.radius();
    let r2 = c2.radius();
    let (rmax, rmin) = if r1 > to_much * r2 {
        (r1, r2)
    } else if r2 > to_much * r1 {
        (r2, r1)
    } else {
        return None;
    };
    let ax1 = c1.axis();
    let ax2 = c2.axis();
    if (FRAC_PI_2 - ax1.direction().angle(ax2.direction())).abs() > crit_ang {
        return None;
    }
    if line_line_distance(&ax1, &ax2) > rmax / 2.0 {
        return None;
    }
    let defl = 0.001;
    let mut nb_p = 3;
    if rmin * defl > 1.0e-3 {
        let ang = 2.0 * (1.0 - defl).acos();
        nb_p = (2.0 * PI / ang) as i32 + 1;
    }
    Some(PI / (nb_p as f64 - 1.0))
}

fn surf_dom(uv1: CylUv, uv2: CylUv, tol3d: f64, tol2d: f64) -> SurfDom {
    SurfDom {
        u1f: uv1.u0,
        u1l: uv1.u1,
        u2f: uv2.u0,
        u2l: uv2.u1,
        v1f: uv1.v0,
        v1l: uv1.v1,
        v2f: uv2.v0,
        v2l: uv2.v1,
        period: PERIOD,
        tol2d,
        tol3d,
    }
}

/// `CyCyNoGeometric` plus `IntCyCy` range-sum reverse.
pub(crate) fn cy_cy_no_geometric(
    c1: &GpCylinder,
    c2: &GpCylinder,
    uv1: CylUv,
    uv2: CylUv,
    tol3d: f64,
    tol2d: f64,
) -> Result<PairOutcome, bool> {
    let Some(coeffs1) = coeffs_from_cylinders(c1, c2) else {
        return Err(false);
    };
    let Some(coeffs2) = coeffs_from_cylinders(c2, c1) else {
        return Err(false);
    };
    let Some(r1) = bounds::boundaries_computing(&coeffs1, PERIOD) else {
        return Ok(PairOutcome::Empty);
    };
    let Some(r2) = bounds::boundaries_computing(&coeffs2, PERIOD) else {
        return Ok(PairOutcome::Empty);
    };
    let sum1 = bounds::sum_inscribed_u(r1, uv1.u0, uv1.u1, tol2d);
    let sum2 = bounds::sum_inscribed_u(r2, uv2.u0, uv2.u1, tol2d);
    if sum2 > sum1 {
        walk_one(c2, c1, uv2, uv1, coeffs2, r2, true, tol3d, tol2d)
    } else {
        walk_one(c1, c2, uv1, uv2, coeffs1, r1, false, tol3d, tol2d)
    }
}

fn walk_one(
    c1: &GpCylinder,
    c2: &GpCylinder,
    uv1: CylUv,
    uv2: CylUv,
    coeffs: Coeffs,
    mut ranges: [bounds::URange; 2],
    reversed: bool,
    tol3d: f64,
    tol2d: f64,
) -> Result<PairOutcome, bool> {
    let Some((dv1, dv2)) = boundary_estimation(c1, c2, uv1, uv2) else {
        return Ok(PairOutcome::Empty);
    };
    if dv1 > 1.0e5 * c1.radius() || dv2 > 1.0e5 * c2.radius() {
        return Err(true);
    }
    let good_du = good_intersection(c1, c2);
    let is_good = good_du.is_some();
    let (nb_max, nb_min, du) = if let Some(optdu) = good_du {
        (200, 50, optdu)
    } else {
        (1000, 200, PERIOD / 1000.0)
    };
    let u1f = uv1.u0;
    let u1l = uv1.u1;
    let nb_pts = (((u1l - u1f) / du).floor() as i32 + 1).min((20.0 * c1.radius()) as i32);
    let nb_points = nb_pts.clamp(nb_min, nb_max);
    let step_min = tol2d.max(occt_core::precision::PCONFUSION);
    let step_max = if u1l - u1f > PI / 100.0 {
        (u1l - u1f) / nb_points as f64
    } else {
        u1l - u1f
    };
    let dom = surf_dom(uv1, uv2, tol3d, tol2d);
    for r in &mut ranges {
        if !r.is_void() {
            let _ = bounds::inscribe_interval(dom.u1f, dom.u1l, r, tol2d, PERIOD);
        }
    }
    bounds::merge_ranges(&mut ranges);
    walk::cy_cy_walk(
        c1,
        c2,
        coeffs,
        reversed,
        dom,
        ranges,
        dv1,
        dv2,
        is_good,
        step_min,
        step_max,
        nb_points,
        nb_min,
        nb_max,
    )
}
