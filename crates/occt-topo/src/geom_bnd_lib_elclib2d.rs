//! 2D `ElCLib` evaluators used by `GeomBndLib_*2d`.
//!
//! Source: `ElCLib.hxx` (`Value` / `D1` / `InPeriod` for `gp_Lin2d`,
//! `gp_Circ2d`, `gp_Elips2d`, `gp_Hypr2d`, `gp_Parab2d`). PerformAreas UV
//! boxes call these through `GeomBndLib_Line2d` / `Circle2d` / `Ellipse2d` /
//! `Hyperbola2d` / `Parabola2d` rather than sampling.

use occt_core::gp::{
    GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d,
};
use occt_core::precision::{Precision, RESOLUTION};

/// `Epsilon(x)` used by `ElCLib::InPeriod` / `AdjustPeriodic`.
fn epsilon_of(x: f64) -> f64 {
    x.abs() * f64::EPSILON
}

/// `ElCLib::InPeriod` (`ElCLib.cxx:95`).
///
/// Infinite inputs are returned unchanged. A degenerate period
/// (`period < Epsilon(ULast)`) returns `u` unchanged. Otherwise
/// `max(UFirst, U + period * ceil((UFirst - U) / period))`.
pub fn in_period(the_u: f64, the_u_first: f64, the_u_last: f64) -> f64 {
    if Precision::is_infinite(the_u)
        || Precision::is_infinite(the_u_first)
        || Precision::is_infinite(the_u_last)
    {
        return the_u;
    }
    let a_period = the_u_last - the_u_first;
    if a_period < epsilon_of(the_u_last) {
        return the_u;
    }
    the_u_first.max(the_u + a_period * ((the_u_first - the_u) / a_period).ceil())
}

/// `ElCLib::AdjustPeriodic` (`ElCLib.cxx:115`).
///
/// Maps `[u1, u2]` into the periodic window `[u_first, u_last]` so that `u1`
/// sits in `[u_first, u_last)` (with a `preci` snap at the closing seam) and
/// `u2 - u1` is a positive length of at most one period.
pub fn adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if Precision::is_infinite(u_first) || Precision::is_infinite(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < epsilon_of(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

/// `ElCLib::Value(U, gp_Lin2d)` — `P = Location + U * Direction`.
pub fn line_value(u: f64, lin: &GpLin2d) -> GpPnt2d {
    let loc = lin.location();
    let d = lin.direction();
    GpPnt2d::new(loc.x() + u * d.x(), loc.y() + u * d.y())
}

/// `ElCLib::D1(U, gp_Lin2d)` — point plus constant tangent.
pub fn line_d1(u: f64, lin: &GpLin2d) -> (GpPnt2d, GpVec2d) {
    let p = line_value(u, lin);
    let d = lin.direction();
    (p, GpVec2d::new(d.x(), d.y()))
}

fn circ_axes(c: &GpCirc2d) -> (GpPnt2d, GpDir2d, GpDir2d, f64) {
    (
        c.location(),
        *c.position().x_direction(),
        *c.position().y_direction(),
        c.radius(),
    )
}

/// `ElCLib::Value(U, gp_Circ2d)` — `O + R*(cos U * Xd + sin U * Yd)`.
pub fn circ_value(u: f64, circ: &GpCirc2d) -> GpPnt2d {
    let (o, xd, yd, r) = circ_axes(circ);
    let cu = u.cos();
    let su = u.sin();
    GpPnt2d::new(
        o.x() + r * cu * xd.x() + r * su * yd.x(),
        o.y() + r * cu * xd.y() + r * su * yd.y(),
    )
}

/// `ElCLib::D1(U, gp_Circ2d)` — first derivative `R*(-sin U * Xd + cos U * Yd)`.
pub fn circ_d1(u: f64, circ: &GpCirc2d) -> (GpPnt2d, GpVec2d) {
    let (o, xd, yd, r) = circ_axes(circ);
    let cu = u.cos();
    let su = u.sin();
    let p = GpPnt2d::new(
        o.x() + r * cu * xd.x() + r * su * yd.x(),
        o.y() + r * cu * xd.y() + r * su * yd.y(),
    );
    let v = GpVec2d::new(
        r * (-su) * xd.x() + r * cu * yd.x(),
        r * (-su) * xd.y() + r * cu * yd.y(),
    );
    (p, v)
}

fn elips_axes(e: &GpElips2d) -> (GpPnt2d, GpDir2d, GpDir2d, f64, f64) {
    (
        e.pos.point,
        e.pos.vxdir,
        e.pos.vydir,
        e.major_radius,
        e.minor_radius,
    )
}

/// `ElCLib::Value(U, gp_Elips2d)` — `O + Major*cos U * Xd + Minor*sin U * Yd`.
pub fn elips_value(u: f64, el: &GpElips2d) -> GpPnt2d {
    let (o, xd, yd, a, b) = elips_axes(el);
    let cu = u.cos();
    let su = u.sin();
    GpPnt2d::new(
        o.x() + a * cu * xd.x() + b * su * yd.x(),
        o.y() + a * cu * xd.y() + b * su * yd.y(),
    )
}

/// `ElCLib::D1(U, gp_Elips2d)`.
pub fn elips_d1(u: f64, el: &GpElips2d) -> (GpPnt2d, GpVec2d) {
    let (o, xd, yd, a, b) = elips_axes(el);
    let cu = u.cos();
    let su = u.sin();
    let p = GpPnt2d::new(
        o.x() + a * cu * xd.x() + b * su * yd.x(),
        o.y() + a * cu * xd.y() + b * su * yd.y(),
    );
    let v = GpVec2d::new(
        -a * su * xd.x() + b * cu * yd.x(),
        -a * su * xd.y() + b * cu * yd.y(),
    );
    (p, v)
}

fn hypr_axes(h: &GpHypr2d) -> (GpPnt2d, GpDir2d, GpDir2d, f64, f64) {
    (
        h.pos.point,
        h.pos.vxdir,
        h.pos.vydir,
        h.major_radius,
        h.minor_radius,
    )
}

/// `ElCLib::Value(U, gp_Hypr2d)` — `O + Major*cosh U * Xd + Minor*sinh U * Yd`.
pub fn hypr_value(u: f64, hy: &GpHypr2d) -> GpPnt2d {
    let (o, xd, yd, a, b) = hypr_axes(hy);
    let ch = u.cosh();
    let sh = u.sinh();
    GpPnt2d::new(
        o.x() + a * ch * xd.x() + b * sh * yd.x(),
        o.y() + a * ch * xd.y() + b * sh * yd.y(),
    )
}

/// `ElCLib::D1(U, gp_Hypr2d)`.
pub fn hypr_d1(u: f64, hy: &GpHypr2d) -> (GpPnt2d, GpVec2d) {
    let (o, xd, yd, a, b) = hypr_axes(hy);
    let ch = u.cosh();
    let sh = u.sinh();
    let p = GpPnt2d::new(
        o.x() + a * ch * xd.x() + b * sh * yd.x(),
        o.y() + a * ch * xd.y() + b * sh * yd.y(),
    );
    let v = GpVec2d::new(
        a * sh * xd.x() + b * ch * yd.x(),
        a * sh * xd.y() + b * ch * yd.y(),
    );
    (p, v)
}

fn parab_axes(p: &GpParab2d) -> (GpPnt2d, GpDir2d, GpDir2d, f64) {
    (p.pos.point, p.pos.vxdir, p.pos.vydir, p.focal)
}

/// `ElCLib::Value(U, gp_Parab2d)` — `O + (U^2 / (4 f)) * Xd + U * Yd`.
pub fn parab_value(u: f64, pa: &GpParab2d) -> GpPnt2d {
    let (o, xd, yd, f) = parab_axes(pa);
    let x = if f.abs() > RESOLUTION {
        u * u / (4.0 * f)
    } else {
        0.0
    };
    GpPnt2d::new(
        o.x() + x * xd.x() + u * yd.x(),
        o.y() + x * xd.y() + u * yd.y(),
    )
}

/// `ElCLib::D1(U, gp_Parab2d)` — `(U / (2 f)) * Xd + Yd`.
pub fn parab_d1(u: f64, pa: &GpParab2d) -> (GpPnt2d, GpVec2d) {
    let (o, xd, yd, f) = parab_axes(pa);
    let p = parab_value(u, pa);
    let s = if f.abs() > RESOLUTION {
        u / (2.0 * f)
    } else {
        0.0
    };
    let v = GpVec2d::new(s * xd.x() + yd.x(), s * xd.y() + yd.y());
    let _ = o;
    (p, v)
}

/// Coordinate of a `gp_XY` / `gp_Dir2d` by 1-based index (`Coord(k)`).
pub fn coord2(x: f64, y: f64, k: i32) -> f64 {
    if k == 1 {
        x
    } else {
        y
    }
}
