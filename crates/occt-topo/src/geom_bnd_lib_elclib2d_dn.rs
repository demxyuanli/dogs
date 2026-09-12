//! 2D `ElCLib` N-th derivatives (`LineDN` / `CircleDN` / `EllipseDN` /
//! `HyperbolaDN` / `ParabolaDN`).
//!
//! Source: `ElCLib.cxx` 1049–1188. Circle/ellipse cycle every 4 derivatives;
//! hyperbola uses odd/even sinh/cosh; parabola is zero for `N > 2`.

use occt_core::gp::{
    GpAx2d, GpAx22d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpVec2d, GpXY,
};
use occt_core::precision::RESOLUTION;

fn xdir(pos: &GpAx22d) -> GpXY {
    GpXY::new(pos.x_direction().x(), pos.x_direction().y())
}

fn ydir(pos: &GpAx22d) -> GpXY {
    GpXY::new(pos.y_direction().x(), pos.y_direction().y())
}

/// `ElCLib::LineDN(U, Pos, N)` — `N == 1` returns the direction, else zero.
pub fn line_dn_ax(_u: f64, pos: &GpAx2d, n: i32) -> GpVec2d {
    if n == 1 {
        GpVec2d::from_dir2d(pos.direction())
    } else {
        GpVec2d::zero()
    }
}

/// `ElCLib::CircleDN(U, Pos, Radius, N)`.
///
/// * `N == 1` or `(N-1) % 4 == 0`: `(-R sin U, R cos U)` in the (Xd, Yd) basis
/// * `(N+2) % 4 == 0`: `(-R cos U, -R sin U)`
/// * `(N+1) % 4 == 0`: `(R sin U, -R cos U)`
/// * `N % 4 == 0`: `(R cos U, R sin U)`
pub fn circle_dn_ax(u: f64, pos: &GpAx22d, radius: f64, n: i32) -> GpVec2d {
    let (mut xc, mut yc) = (0.0, 0.0);
    if n == 1 {
        xc = radius * (-u.sin());
        yc = radius * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = radius * (-u.cos());
        yc = radius * (-u.sin());
    } else if (n + 1) % 4 == 0 {
        xc = radius * u.sin();
        yc = radius * (-u.cos());
    } else if n % 4 == 0 {
        xc = radius * u.cos();
        yc = radius * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = radius * (-u.sin());
        yc = radius * u.cos();
    }
    let xd = xdir(pos);
    let yd = ydir(pos);
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc, &xd, yc, &yd);
    GpVec2d::from_xy(vxy)
}

/// `ElCLib::EllipseDN(U, Pos, Major, Minor, N)` — same 4-cycle as circle,
/// with Major/Minor replacing Radius on the X/Y coefficients.
pub fn ellipse_dn_ax(u: f64, pos: &GpAx22d, major: f64, minor: f64, n: i32) -> GpVec2d {
    let (mut xc, mut yc) = (0.0, 0.0);
    if n == 1 {
        xc = major * (-u.sin());
        yc = minor * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = major * (-u.cos());
        yc = minor * (-u.sin());
    } else if (n + 1) % 4 == 0 {
        xc = major * u.sin();
        yc = minor * (-u.cos());
    } else if n % 4 == 0 {
        xc = major * u.cos();
        yc = minor * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = major * (-u.sin());
        yc = minor * u.cos();
    }
    let xd = xdir(pos);
    let yd = ydir(pos);
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc, &xd, yc, &yd);
    GpVec2d::from_xy(vxy)
}

/// `ElCLib::HyperbolaDN(U, Pos, Major, Minor, N)`.
///
/// Odd `N`: `(Major sinh U, Minor cosh U)`; even: `(Major cosh U, Minor sinh U)`.
pub fn hyperbola_dn_ax(u: f64, pos: &GpAx22d, major: f64, minor: f64, n: i32) -> GpVec2d {
    let (xc, yc) = if n % 2 != 0 {
        (major * u.sinh(), minor * u.cosh())
    } else {
        (major * u.cosh(), minor * u.sinh())
    };
    let xd = xdir(pos);
    let yd = ydir(pos);
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc, &xd, yc, &yd);
    GpVec2d::from_xy(vxy)
}

/// `ElCLib::ParabolaDN(U, Pos, Focal, N)`.
///
/// `N > 2` or `N <= 0` → zero. `N == 1`: D1 tangent. `N == 2`: `Xd / (2 Focal)`
/// (or zero when degenerate).
pub fn parabola_dn_ax(u: f64, pos: &GpAx22d, focal: f64, n: i32) -> GpVec2d {
    if n > 2 || n <= 0 {
        return GpVec2d::zero();
    }
    let xd = xdir(pos);
    if n == 1 {
        if focal.abs() <= RESOLUTION {
            return GpVec2d::from_xy(xd);
        }
        let yd = ydir(pos);
        let mut vxy = GpXY::default();
        vxy.set_linear_form_add_scaled(u / (2.0 * focal), &xd, &yd);
        return GpVec2d::from_xy(vxy);
    }
    if focal.abs() <= RESOLUTION {
        return GpVec2d::zero();
    }
    GpVec2d::from_xy(xd.multiplied_scalar(1.0 / (2.0 * focal)))
}

/// `ElCLib::DN(U, gp_Lin2d, N)`.
pub fn line_dn(u: f64, lin: &GpLin2d, n: i32) -> GpVec2d {
    line_dn_ax(u, lin.position(), n)
}

/// `ElCLib::DN(U, gp_Circ2d, N)`.
pub fn circ_dn(u: f64, circ: &GpCirc2d, n: i32) -> GpVec2d {
    circle_dn_ax(u, circ.position(), circ.radius(), n)
}

/// `ElCLib::DN(U, gp_Elips2d, N)`.
pub fn elips_dn(u: f64, el: &GpElips2d, n: i32) -> GpVec2d {
    ellipse_dn_ax(u, &el.pos, el.major_radius, el.minor_radius, n)
}

/// `ElCLib::DN(U, gp_Hypr2d, N)`.
pub fn hypr_dn(u: f64, hy: &GpHypr2d, n: i32) -> GpVec2d {
    hyperbola_dn_ax(u, &hy.pos, hy.major_radius, hy.minor_radius, n)
}

/// `ElCLib::DN(U, gp_Parab2d, N)`.
pub fn parab_dn(u: f64, pa: &GpParab2d, n: i32) -> GpVec2d {
    parabola_dn_ax(u, &pa.pos, pa.focal, n)
}
