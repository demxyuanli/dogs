//! 2D `ElCLib` second and third derivatives.
//!
//! Source: `ElCLib.cxx` `CircleD2` (694), `EllipseD2` (720), `HyperbolaD2`
//! (750), `ParabolaD2` (779), `CircleD3` (809), `EllipseD3` (843),
//! `HyperbolaD3` (878). The Ax22d overloads are the ones
//! `GeomBndLib_*2d` extrema would use if they asked for D2 (circle/ellipse
//! UV boxes currently only need Value + InPeriod; D2 is the same evaluator).

use occt_core::gp::{
    GpAx22d, GpCirc2d, GpElips2d, GpHypr2d, GpParab2d, GpPnt2d, GpVec2d, GpXY,
};
use occt_core::precision::RESOLUTION;

fn xdir(pos: &GpAx22d) -> GpXY {
    GpXY::new(pos.x_direction().x(), pos.x_direction().y())
}

fn ydir(pos: &GpAx22d) -> GpXY {
    GpXY::new(pos.y_direction().x(), pos.y_direction().y())
}

fn loc(pos: &GpAx22d) -> GpXY {
    *pos.location().xy()
}

fn pnt_from_xy(xy: GpXY) -> GpPnt2d {
    GpPnt2d::from_xy(xy)
}

fn vec_from_xy(xy: GpXY) -> GpVec2d {
    GpVec2d::from_xy(xy)
}

/// `ElCLib::CircleD2(U, Pos, Radius, P, V1, V2)`.
///
/// `V2 = -(R cos U * Xd + R sin U * Yd)`, `P = -V2 + Location`,
/// `V1 = -R sin U * Xd + R cos U * Yd`.
pub fn circle_d2_ax(
    u: f64,
    pos: &GpAx22d,
    radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let yd = ydir(pos);
    let xc = radius * u.cos();
    let yc = radius * u.sin();
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc, &xd, yc, &yd);
    let mut v2 = vec_from_xy(vxy);
    v2.reverse();
    vxy.add(&loc(pos));
    let p = pnt_from_xy(vxy);
    vxy.set_linear_form_2(-yc, &xd, xc, &yd);
    let v1 = vec_from_xy(vxy);
    (p, v1, v2)
}

/// `ElCLib::EllipseD2(U, Pos, Major, Minor, P, V1, V2)`.
pub fn ellipse_d2_ax(
    u: f64,
    pos: &GpAx22d,
    major: f64,
    minor: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let yd = ydir(pos);
    let xc = u.cos();
    let yc = u.sin();
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc * major, &xd, yc * minor, &yd);
    let mut v2 = vec_from_xy(vxy);
    v2.reverse();
    vxy.add(&loc(pos));
    let p = pnt_from_xy(vxy);
    vxy.set_linear_form_2(-yc * major, &xd, xc * minor, &yd);
    let v1 = vec_from_xy(vxy);
    (p, v1, v2)
}

/// `ElCLib::HyperbolaD2(U, Pos, Major, Minor, P, V1, V2)`.
///
/// Unlike circle/ellipse, V2 is not reversed: `V2 = Major cosh U * Xd +
/// Minor sinh U * Yd`.
pub fn hyperbola_d2_ax(
    u: f64,
    pos: &GpAx22d,
    major: f64,
    minor: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let yd = ydir(pos);
    let xc = u.cosh();
    let yc = u.sinh();
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc * major, &xd, yc * minor, &yd);
    let v2 = vec_from_xy(vxy);
    vxy.add(&loc(pos));
    let p = pnt_from_xy(vxy);
    vxy.set_linear_form_2(yc * major, &xd, xc * minor, &yd);
    let v1 = vec_from_xy(vxy);
    (p, v1, v2)
}

/// `ElCLib::ParabolaD2(U, Pos, Focal, P, V1, V2)`.
///
/// Degenerate (`|Focal| <= gp::Resolution`): `V2 = 0`, `V1 = Xd`,
/// `P = Location + U * Xd`. Otherwise `V2 = Xd / (2 Focal)`,
/// `V1 = U * V2 + Yd`, `P = (U^2 / (4 Focal)) Xd + U Yd + Location`.
pub fn parabola_d2_ax(u: f64, pos: &GpAx22d, focal: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let mut vxy = GpXY::default();
    if focal.abs() <= RESOLUTION {
        let v2 = GpVec2d::new(0.0, 0.0);
        let v1 = vec_from_xy(xd);
        vxy.set_linear_form_add_scaled(u, &xd, &loc(pos));
        return (pnt_from_xy(vxy), v1, v2);
    }
    let yd = ydir(pos);
    let v2_xy = xd.multiplied_scalar(1.0 / (2.0 * focal));
    let v2 = vec_from_xy(v2_xy);
    let mut vxy = GpXY::default();
    vxy.set_linear_form_add_scaled(u, &v2_xy, &yd);
    let v1 = vec_from_xy(vxy);
    vxy.set_linear_form_2(u * u / (4.0 * focal), &xd, u, &yd);
    vxy.add(&loc(pos));
    (pnt_from_xy(vxy), v1, v2)
}

/// `ElCLib::CircleD3(U, Pos, Radius, P, V1, V2, V3)`.
///
/// `V3 = -V1`.
pub fn circle_d3_ax(
    u: f64,
    pos: &GpAx22d,
    radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let (p, v1, v2) = circle_d2_ax(u, pos, radius);
    let v3 = v1.reversed();
    (p, v1, v2, v3)
}

/// `ElCLib::EllipseD3(U, Pos, Major, Minor, P, V1, V2, V3)`.
///
/// `V3 = -V1`.
pub fn ellipse_d3_ax(
    u: f64,
    pos: &GpAx22d,
    major: f64,
    minor: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let yd = ydir(pos);
    let xc = u.cos();
    let yc = u.sin();
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc * major, &xd, yc * minor, &yd);
    let mut v2 = vec_from_xy(vxy);
    v2.reverse();
    vxy.add(&loc(pos));
    let p = pnt_from_xy(vxy);
    vxy.set_linear_form_2(-yc * major, &xd, xc * minor, &yd);
    let v1 = vec_from_xy(vxy);
    let v3 = v1.reversed();
    (p, v1, v2, v3)
}

/// `ElCLib::HyperbolaD3(U, Pos, Major, Minor, P, V1, V2, V3)`.
///
/// `V3 = V1` (not reversed — `ElCLib.cxx:906`).
pub fn hyperbola_d3_ax(
    u: f64,
    pos: &GpAx22d,
    major: f64,
    minor: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let xd = xdir(pos);
    let yd = ydir(pos);
    let xc = u.cosh();
    let yc = u.sinh();
    let mut vxy = GpXY::default();
    vxy.set_linear_form_2(xc * major, &xd, yc * minor, &yd);
    let v2 = vec_from_xy(vxy);
    vxy.add(&loc(pos));
    let p = pnt_from_xy(vxy);
    vxy.set_linear_form_2(yc * major, &xd, xc * minor, &yd);
    let v1 = vec_from_xy(vxy);
    let v3 = v1;
    (p, v1, v2, v3)
}

/// `ElCLib::D2(U, gp_Circ2d)`.
pub fn circ_d2(u: f64, circ: &GpCirc2d) -> (GpPnt2d, GpVec2d, GpVec2d) {
    circle_d2_ax(u, circ.position(), circ.radius())
}

/// `ElCLib::D3(U, gp_Circ2d)`.
pub fn circ_d3(u: f64, circ: &GpCirc2d) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    circle_d3_ax(u, circ.position(), circ.radius())
}

/// `ElCLib::D2(U, gp_Elips2d)`.
pub fn elips_d2(u: f64, el: &GpElips2d) -> (GpPnt2d, GpVec2d, GpVec2d) {
    ellipse_d2_ax(u, &el.pos, el.major_radius, el.minor_radius)
}

/// `ElCLib::D3(U, gp_Elips2d)`.
pub fn elips_d3(u: f64, el: &GpElips2d) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    ellipse_d3_ax(u, &el.pos, el.major_radius, el.minor_radius)
}

/// `ElCLib::D2(U, gp_Hypr2d)`.
pub fn hypr_d2(u: f64, hy: &GpHypr2d) -> (GpPnt2d, GpVec2d, GpVec2d) {
    hyperbola_d2_ax(u, &hy.pos, hy.major_radius, hy.minor_radius)
}

/// `ElCLib::D3(U, gp_Hypr2d)`.
pub fn hypr_d3(u: f64, hy: &GpHypr2d) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    hyperbola_d3_ax(u, &hy.pos, hy.major_radius, hy.minor_radius)
}

/// `ElCLib::D2(U, gp_Parab2d)`.
pub fn parab_d2(u: f64, pa: &GpParab2d) -> (GpPnt2d, GpVec2d, GpVec2d) {
    parabola_d2_ax(u, &pa.pos, pa.focal)
}
