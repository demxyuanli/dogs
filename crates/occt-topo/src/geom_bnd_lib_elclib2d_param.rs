//! 2D `ElCLib` inverse parameterisation and angle wrap.
//!
//! Source: `ElCLib.cxx` `normalizeAngle` (56), `LineParameter` (1276),
//! `CircleParameter` (1285), `EllipseParameter` (1295),
//! `HyperbolaParameter` (1315), `ParabolaParameter` (1331).

use occt_core::gp::{
    GpAx2d, GpAx22d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d, GpXY,
};
use occt_core::precision::{COMPUTATIONAL, RESOLUTION};

const PIPI: f64 = 2.0 * std::f64::consts::PI;
const NEGATIVE_RESOLUTION: f64 = -COMPUTATIONAL;

/// `normalizeAngle` — wrap into `[0, 2*PI]`, keep the closing seam at `2*PI`.
pub fn normalize_angle(the_angle: &mut f64) {
    while *the_angle < NEGATIVE_RESOLUTION {
        *the_angle += PIPI;
    }
    while *the_angle > PIPI * (1.0 + RESOLUTION) {
        *the_angle -= PIPI;
    }
    if *the_angle < 0.0 {
        *the_angle = 0.0;
    }
}

/// Signed angle of `gp_Dir2d` vs `gp_Vec2d` (`gp_Dir2d::Angle`).
fn dir_angle_vec(d_x: f64, d_y: f64, v: &GpVec2d) -> f64 {
    let cross = d_x * v.y() - d_y * v.x();
    let dot = d_x * v.x() + d_y * v.y();
    cross.atan2(dot)
}

/// `ElCLib::LineParameter(gp_Ax2d, gp_Pnt2d)` — `(P - Loc) · Dir`.
pub fn line_parameter_ax(axis: &GpAx2d, p: &GpPnt2d) -> f64 {
    let mut coord = *p.xy();
    coord.subtract(axis.location().xy());
    let d = axis.direction();
    coord.dot(&GpXY::new(d.x(), d.y()))
}

/// `ElCLib::CircleParameter(gp_Ax22d, gp_Pnt2d)`.
///
/// Angle of `XDirection` with `(Location → P)`, flipped when
/// `XDirection ^ YDirection < 0`, then [`normalize_angle`].
pub fn circle_parameter_ax(pos: &GpAx22d, p: &GpPnt2d) -> f64 {
    let loc = pos.location();
    let v = GpVec2d::new(p.x() - loc.x(), p.y() - loc.y());
    let xd = pos.x_direction();
    let yd = pos.y_direction();
    let mut teta = dir_angle_vec(xd.x(), xd.y(), &v);
    let cross = xd.x() * yd.y() - xd.y() * yd.x();
    if cross < 0.0 {
        teta = -teta;
    }
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::EllipseParameter(gp_Ax22d, Major, Minor, P)`.
///
/// Project `OP` onto the axes, scale the Y component by `Major/Minor`, then
/// take the signed angle of that reconstructed vector vs `XDirection`.
pub fn ellipse_parameter_ax(pos: &GpAx22d, major: f64, minor: f64, p: &GpPnt2d) -> f64 {
    let mut op = *p.xy();
    op.subtract(pos.location().xy());
    let xaxis = GpXY::new(pos.x_direction().x(), pos.x_direction().y());
    let mut yaxis = GpXY::new(pos.y_direction().x(), pos.y_direction().y());
    let mut om = xaxis.multiplied_scalar(op.dot(&xaxis));
    let scale = if minor.abs() > RESOLUTION {
        major / minor
    } else {
        0.0
    };
    yaxis.multiply_scalar(op.dot(&yaxis) * scale);
    om.add(&yaxis);
    let mut teta = dir_angle_vec(xaxis.x(), xaxis.y(), &GpVec2d::from_xy(om));
    let xd = pos.x_direction();
    let yd = pos.y_direction();
    let cross = xd.x() * yd.y() - xd.y() * yd.x();
    if cross < 0.0 {
        teta = -teta;
    }
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::HyperbolaParameter(gp_Ax22d, Major, Minor, P)`.
///
/// `asinh( ((P-Loc) · Yd) / Minor )`. The unused Major argument is kept
/// to match the OCCT signature.
pub fn hyperbola_parameter_ax(pos: &GpAx22d, _major: f64, minor: f64, p: &GpPnt2d) -> f64 {
    let loc = pos.location();
    let v = GpVec2d::new(p.x() - loc.x(), p.y() - loc.y());
    let yd = pos.y_direction();
    let sht = v.dot(&GpVec2d::from_dir2d(yd)) / minor;
    sht.asinh()
}

/// `ElCLib::ParabolaParameter(gp_Ax22d, P)` — `(P-Loc) · Yd`.
pub fn parabola_parameter_ax(pos: &GpAx22d, p: &GpPnt2d) -> f64 {
    let loc = pos.location();
    let v = GpVec2d::new(p.x() - loc.x(), p.y() - loc.y());
    v.dot(&GpVec2d::from_dir2d(pos.y_direction()))
}

/// `ElCLib::Parameter(gp_Lin2d, P)`.
pub fn line_parameter(lin: &GpLin2d, p: &GpPnt2d) -> f64 {
    line_parameter_ax(lin.position(), p)
}

/// `ElCLib::Parameter(gp_Circ2d, P)`.
pub fn circ_parameter(circ: &GpCirc2d, p: &GpPnt2d) -> f64 {
    circle_parameter_ax(circ.position(), p)
}

/// `ElCLib::Parameter(gp_Elips2d, P)`.
pub fn elips_parameter(el: &GpElips2d, p: &GpPnt2d) -> f64 {
    ellipse_parameter_ax(&el.pos, el.major_radius, el.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Hypr2d, P)`.
pub fn hypr_parameter(hy: &GpHypr2d, p: &GpPnt2d) -> f64 {
    hyperbola_parameter_ax(&hy.pos, hy.major_radius, hy.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Parab2d, P)`.
pub fn parab_parameter(pa: &GpParab2d, p: &GpPnt2d) -> f64 {
    parabola_parameter_ax(&pa.pos, p)
}
