//! Elementary Curves Library. Source: `ElCLib.hxx`
use crate::gp::{
    GpAx2, GpAx22d, GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr, GpLin, GpLin2d,
    GpParab, GpPnt, GpPnt2d, GpVec, GpVec2d, GpXY,
};
use crate::precision::{COMPUTATIONAL, RESOLUTION};

// Helper: GpVec from scaled & added GpXyz refs
fn vec_add(a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpVec {
    let r = a.multiplied(sa).added(&b.multiplied(sb));
    GpVec::from_xyz(&r)
}

fn pt_add(org: &crate::gp::GpXyz, a: &crate::gp::GpXyz, sa: f64, b: &crate::gp::GpXyz, sb: f64) -> GpPnt {
    let r = org.added(&a.multiplied(sa)).added(&b.multiplied(sb));
    GpPnt::from_xyz(&r)
}

pub fn line_value(l: &GpLin, u: f64) -> GpPnt {
    GpPnt::from_xyz(&l.pos.loc.coord.added(&l.pos.vdir.xyz().multiplied(u)))
}
pub fn line_d1(l: &GpLin, u: f64) -> (GpPnt, GpVec) {
    (line_value(l, u), GpVec::from_xyz(l.pos.vdir.xyz()))
}
pub fn line_d2(l: &GpLin, u: f64) -> (GpPnt, GpVec, GpVec) {
    let (p, d1) = line_d1(l, u); (p, d1, GpVec::zero())
}

pub fn circle_value(c: &GpCirc, u: f64) -> GpPnt {
    let r = c.radius;
    pt_add(&c.pos.location().coord, c.pos.x_direction().xyz(), r*u.cos(), c.pos.y_direction().xyz(), r*u.sin())
}
pub fn circle_d1(c: &GpCirc, u: f64) -> (GpPnt, GpVec) {
    let r = c.radius;
    let p = circle_value(c, u);
    (p, vec_add(c.pos.x_direction().xyz(), -r*u.sin(), c.pos.y_direction().xyz(), r*u.cos()))
}
pub fn circle_d2(c: &GpCirc, u: f64) -> (GpPnt, GpVec, GpVec) {
    let r = c.radius;
    let p = circle_value(c, u);
    let d1 = vec_add(c.pos.x_direction().xyz(), -r*u.sin(), c.pos.y_direction().xyz(), r*u.cos());
    let d2 = vec_add(c.pos.x_direction().xyz(), -r*u.cos(), c.pos.y_direction().xyz(), -r*u.sin());
    (p, d1, d2)
}

/// `ElCLib::EllipseValue(U, gp_Ax2, Major, Minor)` (`cxx:176-189`):
/// `P = Loc + Major*cos(U)*XDir + Minor*sin(U)*YDir`.
///
/// The minor-axis term is `+Minor*sin(U)`; the previous `-Minor*sin(U)` mirrored
/// the parameterisation about the major axis. A full ellipse is unchanged (its
/// point set is symmetric), but a trimmed arc is not: the vertex parameters
/// coming out of `ShapeAnalysis_Curve::Project` land on the mirrored angles, so
/// the forward arc mandated by `StepToTopoDS_GeometricTool::UpdateParam3d`
/// (`GeometricTool.cxx:270-279`) became the complement of OCCT's arc and the
/// edge left its faces (ATU01038 faces 130/140/156/216/222).
pub fn ellipse_value(e: &GpElips, u: f64) -> GpPnt {
    let a = e.major_radius; let b = e.minor_radius;
    pt_add(&e.location().coord, e.pos.x_direction().xyz(), a*u.cos(), e.pos.y_direction().xyz(), b*u.sin())
}
/// `ElCLib::EllipseD1(U, gp_Ax2, Major, Minor)` (`cxx:256-274`):
/// `V1 = -Major*sin(U)*XDir + Minor*cos(U)*YDir`.
pub fn ellipse_d1(e: &GpElips, u: f64) -> (GpPnt, GpVec) {
    let a = e.major_radius; let b = e.minor_radius;
    let p = ellipse_value(e, u);
    (p, vec_add(e.pos.x_direction().xyz(), -a*u.sin(), e.pos.y_direction().xyz(), b*u.cos()))
}
/// `ElCLib::EllipseD2(U, gp_Ax2, Major, Minor)` (`cxx:352-374`):
/// `V2 = -Major*cos(U)*XDir - Minor*sin(U)*YDir`.
pub fn ellipse_d2(e: &GpElips, u: f64) -> (GpPnt, GpVec, GpVec) {
    let a = e.major_radius; let b = e.minor_radius;
    let p = ellipse_value(e, u);
    let d1 = vec_add(e.pos.x_direction().xyz(), -a*u.sin(), e.pos.y_direction().xyz(), b*u.cos());
    let d2 = vec_add(e.pos.x_direction().xyz(), -a*u.cos(), e.pos.y_direction().xyz(), -b*u.sin());
    (p, d1, d2)
}

pub fn hyperbola_value(h: &GpHypr, u: f64) -> GpPnt {
    let a = h.major_radius; let b = h.minor_radius;
    pt_add(&h.location().coord, h.pos.x_direction().xyz(), a*u.cosh(), h.pos.y_direction().xyz(), b*u.sinh())
}
pub fn hyperbola_d1(h: &GpHypr, u: f64) -> (GpPnt, GpVec) {
    let a = h.major_radius; let b = h.minor_radius;
    let p = hyperbola_value(h, u);
    (p, vec_add(h.pos.x_direction().xyz(), a*u.sinh(), h.pos.y_direction().xyz(), b*u.cosh()))
}

pub fn parabola_value(p: &GpParab, u: f64) -> GpPnt {
    let f = p.focal;
    pt_add(&p.location().coord, p.pos.x_direction().xyz(), u*u/(4.0*f), p.pos.y_direction().xyz(), u)
}
pub fn parabola_d1(p: &GpParab, u: f64) -> (GpPnt, GpVec) {
    let f = p.focal;
    let pt = parabola_value(p, u);
    (pt, vec_add(p.pos.x_direction().xyz(), u/(2.0*f), p.pos.y_direction().xyz(), 1.0))
}

// 2D
pub fn line2d_value(l: &GpLin2d, u: f64) -> GpPnt2d { GpPnt2d::new(l.pos.loc.x()+u*l.pos.vdir.x, l.pos.loc.y()+u*l.pos.vdir.y) }
/// `ElCLib::CircleValue(U, gp_Ax22d, Radius)` (`cxx:530-539`):
/// `O + R*(cos U * Xd + sin U * Yd)`.
pub fn circle2d_value(c: &GpCirc2d, u: f64) -> GpPnt2d {
    let r = c.radius;
    let cx = c.pos.point.x();
    let cy = c.pos.point.y();
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    GpPnt2d::new(a1 * xd.x + a2 * yd.x + cx, a1 * xd.y + a2 * yd.y + cy)
}

/// `ElCLib::CircleD1` 2d (`cxx:602-619`).
pub fn circle2d_d1(c: &GpCirc2d, u: f64) -> (GpPnt2d, GpVec2d) {
    let r = c.radius;
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    let xc = r * u.cos();
    let yc = r * u.sin();
    let p = GpPnt2d::new(
        xc * xd.x + yc * yd.x + c.pos.point.x(),
        xc * xd.y + yc * yd.y + c.pos.point.y(),
    );
    let v = GpVec2d::new(-yc * xd.x + xc * yd.x, -yc * xd.y + xc * yd.y);
    (p, v)
}

/// `ElCLib::CircleD2` 2d (`cxx:694-716`).
pub fn circle2d_d2(c: &GpCirc2d, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let (p, v1) = circle2d_d1(c, u);
    let r = c.radius;
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    let xc = r * u.cos();
    let yc = r * u.sin();
    let v2 = GpVec2d::new(-(xc * xd.x + yc * yd.x), -(xc * xd.y + yc * yd.y));
    (p, v1, v2)
}
/// `ElCLib::EllipseValue(U, gp_Ax22d, Major, Minor)` (`cxx:543-555`):
/// `P = Loc + Major*cos(U)*Xd + Minor*sin(U)*Yd` (same `+Minor*sin(U)` as the
/// 3D `gp_Ax2` overload, so the 3D and 2D parameterisations agree).
pub fn ellipse2d_value(e: &GpElips2d, u: f64) -> GpPnt2d {
    let a=e.major_radius; let b=e.minor_radius; let cx=e.pos.point.x(); let cy=e.pos.point.y(); let xd=e.pos.vxdir; let yd=e.pos.vydir;
    GpPnt2d::new(cx+a*u.cos()*xd.x+b*u.sin()*yd.x, cy+a*u.cos()*xd.y+b*u.sin()*yd.y)
}

const PIPI: f64 = 2.0 * std::f64::consts::PI;
const NEGATIVE_RESOLUTION: f64 = -COMPUTATIONAL;

/// `gp_Dir2d::Angle(gp_Dir2d)` (`gp_Dir2d.cxx:26-62`). Both arguments must be
/// unit vectors; the result is signed and lies in `[-PI, PI]`.
pub(super) fn dir2d_angle(a: &GpDir2d, b: &GpDir2d) -> f64 {
    let cosinus = a.dot(b);
    let sinus = a.crossed(b);
    const COS_45: f64 = std::f64::consts::FRAC_1_SQRT_2;
    if cosinus > -COS_45 && cosinus < COS_45 {
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        -std::f64::consts::PI - sinus.asin()
    }
}

/// `ElCLib::LineParameter(gp_Ax2d, gp_Pnt2d)` (`ElCLib.cxx:1276-1281`).
pub fn line2d_parameter(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let mut coord = *p.xy();
    coord.subtract(l.pos.loc.xy());
    coord.dot(&GpXY::new(l.pos.vdir.x, l.pos.vdir.y))
}

/// `ElCLib::CircleParameter(gp_Ax22d, gp_Pnt2d)` (`ElCLib.cxx:1285-1291`).
///
/// The OCCT body is `Pos.XDirection().Angle(gp_Vec2d(Pos.Location(), P))`, which
/// relies on the implicit `gp_Dir2d(const gp_Vec2d&)` conversion
/// (`gp_Dir2d.hxx:67-69,290`) that normalizes the vector. A null vector makes
/// OCCT raise `Standard_ConstructionError`; like the 3D `circle_parameter`
/// below, the degenerate case returns `0.0` instead of aborting.
pub fn circle2d_parameter(pos: &GpAx22d, p: &GpPnt2d) -> f64 {
    let vx = p.x() - pos.point.x();
    let vy = p.y() - pos.point.y();
    let Ok(vdir) = GpDir2d::new(vx, vy) else {
        return 0.0;
    };
    let mut teta = dir2d_angle(&pos.vxdir, &vdir);
    if pos.vxdir.crossed(&pos.vydir) < 0.0 {
        teta = -teta;
    }
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::normalizeAngle` (`ElCLib.cxx:56-72`): wrap into `[0, 2*PI]`, keep
/// the closing seam at exactly `2*PI`.
pub(super) fn normalize_angle(the_angle: &mut f64) {
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

/// `ElCLib::CircleParameter(gp_Ax2, gp_Pnt)` (`ElCLib.cxx:1199-1222`).
pub fn circle_parameter(pos: &GpAx2, p: &GpPnt) -> f64 {
    let a_vec = GpVec::from_pnts(&pos.location(), p);
    if a_vec.square_magnitude() < RESOLUTION {
        return 0.0;
    }
    let dir = pos.direction();
    let a_v_proj = dir.xyz().cross_crossed(a_vec.xyz(), dir.xyz());
    if a_v_proj.square_modulus() < RESOLUTION {
        return 0.0;
    }
    let Ok(proj_dir) = GpDir::from_xyz(&a_v_proj) else {
        return 0.0;
    };
    let mut teta = pos.x_direction().angle_with_ref(&proj_dir, &dir);
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::EllipseParameter(gp_Ax2, Major, Minor, gp_Pnt)` (`ElCLib.cxx:1226-1249`).
///
/// `gp_Vec(xaxis).AngleWithRef(gp_Vec(Om), gp_Vec(Pos.Direction()))` normalizes
/// all three vectors, so the `gp_Dir` form is used here; a degenerate `Om`
/// (which makes OCCT raise `Standard_ConstructionError`) returns `0.0`, matching
/// the convention of [`circle_parameter`] above.
pub fn ellipse_parameter(pos: &GpAx2, major_radius: f64, minor_radius: f64, p: &GpPnt) -> f64 {
    let op = p.xyz().subtracted(pos.location().xyz());
    let xaxis = *pos.x_direction().xyz();
    let yaxis = *pos.y_direction().xyz();
    let ny = op.dot(&yaxis);
    let nx = op.dot(&xaxis);

    if nx.abs() <= RESOLUTION && ny.abs() <= RESOLUTION {
        // The point P is on the axis of the ellipse.
        return 0.0;
    }

    let yaxis = yaxis.multiply_scalar(ny * (major_radius / minor_radius));
    let om = xaxis.multiplied(nx).added(&yaxis);

    let (Ok(dx), Ok(dom)) = (GpDir::from_xyz(&xaxis), GpDir::from_xyz(&om)) else {
        return 0.0;
    };
    let mut teta = dx.angle_with_ref(&dom, &pos.direction());
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::HyperbolaParameter(gp_Ax2, Major /*unused*/, Minor, gp_Pnt)`
/// (`ElCLib.cxx:1253-1265`).
pub fn hyperbola_parameter(pos: &GpAx2, _major_radius: f64, minor_radius: f64, p: &GpPnt) -> f64 {
    let sht = GpVec::from_pnts(&pos.location(), p)
        .xyz()
        .dot(&pos.y_direction().xyz())
        / minor_radius;
    sht.asinh()
}

/// `ElCLib::ParabolaParameter(gp_Ax2, gp_Pnt)` (`ElCLib.cxx:1269-1272`).
pub fn parabola_parameter(pos: &GpAx2, p: &GpPnt) -> f64 {
    GpVec::from_pnts(&pos.location(), p)
        .xyz()
        .dot(&pos.y_direction().xyz())
}

/// `ElCLib::Parameter(gp_Elips, gp_Pnt)` (`ElCLib.lxx:335-339`).
pub fn parameter_elips(e: &GpElips, p: &GpPnt) -> f64 {
    ellipse_parameter(&e.pos, e.major_radius, e.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Hypr, gp_Pnt)` (`ElCLib.lxx:341-345`).
pub fn parameter_hypr(h: &GpHypr, p: &GpPnt) -> f64 {
    hyperbola_parameter(&h.pos, h.major_radius, h.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Parab, gp_Pnt)` (`ElCLib.lxx:347-351`).
pub fn parameter_parab(prb: &GpParab, p: &GpPnt) -> f64 {
    parabola_parameter(&prb.pos, p)
}
