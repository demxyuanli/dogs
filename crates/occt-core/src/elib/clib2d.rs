//! Elementary Curves Library, 2D overloads. Source: `ElCLib/ElCLib.cxx`.
//!
//! Every function here is the 2D overload of an `ElCLib` static member that
//! takes a `gp_Ax2d` / `gp_Ax22d` axis, or which only exists in 2D. The name
//! suffix `_ax2d` / `_ax22d` records which axis object OCCT passes, so the
//! `gp`-struct-taking overloads that already live in [`super::clib`]
//! (`line2d_value(&GpLin2d, u)`, `circle2d_d1(&GpCirc2d, u)`, ...) stay
//! unambiguous: different OCCT members, different Rust names.
//!
//! Line-number citations are against OCCT 8.0.0.
use super::clib::{dir2d_angle, normalize_angle};
use crate::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d, GpXY};
use crate::precision::{epsilon, INFINITE, REAL_SMALL};

/// `Precision::IsInfinite` (`Precision.hxx`).
fn is_infinite(r: f64) -> bool {
    r.abs() >= 0.5 * INFINITE
}

fn is_odd(n: i32) -> bool {
    n % 2 != 0
}

// ---------------------------------------------------------------------------
// Period helpers
// ---------------------------------------------------------------------------

/// `ElCLib::InPeriod` (`cxx:95-111`).
pub fn in_period(the_u: f64, the_u_first: f64, the_u_last: f64) -> f64 {
    if is_infinite(the_u) || is_infinite(the_u_first) || is_infinite(the_u_last) {
        return the_u;
    }
    let a_period = the_u_last - the_u_first;
    if a_period < epsilon(the_u_last) {
        return the_u;
    }
    the_u_first.max(the_u + a_period * ((the_u_first - the_u) / a_period).ceil())
}

/// `ElCLib::AdjustPeriodic` (`cxx:115-148`).
pub fn adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if is_infinite(u_first) || is_infinite(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < epsilon(u_last) {
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

// ---------------------------------------------------------------------------
// Value
// ---------------------------------------------------------------------------

/// `ElCLib::LineValue(U, gp_Ax2d)` (`cxx:519-523`).
pub fn line_value_ax2d(u: f64, pos: &GpAx2d) -> GpPnt2d {
    GpPnt2d::new(
        u * pos.vdir.x + pos.loc.x(),
        u * pos.vdir.y + pos.loc.y(),
    )
}

/// `ElCLib::CircleValue(U, gp_Ax22d, Radius)` (`cxx:530-539`).
pub fn circle_value_ax22d(u: f64, pos: &GpAx22d, radius: f64) -> GpPnt2d {
    let a1 = radius * u.cos();
    let a2 = radius * u.sin();
    GpPnt2d::new(
        a1 * pos.vxdir.x + a2 * pos.vydir.x + pos.point.x(),
        a1 * pos.vxdir.y + a2 * pos.vydir.y + pos.point.y(),
    )
}

/// `ElCLib::EllipseValue(U, gp_Ax22d, Major, Minor)` (`cxx:543-555`).
pub fn ellipse_value_ax22d(u: f64, pos: &GpAx22d, major_radius: f64, minor_radius: f64) -> GpPnt2d {
    let a1 = major_radius * u.cos();
    let a2 = minor_radius * u.sin();
    GpPnt2d::new(
        a1 * pos.vxdir.x + a2 * pos.vydir.x + pos.point.x(),
        a1 * pos.vxdir.y + a2 * pos.vydir.y + pos.point.y(),
    )
}

/// `ElCLib::HyperbolaValue(U, gp_Ax22d, Major, Minor)` (`cxx:559-572`).
pub fn hyperbola_value_ax22d(u: f64, pos: &GpAx22d, major_radius: f64, minor_radius: f64) -> GpPnt2d {
    let a1 = major_radius * u.cosh();
    let a2 = minor_radius * u.sinh();
    GpPnt2d::new(
        a1 * pos.vxdir.x + a2 * pos.vydir.x + pos.point.x(),
        a1 * pos.vxdir.y + a2 * pos.vydir.y + pos.point.y(),
    )
}

/// `ElCLib::ParabolaValue(U, gp_Ax22d, Focal)` (`cxx:575-589`).
pub fn parabola_value_ax22d(u: f64, pos: &GpAx22d, focal: f64) -> GpPnt2d {
    if focal.abs() <= REAL_SMALL {
        // `cxx:580-583`: degenerate parabola collapses onto its X axis.
        return GpPnt2d::new(
            u * pos.vxdir.x + pos.point.x(),
            u * pos.vxdir.y + pos.point.y(),
        );
    }
    let a1 = u * u / (4.0 * focal);
    GpPnt2d::new(
        a1 * pos.vxdir.x + u * pos.vydir.x + pos.point.x(),
        a1 * pos.vxdir.y + u * pos.vydir.y + pos.point.y(),
    )
}

// ---------------------------------------------------------------------------
// D1
// ---------------------------------------------------------------------------

/// `ElCLib::LineD1(U, gp_Ax2d)` (`cxx:592-598`).
pub fn line_d1_ax2d(u: f64, pos: &GpAx2d) -> (GpPnt2d, GpVec2d) {
    let coord = GpXY::new(pos.vdir.x, pos.vdir.y);
    let v1 = GpVec2d::from_xy(coord);
    let p = GpPnt2d::new(u * coord.x + pos.loc.x(), u * coord.y + pos.loc.y());
    (p, v1)
}

/// `ElCLib::CircleD1(U, gp_Ax22d, Radius)` (`cxx:602-619`).
pub fn circle_d1_ax22d(u: f64, pos: &GpAx22d, radius: f64) -> (GpPnt2d, GpVec2d) {
    let xc = radius * u.cos();
    let yc = radius * u.sin();
    let p = GpPnt2d::new(
        xc * pos.vxdir.x + yc * pos.vydir.x + pos.point.x(),
        xc * pos.vxdir.y + yc * pos.vydir.y + pos.point.y(),
    );
    let v1 = GpVec2d::new(
        -yc * pos.vxdir.x + xc * pos.vydir.x,
        -yc * pos.vxdir.y + xc * pos.vydir.y,
    );
    (p, v1)
}

/// `ElCLib::EllipseD1(U, gp_Ax22d, Major, Minor)` (`cxx:623-643`).
pub fn ellipse_d1_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d) {
    let xc = u.cos();
    let yc = u.sin();
    let p = GpPnt2d::new(
        xc * major_radius * pos.vxdir.x + yc * minor_radius * pos.vydir.x + pos.point.x(),
        xc * major_radius * pos.vxdir.y + yc * minor_radius * pos.vydir.y + pos.point.y(),
    );
    let v1 = GpVec2d::new(
        -yc * major_radius * pos.vxdir.x + xc * minor_radius * pos.vydir.x,
        -yc * major_radius * pos.vxdir.y + xc * minor_radius * pos.vydir.y,
    );
    (p, v1)
}

/// `ElCLib::HyperbolaD1(U, gp_Ax22d, Major, Minor)` (`cxx:646-666`).
pub fn hyperbola_d1_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d) {
    let xc = u.cosh();
    let yc = u.sinh();
    let p = GpPnt2d::new(
        xc * major_radius * pos.vxdir.x + yc * minor_radius * pos.vydir.x + pos.point.x(),
        xc * major_radius * pos.vxdir.y + yc * minor_radius * pos.vydir.y + pos.point.y(),
    );
    let v1 = GpVec2d::new(
        yc * major_radius * pos.vxdir.x + xc * minor_radius * pos.vydir.x,
        yc * major_radius * pos.vxdir.y + xc * minor_radius * pos.vydir.y,
    );
    (p, v1)
}

/// `ElCLib::ParabolaD1(U, gp_Ax22d, Focal)` (`cxx:669-691`).
pub fn parabola_d1_ax22d(u: f64, pos: &GpAx22d, focal: f64) -> (GpPnt2d, GpVec2d) {
    let xdir = GpXY::new(pos.vxdir.x, pos.vxdir.y);
    if focal.abs() <= REAL_SMALL {
        let v1 = GpVec2d::from_xy(xdir);
        let p = GpPnt2d::new(u * xdir.x + pos.point.x(), u * xdir.y + pos.point.y());
        (p, v1)
    } else {
        let ydir = GpXY::new(pos.vydir.x, pos.vydir.y);
        let v1 = GpVec2d::new(
            u / (2.0 * focal) * xdir.x + ydir.x,
            u / (2.0 * focal) * xdir.y + ydir.y,
        );
        let p = GpPnt2d::new(
            (u * u) / (4.0 * focal) * xdir.x + u * ydir.x + pos.point.x(),
            (u * u) / (4.0 * focal) * xdir.y + u * ydir.y + pos.point.y(),
        );
        (p, v1)
    }
}

// ---------------------------------------------------------------------------
// D2
// ---------------------------------------------------------------------------

/// `ElCLib::CircleD2(U, gp_Ax22d, Radius)` (`cxx:694-716`): `V2 = -O + P`.
pub fn circle_d2_ax22d(u: f64, pos: &GpAx22d, radius: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let (p, v1) = circle_d1_ax22d(u, pos, radius);
    let v2 = GpVec2d::new(
        -(radius * u.cos() * pos.vxdir.x + radius * u.sin() * pos.vydir.x),
        -(radius * u.cos() * pos.vxdir.y + radius * u.sin() * pos.vydir.y),
    );
    (p, v1, v2)
}

/// `ElCLib::EllipseD2(U, gp_Ax22d, Major, Minor)` (`cxx:720-747`): `V2` is the
/// negative of the point's offset from the center.
pub fn ellipse_d2_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xc = u.cos();
    let yc = u.sin();
    let v2_base = GpXY::new(
        xc * major_radius * pos.vxdir.x + yc * minor_radius * pos.vydir.x,
        xc * major_radius * pos.vxdir.y + yc * minor_radius * pos.vydir.y,
    );
    let v2 = GpVec2d::from_xy(v2_base.reversed());
    let p = GpPnt2d::new(v2_base.x + pos.point.x(), v2_base.y + pos.point.y());
    let v1 = GpVec2d::new(
        -yc * major_radius * pos.vxdir.x + xc * minor_radius * pos.vydir.x,
        -yc * major_radius * pos.vxdir.y + xc * minor_radius * pos.vydir.y,
    );
    (p, v1, v2)
}

/// `ElCLib::HyperbolaD2(U, gp_Ax22d, Major, Minor)` (`cxx:750-776`). `V2` is
/// *not* reversed here, unlike the circle and ellipse overloads.
pub fn hyperbola_d2_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xc = u.cosh();
    let yc = u.sinh();
    let v2_base = GpXY::new(
        xc * major_radius * pos.vxdir.x + yc * minor_radius * pos.vydir.x,
        xc * major_radius * pos.vxdir.y + yc * minor_radius * pos.vydir.y,
    );
    let v2 = GpVec2d::from_xy(v2_base);
    let p = GpPnt2d::new(v2_base.x + pos.point.x(), v2_base.y + pos.point.y());
    let v1 = GpVec2d::new(
        yc * major_radius * pos.vxdir.x + xc * minor_radius * pos.vydir.x,
        yc * major_radius * pos.vxdir.y + xc * minor_radius * pos.vydir.y,
    );
    (p, v1, v2)
}

/// `ElCLib::ParabolaD2(U, gp_Ax22d, Focal)` (`cxx:779-806`).
pub fn parabola_d2_ax22d(u: f64, pos: &GpAx22d, focal: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    let xdir = GpXY::new(pos.vxdir.x, pos.vxdir.y);
    if focal.abs() <= REAL_SMALL {
        let v2 = GpVec2d::new(0.0, 0.0);
        let v1 = GpVec2d::from_xy(xdir);
        let p = GpPnt2d::new(u * xdir.x + pos.point.x(), u * xdir.y + pos.point.y());
        (p, v1, v2)
    } else {
        let ydir = GpXY::new(pos.vydir.x, pos.vydir.y);
        let v2 = GpVec2d::new(xdir.x / (2.0 * focal), xdir.y / (2.0 * focal));
        let v1 = GpVec2d::new(u * v2.x() + ydir.x, u * v2.y() + ydir.y);
        let p = GpPnt2d::new(
            u * u / (4.0 * focal) * xdir.x + u * ydir.x + pos.point.x(),
            u * u / (4.0 * focal) * xdir.y + u * ydir.y + pos.point.y(),
        );
        (p, v1, v2)
    }
}

// ---------------------------------------------------------------------------
// D3
// ---------------------------------------------------------------------------

/// `ElCLib::CircleD3(U, gp_Ax22d, Radius)` (`cxx:809-839`): `V3 = -V1`.
pub fn circle_d3_ax22d(u: f64, pos: &GpAx22d, radius: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let (p, v1, v2) = circle_d2_ax22d(u, pos, radius);
    let v3 = v1.reversed();
    (p, v1, v2, v3)
}

/// `ElCLib::EllipseD3(U, gp_Ax22d, Major, Minor)` (`cxx:843-874`): `V3 = -V1`.
pub fn ellipse_d3_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let (p, v1, v2) = ellipse_d2_ax22d(u, pos, major_radius, minor_radius);
    let v3 = v1.reversed();
    (p, v1, v2, v3)
}

/// `ElCLib::HyperbolaD3(U, gp_Ax22d, Major, Minor)` (`cxx:878-907`): `V3 = V1`,
/// i.e. deliberately *not* reversed.
pub fn hyperbola_d3_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    let (p, v1, v2) = hyperbola_d2_ax22d(u, pos, major_radius, minor_radius);
    (p, v1, v2, v1)
}

// ---------------------------------------------------------------------------
// DN
// ---------------------------------------------------------------------------

/// `ElCLib::LineDN(U, gp_Ax2d, N)` (`cxx:1049-1055`).
pub fn line_dn_ax2d(_u: f64, pos: &GpAx2d, n: i32) -> GpVec2d {
    if n == 1 {
        GpVec2d::new(pos.vdir.x, pos.vdir.y)
    } else {
        GpVec2d::new(0.0, 0.0)
    }
}

/// `ElCLib::CircleDN(U, gp_Ax22d, Radius, N)` (`cxx:1060-1092`).
pub fn circle_dn_ax22d(u: f64, pos: &GpAx22d, radius: f64, n: i32) -> GpVec2d {
    let (mut xc, mut yc) = (0.0, 0.0);
    if n == 1 {
        xc = radius * -u.sin();
        yc = radius * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = radius * -u.cos();
        yc = radius * -u.sin();
    } else if (n + 1) % 4 == 0 {
        xc = radius * u.sin();
        yc = radius * -u.cos();
    } else if n % 4 == 0 {
        xc = radius * u.cos();
        yc = radius * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = radius * -u.sin();
        yc = radius * u.cos();
    }
    GpVec2d::new(
        xc * pos.vxdir.x + yc * pos.vydir.x,
        xc * pos.vxdir.y + yc * pos.vydir.y,
    )
}

/// `ElCLib::EllipseDN(U, gp_Ax22d, Major, Minor, N)` (`cxx:1095-1133`).
pub fn ellipse_dn_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
    n: i32,
) -> GpVec2d {
    let (mut xc, mut yc) = (0.0, 0.0);
    if n == 1 {
        xc = major_radius * -u.sin();
        yc = minor_radius * u.cos();
    } else if (n + 2) % 4 == 0 {
        xc = major_radius * -u.cos();
        yc = minor_radius * -u.sin();
    } else if (n + 1) % 4 == 0 {
        xc = major_radius * u.sin();
        yc = minor_radius * -u.cos();
    } else if n % 4 == 0 {
        xc = major_radius * u.cos();
        yc = minor_radius * u.sin();
    } else if (n - 1) % 4 == 0 {
        xc = major_radius * -u.sin();
        yc = minor_radius * u.cos();
    }
    GpVec2d::new(
        xc * pos.vxdir.x + yc * pos.vydir.x,
        xc * pos.vxdir.y + yc * pos.vydir.y,
    )
}

/// `ElCLib::HyperbolaDN(U, gp_Ax22d, Major, Minor, N)` (`cxx:1137-1158`).
pub fn hyperbola_dn_ax22d(
    u: f64,
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
    n: i32,
) -> GpVec2d {
    // `cxx:1145-1153`: `if (IsOdd(N)) ... else if (IsEven(N)) ...`. The two
    // predicates are exhaustive and mutually exclusive, so OCCT's
    // `double Xc = 0, Yc = 0;` initialisers are dead and are folded away here.
    let (xc, yc) = if is_odd(n) {
        (major_radius * u.sinh(), minor_radius * u.cosh())
    } else {
        (major_radius * u.cosh(), minor_radius * u.sinh())
    };
    GpVec2d::new(
        xc * pos.vxdir.x + yc * pos.vydir.x,
        xc * pos.vxdir.y + yc * pos.vydir.y,
    )
}

/// `ElCLib::ParabolaDN(U, gp_Ax22d, Focal, N)` (`cxx:1162-1188`).
pub fn parabola_dn_ax22d(u: f64, pos: &GpAx22d, focal: f64, n: i32) -> GpVec2d {
    if !(1..=2).contains(&n) {
        return GpVec2d::new(0.0, 0.0);
    }
    let xdir = GpXY::new(pos.vxdir.x, pos.vxdir.y);
    if n == 1 {
        if focal.abs() <= REAL_SMALL {
            return GpVec2d::from_xy(xdir);
        }
        let ydir = GpXY::new(pos.vydir.x, pos.vydir.y);
        return GpVec2d::new(u / (2.0 * focal) * xdir.x + ydir.x, u / (2.0 * focal) * xdir.y + ydir.y);
    }
    if focal.abs() <= REAL_SMALL {
        return GpVec2d::new(0.0, 0.0);
    }
    GpVec2d::new(xdir.x / (2.0 * focal), xdir.y / (2.0 * focal))
}

// ---------------------------------------------------------------------------
// Parameter
// ---------------------------------------------------------------------------

/// `ElCLib::LineParameter(gp_Ax2d, gp_Pnt2d)` (`cxx:1276-1281`).
pub fn line_parameter_ax2d(pos: &GpAx2d, p: &GpPnt2d) -> f64 {
    let mut coord = *p.xy();
    coord.subtract(pos.loc.xy());
    coord.dot(&GpXY::new(pos.vdir.x, pos.vdir.y))
}

/// `ElCLib::CircleParameter(gp_Ax22d, gp_Pnt2d)` (`cxx:1285-1291`).
///
/// The OCCT body is `Pos.XDirection().Angle(gp_Vec2d(Pos.Location(), P))`, which
/// relies on the implicit `gp_Dir2d(const gp_Vec2d&)` conversion
/// (`gp_Dir2d.hxx:67-69,290`) that normalizes the vector. A null vector makes
/// OCCT raise `Standard_ConstructionError`; the degenerate case returns `0.0`
/// here instead, matching the 3D `circle_parameter` in [`super::clib`].
pub fn circle_parameter_ax22d(pos: &GpAx22d, p: &GpPnt2d) -> f64 {
    let Ok(vdir) = GpDir2d::new(p.x() - pos.point.x(), p.y() - pos.point.y()) else {
        return 0.0;
    };
    let mut teta = dir2d_angle(&pos.vxdir, &vdir);
    if pos.vxdir.crossed(&pos.vydir) < 0.0 {
        teta = -teta;
    }
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::EllipseParameter(gp_Ax22d, Major, Minor, gp_Pnt2d)`
/// (`cxx:1294-1311`).
pub fn ellipse_parameter_ax22d(
    pos: &GpAx22d,
    major_radius: f64,
    minor_radius: f64,
    p: &GpPnt2d,
) -> f64 {
    let op = GpXY::new(p.x() - pos.point.x(), p.y() - pos.point.y());
    let xaxis = GpXY::new(pos.vxdir.x, pos.vxdir.y);
    let mut yaxis = GpXY::new(pos.vydir.x, pos.vydir.y);
    let mut om = xaxis.multiplied(op.dot(&xaxis));
    yaxis.multiply_scalar(op.dot(&yaxis) * (major_radius / minor_radius));
    om.add(&yaxis);
    let mut teta = dir2d_angle_vec(&xaxis, &om);
    if pos.vxdir.crossed(&pos.vydir) < 0.0 {
        teta = -teta;
    }
    normalize_angle(&mut teta);
    teta
}

/// `ElCLib::HyperbolaParameter(gp_Ax22d, Major /*ignored*/, Minor, gp_Pnt2d)`
/// (`cxx:1314-1327`). OCCT keeps the unused `MajorRadius` slot so the four
/// signatures stay uniform; it is unnamed there and unused here.
pub fn hyperbola_parameter_ax22d(
    pos: &GpAx22d,
    _major_radius: f64,
    minor_radius: f64,
    p: &GpPnt2d,
) -> f64 {
    let v = GpXY::new(pos.vydir.x, pos.vydir.y);
    let sht = GpXY::new(p.x() - pos.point.x(), p.y() - pos.point.y()).dot(&v) / minor_radius;
    sht.asinh()
}

/// `ElCLib::ParabolaParameter(gp_Ax22d, gp_Pnt2d)` (`cxx:1331-1335`).
pub fn parabola_parameter_ax22d(pos: &GpAx22d, p: &GpPnt2d) -> f64 {
    let directrix = GpXY::new(pos.vydir.x, pos.vydir.y);
    GpXY::new(p.x() - pos.point.x(), p.y() - pos.point.y()).dot(&directrix)
}

// ---------------------------------------------------------------------------
// Parameter of a `gp` primitive (ElCLib.lxx:353-381)
// ---------------------------------------------------------------------------

/// `ElCLib::Parameter(gp_Lin2d, gp_Pnt2d)` (`ElCLib.lxx:353-357`).
pub fn parameter_lin2d(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    line_parameter_ax2d(&l.pos, p)
}

/// `ElCLib::Parameter(gp_Circ2d, gp_Pnt2d)` (`ElCLib.lxx:359-363`).
pub fn parameter_circ2d(c: &GpCirc2d, p: &GpPnt2d) -> f64 {
    circle_parameter_ax22d(&c.pos, p)
}

/// `ElCLib::Parameter(gp_Elips2d, gp_Pnt2d)` (`ElCLib.lxx:365-369`).
pub fn parameter_elips2d(e: &GpElips2d, p: &GpPnt2d) -> f64 {
    ellipse_parameter_ax22d(&e.pos, e.major_radius, e.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Hypr2d, gp_Pnt2d)` (`ElCLib.lxx:371-375`).
pub fn parameter_hypr2d(h: &GpHypr2d, p: &GpPnt2d) -> f64 {
    hyperbola_parameter_ax22d(&h.pos, h.major_radius, h.minor_radius, p)
}

/// `ElCLib::Parameter(gp_Parab2d, gp_Pnt2d)` (`ElCLib.lxx:377-381`).
pub fn parameter_parab2d(prb: &GpParab2d, p: &GpPnt2d) -> f64 {
    parabola_parameter_ax22d(&prb.pos, p)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// `gp_Vec2d(xaxis).Angle(gp_Vec2d(om))` (`ElCLib.cxx:1307`): both operands go
/// through the implicit `gp_Dir2d(const gp_Vec2d&)` normalization
/// (`gp_Dir2d.hxx:67-69,290`). A null vector makes OCCT raise
/// `Standard_ConstructionError`; the degenerate case yields `0.0` here, which
/// the caller then normalizes in exactly the same way.
fn dir2d_angle_vec(a: &GpXY, b: &GpXY) -> f64 {
    let (Ok(a_dir), Ok(b_dir)) = (GpDir2d::new(a.x, a.y), GpDir2d::new(b.x, b.y)) else {
        return 0.0;
    };
    dir2d_angle(&a_dir, &b_dir)
}
