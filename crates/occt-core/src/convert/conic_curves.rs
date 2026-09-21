//! Ports of the conic → B-spline converters (OCCT `Convert/`, TKMath):
//! `Convert_CircleToBSplineCurve` (`Convert_CircleToBSplineCurve.cxx:47-172`),
//! `Convert_EllipseToBSplineCurve` (`Convert_EllipseToBSplineCurve.cxx:47-174`),
//! `Convert_HyperbolaToBSplineCurve` (`Convert_HyperbolaToBSplineCurve.cxx:32-79`)
//! and `Convert_ParabolaToBSplineCurve` (`Convert_ParabolaToBSplineCurve.cxx:32-69`).
//!
//! Each converter derives from `Convert_ConicToBSplineCurve` (ported in
//! `convert/conic_to_bspline.rs`), asks it for the `cos`/`sin` tables of the
//! requested parameterisation, scales them by the conic's radii (with the sign
//! of the conic frame's orientation), and finally moves the poles from the
//! conic's own frame into the world plane with
//! `gp_Trsf2d::SetTransformation(conic.XAxis(), gp::OX2d())`.
//!
//! The circle and the ellipse have a **periodic** constructor
//! (`Convert_CircleToBSplineCurve(C, Parameterisation)`, `cxx:47-110`) and a
//! **ranged** one (`(C, UFirst, ULast, Parameterisation)`, `cxx:117-172`);
//! only `Convert_TgtThetaOver2` and `Convert_RationalC1` keep the periodicity
//! (`cxx:59-84`). The hyperbola and the parabola only have the ranged form and
//! always produce a degree-2 rational (parabola: polynomial, all weights 1)
//! arc of 3 poles / 2 knots of multiplicity 3.
//!
//! Consumers: `GeomConvert::CurveToBSplineCurve`
//! (`GeomConvert_CurveToBSplineCurve` arms `GeomConvert.cxx:208-298`, `:363-408`)
//! and, through it, `GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)`
//! (`GeomToIGES_GeomCurve.cxx:620-645`, board card R2-20).

use crate::convert::conic_to_bspline::{
    build_cos_and_sin, build_cos_and_sin_periodic, ConicToBSplineCurve, ConvertError,
    ParameterisationType,
};
use crate::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpParab2d, GpPnt2d, GpTrsf2d};
use crate::precision;

/// `gp::OX2d()` (`gp.hxx`).
fn ox2d() -> GpAx2d {
    GpAx2d::new(GpPnt2d::zero(), GpDir2d::default())
}

/// The sign OCCT gives the second radius: `+r` when the conic frame is
/// right-handed (`Ox.X() * Oy.Y() - Ox.Y() * Oy.X() > 0`), `-r` otherwise
/// (`Convert_CircleToBSplineCurve.cxx:88-99`,
/// `Convert_EllipseToBSplineCurve.cxx:90-101`).
fn signed_radius(pos: &GpAx22d, r: f64) -> f64 {
    let ox = pos.x_direction();
    let oy = pos.y_direction();
    if ox.x * oy.y - ox.y * oy.x > 0.0 {
        r
    } else {
        -r
    }
}

/// `gp_Trsf2d::SetTransformation(conic.XAxis(), gp::OX2d())` +
/// `myPoles(ii).SetCoord(1, R * CosNumerator(ii))`,
/// `myPoles(ii).SetCoord(2, value * SinNumerator(ii))`,
/// `myPoles(ii).Transform(Trsf)` (`Convert_CircleToBSplineCurve.cxx:101-109`,
/// `Convert_EllipseToBSplineCurve.cxx:103-111`).
fn scaled_poles(
    cos_numerator: &[f64],
    sin_numerator: &[f64],
    rx: f64,
    ry: f64,
    x_axis: &GpAx2d,
) -> Vec<GpPnt2d> {
    let mut trsf = GpTrsf2d::identity();
    trsf.set_transformation(x_axis, &ox2d());
    cos_numerator
        .iter()
        .zip(sin_numerator.iter())
        .map(|(&c, &s)| {
            let mut p = GpPnt2d::new(rx * c, ry * s);
            p.transform(&trsf);
            p
        })
        .collect()
}

/// `Convert_CircleToBSplineCurve(C, Parameterisation)` (`cxx:47-110`): the
/// periodic circle (`myIsPeriodic = true` only for `Convert_TgtThetaOver2` and
/// `Convert_RationalC1`).
pub fn circle_to_bspline_curve(
    c: &GpCirc2d,
    parameterisation: ParameterisationType,
) -> Result<ConicToBSplineCurve, ConvertError> {
    let r = c.radius();
    let periodic = parameterisation == ParameterisationType::TgtThetaOver2
        || parameterisation == ParameterisationType::RationalC1;
    let cs = if periodic {
        // `cxx:74-84`.
        build_cos_and_sin_periodic(parameterisation)?
    } else {
        // `cxx:59-73`: trim on 0, 2*PI.
        build_cos_and_sin(parameterisation, 0.0, 2.0 * std::f64::consts::PI)?
    };
    let ry = signed_radius(c.position(), r);
    let poles = scaled_poles(
        &cs.cos_numerator,
        &cs.sin_numerator,
        r,
        ry,
        &c.x_axis(),
    );
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        cs.denominator,
        cs.knots,
        cs.mults,
        cs.degree,
        periodic,
    ))
}

/// `Convert_CircleToBSplineCurve(C, UFirst, ULast, Parameterisation)`
/// (`cxx:117-172`).
pub fn circle_to_bspline_curve_range(
    c: &GpCirc2d,
    u_first: f64,
    u_last: f64,
    parameterisation: ParameterisationType,
) -> Result<ConicToBSplineCurve, ConvertError> {
    let delta = u_last - u_first;
    // `Standard_DomainError_Raise_if((delta > (2 * M_PI + Eps)) || (delta <= 0.),
    //  "Convert_CircleToBSplineCurve")` (`cxx:127-130`).
    if delta > 2.0 * std::f64::consts::PI + precision::PCONFUSION || delta <= 0.0 {
        return Err(ConvertError::DomainError);
    }
    let r = c.radius();
    let cs = build_cos_and_sin(parameterisation, u_first, u_last)?;
    let ry = signed_radius(c.position(), r);
    let poles = scaled_poles(
        &cs.cos_numerator,
        &cs.sin_numerator,
        r,
        ry,
        &c.x_axis(),
    );
    // `myIsPeriodic = false` (`cxx:137`).
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        cs.denominator,
        cs.knots,
        cs.mults,
        cs.degree,
        false,
    ))
}

/// `Convert_EllipseToBSplineCurve(E, Parameterisation)` (`cxx:47-112`): the
/// periodic ellipse.
pub fn ellipse_to_bspline_curve(
    e: &GpElips2d,
    parameterisation: ParameterisationType,
) -> Result<ConicToBSplineCurve, ConvertError> {
    let major = e.major_radius;
    let minor = e.minor_radius;
    let periodic = parameterisation == ParameterisationType::TgtThetaOver2
        || parameterisation == ParameterisationType::RationalC1;
    let cs = if periodic {
        // `cxx:76-86`.
        build_cos_and_sin_periodic(parameterisation)?
    } else {
        // `cxx:61-75`: trim on 0, 2*PI.
        build_cos_and_sin(parameterisation, 0.0, 2.0 * std::f64::consts::PI)?
    };
    let ry = signed_radius(e.axis(), minor);
    let poles = scaled_poles(
        &cs.cos_numerator,
        &cs.sin_numerator,
        major,
        ry,
        &e.x_axis(),
    );
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        cs.denominator,
        cs.knots,
        cs.mults,
        cs.degree,
        periodic,
    ))
}

/// `Convert_EllipseToBSplineCurve(E, UFirst, ULast, Parameterisation)`
/// (`cxx:119-174`).
pub fn ellipse_to_bspline_curve_range(
    e: &GpElips2d,
    u_first: f64,
    u_last: f64,
    parameterisation: ParameterisationType,
) -> Result<ConicToBSplineCurve, ConvertError> {
    let delta = u_last - u_first;
    // `Standard_DomainError_Raise_if((delta > (2 * M_PI + Tol)) || (delta <= 0.),
    //  "Convert_EllipseToBSplineCurve")` (`cxx:127-131`).
    if delta > 2.0 * std::f64::consts::PI + precision::PCONFUSION || delta <= 0.0 {
        return Err(ConvertError::DomainError);
    }
    let major = e.major_radius;
    let minor = e.minor_radius;
    let cs = build_cos_and_sin(parameterisation, u_first, u_last)?;
    let ry = signed_radius(e.axis(), minor);
    let poles = scaled_poles(
        &cs.cos_numerator,
        &cs.sin_numerator,
        major,
        ry,
        &e.x_axis(),
    );
    // `myIsPeriodic = false` (`cxx:139`).
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        cs.denominator,
        cs.knots,
        cs.mults,
        cs.degree,
        false,
    ))
}

/// `Convert_HyperbolaToBSplineCurve(H, U1, U2)` (`cxx:32-79`): 3 poles, 2 knots
/// of multiplicity 3, degree 2; the middle weight is `cosh((UL - UF) / 2)`.
pub fn hyperbola_to_bspline_curve(
    h: &GpHypr2d,
    u1: f64,
    u2: f64,
) -> Result<ConicToBSplineCurve, ConvertError> {
    // `Standard_DomainError_Raise_if(Abs(U2 - U1) < Epsilon(0.),
    //  "Convert_HyperbolaToBSplineCurve")` (`cxx:38`). `Epsilon(0.)` is the
    // smallest positive double (`Standard_Real.hxx:242-246`).
    if (u2 - u1).abs() < precision::epsilon(0.0) {
        return Err(ConvertError::DomainError);
    }
    let uf = u1.min(u2);
    let ul = u1.max(u2);

    let major = h.major_radius;
    let minor = h.minor_radius;
    let pos = h.axis();
    let ox = pos.x_direction();
    let oy = pos.y_direction();
    let s = if ox.x * oy.y - ox.y * oy.x > 0.0 { 1.0 } else { -1.0 };

    // `cxx:62-71`.
    let weights = vec![1.0, ((ul - uf) / 2.0).cosh(), 1.0];
    let delta = (ul - uf).sinh();
    let x = major * (ul.sinh() - uf.sinh()) / delta;
    let y = s * minor * (ul.cosh() - uf.cosh()) / delta;
    let mut poles = vec![
        GpPnt2d::new(major * uf.cosh(), s * minor * uf.sinh()),
        GpPnt2d::new(x, y),
        GpPnt2d::new(major * ul.cosh(), s * minor * ul.sinh()),
    ];
    let mut trsf = GpTrsf2d::identity();
    trsf.set_transformation(&pos.x_axis(), &ox2d());
    for p in poles.iter_mut() {
        p.transform(&trsf);
    }
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        weights,
        vec![uf, ul],
        vec![3, 3],
        2,
        // `myIsPeriodic = false` (`cxx:43`).
        false,
    ))
}

/// `Convert_ParabolaToBSplineCurve(Prb, U1, U2)` (`cxx:32-69`): 3 poles, 2 knots
/// of multiplicity 3, degree 2, polynomial (all weights 1). The axis is `Oy`
/// and the parabola is `y^2 = 2*p*x` in its own frame, so `Prb.Parameter()`
/// (`gp_Parab2d::Parameter() = 2 * Focal()`) is the `2*p` divisor.
pub fn parabola_to_bspline_curve(
    prb: &GpParab2d,
    u1: f64,
    u2: f64,
) -> Result<ConicToBSplineCurve, ConvertError> {
    // `Standard_DomainError_Raise_if(Abs(U2 - U1) < Epsilon(0.),
    //  "Convert_ParabolaToBSplineCurve")` (`cxx:37`).
    if (u2 - u1).abs() < precision::epsilon(0.0) {
        return Err(ConvertError::DomainError);
    }
    let uf = u1.min(u2);
    let ul = u1.max(u2);

    let p = prb.parameter();
    let pos = prb.axis();
    let ox = pos.x_direction();
    let oy = pos.y_direction();
    let s = if ox.x * oy.y - ox.y * oy.x > 0.0 { 1.0 } else { -1.0 };

    // `cxx:50-61`.
    let mut poles = vec![
        GpPnt2d::new((uf * uf) / (2.0 * p), s * uf),
        GpPnt2d::new((uf * ul) / (2.0 * p), s * (uf + ul) / 2.0),
        GpPnt2d::new((ul * ul) / (2.0 * p), s * ul),
    ];
    let mut trsf = GpTrsf2d::identity();
    trsf.set_transformation(&pos.x_axis(), &ox2d());
    for pnt in poles.iter_mut() {
        pnt.transform(&trsf);
    }
    Ok(ConicToBSplineCurve::from_parts(
        poles,
        vec![1.0, 1.0, 1.0],
        vec![uf, ul],
        vec![3, 3],
        2,
        // `myIsPeriodic = false` (`cxx:44`).
        false,
    ))
}
