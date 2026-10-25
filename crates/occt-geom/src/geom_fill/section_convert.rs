//! The static services of `GeomFill.hxx` that do not need any state: the
//! parameterisation bookkeeping (`GetShape`, `Knots`, `Mults`), the circle
//! section builders (`GetCircle`, three overloads) and the two tolerance
//! helpers (`GetMinimalWeights`, `GetTolerance`).
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill.cxx:186-824` (OCCT 8.0.0).
//!
//! `GeomFill::Surface` (`GeomFill.cxx:42-182`) lives in
//! `geom_fill/utilities/mod.rs` (`GeomFill::surface`), next to the other
//! static services.
//!
//! Translation notes:
//! - `Convert_ParameterisationType` is
//!   [`occt_core::convert::ParameterisationType`].
//! - `NCollection_Array1` out-parameters (`Poles`, `Weights`, `DPoles`, …) are
//!   slices; `Poles.Length()` is `poles.len()` and `Lower()` is `0`, so OCCT's
//!   `jj = low + 2` loops start at index `2`.
//! - `GeomConvert::CurveToBSplineCurve(Sect, TConv)` (`cxx:248`, `:253`,
//!   `:337`) is reached through [`circle_arc_bspline`], which feeds the
//!   already-ported `GeomConvert` trimmed-circle arm
//!   (`convert_bspl/curve_convert.rs:157-206`, `GeomConvert.cxx:171-244`)
//!   through a trimmed-circle adapter, so the weights and poles are the same
//!   values OCCT produces.
//! - `GeomFill::GetCircle` returns `void` / `bool` in OCCT and can only fail by
//!   throwing from `GeomConvert`; the port returns `Result` there, but the
//!   `GetCircle` overloads themselves never fail (all three converter branches
//!   are total), so they stay non-`Result`.

use std::sync::Arc;

use occt_core::convert::{ConvertError, ParameterisationType};
use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpPnt, GpVec};
use occt_core::precision::PCONFUSION;

use crate::bspline_curve::GeomBSplineCurve;
use crate::circle::GeomCircle;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::geom_fill::polynomial_convertor::PolynomialConvertor;
use crate::geom_fill::quasi_angular_convertor::QuasiAngularConvertor;

/// A `Geom_TrimmedCurve` whose parameter is the *basis* parameter, which is
/// what `GeomConvert::CurveToBSplineCurve` expects (`GeomConvert.cxx:171-175`).
/// The port's own `GeomTrimmedCurve` normalises its parameter to `[0, 1]`, so it
/// cannot be used here.
struct CircleArc {
    circ: GpCirc,
    first: f64,
    last: f64,
}

impl Curve for CircleArc {
    fn d0(&self, u: f64) -> GpPnt {
        clib::circle_value(&self.circ, u)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        clib::circle_d1(&self.circ, u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        clib::circle_d2(&self.circ, u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        false
    }
    fn is_geom_trimmed(&self) -> bool {
        true
    }
    fn untrimmed_basis(&self) -> Option<(Arc<dyn Curve>, f64, f64)> {
        Some((Arc::new(GeomCircle::new(self.circ.clone())), self.first, self.last))
    }
    fn continuity(&self) -> u8 {
        6
    }
    fn transform(&mut self, t: &occt_core::gp::GpTrsf) {
        self.circ.transform(t);
    }
    fn reverse(&mut self) {
        let z = self.circ.pos.direction().reversed();
        self.circ.pos.set_direction(z);
        let two_pi = 2.0 * std::f64::consts::PI;
        let first = self.first;
        self.first = two_pi - self.last;
        self.last = two_pi - first;
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(Self { circ: self.circ.clone(), first: self.first, last: self.last })
    }
}

/// `GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(new
/// Geom_Circle(C), U1, U2), TConv)` for the circle `C`.
pub(crate) fn circle_arc_bspline_of(
    circ: &GpCirc,
    u1: f64,
    u2: f64,
    parameterisation: ParameterisationType,
) -> Result<GeomBSplineCurve, ConvertError> {
    let arc = CircleArc { circ: circ.clone(), first: u1, last: u2 };
    curve_to_bspline_curve(&arc, parameterisation)
}

/// The standard-position circle of radius `radius`, see
/// [`circle_arc_bspline_of`].
fn circle_arc_bspline(
    radius: f64,
    u1: f64,
    u2: f64,
    parameterisation: ParameterisationType,
) -> Result<GeomBSplineCurve, ConvertError> {
    circle_arc_bspline_of(&GpCirc::new(GpAx2::standard(), radius), u1, u2, parameterisation)
}

/// `GeomFill::GetShape(MaxAng, NbPoles, NbKnots, Degree, TConv)`
/// (`GeomFill.cxx:186-226`). `TConv` is an in-out parameter: the `default`
/// arm rewrites it to `Convert_TgtThetaOver2_1`..`_3` for one to three spans.
pub fn get_shape(
    max_ang: f64,
    nb_poles: &mut i32,
    nb_knots: &mut i32,
    degree: &mut i32,
    tconv: &mut ParameterisationType,
) {
    match *tconv {
        ParameterisationType::QuasiAngular => {
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
        }
        ParameterisationType::Polynomial => {
            *nb_poles = 8;
            *nb_knots = 2;
            *degree = 7;
        }
        _ => {
            let nb_span = (3.0 * max_ang.abs() / 2.0 / std::f64::consts::PI).ceil() as i32;
            *nb_poles = 2 * nb_span + 1;
            *nb_knots = nb_span + 1;
            *degree = 2;
            if nb_span == 1 {
                *tconv = ParameterisationType::TgtThetaOver2_1;
            } else if nb_span == 2 {
                *tconv = ParameterisationType::TgtThetaOver2_2;
            } else if nb_span == 3 {
                *tconv = ParameterisationType::TgtThetaOver2_3;
            }
        }
    }
}

/// `GeomFill::GetMinimalWeights(TConv, MinAng, MaxAng, Weights)`
/// (`GeomFill.cxx:236-263`). OCCT throws out of `GeomConvert` when the angles
/// are out of range; the port reports that through [`ConvertError`].
pub fn get_minimal_weights(
    tconv: ParameterisationType,
    min_ang: f64,
    max_ang: f64,
    weights: &mut [f64],
) -> Result<(), ConvertError> {
    if tconv == ParameterisationType::Polynomial {
        weights.fill(1.0);
        return Ok(());
    }
    let radius = 1.0;
    let sect1 = circle_arc_bspline(radius, 0.0, max_ang, tconv)?;
    let first = sect1.weights.clone().unwrap_or_else(|| vec![1.0; sect1.poles.len()]);
    let count = weights.len().min(first.len());
    weights[..count].copy_from_slice(&first[..count]);

    let angle_min = PCONFUSION.max(min_ang);
    let sect2 = circle_arc_bspline(radius, 0.0, angle_min, tconv)?;
    let poids = sect2.weights.clone().unwrap_or_else(|| vec![1.0; sect2.poles.len()]);
    let count = weights.len().min(poids.len());
    for ii in 0..count {
        if poids[ii] < weights[ii] {
            weights[ii] = poids[ii];
        }
    }
    Ok(())
}

/// `GeomFill::Knots(TConv, TKnots)` (`GeomFill.cxx:268-284`).
pub fn knots(tconv: ParameterisationType, t_knots: &mut [f64]) {
    if tconv != ParameterisationType::QuasiAngular && tconv != ParameterisationType::Polynomial {
        let mut val = 0.0;
        for knot in t_knots.iter_mut() {
            *knot = val;
            val += 1.0;
        }
    } else if t_knots.len() >= 2 {
        t_knots[0] = 0.0;
        t_knots[1] = 1.0;
    }
}

/// `GeomFill::Mults(TConv, TMults)` (`GeomFill.cxx:289-315`).
pub fn mults(tconv: ParameterisationType, t_mults: &mut [i32]) {
    match tconv {
        ParameterisationType::QuasiAngular => {
            if t_mults.len() >= 2 {
                t_mults[0] = 7;
                t_mults[1] = 7;
            }
        }
        ParameterisationType::Polynomial => {
            if t_mults.len() >= 2 {
                t_mults[0] = 8;
                t_mults[1] = 8;
            }
        }
        _ => {
            let n = t_mults.len();
            if n == 0 {
                return;
            }
            t_mults[0] = 3;
            for m in t_mults.iter_mut().take(n - 1).skip(1) {
                *m = 2;
            }
            t_mults[n - 1] = 3;
        }
    }
}

/// `GeomFill::GetTolerance(TConv, AngleMin, Radius, AngularTol, SpatialTol)`
/// (`GeomFill.cxx:323-340`): the 3D tolerance that respects the G1 continuity
/// tolerance. `0.02` is OCCT's own floor on the section angle.
pub fn get_tolerance(
    tconv: ParameterisationType,
    angle_min: f64,
    radius: f64,
    angular_tol: f64,
    spatial_tol: f64,
) -> Result<f64, ConvertError> {
    let sect = circle_arc_bspline(radius, 0.0, angle_min.max(0.02), tconv)?;
    let dist = if sect.poles.len() >= 2 {
        sect.poles[0].distance(&sect.poles[1]) + spatial_tol
    } else {
        spatial_tol
    };
    Ok(dist * angular_tol / 2.0)
}

/// The `Cosa` / `Sina` / `Angle` preamble shared by the three `GetCircle`
/// overloads (`GeomFill.cxx:369-395`, `:473-499`, `:644-668`), including the
/// recadrage of the angle into `]-pi/2, 3pi/2]`.
fn section_angle(ns1: &GpVec, ns2: &GpVec, nplan: &GpVec) -> (f64, f64, f64) {
    let mut cosa = ns1.dot(ns2);
    let mut sina = nplan.dot(&ns1.crossed(ns2));
    if cosa < -1.0 {
        cosa = -1.0;
        sina = 0.0;
    }
    if cosa > 1.0 {
        cosa = 1.0;
        sina = 0.0;
    }
    let mut angle = cosa.acos();
    if sina < 0.0 {
        if cosa > 0.0 {
            angle = -angle;
        } else {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
    }
    (cosa, sina, angle)
}

/// `GeomFill::GetCircle(TConv, ns1, ns2, nplan, pts1, pts2, Rayon, Center,
/// Poles, Weights)` (`GeomFill.cxx:350-447`).
#[allow(clippy::too_many_arguments)]
pub fn get_circle(
    tconv: ParameterisationType,
    ns1: &GpVec,
    ns2: &GpVec,
    nplan: &GpVec,
    pts1: &GpPnt,
    pts2: &GpPnt,
    rayon: f64,
    center: &GpPnt,
    poles: &mut [GpPnt],
    weights: &mut [f64],
) {
    let (_, _, angle) = section_angle(ns1, ns2, nplan);
    match tconv {
        ParameterisationType::QuasiAngular => {
            let mut convertor = QuasiAngularConvertor::new();
            convertor.init();
            let (p, w) = convertor.section(pts1, center, nplan, angle);
            copy_poles(poles, &p);
            copy_values(weights, &w);
        }
        ParameterisationType::Polynomial => {
            let mut convertor = PolynomialConvertor::new();
            convertor.init();
            let p = convertor.section(pts1, center, nplan, angle);
            copy_poles(poles, &p);
            weights.fill(1.0);
        }
        _ => {
            // Cas rational: expression directe, plus performante que GeomConvert.
            if poles.len() < 2 || weights.len() < 2 {
                return;
            }
            let low = 0usize;
            let upp = poles.len() - 1;
            let nb_span = (poles.len() - 1) / 2;
            poles[low] = pts1.clone();
            poles[upp] = pts2.clone();
            weights[low] = 1.0;
            weights[upp] = 1.0;

            let np2 = nplan.crossed(ns1);
            let alpha = angle / nb_span as f64;
            let cosas2 = (alpha / 2.0).cos();

            let mut i = 1usize;
            let mut jj = low + 2;
            while i <= nb_span.saturating_sub(1) {
                let lambda = i as f64 * alpha;
                let cosa = lambda.cos();
                let sina = lambda.sin();
                let temp = GpVec::zero().set_linear_form_2(cosa - 1.0, ns1, sina, &np2);
                if jj < poles.len() {
                    poles[jj] = GpPnt::from_xyz(&pts1.coord.added(&temp.coord.multiply(rayon)));
                    weights[jj] = 1.0;
                }
                i += 1;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            let mut i = 1usize;
            let mut jj = low + 1;
            while i <= nb_span {
                if jj >= 1 && jj + 1 < poles.len() {
                    let temp = poles[jj - 1]
                        .coord
                        .added(&poles[jj + 1].coord)
                        .subtracted(&center.coord.multiply(2.0));
                    poles[jj] = GpPnt::from_xyz(&center.coord.added(&temp.multiply(lambda)));
                    weights[jj] = cosas2;
                }
                i += 1;
                jj += 2;
            }
        }
    }
}

/// `GeomFill::GetCircle(TConv, ns1, ns2, dn1w, dn2w, nplan, dnplan, pts1,
/// pts2, tang1, tang2, Rayon, DRayon, Center, DCenter, Poles, DPoles, Weights,
/// DWeights)` (`GeomFill.cxx:445-602`).
#[allow(clippy::too_many_arguments)]
pub fn get_circle_d1(
    tconv: ParameterisationType,
    ns1: &GpVec,
    ns2: &GpVec,
    dn1w: &GpVec,
    dn2w: &GpVec,
    nplan: &GpVec,
    dnplan: &GpVec,
    pts1: &GpPnt,
    pts2: &GpPnt,
    tang1: &GpVec,
    tang2: &GpVec,
    rayon: f64,
    d_rayon: f64,
    center: &GpPnt,
    d_center: &GpVec,
    poles: &mut [GpPnt],
    d_poles: &mut [GpVec],
    weights: &mut [f64],
    d_weights: &mut [f64],
) {
    let (cosa, sina, angle) = section_angle(ns1, ns2, nplan);
    let d_angle = if sina.abs() > cosa.abs() {
        -(dn1w.dot(ns2) + ns1.dot(dn2w)) / sina
    } else {
        (dnplan.dot(&ns1.crossed(ns2)) + nplan.dot(&dn1w.crossed(ns2).added(&ns1.crossed(dn2w))))
            / cosa
    };

    if poles.len() < 2 || weights.len() < 2 || d_poles.len() < 2 || d_weights.len() < 2 {
        return;
    }
    // Aux extremites.
    let low = 0usize;
    let upp = poles.len() - 1;
    poles[low] = pts1.clone();
    poles[upp] = pts2.clone();
    weights[low] = 1.0;
    weights[upp] = 1.0;
    d_poles[low] = tang1.clone();
    d_poles[upp] = tang2.clone();
    d_weights[low] = 0.0;
    d_weights[upp] = 0.0;

    match tconv {
        ParameterisationType::QuasiAngular => {
            let mut convertor = QuasiAngularConvertor::new();
            convertor.init();
            let (p, dp, w, dw) = convertor.section_d1(
                pts1, tang1, center, d_center, nplan, dnplan, angle, d_angle,
            );
            copy_poles(poles, &p);
            copy_vectors(d_poles, &dp);
            copy_values(weights, &w);
            copy_values(d_weights, &dw);
        }
        ParameterisationType::Polynomial => {
            let mut convertor = PolynomialConvertor::new();
            convertor.init();
            let (p, dp) = convertor.section_d1(
                pts1, tang1, center, d_center, nplan, dnplan, angle, d_angle,
            );
            copy_poles(poles, &p);
            copy_vectors(d_poles, &dp);
            weights.fill(1.0);
            d_weights.fill(0.0);
        }
        _ => {
            let nb_span = (poles.len() - 1) / 2;
            let np2 = nplan.crossed(ns1);
            let dnp2 = dnplan.crossed(ns1).added(&nplan.crossed(dn1w));
            let alpha = angle / nb_span as f64;
            let cosas2 = (alpha / 2.0).cos();
            let sinas2 = (alpha / 2.0).sin();

            let mut i = 1usize;
            let mut jj = low + 2;
            while i <= nb_span.saturating_sub(1) {
                let lambda = i as f64 * alpha;
                let cosa_i = lambda.cos();
                let sina_i = lambda.sin();
                let mut temp = GpVec::zero().set_linear_form_2(cosa_i - 1.0, ns1, sina_i, &np2);
                if jj < poles.len() {
                    poles[jj] = GpPnt::from_xyz(&pts1.coord.added(&temp.coord.multiply(rayon)));

                    let mut dp = GpVec::zero().set_linear_form_2(d_rayon, &temp, 1.0, tang1);
                    temp = GpVec::zero().set_linear_form_2(-sina_i, ns1, cosa_i, &np2);
                    temp = temp.multiply_scalar(i as f64 / nb_span as f64 * d_angle);
                    temp = temp.add(
                        &GpVec::zero().set_linear_form_2(cosa_i - 1.0, dn1w, sina_i, &dnp2),
                    );
                    dp = dp.added(&temp.multiplied_scalar(rayon));
                    if jj < d_poles.len() {
                        d_poles[jj] = dp;
                    }
                }
                i += 1;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            let dlambda = (lambda * sinas2 * d_angle) / (cosas2 * nb_span as f64);
            let mut i = 1usize;
            let mut jj = low;
            while i <= nb_span {
                if jj + 2 < poles.len() {
                    let temp = poles[jj]
                        .coord
                        .added(&poles[jj + 2].coord)
                        .subtracted(&center.coord.multiply(2.0));
                    poles[jj + 1] = GpPnt::from_xyz(&center.coord.added(&temp.multiply(lambda)));
                    if jj + 1 < d_poles.len() && jj + 2 < d_poles.len() {
                        d_poles[jj + 1] = GpVec::zero().set_linear_form_3(
                            dlambda,
                            &GpVec::from_xyz(&temp),
                            1.0 - 2.0 * lambda,
                            d_center,
                            lambda,
                            &d_poles[jj].added(&d_poles[jj + 2]),
                        );
                    }
                }
                i += 1;
                jj += 2;
            }

            // Les poids.
            let dlambda = -sinas2 * d_angle / (2.0 * nb_span as f64);
            let mut i = low;
            while i < upp {
                weights[i] = 1.0;
                weights[i + 1] = cosas2;
                d_weights[i] = 0.0;
                d_weights[i + 1] = dlambda;
                i += 2;
            }
        }
    }
}

/// `GeomFill::GetCircle(TConv, ns1, ns2, dn1w, dn2w, d2n1w, d2n2w, nplan,
/// dnplan, d2nplan, pts1, pts2, tang1, tang2, Dtang1, Dtang2, Rayon, DRayon,
/// D2Rayon, Center, DCenter, D2Center, Poles, DPoles, D2Poles, Weights,
/// DWeights, D2Weights)` (`GeomFill.cxx:604-824`).
#[allow(clippy::too_many_arguments)]
pub fn get_circle_d2(
    tconv: ParameterisationType,
    ns1: &GpVec,
    ns2: &GpVec,
    dn1w: &GpVec,
    dn2w: &GpVec,
    d2n1w: &GpVec,
    d2n2w: &GpVec,
    nplan: &GpVec,
    dnplan: &GpVec,
    d2nplan: &GpVec,
    pts1: &GpPnt,
    pts2: &GpPnt,
    tang1: &GpVec,
    tang2: &GpVec,
    d_tang1: &GpVec,
    d_tang2: &GpVec,
    rayon: f64,
    d_rayon: f64,
    d2_rayon: f64,
    center: &GpPnt,
    d_center: &GpVec,
    d2_center: &GpVec,
    poles: &mut [GpPnt],
    d_poles: &mut [GpVec],
    d2_poles: &mut [GpVec],
    weights: &mut [f64],
    d_weights: &mut [f64],
    d2_weights: &mut [f64],
) {
    let (cosa, sina, angle) = section_angle(ns1, ns2, nplan);
    let (d_angle, d2_angle);
    if sina.abs() > cosa.abs() {
        let aux = dn1w.dot(ns2) + ns1.dot(dn2w);
        d_angle = -aux / sina;
        d2_angle = -(d2n1w.dot(ns2) + 2.0 * dn1w.dot(dn2w) + ns1.dot(d2n2w)) / sina
            + aux * (dnplan.dot(&ns1.crossed(ns2))
                + nplan.dot(&dn1w.crossed(ns2).added(&ns1.crossed(dn2w))))
                / (sina * sina);
    } else {
        let temp = dn1w.crossed(ns2).added(&ns1.crossed(dn2w));
        d_angle = (dnplan.dot(&ns1.crossed(ns2)) + nplan.dot(&temp)) / cosa;
        d2_angle = (d2nplan.dot(&ns1.crossed(ns2)) + 2.0 * dnplan.dot(&temp)
            + nplan.dot(
                &d2n1w
                    .crossed(ns2)
                    .added(&dn1w.crossed(dn2w).multiplied_scalar(2.0))
                    .added(&ns1.crossed(d2n2w)),
            ))
            / cosa
            - (dn1w.dot(ns2) + ns1.dot(dn2w))
                * (dnplan.dot(&ns1.crossed(ns2)) + nplan.dot(&temp))
                / (cosa * cosa);
    }

    if poles.len() < 2
        || weights.len() < 2
        || d_poles.len() < 2
        || d2_poles.len() < 2
        || d_weights.len() < 2
        || d2_weights.len() < 2
    {
        return;
    }
    // Aux extremites.
    let low = 0usize;
    let upp = poles.len() - 1;
    poles[low] = pts1.clone();
    poles[upp] = pts2.clone();
    weights[low] = 1.0;
    weights[upp] = 1.0;
    d_poles[low] = tang1.clone();
    d_poles[upp] = tang2.clone();
    d_weights[low] = 0.0;
    d_weights[upp] = 0.0;
    d2_poles[low] = d_tang1.clone();
    d2_poles[upp] = d_tang2.clone();
    d2_weights[low] = 0.0;
    d2_weights[upp] = 0.0;

    match tconv {
        ParameterisationType::QuasiAngular => {
            let mut convertor = QuasiAngularConvertor::new();
            convertor.init();
            let (p, dp, d2p, w, dw, d2w) = convertor.section_d2(
                pts1, tang1, d_tang1, center, d_center, d2_center, nplan, dnplan, d2nplan, angle,
                d_angle, d2_angle,
            );
            copy_poles(poles, &p);
            copy_vectors(d_poles, &dp);
            copy_vectors(d2_poles, &d2p);
            copy_values(weights, &w);
            copy_values(d_weights, &dw);
            copy_values(d2_weights, &d2w);
        }
        ParameterisationType::Polynomial => {
            let mut convertor = PolynomialConvertor::new();
            convertor.init();
            let (p, dp, d2p) = convertor.section_d2(
                pts1, tang1, d_tang1, center, d_center, d2_center, nplan, dnplan, d2nplan, angle,
                d_angle, d2_angle,
            );
            copy_poles(poles, &p);
            copy_vectors(d_poles, &dp);
            copy_vectors(d2_poles, &d2p);
            weights.fill(1.0);
            d_weights.fill(0.0);
            d2_weights.fill(0.0);
        }
        _ => {
            let nb_span = (poles.len() - 1) / 2;
            let np2 = nplan.crossed(ns1);
            let dnp2 = dnplan.crossed(ns1).added(&nplan.crossed(dn1w));
            let d2np2 = d2nplan
                .crossed(ns1)
                .added(&nplan.crossed(dn2w))
                .added(&dnplan.crossed(dn1w).multiplied_scalar(2.0));
            let alpha = angle / nb_span as f64;
            let cosas2 = (alpha / 2.0).cos();
            let sinas2 = (alpha / 2.0).sin();

            let mut i = 1usize;
            let mut jj = low + 2;
            while i <= nb_span.saturating_sub(1) {
                let lambda = i as f64 * alpha;
                let cosa_i = lambda.cos();
                let sina_i = lambda.sin();
                let temp = GpVec::zero().set_linear_form_2(cosa_i - 1.0, ns1, sina_i, &np2);
                if jj < poles.len() {
                    poles[jj] = GpPnt::from_xyz(&pts1.coord.added(&temp.coord.multiply(rayon)));

                    let mut dp = GpVec::zero().set_linear_form_2(d_rayon, &temp, 1.0, tang1);
                    let mut dtemp = GpVec::zero().set_linear_form_2(-sina_i, ns1, cosa_i, &np2);
                    let aux = i as f64 / nb_span as f64;
                    dtemp = dtemp.multiply_scalar(aux * d_angle);
                    dtemp = dtemp.add(
                        &GpVec::zero().set_linear_form_2(cosa_i - 1.0, dn1w, sina_i, &dnp2),
                    );
                    dp = dp.added(&dtemp.multiplied_scalar(rayon));
                    if jj < d_poles.len() {
                        d_poles[jj] = dp.clone();
                    }

                    if jj < d2_poles.len() {
                        let mut d2p = GpVec::zero().set_linear_form_3(
                            d2_rayon,
                            &temp,
                            2.0 * d_rayon,
                            &dtemp,
                            1.0,
                            &d_tang1,
                        );
                        let mut t2 =
                            GpVec::zero().set_linear_form_2(cosa_i - 1.0, dn2w, sina_i, &d2np2);
                        let mut dt2 = GpVec::zero().set_linear_form_2(-sina_i, ns1, cosa_i, &np2);
                        t2 = t2.added(&dt2.multiplied_scalar(aux * aux * d2_angle));
                        dt2 = GpVec::zero().set_linear_form_3(
                            -sina_i,
                            &dn1w.added(&np2),
                            cosa_i,
                            &dnp2,
                            -cosa_i,
                            ns1,
                        );
                        t2 = t2.added(&dt2.multiplied_scalar(aux * d_angle));
                        d2p = d2p.added(&t2.multiplied_scalar(rayon));
                        d2_poles[jj] = d2p;
                    }
                }
                i += 1;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            let dlambda = (lambda * sinas2 * d_angle) / (cosas2 * nb_span as f64);
            let aux = sinas2 / cosas2;
            let d2lambda = (dlambda * aux * d_angle
                + d2_angle * aux * lambda
                + (1.0 + aux * aux) * (d_angle / (2.0 * nb_span as f64)) * d_angle * lambda)
                / nb_span as f64;
            let mut i = 1usize;
            let mut jj = low;
            while i <= nb_span {
                if jj + 2 < poles.len() {
                    let temp = poles[jj]
                        .coord
                        .added(&poles[jj + 2].coord)
                        .subtracted(&center.coord.multiply(2.0));
                    poles[jj + 1] = GpPnt::from_xyz(&center.coord.added(&temp.multiply(lambda)));
                    if jj + 1 < d_poles.len() && jj + 2 < d_poles.len() {
                        let dtemp = d_poles[jj]
                            .coord
                            .added(&d_poles[jj + 2].coord)
                            .subtracted(&d_center.coord.multiply(2.0));
                        d_poles[jj + 1] = GpVec::zero().set_linear_form_3(
                            dlambda,
                            &GpVec::from_xyz(&temp),
                            lambda,
                            &GpVec::from_xyz(&dtemp),
                            1.0,
                            d_center,
                        );
                        if jj + 1 < d2_poles.len() && jj + 2 < d2_poles.len() {
                            let mut d2p = GpVec::zero().set_linear_form_3(
                                d2lambda,
                                &GpVec::from_xyz(&temp),
                                2.0 * dlambda,
                                &GpVec::from_xyz(&dtemp),
                                lambda,
                                &d2_poles[jj].added(&d2_poles[jj + 2]),
                            );
                            d2p = d2p.added(&d2_center.multiplied_scalar(1.0 - 2.0 * lambda));
                            d2_poles[jj + 1] = d2p;
                        }
                    }
                }
                i += 1;
                jj += 2;
            }

            // Les poids.
            let dlambda = -sinas2 * d_angle / (2.0 * nb_span as f64);
            let d2lambda = -sinas2 * d2_angle / (2.0 * nb_span as f64)
                - cosas2 * (d_angle / (2.0 * nb_span as f64)).powi(2);
            let mut i = low;
            while i < upp {
                weights[i] = 1.0;
                weights[i + 1] = cosas2;
                d_weights[i] = 0.0;
                d_weights[i + 1] = dlambda;
                d2_weights[i] = 0.0;
                d2_weights[i + 1] = d2lambda;
                i += 2;
            }
        }
    }
}

fn copy_poles(dst: &mut [GpPnt], src: &[GpPnt]) {
    let n = dst.len().min(src.len());
    dst[..n].clone_from_slice(&src[..n]);
}

fn copy_vectors(dst: &mut [GpVec], src: &[GpVec]) {
    let n = dst.len().min(src.len());
    dst[..n].clone_from_slice(&src[..n]);
}

fn copy_values(dst: &mut [f64], src: &[f64]) {
    let n = dst.len().min(src.len());
    dst[..n].copy_from_slice(&src[..n]);
}
