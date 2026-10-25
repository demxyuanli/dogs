//! `GeomFill::GetCircle` — the three overloads that build the arc pole /
//! weight array of a section, with D1 and D2.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill.cxx:350-827` (OCCT 8.0.0).

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpVec};

use crate::geom_fill::polynomial_convertor::PolynomialConvertor;
use crate::geom_fill::quasi_angular_convertor::QuasiAngularConvertor;
use crate::geom_fill::utilities::GeomFill;

impl GeomFill {
    /// `GeomFill::GetCircle(TConv, ns1, ns2, nplan, pts1, pts2, Rayon, Center,
    /// Poles, Weights)` (`cxx:350-446`).
    ///
    /// `ns1` / `ns2` are the inward normals at the two arc ends, `nplan` the
    /// plane normal; `Poles` and `Weights` are sized from `GetShape`.
    pub fn get_circle(
        t_conv: ParameterisationType,
        ns1: &GpVec,
        ns2: &GpVec,
        nplan: &GpVec,
        pts1: &GpPnt,
        pts2: &GpPnt,
        rayon: f64,
        center: &GpPnt,
        poles: &mut [GpPnt],
        weights: &mut [f64],
    ) -> Result<(), String> {
        let (angle, _sina, _cosa) = arc_angle(ns1, ns2, nplan);
        match t_conv {
            ParameterisationType::QuasiAngular => {
                let mut q = QuasiAngularConvertor::new();
                let (p, w) = q.section(pts1, center, nplan, angle);
                copy_poles_weights(poles, weights, &p, &w)
            }
            ParameterisationType::Polynomial => {
                let mut pc = PolynomialConvertor::new();
                let p = pc.section(pts1, center, nplan, angle);
                if poles.len() != p.len() {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                poles.copy_from_slice(&p);
                weights.fill(1.0);
                Ok(())
            }
            _ => {
                // `cxx:412-443`: the rational expression, far cheaper than
                // `GeomConvert` and consistent with `GetCircle` D1 / D2.
                let nb_span = (poles.len() - 1) / 2;
                if nb_span == 0 {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                let upp = poles.len() - 1;
                poles[0] = *pts1;
                poles[upp] = *pts2;
                weights[0] = 1.0;
                weights[upp] = 1.0;

                let np2 = nplan.crossed(ns1);
                let alpha = angle / nb_span as f64;
                let cosas2 = (alpha / 2.0).cos();

                for i in 1..=(nb_span - 1) {
                    let jj = 2 * i;
                    let lambda = i as f64 * alpha;
                    let cosa = lambda.cos();
                    let sina = lambda.sin();
                    let temp = ns1
                        .multiplied_scalar(cosa - 1.0)
                        .added(&np2.multiplied_scalar(sina));
                    poles[jj] = pts1.translated_vec(&temp.multiplied_scalar(rayon));
                    weights[jj] = 1.0;
                }

                let lambda = 1.0 / (2.0 * cosas2 * cosas2);
                for i in 1..=nb_span {
                    let jj = 2 * i - 1;
                    let temp = GpVec::from_pnts(
                        center,
                        &poles[jj - 1].translated_vec(&GpVec::from_pnts(center, &poles[jj + 1])),
                    );
                    poles[jj] = center.translated_vec(&temp.multiplied_scalar(lambda));
                    weights[jj] = cosas2;
                }
                Ok(())
            }
        }
    }

    /// `GeomFill::GetCircle(TConv, ns1, ns2, dn1w, dn2w, nplan, dnplan, pts1,
    /// pts2, tang1, tang2, Rayon, DRayon, Center, DCenter, Poles, DPoles,
    /// Weights, DWeights)` (`cxx:448-602`).
    #[allow(clippy::too_many_arguments)]
    pub fn get_circle_d1(
        t_conv: ParameterisationType,
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
    ) -> Result<(), String> {
        let (angle, sina, cosa) = arc_angle(ns1, ns2, nplan);
        let d_angle = if sina.abs() > cosa.abs() {
            -(dn1w.dot(ns2) + ns1.dot(dn2w)) / sina
        } else {
            (dnplan.dot(&ns1.crossed(ns2))
                + nplan.dot(&dn1w.crossed(ns2).added(&ns1.crossed(dn2w))))
                / cosa
        };

        let low = 0;
        let upp = poles.len() - 1;
        poles[low] = *pts1;
        poles[upp] = *pts2;
        weights[low] = 1.0;
        weights[upp] = 1.0;
        d_poles[low] = *tang1;
        d_poles[upp] = *tang2;
        d_weights[low] = 0.0;
        d_weights[upp] = 0.0;

        match t_conv {
            ParameterisationType::QuasiAngular => {
                let mut q = QuasiAngularConvertor::new();
                let (p, dp, w, dw) = q.section_d1(
                    pts1, tang1, center, d_center, nplan, dnplan, angle, d_angle,
                );
                copy_poles_weights(poles, weights, &p, &w)?;
                if d_poles.len() != dp.len() || d_weights.len() != dw.len() {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                d_poles.copy_from_slice(&dp);
                d_weights.copy_from_slice(&dw);
                Ok(())
            }
            ParameterisationType::Polynomial => {
                let mut pc = PolynomialConvertor::new();
                let (p, dp) = pc.section_d1(
                    pts1, tang1, center, d_center, nplan, dnplan, angle, d_angle,
                );
                if poles.len() != p.len() {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                poles.copy_from_slice(&p);
                d_poles.copy_from_slice(&dp);
                weights.fill(1.0);
                d_weights.fill(0.0);
                Ok(())
            }
            _ => {
                // `cxx:553-600`: the rational expression with derivatives.
                let nb_span = (poles.len() - 1) / 2;
                if nb_span == 0 {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                let np2 = nplan.crossed(ns1);
                let dnp2 = dnplan.crossed(ns1).added(&nplan.crossed(dn1w));

                let alpha = angle / nb_span as f64;
                let cosas2 = (alpha / 2.0).cos();
                let sinas2 = (alpha / 2.0).sin();

                for i in 1..=(nb_span - 1) {
                    let jj = 2 * i;
                    let lambda = i as f64 * alpha;
                    let cosa = lambda.cos();
                    let sina = lambda.sin();
                    let temp = ns1
                        .multiplied_scalar(cosa - 1.0)
                        .added(&np2.multiplied_scalar(sina));
                    poles[jj] = pts1.translated_vec(&temp.multiplied_scalar(rayon));

                    d_poles[jj] = temp.multiplied_scalar(d_rayon).added(tang1);
                    let mut t = ns1
                        .multiplied_scalar(-sina)
                        .added(&np2.multiplied_scalar(cosa));
                    t = t.multiplied_scalar(i as f64 / nb_span as f64 * d_angle);
                    t = t.added(
                        &dn1w
                            .multiplied_scalar(cosa - 1.0)
                            .added(&dnp2.multiplied_scalar(sina)),
                    );
                    d_poles[jj] = d_poles[jj].added(&t.multiplied_scalar(rayon));
                }

                let lambda = 1.0 / (2.0 * cosas2 * cosas2);
                let dlambda = (lambda * sinas2 * d_angle) / (cosas2 * nb_span as f64);

                for i in 1..=nb_span {
                    let jj = 2 * (i - 1);
                    let temp = GpVec::from_pnts(
                        center,
                        &poles[jj].translated_vec(&GpVec::from_pnts(center, &poles[jj + 2])),
                    );
                    poles[jj + 1] = center.translated_vec(&temp.multiplied_scalar(lambda));

                    let dtemp = d_poles[jj]
                        .added(&d_poles[jj + 2])
                        .subtracted(&d_center.multiplied_scalar(2.0));
                    d_poles[jj + 1] = temp
                        .multiplied_scalar(dlambda)
                        .added(&d_center.multiplied_scalar(1.0 - 2.0 * lambda))
                        .added(&dtemp.multiplied_scalar(lambda));
                }

                let dl = -sinas2 * d_angle / (2.0 * nb_span as f64);
                let mut i = 0;
                while i < upp {
                    weights[i] = 1.0;
                    weights[i + 1] = cosas2;
                    d_weights[i] = 0.0;
                    d_weights[i + 1] = dl;
                    i += 2;
                }
                Ok(())
            }
        }
    }

    /// `GeomFill::GetCircle(TConv, ns1, ns2, dn1w, dn2w, d2n1w, d2n2w, nplan,
    /// dnplan, d2nplan, pts1, pts2, tang1, tang2, Dtang1, Dtang2, Rayon, DRayon,
    /// D2Rayon, Center, DCenter, D2Center, Poles, DPoles, D2Poles, Weights,
    /// DWeights, D2Weights)` (`cxx:604-827`).
    #[allow(clippy::too_many_arguments)]
    pub fn get_circle_d2(
        t_conv: ParameterisationType,
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
    ) -> Result<(), String> {
        let (angle, sina, cosa) = arc_angle(ns1, ns2, nplan);
        let (d_angle, d2_angle) = if sina.abs() > cosa.abs() {
            let aux = dn1w.dot(ns2) + ns1.dot(dn2w);
            let da = -aux / sina;
            let d2a = -(d2n1w.dot(ns2) + 2.0 * dn1w.dot(dn2w) + ns1.dot(d2n2w)) / sina
                + aux * (dnplan.dot(&ns1.crossed(ns2))
                    + nplan.dot(&dn1w.crossed(ns2).added(&ns1.crossed(dn2w))))
                    / (sina * sina);
            (da, d2a)
        } else {
            let temp = dn1w.crossed(ns2).added(&ns1.crossed(dn2w));
            let da =
                (dnplan.dot(&ns1.crossed(ns2)) + nplan.dot(&temp)) / cosa;
            let d2a = (d2nplan.dot(&ns1.crossed(ns2))
                + 2.0 * dnplan.dot(&temp)
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
            (da, d2a)
        };

        let low = 0;
        let upp = poles.len() - 1;
        poles[low] = *pts1;
        poles[upp] = *pts2;
        weights[low] = 1.0;
        weights[upp] = 1.0;
        d_poles[low] = *tang1;
        d_poles[upp] = *tang2;
        d_weights[low] = 0.0;
        d_weights[upp] = 0.0;
        d2_poles[low] = *d_tang1;
        d2_poles[upp] = *d_tang2;
        d2_weights[low] = 0.0;
        d2_weights[upp] = 0.0;

        match t_conv {
            ParameterisationType::QuasiAngular => {
                let mut q = QuasiAngularConvertor::new();
                let (p, dp, d2p, w, dw, d2w) = q.section_d2(
                    pts1, tang1, d_tang1, center, d_center, d2_center, nplan, dnplan, d2nplan,
                    angle, d_angle, d2_angle,
                );
                copy_poles_weights(poles, weights, &p, &w)?;
                if d_poles.len() != dp.len()
                    || d2_poles.len() != d2p.len()
                    || d_weights.len() != dw.len()
                    || d2_weights.len() != d2w.len()
                {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                d_poles.copy_from_slice(&dp);
                d2_poles.copy_from_slice(&d2p);
                d_weights.copy_from_slice(&dw);
                d2_weights.copy_from_slice(&d2w);
                Ok(())
            }
            ParameterisationType::Polynomial => {
                let mut pc = PolynomialConvertor::new();
                let (p, dp, d2p) = pc.section_d2(
                    pts1, tang1, d_tang1, center, d_center, d2_center, nplan, dnplan, d2nplan,
                    angle, d_angle, d2_angle,
                );
                if poles.len() != p.len() {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                poles.copy_from_slice(&p);
                d_poles.copy_from_slice(&dp);
                d2_poles.copy_from_slice(&d2p);
                weights.fill(1.0);
                d_weights.fill(0.0);
                d2_weights.fill(0.0);
                Ok(())
            }
            _ => {
                // `cxx:754-825`: the rational expression with second
                // derivatives.
                let nb_span = (poles.len() - 1) / 2;
                if nb_span == 0 {
                    return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
                }
                let np2 = nplan.crossed(ns1);
                let dnp2 = dnplan.crossed(ns1).added(&nplan.crossed(dn1w));
                // `cxx:758` verbatim: `d2np2 = d2nplan x ns1 + nplan x dn2w
                // + 2 * dnplan x dn1w` (OCCT uses `dn2w`, the first derivative
                // of the second normal, not `d2n1w`).
                let d2np2 = d2nplan
                    .crossed(ns1)
                    .added(&nplan.crossed(dn2w))
                    .added(&dnplan.crossed(dn1w).multiplied_scalar(2.0));

                let alpha = angle / nb_span as f64;
                let cosas2 = (alpha / 2.0).cos();
                let sinas2 = (alpha / 2.0).sin();

                for i in 1..=(nb_span - 1) {
                    let jj = 2 * i;
                    let lambda = i as f64 * alpha;
                    let cosa = lambda.cos();
                    let sina = lambda.sin();
                    let temp = ns1
                        .multiplied_scalar(cosa - 1.0)
                        .added(&np2.multiplied_scalar(sina));
                    poles[jj] = pts1.translated_vec(&temp.multiplied_scalar(rayon));

                    d_poles[jj] = temp.multiplied_scalar(d_rayon).added(tang1);
                    let mut dtemp = ns1
                        .multiplied_scalar(-sina)
                        .added(&np2.multiplied_scalar(cosa));
                    let aux = i as f64 / nb_span as f64;
                    dtemp = dtemp.multiplied_scalar(aux * d_angle);
                    dtemp = dtemp.added(
                        &dn1w
                            .multiplied_scalar(cosa - 1.0)
                            .added(&dnp2.multiplied_scalar(sina)),
                    );
                    d_poles[jj] = d_poles[jj].added(&dtemp.multiplied_scalar(rayon));

                    let mut t = d_tang1
                        .multiplied_scalar(1.0)
                        .added(&temp.multiplied_scalar(d2_rayon))
                        .added(&dtemp.multiplied_scalar(2.0 * d_rayon));
                    let mut u = ns1
                        .multiplied_scalar(cosa - 1.0)
                        .added(&d2np2.multiplied_scalar(sina));
                    let mut v = ns1
                        .multiplied_scalar(-sina)
                        .added(&np2.multiplied_scalar(cosa));
                    u = u.added(&v.multiplied_scalar(aux * aux * d2_angle));
                    v = dn1w
                        .added(&np2)
                        .multiplied_scalar(-sina)
                        .added(&dnp2.multiplied_scalar(cosa))
                        .added(&ns1.multiplied_scalar(-cosa));
                    u = u.added(&v.multiplied_scalar(aux * d_angle));
                    t = t.added(&u.multiplied_scalar(rayon));
                    d2_poles[jj] = t;
                }

                let lambda = 1.0 / (2.0 * cosas2 * cosas2);
                let dlambda = (lambda * sinas2 * d_angle) / (cosas2 * nb_span as f64);
                let aux = sinas2 / cosas2;
                let d2lambda = (dlambda * aux * d_angle
                    + d2_angle * aux * lambda
                    + (1.0 + aux * aux) * (d_angle / (2.0 * nb_span as f64)) * d_angle * lambda)
                    / nb_span as f64;

                for i in 1..=nb_span {
                    let jj = 2 * (i - 1);
                    let temp = GpVec::from_pnts(
                        center,
                        &poles[jj].translated_vec(&GpVec::from_pnts(center, &poles[jj + 2])),
                    );
                    poles[jj + 1] = center.translated_vec(&temp.multiplied_scalar(lambda));

                    let dtemp = d_poles[jj]
                        .added(&d_poles[jj + 2])
                        .subtracted(&d_center.multiplied_scalar(2.0));
                    d_poles[jj + 1] = temp
                        .multiplied_scalar(dlambda)
                        .added(&dtemp.multiplied_scalar(lambda))
                        .added(d_center);

                    d2_poles[jj + 1] = d2_poles[jj]
                        .added(&d2_poles[jj + 2])
                        .multiplied_scalar(lambda)
                        .added(&temp.multiplied_scalar(d2lambda))
                        .added(&dtemp.multiplied_scalar(2.0 * dlambda))
                        .added(&d2_center.multiplied_scalar(1.0 - 2.0 * lambda));
                }

                let dl = -sinas2 * d_angle / (2.0 * nb_span as f64);
                let d2l =
                    -sinas2 * d2_angle / (2.0 * nb_span as f64) - cosas2 * (d_angle / (2.0 * nb_span as f64)).powi(2);
                let mut i = 0;
                while i < upp {
                    weights[i] = 1.0;
                    weights[i + 1] = cosas2;
                    d_weights[i] = 0.0;
                    d_weights[i + 1] = dl;
                    d2_weights[i] = 0.0;
                    d2_weights[i + 1] = d2l;
                    i += 2;
                }
                Ok(())
            }
        }
    }
}

/// `Cosa` / `Sina` / `Angle` of `GeomFill::GetCircle` (`cxx:369-398`).
///
/// Returns `(Angle, Sina, Cosa)` after the `[-1, 1]` clamp and the
/// `]-pi/2, 3pi/2]` re-framing.
fn arc_angle(ns1: &GpVec, ns2: &GpVec, nplan: &GpVec) -> (f64, f64, f64) {
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
    (angle, sina, cosa)
}

/// Copy the convertor's pole / weight vectors into the caller's arrays.
fn copy_poles_weights(
    poles: &mut [GpPnt],
    weights: &mut [f64],
    p: &[GpPnt],
    w: &[f64],
) -> Result<(), String> {
    if poles.len() != p.len() || weights.len() != w.len() {
        return Err("Standard_OutOfRange: GeomFill::GetCircle".into());
    }
    poles.copy_from_slice(p);
    weights.copy_from_slice(w);
    Ok(())
}
