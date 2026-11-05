//! `MAT2d_Tool2d::TrimBisector` (both overloads, `MAT2d_Tool2d.cxx:490-557`),
//! `BisecFusion` (`cxx:1186-1249`), the static `SetTrim` (`cxx:1289-1327`) and
//! `Domain` (`cxx:1332-1380`).
//!
//! A bisector is parameterised by the distance of its points to the two items
//! it separates, so cutting it by a parameter, by a point of the locus or by
//! the first crossing of a connexion half-line are all the same trim. `Domain`
//! bounds the semi-infinite bisectors (the half-lines and the conic branches
//! whose parameter range is `Precision::Infinite()`) so that they can be
//! intersected.

use std::sync::Arc;

use occt_core::gp::GpVec2d;
use occt_core::intres2d::IntRes2dDomain;
use occt_core::precision::INFINITE;
use occt_geom2d::bisector::{bisector_parameter, BisectorBisec};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use super::{is_infinite, Mat2dTool2d, TOLCONF};
use crate::mat::MatBisector;

impl Mat2dTool2d {
    /// `TrimBisector(abisector)` (`cxx:492-515`): the bisector is cut at the
    /// `FirstParameter` of the source `MAT_Bisector`.
    pub fn trim_bisector(&mut self, abisector: &MatBisector) -> bool {
        let mut param = abisector.first_parameter();
        let index = abisector.bisector_number();

        let bisector = self
            .change_geom_bis(index)
            .change_value()
            .expect("MAT2d_Tool2d::TrimBisector: null thebisector");

        // `cxx:502-509`.
        if bisector.basis_curve().is_periodic() && is_infinite(param) {
            param = bisector.first_parameter() + 2.0 * std::f64::consts::PI;
        }
        let last = bisector.basis_curve().last_parameter();
        if param > last {
            param = last;
        }
        if bisector.first_parameter() == param {
            return false;
        }

        let first = bisector.first_parameter();
        bisector.set_trim_range(first, param);
        true
    }

    /// `TrimBisector(abisector, apoint)` (`cxx:522-555`): the bisector is cut at
    /// the parameter of the point `apoint` of `theGeomPnts`.
    pub fn trim_bisector_by_point(&mut self, abisector: &MatBisector, apoint: i32) -> bool {
        let point = self.geom_pnt(apoint);
        let index = abisector.bisector_number();

        let bisector = self
            .change_geom_bis(index)
            .change_value()
            .expect("MAT2d_Tool2d::TrimBisector: null thebisector");

        // `cxx:530`: `Bis->Parameter(GeomPnt(apoint))` on the basis curve.
        let basis = Arc::clone(bisector.basis_curve());
        let mut param = bisector_parameter(&*basis, &point);

        // `cxx:532-546`.
        if bisector.basis_curve().is_periodic() && bisector.first_parameter() > param {
            param += 2.0 * std::f64::consts::PI;
        }
        if bisector.first_parameter() >= param {
            return false;
        }
        if bisector.last_parameter() < param {
            return false;
        }

        let first = bisector.first_parameter();
        bisector.set_trim_range(first, param);
        true
    }

    /// `BisecFusion(I1, I2)` (`cxx:1186-1249`): a curve/curve bisector cannot be
    /// trimmed past its start point, so the fusion rebuilds it from the last
    /// point of `I2` and then cuts it at the last point of `I1`.
    pub fn bisec_fusion(&mut self, i1: i32, i2: i32) {
        let bisector1 = self
            .geom_bis(i1)
            .value()
            .expect("MAT2d_Tool2d::BisecFusion: null thebisector")
            .clone();
        let bisector2 = self
            .geom_bis(i2)
            .value()
            .expect("MAT2d_Tool2d::BisecFusion: null thebisector")
            .clone();

        let mut uf1 = bisector1.first_parameter();
        let mut ul1 = bisector1.last_parameter();
        let tolerance = TOLCONF;

        if let Some(bcc1) = bisector1.basis_curve().bisector_bisec_cc() {
            // `cxx:1199-1230`.
            let p2 = bisector2.d0(bisector2.last_parameter());
            let p1 = bisector1.d0(bisector1.last_parameter());
            let vbid = GpVec2d::new(1.0, 0.0);

            let curve2 = Arc::clone(
                bcc1
                    .curve(2)
                    .expect("MAT2d_Tool2d::BisecFusion: missing Curve(2)"),
            );
            let curve1 = Arc::clone(
                bcc1
                    .curve(1)
                    .expect("MAT2d_Tool2d::BisecFusion: missing Curve(1)"),
            );

            let mut bis = BisectorBisec::new();
            bis.perform_curve_curve(
                &*curve2,
                &*curve1,
                &p2,
                &vbid,
                &vbid,
                self.direction(),
                self.join_type(),
                tolerance,
                false,
            )
            .expect("MAT2d_Tool2d::BisecFusion: StdFail_NotDone");

            // `cxx:1226-1229`: the new trim runs from the new bisector's own
            // start to the parameter of `P1` on it.
            {
                let new_bisector = bis
                    .value()
                    .expect("MAT2d_Tool2d::BisecFusion: null thebisector");
                let new_bcc = new_bisector
                    .basis_curve()
                    .bisector_bisec_cc()
                    .expect("MAT2d_Tool2d::BisecFusion: not a Bisector_BisecCC");
                uf1 = new_bcc.first_parameter();
                ul1 = new_bcc.parameter(&p1);
            }
            bis.change_value()
                .expect("MAT2d_Tool2d::BisecFusion: null thebisector")
                .set_trim_range(uf1, ul1);
            self.the_geom_bisectors.insert(i1, bis);
        } else {
            // `cxx:1232-1248`: an analytic bisector keeps its curve, only its
            // parameter range is shifted by the length of `I2`.
            let du = bisector2.last_parameter() - bisector2.first_parameter();
            uf1 -= du;

            if let Some(mut ana) = bisector1.basis_curve().bisector_bisec_ana().cloned() {
                ana.set_trim(uf1, ul1);
                self.change_geom_bis(i1)
                    .change_value()
                    .expect("MAT2d_Tool2d::BisecFusion: null thebisector")
                    .set_basis_curve(Arc::new(ana));
            }
            self.change_geom_bis(i1)
                .change_value()
                .expect("MAT2d_Tool2d::BisecFusion: null thebisector")
                .set_trim_range(uf1, ul1);
        }
    }
}

/// `SetTrim(Bis, Line1)` (`cxx:1289-1327`): the bisector is cut at its first
/// crossing of the connexion half-line `Line1`.
pub(crate) fn set_trim(bis: &mut BisectorBisec, line: &Geom2dTrimmedCurve) {
    let tolerance = TOLCONF;
    let bisector = bis
        .change_value()
        .expect("MAT2d_Tool2d::SetTrim: null thebisector");

    // `cxx:1296-1319`.
    let domain1 = domain(bisector, tolerance);
    let ub1 = bisector.first_parameter();
    let ub2 = bisector.last_parameter();
    let first_point = bisector.d0(ub1);
    let mut u_trim = INFINITE;

    let mut intersect = Geom2dIntGInter::new();
    let domain2 = intersect.compute_domain(line, tolerance.max(tolerance));
    intersect.perform(bisector, &domain1, line, &domain2, tolerance, tolerance);

    if intersect.is_done() && intersect.nb_points() > 0 {
        for i in 1..=intersect.nb_points() {
            let point = intersect.point(i);
            let distance = first_point.distance(point.value());
            if distance > 10.0 * tolerance && point.param_on_first() < u_trim {
                u_trim = point.param_on_first();
            }
        }
    }

    // `cxx:1323-1326`.
    if u_trim < ub2 && u_trim > ub1 {
        bisector.set_trim_range(ub1, u_trim);
    }
}

/// `Domain(Bisector1, Tolerance)` (`cxx:1332-1380`): the intersection domain of
/// a bisector, whose last parameter is bounded when it is semi-infinite.
pub(crate) fn domain(bisector: &Geom2dTrimmedCurve, tolerance: f64) -> IntRes2dDomain {
    let param1 = bisector.first_parameter();
    let mut param2 = bisector.last_parameter();

    // `cxx:1335-1368`.
    if param2 > 10000.0 {
        param2 = 10000.0;
        let basis = bisector.basis_curve();
        let basis_curve: &dyn Curve2d = match basis.bisector_bisec_ana() {
            Some(ana) => match ana.geom2d_curve() {
                Some(curve) => &**curve,
                None => &**basis,
            },
            None => &**basis,
        };
        let limit = 50000.0;
        if let Some(parabola) = basis_curve.gp_parab2d() {
            let val1 = (limit * parabola.focal).sqrt();
            let val2 = (limit * limit).sqrt();
            param2 = if val1 <= val2 { val1 } else { val2 };
        } else if let Some(hyperbola) = basis_curve.gp_hypr2d() {
            let valu1 = limit / hyperbola.major_radius;
            let valu2 = limit / hyperbola.minor_radius;
            let val1 = (valu1 + (valu1 * valu1 - 1.0).sqrt()).ln();
            let val2 = (valu2 + (valu2 * valu2 + 1.0).sqrt()).ln();
            param2 = if val1 <= val2 { val1 } else { val2 };
        }
    }

    let mut result = IntRes2dDomain::bounded(        &bisector.d0(param1),
        param1,
        tolerance,
        &bisector.d0(param2),
        param2,
        tolerance,
    );
    if bisector.basis_curve().is_periodic() {
        result.set_equivalent_parameters(0.0, 2.0 * std::f64::consts::PI);
    }
    result
}
