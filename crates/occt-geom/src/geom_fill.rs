//! `GeomFill_Profiler` and `GeomFill_Generator` for a ruled section pair.
//!
//! Source: `GeomFill_Profiler.cxx` (`UnifyByInsertingAllKnots` `cxx:29-74`,
//! `UnifyBySettingMiddleKnots` `cxx:78-109`, `AddCurve` `cxx:125-162`,
//! `Perform` `cxx:166-246`) and `GeomFill_Generator.cxx:32-82`.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::convert_approx_curve::GeomConvertApproxCurve;
use crate::convert_bspl::curve_to_bspline_curve;
use crate::curve::Curve;
use crate::surface::Surface;

/// `Geom_TrimmedCurve(C, U1, U2)` as `BRepFill` builds it (`BRepFill.cxx:244`).
///
/// The parameter is the basis parameter. This is not
/// `occt_geom::GeomTrimmedCurve`, which remaps the span onto `[0, 1]`.
pub struct ParameterTrim {
    basis: Arc<dyn Curve>,
    first: f64,
    last: f64,
}

impl ParameterTrim {
    pub fn new(basis: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { basis, first, last }
    }
}

impl Curve for ParameterTrim {
    fn d0(&self, u: f64) -> occt_core::gp::GpPnt {
        self.basis.d0(u)
    }
    fn d1(&self, u: f64) -> (occt_core::gp::GpPnt, occt_core::gp::GpVec) {
        self.basis.d1(u)
    }
    fn d2(&self, u: f64) -> (occt_core::gp::GpPnt, occt_core::gp::GpVec, occt_core::gp::GpVec) {
        self.basis.d2(u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }
    fn is_geom_trimmed(&self) -> bool {
        true
    }
    fn untrimmed_basis(&self) -> Option<(Arc<dyn Curve>, f64, f64)> {
        match self.basis.untrimmed_basis() {
            Some((inner, _, _)) => Some((inner, self.first, self.last)),
            None => Some((self.basis.clone(), self.first, self.last)),
        }
    }
    fn transform(&mut self, t: &occt_core::gp::GpTrsf) {
        let mut basis = self.basis.clone_dyn();
        basis.transform(t);
        self.basis = Arc::from(basis);
    }
    /// `Geom_TrimmedCurve::Reverse` (`Geom_TrimmedCurve.cxx:88-96`).
    fn reverse(&mut self) {
        let first = self.basis.reversed_parameter(self.last);
        let last = self.basis.reversed_parameter(self.first);
        let mut basis = self.basis.clone_dyn();
        basis.reverse();
        self.basis = Arc::from(basis);
        self.first = first;
        self.last = last;
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(ParameterTrim {
            basis: Arc::from(self.basis.clone_dyn()),
            first: self.first,
            last: self.last,
        })
    }
}

/// `GeomFill_Generator::Perform` (`GeomFill_Generator.cxx:32-82`) after
/// `GeomFill_Profiler::AddCurve` of each section.
pub fn generator_surface(sections: &[Arc<dyn Curve>]) -> Result<Arc<dyn Surface>, String> {
    let mut profiler = Profiler::default();
    for curve in sections {
        profiler.add_curve(curve.as_ref())?;
    }
    profiler.perform(PCONFUSION)?;
    profiler.surface()
}

#[derive(Default)]
struct Profiler {
    curves: Vec<GeomBSplineCurve>,
    /// `myIsPeriodic` starts true (`GeomFill_Profiler.cxx:113-116`).
    periodic: bool,
    done: bool,
}

impl Profiler {
    fn add_curve(&mut self, curve: &dyn Curve) -> Result<(), String> {
        // `GeomFill_Profiler::AddCurve` (`cxx:125-162`).
        let mut converted: Option<GeomBSplineCurve> = None;
        let kind = if curve.is_geom_trimmed() {
            curve.trim_basis()
        } else {
            None
        };
        let kind_curve: &dyn Curve = kind.as_deref().unwrap_or(curve);
        if is_conic(kind_curve) {
            // `GeomConvert_ApproxCurve(Curve, Confusion, C1, 16, 14)` (`cxx:136`).
            // The approximation sees the original argument, including a trim.
            let appr = GeomConvertApproxCurve::new(curve, CONFUSION, Shape::C1, 16, 14);
            if appr.has_result() {
                converted = appr.curve().cloned();
            }
        }
        let bspline = match converted {
            Some(curve) => curve,
            None => curve_to_bspline_curve(curve, ParameterisationType::TgtThetaOver2)
                .map_err(|_| "GeomConvert::CurveToBSplineCurve".to_string())?,
        };
        if self.periodic && !bspline.is_periodic() {
            self.periodic = false;
        }
        self.curves.push(bspline);
        Ok(())
    }

    fn perform(&mut self, ptol: f64) -> Result<(), String> {
        // `GeomFill_Profiler::Perform` (`cxx:166-246`).
        if self.curves.is_empty() {
            return Err("GeomFill_Profiler::Perform".into());
        }
        let mut degree = 0usize;
        let mut span = 0.0;
        let mut u_first = 0.0;
        let mut u_last = 0.0;
        let periodic = self.periodic;
        for curve in &mut self.curves {
            let (u1, u2) = knot_ends(curve);
            if !periodic && curve.is_periodic() {
                curve.set_not_periodic();
                curve
                    .segment(u1, u2, PCONFUSION)
                    .map_err(|_| "Geom_BSplineCurve::Segment")?;
            }
            degree = degree.max(curve.degree());
            if u2 - u1 > span {
                span = u2 - u1;
                u_first = u1;
                u_last = u2;
            }
        }
        for curve in &mut self.curves {
            curve
                .increase_degree(degree)
                .map_err(|_| "Geom_BSplineCurve::IncreaseDegree")?;
            let (mut knots, _) = curve.distinct_knots_and_mults();
            occt_core::bspl::knots::reparametrize(u_first, u_last, &mut knots);
            curve
                .set_knots(&knots)
                .map_err(|_| "Geom_BSplineCurve::SetKnots")?;
        }
        let mut copies: Vec<GeomBSplineCurve> = self.curves.iter().cloned().collect();
        unify_by_inserting_all_knots(&mut copies, ptol)?;
        let n0 = distinct_count(&copies[0]);
        let unified = copies.iter().skip(1).all(|c| distinct_count(c) == n0);
        if unified {
            self.curves = copies;
        } else {
            unify_by_setting_middle_knots(&mut self.curves)?;
        }
        self.done = true;
        Ok(())
    }

    fn surface(&self) -> Result<Arc<dyn Surface>, String> {
        // `GeomFill_Generator::Perform` (`GeomFill_Generator.cxx:38-81`).
        if !self.done || self.curves.is_empty() {
            return Err("GeomFill_Generator::Perform".into());
        }
        let nu = self.curves[0].poles.len();
        let nv = self.curves.len();
        let mut poles = vec![vec![occt_core::gp::GpPnt::zero(); nv]; nu];
        let mut weights = vec![vec![1.0; nv]; nu];
        let mut rational = false;
        let mut v_knots = Vec::with_capacity(nv);
        let mut v_mults = vec![1i32; nv];
        if nv >= 1 {
            v_mults[0] = 2;
            v_mults[nv - 1] = 2;
        }
        for (j, curve) in self.curves.iter().enumerate() {
            v_knots.push(j as f64);
            if curve.poles.len() != nu {
                return Err("GeomFill_Generator::Perform".into());
            }
            for i in 0..nu {
                poles[i][j] = curve.poles[i];
                let w = curve.weight((i + 1) as i32);
                weights[i][j] = w;
                if (w - 1.0).abs() > 0.0 {
                    rational = true;
                }
            }
        }
        let (u_knots, u_mults) = self.curves[0].distinct_knots_and_mults();
        let mut surface = GeomBSplineSurface::from_poles_knots_mults(
            poles,
            u_knots,
            v_knots,
            u_mults,
            v_mults,
            self.curves[0].degree(),
            1,
        )?;
        if rational {
            surface.weights = Some(weights);
        }
        Ok(Arc::new(surface))
    }
}

fn is_conic(curve: &dyn Curve) -> bool {
    curve.gp_circ().is_some()
        || curve.gp_ellipse().is_some()
        || curve.gp_hyperbola().is_some()
        || curve.gp_parabola().is_some()
}

/// `Knot(FirstUKnotIndex)` / `Knot(LastUKnotIndex)` (`Geom_BSplineCurve_1.cxx:334-344`).
fn knot_ends(curve: &GeomBSplineCurve) -> (f64, f64) {
    let (knots, mults) = curve.distinct_knots_and_mults();
    if knots.is_empty() {
        return (curve.first_parameter(), curve.last_parameter());
    }
    let (first, last) = if curve.is_periodic() {
        (1usize, knots.len())
    } else {
        (
            occt_core::bspl::locate::first_u_knot_index(curve.degree() as i32, &mults).max(1)
                as usize,
            occt_core::bspl::locate::last_u_knot_index(curve.degree() as i32, &mults).max(1)
                as usize,
        )
    };
    let first = first.min(knots.len());
    let last = last.min(knots.len()).max(first);
    (knots[first - 1], knots[last - 1])
}

fn distinct_count(curve: &GeomBSplineCurve) -> usize {
    curve.distinct_knots_and_mults().0.len()
}

/// `UnifyByInsertingAllKnots` (`GeomFill_Profiler.cxx:29-74`).
fn unify_by_inserting_all_knots(curves: &mut [GeomBSplineCurve], ptol: f64) -> Result<(), String> {
    if curves.is_empty() {
        return Ok(());
    }
    for i in 1..curves.len() {
        let (knots, mults) = curves[i].distinct_knots_and_mults();
        curves[0]
            .insert_knots(&knots, Some(&mults), ptol, false)
            .map_err(|_| "Geom_BSplineCurve::InsertKnots")?;
    }
    let (new_knots, new_mults) = curves[0].distinct_knots_and_mults();
    for curve in curves.iter_mut().skip(1) {
        curve
            .insert_knots(&new_knots, Some(&new_mults), ptol, false)
            .map_err(|_| "Geom_BSplineCurve::InsertKnots")?;
    }
    for curve in curves.iter_mut() {
        if !curve.is_rational() {
            continue;
        }
        let np = curve.poles.len();
        if np == 0 {
            continue;
        }
        let mut sigma = 0.0;
        for j in 1..=np {
            sigma += curve.weight(j as i32);
        }
        sigma /= np as f64;
        for j in 1..=np {
            let w = curve.weight(j as i32) / sigma;
            curve.set_weight(j - 1, w);
        }
    }
    Ok(())
}

/// `UnifyBySettingMiddleKnots` (`GeomFill_Profiler.cxx:78-109`).
fn unify_by_setting_middle_knots(curves: &mut [GeomBSplineCurve]) -> Result<(), String> {
    if curves.is_empty() {
        return Ok(());
    }
    let (knots, _) = curves[0].distinct_knots_and_mults();
    let nb = knots.len();
    if nb == 0 {
        return Err("GeomFill_Profiler::KnotsAndMults".into());
    }
    let mut middle = knots;
    for j in 1..nb - 1 {
        let mut sum = 0.0;
        for curve in curves.iter() {
            let (ck, _) = curve.distinct_knots_and_mults();
            if j >= ck.len() {
                return Err("GeomFill_Profiler::KnotsAndMults".into());
            }
            sum += ck[j];
        }
        middle[j] = sum / curves.len() as f64;
    }
    for curve in curves.iter_mut() {
        curve
            .set_knots(&middle)
            .map_err(|_| "Geom_BSplineCurve::SetKnots")?;
    }
    Ok(())
}
