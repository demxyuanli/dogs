//! `Convert_CompBezierCurvesToBSplineCurve` /
//! `Convert_CompBezierCurves2dToBSplineCurve2d`.
//!
//! Source: `FoundationClasses/TKMath/Convert/Convert_CompBezierCurvesToBSplineCurveBase.hxx:28-220`
//! (both public classes are thin wrappers of this template,
//! `Convert_CompBezierCurvesToBSplineCurve.cxx:18-25` and
//! `Convert_CompBezierCurves2dToBSplineCurve2d.cxx:18-25`).
//!
//! Translation notes:
//! - `AddCurve` + `Perform` + accessors are kept as the same state machine; the
//!   OCCT `NCollection_Sequence` members become `Vec`.
//! - `BSplCLib::IncreaseDegree(deg, poles, NoWeights, points, NoWeights)` is the
//!   Bezier syntax (`BSplCLib_CurveComputation.pxx:2080-2098`), i.e. flat
//!   Bezier knots `{0, 1}` with multiplicities `{deg+1, deg+1}`.
//! - `Epsilon(X)` is `RealEpsilon() * |X|`. Only the 3d junction arm carries the
//!   `aCurveKnVals(i-1) * aLambda > 10 * Epsilon(aDet)` guard
//!   (`Base.hxx:125-127`); the 2d arm (`Base.hxx:135-140`) has none.
//! - The template's `PointType` / `VecType` split is kept by the `is_parallel`
//!   closure the caller supplies (`gp_Vec` vs `gp_Vec2d`).

use crate::bspl::increase_degree::increase_degree;
use crate::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_foundation::precision::{REAL_EPSILON, RESOLUTION};

/// `BSplCLib::IncreaseDegree` in the Bezier syntax, 3d
/// (`BSplCLib_CurveComputation.pxx:2080-2098`).
pub fn increase_degree_bezier(new_degree: i32, poles: &[GpPnt]) -> Vec<GpPnt> {
    let degree = poles.len() as i32 - 1;
    if new_degree <= degree || poles.is_empty() {
        return poles.to_vec();
    }
    let knots = [0.0f64, 1.0f64];
    let mults = [degree + 1, degree + 1];
    increase_degree(
        degree,
        new_degree,
        false,
        poles,
        None,
        &knots,
        &mults,
        (new_degree + 1) as usize,
        2,
    )
    .poles
}

/// `BSplCLib::IncreaseDegree` in the Bezier syntax, 2d.
///
/// `BSplCLib` is not instantiated for `gp_Pnt2d` in this port; degree elevation
/// is linear in the coordinates, so the poles are lifted to 3d with `z == 0`,
/// elevated, and projected back.
pub fn increase_degree_bezier_2d(new_degree: i32, poles: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let degree = poles.len() as i32 - 1;
    if new_degree <= degree || poles.is_empty() {
        return poles.to_vec();
    }
    let lifted: Vec<GpPnt> = poles
        .iter()
        .map(|p| GpPnt::new(p.x(), p.y(), 0.0))
        .collect();
    increase_degree_bezier(new_degree, &lifted)
        .iter()
        .map(|p| GpPnt2d::new(p.x(), p.y()))
        .collect()
}

/// Result of `Perform()` (`Base.hxx:50-174`).
#[derive(Debug, Clone)]
pub struct CompBezierCurvesResult<P> {
    /// `Degree()`.
    pub degree: i32,
    /// `Poles()`.
    pub poles: Vec<P>,
    /// `Knots()`.
    pub knots: Vec<f64>,
    /// `KnotsAndMults()` multiplicities.
    pub mults: Vec<i32>,
}

/// Shared `Perform()` body (`Base.hxx:50-174`).
///
/// `is_3d` selects the `gp_Pnt` arm of the `if constexpr` at `Base.hxx:121`.
fn perform<F>(
    sequence: &[Vec<F>],
    is_3d: bool,
    to_pnt: impl Fn(&F) -> GpPnt,
    from_pnt: impl Fn(&[f64; 3]) -> F,
    elevate: impl Fn(i32, &[GpPnt]) -> Vec<GpPnt>,
    is_parallel: impl Fn([f64; 3], [f64; 3]) -> bool,
) -> CompBezierCurvesResult<F>
where
    F: Clone,
{
    let mut out: CompBezierCurvesResult<F> = CompBezierCurvesResult {
        degree: 0,
        poles: Vec::new(),
        knots: Vec::new(),
        mults: Vec::new(),
    };
    if sequence.is_empty() {
        return out;
    }
    let nbr_curv = sequence.len();
    let mut curve_kn_vals = vec![0.0f64; nbr_curv];

    let mut degree = 0i32;
    for curve in sequence {
        degree = degree.max(curve.len() as i32 - 1);
    }
    out.degree = degree;
    let max_degree = degree;

    let mut det = 0.0f64;
    let mut p1 = [0.0f64; 3];
    let mut points: Vec<GpPnt> = Vec::new();

    for (i, curve) in sequence.iter().enumerate() {
        // 1- Raise the Bezier curve to the maximum degree (`Base.hxx:77-91`).
        let src: Vec<GpPnt> = curve.iter().map(|p| to_pnt(p)).collect();
        points = elevate(max_degree, &src);

        // 2- Process the node of junction between 2 Bezier curves.
        if i == 0 {
            for j in 0..max_degree as usize {
                out.poles
                    .push(from_pnt(&[points[j].x(), points[j].y(), points[j].z()]));
            }
            curve_kn_vals[0] = 1.0;
            out.mults.push(max_degree + 1);
            det = 1.0;
        }

        if i != 0 {
            let a_p2 = [points[0].x(), points[0].y(), points[0].z()];
            let a_p3 = [points[1].x(), points[1].y(), points[1].z()];
            let v1 = [a_p2[0] - p1[0], a_p2[1] - p1[1], a_p2[2] - p1[2]];
            let v2 = [a_p3[0] - a_p2[0], a_p3[1] - a_p2[1], a_p3[2] - a_p2[2]];
            let d1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
            let d2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
            if max_degree > 1 && d1 > RESOLUTION && d2 > RESOLUTION && is_parallel(v1, v2) {
                let lambda = (d2 / d1).sqrt();
                let guarded = if is_3d {
                    curve_kn_vals[i - 1] * lambda > 10.0 * (REAL_EPSILON * det.abs())
                } else {
                    true
                };
                if guarded {
                    out.mults.push(max_degree - 1);
                    curve_kn_vals[i] = curve_kn_vals[i - 1] * lambda;
                } else {
                    out.poles
                        .push(from_pnt(&[points[0].x(), points[0].y(), points[0].z()]));
                    out.mults.push(max_degree);
                    curve_kn_vals[i] = 1.0;
                }
            } else {
                out.poles
                    .push(from_pnt(&[points[0].x(), points[0].y(), points[0].z()]));
                out.mults.push(max_degree);
                curve_kn_vals[i] = 1.0;
            }
            det += curve_kn_vals[i];

            // Store the poles (`Base.hxx:150-154`).
            for j in 1..max_degree as usize {
                out.poles
                    .push(from_pnt(&[points[j].x(), points[j].y(), points[j].z()]));
            }
        }

        if i == nbr_curv - 1 {
            let last = &points[max_degree as usize];
            out.poles.push(from_pnt(&[last.x(), last.y(), last.z()]));
            out.mults.push(max_degree + 1);
        }
        let pm = &points[max_degree as usize];
        p1 = [pm.x(), pm.y(), pm.z()];
    }

    // Correct nodal values to make them vary within [0.,1.] (`Base.hxx:166-173`).
    out.knots.push(0.0);
    for i in 1..nbr_curv {
        let prev = out.knots[i - 1];
        out.knots.push(prev + curve_kn_vals[i - 1] / det);
    }
    out.knots.push(1.0);
    out
}

/// `Convert_CompBezierCurvesToBSplineCurve` (`Base.hxx` instantiated with
/// `gp_Pnt` / `gp_Vec`).
#[derive(Debug, Clone)]
pub struct CompBezierCurves {
    sequence: Vec<Vec<GpPnt>>,
    angular: f64,
    result: CompBezierCurvesResult<GpPnt>,
}

impl CompBezierCurves {
    /// `Convert_CompBezierCurvesToBSplineCurve(theAngularTolerance = 1.0e-4)`.
    pub fn new(angular_tolerance: f64) -> Self {
        Self {
            sequence: Vec::new(),
            angular: angular_tolerance,
            result: CompBezierCurvesResult {
                degree: 0,
                poles: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
            },
        }
    }

    /// `AddCurve(thePoles)` (`Base.hxx:44-47`).
    pub fn add_curve(&mut self, poles: &[GpPnt]) {
        self.sequence.push(poles.to_vec());
    }

    /// `Perform()` (`Base.hxx:50-174`).
    pub fn perform(&mut self) {
        let angular = self.angular;
        let seq = self.sequence.clone();
        self.result = perform(
            &seq,
            true,
            |p| *p,
            |c| GpPnt::new(c[0], c[1], c[2]),
            increase_degree_bezier,
            |v1, v2| {
                GpVec::new(v1[0], v1[1], v1[2])
                    .is_parallel_ang(&GpVec::new(v2[0], v2[1], v2[2]), angular)
            },
        );
    }

    /// `Degree()`.
    pub fn degree(&self) -> i32 {
        self.result.degree
    }

    /// `NbPoles()`.
    pub fn nb_poles(&self) -> usize {
        self.result.poles.len()
    }

    /// `Poles()`.
    pub fn poles(&self) -> &[GpPnt] {
        &self.result.poles
    }

    /// `NbKnots()`.
    pub fn nb_knots(&self) -> usize {
        self.result.knots.len()
    }

    /// `KnotsAndMults()`.
    pub fn knots_and_mults(&self) -> (&[f64], &[i32]) {
        (&self.result.knots, &self.result.mults)
    }
}

/// `Convert_CompBezierCurves2dToBSplineCurve2d` (`Base.hxx` instantiated with
/// `gp_Pnt2d` / `gp_Vec2d`).
#[derive(Debug, Clone)]
pub struct CompBezierCurves2d {
    sequence: Vec<Vec<GpPnt2d>>,
    angular: f64,
    result: CompBezierCurvesResult<GpPnt2d>,
}

impl CompBezierCurves2d {
    /// `Convert_CompBezierCurves2dToBSplineCurve2d(theAngularTolerance = 1.0e-4)`.
    pub fn new(angular_tolerance: f64) -> Self {
        Self {
            sequence: Vec::new(),
            angular: angular_tolerance,
            result: CompBezierCurvesResult {
                degree: 0,
                poles: Vec::new(),
                knots: Vec::new(),
                mults: Vec::new(),
            },
        }
    }

    /// `AddCurve(thePoles)` (`Base.hxx:44-47`).
    pub fn add_curve(&mut self, poles: &[GpPnt2d]) {
        self.sequence.push(poles.to_vec());
    }

    /// `Perform()` (`Base.hxx:50-174`), 2d arm (`Base.hxx:135-140`).
    pub fn perform(&mut self) {
        let angular = self.angular;
        let seq = self.sequence.clone();
        self.result = perform(
            &seq,
            false,
            |p: &GpPnt2d| GpPnt::new(p.x(), p.y(), 0.0),
            |c| GpPnt2d::new(c[0], c[1]),
            increase_degree_bezier,
            |v1, v2| GpVec2d::new(v1[0], v1[1]).is_parallel(&GpVec2d::new(v2[0], v2[1]), angular),
        );
    }

    /// `Degree()`.
    pub fn degree(&self) -> i32 {
        self.result.degree
    }

    /// `NbPoles()`.
    pub fn nb_poles(&self) -> usize {
        self.result.poles.len()
    }

    /// `Poles()`.
    pub fn poles(&self) -> &[GpPnt2d] {
        &self.result.poles
    }

    /// `NbKnots()`.
    pub fn nb_knots(&self) -> usize {
        self.result.knots.len()
    }

    /// `KnotsAndMults()`.
    pub fn knots_and_mults(&self) -> (&[f64], &[i32]) {
        (&self.result.knots, &self.result.mults)
    }
}
