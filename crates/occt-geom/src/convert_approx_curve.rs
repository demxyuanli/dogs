//! `GeomConvert_ApproxCurve` (TKGeomBase): approximation of a 3D curve by a
//! B-spline within a tolerance.
//!
//! Source: `GeomConvert_ApproxCurve.cxx` (`GeomConvert_ApproxCurve_Eval`
//! `cxx:31-106`, constructors `cxx:108-125`, `Approximate` `cxx:127-182`,
//! accessors `cxx:184-207`) and `GeomConvert_ApproxCurve.hxx:29-94`.
//!
//! The curve is fed to `AdvApprox_ApproxAFunction` as
//! `Num1DSS = 0, Num2DSS = 0, Num3DSS = 1` (`cxx:135-138`), with the C2 and
//! C3 break points of the curve as the recommended / preferred cutting points
//! of `AdvApprox_PrefAndRec` (`cxx:143-150`). `GeomConvert_ApproxCurve`
//! exposes no `AverageError()` in OCCT 8.0.0 (see the header); the underlying
//! `AdvApprox_ApproxAFunction` average error is available on
//! [`crate::adv_approx::ApproxAFunction3d::average_error`].
//!
//! Numerical correctness depends on the `occt-core` `PLib_JacobiPolynomial`
//! transformation table (`plib_jacobi_data.pxx`) being read exactly as OCCT
//! reads it (`PLib_JacobiPolynomial.cxx:316-372`); see the R2-6 report for the
//! out-of-scope reader defect observed while validating this port.

use crate::adv_approx::{ApproxAFunction3d, PrefAndRec};
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use occt_core::kernel::geomabs::Shape;

/// `GeomConvert_ApproxCurve` (`GeomConvert_ApproxCurve.hxx:29-94`).
pub struct GeomConvertApproxCurve {
    is_done: bool,
    has_result: bool,
    curve: Option<GeomBSplineCurve>,
    max_error: f64,
}

impl GeomConvertApproxCurve {
    /// `GeomConvert_ApproxCurve(const Handle(Geom_Curve)&, Tol3d, Order,
    /// MaxSegments, MaxDegree)` (`cxx:108-116`). OCCT wraps the curve in a
    /// `GeomAdaptor_Curve`; here the `Curve` trait plays that role.
    pub fn new(
        curve: &dyn Curve,
        tol3d: f64,
        order: Shape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        Self::approximate(curve, tol3d, order, max_segments, max_degree)
    }

    /// `GeomConvert_ApproxCurve::Approximate` (`cxx:127-182`).
    fn approximate(
        the_curve: &dyn Curve,
        the_tol3d: f64,
        the_order: Shape,
        the_max_segments: i32,
        the_max_degree: i32,
    ) -> Self {
        // `cxx:135-138`: Num1DSS = Num2DSS = 0, Num3DSS = 1, ThreeDTol = {Tol3d}.
        // `cxx:140-141`.
        let first = the_curve.first_parameter();
        let last = the_curve.last_parameter();

        // `cxx:143-148`: `Intervals(GeomAbs_C2)` / `Intervals(GeomAbs_C3)`.
        // `GeomAbs_Shape` codes: C0=0, C1=2, C2=4, C3=5, CN=6.
        let cut_pnts_c2 = the_curve.parameter_intervals(4);
        let cut_pnts_c3 = the_curve.parameter_intervals(5);

        // `cxx:150`: `AdvApprox_PrefAndRec CutTool(CutPnts_C2, CutPnts_C3)`
        // (default `Weight = 5`).
        let cut_tool = match PrefAndRec::new(cut_pnts_c2, cut_pnts_c3, 5.0) {
            Ok(cut_tool) => cut_tool,
            Err(()) => return Self::no_result(),
        };

        // `GeomConvert_ApproxCurve_Eval::Evaluate` (`cxx:55-106`): order 0 is
        // the point, 1 the first derivative, 2 the second; any other order has
        // error code 3 with a zero result (`cxx:101-104`).
        let eval = |u: f64, order: i32, out: &mut [f64]| -> i32 {
            match order {
                0 => {
                    let p = the_curve.d0(u);
                    out[0] = p.x();
                    out[1] = p.y();
                    out[2] = p.z();
                }
                1 => {
                    let (_, v) = the_curve.d1(u);
                    out[0] = v.x();
                    out[1] = v.y();
                    out[2] = v.z();
                }
                2 => {
                    let (_, _, v) = the_curve.d2(u);
                    out[0] = v.x();
                    out[1] = v.y();
                    out[2] = v.z();
                }
                _ => {
                    out[0] = 0.0;
                    out[1] = 0.0;
                    out[2] = 0.0;
                    return 3;
                }
            }
            0
        };

        // `cxx:155-167`: `AdvApprox_ApproxAFunction(Num1DSS, Num2DSS, Num3DSS,
        // OneDTol, TwoDTol, ThreeDTol, First, Last, theOrder, theMaxDegree,
        // theMaxSegments, ev, CutTool)`.
        let a_approx = ApproxAFunction3d::approx(
            first,
            last,
            the_order,
            the_max_degree,
            the_max_segments,
            the_tol3d,
            &cut_tool,
            &eval,
        );

        match a_approx {
            Ok(a_approx) => {
                // `cxx:169-170`.
                let is_done = a_approx.done;
                let has_result = a_approx.has_result;
                // `cxx:172-181`: the poles are packed into a
                // `Geom_BSplineCurve(Poles, Knots, Mults, Degree)` and
                // `myMaxError = aApprox.MaxError(3, 1)`.
                let (curve, max_error) = if has_result {
                    let curve = GeomBSplineCurve::from_poles_knots_mults(
                        a_approx.poles,
                        a_approx.knots,
                        a_approx.mults,
                        a_approx.degree.max(0) as usize,
                    )
                    .ok();
                    (curve, a_approx.max_error)
                } else {
                    (None, 0.0)
                };
                Self {
                    is_done,
                    has_result,
                    curve,
                    max_error,
                }
            }
            Err(()) => Self::no_result(),
        }
    }

    fn no_result() -> Self {
        Self {
            is_done: false,
            has_result: false,
            curve: None,
            max_error: 0.0,
        }
    }

    /// `GeomConvert_ApproxCurve::Curve` (`cxx:184-187`): the B-spline
    /// resulting from the approximation, or `None` when OCCT's handle is null.
    pub fn curve(&self) -> Option<&GeomBSplineCurve> {
        self.curve.as_ref()
    }

    /// `GeomConvert_ApproxCurve::IsDone` (`cxx:189-192`): true when the
    /// approximation was done **within** the required tolerance.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `GeomConvert_ApproxCurve::HasResult` (`cxx:194-197`): true when a
    /// result was produced, not necessarily within tolerance.
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `GeomConvert_ApproxCurve::MaxError` (`cxx:199-202`): greatest distance
    /// between a point of the source curve and the B-spline.
    pub fn max_error(&self) -> f64 {
        self.max_error
    }
}
