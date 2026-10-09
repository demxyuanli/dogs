//! `Geom2dConvert_ApproxCurve` (TKGeomBase): approximation of a 2D curve by a
//! B-spline within a tolerance.
//!
//! Source: `Geom2dConvert_ApproxCurve.cxx` (`Geom2dConvert_ApproxCurve_Eval`
//! `cxx:31-107`, constructors `cxx:111-125`, `Approximate` `cxx:127-183`,
//! accessors `cxx:185-203`) and `Geom2dConvert_ApproxCurve.hxx:29-92`.
//!
//! The curve is fed to `AdvApprox_ApproxAFunction` as
//! `Num1DSS = 0, Num2DSS = 1, Num3DSS = 0` (`cxx:137`), with the C2 and C3
//! break points of the curve as the recommended / preferred cutting points of
//! `AdvApprox_PrefAndRec` (`cxx:145-151`). The break points are read through
//! [`adaptor_intervals`], the local `Geom2dAdaptor_Curve::NbIntervals` /
//! `Intervals` resolution (trimmed unwrap + offset `BaseS` mapping).

use occt_core::adv_approx::{ApproxAFunction2d, PrefAndRec};
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::gp::GpPnt2d;
use occt_core::kernel::geomabs::Shape;

use crate::bspline_curve::Geom2dBSplineCurve;
use crate::curve::Curve2d;

/// `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:285-288`): a
/// `Geom2d_TrimmedCurve` stores its basis (recursively) while
/// `FirstParameter`/`LastParameter` stay the trimmed range.
fn stored_curve(c: &dyn Curve2d) -> &dyn Curve2d {
    let mut cur = c;
    while let Some(basis) = cur.trimmed_basis() {
        cur = basis;
    }
    cur
}

/// `Geom2dAdaptor_Curve::NbIntervals` / `Intervals` (`Geom2dAdaptor_Curve.cxx`
/// `:409-486`, `:490-573`) for the requested shape: the stored curve is the
/// trimmed-unwrapped one, and an offset curve resolves its break points on a
/// basis adaptor with `BaseS = GeomAbs_C0 -> C1, C1 -> C2, C2 -> C3, else CN`
/// (`:453-480`, `:536-566`). OCCT then forces the two end points to the
/// offset's own range (`:564-565`); the port's break-point arrays are a
/// superset, which `AdvApprox_PrefAndRec` ignores when outside `[First, Last]`.
fn adaptor_intervals(c: &dyn Curve2d, shape: u8) -> Vec<f64> {
    let stored = stored_curve(c);
    match stored.offset_basis() {
        Some(basis) => {
            // `GeomAbs_Shape` codes: C0=0, C1=2, C2=4, C3=5, CN=6.
            let base_s = match shape {
                0 => 2,
                2 => 4,
                4 => 5,
                _ => 6,
            };
            adaptor_intervals(basis, base_s)
        }
        None => stored.parameter_intervals(shape),
    }
}

/// `Geom2dConvert_ApproxCurve` (`Geom2dConvert_ApproxCurve.hxx:29-92`).
pub struct Geom2dConvertApproxCurve {
    is_done: bool,
    has_result: bool,
    curve: Option<Geom2dBSplineCurve>,
    max_error: f64,
}

impl Geom2dConvertApproxCurve {
    /// `Geom2dConvert_ApproxCurve(const Handle(Geom2d_Curve)&, Tol2d, Order,
    /// MaxSegments, MaxDegree)` (`cxx:111-118`). OCCT wraps the curve in a
    /// `Geom2dAdaptor_Curve`; here the `Curve2d` trait plays that role.
    pub fn new(
        curve: &dyn Curve2d,
        tol2d: f64,
        order: Shape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        Self::approximate(curve, tol2d, order, max_segments, max_degree)
    }

    /// `Geom2dConvert_ApproxCurve::Approximate` (`cxx:127-183`).
    fn approximate(
        the_curve: &dyn Curve2d,
        the_tol2d: f64,
        the_order: Shape,
        the_max_segments: i32,
        the_max_degree: i32,
    ) -> Self {
        // `cxx:137-139`: Num1DSS = Num3DSS = 0, Num2DSS = 1, TwoDTol = {Tol2d}.
        // `cxx:141-142`.
        let first = the_curve.first_parameter();
        let last = the_curve.last_parameter();

        // `cxx:145-150`: `NbIntervals(GeomAbs_C2)` / `NbIntervals(GeomAbs_C3)`
        // then `Intervals(CutPnts_C2, GeomAbs_C2)` / `Intervals(CutPnts_C3,
        // GeomAbs_C3)`.
        let cut_pnts_c2 = adaptor_intervals(the_curve, 4);
        let cut_pnts_c3 = adaptor_intervals(the_curve, 5);

        // `cxx:151`: `AdvApprox_PrefAndRec CutTool(CutPnts_C2, CutPnts_C3)`
        // (default `Weight = 5`).
        let cut_tool = match PrefAndRec::new(cut_pnts_c2, cut_pnts_c3, 5.0) {
            Ok(cut_tool) => cut_tool,
            Err(()) => return Self::no_result(),
        };

        // `Geom2dConvert_ApproxCurve_Eval::Evaluate` (`cxx:55-107`): order 0 is
        // the point, 1 the first derivative, 2 the second; any other order has
        // error code 3 with a zero result (`cxx:101-104`).
        let eval = |u: f64, order: i32, out: &mut [f64]| -> i32 {
            match order {
                0 => {
                    let p = the_curve.d0(u);
                    out[0] = p.x();
                    out[1] = p.y();
                }
                1 => {
                    let (_, v) = the_curve.d1(u);
                    out[0] = v.x();
                    out[1] = v.y();
                }
                2 => {
                    let (_, _, v) = the_curve.d2(u);
                    out[0] = v.x();
                    out[1] = v.y();
                }
                _ => {
                    out[0] = 0.0;
                    out[1] = 0.0;
                    return 3;
                }
            }
            0
        };

        // `cxx:155-168`: `AdvApprox_ApproxAFunction(Num1DSS, Num2DSS, Num3DSS,
        // OneDTol, TwoDTol, ThreeDTol, First, Last, theOrder, theMaxDegree,
        // theMaxSegments, ev, CutTool)`.
        let a_approx = ApproxAFunction2d::approx(
            first,
            last,
            the_order,
            the_max_degree,
            the_max_segments,
            the_tol2d,
            &cut_tool,
            &eval,
        );

        match a_approx {
            Ok(a_approx) => {
                // `cxx:170-171`.
                let is_done = a_approx.done;
                let has_result = a_approx.has_result;
                // `cxx:172-183`: the `Poles2d` are packed into a
                // `Geom2d_BSplineCurve(Poles, Knots, Mults, Degree)` and
                // `myMaxError = aApprox.MaxError(2, 1)`.
                let (curve, max_error) = if has_result {
                    let curve = bspline_from_poles_knots_mults(
                        &a_approx.poles,
                        &a_approx.knots,
                        &a_approx.mults,
                        a_approx.degree,
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

    /// `Geom2dConvert_ApproxCurve::Curve` (`cxx:185-188`): the B-spline
    /// resulting from the approximation, or `None` when OCCT's handle is null.
    pub fn curve(&self) -> Option<&Geom2dBSplineCurve> {
        self.curve.as_ref()
    }

    /// `Geom2dConvert_ApproxCurve::IsDone` (`cxx:190-193`): true when the
    /// approximation was done **within** the required tolerance.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Geom2dConvert_ApproxCurve::HasResult` (`cxx:195-198`): true when a
    /// result was produced, not necessarily within tolerance.
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `Geom2dConvert_ApproxCurve::MaxError` (`cxx:200-203`): greatest distance
    /// between a point of the source curve and the B-spline.
    pub fn max_error(&self) -> f64 {
        self.max_error
    }
}

/// `Geom2d_BSplineCurve(Poles, Knots, Mults, Degree)`
/// (`Geom2dConvert_ApproxCurve.cxx:180`): the distinct knots/multiplicities are
/// expanded to the flat knot sequence the port stores.
fn bspline_from_poles_knots_mults(
    poles: &[GpPnt2d],
    knots: &[f64],
    mults: &[i32],
    degree: i32,
) -> Result<Geom2dBSplineCurve, &'static str> {
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    let flat = knot_sequence(knots, mults, degree);
    Geom2dBSplineCurve::from_flat(xs, ys, None, flat, degree.max(0) as usize, false)
}
