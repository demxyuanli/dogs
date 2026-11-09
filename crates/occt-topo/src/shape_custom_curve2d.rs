//! `ShapeCustom_Curve2d` (`ShapeCustom_Curve2d.hxx` /
//! `ShapeCustom_Curve2d.cxx`, TKShHealing).
//!
//! Only the two 2D helpers `BRepOffset_Tool::CheckBounds` needs are ported:
//! `IsLinear` (`cxx:51-105`) and `ConvertToLine2d` (`cxx:109-152`).
//! `SimplifyBSpline2d` (`cxx:158-196`) is still UNPORTED - its only reader in
//! OCCT is `ShapeCustom_BSplineRestriction`, which is not in this port either.

use occt_core::elib::clib2d::{line_value_ax2d, parameter_lin2d};
use occt_core::gp::{GpAx2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::{Curve2d, Geom2dLine};

/// Static helpers of `ShapeCustom_Curve2d`.
pub struct ShapeCustomCurve2d;

impl ShapeCustomCurve2d {
    /// `ShapeCustom_Curve2d::IsLinear(thePoles, tolerance, Deviation)`
    /// (`cxx:51-105`).
    ///
    /// `Some(deviation)` is OCCT's `true` with the `Deviation` out parameter
    /// written; `None` is `false`, in which case OCCT leaves `Deviation`
    /// untouched.
    ///
    /// The pole array is `NCollection_Array1<gp_Pnt2d>` (1-based) in OCCT; the
    /// port takes the same sequence as a slice.
    pub fn is_linear(poles: &[GpPnt2d], tolerance: f64) -> Option<f64> {
        let nb_poles = poles.len();
        if nb_poles < 2 {
            return None;
        }

        // `cxx:60-77`: the two farthest poles.
        let mut d_max = 0.0;
        let mut i_max1 = 0usize;
        let mut i_max2 = 0usize;
        for i in 0..nb_poles - 1 {
            for j in i + 1..nb_poles {
                let dist = poles[i].square_distance(&poles[j]);
                if dist > d_max {
                    d_max = dist;
                    i_max1 = i;
                    i_max2 = j;
                }
            }
        }

        // `cxx:77-82`: `dPreci = Precision::PConfusion()^2`.
        let d_preci = PCONFUSION * PCONFUSION;
        if d_max < d_preci {
            return None;
        }

        // `cxx:84-104`: every pole within `tol` of the line through the two
        // farthest ones.
        let tol2 = tolerance * tolerance;
        let a_vec = GpVec2d::new(
            poles[i_max2].x() - poles[i_max1].x(),
            poles[i_max2].y() - poles[i_max1].y(),
        );
        let a_dir = GpDir2d::from_vec2d(&a_vec).ok()?;
        let a_lin = GpLin2d::from_pnt_dir(poles[i_max1], a_dir);

        let mut a_max = 0.0;
        for pole in poles {
            let dist = a_lin.square_distance(pole);
            if dist > tol2 {
                return None;
            }
            if dist > a_max {
                a_max = dist;
            }
        }
        Some(a_max.sqrt())
    }

    /// `ShapeCustom_Curve2d::ConvertToLine2d(theCurve, c1, c2, theTolerance,
    /// cf, cl, theDeviation)` (`cxx:109-152`).
    ///
    /// `Some((line, cf, cl, deviation))` is a non-null `Geom2d_Line`; `None` is
    /// OCCT's null handle ("it is not a line"). Only the two `down_cast` arms
    /// of OCCT are reachable: a `Geom2d_BSplineCurve` then a
    /// `Geom2d_BezierCurve`, each fed to [`Self::is_linear`]; every other curve
    /// type answers the null handle exactly like `cxx:150-151`.
    pub fn convert_to_line2d(
        curve: &dyn Curve2d,
        c1: f64,
        c2: f64,
        tolerance: f64,
    ) -> Option<(Geom2dLine, f64, f64, f64)> {
        let p1 = curve.d0(c1);
        let p2 = curve.d0(c2);
        let d_preci = tolerance * tolerance;
        if p1.square_distance(&p2) < d_preci {
            return None; // it is not a line
        }

        // `cxx:128-147`: the pole array of the spline/Bezier, whichever the
        // dynamic type is.
        let poles = if curve.is_bspline2d() {
            curve.poles2d()?
        } else if curve.is_bezier2d() {
            curve.poles2d()?
        } else {
            return None;
        };
        let deviation = Self::is_linear(&poles, tolerance)?;

        // `GetLine` (`cxx:34-48`).
        let (a_lin, cf, cl) = Self::get_line(&p1, &p2, c1);
        Some((
            Geom2dLine::from_pnt_dir(a_lin.location(), *a_lin.direction()),
            cf,
            cl,
            deviation,
        ))
    }

    /// `GetLine(P1, P2, c1, cf, cl)` (`cxx:34-48`): the line through `P1` and
    /// `P2` relocated to `ElCLib::Value(c1, alin)`, with `cf`/`cl` the
    /// parameters of `P1`/`P2` on it.
    fn get_line(p1: &GpPnt2d, p2: &GpPnt2d, c1: f64) -> (GpLin2d, f64, f64) {
        let a_vec = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        let a_dir = GpDir2d::from_vec2d(&a_vec).unwrap_or_default();
        let mut a_lin = GpLin2d::from_pnt_dir(*p1, a_dir);
        let pos = *a_lin.position();
        a_lin.set_position(GpAx2d::new(line_value_ax2d(c1, &pos), *a_lin.direction()));
        let cf = parameter_lin2d(&a_lin, p1);
        let cl = parameter_lin2d(&a_lin, p2);
        (a_lin, cf, cl)
    }

    /// `Precision::Confusion()`, the tolerance `BRepOffset_Tool::CheckBounds`
    /// passes to [`Self::convert_to_line2d`] (`BRepOffset_Tool.cxx:3362-3368`).
    pub const CONFUSION: f64 = CONFUSION;
}
