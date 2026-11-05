//! `BRepFill_ApproxSeewing` — the B-spline approximation of a
//! [`MultiLine`](crate::brep_fill_multi_line::MultiLine) used by
//! `BRepFill_TrimSurfaceTool::Project` when the sawing line is not a
//! particular case.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_ApproxSeewing.cxx:17-245`
//! and `.hxx:25-63` (OCCT 8.0.0).
//!
//! The body is `Perform` (`cxx:45-223`) in three steps:
//!
//! 1. `cxx:49-68` — walk `2 * NbPoints` samples of the sawing line at a
//!    constant step of the bisector parameter and accumulate the 3D chord
//!    length, storing `(length, parameter)` pairs in `LP`.
//! 2. `cxx:72-142` — build an `AppDef_MultiLine` of `NbPoints` columns placed
//!    at equal chord length (`Value3dOnF1OnF2` of the interpolated parameter),
//!    then run `AppDef_Compute` on it.
//! 3. `cxx:146-166` / `cxx:169-223` — either the polygon fallback (one
//!    degree-1 B-spline through the same columns) when the approximation
//!    produced no multi-curve, or a B-spline of the approximation's poles,
//!    knots, multiplicities and degree; both feed one 3D curve and two
//!    p-curves (`myCurve` / `myPCurve1` / `myPCurve2`).
//!
//! `AppDef_Compute` is `Approx_ComputeLine.gxx` instantiated on
//! `AppDef_MultiLine` / `AppDef_MyLineTool`, ported as
//! [`ComputeLine`](occt_geom::approx::compute_line::ComputeLine); the
//! constructor `AppDef_Compute(const AppDef_MultiLine& Line)` with its default
//! arguments is `degreemin = 4`, `degreemax = 8`, `Tolerance3d = 1.0e-3`,
//! `Tolerance2d = 1.0e-6`, `NbIterations = 5`, `cutting = true`,
//! `parametrization = Approx_ChordLength`, `Squares = false`
//! (`AppDef_Compute.hxx:56-66`).
//!
//! UNPORTED: nothing of `cxx:45-223` is dropped; the `DUMP_ML` /
//! `OCCT_DEBUG` blocks (`cxx:88-96`, `cxx:113-118`, `cxx:139-143`,
//! `cxx:150-154`) only print.

use std::sync::Arc;

use occt_core::approx::ApproxParametrizationType;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::app_def::multi_line::MultiLine as AppDefMultiLine;
use occt_geom::app_def::MultiPointConstraint;
use occt_geom::approx::compute_line::ComputeLine;
use occt_geom::curve::Curve;
use occt_geom::GeomBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::brep_fill_multi_line::{Error as MultiLineError, MultiLine};

/// `NbPoints` of `Perform` (`cxx:56`): the number of columns of the
/// approximation and, with `2 * NbPoints` samples, the length of `LP`.
const NB_POINTS: usize = 50;

/// Failure modes of [`ApproxSeewing`].
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// One of the `BRepFill_MultiLine` evaluations of `cxx:60`, `cxx:75`,
    /// `cxx:84`, `cxx:153` failed.
    MultiLine(MultiLineError),
    /// The pole / knot arrays of `cxx:196-200` could not form a B-spline
    /// (`Geom_BSplineCurve(Poles, Knots, Mults, Degree)` raised in OCCT).
    SplineConstruction,
    /// `StdFail_NotDone_Raise_if(!myIsDone, ...)` of `Curve` / `CurveOnF1` /
    /// `CurveOnF2` (`cxx:238`, `cxx:248`, `cxx:258`).
    NotDone,
}

impl From<MultiLineError> for Error {
    fn from(e: MultiLineError) -> Self {
        Error::MultiLine(e)
    }
}

/// `BRepFill_ApproxSeewing` (`hxx:33-63`).
#[derive(Clone)]
pub struct ApproxSeewing {
    /// `myML` (`hxx:...`): the sawing line being approximated.
    ml: MultiLine,
    /// `myCurve` (`hxx:...`).
    curve: Option<Arc<dyn Curve>>,
    /// `myPCurve1`.
    pcurve1: Option<Arc<dyn Curve2d>>,
    /// `myPCurve2`.
    pcurve2: Option<Arc<dyn Curve2d>>,
    /// `myIsDone`.
    is_done: bool,
}

impl Default for ApproxSeewing {
    fn default() -> Self {
        Self::new()
    }
}

impl ApproxSeewing {
    /// `BRepFill_ApproxSeewing()` (`cxx:30-33`): `myIsDone = false`, `myML` is
    /// a default `BRepFill_MultiLine`.
    pub fn new() -> Self {
        Self {
            ml: MultiLine::new(),
            curve: None,
            pcurve1: None,
            pcurve2: None,
            is_done: false,
        }
    }

    /// `BRepFill_ApproxSeewing(const BRepFill_MultiLine&)` (`cxx:35-41`):
    /// `myIsDone` starts false and `Perform` runs right away.
    pub fn from_multi_line(ml: &MultiLine) -> Result<Self, Error> {
        let mut this = Self::new();
        this.perform(ml)?;
        Ok(this)
    }

    /// `Perform(const BRepFill_MultiLine&)` (`cxx:45-223`).
    pub fn perform(&mut self, ml: &MultiLine) -> Result<(), Error> {
        self.ml = ml.clone();

        // `cxx:49-68`: chord-length table of `2 * NbPoints` samples.
        let u1 = self.ml.first_parameter();
        let u2 = self.ml.last_parameter();
        let d_u = (u2 - u1) / (2.0 * NB_POINTS as f64 - 1.0);
        let mut length = 0.0;
        let mut lp: Vec<(f64, f64)> = Vec::with_capacity(2 * NB_POINTS);
        let mut a_pnt1 = self.ml.value(u1)?;
        for i in 0..2 * NB_POINTS {
            let a_pnt2 = self.ml.value(u1 + i as f64 * d_u)?;
            length += a_pnt1.distance(&a_pnt2);
            lp.push((length, u1 + i as f64 * d_u));
            a_pnt1 = a_pnt2;
        }

        // `cxx:72-78`: the first column, at `U1`.
        let mut mls = AppDefMultiLine::with_size(NB_POINTS);
        let (p3d, pf1, pf2) = self.ml.value_3d_on_f1_on_f2(u1)?;
        mls.set_value(1, &constraint(p3d, pf1, pf2));

        // `cxx:98-121`: one column every `DCorde`, the parameter interpolated
        // in the `LP` segment that brackets that arc length.
        let d_corde = length / (NB_POINTS as f64 - 1.0);
        let mut corde = d_corde;
        let mut index = 1usize;
        for i in 2..NB_POINTS {
            while lp[index - 1].0 < corde {
                index += 1;
            }
            let alpha = (corde - lp[index - 2].0) / (lp[index - 1].0 - lp[index - 2].0);
            let u = lp[index - 2].1 + alpha * (lp[index - 1].1 - lp[index - 2].1);
            let (p3d, pf1, pf2) = self.ml.value_3d_on_f1_on_f2(u)?;
            mls.set_value(i as i32, &constraint(p3d, pf1, pf2));
            corde = i as f64 * d_corde;
        }

        // `cxx:130-142`: the last column, at `U2`.
        let (p3d, pf1, pf2) = self.ml.value_3d_on_f1_on_f2(u2)?;
        mls.set_value(NB_POINTS as i32, &constraint(p3d, pf1, pf2));

        // `AppDef_Compute Fit(MLS);` (`cxx:144`), then `Fit.NbMultiCurves()`
        // (`cxx:146`).
        let mut fit = ComputeLine::from_line(
            &mls,
            4,
            8,
            1.0e-3,
            1.0e-6,
            5,
            true,
            ApproxParametrizationType::ChordLength,
            false,
        );
        if fit.nb_multi_curves() == 0 {
            return self.polygon_fallback(&mls);
        }

        // `cxx:169-223`: the approximation succeeded.
        let mbsp = fit.spline_value();
        let base = mbsp.base();
        let poles = base.curve3d(1);
        let poles2d1 = base.curve2d(2);
        let poles2d2 = base.curve2d(3);
        let knots = mbsp.knots().to_vec();
        let mults = mbsp.multiplicities().to_vec();
        let degree = mbsp.degree() as usize;

        let curve = GeomBSplineCurve::from_poles_knots_mults(poles, knots.clone(), mults.clone(), degree)
            .map_err(|_| Error::SplineConstruction)?;
        // `new Geom2d_BSplineCurve(Poles2d1, Knots, Mults, Degree)`
        // (`cxx:219-220`): the same distinct knots / multiplicities expanded
        // to the flat sequence `Geom2d_BSplineCurve` stores.
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree as i32);
        let (xs1, ys1) = split_pnt2d(&poles2d1);
        let (xs2, ys2) = split_pnt2d(&poles2d2);
        let p1 = Geom2dBSplineCurve::new(xs1, ys1, flat.clone(), degree)
            .map_err(|_| Error::SplineConstruction)?;
        let p2 = Geom2dBSplineCurve::new(xs2, ys2, flat, degree)
            .map_err(|_| Error::SplineConstruction)?;

        // `cxx:222-223`.
        self.curve = Some(Arc::new(curve));
        self.pcurve1 = Some(Arc::new(p1));
        self.pcurve2 = Some(Arc::new(p2));
        self.is_done = true;
        Ok(())
    }

    /// `cxx:146-166`: `NbCurves == 0`, so OCCT builds one degree-1 B-spline
    /// through the columns, with the knots at the sampling parameters and
    /// multiplicity 2 at both ends.
    fn polygon_fallback(&mut self, mls: &AppDefMultiLine) -> Result<(), Error> {
        let mut knots = vec![0.0; NB_POINTS];
        let mut mults = vec![1i32; NB_POINTS];
        mults[0] = 2;
        mults[NB_POINTS - 1] = 2;
        let mut p: Vec<GpPnt> = Vec::with_capacity(NB_POINTS);
        let mut p1: Vec<GpPnt2d> = Vec::with_capacity(NB_POINTS);
        let mut p2: Vec<GpPnt2d> = Vec::with_capacity(NB_POINTS);

        let uf = self.ml.first_parameter();
        let ul = self.ml.last_parameter();
        let d_ulf = (ul - uf) / (NB_POINTS as f64 - 1.0);
        // `for (i = 1; i <= NbPoints - 1; i++)` (`cxx:156`): the last column is
        // read separately right after (`cxx:160-165`).
        for i in 1..NB_POINTS {
            let mpc = mls.value(i as i32);
            p.push(*mpc.point(1));
            p1.push(*mpc.point2d(2));
            p2.push(*mpc.point2d(3));
            knots[i - 1] = uf + (i - 1) as f64 * d_ulf;
        }
        let mpc = mls.value(NB_POINTS as i32);
        p.push(*mpc.point(1));
        p1.push(*mpc.point2d(2));
        p2.push(*mpc.point2d(3));
        knots[NB_POINTS - 1] = ul;

        let curve = GeomBSplineCurve::from_poles_knots_mults(p, knots.clone(), mults.clone(), 1)
            .map_err(|_| Error::SplineConstruction)?;
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, 1);
        let (xs1, ys1) = split_pnt2d(&p1);
        let (xs2, ys2) = split_pnt2d(&p2);
        let c1 = Geom2dBSplineCurve::new(xs1, ys1, flat.clone(), 1)
            .map_err(|_| Error::SplineConstruction)?;
        let c2 =
            Geom2dBSplineCurve::new(xs2, ys2, flat, 1).map_err(|_| Error::SplineConstruction)?;

        self.curve = Some(Arc::new(curve));
        self.pcurve1 = Some(Arc::new(c1));
        self.pcurve2 = Some(Arc::new(c2));
        self.is_done = true;
        Ok(())
    }

    /// `IsDone()` (`cxx:227-230`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Curve()` (`cxx:236-240`): `StdFail_NotDone` when `Perform` failed.
    pub fn curve(&self) -> Result<&Arc<dyn Curve>, Error> {
        if !self.is_done {
            return Err(Error::NotDone);
        }
        self.curve.as_ref().ok_or(Error::SplineConstruction)
    }

    /// `CurveOnF1()` (`cxx:246-250`).
    pub fn curve_on_f1(&self) -> Result<&Arc<dyn Curve2d>, Error> {
        if !self.is_done {
            return Err(Error::NotDone);
        }
        self.pcurve1.as_ref().ok_or(Error::SplineConstruction)
    }

    /// `CurveOnF2()` (`cxx:256-260`).
    pub fn curve_on_f2(&self) -> Result<&Arc<dyn Curve2d>, Error> {
        if !self.is_done {
            return Err(Error::NotDone);
        }
        self.pcurve2.as_ref().ok_or(Error::SplineConstruction)
    }
}

/// One `AppDef_MultiPointConstraint(1, 2)` column: the 3D point at 1 and the
/// two p-curve points at 2 / 3 (`cxx:82-86`).
fn constraint(p3d: GpPnt, pf1: GpPnt2d, pf2: GpPnt2d) -> MultiPointConstraint {
    let mut mp = MultiPointConstraint::with_size(1, 2);
    mp.set_point(1, &p3d);
    mp.set_point2d(2, &pf1);
    mp.set_point2d(3, &pf2);
    mp
}

/// The `(x, y)` pole arrays of a `Geom2d_BSplineCurve` constructor.
fn split_pnt2d(poles: &[GpPnt2d]) -> (Vec<f64>, Vec<f64>) {
    (
        poles.iter().map(|p| p.x()).collect(),
        poles.iter().map(|p| p.y()).collect(),
    )
}
