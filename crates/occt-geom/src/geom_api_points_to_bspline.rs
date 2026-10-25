//! `GeomAPI_PointsToBSpline` — approximate a point sequence by a B-spline.
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomAPI/GeomAPI_PointsToBSpline.cxx`
//! (`:41-357`), declaration `GeomAPI_PointsToBSpline.hxx`.
//!
//! **Ported (faithful)**:
//! - default constructor (`cxx:43-47`);
//! - `Init(Points, ParType, DegMin, DegMax, Continuity, Tol3D)` (`cxx:117-166`)
//!   and the convenience overloads that forward to it (`cxx:52-112`);
//! - `Init(Points, Params, ...)` (`cxx:170-238`);
//! - `Curve` / `IsDone` (`cxx:335-357`).
//!
//! The fit itself is `AppDef_BSplineCompute`, ported as
//! [`BsplComputeLine`](crate::approx::bspl_compute_line::BsplComputeLine); the
//! point sequence becomes an `AppDef_MultiLine` through
//! [`MultiLine::from_3d`](crate::app_def::multi_line::MultiLine::from_3d)
//! (`AppDef_MultiLine(TabP3d)`).
//!
//! **UNPORTED**: `Init(Points, W1, W2, W3, DegMax, Continuity, Tol3D)`
//! (`cxx:243-329`), which runs `AppDef_Variational`. That class has no port in
//! this crate, so the overload is left out rather than emulated by the
//! least-squares computer.
//!
//! Differences in shape only: `myIsDone` / `myCurve` become
//! [`PointsToBSpline::is_done`] / [`PointsToBSpline::curve`]; the
//! `StdFail_NotDone` throw of `Curve()` (`cxx:339-343`) is an `Err`.

use occt_core::approx::ApproxParametrizationType;
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::REAL_SMALL;
use occt_math::MathVector;

use crate::app_def::multi_line::MultiLine;
use crate::approx::bspl_compute_line::BsplComputeLine;
use crate::bspline_curve::GeomBSplineCurve;

/// `GeomAPI_PointsToBSpline` (`hxx:53-...`): `myCurve` / `myIsDone`
/// (`hxx:186-188`).
#[derive(Clone, Default)]
pub struct PointsToBSpline {
    curve: Option<GeomBSplineCurve>,
    is_done: bool,
}

impl PointsToBSpline {
    /// `GeomAPI_PointsToBSpline()` (`cxx:43-47`).
    pub fn new() -> Self {
        Self {
            curve: None,
            is_done: false,
        }
    }

    /// `GeomAPI_PointsToBSpline(Points, DegMin, DegMax, Continuity, Tol3D)`
    /// (`cxx:52-59`): forwards to the `Approx_ChordLength` overload.
    pub fn with_degrees(
        points: &[GpPnt],
        deg_min: i32,
        deg_max: i32,
        continuity: Shape,
        tol3d: f64,
    ) -> Self {
        Self::with_params_type(
            points,
            ApproxParametrizationType::ChordLength,
            deg_min,
            deg_max,
            continuity,
            tol3d,
        )
    }

    /// `GeomAPI_PointsToBSpline(Points, ParType, DegMin, DegMax, Continuity,
    /// Tol3D)` (`cxx:64-72`).
    pub fn with_params_type(
        points: &[GpPnt],
        par_type: ApproxParametrizationType,
        deg_min: i32,
        deg_max: i32,
        continuity: Shape,
        tol3d: f64,
    ) -> Self {
        let mut this = Self::new();
        this.init(points, par_type, deg_min, deg_max, continuity, tol3d);
        this
    }

    /// `Init(Points, DegMin, DegMax, Continuity, Tol3D)` (`cxx:104-112`): the
    /// parametrisation is `Approx_ChordLength`.
    pub fn init(
        &mut self,
        points: &[GpPnt],
        par_type: ApproxParametrizationType,
        deg_min: i32,
        deg_max: i32,
        continuity: Shape,
        tol3d: f64,
    ) {
        self.is_done = false;
        // `double Tol2D = 0.;` (`cxx:122`), `int nbit = 2` (`cxx:124`).
        let tol2d = 0.0;
        let nbit = 2;
        // `if (Tol3D <= 1.e-3) UseSquares = true;` (`cxx:126-131`).
        let use_squares = tol3d <= 1.0e-3;

        let mut computer = BsplComputeLine::with_settings(
            deg_min, deg_max, tol3d, tol2d, nbit, true, par_type, use_squares,
        );
        computer.set_continuity(continuity_index(continuity));

        let line = MultiLine::from_3d(points);
        computer.perform(&line);

        let value = computer.value();
        let poles = value.base().curve3d(1);
        self.curve = GeomBSplineCurve::from_poles_knots_mults(
            poles,
            value.knots().to_vec(),
            value.multiplicities().to_vec(),
            value.degree() as usize,
        )
        .ok();
        // `myIsDone = true;` (`cxx:165`) is unconditional in OCCT.
        self.is_done = true;
    }

    /// `GeomAPI_PointsToBSpline(Points, Params, DegMin, DegMax, Continuity,
    /// Tol3D)` (`cxx:76-85`).
    pub fn with_parameters(
        points: &[GpPnt],
        params: &[f64],
        deg_min: i32,
        deg_max: i32,
        continuity: Shape,
        tol3d: f64,
    ) -> Result<Self, String> {
        let mut this = Self::new();
        this.init_with_parameters(points, params, deg_min, deg_max, continuity, tol3d)?;
        Ok(this)
    }

    /// `Init(Points, Params, ...)` (`cxx:170-238`).
    ///
    /// `Err` is the `Standard_OutOfRange` raised when the two arrays differ in
    /// length (`cxx:177-181`). A degenerate parameter span leaves the object
    /// not done (`cxx:190-193`).
    pub fn init_with_parameters(
        &mut self,
        points: &[GpPnt],
        params: &[f64],
        deg_min: i32,
        deg_max: i32,
        continuity: Shape,
        tol3d: f64,
    ) -> Result<(), String> {
        self.is_done = false;
        if params.len() != points.len() {
            return Err("GeomAPI_PointsToBSpline::Init() - invalid input".into());
        }

        // `double Tol2D = 0.; int Nbp = Params.Length();` (`cxx:182-184`).
        let tol2d = 0.0;
        let nbp = params.len();
        let mut the_params = MathVector::new(1, nbp);
        the_params.set_value(1, 0.0);
        the_params.set_value(nbp, 1.0);

        // `double Uf = Params(Params.Lower());`
        // `double Ul = Params(Params.Upper()) - Uf;` (`cxx:188-189`).
        let uf = match params.first() {
            Some(v) => *v,
            None => return Ok(()),
        };
        let ul = match params.last() {
            Some(v) => *v - uf,
            None => return Ok(()),
        };
        // `hasMeaningfulSpan(Ul)` (`cxx:34-39`, `cxx:190-193`): the bound is
        // `gp::Resolution()` = `RealSmall()` (`gp.hxx:60`), not
        // `Precision::Angular()`.
        if !(ul * ul > REAL_SMALL * REAL_SMALL) {
            return Ok(());
        }

        // `for (int i = 2; i < Nbp; i++)` (`cxx:195-198`): the OCCT loop stops
        // before `Nbp`, whose entry keeps the 1. written above.
        for i in 2..nbp {
            the_params.set_value(i, (params[i - 1] - uf) / ul);
        }

        // `AppDef_BSplineCompute TheComputer(DegMin, DegMax, Tol3D, Tol2D, 0,
        // true, Approx_IsoParametric, true);` (`cxx:200-201`).
        let mut computer = BsplComputeLine::with_settings(
            deg_min,
            deg_max,
            tol3d,
            tol2d,
            0,
            true,
            ApproxParametrizationType::IsoParametric,
            true,
        );
        computer.set_parameters(&the_params);
        computer.set_continuity(continuity_index(continuity));

        let line = MultiLine::from_3d(points);
        computer.perform(&line);

        let value = computer.value();
        let poles = value.base().curve3d(1);
        // `BSplCLib::Reparametrize(Params(Lower), Params(Upper), Knots);`
        // (`cxx:235`).
        let mut knots = value.knots().to_vec();
        reparametrize(uf, uf + ul, &mut knots);
        self.curve = GeomBSplineCurve::from_poles_knots_mults(
            poles,
            knots,
            value.multiplicities().to_vec(),
            value.degree() as usize,
        )
        .ok();
        self.is_done = true;
        Ok(())
    }

    /// `Curve()` (`cxx:335-343`): the `StdFail_NotDone` throw is the `Err`.
    pub fn curve(&self) -> Result<&GeomBSplineCurve, String> {
        if !self.is_done {
            return Err("GeomAPI_PointsToBSpline::Curve ".into());
        }
        match self.curve.as_ref() {
            Some(c) => Ok(c),
            None => Err("GeomAPI_PointsToBSpline::Curve ".into()),
        }
    }

    /// `IsDone()` (`cxx:353-356`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }
}

/// `switch (Continuity)` (`cxx:133-153`, `cxx:203-224`): `SetContinuity` order
/// of the computer. `GeomAbs_G1` shares the `GeomAbs_C1` arm and `GeomAbs_G2`
/// the `GeomAbs_C2` one; anything above `C2` becomes 3.
fn continuity_index(continuity: Shape) -> i32 {
    if continuity == Shape::C0 {
        0
    } else if continuity == Shape::C1 || continuity == Shape::G1 {
        1
    } else if continuity == Shape::C2 || continuity == Shape::G2 {
        2
    } else {
        3
    }
}
