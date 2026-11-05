//! `Bisector` package (`TKTopAlgo/Bisector`): the analytic bisecting curves
//! used to build offset wires.
//!
//! | OCCT | lines |
//! |---|---|
//! | `Bisector_Curve` | `Bisector_Curve.hxx:32-51`, `Bisector_Curve.cxx:1-21` |
//! | `Bisector_BisecAna` | `Bisector_BisecAna.cxx:17-1825` |
//! | `Bisector_Bisec` | `Bisector_Bisec.cxx:49-681` |
//! | `Bisector_BisecCC` | `Bisector_BisecCC.cxx:60-1960` |
//! | `Bisector_BisecPC` | `Bisector_BisecPC.cxx` |
//! | `Bisector_FunctionH` | `Bisector_FunctionH.cxx:19-85` |
//! | `Bisector_FunctionInter` | `Bisector_FunctionInter.cxx:19-106` |
//! | `Bisector_Inter` | `Bisector_Inter.cxx:19-463` |
//!
//! `Bisector_BisecAna` wraps a `Geom2d_TrimmedCurve` (`thebisector`) whose basis
//! is one of the elementary conics produced by the `GccAna_*Bisec` solvers; it
//! does not own any geometry of its own. The four `Perform` overloads build that
//! trimmed curve from the `GccAna` solutions, and `Distance` picks the sector of
//! the candidate bisector that matches the input tangents.
//!
//! Not ported: none of the `Bisector` package remains.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::CONFUSION;

use crate::curve::Curve2d;
use crate::trimmed::Geom2dTrimmedCurve;

mod bisec;
mod bisec_cc;
mod bisec_pc;
mod degenerate;
mod distance;
mod function_h;
mod function_inter;
mod inter;
mod null_curve;
mod perform_curve_curve;
mod perform_point;
mod point_on_bis;
mod poly_bis;

pub use bisec::{BisectorBisec, BisectorBisecError};
pub use bisec_cc::BisectorBisecCC;
pub use bisec_pc::BisectorBisecPC;
pub use degenerate::degenerate;
pub use function_h::BisectorFunctionH;
pub use function_inter::BisectorFunctionInter;
pub use inter::BisectorInter;
pub use null_curve::{null_curve, NullCurve};
pub use point_on_bis::BisectorPointOnBis;
pub use poly_bis::BisectorPolyBis;

/// `Bisector::IsConvex` (`Bisector.cxx:26-34`): true when the curvature of `Cu`
/// at its mid parameter has the sign of `Sign`.
pub fn is_convex(cu: &dyn Curve2d, sign: f64) -> bool {
    let u1 = (cu.last_parameter() + cu.first_parameter()) / 2.0;
    let (_, v1, v2) = cu.d2(u1);
    const TOL: f64 = 1.0e-5;
    sign * v1.crossed(&v2) < TOL
}

/// `Geom2d_Curve::IsClosed()` for the curve kinds this package can receive.
///
/// Circles and ellipses are closed (`Geom2d_Circle.cxx:182-185`,
/// `Geom2d_Ellipse.cxx:191-194`), a trimmed curve is closed when its range
/// covers a whole period (`Geom2d_TrimmedCurve.cxx`), and the offset, B-spline
/// and Bezier forms compare their end points against
/// `Precision::Computational()` (`Geom2d_OffsetCurve.cxx`,
/// `Geom2d_BSplineCurve_1.cxx`).
///
/// UNPORTED: `Geom2d_BezierCurve::IsClosed` returns the stored `myClosed` flag,
/// which the port's `Geom2dBezierCurve` does not keep; it therefore takes the
/// end-point comparison.
pub(crate) fn curve_is_closed(c: &dyn Curve2d) -> bool {
    if c.is_geom2d_circle() || c.gp_elips2d().is_some() {
        return true;
    }
    if let Some(basis) = c.trimmed_basis() {
        if basis.is_periodic() {
            let period = basis.period();
            let length = c.last_parameter() - c.first_parameter();
            if length > occt_core::precision::PCONFUSION {
                let rem = length - (length / period).round() * period;
                if rem.abs() <= occt_core::precision::PCONFUSION {
                    return true;
                }
            }
        }
    }
    if c.offset_basis().is_some() || c.is_bspline2d() || c.is_bezier2d() {
        return c
            .d0(c.first_parameter())
            .square_distance(&c.d0(c.last_parameter()))
            <= occt_core::precision::COMPUTATIONAL;
    }
    false
}

/// `class Bisector_Curve` (`Bisector_Curve.hxx:30-51`): the extra queries the
/// bisector family adds on top of `Geom2d_Curve`.
pub trait BisectorCurve: Curve2d {
    /// `Bisector_Curve::IsExtendAtStart`.
    fn is_extend_at_start(&self) -> bool;
    /// `Bisector_Curve::IsExtendAtEnd`.
    fn is_extend_at_end(&self) -> bool;
    /// `Bisector_Curve::NbIntervals`.
    fn nb_intervals(&self) -> i32;
    /// `Bisector_Curve::IntervalFirst`.
    fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError>;
    /// `Bisector_Curve::IntervalLast`.
    fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError>;
}

/// Exceptions raised by `Bisector_BisecAna`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BisectorBisecAnaError {
    /// `Standard_OutOfRange` - `IntervalFirst` / `IntervalLast` with `I != 1`.
    OutOfRange,
    /// `StdFail_NotDone` - the `switch (cas)` default arm of
    /// `Bisector_BisecAna.cxx:957-961`.
    NotDone,
    /// OCCT dereferences a null `thebisector` handle here; the Rust port reports
    /// it instead.
    NullBisector,
}

impl std::fmt::Display for BisectorBisecAnaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::OutOfRange => "Standard_OutOfRange",
            Self::NotDone => "StdFail_NotDone",
            Self::NullBisector => "null thebisector",
        };
        write!(f, "{s}")
    }
}

impl std::error::Error for BisectorBisecAnaError {}

impl From<crate::gcc_ana_bisec::GccAnaBisecError> for BisectorBisecAnaError {
    /// A `GccAna_*Bisec` failure surfaces as `StdFail_NotDone` on the bisector
    /// side: the callers of `ThisSolution` in `Bisector_BisecAna.cxx` are all
    /// behind an `IsDone()` guard, so only the `NbrSol` guard can raise, and it
    /// raises `Standard_OutOfRange`.
    fn from(e: crate::gcc_ana_bisec::GccAnaBisecError) -> Self {
        use crate::gcc_ana_bisec::GccAnaBisecError;
        match e {
            GccAnaBisecError::OutOfRange => Self::OutOfRange,
            GccAnaBisecError::NotDone | GccAnaBisecError::NoSolution => Self::NotDone,
        }
    }
}

/// `class Bisector_BisecAna`.
#[derive(Clone, Default)]
pub struct BisectorBisecAna {
    /// `occ::handle<Geom2d_TrimmedCurve> thebisector` (`BisecAna.hxx:180`). The
    /// handle is null after the default constructor and after a `Perform` that
    /// found no solution, exactly as in OCCT.
    thebisector: Option<Geom2dTrimmedCurve>,
}

impl BisectorBisecAna {
    /// `Bisector_BisecAna::Bisector_BisecAna` (`cxx:53`): leaves the handle null.
    pub fn new() -> Self {
        Self { thebisector: None }
    }

    /// `Bisector_BisecAna::Init` (`cxx:1725-1731`).
    pub fn init(&mut self, bis: Geom2dTrimmedCurve) {
        self.thebisector = Some(bis);
    }

    /// The wrapped curve, `None` while the internal handle is null.
    pub fn thebisector(&self) -> Option<&Geom2dTrimmedCurve> {
        self.thebisector.as_ref()
    }

    /// `Bisector_BisecAna::IsExtendAtStart` (`cxx:1305-1311`): always false.
    pub fn is_extend_at_start(&self) -> bool {
        false
    }

    /// `Bisector_BisecAna::IsExtendAtEnd` (`cxx:1313-1318`): always false.
    pub fn is_extend_at_end(&self) -> bool {
        false
    }

    /// `Bisector_BisecAna::SetTrim(const occ::handle<Geom2d_Curve>&)`
    /// (`cxx:1328-1516`): the whole body is commented out in OCCT, so this
    /// overload is a no-op there and here.
    pub fn set_trim_on_curve(&mut self, _cu: &dyn Curve2d) {}

    /// `Bisector_BisecAna::SetTrim(uf, ul)` (`cxx:1518-1522`).
    pub fn set_trim(&mut self, uf: f64, ul: f64) {
        if let Some(b) = self.thebisector.as_mut() {
            b.set_trim_range(uf, ul);
        }
    }

    /// `Bisector_BisecAna::Reverse` (`cxx:1526-1529`).
    pub fn reverse_curve(&mut self) {
        if let Some(b) = self.thebisector.as_mut() {
            b.reverse();
        }
    }

    /// `Bisector_BisecAna::ReversedParameter` (`cxx:1535-1538`).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        match self.thebisector.as_ref() {
            Some(b) => b.reversed_parameter(u),
            None => u,
        }
    }

    /// `Bisector_BisecAna::IsCN` (`cxx:1542-1546`): `thebisector->IsCN(N)`,
    /// which for `Geom2d_TrimmedCurve` is `BasisCurve()->IsCN(N)` - i.e.
    /// `Continuity() >= N` on the `GeomAbs_Shape` scale.
    pub fn is_cn(&self, n: i32) -> bool {
        match self.thebisector.as_ref() {
            Some(b) => b.continuity() as i32 >= n,
            None => false,
        }
    }

    /// `Bisector_BisecAna::Copy` (`cxx:1550-1556`): copies the wrapped trimmed
    /// curve into a fresh `Bisector_BisecAna`.
    pub fn copy_bisector(&self) -> Self {
        match self.thebisector.as_ref() {
            Some(b) => Self {
                thebisector: Some(b.clone()),
            },
            None => Self::new(),
        }
    }

    /// `Bisector_BisecAna::Transform` (`cxx:1560-1563`).
    pub fn transform(&mut self, t: &GpTrsf2d) {
        if let Some(b) = self.thebisector.as_mut() {
            b.transform(t);
        }
    }

    /// `Bisector_BisecAna::FirstParameter` (`cxx:1566-1569`): the *trimmed*
    /// first parameter (the `BasisCurve()->FirstParameter()` form is commented
    /// out at `cxx:1564-1567`).
    pub fn first_parameter(&self) -> f64 {
        match self.thebisector.as_ref() {
            Some(b) => b.first_parameter(),
            None => 0.0,
        }
    }

    /// `Bisector_BisecAna::LastParameter` (`cxx:1574-1577`).
    pub fn last_parameter(&self) -> f64 {
        match self.thebisector.as_ref() {
            Some(b) => b.last_parameter(),
            None => 0.0,
        }
    }

    /// `Bisector_BisecAna::IsClosed` (`cxx:1581-1584`): forwarded to the
    /// **basis** curve.
    pub fn is_closed(&self) -> bool {
        match self.thebisector.as_ref() {
            Some(b) => {
                let bp = b.basis_curve();
                bp.first_parameter() == bp.last_parameter()
            }
            None => false,
        }
    }

    /// `Bisector_BisecAna::IsPeriodic` (`cxx:1588-1591`): forwarded to the
    /// **basis** curve.
    pub fn is_periodic(&self) -> bool {
        match self.thebisector.as_ref() {
            Some(b) => b.basis_curve().is_periodic(),
            None => false,
        }
    }

    /// `Bisector_BisecAna::Continuity` (`cxx:1595-1598`): the *trimmed*
    /// continuity.
    pub fn continuity(&self) -> u8 {
        match self.thebisector.as_ref() {
            Some(b) => b.continuity(),
            None => 0,
        }
    }

    /// The basis curve behind the trim, `None` while the handle is null.
    pub fn geom2d_curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.thebisector.as_ref().map(|b| b.basis_curve())
    }

    /// `Bisector_BisecAna::ParameterOfStartPoint` (`cxx:1644-1647`).
    pub fn parameter_of_start_point(&self) -> f64 {
        self.first_parameter()
    }

    /// `Bisector_BisecAna::ParameterOfEndPoint` (`cxx:1652-1655`).
    pub fn parameter_of_end_point(&self) -> f64 {
        self.last_parameter()
    }

    /// `Bisector_BisecAna::NbIntervals` (`cxx:1697-1700`): always one interval.
    pub fn nb_intervals(&self) -> i32 {
        1
    }

    /// `Bisector_BisecAna::IntervalFirst` (`cxx:1705-1712`).
    pub fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        if index != 1 {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.first_parameter())
    }

    /// `Bisector_BisecAna::IntervalLast` (`cxx:1717-1723`).
    pub fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        if index != 1 {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.last_parameter())
    }

    /// True while `thebisector` is the null handle.
    pub fn is_null(&self) -> bool {
        self.thebisector.is_none()
    }
}

impl Curve2d for BisectorBisecAna {
    /// `Bisector_BisecAna::EvalD0` (`cxx:1602-1605`): evaluated on the basis
    /// curve, not on the trimmed curve.
    fn d0(&self, u: f64) -> GpPnt2d {
        match self.geom2d_curve() {
            Some(b) => b.d0(u),
            None => GpPnt2d::zero(),
        }
    }

    /// `Bisector_BisecAna::EvalD1` (`cxx:1609-1612`).
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        match self.geom2d_curve() {
            Some(b) => b.d1(u),
            None => (GpPnt2d::zero(), GpVec2d::zero()),
        }
    }

    /// `Bisector_BisecAna::EvalD2` (`cxx:1616-1619`).
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        match self.geom2d_curve() {
            Some(b) => b.d2(u),
            None => (GpPnt2d::zero(), GpVec2d::zero(), GpVec2d::zero()),
        }
    }

    /// `Bisector_BisecAna::EvalD3` (`cxx:1623-1626`).
    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        match self.geom2d_curve() {
            Some(b) => b.d3(u),
            None => (
                GpPnt2d::zero(),
                GpVec2d::zero(),
                GpVec2d::zero(),
                GpVec2d::zero(),
            ),
        }
    }

    /// `Bisector_BisecAna::EvalDN` (`cxx:1630-1633`).
    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        match self.geom2d_curve() {
            Some(b) => b.eval_dn(u, n),
            None => GpVec2d::zero(),
        }
    }

    fn first_parameter(&self) -> f64 {
        BisectorBisecAna::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecAna::last_parameter(self)
    }

    fn is_periodic(&self) -> bool {
        BisectorBisecAna::is_periodic(self)
    }

    fn continuity(&self) -> u8 {
        BisectorBisecAna::continuity(self)
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        BisectorBisecAna::transform(self, t);
    }

    fn reverse(&mut self) {
        BisectorBisecAna::reverse_curve(self);
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecAna::reversed_parameter(self, u)
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }

    /// `STANDARD_TYPE(Bisector_BisecAna)` downcast (`Bisector_Inter.cxx:239`).
    fn bisector_bisec_ana(&self) -> Option<&BisectorBisecAna> {
        Some(self)
    }
}

impl BisectorCurve for BisectorBisecAna {
    fn is_extend_at_start(&self) -> bool {
        BisectorBisecAna::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecAna::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecAna::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecAna::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecAna::interval_last(self, index)
    }
}

/// `Precision::Confusion()` re-export so the `Perform` overloads and their
/// callers agree on the tolerance the OCCT code hard-codes as `PreConf`.
pub const PRE_CONF: f64 = CONFUSION;

/// `Bisector_Curve::Parameter(const gp_Pnt2d&)` (`Bisector_Curve.hxx`): the
/// parameter of `P` on a bisector, dispatched on the concrete family.
///
/// `MAT2d_Tool2d::TrimBisector` (`MAT2d_Tool2d.cxx:530`) and
/// `MAT2d_Tool2d::BisecFusion` (`cxx:1228`) reach this through the virtual
/// `Bisector_Curve` base of `Geom2d_TrimmedCurve::BasisCurve()`.
pub fn bisector_parameter(c: &dyn Curve2d, p: &GpPnt2d) -> f64 {
    if let Some(cc) = c.bisector_bisec_cc() {
        return cc.parameter(p);
    }
    if let Some(pc) = c.bisector_bisec_pc() {
        return pc.parameter(p);
    }
    if let Some(ana) = c.bisector_bisec_ana() {
        return ana.parameter(p);
    }
    panic!("Bisector_Curve::Parameter: not a bisector curve");
}
