//! `BRepFill_TrimEdgeTool` — the trimming of an offset edge pair by a 2D
//! bisector.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_TrimEdgeTool.cxx:17-781`
//! and `.hxx:30-80` (OCCT 8.0.0).
//!
//! The tool holds the two supports of the offset (`S1` / `S2`, each a
//! `Geom2d_Geometry`: a point or a curve), the offset value and the 2D bisector
//! between them. `BRepFill_OffsetWire` uses it in three ways:
//!
//! * [`TrimEdgeTool::intersect_with`] (`cxx:272-628`): intersect each of the two
//!   edges' p-curves with the bisector and answer the
//!   `(bisector, edge1, edge2)` parameter triples.
//! * [`TrimEdgeTool::add_or_confuse`] (`cxx:636-734`): when the start / end of
//!   the arc was not found among the intersections, project the bisector
//!   extremity onto both parallels and prepend / append it.
//! * [`TrimEdgeTool::is_inside`] (`cxx:737-781`): is a point closer to `S1`
//!   than the offset?
//!
//! Ported layout:
//! * this file — the type, [`Error`], [`Geometry2d`], the constructors and
//!   `IsInside`,
//! * [`eval_parameters`] — the static `Bubble` / `EvalParameters` /
//!   `EvalParametersBis` of `cxx:104-269`,
//! * [`intersect_with`] — `IntersectWith` (`cxx:272-628`),
//! * [`add_or_confuse`] — `AddOrConfuse` (`cxx:636-734`).
//!
//! Deviations (all of them handle / exception correspondence, none of them a
//! missing branch):
//! * `Bisector_Bisec` and `Bisector_BisecAna` have no port (the bisector is an
//!   `Arc<dyn Curve2d>` everywhere else in this crate), so the constructor takes
//!   the curve directly instead of `Bisector_Bisec`. `SimpleExpression`
//!   (`cxx:45-60`) only re-wraps a `Bisector_BisecAna` held by a
//!   `Geom2d_TrimmedCurve`; `Geom2dAdaptor_Curve::load` already replaces a
//!   trimmed curve by its basis while keeping the trimmed range
//!   (`Geom2dAdaptor_Curve.cxx:285-288`), so the adaptor the port builds from
//!   the curve equals the one `SimpleExpression` produces.
//! * `handle<Geom2d_Geometry>` is [`Geometry2d`]; a null `myC1` / `myC2` is the
//!   `Curve`-less state of a point support, which is what `isPoint1` /
//!   `isPoint2` guard in OCCT.
//! * `TopoDS_Vertex End1 / End2` may be a null handle at the `BRepFill_OffsetWire`
//!   call site (`BRepFill_OffsetWire.cxx:909-920`); the port takes
//!   `Option<&Vertex>` and `V.IsSame(null)` answers `false`, as in OCCT.
//! * `myBisec` (the `Bisector_Bisec` member, `hxx:75`) is not stored: only
//!   `myBis` (`hxx:76`) is read by the four methods.
//! * A null `BRep_Tool::CurveOnSurface` handle (`cxx:289-298`, `cxx:650-652`)
//!   raises `Standard_NullObject` in OCCT and answers [`Error::NullCurve`] here;
//!   `StdFail_NotDone` (`cxx:148`, `cxx:219`) answers
//!   [`Error::IntersectionNotDone`]. `TopoDS::Edge(InitShape1)` on a non-edge
//!   (`cxx:547`) raises `Standard_TypeMismatch` and answers
//!   [`Error::TypeMismatch`].
//! * The default constructor (`cxx:64`) is `= default`, so OCCT leaves
//!   `isPoint1` / `isPoint2` indeterminate; the port initialises them, since the
//!   only caller (`BRepFill_OffsetWire`) always uses the four-argument form.
//! * The `OCCT_DEBUG` printouts (`cxx:378-402`, `:500-511`, `:684-719`) are not
//!   ported.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

pub mod add_or_confuse;
pub mod eval_parameters;
pub mod intersect_with;

pub use eval_parameters::{bubble, eval_parameters, eval_parameters_bis};

/// `Geom2d_Geometry` as `BRepFill_TrimEdgeTool` reads it (`cxx:77-95`): a
/// `Geom2d_CartesianPoint` (`isPoint*`) or a `Geom2d_Curve` (`myC1` / `myC2`).
#[derive(Clone)]
pub enum Geometry2d {
    /// `Geom2d_CartesianPoint` — `myP1` / `myP2` (`hxx:70-71`).
    Point(GpPnt2d),
    /// `Geom2d_Curve` — `myC1` / `myC2` (`hxx:72-73`).
    Curve(Arc<dyn Curve2d>),
}

/// Failure modes of [`TrimEdgeTool`].
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// A null `BRep_Tool::CurveOnSurface` handle (`cxx:289-298`, `cxx:650-652`):
    /// OCCT dereferences it and raises `Standard_NullObject`.
    NullCurve,
    /// `StdFail_NotDone("BRepFill_TrimSurfaceTool::IntersectWith")`
    /// (`cxx:148`, `cxx:219`).
    IntersectionNotDone,
    /// `TopoDS::Edge(InitShape1)` on a shape that is not an edge (`cxx:547`):
    /// `Standard_TypeMismatch`.
    TypeMismatch,
}

/// `SimpleExpression` (`cxx:45-60`): a trimmed curve whose basis is a
/// `Bisector_BisecAna` is rebuilt as the trim of `Bisector_BisecAna::Geom2dCurve()`,
/// so that the adaptor type test sees the underlying analytic curve (line /
/// circle) instead of the `Bisector_BisecAna` wrapper. Any other curve is
/// answered unchanged.
fn simple_expression(bis: Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    let Some(basis) = bis.trimmed_basis() else {
        return bis;
    };
    let Some(ana) = basis.bisector_bisec_ana() else {
        return bis;
    };
    let Some(inner) = ana.geom2d_curve() else {
        return bis;
    };
    let (first, last) = (bis.first_parameter(), bis.last_parameter());
    Arc::new(Geom2dTrimmedCurve::new(Arc::clone(inner), first, last))
}

/// `BRepFill_TrimEdgeTool` (`hxx:37-77`).
#[derive(Clone)]
pub struct TrimEdgeTool {
    /// `isPoint1` (`hxx:68`).
    is_point1: bool,
    /// `isPoint2` (`hxx:69`).
    is_point2: bool,
    /// `myP1` (`hxx:70`).
    p1: GpPnt2d,
    /// `myP2` (`hxx:71`).
    p2: GpPnt2d,
    /// `myC1` (`hxx:72`), `None` on a point support.
    c1: Option<Arc<dyn Curve2d>>,
    /// `myC2` (`hxx:73`), `None` on a point support. OCCT never reads `myC2`
    /// either — only `myC1` is used, by `IsInside` (`cxx:752`) — but the member
    /// is kept so the two supports stay symmetric.
    #[allow(dead_code)]
    c2: Option<Arc<dyn Curve2d>>,
    /// `myOffset` (`hxx:74`).
    offset: f64,
    /// `myBis` (`hxx:76`), the `Geom2dAdaptor_Curve` over `SimpleExpression`
    /// of `myBisec`.
    bis: Geom2dAdaptorCurve,
}

impl Default for TrimEdgeTool {
    /// `BRepFill_TrimEdgeTool()` (`cxx:64`): OCCT's `= default` leaves the
    /// members uninitialised; the port answers the same shape with zero state.
    fn default() -> Self {
        Self {
            is_point1: false,
            is_point2: false,
            p1: GpPnt2d::zero(),
            p2: GpPnt2d::zero(),
            c1: None,
            c2: None,
            offset: 0.0,
            bis: Geom2dAdaptorCurve::new(),
        }
    }
}

impl TrimEdgeTool {
    /// `BRepFill_TrimEdgeTool(Bisec, S1, S2, Offset)` (`cxx:68-100`).
    ///
    /// `bisec` is the `Bisector_Bisec` value (see the module note); its
    /// `SimpleExpression` (`cxx:45-60`) collapses into
    /// [`Geom2dAdaptorCurve::from_curve`], which keeps the trimmed range while
    /// unwrapping the basis (`Geom2dAdaptor_Curve.cxx:285-288`).
    pub fn new(bisec: Arc<dyn Curve2d>, s1: &Geometry2d, s2: &Geometry2d, offset: f64) -> Self {
        let is_point1 = matches!(s1, Geometry2d::Point(_));
        let is_point2 = matches!(s2, Geometry2d::Point(_));

        // `cxx:80-95`: geometries of the supports.
        let (p1, c1) = match s1 {
            Geometry2d::Point(p) => (*p, None),
            Geometry2d::Curve(c) => (GpPnt2d::zero(), Some(Arc::clone(c))),
        };
        let (p2, c2) = match s2 {
            Geometry2d::Point(p) => (*p, None),
            Geometry2d::Curve(c) => (GpPnt2d::zero(), Some(Arc::clone(c))),
        };

        // `cxx:97-99`: `SimpleExpression(myBisec, Bis); myBis =
        // Geom2dAdaptor_Curve(Bis);`.
        let bis = Geom2dAdaptorCurve::from_curve(simple_expression(bisec));

        Self {
            is_point1,
            is_point2,
            p1,
            p2,
            c1,
            c2,
            offset,
            bis,
        }
    }

    /// `myOffset` (`hxx:74`).
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// `myBis` (`hxx:76`).
    pub fn bisector(&self) -> &Geom2dAdaptorCurve {
        &self.bis
    }

    /// `IsInside(P)` (`cxx:737-781`): `|P - S1| < |Offset| - Precision::Confusion()`,
    /// where the distance to a curve support is the projection distance,
    /// further reduced by the distance to the curve's two endpoints
    /// (`cxx:764-775`).
    pub fn is_inside(&self, p: &GpPnt2d) -> bool {
        let mut dist = f64::MAX;
        if self.is_point1 {
            // `cxx:743-746`.
            dist = p.distance(&self.p1);
        } else if self.is_point2 {
            // `cxx:747-750`.
            dist = p.distance(&self.p2);
        } else if let Some(c1) = self.c1.as_ref() {
            // `cxx:752-776`.
            let projector = ProjectPointOnCurve2d::with_curve(p, c1);
            if projector.nb_points() > 0 {
                dist = projector.lower_distance().unwrap_or(f64::MAX);
            }
            // `cxx:766-775`: the endpoint check runs even when the projection
            // succeeded.
            let pf = c1.d0(c1.first_parameter());
            let pl = c1.d0(c1.last_parameter());
            let dist_min = p.distance(&pf).min(p.distance(&pl));
            if dist > dist_min {
                dist = dist_min;
            }
        }
        // `cxx:781`: note the `- Precision::Confusion()`, not `+`.
        dist < self.offset.abs() - PCONFUSION
    }
}
