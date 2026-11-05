//! `Geom2dAdaptor_Curve` queries that are not a plain delegation.
//!
//! Sources: `Geom2dAdaptor_Curve.cxx` (the `load` type table `:279-346`,
//! `GetType` `:96-120`). These are the queries an `Adaptor2d_Curve2d` facade is
//! expected to answer on the *stored* curve, so they live on their own rather
//! than inside one consumer (`Extrema_Curve2dTool`, `GeomBndLib_Curve2d` and
//! `BndLib_Add2dCurve` all read the same table).

use std::sync::Arc;

use occt_core::gp::{GpCirc2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::CurveType;

use crate::curve::Curve2d;

/// Unwrap `Geom2d_TrimmedCurve`s to their basis:
/// `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:285-288`) replaces a
/// trimmed curve with its basis and keeps `myFirst`/`myLast` at the trimmed
/// range, so the *type* resolves on the basis while the range stays trimmed.
pub fn basis2d(mut c: &dyn Curve2d) -> &dyn Curve2d {
    while let Some(b) = c.trimmed_basis() {
        c = b;
    }
    c
}

/// `Geom2dAdaptor_Curve::GetType` (`Geom2dAdaptor_Curve.cxx:96-120`), i.e. the
/// `load` type table (`:279-346`) resolved on the basis of a trimmed curve.
///
/// The type queries are the port's stand-ins for OCCT's
/// `DynamicType() == STANDARD_TYPE(Geom2d_*)` branches; each is answered only
/// by the matching concrete curve, and a `Geom2d_TrimmedCurve` forwards to its
/// basis because `load` recurses into it (`:285-288`).
pub fn type2d(c: &dyn Curve2d) -> CurveType {
    let base = basis2d(c);
    if base.is_line() {
        CurveType::Line // `:292-296`
    } else if base.gp_circ2d().is_some() {
        CurveType::Circle // `:289-292`
    } else if base.gp_elips2d().is_some() {
        CurveType::Ellipse // `:299-303`
    } else if base.gp_parab2d().is_some() {
        CurveType::Parabola // `:304-308`
    } else if base.gp_hypr2d().is_some() {
        CurveType::Hyperbola // `:309-313`
    } else if base.offset_basis().is_some() {
        CurveType::OffsetCurve // `:330-342`
    } else if base.bezier_nb_poles().is_some() {
        CurveType::BezierCurve // `:314-321`
    } else if base.bspline_nb_knots().is_some() || base.bspline_degree().is_some() {
        CurveType::BSplineCurve // `:322-329`
    } else {
        CurveType::OtherCurve // `:343-345`
    }
}

/// `Geom2dAdaptor_Curve` (`Geom2dAdaptor_Curve.hxx:88-260`), partial port: the
/// `load` state (`myCurve` / `myFirst` / `myLast` / `myTypeCurve`) plus the
/// queries the sweep tools ask of it.
///
/// Ported entries:
/// - `Geom2dAdaptor_Curve(C, UFirst, ULast)` (`cxx:250-258`) and `load`
///   (`cxx:274-346`): the range is stored as given and a `Geom2d_TrimmedCurve`
///   is replaced by its basis (`cxx:285-288`), recursively, while keeping the
///   caller's range.
/// - `Load(C)` (`cxx:250-258` pair): re-load a curve into an existing adaptor
///   without touching `myFirst` / `myLast`.
/// - `Curve` (`hxx:136`), `FirstParameter` / `LastParameter` (`hxx:138-140`),
///   `Value` (`hxx:172`), `DN` (`hxx:204`), `GetType` (`hxx:209`), `Circle`
///   (`hxx:213`), `Line` (`hxx:...`).
///
/// UNPORTED (`Geom2dAdaptor_Curve.cxx`): the `EvalRepresentation` / cache
/// variants (`BezierData` / `BSplineData` / `OffsetData`, `cxx:315-342`), the
/// `ShallowCopy` (`cxx:89-134`), `Resolution` (`cxx:1186-1219`), `Intervals` /
/// `NbIntervals` and `Trim` — none of them is reached by `BRepFill_MultiLine`.
///
/// A null `Handle(Geom2d_Curve)` is a legitimate state in OCCT (a null curve
/// still carries `myFirst` / `myLast`); [`Option`] models it here.
#[derive(Clone)]
pub struct Geom2dAdaptorCurve {
    /// `myCurve` (`cxx:280`): the handle after `load` unwrapped the trimmed
    /// curves.
    curve: Option<Arc<dyn Curve2d>>,
    /// `myFirst` (`hxx:...`).
    first: f64,
    /// `myLast`.
    last: f64,
}

impl Default for Geom2dAdaptorCurve {
    /// `Geom2dAdaptor_Curve()` (`cxx:238-244`): `myFirst = myLast = 0.0`,
    /// no curve, type `GeomAbs_OtherCurve`.
    fn default() -> Self {
        Self::new()
    }
}

impl Geom2dAdaptorCurve {
    /// `Geom2dAdaptor_Curve()` (`cxx:238-244`).
    pub fn new() -> Self {
        Self {
            curve: None,
            first: 0.0,
            last: 0.0,
        }
    }

    /// `Geom2dAdaptor_Curve(C)` (`cxx:245-248`): range from the curve itself.
    pub fn from_curve(curve: Arc<dyn Curve2d>) -> Self {
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        Self::with_range(curve, first, last)
    }

    /// `Geom2dAdaptor_Curve(C, UFirst, ULast)` (`cxx:250-258`): `myFirst` /
    /// `myLast` are the arguments (`cxx:254-255`) and `load` (`cxx:274-346`)
    /// keeps them while unwrapping a `Geom2d_TrimmedCurve` (`cxx:285-288`).
    pub fn with_range(curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Self {
        Self {
            curve: Some(Self::unwrap_trimmed(curve)),
            first,
            last,
        }
    }

    /// `load` (`cxx:285-288`) / `Load` (`hxx:118-130`): a trimmed curve is
    /// replaced by its basis; the recursion stops on the first non-trimmed
    /// curve.
    fn unwrap_trimmed(mut curve: Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
        loop {
            let basis = curve
                .trimmed_basis()
                .map(|b| Arc::from(b.clone_dyn()) as Arc<dyn Curve2d>);
            match basis {
                Some(b) => curve = b,
                None => return curve,
            }
        }
    }

    /// `Load(C)` (`Geom2dAdaptor_Curve.cxx:250-258` overload pair): the range
    /// (`myFirst` / `myLast`) is kept and the curve is re-loaded, so a trimmed
    /// curve is unwrapped to its basis again (`load`, `cxx:274-346`).
    pub fn load(&mut self, curve: Arc<dyn Curve2d>) {
        self.curve = Some(Self::unwrap_trimmed(curve));
    }

    /// `Curve()` (`hxx:136`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.curve.as_ref()
    }

    /// `FirstParameter()` (`hxx:138`).
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `LastParameter()` (`hxx:140`).
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Value(U)` (`hxx:172`): `D0` on the loaded curve.
    pub fn value(&self, u: f64) -> Option<GpPnt2d> {
        self.curve.as_ref().map(|c| c.d0(u))
    }

    /// `DN(U, N)` (`hxx:204`).
    pub fn dn(&self, u: f64, n: i32) -> Option<GpVec2d> {
        self.curve.as_ref().map(|c| c.eval_dn(u, n))
    }

    /// `GetType()` (`hxx:209`): the type stored by `load` (`cxx:288-346`),
    /// which is the type of the unwrapped basis curve.
    pub fn get_type(&self) -> Option<CurveType> {
        self.curve.as_ref().map(|c| type2d(c.as_ref()))
    }

    /// `Circle()` (`hxx:213`): the `gp_Circ2d` of a loaded circle.
    pub fn circle(&self) -> Option<GpCirc2d> {
        self.curve.as_ref().and_then(|c| basis2d(c.as_ref()).gp_circ2d())
    }

    /// `Line()` (`Geom2dAdaptor_Curve::Line`): the `gp_Lin2d` of a loaded line.
    pub fn line(&self) -> Option<occt_core::gp::GpLin2d> {
        self.curve.as_ref().and_then(|c| basis2d(c.as_ref()).gp_lin2d())
    }
}
