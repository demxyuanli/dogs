//! `Geom2dGcc_QCurve` and `Geom2dGcc_QualifiedCurve`: a 2D curve plus the
//! position the requested solution must take with respect to it.
//!
//! | OCCT | lines |
//! |---|---|
//! | `Geom2dGcc_QualifiedCurve` | `Geom2dGcc_QualifiedCurve.cxx:17-53` |
//! | `Geom2dGcc_QCurve` | `Geom2dGcc_QCurve.cxx:18-53` |
//!
//! The two classes are structurally identical in OCCT; they differ only in the
//! adaptor they wrap (`Geom2dAdaptor_Curve` in both cases, but
//! `QualifiedCurve` is handed to the solvers that take an unbounded curve and
//! `QCurve` to the ones that take a curve with a range). The port stores
//! `Handle(Geom2d_Curve)` as `Arc<dyn Curve2d>` in
//! [`Geom2dGccQualifiedCurve`], and keeps the [`Geom2dAdaptorCurve`] itself in
//! [`Geom2dGccQCurve`], whose consumers (`Geom2dGcc_CurveTool`) need the range
//! and the parametric resolution.

use std::sync::Arc;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::curve::Curve2d;
use crate::gcc::gcc_ent::GccEntPosition;

/// `Handle(Geom2dAdaptor_Curve)` (`Geom2dAdaptor_Curve.hxx`): the adaptor owns
/// the curve and answers `GetType()` through
/// [`crate::adaptor_curve::type2d`], which is why the port keeps the raw curve.
pub type QualifiedCurve2d = Arc<dyn Curve2d>;

/// `Geom2dGcc_QualifiedCurve` (`cxx:17-53`) / `Geom2dGcc_QCurve`
/// (`cxx:18-53`).
#[derive(Clone)]
pub struct Geom2dGccQualifiedCurve {
    the_qualified: QualifiedCurve2d,
    the_qualifier: GccEntPosition,
}

impl Geom2dGccQualifiedCurve {
    /// `Geom2dGcc_QualifiedCurve(Curve, Qualifier)` (`cxx:17-21`).
    pub fn new(curve: QualifiedCurve2d, qualifier: GccEntPosition) -> Self {
        Self { the_qualified: curve, the_qualifier: qualifier }
    }

    /// `Qualified()` (`cxx:23-26`).
    pub fn qualified(&self) -> &QualifiedCurve2d {
        &self.the_qualified
    }

    /// `Qualifier()` (`cxx:28-31`).
    pub fn qualifier(&self) -> GccEntPosition {
        self.the_qualifier
    }

    /// `IsUnqualified()` (`cxx:33-36`).
    pub fn is_unqualified(&self) -> bool {
        self.the_qualifier == GccEntPosition::Unqualified
    }

    /// `IsEnclosing()` (`cxx:38-41`).
    pub fn is_enclosing(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosing
    }

    /// `IsEnclosed()` (`cxx:43-46`).
    pub fn is_enclosed(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosed
    }

    /// `IsOutside()` (`cxx:48-51`).
    pub fn is_outside(&self) -> bool {
        self.the_qualifier == GccEntPosition::Outside
    }
}

/// `Geom2dGcc_QCurve` (`Geom2dGcc_QCurve.{hxx:28-60,cxx:18-53}`): a
/// `Geom2dAdaptor_Curve` plus the position the requested solution must take
/// with respect to it.
///
/// Unlike [`Geom2dGccQualifiedCurve`] this keeps the adaptor, because the
/// `Geom2dGcc_Circ2d2TanOnIter` solvers read `FirstParameter` /
/// `LastParameter` / `Resolution` off `Qualified()` (`Geom2dGcc_CurveTool`).
#[derive(Clone)]
pub struct Geom2dGccQCurve {
    the_qualified: Geom2dAdaptorCurve,
    the_qualifier: GccEntPosition,
}

impl Geom2dGccQCurve {
    /// `Geom2dGcc_QCurve(Curve, Qualifier)` (`cxx:18-21`).
    pub fn new(curve: Geom2dAdaptorCurve, qualifier: GccEntPosition) -> Self {
        Self { the_qualified: curve, the_qualifier: qualifier }
    }

    /// `Qualified()` (`cxx:23-26`).
    pub fn qualified(&self) -> &Geom2dAdaptorCurve {
        &self.the_qualified
    }

    /// `Qualifier()` (`cxx:28-31`).
    pub fn qualifier(&self) -> GccEntPosition {
        self.the_qualifier
    }

    /// `IsUnqualified()` (`cxx:33-36`).
    pub fn is_unqualified(&self) -> bool {
        self.the_qualifier == GccEntPosition::Unqualified
    }

    /// `IsEnclosing()` (`cxx:38-41`).
    pub fn is_enclosing(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosing
    }

    /// `IsEnclosed()` (`cxx:43-46`).
    pub fn is_enclosed(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosed
    }

    /// `IsOutside()` (`cxx:48-51`).
    pub fn is_outside(&self) -> bool {
        self.the_qualifier == GccEntPosition::Outside
    }
}
