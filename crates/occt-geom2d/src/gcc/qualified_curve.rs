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
//! `Handle(Geom2d_Curve)` as `Arc<dyn Curve2d>`, so a single wrapper serves
//! both roles and the type alias keeps the OCCT names at the call sites.

use std::sync::Arc;

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
