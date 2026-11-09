//! `Geom2dGcc` (`Geom2dGcc.{hxx:40-114,cxx:19-37}`): the package facade with
//! the four "qualified curve" builders.
//!
//! Every builder returns a `Geom2dGcc_QualifiedCurve` that carries the adaptor
//! and fixed qualifier; the solvers that take a `Geom2dGcc_QualifiedCurve`
//! (`Geom2dGcc_Circ2d2TanRad`, `Geom2dGcc_Circ2dTanCen`, `Geom2dGcc_Lin2d2Tan`,
//! `Geom2dGcc_Lin2dTanObl`) are then handed it directly.
//!
//! Translation notes:
//! - OCCT keeps the `Geom2dAdaptor_Curve` by reference in the returned
//!   `Geom2dGcc_QualifiedCurve`; the port stores the curve handle
//!   (`Arc<dyn Curve2d>`) exactly like [`Geom2dGccQualifiedCurve::new`], so a
//!   caller that mutates its adaptor after the call does not affect the
//!   qualified curve in Rust. No OCCT control flow depends on that aliasing.
//!
//! UNPORTED: none.

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc::qualified_curve::Geom2dGccQualifiedCurve;

/// `Geom2dGcc` (`Geom2dGcc.hxx:40-114`).
pub struct Geom2dGcc;

impl Geom2dGcc {
    /// `Geom2dGcc::Unqualified(Obj)` (`cxx:19-22`).
    pub fn unqualified(curve: &Geom2dAdaptorCurve) -> Geom2dGccQualifiedCurve {
        qualified(curve, GccEntPosition::Unqualified)
    }

    /// `Geom2dGcc::Enclosing(Obj)` (`cxx:24-27`).
    pub fn enclosing(curve: &Geom2dAdaptorCurve) -> Geom2dGccQualifiedCurve {
        qualified(curve, GccEntPosition::Enclosing)
    }

    /// `Geom2dGcc::Enclosed(Obj)` (`cxx:29-32`).
    pub fn enclosed(curve: &Geom2dAdaptorCurve) -> Geom2dGccQualifiedCurve {
        qualified(curve, GccEntPosition::Enclosed)
    }

    /// `Geom2dGcc::Outside(Obj)` (`cxx:34-37`).
    pub fn outside(curve: &Geom2dAdaptorCurve) -> Geom2dGccQualifiedCurve {
        qualified(curve, GccEntPosition::Outside)
    }
}

/// `Geom2dGcc_QualifiedCurve(Obj, Qualifier)` (`cxx:21`, `:26`, `:31`, `:36`).
fn qualified(curve: &Geom2dAdaptorCurve, position: GccEntPosition) -> Geom2dGccQualifiedCurve {
    Geom2dGccQualifiedCurve::new(
        curve
            .curve()
            .expect("Geom2dAdaptor_Curve::Curve on the qualified curve")
            .clone(),
        position,
    )
}
