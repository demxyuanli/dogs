//! Helper bodies shared by the `Geom2dGcc_*Geo` solvers whose centre locus is
//! an operand offset by `-+Radius`.
//!
//! OCCT has no shared class for these blocks: `Geom2dGcc_Circ2dTanOnRadGeo.cxx`
//! (`:136-143`, `:252-267`, `:391-396`, `:675-692`) and
//! `Geom2dGcc_Circ2d2TanRadGeo.cxx` (`:199-206`, `:408-421`, `:571-580`) repeat
//! the axis construction, the offset-curve wrapper and the two
//! `IntRes2d_Domain` builders verbatim. The port keeps one copy here so the two
//! solvers cannot drift apart; every item names the OCCT lines it comes from.
//!
//! UNPORTED: none.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::clib2d::{circle_value_ax22d, line_value_ax2d};
use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intres2d::IntRes2dDomain;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::curve::Curve2d;
use crate::offset::Geom2dOffsetCurve;

/// `thefirst = -100000.` (`Circ2dTanOnRadGeo.cxx:88`, `Circ2d2TanRadGeo.cxx:79`).
pub const THEFIRST: f64 = -100000.0;

/// `thelast = 100000.` (`Circ2dTanOnRadGeo.cxx:89`, `Circ2d2TanRadGeo.cxx:80`).
pub const THELAST: f64 = 100000.0;

/// `gp_Ax2d(Center, gp_Dir2d::D::X)`.
pub fn axis_at(center: GpPnt2d) -> GpAx22d {
    GpAx22d::from_xdir(center, GpDir2d::new(1.0, 0.0).expect("gp::DX2d()"))
}

/// `gp_Circ2d(gp_Ax2d(Center, dirx), Radius)` (`Cir2dTanOnRadGeo.cxx:155`,
/// `Circ2d2TanRadGeo.cxx:218`).
pub fn circle_at(center: GpPnt2d, radius: f64) -> GpCirc2d {
    GpCirc2d::new(axis_at(center), radius)
}

/// `Adaptor2d_OffsetCurve(HCu, Offset)` where the basis is the adaptor's curve
/// (`Circ2dTanOnRadGeo.cxx:136`, `Circ2d2TanRadGeo.cxx:200`).
pub fn offset_of(adaptor: &Geom2dAdaptorCurve, offset: f64) -> Arc<dyn Curve2d> {
    let basis = adaptor
        .curve()
        .expect("Geom2dAdaptor_Curve::Curve on the operand")
        .clone();
    Arc::new(Geom2dOffsetCurve::new(basis, offset))
}

/// `IntRes2d_Domain D2(C2.Value(firstparam), firstparam, Tol,
/// C2.Value(lastparam), lastparam, Tol)` with
/// `firstparam = max(C2.FirstParameter(), thefirst)` and
/// `lastparam = min(C2.LastParameter(), thelast)`
/// (`Circ2dTanOnRadGeo.cxx:137-143`, `Circ2d2TanRadGeo.cxx:201-207`).
pub fn curve_domain(curve: &dyn Curve2d, tol: f64) -> IntRes2dDomain {
    let first = curve.first_parameter().max(THEFIRST);
    let last = curve.last_parameter().min(THELAST);
    let p_first = curve.d0(first);
    let p_last = curve.d0(last);
    IntRes2dDomain::bounded(&p_first, first, tol, &p_last, last, tol)
}

/// `IntRes2d_Domain D1(ElCLib::Value(0., Circ), 0., Tol, ElCLib::Value(2.*M_PI,
/// Circ), 2.*M_PI, Tol); D1.SetEquivalentParameters(0., 2.*M_PI);`
/// (`Circ2dTanOnRadGeo.cxx:254-260`, `Circ2d2TanRadGeo.cxx:408-410`).
pub fn circle_domain(circle: &GpCirc2d, tol: f64) -> IntRes2dDomain {
    let first = circle_value_ax22d(0.0, circle.position(), circle.radius());
    let last = circle_value_ax22d(2.0 * PI, circle.position(), circle.radius());
    let mut domain = IntRes2dDomain::bounded(&first, 0.0, tol, &last, 2.0 * PI, tol);
    domain.set_equivalent_parameters(0.0, 2.0 * PI);
    domain
}

/// `ElCLib::Value(Par, C1)` for a circle (`Circ2dTanOnRadGeo.cxx:437`).
pub fn value_on_circ(circle: &GpCirc2d, par: f64) -> GpPnt2d {
    circle_value_ax22d(par, circle.position(), circle.radius())
}

/// `ElCLib::Value(Par, L1)` for a line (`Circ2d2TanRadGeo.cxx:238`).
pub fn value_on_line(line: &GpLin2d, par: f64) -> GpPnt2d {
    line_value_ax2d(par, line.position())
}
