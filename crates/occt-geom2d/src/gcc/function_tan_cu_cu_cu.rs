//! `Geom2dGcc_FunctionTanCuCuCu`
//! (`Geom2dGcc_FunctionTanCuCuCu.{hxx:34-114,cxx:27-277}`): the three-equation
//! system of the general-curve branches of `Geom2dGcc_Circ2d3TanIter`.
//!
//! Unknowns are `X = (u1, u2, u3)`, the parameters of the three operands. The
//! residual of each pair is the sum of the two signed normal projections of the
//! chord on the two operands' normals:
//!
//! ```text
//! F1 = s1 * nor1.p1p2 + s2 * nor2.p1p2
//! F2 = s2 * nor2.p2p3 + s3 * nor3.p2p3
//! F3 = s3 * nor3.p3p1 + s1 * nor1.p3p1
//! ```
//!
//! where `nor` is the unit left normal of the operand tangent and `s` its sign
//! relative to the centroid.
//!
//! Translation notes:
//! - `Derivatives` and `Values` are declared but return `true` **without
//!   touching `D`** (`cxx:270-277`); the port keeps the empty bodies,
//! - the three point/point/point overloads declared at `hxx:63-74` have **no
//!   definition** in the `.cxx`, so they cannot be linked and are not ported
//!   (see UNPORTED),
//! - `gp::Resolution()` is `REAL_SMALL` (`gp.hxx:60`), not the repository's
//!   `precision::RESOLUTION`,
//! - `InitDerivative`'s `default:` arm throws `Standard_ConstructionError`
//!   (`cxx:81-83`); the port panics with that name,
//! - the zero-normal and zero-chord guards set the whole vector to
//!   `gp_XY(0., 0.)` (`cxx:179-201`, `:216-237`), which leaves the paired
//!   residual equal to whatever the other normal contributes.
//!
//! UNPORTED: the three `gp_Pnt2d` overloads declared at `hxx:63-74`
//! (`(gp_Circ2d, Adaptor, gp_Pnt2d)`, `(gp_Lin2d, Adaptor, gp_Pnt2d)`,
//! `(Adaptor, gp_Pnt2d, gp_Pnt2d)`) have no `.cxx` definition, so no OCCT
//! object can be constructed through them and there is no control flow to port.

use occt_core::elib::clib2d::{circle_d2_ax22d, line_d1_ax2d};
use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d, GpVec2d, GpXY};
use occt_core::precision::REAL_SMALL;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::geom2d_gcc_type::Geom2dGccType1;

/// `InitDerivative`'s nine outputs (`cxx:27-38`).
struct DerivativeState {
    point1: GpPnt2d,
    point2: GpPnt2d,
    point3: GpPnt2d,
    tan1: GpVec2d,
    tan2: GpVec2d,
    tan3: GpVec2d,
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `if (Nor >= gp::Resolution()) { out = v / Nor; } else { out = gp_XY(0., 0.); }`
/// (`cxx:179-201`, `:216-237`).
fn normalized_or_zero(v: &GpXY) -> GpXY {
    let n = v.modulus();
    if n >= REAL_SMALL {
        v.divided(n)
    } else {
        GpXY::zero()
    }
}

/// `Geom2dGcc_FunctionTanCuCuCu` (`Geom2dGcc_FunctionTanCuCuCu.hxx:34-114`).
pub struct Geom2dGccFunctionTanCuCuCu {
    /// `Curv1` (`hxx:103`).
    curv1: Geom2dAdaptorCurve,
    /// `Curv2` (`hxx:104`).
    curv2: Geom2dAdaptorCurve,
    /// `Curv3` (`hxx:105`).
    curv3: Geom2dAdaptorCurve,
    /// `Circ1` (`hxx:106`).
    circ1: GpCirc2d,
    /// `Circ2` (`hxx:107`).
    circ2: GpCirc2d,
    /// `Lin1` (`hxx:108`).
    lin1: GpLin2d,
    /// `Lin2` (`hxx:109`).
    lin2: GpLin2d,
    /// `TheType` (`hxx:110`).
    the_type: Geom2dGccType1,
}

impl Geom2dGccFunctionTanCuCuCu {
    /// `(C1, C2, C3)` (`cxx:86-95`).
    pub fn new_cu_cu_cu(
        c1: &Geom2dAdaptorCurve,
        c2: &Geom2dAdaptorCurve,
        c3: &Geom2dAdaptorCurve,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: c2.clone(),
            curv3: c3.clone(),
            circ1: GpCirc2d::default(),
            circ2: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            lin2: GpLin2d::default(),
            the_type: Geom2dGccType1::CuCuCu,
        }
    }

    /// `(gp_Circ2d C1, C2, C3)` (`cxx:97-104`).
    pub fn new_ci_cu_cu(
        c1: &GpCirc2d,
        c2: &Geom2dAdaptorCurve,
        c3: &Geom2dAdaptorCurve,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            curv3: c3.clone(),
            circ1: *c1,
            circ2: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            lin2: GpLin2d::default(),
            the_type: Geom2dGccType1::CiCuCu,
        }
    }

    /// `(gp_Circ2d C1, gp_Circ2d C2, C3)` (`cxx:106-114`).
    pub fn new_ci_ci_cu(c1: &GpCirc2d, c2: &GpCirc2d, c3: &Geom2dAdaptorCurve) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: Geom2dAdaptorCurve::new(),
            curv3: c3.clone(),
            circ1: *c1,
            circ2: *c2,
            lin1: GpLin2d::default(),
            lin2: GpLin2d::default(),
            the_type: Geom2dGccType1::CiCiCu,
        }
    }

    /// `(gp_Circ2d C1, gp_Lin2d L2, C3)` (`cxx:116-124`): the line goes to
    /// `Lin2`, which is the field the `CiLiCu` arm of `InitDerivative` reads
    /// (`cxx:53`).
    pub fn new_ci_li_cu(c1: &GpCirc2d, l2: &GpLin2d, c3: &Geom2dAdaptorCurve) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: Geom2dAdaptorCurve::new(),
            curv3: c3.clone(),
            circ1: *c1,
            circ2: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            lin2: *l2,
            the_type: Geom2dGccType1::CiLiCu,
        }
    }

    /// `(gp_Lin2d L1, gp_Lin2d L2, C3)` (`cxx:126-134`).
    pub fn new_li_li_cu(l1: &GpLin2d, l2: &GpLin2d, c3: &Geom2dAdaptorCurve) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: Geom2dAdaptorCurve::new(),
            curv3: c3.clone(),
            circ1: GpCirc2d::default(),
            circ2: GpCirc2d::default(),
            lin1: *l1,
            lin2: *l2,
            the_type: Geom2dGccType1::LiLiCu,
        }
    }

    /// `(gp_Lin2d L1, C2, C3)` (`cxx:136-144`).
    pub fn new_li_cu_cu(l1: &GpLin2d, c2: &Geom2dAdaptorCurve, c3: &Geom2dAdaptorCurve) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            curv3: c3.clone(),
            circ1: GpCirc2d::default(),
            circ2: GpCirc2d::default(),
            lin1: *l1,
            lin2: GpLin2d::default(),
            the_type: Geom2dGccType1::LiCuCu,
        }
    }

    /// `TheType` (`hxx:110`).
    pub fn the_type(&self) -> Geom2dGccType1 {
        self.the_type
    }

    /// `InitDerivative(X, Point1, Point2, Point3, Tan1, Tan2, Tan3, D21, D22,
    /// D23)` (`cxx:27-84`): the six branches that have a constructor.
    fn init_derivative(&self, x: &MathVector) -> DerivativeState {
        let (u1, u2, u3) = (x.value(1), x.value(2), x.value(3));
        let mut st = DerivativeState {
            point1: GpPnt2d::zero(),
            point2: GpPnt2d::zero(),
            point3: GpPnt2d::zero(),
            tan1: GpVec2d::zero(),
            tan2: GpVec2d::zero(),
            tan3: GpVec2d::zero(),
        };
        match self.the_type {
            Geom2dGccType1::CiCuCu => {
                // `cxx:41-45`: the circle also yields a second derivative,
                // which no residual reads.
                let (p, t, _) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                let (p, t, _) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
            Geom2dGccType1::CiCiCu => {
                // `cxx:46-51`.
                let (p, t, _) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                let (p, t, _) = circle_d2_ax22d(u2, self.circ2.position(), self.circ2.radius());
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
            Geom2dGccType1::CiLiCu => {
                // `cxx:52-58`.
                let (p, t, _) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                let (p, t) = line_d1_ax2d(u2, self.lin2.position());
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
            Geom2dGccType1::LiCuCu => {
                // `cxx:59-65`.
                let (p, t) = line_d1_ax2d(u1, self.lin1.position());
                st.point1 = p;
                st.tan1 = t;
                let (p, t, _) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
            Geom2dGccType1::LiLiCu => {
                // `cxx:66-73`.
                let (p, t) = line_d1_ax2d(u1, self.lin1.position());
                st.point1 = p;
                st.tan1 = t;
                let (p, t) = line_d1_ax2d(u2, self.lin2.position());
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
            Geom2dGccType1::CuCuCu => {
                // `cxx:74-80`.
                let (p, t, _) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                let (p, t, _) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                let (p, t, _) = curve_tool::d2(&self.curv3, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
            }
        }
        st
    }

    /// The three residuals of `Value` (`cxx:171-266`).
    fn equations(&self, st: &DerivativeState) -> [f64; 3] {
        let p1p2 = vec_between(&st.point1, &st.point2).xy();
        let p2p3 = vec_between(&st.point2, &st.point3).xy();
        let p3p1 = vec_between(&st.point3, &st.point1).xy();
        let p1p2 = normalized_or_zero(&p1p2);
        let p2p3 = normalized_or_zero(&p2p3);
        let p3p1 = normalized_or_zero(&p3p1);
        // `Nor1(-Deriv1u1.Y(), Deriv1u1.X())` etc. (`cxx:211-213`).
        let nor1 = normalized_or_zero(&GpXY::new(-st.tan1.y(), st.tan1.x()));
        let nor2 = normalized_or_zero(&GpXY::new(-st.tan2.y(), st.tan2.x()));
        let nor3 = normalized_or_zero(&GpXY::new(-st.tan3.y(), st.tan3.x()));
        // `Pcenter` is the centroid of the three points (`cxx:244`).
        let center = GpXY::new(
            (st.point1.x() + st.point2.x() + st.point3.x()) / 3.0,
            (st.point1.y() + st.point2.y() + st.point3.y()) / 3.0,
        );
        let fic1 = center.subtracted(st.point1.xy());
        let fic2 = center.subtracted(st.point2.xy());
        let fic3 = center.subtracted(st.point3.xy());
        let signe1 = if nor1.dot(&fic1) <= 0.0 { -1.0 } else { 1.0 };
        let signe2 = if nor2.dot(&fic2) <= 0.0 { -1.0 } else { 1.0 };
        let signe3 = if nor3.dot(&fic3) <= 0.0 { -1.0 } else { 1.0 };
        [
            signe1 * nor1.dot(&p1p2) + signe2 * nor2.dot(&p1p2),
            signe2 * nor2.dot(&p2p3) + signe3 * nor3.dot(&p2p3),
            signe3 * nor3.dot(&p3p1) + signe1 * nor1.dot(&p3p1),
        ]
    }
}

impl MathFunctionSetWithDerivatives for Geom2dGccFunctionTanCuCuCu {
    /// `NbVariables()` (`cxx:149-152`): 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:154-157`): 3.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, Fval)` (`cxx:159-268`): always `true`.
    fn value(&mut self, x: &MathVector, fval: &mut MathVector) -> bool {
        let st = self.init_derivative(x);
        for (i, v) in self.equations(&st).iter().enumerate() {
            fval.set_value(i + 1, *v);
        }
        true
    }

    /// `Derivatives(X, D)` (`cxx:270-273`): returns `true` and leaves `D`
    /// untouched.
    fn derivatives(&mut self, _x: &MathVector, _d: &mut MathMatrix) -> bool {
        true
    }

    /// `Values(X, F, D)` (`cxx:275-278`): returns `true` and leaves both
    /// outputs untouched.
    fn values(&mut self, _x: &MathVector, _f: &mut MathVector, _d: &mut MathMatrix) -> bool {
        true
    }
}
