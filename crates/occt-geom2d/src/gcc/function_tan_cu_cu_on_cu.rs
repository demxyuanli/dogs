//! `Geom2dGcc_FunctionTanCuCuOnCu`
//! (`Geom2dGcc_FunctionTanCuCuOnCu.{hxx:31-140,cxx:26-405}`): the 4-equation /
//! 4-variable system of `Geom2dGcc_Circ2d2TanOnIter`.
//!
//! Unknowns are `X = (u1, u2, u3, R)`:
//!
//! ```text
//! F1 = (P3P1.P3P1 - R*R) / Rad^2
//! F2 = (P2P3.P2P3 - R*R) / Rad^2
//! F3 = P3P1.Tan1 / (|Tan1| * Rad)
//! F4 = P2P3.Tan2 / (|Tan2| * Rad)
//! ```
//!
//! `X(1..=3)` parameterize the two tangent operands and the carrier curve,
//! `X(4)` is the solution radius. `FirstRad` (`Rad`) is the radius estimate the
//! caller starts from and only scales the equations.
//!
//! Translation notes:
//! - The twelve OCCT constructors (`cxx:131-275`) are named
//!   `new_<arg1>_<arg2>_on_<carrier>`, the three suffix forms being
//!   `cu` (curve), `ci` (circle), `li` (line), `pt` (point);
//! - `cxx:240-248` (the `(gp_Circ2d, curve, circle)` constructor) sets
//!   `TheType = Geom2dGcc_CuCuOnCi` while storing its circle in `Circ1`. OCCT
//!   therefore leaves `Curv1` null in `InitDerivative`, which dereferences it at
//!   `cxx:69`; the port keeps the same `TheType` assignment
//!   ([`Geom2dGccFunctionTanCuCuOnCu::new_ci_cu_on_ci`]) and the same
//!   consequence;
//! - `p1p2` / `p2p3` / `p3p1` (`cxx:294-298`, `:323-327`, `:366-370`) are
//!   computed and never read, so they are not written here (dead stores);
//! - `Deriv(4,3)` divides by `nnor1` (`cxx:350`), not by `nnor2`; the port
//!   keeps the OCCT expression;
//! - `InitDerivative`'s `default:` arm throws `Standard_ConstructionError`
//!   (`cxx:125-127`); the port panics with that name.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::{circle_d2_ax22d, line_d1_ax2d};
use occt_core::gp::{GpCirc2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::geom2d_gcc_type::Geom2dGccType2;

/// `InitDerivative`'s nine outputs (`cxx:26-38`).
struct DerivativeState {
    point1: GpPnt2d,
    point2: GpPnt2d,
    point3: GpPnt2d,
    tan1: GpVec2d,
    tan2: GpVec2d,
    tan3: GpVec2d,
    d21: GpVec2d,
    d22: GpVec2d,
    d23: GpVec2d,
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Geom2dGcc_FunctionTanCuCuOnCu`
/// (`Geom2dGcc_FunctionTanCuCuOnCu.hxx:31-139`).
pub struct Geom2dGccFunctionTanCuCuOnCu {
    /// `Curv1` (`hxx:128`).
    curv1: Geom2dAdaptorCurve,
    /// `Curv2` (`hxx:129`).
    curv2: Geom2dAdaptorCurve,
    /// `Circ1` (`hxx:130`).
    circ1: GpCirc2d,
    /// `Lin1` (`hxx:131`).
    lin1: GpLin2d,
    /// `Pnt2` (`hxx:132`).
    pnt2: GpPnt2d,
    /// `Circon` (`hxx:133`).
    circon: GpCirc2d,
    /// `Linon` (`hxx:134`).
    linon: GpLin2d,
    /// `Curvon` (`hxx:135`).
    curvon: Geom2dAdaptorCurve,
    /// `FirstRad` (`hxx:136`).
    first_rad: f64,
    /// `TheType` (`hxx:137`).
    the_type: Geom2dGccType2,
}

impl Geom2dGccFunctionTanCuCuOnCu {
    /// `(Curve, Curve, Circle, Rad)` (`cxx:131-141`).
    pub fn new_cu_cu_on_ci(
        c1: &Geom2dAdaptorCurve,
        c2: &Geom2dAdaptorCurve,
        on_ci: &GpCirc2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: *on_ci,
            linon: GpLin2d::default(),
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuCuOnCi,
        }
    }

    /// `(Circle, Curve, Circle, Rad)` (`cxx:238-248`).
    ///
    /// OCCT stores `TheType = Geom2dGcc_CuCuOnCi` here (not `CiCuOnCi`), so
    /// `InitDerivative` takes the `CuCuOnCi` arm and reads the null `Curv1`;
    /// the port reproduces the assignment.
    pub fn new_ci_cu_on_ci(
        c1: &GpCirc2d,
        c2: &Geom2dAdaptorCurve,
        on_ci: &GpCirc2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: *c1,
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: *on_ci,
            linon: GpLin2d::default(),
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuCuOnCi,
        }
    }

    /// `(Line, Curve, Circle, Rad)` (`cxx:249-259`).
    pub fn new_li_cu_on_ci(
        l1: &GpLin2d,
        c2: &Geom2dAdaptorCurve,
        on_ci: &GpCirc2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: *l1,
            pnt2: GpPnt2d::zero(),
            circon: *on_ci,
            linon: GpLin2d::default(),
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::LiCuOnCi,
        }
    }

    /// `(Curve, Point, Circle, Rad)` (`cxx:262-272`).
    pub fn new_cu_pt_on_ci(
        c1: &Geom2dAdaptorCurve,
        p2: &GpPnt2d,
        on_ci: &GpCirc2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: Geom2dAdaptorCurve::new(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: *p2,
            circon: *on_ci,
            linon: GpLin2d::default(),
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuPtOnCi,
        }
    }

    /// `(Curve, Curve, Line, Rad)` (`cxx:178-188`).
    pub fn new_cu_cu_on_li(
        c1: &Geom2dAdaptorCurve,
        c2: &Geom2dAdaptorCurve,
        on_li: &GpLin2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: *on_li,
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuCuOnLi,
        }
    }

    /// `(Circle, Curve, Line, Rad)` (`cxx:191-200`).
    pub fn new_ci_cu_on_li(
        c1: &GpCirc2d,
        c2: &Geom2dAdaptorCurve,
        on_li: &GpLin2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: *c1,
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: *on_li,
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CiCuOnLi,
        }
    }

    /// `(Line, Curve, Line, Rad)` (`cxx:201-211`).
    pub fn new_li_cu_on_li(
        l1: &GpLin2d,
        c2: &Geom2dAdaptorCurve,
        on_li: &GpLin2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: *l1,
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: *on_li,
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::LiCuOnLi,
        }
    }

    /// `(Curve, Point, Line, Rad)` (`cxx:214-224`).
    pub fn new_cu_pt_on_li(
        c1: &Geom2dAdaptorCurve,
        p2: &GpPnt2d,
        on_li: &GpLin2d,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: Geom2dAdaptorCurve::new(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: *p2,
            circon: GpCirc2d::default(),
            linon: *on_li,
            curvon: Geom2dAdaptorCurve::new(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuPtOnLi,
        }
    }

    /// `(Curve, Curve, Curve, Rad)` (`cxx:130-141`).
    pub fn new_cu_cu_on_cu(
        c1: &Geom2dAdaptorCurve,
        c2: &Geom2dAdaptorCurve,
        c3: &Geom2dAdaptorCurve,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: GpLin2d::default(),
            curvon: c3.clone(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuCuOnCu,
        }
    }

    /// `(Circle, Curve, Curve, Rad)` (`cxx:143-152`).
    pub fn new_ci_cu_on_cu(
        c1: &GpCirc2d,
        c2: &Geom2dAdaptorCurve,
        c3: &Geom2dAdaptorCurve,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: *c1,
            lin1: GpLin2d::default(),
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: GpLin2d::default(),
            curvon: c3.clone(),
            first_rad: rad,
            the_type: Geom2dGccType2::CiCuOnCu,
        }
    }

    /// `(Line, Curve, Curve, Rad)` (`cxx:154-164`).
    pub fn new_li_cu_on_cu(
        l1: &GpLin2d,
        c2: &Geom2dAdaptorCurve,
        c3: &Geom2dAdaptorCurve,
        rad: f64,
    ) -> Self {
        Self {
            curv1: Geom2dAdaptorCurve::new(),
            curv2: c2.clone(),
            circ1: GpCirc2d::default(),
            lin1: *l1,
            pnt2: GpPnt2d::zero(),
            circon: GpCirc2d::default(),
            linon: GpLin2d::default(),
            curvon: c3.clone(),
            first_rad: rad,
            the_type: Geom2dGccType2::LiCuOnCu,
        }
    }

    /// `(Curve, Point, Curve, Rad)` (`cxx:165-176`).
    pub fn new_cu_pt_on_cu(
        c1: &Geom2dAdaptorCurve,
        p2: &GpPnt2d,
        c3: &Geom2dAdaptorCurve,
        rad: f64,
    ) -> Self {
        Self {
            curv1: c1.clone(),
            curv2: Geom2dAdaptorCurve::new(),
            circ1: GpCirc2d::default(),
            lin1: GpLin2d::default(),
            pnt2: *p2,
            circon: GpCirc2d::default(),
            linon: GpLin2d::default(),
            curvon: c3.clone(),
            first_rad: rad,
            the_type: Geom2dGccType2::CuPtOnCu,
        }
    }

    /// `FirstRad` (`hxx:136`).
    pub fn first_rad(&self) -> f64 {
        self.first_rad
    }

    /// `TheType` (`hxx:137`).
    pub fn the_type(&self) -> Geom2dGccType2 {
        self.the_type
    }

    /// `InitDerivative(X, Point1, Point2, Point3, Tan1, Tan2, Tan3, D21, D22,
    /// D23)` (`cxx:26-129`): evaluates the three points with their first and
    /// second derivatives according to `TheType`.
    fn init_derivative(&self, x: &MathVector) -> DerivativeState {
        let (u1, u2, u3) = (x.value(1), x.value(2), x.value(3));
        let zero = GpVec2d::zero();
        let mut st = DerivativeState {
            point1: GpPnt2d::zero(),
            point2: GpPnt2d::zero(),
            point3: GpPnt2d::zero(),
            tan1: zero,
            tan2: zero,
            tan3: zero,
            d21: zero,
            d22: zero,
            d23: zero,
        };
        match self.the_type {
            Geom2dGccType2::CuCuOnCu => {
                // `cxx:40-44`.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = curve_tool::d2(&self.curvon, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::CiCuOnCu => {
                // `cxx:45-50`.
                let (p, t, n) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = curve_tool::d2(&self.curvon, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::LiCuOnCu => {
                // `cxx:51-57`: `D21` is the zero vector.
                let (p, t) = line_d1_ax2d(u1, self.lin1.position());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = zero;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = curve_tool::d2(&self.curvon, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::CuPtOnCu => {
                // `cxx:58-66`: `Point2` is the constant `Pnt2`, `Tan2` and
                // `D22` are zero.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curvon, u3).expect("undefined curve value");
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
                st.point2 = self.pnt2;
                st.tan2 = zero;
                st.d22 = zero;
            }
            Geom2dGccType2::CuCuOnCi => {
                // `cxx:67-70`.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = circle_d2_ax22d(u3, self.circon.position(), self.circon.radius());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::CiCuOnCi => {
                // `cxx:71-76`.
                let (p, t, n) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = circle_d2_ax22d(u3, self.circon.position(), self.circon.radius());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::LiCuOnCi => {
                // `cxx:77-85`: `D21` is the zero vector.
                let (p, t) = line_d1_ax2d(u1, self.lin1.position());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = zero;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t, n) = circle_d2_ax22d(u3, self.circon.position(), self.circon.radius());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::CuPtOnCi => {
                // `cxx:86-92`.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                st.point2 = self.pnt2;
                st.tan2 = zero;
                st.d22 = zero;
                let (p, t, n) = circle_d2_ax22d(u3, self.circon.position(), self.circon.radius());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = n;
            }
            Geom2dGccType2::CuCuOnLi => {
                // `cxx:93-99`: `D23` is the zero vector.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t) = line_d1_ax2d(u3, self.linon.position());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = zero;
            }
            Geom2dGccType2::CiCuOnLi => {
                // `cxx:100-107`.
                let (p, t, n) = circle_d2_ax22d(u1, self.circ1.position(), self.circ1.radius());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t) = line_d1_ax2d(u3, self.linon.position());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = zero;
            }
            Geom2dGccType2::LiCuOnLi => {
                // `cxx:108-115`: `D21` and `D23` are zero vectors.
                let (p, t) = line_d1_ax2d(u1, self.lin1.position());
                st.point1 = p;
                st.tan1 = t;
                st.d21 = zero;
                let (p, t, n) = curve_tool::d2(&self.curv2, u2).expect("undefined curve value");
                st.point2 = p;
                st.tan2 = t;
                st.d22 = n;
                let (p, t) = line_d1_ax2d(u3, self.linon.position());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = zero;
            }
            Geom2dGccType2::CuPtOnLi => {
                // `cxx:116-124`.
                let (p, t, n) = curve_tool::d2(&self.curv1, u1).expect("undefined curve value");
                st.point1 = p;
                st.tan1 = t;
                st.d21 = n;
                st.point2 = self.pnt2;
                st.tan2 = zero;
                st.d22 = zero;
                let (p, t) = line_d1_ax2d(u3, self.linon.position());
                st.point3 = p;
                st.tan3 = t;
                st.d23 = zero;
            }
        }
        st
    }

    /// The four residuals of `Value` / `Values` (`cxx:304-308`, `:377-380`).
    fn equations(&self, x: &MathVector, st: &DerivativeState) -> [f64; 4] {
        let p2p3 = vec_between(&st.point2, &st.point3);
        let p3p1 = vec_between(&st.point3, &st.point1);
        let nnor1 = st.tan1.magnitude();
        let nnor2 = st.tan2.magnitude();
        let rad2 = self.first_rad * self.first_rad;
        let r = x.value(4);
        let rad = self.first_rad;
        [
            (p3p1.dot(&p3p1) - r * r) / rad2,
            (p2p3.dot(&p2p3) - r * r) / rad2,
            p3p1.dot(&st.tan1) / (nnor1 * rad),
            p2p3.dot(&st.tan2) / (nnor2 * rad),
        ]
    }

    /// The sixteen partials of `Derivatives` / `Values` (`cxx:330-352`,
    /// `:383-403`), row-major `D[i][j]` for `i, j` in `1..=4`.
    fn jacobian(&self, x: &MathVector, st: &DerivativeState) -> [[f64; 4]; 4] {
        let p2p3 = vec_between(&st.point2, &st.point3);
        let p3p1 = vec_between(&st.point3, &st.point1);
        let nnor1 = st.tan1.magnitude();
        let nnor2 = st.tan2.magnitude();
        let rad2 = self.first_rad * self.first_rad;
        let rad = self.first_rad;
        let r = x.value(4);
        let mut d = [[0.0_f64; 4]; 4];
        // Row 1.
        d[0][0] = 2. * st.tan1.dot(&p3p1) / rad2;
        d[0][1] = 0.;
        d[0][2] = -2. * st.tan3.dot(&p3p1) / rad2;
        d[0][3] = -2. * r / rad2;
        // Row 2.
        d[1][0] = 0.;
        d[1][1] = -2. * st.tan2.dot(&p2p3) / rad2;
        d[1][2] = 2. * st.tan3.dot(&p2p3) / rad2;
        d[1][3] = -2. * r / rad2;
        // Row 3.
        d[2][0] = (p3p1.dot(&st.d21) + st.tan1.dot(&st.tan1)) / (rad * nnor1)
            - (p3p1.dot(&st.tan1) * st.d21.dot(&st.tan1)) / (rad * nnor1 * nnor1 * nnor1);
        d[2][1] = 0.;
        d[2][2] = -(st.tan3.dot(&st.tan1)) / (rad * nnor1);
        d[2][3] = 0.;
        // Row 4. OCCT divides the last partial by `nnor1` (`cxx:350`).
        d[3][0] = 0.;
        d[3][1] = (p2p3.dot(&st.d22) - st.tan2.dot(&st.tan2)) / (rad * nnor2)
            - p2p3.dot(&st.tan2) * st.tan2.dot(&st.d22) / (rad * nnor2 * nnor2 * nnor2);
        d[3][2] = st.tan3.dot(&st.tan2) / (rad * nnor1);
        d[3][3] = 0.;
        d
    }
}

impl MathFunctionSetWithDerivatives for Geom2dGccFunctionTanCuCuOnCu {
    /// `NbVariables()` (`cxx:276-279`): 4.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:281-284`): 4.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, Fval)` (`cxx:286-309`): always `true`.
    fn value(&mut self, x: &MathVector, fval: &mut MathVector) -> bool {
        let st = self.init_derivative(x);
        for (i, v) in self.equations(x, &st).iter().enumerate() {
            fval.set_value(i + 1, *v);
        }
        true
    }

    /// `Derivatives(X, Deriv)` (`cxx:311-353`): always `true`.
    fn derivatives(&mut self, x: &MathVector, deriv: &mut MathMatrix) -> bool {
        let st = self.init_derivative(x);
        let jac = self.jacobian(x, &st);
        for i in 0..4 {
            for j in 0..4 {
                deriv.set_value(i + 1, j + 1, jac[i][j]);
            }
        }
        true
    }

    /// `Values(X, Fval, Deriv)` (`cxx:355-405`): `Value` and `Derivatives` in
    /// a single evaluation of `InitDerivative`.
    fn values(&mut self, x: &MathVector, fval: &mut MathVector, deriv: &mut MathMatrix) -> bool {
        let st = self.init_derivative(x);
        for (i, v) in self.equations(x, &st).iter().enumerate() {
            fval.set_value(i + 1, *v);
        }
        let jac = self.jacobian(x, &st);
        for i in 0..4 {
            for j in 0..4 {
                deriv.set_value(i + 1, j + 1, jac[i][j]);
            }
        }
        true
    }
}
