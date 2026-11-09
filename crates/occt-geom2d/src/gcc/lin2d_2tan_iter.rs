//! `Geom2dGcc_Lin2d2TanIter`
//! (`Geom2dGcc_Lin2d2TanIter.{hxx:38-103,cxx:40-338}`): one line tangent to a
//! circle/curve pair, a curve pair or a curve/point pair, solved by a local
//! root search from the caller's parametric guess.
//!
//! The three constructors back the `Geom2dGcc_Lin2d2Tan` sampling loops; each
//! one either accepts or rejects the single root it finds, so the object holds
//! at most one solution.
//!
//! | OCCT | arguments | `cxx` |
//! |---|---|---|
//! | `(QualifiedCirc, QCurve, Param2, Tolang)` | `:40-135` |
//! | `(QCurve, QCurve, Param1, Param2, Tolang)` | `:137-217` |
//! | `(QCurve, gp_Pnt2d, Param1, Tolang)` | `:219-282` |
//! | `IsDone`..`Tangency2` | `:284-338` |
//!
//! Translation notes:
//! - the circle/curve constructor throws `GccEnt_BadQualifier` for an enclosed
//!   circle (`cxx:55-58`), kept as a panic,
//! - `math_FunctionRoot(func, Param2, EpsX, U1, U2, 100)` (`cxx:65-69`) is
//!   [`MathFunctionRoot::new_with_bounds`] with `EpsX = CurveTool::EpsX(Cu,
//!   |Tolang|)`; the root is then re-evaluated and only accepted when
//!   `|Norm| < Tolang` with the **signed** `Tolang` (`cxx:77`), kept,
//! - the curve/curve constructor builds `math_FunctionSetRoot Root(Func, tol)`
//!   without `Perform`'s `WithSigma` flag and tolerances `EpsX(Cu1, |Tolang|)` /
//!   `EpsX(Cu2, |Tolang|)` (`cxx:171-175`); the port calls
//!   `perform_with_bounds(..., false)`,
//! - `Root.Root(Ufirst)` overwrites `Ufirst` in place (`cxx:178`), so the
//!   subsequent `Func.Value(Ufirst, Norm)` and every `Ufirst(i)` read use the
//!   root,
//! - `Vec.Angle(Vect12)` (`cxx:189`) angles the chord against the *second*
//!   derivative of curve 1, and `Vec.Angle(Vect22)` (`cxx:196`) the same chord
//!   against the second derivative of curve 2; kept verbatim,
//! - `Vec1.Angle(Vdir)` in the circle/curve constructor (`cxx:97-98`) and
//!   `gp_Vec2d::Angle` elsewhere return `0.` for a null operand here
//!   (`gp_Vec2d.cxx:57` throws in OCCT),
//! - `gp_Dir2d(...)` built from a vector that may be null throws
//!   `gp_VectorWithNullMagnitude` in OCCT (`cxx:98`, `:209`, `:271`); the port
//!   panics through `expect`,
//! - `pararg2 = ThePoint.Distance(Origine)` and `par2sol = 0.` in the
//!   curve/point constructor (`cxx:277-278`) are not the usual pairing (the
//!   other two constructors put the distance in `par2sol`); kept.
//!
//! UNPORTED: none - the three constructors and all five accessors.

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d, GpXY};
use occt_math::{
    FunctionWithDerivative, MathFunctionRoot, MathFunctionSetRoot,
    MathFunctionSetWithDerivatives, MathVector,
};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_cir_cu::Geom2dGccFunctionTanCirCu;
use crate::gcc::function_tan_cu_cu::Geom2dGccFunctionTanCuCu;
use crate::gcc::function_tan_cu_pnt::Geom2dGccFunctionTanCuPnt;
use crate::gcc::gcc_ent::{GccEntPosition, GccEntQualifiedCirc};
use crate::gcc::qualified_curve::Geom2dGccQCurve;

/// `nbIterations = 100` of every `math_FunctionRoot` / `math_FunctionSetRoot`
/// call in the three constructors (`cxx:69`, `:175`, `:247`).
const NB_ITERATIONS: i32 = 100;

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `gp_Vec2d(V.XY, W.XY)` (`gp_Vec2d.cxx:66-73`): `W - V`.
fn vec_between_xy(from: &GpXY, to: &GpXY) -> GpVec2d {
    GpVec2d::new(to.x - from.x, to.y - from.y)
}

/// `Geom2dGcc_Lin2d2TanIter` (`hxx:38-103`).
#[derive(Clone)]
pub struct Geom2dGccLin2d2TanIter {
    /// `WellDone` (`hxx:97`).
    well_done: bool,
    /// `qualifier1` (`hxx:98`).
    qualifier1: GccEntPosition,
    /// `qualifier2` (`hxx:99`).
    qualifier2: GccEntPosition,
    /// `linsol` (`hxx:100`).
    linsol: GpLin2d,
    /// `pnttg1sol` (`hxx:101`).
    pnttg1sol: GpPnt2d,
    /// `pnttg2sol` (`hxx:102`).
    pnttg2sol: GpPnt2d,
    /// `par1sol` (`hxx:103`).
    par1sol: f64,
    /// `pararg1` (`hxx:104`).
    pararg1: f64,
    /// `par2sol` (`hxx:105`).
    par2sol: f64,
    /// `pararg2` (`hxx:106`).
    pararg2: f64,
}

impl Geom2dGccLin2d2TanIter {
    /// `par1sol = 0.; pararg1 = 0.; par2sol = 0.0; pararg2 = 0.0;
    /// WellDone = false; qualifier1 = GccEnt_noqualifier;
    /// qualifier2 = GccEnt_noqualifier;` (`cxx:46-53`, `:143-149`, `:224-230`).
    fn not_done() -> Self {
        Self {
            well_done: false,
            qualifier1: GccEntPosition::NoQualifier,
            qualifier2: GccEntPosition::NoQualifier,
            linsol: GpLin2d::default(),
            pnttg1sol: GpPnt2d::zero(),
            pnttg2sol: GpPnt2d::zero(),
            par1sol: 0.0,
            pararg1: 0.0,
            par2sol: 0.0,
            pararg2: 0.0,
        }
    }

    /// `(GccEnt_QualifiedCirc, Geom2dGcc_QCurve, Param2, Tolang)`
    /// (`cxx:40-135`).
    pub fn new_circ_curve(
        qualified1: &GccEntQualifiedCirc,
        qualified2: &Geom2dGccQCurve,
        param2: f64,
        tolang: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if qualified1.is_enclosed() {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Lin2d2TanIter");
        }
        let c1 = qualified1.qualified();
        let cu2 = qualified2.qualified();
        let u1 = curve_tool::first_parameter(cu2);
        let u2 = curve_tool::last_parameter(cu2);
        let mut func = Geom2dGccFunctionTanCirCu::new(&c1, cu2);
        let sol = MathFunctionRoot::new_with_bounds(
            &mut func,
            param2,
            curve_tool::eps_x(cu2, tolang.abs()),
            u1,
            u2,
            NB_ITERATIONS,
        );
        if sol.is_done() {
            let usol = sol.root();
            let mut norm = 0.0f64;
            func.value(usol, &mut norm);
            if norm.abs() < tolang {
                let (origine, mut vect1, vect2) =
                    curve_tool::d2(cu2, usol).expect("Geom2dAdaptor_Curve::D2");
                let vdir = vec_between_xy(&origine.xy(), &c1.location().xy());
                let mut sign1 = vect1.dot(&vdir);
                if sign1 <= 0.0 {
                    vect1.reverse();
                }
                let sign2 = vect2.crossed(&vect1);
                if qualified2.is_unqualified()
                    || (qualified2.is_enclosing() && sign2 <= 0.0)
                    || (qualified2.is_outside() && sign1 <= 0.0 && sign2 >= 0.0)
                    || (qualified2.is_enclosed() && sign1 >= 0.0 && sign2 >= 0.0)
                {
                    if qualified1.is_unqualified()
                        || (qualified1.is_outside() && vect1.angle(&vdir) <= 0.0)
                        || (qualified1.is_enclosing() && vect1.angle(&vdir) >= 0.0)
                    {
                        let direc = GpDir2d::new(vect1.x(), vect1.y())
                            .expect("gp_VectorWithNullMagnitude in Geom2dGcc_Lin2d2TanIter");
                        let r1 = c1.radius();
                        let normal = GpXY::new(-r1 * direc.y(), r1 * direc.x());
                        let center = c1.location();
                        sign1 = vect1.crossed(&vdir);
                        if qualified1.is_enclosing() {
                            s.pnttg1sol = GpPnt2d::new(
                                center.x() - normal.x,
                                center.y() - normal.y,
                            );
                        } else if qualified1.is_outside() {
                            s.pnttg1sol = GpPnt2d::new(
                                center.x() + normal.x,
                                center.y() + normal.y,
                            );
                        } else if sign1 >= 0.0 {
                            s.pnttg1sol = GpPnt2d::new(
                                center.x() - normal.x,
                                center.y() - normal.y,
                            );
                        } else {
                            s.pnttg1sol = GpPnt2d::new(
                                center.x() + normal.x,
                                center.y() + normal.y,
                            );
                        }
                        s.pnttg2sol = origine;
                        s.linsol = GpLin2d::from_pnt_dir(s.pnttg1sol, direc);
                        s.well_done = true;
                        s.qualifier1 = qualified1.qualifier();
                        s.qualifier2 = qualified2.qualifier();
                        s.pararg2 = usol;
                        s.par1sol = 0.0;
                        s.par2sol = s.pnttg2sol.distance(&s.pnttg1sol);
                        s.pararg1 = 0.0;
                    }
                }
            }
        }
        s
    }

    /// `(Geom2dGcc_QCurve, Geom2dGcc_QCurve, Param1, Param2, Tolang)`
    /// (`cxx:137-217`).
    pub fn new_curve_curve(
        qualified1: &Geom2dGccQCurve,
        qualified2: &Geom2dGccQCurve,
        param1: f64,
        param2: f64,
        tolang: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
            || !(qualified2.is_enclosed()
                || qualified2.is_enclosing()
                || qualified2.is_outside()
                || qualified2.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Lin2d2TanIter");
        }
        let cu1 = qualified1.qualified();
        let cu2 = qualified2.qualified();
        let mut func = Geom2dGccFunctionTanCuCu::new_two_curves(cu1, cu2);
        let mut umin = MathVector::with_init(1, 2, 0.0);
        let mut umax = MathVector::with_init(1, 2, 0.0);
        let mut ufirst = MathVector::with_init(1, 2, 0.0);
        let mut tol = MathVector::with_init(1, 2, 0.0);
        umin.set_value(1, curve_tool::first_parameter(cu1));
        umin.set_value(2, curve_tool::first_parameter(cu2));
        umax.set_value(1, curve_tool::last_parameter(cu1));
        umax.set_value(2, curve_tool::last_parameter(cu2));
        ufirst.set_value(1, param1);
        ufirst.set_value(2, param2);
        tol.set_value(1, curve_tool::eps_x(cu1, tolang.abs()));
        tol.set_value(2, curve_tool::eps_x(cu2, tolang.abs()));
        let mut root = MathFunctionSetRoot::new(&func, &tol, NB_ITERATIONS);
        root.perform_with_bounds(&mut func, &ufirst, &umin, &umax, false);
        if root.is_done() {
            let x = root.root();
            let mut norm = MathVector::with_init(1, 2, 0.0);
            func.value(&x, &mut norm);
            if norm.value(1).abs() < tolang && norm.value(2).abs() < tolang {
                let (point1, vect11, vect12) =
                    curve_tool::d2(cu1, x.value(1)).expect("Geom2dAdaptor_Curve::D2");
                let (point2, vect21, vect22) =
                    curve_tool::d2(cu2, x.value(2)).expect("Geom2dAdaptor_Curve::D2");
                let vec = vec_between(&point1, &point2);
                let mut angle1 = vec.angle(&vect12);
                let mut sign1 = vect11.dot(&vec);
                if qualified1.is_unqualified()
                    || (qualified1.is_enclosing() && angle1 >= 0.0)
                    || (qualified1.is_outside() && angle1 <= 0.0 && sign1 <= 0.0)
                    || (qualified1.is_enclosed() && angle1 <= 0.0 && sign1 >= 0.0)
                {
                    angle1 = vec.angle(&vect22);
                    sign1 = vect21.dot(&vec);
                    if qualified2.is_unqualified()
                        || (qualified2.is_enclosing() && angle1 >= 0.0)
                        || (qualified2.is_outside() && angle1 <= 0.0 && sign1 <= 0.0)
                        || (qualified2.is_enclosed() && angle1 <= 0.0 && sign1 >= 0.0)
                    {
                        s.qualifier1 = qualified1.qualifier();
                        s.qualifier2 = qualified2.qualifier();
                        s.pararg1 = x.value(1);
                        s.par1sol = 0.0;
                        s.pnttg1sol = point1;
                        s.pararg2 = x.value(2);
                        s.pnttg2sol = point2;
                        s.par2sol = s.pnttg2sol.distance(&s.pnttg1sol);
                        let dir = GpDir2d::new(
                            s.pnttg2sol.x() - s.pnttg1sol.x(),
                            s.pnttg2sol.y() - s.pnttg1sol.y(),
                        )
                        .expect("gp_VectorWithNullMagnitude in Geom2dGcc_Lin2d2TanIter");
                        s.linsol = GpLin2d::from_pnt_dir(s.pnttg1sol, dir);
                        s.well_done = true;
                    }
                }
            }
        }
        s
    }

    /// `(Geom2dGcc_QCurve, gp_Pnt2d, Param1, Tolang)` (`cxx:219-282`).
    pub fn new_curve_point(
        qualified1: &Geom2dGccQCurve,
        the_point: &GpPnt2d,
        param1: f64,
        tolang: f64,
    ) -> Self {
        let mut s = Self::not_done();
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Lin2d2TanIter");
        }
        let cu1 = qualified1.qualified();
        let u1 = curve_tool::first_parameter(cu1);
        let u2 = curve_tool::last_parameter(cu1);
        let mut func = Geom2dGccFunctionTanCuPnt::new(cu1, the_point);
        let sol = MathFunctionRoot::new_with_bounds(
            &mut func,
            param1,
            curve_tool::eps_x(cu1, tolang.abs()),
            u1,
            u2,
            NB_ITERATIONS,
        );
        if sol.is_done() {
            let usol = sol.root();
            let mut norm = 0.0f64;
            func.value(usol, &mut norm);
            if norm.abs() < tolang {
                let (origine, vect1, vect2) =
                    curve_tool::d2(cu1, usol).expect("Geom2dAdaptor_Curve::D2");
                let vdir = vec_between_xy(&origine.xy(), &the_point.xy());
                let sign1 = vect1.dot(&vdir);
                let sign2 = vect2.crossed(&vdir);
                if qualified1.is_unqualified()
                    || (qualified1.is_enclosing()
                        && ((sign1 >= 0.0 && sign2 <= 0.0) || (sign1 <= 0.0 && sign2 <= 0.0)))
                    || (qualified1.is_outside() && sign1 <= 0.0 && sign2 >= 0.0)
                    || (qualified1.is_enclosed() && sign1 >= 0.0 && sign2 >= 0.0)
                {
                    s.well_done = true;
                    s.linsol = GpLin2d::from_pnt_dir(
                        origine,
                        GpDir2d::new(vdir.x(), vdir.y()).expect(
                            "gp_VectorWithNullMagnitude in Geom2dGcc_Lin2d2TanIter",
                        ),
                    );
                    s.qualifier1 = qualified1.qualifier();
                    s.qualifier2 = GccEntPosition::NoQualifier;
                    s.pnttg1sol = origine;
                    s.pnttg2sol = *the_point;
                    s.pararg1 = usol;
                    s.par1sol = 0.0;
                    s.pararg2 = the_point.distance(&origine);
                    s.par2sol = 0.0;
                }
            }
        }
        s
    }

    /// `IsDone()` (`cxx:284-287`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `ThisSolution()` (`cxx:289-296`).
    pub fn this_solution(&self) -> GpLin2d {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2d2TanIter");
        self.linsol
    }

    /// `WhichQualifier(Qualif1, Qualif2)` (`cxx:298-310`).
    pub fn which_qualifier(&self) -> (GccEntPosition, GccEntPosition) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2d2TanIter");
        (self.qualifier1, self.qualifier2)
    }

    /// `Tangency1(ParSol, ParArg, Pnt)` (`cxx:312-326`).
    pub fn tangency1(&self) -> (f64, f64, GpPnt2d) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2d2TanIter");
        (self.par1sol, self.pararg1, self.pnttg1sol)
    }

    /// `Tangency2(ParSol, ParArg, Pnt)` (`cxx:328-338`).
    pub fn tangency2(&self) -> (f64, f64, GpPnt2d) {
        assert!(self.well_done, "StdFail_NotDone in Geom2dGcc_Lin2d2TanIter");
        (self.par2sol, self.pararg2, self.pnttg2sol)
    }
}

/// `Geom2dAdaptor_Curve()` (`cxx:180`): the null adaptor the point overload of
/// `Add` passes as its second curve.
pub fn null_curve() -> Geom2dAdaptorCurve {
    Geom2dAdaptorCurve::default()
}
