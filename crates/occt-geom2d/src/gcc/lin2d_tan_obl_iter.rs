//! `Geom2dGcc_Lin2dTanOblIter`
//! (`Geom2dGcc_Lin2dTanOblIter.{hxx:34-77,cxx:38-195}`): one line tangent to
//! a curve at the angle `Angle` with a reference line, solved by a local root
//! search from the caller's parametric guess.
//!
//! | OCCT | `cxx` |
//! |---|---|
//! | `(QCurve, gp_Lin2d, Param1, TolAng, Angle)` | `:38-128` |
//! | `IsDone`..`Intersection2` | `:130-195` |
//!
//! Translation notes:
//! - `gp_Dir2d TheDirection(Dir)` starts as the reference line direction and is
//!   replaced by `Dir.Reversed()` when `| |Angle| - PI | <= |TolAng|`
//!   (`cxx:64-70`), by `gp_Dir2d(-B, A)` when `|Angle - PI/2| <= |TolAng|`
//!   (`cxx:71-74`), by `gp_Dir2d(B, -A)` when `|Angle + PI/2| <= |TolAng|`
//!   (`cxx:75-78`) and by the rotated `(A cos - B sin, A sin + B cos)` otherwise
//!   (`cxx:79-82`); `Paral2` is only raised by the reversed and the `|Angle| <=
//!   |TolAng|` (`cxx:83-87`) cases,
//! - `math_FunctionRoot sol(func, Param1, EpsX(Cu1, |TolAng|), U1, U2, 100)`
//!   (`cxx:88-92`) is [`MathFunctionRoot::new_with_bounds`]; the root search
//!   itself does not filter by the function value, so every accepted parameter
//!   only has to satisfy the sign conditions of `cxx:104-108`,
//! - `sign1 = Vect1 . TheDirection` and `sign2 = Vect2 x TheDirection`
//!   (`cxx:101-102`) are evaluated before the qualifier test,
//! - `par2sol` / `pararg2` are the projections of the intersection point on the
//!   respective line directions (`cxx:121-122`), computed only when the
//!   analytic intersection reports exactly one point,
//! - `Intersection2` throws `Geom2dGcc_IsParallel` when `Paral2` is set
//!   (`cxx:184-187`); the port panics with the same name,
//! - `IsParallel2` (`cxx:157-160`) is the only accessor that does not check
//!   `WellDone`,
//! - `StdFail_NotDone` is a panic.
//!
//! UNPORTED: none - the constructor and all seven accessors.

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intana2d::ana_intersection::IntAna2dAnaIntersection;
use occt_math::MathFunctionRoot;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::function_tan_obl::Geom2dGccFunctionTanObl;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc::qualified_curve::Geom2dGccQCurve;

/// `100` of `math_FunctionRoot` (`cxx:92`).
const NB_ITERATIONS: i32 = 100;

/// `Geom2dGcc_Lin2dTanOblIter` (`hxx:34-77`).
#[derive(Clone)]
pub struct Geom2dGccLin2dTanOblIter {
    /// `WellDone` (`hxx:65`).
    well_done: bool,
    /// `Paral2` (`hxx:66`).
    paral2: bool,
    /// `linsol` (`hxx:67`).
    linsol: GpLin2d,
    /// `qualifier1` (`hxx:68`).
    qualifier1: GccEntPosition,
    /// `pnttg1sol` (`hxx:69`).
    pnttg1sol: GpPnt2d,
    /// `pntint2sol` (`hxx:70`).
    pntint2sol: GpPnt2d,
    /// `par1sol` (`hxx:71`).
    par1sol: f64,
    /// `par2sol` (`hxx:72`).
    par2sol: f64,
    /// `pararg1` (`hxx:73`).
    pararg1: f64,
    /// `pararg2` (`hxx:74`).
    pararg2: f64,
}

impl Geom2dGccLin2dTanOblIter {
    /// `(QCurve, gp_Lin2d, Param1, TolAng, Angle)` (`cxx:38-128`).
    pub fn new(
        qualified1: &Geom2dGccQCurve,
        the_lin: &GpLin2d,
        param1: f64,
        tol_ang: f64,
        angle: f64,
    ) -> Self {
        // `: par2sol(0.0), pararg2(0.0)` (`cxx:41-42`), `par1sol = 0.;
        // pararg1 = 0.; WellDone = false;` (`cxx:47-49`); `Paral2` is set at
        // `cxx:53`.
        let mut s = Self {
            well_done: false,
            paral2: false,
            linsol: GpLin2d::default(),
            qualifier1: GccEntPosition::NoQualifier,
            pnttg1sol: GpPnt2d::zero(),
            pntint2sol: GpPnt2d::zero(),
            par1sol: 0.0,
            par2sol: 0.0,
            pararg1: 0.0,
            pararg2: 0.0,
        };
        // `cxx:49-53`.
        if !(qualified1.is_enclosed()
            || qualified1.is_enclosing()
            || qualified1.is_outside()
            || qualified1.is_unqualified())
        {
            panic!("GccEnt_BadQualifier in Geom2dGcc_Lin2dTanOblIter");
        }
        s.paral2 = false;
        let cu1: Geom2dAdaptorCurve = qualified1.qualified().clone();
        let u1 = curve_tool::first_parameter(&cu1);
        let u2 = curve_tool::last_parameter(&cu1);
        let dir = *the_lin.direction();
        let a = dir.x();
        let b = dir.y();
        let mut the_direction = dir;
        if angle.abs() > tol_ang.abs() {
            if (angle.abs() - std::f64::consts::PI).abs() <= tol_ang.abs() {
                // `cxx:68-70`.
                s.paral2 = true;
                the_direction = dir.reversed();
            } else if (angle - std::f64::consts::FRAC_PI_2).abs() <= tol_ang.abs() {
                // `cxx:71-74`.
                the_direction = GpDir2d::new(-b, a).expect("gp_Dir2d");
            } else if (angle + std::f64::consts::FRAC_PI_2).abs() <= tol_ang.abs() {
                // `cxx:75-78`.
                the_direction = GpDir2d::new(b, -a).expect("gp_Dir2d");
            } else {
                // `cxx:79-82`.
                the_direction = GpDir2d::new(
                    a * angle.cos() - b * angle.sin(),
                    a * angle.sin() + b * angle.cos(),
                )
                .expect("gp_Dir2d");
            }
        } else {
            // `cxx:83-87`.
            s.paral2 = true;
        }
        let mut func = Geom2dGccFunctionTanObl::new(&cu1, &the_direction);
        let sol = MathFunctionRoot::new_with_bounds(
            &mut func,
            param1,
            curve_tool::eps_x(&cu1, tol_ang.abs()),
            u1,
            u2,
            NB_ITERATIONS,
        );
        if sol.is_done() {
            let usol = sol.root();
            let (origine, vect1, vect2) =
                curve_tool::d2(&cu1, usol).expect("Geom2dAdaptor_Curve::D2");
            // `cxx:101-102`.
            let sign1 = vect1.xy().dot(&dir_xy(&the_direction));
            let sign2 = vect2.xy().crossed(&dir_xy(&the_direction));
            // `cxx:104-108`.
            if qualified1.is_unqualified()
                || (qualified1.is_enclosing() && sign2 <= 0.0)
                || (qualified1.is_outside() && sign1 <= 0.0 && sign2 >= 0.0)
                || (qualified1.is_enclosed() && sign1 >= 0.0 && sign2 >= 0.0)
            {
                s.well_done = true;
                s.linsol = GpLin2d::from_pnt_dir(origine, the_direction);
                s.pnttg1sol = origine;
                s.qualifier1 = qualified1.qualifier();
                s.pararg1 = usol;
                s.par1sol = 0.0;
                if !s.paral2 {
                    // `cxx:115-125`.
                    let intp = intersect_lin_lin(&s.linsol, the_lin);
                    if intp.is_done() && !intp.is_empty() && intp.nb_points() == 1 {
                        s.pntint2sol = *intp.point(1).value();
                        s.par2sol = projection(&s.linsol, &s.pntint2sol);
                        s.pararg2 = projection(the_lin, &s.pntint2sol);
                    }
                }
            }
        }
        s
    }

    /// `IsDone()` (`cxx:130-133`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `ThisSolution()` (`cxx:135-142`).
    pub fn this_solution(&self) -> GpLin2d {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Lin2dTanOblIter"
        );
        self.linsol
    }

    /// `WhichQualifier(Qualif1)` (`cxx:144-155`).
    pub fn which_qualifier(&self) -> GccEntPosition {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Lin2dTanOblIter"
        );
        self.qualifier1
    }

    /// `IsParallel2()` (`cxx:157-160`): no `WellDone` check.
    pub fn is_parallel2(&self) -> bool {
        self.paral2
    }

    /// `Tangency1(ParSol, ParArg, PntSol)` (`cxx:162-174`).
    pub fn tangency1(&self) -> (f64, f64, GpPnt2d) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Lin2dTanOblIter"
        );
        (self.par1sol, self.pararg1, self.pnttg1sol)
    }

    /// `Intersection2(ParSol, ParArg, PntSol)` (`cxx:176-195`).
    pub fn intersection2(&self) -> (f64, f64, GpPnt2d) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Lin2dTanOblIter"
        );
        assert!(
            !self.paral2,
            "Geom2dGcc_IsParallel in Geom2dGcc_Lin2dTanOblIter"
        );
        (self.par2sol, self.pararg2, self.pntint2sol)
    }
}

/// `gp_Dir2d::XY()` (`gp_Dir2d.lxx:24`).
fn dir_xy(dir: &GpDir2d) -> occt_core::gp::GpXY {
    occt_core::gp::GpXY::new(dir.x(), dir.y())
}

/// `IntAna2d_AnaIntersection(l1, l2)` (`IntAna2d_AnaIntersection_1.cxx`).
fn intersect_lin_lin(l1: &GpLin2d, l2: &GpLin2d) -> IntAna2dAnaIntersection {
    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(l1, l2);
    intp
}

/// `gp_Vec2d(L.Direction()).Dot(gp_Vec2d(L.Location(), P))` (`cxx:121-122`):
/// the signed distance of `P` along the line direction, measured from the line
/// location.
fn projection(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let loc = l.location();
    let dir = l.direction();
    dir.x() * (p.x() - loc.x()) + dir.y() * (p.y() - loc.y())
}
