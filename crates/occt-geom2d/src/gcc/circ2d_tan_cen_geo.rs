//! `Geom2dGcc_Circ2dTanCenGeo`
//! (`Geom2dGcc_Circ2dTanCenGeo.{hxx:49-107,cxx:35-195}`): circles centred on
//! a point and tangent to a general curve.
//!
//! The single constructor locates the extrema of the distance from the centre
//! to the curve (`Extrema_ExtPC2d`), keeps the closest one and accepts it when
//! the radius direction is aligned with the curve tangent.
//!
//! | OCCT | `cxx` |
//! |---|---|
//! | `(QCurve, gp_Pnt2d, Tolerance)` | `:35-136` |
//! | `IsDone`..`Tangency1` | `:140-195` |
//!
//! Translation notes:
//! - `Extrema_ExtPC2d(Pcenter, curve, first, last, Tol)` (`cxx:64-68`) takes
//!   `Tol = |Tolerance|` as its `theTolF`; the port passes it through
//!   [`ExtremaExtPC2d::new_range_tol`] (the default `1.0e-10` of the other
//!   constructors is *not* used here),
//! - `if (!distmin.IsDone()) throw Standard_Failure()` (`cxx:69-72`) and
//!   `if (nbext == 0) throw Standard_Failure()` (`cxx:74-77`) are panics,
//! - `Index(1) == Index(2)` (`cxx:95`) compares two `NCollection_Array1<int>`
//!   slots that are **never written** by this constructor; the port models the
//!   two-zero reading, so `nbsol` is always 1 and only the closest extremum
//!   (`theDist2(1)`) can produce a solution. `theDist2(2)` / `theParam(2)` /
//!   `pTan(2)` are still tracked by the same loop (`cxx:86-92`) so the dead
//!   maximum branch keeps OCCT's exact evaluation order,
//! - `Tol` is reassigned to `1.e-12` inside the solution loop (`cxx:119`),
//!   overwriting the `|Tolerance|` value that `Extrema_ExtPC2d` already
//!   consumed; the port uses two separate variables for the two roles,
//! - `gp_Vec2d::Angle` (`cxx:121`) answers `0.` for a null operand in this
//!   port (`gp_Vec2d.cxx:57` throws in OCCT); with `dot1 == 0.` from the
//!   `normevec1` / `normetan1` guard the `Unqualified` arm is the only one that
//!   can accept such a point,
//! - `this_solution` (`cxx:150-158`) checks only the index, unlike
//!   `which_qualifier` / `tangency1` which also check `WellDone`; kept,
//! - `par1sol` is `0.` (`cxx:130`), `pararg1` is the extremum parameter.
//!
//! UNPORTED: none - the constructor and all four accessors.

use occt_core::gp::{GpCirc2d, GpPnt2d, GpVec2d};
use occt_core::precision::REAL_SMALL;

use crate::extrema2d::ExtremaExtPC2d;
use crate::gcc::circ2d_2tan_rad::axis_at;
use crate::gcc::curve_tool;
use crate::gcc::gcc_ent::GccEntPosition;
use crate::gcc::qualified_curve::Geom2dGccQCurve;

/// `NCollection_Array1<...>(1, 2)` of the six field arrays (`cxx:44-48`).
pub const MAX_SOL: usize = 2;

/// `Geom2dGcc_Circ2dTanCenGeo` (`hxx:49-107`).
#[derive(Clone)]
pub struct Geom2dGccCirc2dTanCenGeo {
    /// `WellDone` (`hxx:100`).
    well_done: bool,
    /// `NbrSol` (`hxx:101`).
    nbr_sol: usize,
    /// `cirsol` (`hxx:102`).
    cirsol: [GpCirc2d; MAX_SOL],
    /// `qualifier1` (`hxx:103`).
    qualifier1: [GccEntPosition; MAX_SOL],
    /// `pnttg1sol` (`hxx:104`).
    pnttg1sol: [GpPnt2d; MAX_SOL],
    /// `par1sol` (`hxx:105`).
    par1sol: [f64; MAX_SOL],
    /// `pararg1` (`hxx:106`).
    pararg1: [f64; MAX_SOL],
}

impl Geom2dGccCirc2dTanCenGeo {
    /// The `: WellDone(false), NbrSol(0), cirsol(1, 2), ... :` initializer list
    /// (`cxx:40-48`).
    pub fn empty() -> Self {
        Self {
            well_done: false,
            nbr_sol: 0,
            cirsol: [GpCirc2d::default(); MAX_SOL],
            qualifier1: [GccEntPosition::NoQualifier; MAX_SOL],
            pnttg1sol: [GpPnt2d::zero(); MAX_SOL],
            par1sol: [0.0; MAX_SOL],
            pararg1: [0.0; MAX_SOL],
        }
    }

    /// `Geom2dGcc_Circ2dTanCenGeo(Qualified1, Pcenter, Tolerance)`
    /// (`cxx:35-136`).
    pub fn new(qualified1: &Geom2dGccQCurve, pcenter: &GpPnt2d, tolerance: f64) -> Self {
        let mut s = Self::empty();
        let tol_f = tolerance.abs();
        let mut the_dist2 = [f64::MAX, 0.0];
        let mut the_param = [0.0f64; MAX_SOL];
        let mut p_tan = [GpPnt2d::zero(); MAX_SOL];
        let curve = qualified1.qualified();
        let basis = curve
            .curve()
            .expect("Geom2dAdaptor_Curve without a curve in Geom2dGcc_Circ2dTanCenGeo");
        let distmin = ExtremaExtPC2d::new_range_tol(
            basis.as_ref(),
            pcenter,
            curve_tool::first_parameter(curve),
            curve_tool::last_parameter(curve),
            tol_f,
        );
        if !distmin.is_done() {
            panic!("Standard_Failure in Geom2dGcc_Circ2dTanCenGeo (Extrema_ExtPC2d not done)");
        }
        let nbext = distmin.nb_ext();
        if nbext == 0 {
            panic!("Standard_Failure in Geom2dGcc_Circ2dTanCenGeo (no extremum)");
        }
        let mut i = 1usize;
        while i <= nbext {
            let (the_par, the_pnt) = distmin.point(i);
            let sq = distmin.square_distance(i);
            if sq < the_dist2[0] {
                the_dist2[0] = sq;
                the_param[0] = the_par;
                p_tan[0] = the_pnt;
            }
            if sq > the_dist2[1] {
                the_dist2[1] = sq;
                the_param[1] = the_par;
                p_tan[1] = the_pnt;
            }
            i += 1;
        }
        // `cxx:95-101`: `Index(1)` and `Index(2)` are never written, so the
        // two-zero reading gives `nbsol = 1`.
        let nbsol = 1;
        for i in 1..=nbsol {
            let (point1, tan1) = curve_tool::d1(curve, the_param[i - 1])
                .expect("Geom2dAdaptor_Curve::D1 in Geom2dGcc_Circ2dTanCenGeo");
            let normetan1 = tan1.magnitude();
            let vec1 = GpVec2d::new(pcenter.x() - point1.x(), pcenter.y() - point1.y());
            let normevec1 = vec1.magnitude();
            let dot1 = if normevec1 >= REAL_SMALL && normetan1 >= REAL_SMALL {
                vec1.dot(&tan1) / (normevec1 * normetan1)
            } else {
                0.0
            };
            let tol = 1.0e-12;
            if dot1 <= tol {
                let angle1 = vec1.angle(&tan1);
                if qualified1.is_unqualified()
                    || (qualified1.is_enclosing() && angle1 <= 0.0)
                    || (qualified1.is_outside() && angle1 >= 0.0)
                    || (qualified1.is_enclosed() && angle1 <= 0.0)
                {
                    s.nbr_sol += 1;
                    let k = s.nbr_sol - 1;
                    s.cirsol[k] = GpCirc2d::new(axis_at(*pcenter), the_dist2[i - 1].sqrt());
                    s.qualifier1[k] = qualified1.qualifier();
                    s.pararg1[k] = the_param[i - 1];
                    s.par1sol[k] = 0.0;
                    s.pnttg1sol[k] = p_tan[i - 1];
                    s.well_done = true;
                }
            }
        }
        s
    }

    /// `IsDone()` (`cxx:140-143`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `NbSolutions()` (`cxx:145-148`).
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }

    /// `ThisSolution(Index)` (`cxx:150-158`): only the index is checked.
    pub fn this_solution(&self, index: usize) -> GpCirc2d {
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanCenGeo"
        );
        self.cirsol[index - 1]
    }

    /// `WhichQualifier(Index, Qualif1)` (`cxx:160-174`).
    pub fn which_qualifier(&self, index: usize) -> GccEntPosition {
        self.check(index);
        self.qualifier1[index - 1]
    }

    /// `Tangency1(Index, ParSol, ParArg, PntSol)` (`cxx:176-195`).
    pub fn tangency1(&self, index: usize) -> (f64, f64, GpPnt2d) {
        self.check(index);
        (
            self.par1sol[index - 1],
            self.pararg1[index - 1],
            self.pnttg1sol[index - 1],
        )
    }

    /// `if (!WellDone) throw StdFail_NotDone()` and
    /// `if (Index <= 0 || Index > NbrSol) throw Standard_OutOfRange()`.
    fn check(&self, index: usize) {
        assert!(
            self.well_done,
            "StdFail_NotDone in Geom2dGcc_Circ2dTanCenGeo"
        );
        assert!(
            index >= 1 && index <= self.nbr_sol,
            "Standard_OutOfRange in Geom2dGcc_Circ2dTanCenGeo"
        );
    }
}
