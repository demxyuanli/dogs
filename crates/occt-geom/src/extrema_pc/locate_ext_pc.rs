//! Port of `Extrema_GLocateExtPC` (`Extrema_GLocateExtPC.hxx:40-399`)
//! monomorphised onto the 3D instantiation (`Extrema_LocateExtPC.hxx:29-36`):
//! `TheCurve=Adaptor3d_Curve`, `TheCurveTool=Extrema_CurveTool`,
//! `ThePoint=gp_Pnt`, `TheVector=gp_Vec`, `ThePOnC=Extrema_POnCurv`,
//! `TheELPC=Extrema_ExtPC`, `TheLocEPC=Extrema_GenLocateExtPC`.
//!
//! It is the class OCCT consumers build when they need "the extremum of the
//! distance from a point to a curve near a known parameter": `BRepLib.cxx:2343`
//! (`ContinuityOfFaces`), `IntTools_Context.cxx:829,890` (`IsVertexOnLine`),
//! `BRepLib_FuseEdges.cxx:185`, `GeomLib.cxx:511`.
//!
//! `Initialize` (`hxx:110-122`) switches on `TheCurveTool::GetType` between the
//! local root finder (`Extrema_GenLocateExtPC`, for Bezier / BSpline / Offset /
//! Other) and the global analytic solver (`Extrema_ExtPC`, for the elementary
//! types); `Perform` (`hxx:125-306`) walks the `GeomAbs_C2` intervals for the
//! free-form arm. Both dependencies are already ported:
//! [`super::gen_locate_ext_pc::GenLocateExtPC`] and
//! [`super::general_extrema_pc::extrema_ext_pc_range`] (= `Extrema_ExtPC`).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, Precision};

use crate::curve::Curve;
use crate::extrema_cc::POnCurv;

use super::curve_tool::{self, GeomAbsCurveType};
use super::gen_locate_ext_pc::GenLocateExtPC;
use super::general_extrema_pc::{extrema_ext_pc_range_done, ExtPcSolution};

/// `GeomAbs_C2` (`GeomAbs_Shape.hxx:49-55`: C0, G1, C1, G2, C2, C3, CN).
const GEOMABS_C2: u8 = 4;

/// `Extrema_LocateExtPC` (`Extrema_LocateExtPC.hxx:29-36`), 3D instantiation.
pub struct LocateExtPC<'a> {
    my_c: Option<&'a dyn Curve>,
    my_dist2: f64,
    my_ismin: bool,
    my_done: bool,
    my_umin: f64,
    my_usup: f64,
    mytol: f64,
    type_: GeomAbsCurveType,
    numberext: usize,
    mypp: POnCurv,
    my_loc_ext_pc: GenLocateExtPC<'a>,
    my_ext_pc_sols: Vec<ExtPcSolution>,
}

impl<'a> Default for LocateExtPC<'a> {
    /// `Extrema_GLocateExtPC()` (`hxx:55-66`).
    fn default() -> Self {
        Self {
            my_c: None,
            my_dist2: 0.0,
            my_ismin: false,
            my_done: false,
            my_umin: 0.0,
            my_usup: 0.0,
            mytol: 0.0,
            type_: GeomAbsCurveType::OtherCurve,
            numberext: 0,
            mypp: POnCurv::new(0.0, GpPnt::zero()),
            my_loc_ext_pc: GenLocateExtPC::default(),
            my_ext_pc_sols: Vec::new(),
        }
    }
}

impl<'a> LocateExtPC<'a> {
    /// `Extrema_GLocateExtPC(P, C, U0, TolF)` (`hxx:73-82`): the search window
    /// is the curve's own `[FirstParameter, LastParameter]`.
    pub fn new(p: &GpPnt, c: &'a dyn Curve, u0: f64, tol_f: f64) -> Self {
        let mut r = Self::default();
        r.initialize(
            c,
            curve_tool::extrema_first_parameter(c),
            curve_tool::extrema_last_parameter(c),
            tol_f,
        );
        r.perform(p, u0);
        r
    }

    /// `Extrema_GLocateExtPC(P, C, U0, Umin, Usup, TolF)` (`hxx:90-97`).
    pub fn new_range(
        p: &GpPnt,
        c: &'a dyn Curve,
        u0: f64,
        u_min: f64,
        u_sup: f64,
        tol_f: f64,
    ) -> Self {
        let mut r = Self::default();
        r.initialize(c, u_min, u_sup, tol_f);
        r.perform(p, u0);
        r
    }

    /// `Initialize(theC, Umin, Usup, TolF)` (`hxx:110-122`).
    pub fn initialize(&mut self, c: &'a dyn Curve, u_min: f64, u_sup: f64, tol_f: f64) {
        self.my_c = Some(c);
        self.mytol = tol_f;
        self.my_umin = u_min;
        self.my_usup = u_sup;
        self.type_ = curve_tool::extrema_curve_type(c);
        let tolu = curve_tool::resolution(c, CONFUSION);
        match self.type_ {
            GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::BezierCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::OtherCurve => {
                self.my_loc_ext_pc.initialize(c, u_min, u_sup, tolu);
            }
            _ => {
                // `myExtremPC.Initialize(theC, theUmin, theUsup, tolu)`
                // (`hxx:120`). The window is stored and replayed in `Perform`,
                // where `Extrema_ExtPC::Perform(theP)` runs the analytic
                // solver over exactly those bounds.
                self.my_umin = u_min;
                self.my_usup = u_sup;
            }
        }
    }

    /// `Perform(theP, U0)` (`hxx:125-306`).
    pub fn perform(&mut self, the_p: &GpPnt, the_u0: f64) {
        let c = match self.my_c {
            Some(c) => c,
            None => return,
        };
        match self.type_ {
            GeomAbsCurveType::OtherCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::BSplineCurve => self.perform_freeform(c, the_p, the_u0),
            GeomAbsCurveType::BezierCurve => {
                self.my_loc_ext_pc.perform(the_p, the_u0);
                self.my_done = self.my_loc_ext_pc.is_done();
            }
            _ => self.perform_elementary(c, the_p, the_u0),
        }
    }

    /// `hxx:132-272`: the free-form arm, walking the `GeomAbs_C2` intervals.
    fn perform_freeform(&mut self, c: &'a dyn Curve, the_p: &GpPnt, the_u0: f64) {
        let n = curve_tool::extrema_nb_intervals(c, GEOMABS_C2);
        let the_inter = curve_tool::extrema_intervals(c, GEOMABS_C2);
        // Be gentle with the caller (`hxx:142-150`).
        let mut local_u0 = the_u0;
        if local_u0 < self.my_umin {
            local_u0 = self.my_umin;
        } else if local_u0 > self.my_usup {
            local_u0 = self.my_usup;
        }
        let mut myintuinf = 0.0;
        let mut myintusup = 0.0;
        // Search for the interval containing U0 (`hxx:152-162`).
        let mut found = false;
        let mut inter: i32 = 1;
        while !found && inter <= n {
            myintuinf = the_inter[(inter - 1) as usize].max(self.my_umin);
            myintusup = the_inter[inter as usize].min(self.my_usup);
            if local_u0 >= myintuinf && local_u0 < myintusup {
                found = true;
            }
            inter += 1;
        }
        if found {
            // `hxx:164-165`: inter is increased after found.
            inter -= 1;
        }

        self.my_loc_ext_pc
            .initialize(c, myintuinf, myintusup, self.mytol);
        self.my_loc_ext_pc.perform(the_p, local_u0);
        self.my_done = self.my_loc_ext_pc.is_done();
        if self.my_done {
            self.mypp = self.my_loc_ext_pc.point(1);
            self.my_ismin = self.my_loc_ext_pc.is_min();
            self.my_dist2 = self.my_loc_ext_pc.square_distance();
            return;
        }

        // Try on the neighbouring intervals (`hxx:178-270`).
        let mut k: i32 = 1;
        let mut i1 = inter;
        let mut i2 = inter;
        // `s2inf` / `s1sup` are the dot products at the two ends of the
        // interval that was tried first (`hxx:185-188`).
        let (p1, v1) = curve_tool::extrema_d1(c, myintuinf);
        let mut s2inf = v1.dot(&GpVec::from_pnts(the_p, &p1));
        let (p1, v1) = curve_tool::extrema_d1(c, myintusup);
        let mut s1sup = v1.dot(&GpVec::from_pnts(the_p, &p1));

        while !self.my_done && i2 > 0 && i1 <= n {
            i1 = inter + k;
            i2 = inter - k;
            if i1 <= n {
                myintuinf = the_inter[(i1 - 1) as usize].max(self.my_umin);
                myintusup = the_inter[i1 as usize].min(self.my_usup);
                if myintuinf < myintusup {
                    let (p1, v1) = curve_tool::extrema_d1(c, myintuinf);
                    let s2sup = v1.dot(&GpVec::from_pnts(the_p, &p1));
                    if Precision::is_infinite(s2sup) || Precision::is_infinite(s1sup) {
                        break;
                    }
                    if s1sup * s2sup <= f64::EPSILON {
                        // Extremum (`hxx:208-215`).
                        self.my_done = true;
                        self.mypp = POnCurv::new(myintuinf, p1);
                        self.my_ismin = s1sup <= 0.0;
                        self.my_dist2 = the_p.square_distance(&p1);
                        break;
                    }
                    let (p1, v1) = curve_tool::extrema_d1(c, myintusup);
                    s1sup = v1.dot(&GpVec::from_pnts(the_p, &p1));
                    self.my_loc_ext_pc
                        .initialize(c, myintuinf, myintusup, self.mytol);
                    self.my_loc_ext_pc
                        .perform(the_p, (myintuinf + myintusup) * 0.5);
                    self.my_done = self.my_loc_ext_pc.is_done();
                    if self.my_done {
                        self.mypp = self.my_loc_ext_pc.point(1);
                        self.my_ismin = self.my_loc_ext_pc.is_min();
                        self.my_dist2 = self.my_loc_ext_pc.square_distance();
                        break;
                    }
                }
            }

            if i2 > 0 {
                myintuinf = the_inter[(i2 - 1) as usize].max(self.my_umin);
                myintusup = the_inter[i2 as usize].min(self.my_usup);
                if myintuinf < myintusup {
                    let (p1, v1) = curve_tool::extrema_d1(c, myintusup);
                    let s1inf = v1.dot(&GpVec::from_pnts(the_p, &p1));
                    if Precision::is_infinite(s2inf) || Precision::is_infinite(s1inf) {
                        break;
                    }
                    if s1inf * s2inf <= f64::EPSILON {
                        // Extremum (`hxx:242-249`).
                        self.my_done = true;
                        self.mypp = POnCurv::new(myintusup, p1);
                        self.my_ismin = s1inf <= 0.0;
                        self.my_dist2 = the_p.square_distance(&p1);
                        break;
                    }
                    let (p1, v1) = curve_tool::extrema_d1(c, myintuinf);
                    s2inf = v1.dot(&GpVec::from_pnts(the_p, &p1));
                    self.my_loc_ext_pc
                        .initialize(c, myintuinf, myintusup, self.mytol);
                    self.my_loc_ext_pc
                        .perform(the_p, (myintuinf + myintusup) * 0.5);
                    self.my_done = self.my_loc_ext_pc.is_done();
                    if self.my_done {
                        self.mypp = self.my_loc_ext_pc.point(1);
                        self.my_ismin = self.my_loc_ext_pc.is_min();
                        self.my_dist2 = self.my_loc_ext_pc.square_distance();
                        break;
                    }
                }
            }

            k += 1;
        }
    }

    /// `hxx:281-304`: the elementary arm runs the whole `Extrema_ExtPC` and
    /// keeps the solution whose parameter is closest to U0.
    fn perform_elementary(&mut self, c: &dyn Curve, the_p: &GpPnt, the_u0: f64) {
        let (sols, done) = extrema_ext_pc_range_done(c, the_p, self.my_umin, self.my_usup);
        self.my_ext_pc_sols = sols;
        self.numberext = 0;
        if done {
            // `hxx:285-297`: `for (i = 1; i <= NbExt; i++)`, `valU2` starts at
            // `RealLast()`.
            let mut valu2 = f64::MAX;
            for (idx, s) in self.my_ext_pc_sols.iter().enumerate() {
                let valu = (s.u - the_u0).abs();
                if valu <= valu2 {
                    valu2 = valu;
                    self.numberext = idx + 1;
                    self.my_done = true;
                }
            }
        }
        if self.numberext == 0 {
            self.my_done = false;
        }
    }

    /// `IsDone` (`hxx:309`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `SquareDistance` (`hxx:312-336`).
    ///
    /// OCCT raises `StdFail_NotDone` when `IsDone()` is false; the port leaves
    /// `IsDone` as the caller's precondition (as everywhere else in this
    /// crate) and otherwise follows the same body, including the
    /// `if (numberext != 0)` zero fallback of the elementary arm.
    pub fn square_distance(&self) -> f64 {
        match self.type_ {
            GeomAbsCurveType::BezierCurve => self.my_loc_ext_pc.square_distance(),
            GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::OtherCurve => self.my_dist2,
            _ => {
                if self.numberext != 0 {
                    self.my_ext_pc_sols[self.numberext - 1].sq_dist
                } else {
                    0.0
                }
            }
        }
    }

    /// `IsMin` (`hxx:339-362`); see [`Self::square_distance`] for the
    /// `NotDone` precondition.
    pub fn is_min(&self) -> bool {
        match self.type_ {
            GeomAbsCurveType::BezierCurve => self.my_loc_ext_pc.is_min(),
            GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::OtherCurve => self.my_ismin,
            _ => {
                if self.numberext != 0 {
                    self.my_ext_pc_sols[self.numberext - 1].is_min
                } else {
                    false
                }
            }
        }
    }

    /// `Point` (`hxx:365-382`); see [`Self::square_distance`] for the
    /// `NotDone` precondition.
    pub fn point(&self) -> POnCurv {
        match self.type_ {
            GeomAbsCurveType::BezierCurve => self.my_loc_ext_pc.point(1),
            GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::OffsetCurve
            | GeomAbsCurveType::OtherCurve => self.mypp,
            _ => {
                if self.numberext != 0 {
                    let s = &self.my_ext_pc_sols[self.numberext - 1];
                    POnCurv::new(s.u, s.point)
                } else {
                    POnCurv::default()
                }
            }
        }
    }

    /// `GetType` of the curve this instance was initialized with (`hxx:110`).
    pub fn curve_type(&self) -> GeomAbsCurveType {
        self.type_
    }
}
