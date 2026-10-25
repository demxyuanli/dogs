//! `Extrema_ExtElC` — distances between two elementary curves.
//!
//! Source: `ModelingData/TKGeomBase/Extrema/Extrema_ExtElC.hxx:36-92`,
//! `.cxx:256-1168` (OCCT 8.0.0).
//!
//! Only the `gp_Lin` / `gp_Lin` constructor (`.cxx:268-358`) is ported; the
//! five other constructors (line/circle `.cxx:471-625`, line/ellipse
//! `.cxx:627-754`, line/hyperbola `.cxx:757-859`, line/parabola `.cxx:862-952`,
//! circle/circle `.cxx:955-1114`) are **UNPORTED**. The one ported arm is the
//! only one OCCT needs to build the circle through three points
//! (`gce_MakeCirc.cxx:94`, see [`crate::gce_make_circ`]).
//!
//! `Extrema_POnCurv` (`Extrema_POnCurv.hxx:24-60`) is defined here because this
//! is its only consumer: a parameter on a curve plus the point at that
//! parameter.

use occt_core::elib::clib;
use occt_core::gp::{GpLin, GpPnt};
use occt_core::precision::{Precision, ANGULAR, REAL_SMALL};

/// `RealLast()` (`Standard_Real.hxx:110-113`) — the initial `mySqDist` value
/// (`cxx:261-264`).
const REAL_LAST: f64 = f64::MAX;

/// `Extrema_POnCurv` (`Extrema_POnCurv.hxx:24-60`): `myU` plus `myP`, plus the
/// 1-based readers `Value()` / `Parameter()`.
#[derive(Debug, Clone, Copy, Default)]
pub struct ExtremaPOnCurv {
    /// `myU` (`hxx:58`).
    pub u: f64,
    /// `myP` (`hxx:59`).
    pub p: GpPnt,
}

impl ExtremaPOnCurv {
    /// `Extrema_POnCurv(theU, theP)` (`hxx:38-42`).
    pub fn new(u: f64, p: GpPnt) -> Self {
        Self { u, p }
    }

    /// `Parameter()` (`hxx:56-57`).
    pub fn parameter(&self) -> f64 {
        self.u
    }

    /// `Value()` (`hxx:52-53`).
    pub fn value(&self) -> &GpPnt {
        &self.p
    }

    /// `SetValues(theU, theP)` (`hxx:44-48`).
    pub fn set_values(&mut self, u: f64, p: GpPnt) {
        self.u = u;
        self.p = p;
    }
}

/// `class Extrema_ExtElC` (`hxx:36-92`): the six `myPoint` slots of
/// `hxx:92` are the largest array OCCT declares for any of its constructors;
/// this port only fills slot `0`.
#[derive(Debug, Clone)]
pub struct ExtremaExtElC {
    /// `myDone` (`hxx:88`).
    done: bool,
    /// `myIsPar` (`hxx:89`).
    is_par: bool,
    /// `myNbExt` (`hxx:90`).
    nb_ext: usize,
    /// `mySqDist[6]` (`hxx:91`).
    sq_dist: [f64; 6],
    /// `myPoint[6][2]` (`hxx:92`).
    point: [[ExtremaPOnCurv; 2]; 6],
}

impl ExtremaExtElC {
    /// `Extrema_ExtElC()` (`cxx:256-265`): `mySqDist` starts at `RealLast()`.
    pub fn new() -> Self {
        Self {
            done: false,
            is_par: false,
            nb_ext: 0,
            sq_dist: [REAL_LAST; 6],
            point: [[ExtremaPOnCurv::default(); 2]; 6],
        }
    }

    /// `Extrema_ExtElC(const gp_Lin& C1, const gp_Lin& C2, const double)`
    /// (`cxx:268-358`). The tolerance argument is unnamed and unused in OCCT
    /// (the parallel test uses `gp::Resolution()` and `Precision::Angular()`),
    /// so the port takes no tolerance.
    pub fn lines(c1: &GpLin, c2: &GpLin) -> Self {
        let mut this = Self::new();

        // `const gp_Dir &aD1 = theC1.Position().Direction(), &aD2 = ...;
        // const double aCosA = aD1.Dot(aD2); const double aSqSinA = 1.0 -
        // aCosA * aCosA;` (`cxx:324-327`).
        let d1 = c1.direction();
        let d2 = c2.direction();
        let cos_a = d1.dot(&d2);
        let sq_sin_a = 1.0 - cos_a * cos_a;
        let mut u1 = 0.0f64;
        let mut u2 = 0.0f64;

        // `if (aSqSinA < gp::Resolution() || aD1.IsParallel(aD2,
        // Precision::Angular())) { myIsPar = true; } else { ... }`
        // (`cxx:328-339`).
        if sq_sin_a < REAL_SMALL || d1.is_parallel_tol(&d2, ANGULAR) {
            this.is_par = true;
        } else {
            let l = c2.location().xyz().subtracted(c1.location().xyz());
            let d1l = d1.xyz().dot(&l);
            let d2l = d2.xyz().dot(&l);
            u1 = (d1l - cos_a * d2l) / sq_sin_a;
            u2 = (cos_a * d1l - d2l) / sq_sin_a;
            this.is_par = Precision::is_infinite(u1) || Precision::is_infinite(u2);
        }

        // `if (myIsPar) { mySqDist[0] = theC2.SquareDistance(theC1.Location());
        // myNbExt = 1; myDone = true; return; }` (`cxx:341-348`).
        if this.is_par {
            this.sq_dist[0] = c2.square_distance(&c1.location());
            this.nb_ext = 1;
            this.done = true;
            return this;
        }

        // `const gp_Pnt aP1(ElCLib::Value(aU1, theC1)), aP2(ElCLib::Value(aU2,
        // theC2)); mySqDist[myNbExt] = aP1.SquareDistance(aP2);
        // myPoint[myNbExt][0] = Extrema_POnCurv(aU1, aP1); myPoint[myNbExt][1]
        // = Extrema_POnCurv(aU2, aP2); myNbExt = 1; myDone = true;`
        // (`cxx:352-357`).
        let p1 = clib::line_value(c1, u1);
        let p2 = clib::line_value(c2, u2);
        this.sq_dist[0] = p1.square_distance(&p2);
        this.point[0] = [ExtremaPOnCurv::new(u1, p1), ExtremaPOnCurv::new(u2, p2)];
        this.nb_ext = 1;
        this.done = true;
        this
    }

    /// `IsDone()` (`cxx:1117-1120`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`cxx:1124-1131`): OCCT throws `StdFail_NotDone`.
    pub fn is_parallel(&self) -> Result<bool, String> {
        if !self.is_done() {
            return Err("StdFail_NotDone: Extrema_ExtElC::IsParallel".into());
        }
        Ok(self.is_par)
    }

    /// `NbExt()` (`cxx:1135-1142`): OCCT throws `StdFail_NotDone`.
    pub fn nb_ext(&self) -> Result<usize, String> {
        if !self.is_done() {
            return Err("StdFail_NotDone: Extrema_ExtElC::NbExt".into());
        }
        Ok(self.nb_ext)
    }

    /// `SquareDistance(N)` (`cxx:1146-1154`), 1-based: OCCT throws
    /// `Standard_OutOfRange` outside `1..=NbExt()`.
    pub fn square_distance(&self, n: usize) -> Result<f64, String> {
        let nb = self.nb_ext()?;
        if n < 1 || n > nb {
            return Err("Standard_OutOfRange: Extrema_ExtElC::SquareDistance".into());
        }
        Ok(self.sq_dist[n - 1])
    }

    /// `Points(N, P1, P2)` (`cxx:1158-1166`), 1-based, same range check.
    pub fn points(&self, n: usize) -> Result<(ExtremaPOnCurv, ExtremaPOnCurv), String> {
        let nb = self.nb_ext()?;
        if n < 1 || n > nb {
            return Err("Standard_OutOfRange: Extrema_ExtElC::Points".into());
        }
        Ok((self.point[n - 1][0], self.point[n - 1][1]))
    }
}

impl Default for ExtremaExtElC {
    fn default() -> Self {
        Self::new()
    }
}
