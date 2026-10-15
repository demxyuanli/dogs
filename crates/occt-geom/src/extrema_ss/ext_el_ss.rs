//! `Extrema_ExtElSS` (TKGeomBase). Source: `Extrema_ExtElSS.cxx` / `.hxx`.
//!
//! Only the `(gp_Pln, gp_Pln)` overload is implemented in OCCT
//! (`Extrema_ExtElSS.cxx:40-54`); every other pair sets the flags and throws
//! `Standard_NotImplemented` (`cxx:62-69`, `cxx:77-83`, `cxx:91-97`,
//! `cxx:106-112`, `cxx:121-127`). For non-parallel planes `Perform` leaves
//! `myNbExt = 0`; for parallel planes it records the square distance and
//! **no points** (`myPOnS1` / `myPOnS2` are never allocated).

use occt_core::gp::GpPln;
use occt_core::precision::ANGULAR;

/// `Extrema_ExtElSS` (`Extrema_ExtElSS.hxx:30-100`).
#[derive(Debug, Clone)]
pub struct ExtremaExtElSS {
    /// `myDone`.
    done: bool,
    /// `myIsPar`.
    is_par: bool,
    /// `myNbExt`.
    nb_ext: i32,
    /// `mySqDist` — `Extrema_ExtElSS.cxx:50-52` allocates `(1, 1)` only on the
    /// parallel branch.
    sq_dist: Vec<f64>,
}

impl ExtremaExtElSS {
    /// `Extrema_ExtElSS()` (`Extrema_ExtElSS.cxx:28-33`).
    pub fn new() -> Self {
        ExtremaExtElSS {
            done: false,
            is_par: false,
            nb_ext: 0,
            sq_dist: Vec::new(),
        }
    }

    /// `Extrema_ExtElSS::Perform(const gp_Pln& S1, const gp_Pln& S2)`
    /// (`Extrema_ExtElSS.cxx:40-54`).
    pub fn new_plane_plane(s1: &GpPln, s2: &GpPln) -> Self {
        let mut me = ExtremaExtElSS::new();
        me.perform_plane_plane(s1, s2);
        me
    }

    /// `Extrema_ExtElSS::Perform(const gp_Pln&, const gp_Pln&)`
    /// (`Extrema_ExtElSS.cxx:40-54`).
    ///
    /// `myDone = true; myIsPar = false; myNbExt = 0;` then, when the two axes
    /// are parallel within `Precision::Angular()`, `myIsPar = true`,
    /// `myNbExt = 1` and `mySqDist(1) = S1.SquareDistance(S2)`.
    pub fn perform_plane_plane(&mut self, s1: &GpPln, s2: &GpPln) {
        self.done = true;
        self.is_par = false;
        self.nb_ext = 0;
        self.sq_dist.clear();
        if s1
            .axis()
            .direction()
            .is_parallel_tol(s2.axis().direction(), ANGULAR)
        {
            self.is_par = true;
            self.nb_ext = 1;
            self.sq_dist.push(s1.square_distance_pln(s2));
        }
    }

    /// `Extrema_ExtElSS::IsDone` (`Extrema_ExtElSS.cxx:130-133`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Extrema_ExtElSS::IsParallel` (`Extrema_ExtElSS.cxx:135-142`): throws
    /// `StdFail_NotDone` when `!IsDone()`.
    pub fn is_parallel(&self) -> bool {
        assert!(self.is_done(), "Extrema_ExtElSS::IsParallel: StdFail_NotDone");
        self.is_par
    }

    /// `Extrema_ExtElSS::NbExt` (`Extrema_ExtElSS.cxx:144-151`): throws
    /// `StdFail_NotDone` when `!IsDone()`.
    pub fn nb_ext(&self) -> i32 {
        assert!(self.is_done(), "Extrema_ExtElSS::NbExt: StdFail_NotDone");
        self.nb_ext
    }

    /// `Extrema_ExtElSS::SquareDistance(N)` (`Extrema_ExtElSS.cxx:154-161`):
    /// `Standard_OutOfRange` outside `1..=NbExt()`.
    pub fn square_distance(&self, n: i32) -> f64 {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_ExtElSS::SquareDistance: Standard_OutOfRange"
        );
        self.sq_dist[(n - 1) as usize]
    }
}

impl Default for ExtremaExtElSS {
    fn default() -> Self {
        ExtremaExtElSS::new()
    }
}
