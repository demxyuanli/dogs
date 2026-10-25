//! `TopTrans_CurveTransition` — transition of a curve crossing a curvilinear
//! boundary.
//!
//! Source: `TKGeomAlgo/TopTrans/TopTrans_CurveTransition.cxx` (OCCT 8.0.0,
//! 441 lines) with the declaration in `TopTrans_CurveTransition.hxx:47-110`.
//!
//! The boundary is a set of curve elements, each described locally at the
//! intersection point by its tangent, its normal and its curvature. The first
//! element that agrees with the swept curve becomes both the "first" and the
//! "last" element (`cxx:94-128`); every later element is compared against them
//! through the sign of `myTgt * T` and, on a tie, through
//! [`CurveTransition::is_before`].
//!
//! Kept verbatim from OCCT:
//! - `myNorm` is only written by `Reset(Tgt, Norm, Curv)` (`cxx:42-48`); the
//!   straight-line overload `Reset(Tgt)` (`cxx:56-61`) leaves it untouched, so
//!   a later `IsBefore` reads whatever the curve reset wrote.
//! - `IsBefore` uses `N1 * myNorm` / `N2 * myNorm` even when `myCurv == 0`
//!   (`cxx:369`, `:377`) — the guards `C1 == 0. || myCurv == 0.` are what
//!   decide, not `myCurv` alone.
//! - `StateBefore` / `StateAfter` return `OUT` for an unreachable
//!   orientation (`cxx:293`, `:317`).

use occt_core::gp::GpDir;

use crate::abs::{Orientation, State};

/// `GREATER (1)` (`cxx:20`).
const GREATER: i32 = 1;
/// `SAME (0)` (`cxx:21`).
const SAME: i32 = 0;
/// `LOWER (-1)` (`cxx:22`).
const LOWER: i32 = -1;

/// `TopTrans_CurveTransition` (`hxx:47-110`).
#[derive(Debug, Clone)]
pub struct CurveTransition {
    /// `myTgt`
    my_tgt: GpDir,
    /// `myNorm`
    my_norm: GpDir,
    /// `myCurv`
    my_curv: f64,
    /// `Init`
    init: bool,
    /// `TgtFirst`
    tgt_first: GpDir,
    /// `NormFirst`
    norm_first: GpDir,
    /// `CurvFirst`
    curv_first: f64,
    /// `TranFirst`
    tran_first: Orientation,
    /// `TgtLast`
    tgt_last: GpDir,
    /// `NormLast`
    norm_last: GpDir,
    /// `CurvLast`
    curv_last: f64,
    /// `TranLast`
    tran_last: Orientation,
}

impl Default for CurveTransition {
    fn default() -> Self {
        Self::new()
    }
}

impl CurveTransition {
    /// `TopTrans_CurveTransition()` (`cxx:26-36`): `myCurv = 0`, `Init = false`,
    /// both curvatures 0 and both transitions `FORWARD`.
    pub fn new() -> Self {
        Self {
            my_tgt: GpDir::default_dir(),
            my_norm: GpDir::default_dir(),
            my_curv: 0.0,
            init: false,
            tgt_first: GpDir::default_dir(),
            norm_first: GpDir::default_dir(),
            curv_first: 0.0,
            tran_first: Orientation::Forward,
            tgt_last: GpDir::default_dir(),
            norm_last: GpDir::default_dir(),
            curv_last: 0.0,
            tran_last: Orientation::Forward,
        }
    }

    /// `Reset(Tgt, Norm, Curv)` (`cxx:42-48`).
    pub fn reset(&mut self, tgt: &GpDir, norm: &GpDir, curv: f64) {
        self.my_tgt = *tgt;
        self.my_norm = *norm;
        self.my_curv = curv;
        self.init = true;
    }

    /// `Reset(Tgt)` (`cxx:56-61`): the straight-line initializer; `myCurv` is
    /// set to `0.` and `myNorm` is left as it was.
    pub fn reset_line(&mut self, tgt: &GpDir) {
        self.my_tgt = *tgt;
        self.my_curv = 0.0;
        self.init = true;
    }

    /// `Compare(Tole, T, N, C, St, Or)` (`cxx:69-273`).
    ///
    /// `st` is how the curve crosses the boundary, `or` how the element is
    /// oriented on the boundary.
    pub fn compare(
        &mut self,
        tole: f64,
        t: &GpDir,
        n: &GpDir,
        c: f64,
        st: Orientation,
        or: Orientation,
    ) {
        // S is the transition, how the curve cross the boundary
        // O is the orientation, how the intersection is set on the boundary
        let mut s = st;
        let o = or;

        // adjustment for INTERNAL transition
        if s == Orientation::Internal {
            if t.dot(&self.my_tgt) < 0.0 {
                s = o.reversed();
            } else {
                s = o;
            }
        }

        // It is the first comparison for this complex transition
        if self.init {
            self.init = false;
            self.tgt_first = *t;
            self.norm_first = *n;
            self.curv_first = c;
            self.tran_first = s;
            self.tgt_last = *t;
            self.norm_last = *n;
            self.curv_last = c;
            self.tran_last = s;
            match o {
                // Interference en fin d'arete il faut inverser la tangente
                Orientation::Reversed => {
                    self.tgt_first.reverse();
                    self.tgt_last.reverse();
                }
                Orientation::Internal => {
                    // Interference en milieu d'arete il faut inverser en fonction de la
                    // position de la tangente de reference
                    if self.my_tgt.dot(t) > 0.0 {
                        self.tgt_first.reverse();
                    } else {
                        self.tgt_last.reverse();
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
        }
        // Compare with the existent first and last transition :
        else {
            let mut first_set = false;
            let mut cos_ang_with_t = self.my_tgt.dot(t);
            match o {
                Orientation::Reversed => {
                    cos_ang_with_t = -cos_ang_with_t;
                }
                Orientation::Internal => {
                    if cos_ang_with_t > 0.0 {
                        cos_ang_with_t = -cos_ang_with_t;
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
            let cos_ang_with1 = self.my_tgt.dot(&self.tgt_first);

            match self.compare_angles(cos_ang_with_t, cos_ang_with1, tole) {
                LOWER => {
                    // If the angle is greater than the first the new become the first
                    first_set = true;
                    self.tgt_first = *t;
                    match o {
                        Orientation::Reversed => self.tgt_first.reverse(),
                        Orientation::Internal => {
                            if self.my_tgt.dot(t) > 0.0 {
                                self.tgt_first.reverse();
                            }
                        }
                        Orientation::Forward | Orientation::External => {}
                    }
                    self.norm_first = *n;
                    self.curv_first = c;
                    self.tran_first = s;
                }
                SAME => {
                    // If same angles we look at the Curvature
                    if self.is_before(tole, cos_ang_with_t, n, c, &self.norm_first, self.curv_first)
                    {
                        first_set = true;
                        self.tgt_first = *t;
                        match o {
                            Orientation::Reversed => self.tgt_first.reverse(),
                            Orientation::Internal => {
                                if self.my_tgt.dot(t) > 0.0 {
                                    self.tgt_first.reverse();
                                }
                            }
                            Orientation::Forward | Orientation::External => {}
                        }
                        self.norm_first = *n;
                        self.curv_first = c;
                        self.tran_first = s;
                    }
                }
                _ => {
                    // GREATER
                }
            }

            if !first_set || o == Orientation::Internal {
                // Dans les cas de tangence le premier peut etre aussi le dernier
                if o == Orientation::Internal {
                    cos_ang_with_t = -cos_ang_with_t;
                }
                let cos_ang_with2 = self.my_tgt.dot(&self.tgt_last);

                match self.compare_angles(cos_ang_with_t, cos_ang_with2, tole) {
                    GREATER => {
                        // If the angle is lower than the last the new become the last
                        self.tgt_last = *t;
                        match o {
                            Orientation::Reversed => self.tgt_last.reverse(),
                            Orientation::Internal => {
                                if self.my_tgt.dot(t) < 0.0 {
                                    self.tgt_last.reverse();
                                }
                            }
                            Orientation::Forward | Orientation::External => {}
                        }
                        self.norm_last = *n;
                        self.curv_last = c;
                        self.tran_last = s;
                    }
                    SAME => {
                        // If the angle is the same we look at the curvature
                        if self.is_before(tole, cos_ang_with_t, n, c, &self.norm_last, self.curv_last)
                        {
                            self.tgt_last = *t;
                            match o {
                                Orientation::Reversed => self.tgt_last.reverse(),
                                Orientation::Internal => {
                                    if self.my_tgt.dot(t) < 0.0 {
                                        self.tgt_last.reverse();
                                    }
                                }
                                Orientation::Forward | Orientation::External => {}
                            }
                            self.norm_last = *n;
                            self.curv_last = c;
                            self.tran_last = s;
                        }
                    }
                    _ => {
                        // LOWER
                    }
                }
            }
        }
    }

    /// `StateBefore()` (`cxx:280-296`).
    pub fn state_before(&self) -> State {
        if self.init {
            return State::Unknown;
        }
        match self.tran_first {
            Orientation::Forward | Orientation::External => State::Out,
            Orientation::Reversed | Orientation::Internal => State::In,
        }
    }

    /// `StateAfter()` (`cxx:304-319`).
    pub fn state_after(&self) -> State {
        if self.init {
            return State::Unknown;
        }
        match self.tran_last {
            Orientation::Forward | Orientation::Internal => State::In,
            Orientation::Reversed | Orientation::External => State::Out,
        }
    }

    /// `IsBefore(Tole, CosAngl, N1, C1, N2, C2)` (`cxx:327-424`): `(N1, C1)` is
    /// before `(N2, C2)` in the edge orientation.
    fn is_before(
        &self,
        tole: f64,
        cos_angl: f64,
        n1: &GpDir,
        c1: f64,
        n2: &GpDir,
        c2: f64,
    ) -> bool {
        let tn1 = self.my_tgt.dot(n1);
        let tn2 = self.my_tgt.dot(n2);
        let mut one_before = false;

        if tn1.abs() <= tole || tn2.abs() <= tole {
            // Tangent : The first is the interference which have the nearest curvature
            //           from the reference.
            if self.my_curv == 0.0 {
                // The reference is straight
                // The first is the interference which have the lowest curvature.
                if c1 < c2 {
                    one_before = true;
                }
                //  Modified by Sergey KHROMOV - Wed Dec 27 17:08:49 2000 Begin
                if cos_angl > 0.0 {
                    one_before = !one_before;
                }
                //  Modified by Sergey KHROMOV - Wed Dec 27 17:08:50 2000 End
            } else {
                // The reference is curv
                // The first is the interference which have the nearest curvature
                // in the direction
                let delta_c1;
                let delta_c2;
                if c1 == 0.0 || self.my_curv == 0.0 {
                    delta_c1 = c1 - self.my_curv;
                } else {
                    delta_c1 = (c1 - self.my_curv) * n1.dot(&self.my_norm);
                }
                if c2 == 0.0 || self.my_curv == 0.0 {
                    delta_c2 = c2 - self.my_curv;
                } else {
                    delta_c2 = (c2 - self.my_curv) * n2.dot(&self.my_norm);
                }
                if delta_c1 < delta_c2 {
                    one_before = true;
                }
                if cos_angl > 0.0 {
                    one_before = !one_before;
                }
            }
        } else if tn1 < 0.0 {
            // Before the first interference we are in the curvature
            if tn2 > 0.0 {
                // Before the second  interference we are out the curvature
                // The first interference is before  /* ->)( */
                one_before = true;
            } else {
                // Before the second interference we are in the curvature
                if c1 > c2 {
                    // We choice the greater curvature
                    // The first interference is before   /* ->)) */
                    one_before = true;
                }
            }
        } else if tn1 > 0.0 {
            // Before the first interference we are out the curvature
            if tn2 > 0.0 {
                // Before the second interference we are out the curvature /* ->(( */
                if c1 < c2 {
                    // We choice the lower curvature
                    // The first interference is before
                    one_before = true;
                }
            }
        }
        one_before
    }

    /// `Compare(Ang1, Ang2, Tole)` (`cxx:428-441`): `GREATER`/`SAME`/`LOWER`.
    fn compare_angles(&self, ang1: f64, ang2: f64, tole: f64) -> i32 {
        let mut res = SAME;
        if ang1 - ang2 > tole {
            res = GREATER;
        } else if ang2 - ang1 > tole {
            res = LOWER;
        }
        res
    }
}
