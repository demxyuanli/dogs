//! `TopTrans_CurveTransition` — transition of a curve crossing a curvilinear
//! boundary, used by `IntPatch_ImpPrmIntersection::ComputeTangency`
//! (`IntPatch_ImpPrmIntersection.cxx:333` Reset, `:362`/`:412` Compare).
//! Source: `TopTrans_CurveTransition.cxx` (OCCT 8.0.0).

use occt_core::gp::{GpDir, GpVec};

use crate::abs::Orientation;
use crate::fclass2d::FaceState;

/// `#define GREATER 1` (`TopTrans_CurveTransition.cxx:20`).
const GREATER: i32 = 1;
/// `#define SAME 0` (`TopTrans_CurveTransition.cxx:21`).
const SAME: i32 = 0;
/// `#define LOWER -1` (`TopTrans_CurveTransition.cxx:22`).
const LOWER: i32 = -1;

/// `TopTrans_CurveTransition` (`TopTrans_CurveTransition.hxx:45-109`).
pub(crate) struct CurveTransition {
    my_tgt: GpDir,
    my_norm: GpDir,
    my_curv: f64,
    init: bool,
    tgt_first: GpDir,
    norm_first: GpDir,
    curv_first: f64,
    tran_first: Orientation,
    tgt_last: GpDir,
    norm_last: GpDir,
    curv_last: f64,
    tran_last: Orientation,
}

impl CurveTransition {
    /// `TopTrans_CurveTransition()` (`TopTrans_CurveTransition.cxx:26-34`):
    /// `Init = false`, both transitions FORWARD, curvatures 0.
    pub(crate) fn new() -> Self {
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

    /// `Reset(Tgt, Norm, Curv)` (`TopTrans_CurveTransition.cxx:42-48`).
    /// `Tgt` is a `gp_Dir` in OCCT, i.e. the normalised `gp_Vec`.
    pub(crate) fn reset(&mut self, tgt: &GpVec, norm: &GpDir, curv: f64) {
        let Ok(tgt) = GpDir::from_vec(tgt) else {
            // OCCT converts `gp_Vec` to `gp_Dir` implicitly; a zero vector
            // raises `Standard_ConstructionError`. Leave the transition
            // uninitialised so the caller rejects the point
            // (`StateBefore`/`StateAfter` return `TopAbs_UNKNOWN`).
            return;
        };
        self.my_tgt = tgt;
        self.my_norm = *norm;
        self.my_curv = curv;
        self.init = true;
    }

    /// `Reset(Tgt)` (`TopTrans_CurveTransition.cxx:54-61`): sets the tangent,
    /// clears the curvature and marks the transition initialised. The normal
    /// is left unchanged, as in OCCT.
    pub(crate) fn reset_tgt(&mut self, tgt: &GpDir) {
        self.my_tgt = *tgt;
        self.my_curv = 0.0;
        self.init = true;
    }

    /// `Compare(Tole, T, N, C, S, O)` (`TopTrans_CurveTransition.cxx:69-273`).
    /// `T` is a `gp_Dir` in OCCT, i.e. the normalised `gp_Vec`.
    pub(crate) fn compare(
        &mut self,
        tole: f64,
        t: &GpVec,
        n: &GpDir,
        c: f64,
        st: Orientation,
        or: Orientation,
    ) {
        let Ok(t) = GpDir::from_vec(t) else {
            return;
        };
        // S is the transition, how the curve crosses the boundary
        // O is the orientation, how the intersection is set on the boundary
        let mut s = st;
        let o = or;

        // adjustment for INTERNAL transition (`:82-92`)
        if s == Orientation::Internal {
            if t.dot(&self.my_tgt) < 0.0 {
                s = o.reversed();
            } else {
                s = o;
            }
        }

        if self.init {
            // It is the first comparison for this complex transition (`:95-129`)
            self.init = false;
            self.tgt_first = t;
            self.norm_first = *n;
            self.curv_first = c;
            self.tran_first = s;
            self.tgt_last = t;
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
                    // Interference en milieu d'arete il faut inverser en
                    // fonction de la position de la tangente de reference
                    if self.my_tgt.dot(&t) > 0.0 {
                        self.tgt_first.reverse();
                    } else {
                        self.tgt_last.reverse();
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
        } else {
            // Compare with the existent first and last transition (`:132-272`)
            let mut first_set = false;
            let mut cos_ang_with_t = self.my_tgt.dot(&t);
            match o {
                Orientation::Reversed => cos_ang_with_t = -cos_ang_with_t,
                Orientation::Internal => {
                    if cos_ang_with_t > 0.0 {
                        cos_ang_with_t = -cos_ang_with_t;
                    }
                }
                Orientation::Forward | Orientation::External => {}
            }
            let cos_ang_with_1 = self.my_tgt.dot(&self.tgt_first);

            match compare_angles(cos_ang_with_t, cos_ang_with_1, tole) {
                LOWER => {
                    // If the angle is greater than the first the new become the first
                    first_set = true;
                    self.tgt_first = t;
                    match o {
                        Orientation::Reversed => self.tgt_first.reverse(),
                        Orientation::Internal => {
                            if self.my_tgt.dot(&t) > 0.0 {
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
                        self.tgt_first = t;
                        match o {
                            Orientation::Reversed => self.tgt_first.reverse(),
                            Orientation::Internal => {
                                if self.my_tgt.dot(&t) > 0.0 {
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
                _ => {}
            }

            if !first_set || o == Orientation::Internal {
                // Dans les cas de tangence le premier peut etre aussi le dernier
                if o == Orientation::Internal {
                    cos_ang_with_t = -cos_ang_with_t;
                }
                let cos_ang_with_2 = self.my_tgt.dot(&self.tgt_last);

                match compare_angles(cos_ang_with_t, cos_ang_with_2, tole) {
                    GREATER => {
                        // If the angle is lower than the last the new become the last
                        self.tgt_last = t;
                        match o {
                            Orientation::Reversed => self.tgt_last.reverse(),
                            Orientation::Internal => {
                                if self.my_tgt.dot(&t) < 0.0 {
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
                        if self.is_before(tole, cos_ang_with_t, &self.norm_last, self.curv_last, n, c)
                        {
                            self.tgt_last = t;
                            match o {
                                Orientation::Reversed => self.tgt_last.reverse(),
                                Orientation::Internal => {
                                    if self.my_tgt.dot(&t) < 0.0 {
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
                    _ => {}
                }
            }
        }
    }

    /// `StateBefore()` (`TopTrans_CurveTransition.cxx:280-296`).
    pub(crate) fn state_before(&self) -> FaceState {
        if self.init {
            return FaceState::Unknown;
        }
        match self.tran_first {
            Orientation::Forward | Orientation::External => FaceState::Out,
            Orientation::Reversed | Orientation::Internal => FaceState::In,
        }
    }

    /// `StateAfter()` (`TopTrans_CurveTransition.cxx:303-319`).
    pub(crate) fn state_after(&self) -> FaceState {
        if self.init {
            return FaceState::Unknown;
        }
        match self.tran_last {
            Orientation::Forward | Orientation::Internal => FaceState::In,
            Orientation::Reversed | Orientation::External => FaceState::Out,
        }
    }

    /// `IsBefore(Tole, CosAngl, N1, C1, N2, C2)`
    /// (`TopTrans_CurveTransition.cxx:327-424`).
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
            // Tangent : The first is the interference which have the nearest
            // curvature from the reference.
            if self.my_curv == 0.0 {
                // The reference is straight
                // The first is the interference which have the lowest curvature.
                if c1 < c2 {
                    one_before = true;
                }
                if cos_angl > 0.0 {
                    one_before = !one_before;
                }
            } else {
                // The reference is curv
                // The first is the interference which have the nearest curvature
                // in the direction
                let delta_c1 = if c1 == 0.0 || self.my_curv == 0.0 {
                    c1 - self.my_curv
                } else {
                    (c1 - self.my_curv) * n1.dot(&self.my_norm)
                };
                let delta_c2 = if c2 == 0.0 || self.my_curv == 0.0 {
                    c2 - self.my_curv
                } else {
                    (c2 - self.my_curv) * n2.dot(&self.my_norm)
                };
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
                // Before the second interference we are out the curvature
                // The first interference is before  /* ->)( */
                one_before = true;
            } else if c1 > c2 {
                // We choice the greater curvature
                // The first interference is before   /* ->)) */
                one_before = true;
            }
        } else if tn1 > 0.0 && tn2 > 0.0 && c1 < c2 {
            // Before the second interference we are out the curvature /* ->(( */
            // We choice the lower curvature
            // The first interference is before
            one_before = true;
        }
        one_before
    }
}

impl Default for CurveTransition {
    fn default() -> Self {
        Self::new()
    }
}

/// `Compare(Ang1, Ang2, Tole)` (`TopTrans_CurveTransition.cxx:428-441`).
fn compare_angles(ang1: f64, ang2: f64, tole: f64) -> i32 {
    if ang1 - ang2 > tole {
        GREATER
    } else if ang2 - ang1 > tole {
        LOWER
    } else {
        SAME
    }
}
