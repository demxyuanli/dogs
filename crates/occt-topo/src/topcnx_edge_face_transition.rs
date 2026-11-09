//! `TopCnx_EdgeFaceTransition`
//! (`ModelingAlgorithms/TKHLR/TopCnx/TopCnx_EdgeFaceTransition.{hxx,cxx}`,
//! `hxx:29-67`, `cxx:22-134`): the cumulated transition of a bundle of
//! interferences sharing one geometric locus on an edge.
//!
//! `AddInterference` (`cxx:46-68`) forwards the local curve data to
//! [`CurveTransition::compare`] and counts the two boundary orientations, so
//! `Transition()` (`cxx:74-111`) is a pure 2x2 table over the curve
//! transition's `(StateBefore, StateAfter)` pair, and `BoundaryTransition()`
//! (`cxx:117-134`) is the majority vote of the counted boundary orientations.
//!
//! The `OCCT_DEBUG` printing of the unprocessed states (`cxx:88-90`,
//! `:102-104`, `:107-109`) has no port counterpart; the three arms that raise it
//! all fall through to the same `TopAbs_INTERNAL` (`cxx:111`), which the port
//! returns directly.

use occt_core::gp::GpDir;

use crate::abs::{Orientation, State};
use crate::top_trans_curve_transition::CurveTransition;

/// `TopCnx_EdgeFaceTransition` (`hxx:29-67`).
#[derive(Debug, Clone)]
pub struct EdgeFaceTransition {
    /// `TopTrans_CurveTransition myCurveTransition` (`hxx:67`).
    my_curve_transition: CurveTransition,
    /// `int nbBoundForward` (`hxx:68`).
    nb_bound_forward: i32,
    /// `int nbBoundReversed` (`hxx:69`).
    nb_bound_reversed: i32,
}

impl Default for EdgeFaceTransition {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeFaceTransition {
    /// `TopCnx_EdgeFaceTransition()` (`cxx:22-26`): both counters zero, the
    /// curve transition left as its own default constructor.
    pub fn new() -> Self {
        Self {
            my_curve_transition: CurveTransition::new(),
            nb_bound_forward: 0,
            nb_bound_reversed: 0,
        }
    }

    /// `Reset(Tgt, Norm, Curv)` (`cxx:28-33`).
    pub fn reset(&mut self, tgt: &GpDir, norm: &GpDir, curv: f64) {
        self.my_curve_transition.reset(tgt, norm, curv);
        self.nb_bound_forward = 0;
        self.nb_bound_reversed = 0;
    }

    /// `Reset(Tgt)` (`cxx:37-42`): the linear-edge initializer.
    pub fn reset_line(&mut self, tgt: &GpDir) {
        self.my_curve_transition.reset_line(tgt);
        self.nb_bound_forward = 0;
        self.nb_bound_reversed = 0;
    }

    /// `AddInterference(Tole, Tang, Norm, Curv, Or, Tr, BTr)` (`cxx:46-68`).
    ///
    /// `or` is the orientation of the interference on the boundary curve, `tr`
    /// its transition and `b_tr` its boundary transition.
    pub fn add_interference(
        &mut self,
        tole: f64,
        tang: &GpDir,
        norm: &GpDir,
        curv: f64,
        or: Orientation,
        tr: Orientation,
        b_tr: Orientation,
    ) {
        self.my_curve_transition
            .compare(tole, tang, norm, curv, tr, or);
        match b_tr {
            Orientation::Forward => self.nb_bound_forward += 1,
            Orientation::Reversed => self.nb_bound_reversed += 1,
            Orientation::Internal | Orientation::External => {}
        }
    }

    /// `Transition()` (`cxx:74-111`).
    pub fn transition(&self) -> Orientation {
        let bef = self.my_curve_transition.state_before();
        let aft = self.my_curve_transition.state_after();
        if bef == State::In {
            if aft == State::In {
                return Orientation::Internal;
            } else if aft == State::Out {
                return Orientation::Reversed;
            }
        } else if bef == State::Out {
            if aft == State::In {
                return Orientation::Forward;
            } else if aft == State::Out {
                return Orientation::External;
            }
        }
        Orientation::Internal
    }

    /// `BoundaryTransition()` (`cxx:117-134`): the three-way vote, where the
    /// tie arm returns `TopAbs_EXTERNAL` in both branches (`cxx:125-133`).
    pub fn boundary_transition(&self) -> Orientation {
        if self.nb_bound_forward > self.nb_bound_reversed {
            Orientation::Forward
        } else if self.nb_bound_forward < self.nb_bound_reversed {
            Orientation::Reversed
        } else {
            Orientation::External
        }
    }
}
