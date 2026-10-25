//! `TopOpeBRepBuild_LoopClassifier` — classifies loops in order to build areas.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_LoopClassifier.hxx:26-40`,
//! `.cxx:19-21`.

use crate::abs::State;

use super::loop_shape::LoopHandle;

/// `class TopOpeBRepBuild_LoopClassifier` (`hxx:29`).
///
/// `Compare` is non-const in OCCT (`hxx:37`) and mutates the classifier's
/// current-pair cache (`PaveClassifier.cxx:325-336`), hence `&mut self`.
pub trait LoopClassifier {
    /// `Compare(L1, L2)` (`hxx:37-38`): state of loop `L1` compared with
    /// loop `L2`.
    fn compare(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State;
}
