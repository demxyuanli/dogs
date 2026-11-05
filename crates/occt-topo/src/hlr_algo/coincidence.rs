//! Information on the "hiding" edge of an interference. Source:
//! `HLRAlgo_Coincidence`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_Coincidence.hxx:38-81`).

use crate::abs::State;

/// The Coincidence class is used in an Interference to store information on
/// the "hiding" edge.
///
/// 2D data: the face/edge index `FE` and the parameter of the intersection
/// point, needed when the intersection is at the extremity of the edge.
/// 3D data: the state of the edge before and after the intersection with the
/// face, needed when the intersection is "ON" the face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HlrAlgoCoincidence {
    fe: i32,
    param: f64,
    st_bef: State,
    st_aft: State,
}

impl Default for HlrAlgoCoincidence {
    /// `HLRAlgo_Coincidence()` (`hxx:43-49`).
    fn default() -> Self {
        Self {
            fe: 0,
            param: 0.0,
            st_bef: State::In,
            st_aft: State::In,
        }
    }
}

impl HlrAlgoCoincidence {
    /// `HLRAlgo_Coincidence()` (`hxx:43-49`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Set2D(FE, Param)` (`hxx:51-55`).
    pub fn set_2d(&mut self, fe: i32, param: f64) {
        self.fe = fe;
        self.param = param;
    }

    /// `SetState3D(stbef, staft)` (`hxx:57-61`).
    pub fn set_state_3d(&mut self, st_bef: State, st_aft: State) {
        self.st_bef = st_bef;
        self.st_aft = st_aft;
    }

    /// `Value2D(FE, Param)` (`hxx:63-67`).
    pub fn value_2d(&self) -> (i32, f64) {
        (self.fe, self.param)
    }

    /// `State3D(stbef, staft)` (`hxx:69-73`).
    pub fn state_3d(&self) -> (State, State) {
        (self.st_bef, self.st_aft)
    }
}
