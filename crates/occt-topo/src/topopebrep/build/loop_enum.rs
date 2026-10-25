//! `TopOpeBRepBuild_LoopEnum` — which loops of a list of loops a comparison
//! test selects.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_LoopEnum.hxx`.

/// `enum TopOpeBRepBuild_LoopEnum`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LoopEnum {
    /// `TopOpeBRepBuild_ANYLOOP`: every loop is tested.
    AnyLoop,
    /// `TopOpeBRepBuild_BOUNDARY`: only the boundary loops are tested.
    Boundary,
    /// `TopOpeBRepBuild_BLOCK`: only the block loops are tested.
    Block,
}
