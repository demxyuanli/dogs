//! `TopOpeBRepBuild_GTool` — the `GTopo` tables of fuse / cut / common.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GTool.hxx:25-57`,
//! `.cxx:20-244`.

use crate::abs::ShapeType;
use crate::topopebrep::ds::kind::Config;

use super::g_topo::GTopo;

/// `class TopOpeBRepBuild_GTool` (`hxx:30`).
pub struct GTool;

impl GTool {
    /// `GFusUnsh(t1, t2)` (`cxx:21-42`).
    pub fn g_fus_unsh(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            false,
            false,
            false,
            false,
            true,
            false,
            true,
            false,
            t1,
            t2,
            Config::UnshGeometry,
            Config::UnshGeometry,
        )
    }

    /// `GFusSame(t1, t2)` (`cxx:46-67`).
    pub fn g_fus_same(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            false,
            false,
            false,
            true,
            true,
            false,
            true,
            false,
            t1,
            t2,
            Config::SameOriented,
            Config::SameOriented,
        )
    }

    /// `GFusDiff(t1, t2)` (`cxx:71-92`).
    pub fn g_fus_diff(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            false,
            false,
            false,
            false,
            true,
            false,
            true,
            false,
            t1,
            t2,
            Config::DiffOriented,
            Config::SameOriented,
        )
    }

    /// `GCutUnsh(t1, t2)` (`cxx:96-117`).
    pub fn g_cut_unsh(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            false,
            false,
            true,
            false,
            false,
            false,
            t1,
            t2,
            Config::UnshGeometry,
            Config::UnshGeometry,
        )
    }

    /// `GCutSame(t1, t2)` (`cxx:121-142`).
    pub fn g_cut_same(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            false,
            false,
            true,
            false,
            false,
            false,
            t1,
            t2,
            Config::SameOriented,
            Config::SameOriented,
        )
    }

    /// `GCutDiff(t1, t2)` (`cxx:146-167`).
    pub fn g_cut_diff(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            false,
            true,
            true,
            false,
            false,
            false,
            t1,
            t2,
            Config::DiffOriented,
            Config::SameOriented,
        )
    }

    /// `GComUnsh(t1, t2)` (`cxx:171-192`).
    pub fn g_com_unsh(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            true,
            false,
            false,
            false,
            false,
            false,
            t1,
            t2,
            Config::UnshGeometry,
            Config::UnshGeometry,
        )
    }

    /// `GComSame(t1, t2)` (`cxx:196-217`).
    pub fn g_com_same(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            true,
            true,
            false,
            false,
            false,
            false,
            t1,
            t2,
            Config::SameOriented,
            Config::SameOriented,
        )
    }

    /// `GComDiff(t1, t2)` (`cxx:221-242`).
    pub fn g_com_diff(t1: ShapeType, t2: ShapeType) -> GTopo {
        GTopo::with_params(
            false,
            true,
            false,
            true,
            false,
            false,
            false,
            false,
            false,
            t1,
            t2,
            Config::DiffOriented,
            Config::SameOriented,
        )
    }

    /// `Dump(OS)` (`cxx:246-293`) — debug output only.
    pub fn dump() {}
}
