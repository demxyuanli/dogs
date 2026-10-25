//! `IntPatch_Line` — the state every intersection line of two surfaces
//! carries: whether it is a line of tangency, its (constant) transition on
//! each surface with the matching situation when the transition is `Touch`,
//! and the four "is an iso-parametric curve" flags.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_Line.hxx:35-137`,
//! `.lxx:17-94`, `.cxx:21-65`.
//!
//! Translation notes:
//! - `typ` (`hxx:124`) is the C++ dynamic type, which the port's
//!   `GeomIntLine` sum type already carries (`arc_type`), so it is not
//!   duplicated here.
//! - `SituationS1` / `SituationS2` (`lxx:51-69`) raise
//!   `Standard_DomainError` when the transition is not `Touch`; the port
//!   returns `Err`.
//! - Every `IntPatch_*Line` constructor chains to one of the three
//!   `IntPatch_Line` constructors (`cxx:21-65`), so the port's line structs
//!   hold this type by value and the three constructors become
//!   [`IntPatchLine::with_transitions`], [`IntPatchLine::with_situations`] and
//!   [`IntPatchLine::undecided`].

use crate::intsurf_transition::{IntSurfSituation, IntSurfTypeTrans};

/// `class IntPatch_Line` (`hxx:35-137`), without the type tag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntPatchLine {
    /// `tg` (`hxx:127`).
    tg: bool,
    /// `tS1` (`hxx:128`).
    t_s1: IntSurfTypeTrans,
    /// `tS2` (`hxx:129`).
    t_s2: IntSurfTypeTrans,
    /// `sit1` (`hxx:130`).
    sit1: IntSurfSituation,
    /// `sit2` (`hxx:131`).
    sit2: IntSurfSituation,
    /// `uS1` / `vS1` (`hxx:132-133`).
    u_s1: bool,
    v_s1: bool,
    /// `uS2` / `vS2` (`hxx:134-135`).
    u_s2: bool,
    v_s2: bool,
}

impl Default for IntPatchLine {
    /// `IntPatch_Line(const bool Tang)` (`cxx:51-65`).
    fn default() -> Self {
        Self::undecided(false)
    }
}

impl IntPatchLine {
    /// `IntPatch_Line(Tang, Trans1, Trans2)` (`cxx:21-34`).
    pub fn with_transitions(tang: bool, trans1: IntSurfTypeTrans, trans2: IntSurfTypeTrans) -> Self {
        Self {
            tg: tang,
            t_s1: trans1,
            t_s2: trans2,
            sit1: IntSurfSituation::Unknown,
            sit2: IntSurfSituation::Unknown,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// `IntPatch_Line(Tang, Situ1, Situ2)` (`cxx:36-49`): both transitions are
    /// forced to `Touch`.
    pub fn with_situations(tang: bool, situ1: IntSurfSituation, situ2: IntSurfSituation) -> Self {
        Self {
            tg: tang,
            t_s1: IntSurfTypeTrans::Touch,
            t_s2: IntSurfTypeTrans::Touch,
            sit1: situ1,
            sit2: situ2,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// `IntPatch_Line(Tang)` (`cxx:51-65`).
    pub fn undecided(tang: bool) -> Self {
        Self {
            tg: tang,
            t_s1: IntSurfTypeTrans::Undecided,
            t_s2: IntSurfTypeTrans::Undecided,
            sit1: IntSurfSituation::Unknown,
            sit2: IntSurfSituation::Unknown,
            u_s1: false,
            v_s1: false,
            u_s2: false,
            v_s2: false,
        }
    }

    /// `SetValue(Uiso1, Viso1, Uiso2, Viso2)` (`lxx:17-26`).
    pub fn set_value(&mut self, u_iso1: bool, v_iso1: bool, u_iso2: bool, v_iso2: bool) {
        self.u_s1 = u_iso1;
        self.v_s1 = v_iso1;
        self.u_s2 = u_iso2;
        self.v_s2 = v_iso2;
    }

    /// `IsTangent()` (`lxx:33-37`).
    pub fn is_tangent(&self) -> bool {
        self.tg
    }

    /// `TransitionOnS1()` (`lxx:39-43`).
    pub fn transition_on_s1(&self) -> IntSurfTypeTrans {
        self.t_s1
    }

    /// `TransitionOnS2()` (`lxx:45-49`).
    pub fn transition_on_s2(&self) -> IntSurfTypeTrans {
        self.t_s2
    }

    /// `SituationS1()` (`lxx:51-59`).
    pub fn situation_s1(&self) -> Result<IntSurfSituation, String> {
        if self.t_s1 != IntSurfTypeTrans::Touch {
            return Err("IntPatch_Line::SituationS1: Standard_DomainError".to_string());
        }
        Ok(self.sit1)
    }

    /// `SituationS2()` (`lxx:61-69`).
    pub fn situation_s2(&self) -> Result<IntSurfSituation, String> {
        if self.t_s2 != IntSurfTypeTrans::Touch {
            return Err("IntPatch_Line::SituationS2: Standard_DomainError".to_string());
        }
        Ok(self.sit2)
    }

    /// `IsUIsoOnS1()` (`lxx:71-75`).
    pub fn is_u_iso_on_s1(&self) -> bool {
        self.u_s1
    }

    /// `IsVIsoOnS1()` (`lxx:77-81`).
    pub fn is_v_iso_on_s1(&self) -> bool {
        self.v_s1
    }

    /// `IsUIsoOnS2()` (`lxx:83-87`).
    pub fn is_u_iso_on_s2(&self) -> bool {
        self.u_s2
    }

    /// `IsVIsoOnS2()` (`lxx:89-93`).
    pub fn is_v_iso_on_s2(&self) -> bool {
        self.v_s2
    }
}
