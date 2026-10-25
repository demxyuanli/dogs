//! The five `ChFiDS` enumerations.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKFillet/ChFiDS/`):
//! `ChFiDS_State.hxx:32-40`, `ChFiDS_TypeOfConcavity.hxx:21-29`,
//! `ChFiDS_ChamfMethod.hxx:20-25`, `ChFiDS_ChamfMode.hxx:19-35`,
//! `ChFiDS_ErrorStatus.hxx:21-29`.

/// `ChFiDS_State` (`ChFiDS_State.hxx:32-40`): the kinds of extremity of a
/// fillet.
///
/// `OnSame`, `OnDiff` and `AllSame` are the particular cases of `BreakPoint`
/// for a corner with three edges and three faces — the header spells the
/// distinction out at `hxx:23-30`.
/// No `Default` impl: OCCT's initial value is a per-member choice
/// (`ChFiDS_Spine.cxx:43-44` starts both ends at `ChFiDS_OnSame`), so the port
/// writes it at the use site rather than guessing a package-wide default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    /// `ChFiDS_OnSame`
    OnSame,
    /// `ChFiDS_OnDiff`
    OnDiff,
    /// `ChFiDS_AllSame`
    AllSame,
    /// `ChFiDS_BreakPoint`
    BreakPoint,
    /// `ChFiDS_FreeBoundary`
    FreeBoundary,
    /// `ChFiDS_Closed`
    Closed,
    /// `ChFiDS_Tangent`
    Tangent,
}

/// `ChFiDS_TypeOfConcavity` (`ChFiDS_TypeOfConcavity.hxx:21-29`).
/// No `Default` impl for the same reason as [`State`]: `ChFiDS_Spine.cxx:43`
/// initializes its own member to `ChFiDS_Other`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeOfConcavity {
    /// `ChFiDS_Concave`
    Concave,
    /// `ChFiDS_Convex`
    Convex,
    /// `ChFiDS_Tangential`
    Tangential,
    /// `ChFiDS_FreeBound`
    FreeBound,
    /// `ChFiDS_Other`
    Other,
    /// `ChFiDS_Mixed`
    Mixed,
}

/// `ChFiDS_ChamfMethod` (`ChFiDS_ChamfMethod.hxx:20-25`).
/// No `Default` impl: `ChFiDS_Spine.cxx:39` initializes its own member to
/// `ChFiDS_ClassicChamfer`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamfMethod {
    /// `ChFiDS_Sym`
    Sym,
    /// `ChFiDS_TwoDist`
    TwoDist,
    /// `ChFiDS_DistAngle`
    DistAngle,
}

/// `ChFiDS_ChamfMode` (`ChFiDS_ChamfMode.hxx:19-35`).
/// No `Default` impl: `ChFiDS_Spine.cxx:39` initializes its own member to
/// `ChFiDS_ClassicChamfer`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamfMode {
    /// `ChFiDS_ClassicChamfer`: constant distance from the spine to one of the
    /// two surfaces (`hxx:20-21`).
    ClassicChamfer,
    /// `ChFiDS_ConstThroatChamfer`: symmetric chamfer with constant throat,
    /// the height of the isosceles triangle in section (`hxx:23-25`).
    ConstThroatChamfer,
    /// `ChFiDS_ConstThroatWithPenetrationChamfer`: right-angled section whose
    /// leg from apex to top is the constant throat (`hxx:27-34`).
    ConstThroatWithPenetrationChamfer,
}

/// `ChFiDS_ErrorStatus` (`ChFiDS_ErrorStatus.hxx:21-29`): the cause of an
/// error.
/// No `Default` impl: `ChFiDS_Spine.cxx:57` initializes its own member to
/// `ChFiDS_Ok`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorStatus {
    /// `ChFiDS_Ok`
    Ok,
    /// `ChFiDS_Error`
    Error,
    /// `ChFiDS_WalkingFailure`
    WalkingFailure,
    /// `ChFiDS_StartsolFailure`
    StartsolFailure,
    /// `ChFiDS_TwistedSurface`
    TwistedSurface,
}
