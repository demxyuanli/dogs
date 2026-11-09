//! `BRepFeat_StatusError` — the failure codes reported by the feature classes.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_StatusError.hxx:21-51`
//! (the whole file is the enum). `BRepFeat::Print` (`BRepFeat.cxx:718-808`)
//! maps every code to its message.

/// `enum BRepFeat_StatusError` (`hxx:21-51`), same order and names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusError {
    /// `BRepFeat_OK`
    Ok,
    /// `BRepFeat_BadDirect`
    BadDirect,
    /// `BRepFeat_BadIntersect`
    BadIntersect,
    /// `BRepFeat_EmptyBaryCurve`
    EmptyBaryCurve,
    /// `BRepFeat_EmptyCutResult`
    EmptyCutResult,
    /// `BRepFeat_FalseSide`
    FalseSide,
    /// `BRepFeat_IncDirection`
    IncDirection,
    /// `BRepFeat_IncSlidFace`
    IncSlidFace,
    /// `BRepFeat_IncParameter`
    IncParameter,
    /// `BRepFeat_IncTypes`
    IncTypes,
    /// `BRepFeat_IntervalOverlap`
    IntervalOverlap,
    /// `BRepFeat_InvFirstShape`
    InvFirstShape,
    /// `BRepFeat_InvOption`
    InvOption,
    /// `BRepFeat_InvShape`
    InvShape,
    /// `BRepFeat_LocOpeNotDone`
    LocOpeNotDone,
    /// `BRepFeat_LocOpeInvNotDone`
    LocOpeInvNotDone,
    /// `BRepFeat_NoExtFace`
    NoExtFace,
    /// `BRepFeat_NoFaceProf`
    NoFaceProf,
    /// `BRepFeat_NoGluer`
    NoGluer,
    /// `BRepFeat_NoIntersectF`
    NoIntersectF,
    /// `BRepFeat_NoIntersectU`
    NoIntersectU,
    /// `BRepFeat_NoParts`
    NoParts,
    /// `BRepFeat_NoProjPt`
    NoProjPt,
    /// `BRepFeat_NotInitialized`
    NotInitialized,
    /// `BRepFeat_NotYetImplemented`
    NotYetImplemented,
    /// `BRepFeat_NullRealTool`
    NullRealTool,
    /// `BRepFeat_NullToolF`
    NullToolF,
    /// `BRepFeat_NullToolU`
    NullToolU,
}

impl Default for StatusError {
    /// The `BRepFeat_Form` constructor's initial value (`lxx:22-35`).
    fn default() -> Self {
        Self::NotInitialized
    }
}

/// `BRepFeat::Print(const BRepFeat_StatusError SE, Standard_OStream& S)`
/// (`BRepFeat.cxx:718-808`): the text each code writes, in the reference's
/// `switch` order.
pub fn print(se: StatusError) -> &'static str {
    match se {
        StatusError::Ok => "No error",
        StatusError::BadDirect => "Directions must be opposite",
        StatusError::BadIntersect => "Intersection failure",
        StatusError::EmptyBaryCurve => "Empty BaryCurve",
        StatusError::EmptyCutResult => "Failure in Cut : Empty resulting shape",
        StatusError::FalseSide => "Verify plane and wire orientation",
        StatusError::IncDirection => "Incoherent Direction for shapes From and Until",
        StatusError::IncSlidFace => "Sliding face not in Base shape",
        StatusError::IncParameter => "Incoherent Parameter : shape Until before shape From",
        StatusError::IncTypes => {
            "Invalid option for faces From and Until : 1 Support and 1 not"
        }
        StatusError::IntervalOverlap => "Shapes From and Until overlap",
        StatusError::InvFirstShape => "Invalid First shape : more than 1 face",
        StatusError::InvOption => "Invalid option",
        StatusError::InvShape => "Invalid shape",
        StatusError::LocOpeNotDone => "Local Operation not done",
        StatusError::LocOpeInvNotDone => "Local Operation : intersection line conflict",
        StatusError::NoExtFace => "No Extreme faces",
        StatusError::NoFaceProf => "No Face Profile",
        StatusError::NoGluer => "Gluer Failure",
        StatusError::NoIntersectF => "No intersection between Feature and shape From",
        StatusError::NoIntersectU => "No intersection between Feature and shape Until",
        StatusError::NoParts => "No parts of tool kept",
        StatusError::NoProjPt => "No projection points",
        StatusError::NotInitialized => "Fields not initialized",
        StatusError::NotYetImplemented => "Not yet implemented",
        StatusError::NullRealTool => "Real Tool : Null DPrism",
        StatusError::NullToolF => "Null Tool : Invalid type for shape Form",
        StatusError::NullToolU => "Null Tool : Invalid type for shape Until",
    }
}
