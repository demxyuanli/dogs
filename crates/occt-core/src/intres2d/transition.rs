//! Port of `IntRes2d_Transition` and its enums
//! (`src/ModelingAlgorithms/TKGeomAlgo/IntRes2d/`):
//! `IntRes2d_Position.hxx`, `IntRes2d_TypeTrans.hxx`,
//! `IntRes2d_Situation.hxx`, `IntRes2d_Transition.hxx/.lxx/.cxx`.
//!
//! The three enum values are part of the result contract read by
//! `ShapeAnalysis_Wire::CheckSelfIntersectingEdge`
//! (`ShapeAnalysis_Wire.cxx:1324-1329`) and
//! `ShapeAnalysis_Wire::CheckIntersectingEdges`
//! (`ShapeAnalysis_Wire.cxx:1490-1493`), which is why a sampled intersector
//! cannot stand in for `Geom2dInt_GInter`.

/// `IntRes2d_Position` (`IntRes2d_Position.hxx:20-25`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IntRes2dPosition {
    /// `IntRes2d_Head`
    Head,
    /// `IntRes2d_Middle`
    Middle,
    /// `IntRes2d_End`
    End,
}

/// `IntRes2d_TypeTrans` (`IntRes2d_TypeTrans.hxx:20-26`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IntRes2dTypeTrans {
    /// `IntRes2d_In`
    In,
    /// `IntRes2d_Out`
    Out,
    /// `IntRes2d_Touch`
    Touch,
    /// `IntRes2d_Undecided`
    Undecided,
}

/// `IntRes2d_Situation` (`IntRes2d_Situation.hxx:20-25`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IntRes2dSituation {
    /// `IntRes2d_Inside`
    Inside,
    /// `IntRes2d_Outside`
    Outside,
    /// `IntRes2d_Unknown`
    Unknown,
}

/// `IntRes2d_Transition` (`IntRes2d_Transition.hxx:37-116`).
#[derive(Clone, Copy, Debug)]
pub struct IntRes2dTransition {
    /// `tangent`
    tangent: bool,
    /// `posit`
    posit: IntRes2dPosition,
    /// `typetra`
    typetra: IntRes2dTypeTrans,
    /// `situat`
    situat: IntRes2dSituation,
    /// `oppos`
    oppos: bool,
}

impl Default for IntRes2dTransition {
    /// `IntRes2d_Transition()` (`IntRes2d_Transition.cxx:20-26`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntRes2dTransition {
    /// Empty constructor (`IntRes2d_Transition.cxx:20-26`).
    pub fn new() -> Self {
        Self {
            tangent: true,
            posit: IntRes2dPosition::Middle,
            typetra: IntRes2dTypeTrans::Undecided,
            situat: IntRes2dSituation::Unknown,
            oppos: false,
        }
    }

    /// `IntRes2d_Transition(Tangent, Pos, Type)` (`Transition.lxx:20-33`):
    /// creates an IN or OUT transition.
    pub fn in_out(tangent: bool, pos: IntRes2dPosition, typ: IntRes2dTypeTrans) -> Self {
        Self {
            tangent,
            posit: pos,
            typetra: typ,
            situat: IntRes2dSituation::Unknown,
            oppos: false,
        }
    }

    /// `IntRes2d_Transition(Tangent, Pos, Situ, Oppos)`
    /// (`Transition.lxx:35-45`): creates a TOUCH transition.
    pub fn touch(
        tangent: bool,
        pos: IntRes2dPosition,
        situ: IntRes2dSituation,
        oppos: bool,
    ) -> Self {
        Self {
            tangent,
            posit: pos,
            typetra: IntRes2dTypeTrans::Touch,
            situat: situ,
            oppos,
        }
    }

    /// `IntRes2d_Transition(Pos)` (`Transition.lxx:47-56`): creates an
    /// UNDECIDED transition.
    pub fn undecided(pos: IntRes2dPosition) -> Self {
        Self {
            tangent: true,
            posit: pos,
            typetra: IntRes2dTypeTrans::Undecided,
            situat: IntRes2dSituation::Unknown,
            oppos: false,
        }
    }

    /// `SetValue(Tangent, Pos, Type)` (`Transition.lxx:58-66`).
    pub fn set_in_out(
        &mut self,
        tangent: bool,
        pos: IntRes2dPosition,
        typ: IntRes2dTypeTrans,
    ) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = typ;
    }

    /// `SetValue(Tangent, Pos, Situ, Oppos)` (`Transition.lxx:68-78`).
    pub fn set_touch(
        &mut self,
        tangent: bool,
        pos: IntRes2dPosition,
        situ: IntRes2dSituation,
        oppos: bool,
    ) {
        self.tangent = tangent;
        self.posit = pos;
        self.typetra = IntRes2dTypeTrans::Touch;
        self.situat = situ;
        self.oppos = oppos;
    }

    /// `SetValue(Pos)` (`Transition.lxx:80-85`): only the position and the
    /// UNDECIDED type are written, as in OCCT.
    pub fn set_undecided(&mut self, pos: IntRes2dPosition) {
        self.posit = pos;
        self.typetra = IntRes2dTypeTrans::Undecided;
    }

    /// `SetPosition(Pos)` (`Transition.lxx:87-91`).
    pub fn set_position(&mut self, pos: IntRes2dPosition) {
        self.posit = pos;
    }

    /// `PositionOnCurve()` (`Transition.lxx:93-97`).
    pub fn position_on_curve(&self) -> IntRes2dPosition {
        self.posit
    }

    /// `TransitionType()` (`Transition.lxx:99-103`).
    pub fn transition_type(&self) -> IntRes2dTypeTrans {
        self.typetra
    }

    /// `IsTangent()` (`Transition.lxx:105-113`). OCCT throws
    /// `Standard_DomainError` when the type is UNDECIDED.
    pub fn is_tangent(&self) -> bool {
        if self.typetra == IntRes2dTypeTrans::Undecided {
            panic!("IntRes2d_Transition::IsTangent: undecided transition");
        }
        self.tangent
    }

    /// `Situation()` (`Transition.lxx:115-123`). OCCT throws
    /// `Standard_DomainError` when the type is not TOUCH.
    pub fn situation(&self) -> IntRes2dSituation {
        if self.typetra != IntRes2dTypeTrans::Touch {
            panic!("IntRes2d_Transition::Situation: not a touch transition");
        }
        self.situat
    }

    /// `IsOpposite()` (`Transition.lxx:125-133`). OCCT throws
    /// `Standard_DomainError` when the type is not TOUCH.
    pub fn is_opposite(&self) -> bool {
        if self.typetra != IntRes2dTypeTrans::Touch {
            panic!("IntRes2d_Transition::IsOpposite: not a touch transition");
        }
        self.oppos
    }
}
