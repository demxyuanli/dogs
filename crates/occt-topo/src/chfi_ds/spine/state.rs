//! The extremity targets, the reference and the status flags of `ChFiDS_Spine`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Spine.lxx:20-170`,
//! `ChFiDS_Spine.cxx:353-500`, `:919-932` (OCCT 8.0.0).

use super::*;

impl Spine {
    /// `SetFirstTgt(W)` (`cxx:353-369`): throws `Standard_Failure` on a
    /// periodic guideline; the port answers `None`. The `hasfirsttgt` flag is
    /// suspended around the `D1` call exactly as in OCCT, so that a `D1` inside
    /// the extension does not stop on the target.
    pub fn set_first_tgt(&mut self, w: f64) -> Option<()> {
        if self.is_periodic() {
            return None;
        }
        self.hasfirsttgt = false;
        let (p, v) = self.d1(w)?;
        self.firstori = p;
        self.firsttgt = v;
        self.hasfirsttgt = true;
        self.firsttgtpar = w;
        Some(())
    }

    /// `SetLastTgt(W)` (`cxx:373-392`).
    pub fn set_last_tgt(&mut self, w: f64) -> Option<()> {
        if self.is_periodic() {
            return None;
        }
        self.haslasttgt = false;
        let (p, v) = self.d1(w)?;
        self.lastori = p;
        self.lasttgt = v;
        self.haslasttgt = true;
        self.lasttgtpar = w;
        Some(())
    }

    /// `HasFirstTgt()` (`cxx:396-399`).
    pub fn has_first_tgt(&self) -> bool {
        self.hasfirsttgt
    }

    /// `HasLastTgt()` (`cxx:403-406`).
    pub fn has_last_tgt(&self) -> bool {
        self.haslasttgt
    }

    /// `SetReference(W)` (`cxx:410-421`).
    pub fn set_reference(&mut self, w: f64) -> Option<()> {
        self.hasref = true;
        let lll = self.abscissa_upper()?;
        self.valref = if self.is_periodic() {
            clib::in_period(w, 0.0, lll)
        } else {
            w
        };
        Some(())
    }

    /// `SetReference(I)` (`cxx:425-437`): the middle of edge `I`, 1-based.
    pub fn set_reference_of(&mut self, i: i32) -> Option<()> {
        self.hasref = true;
        self.valref = if i == 1 {
            self.abscissa_value(1)? * 0.5
        } else {
            (self.abscissa_value(i)? + self.abscissa_value(i - 1)?) * 0.5
        };
        Some(())
    }

    /// `UnsetReference()` (`cxx:497-500`).
    pub fn unset_reference(&mut self) {
        self.hasref = false;
    }

    /// `SetErrorStatus(state)` (`cxx:919-922`).
    pub fn set_error_status(&mut self, state: ErrorStatus) {
        self.errorstate = state;
    }

    /// `ErrorStatus()` (`cxx:927-932`).
    pub fn error_status(&self) -> ErrorStatus {
        self.errorstate
    }

    /// `FirstStatus()` (`lxx:48-51`).
    pub fn first_status(&self) -> State {
        self.first_state
    }

    /// `LastStatus()` (`lxx:55-58`).
    pub fn last_status(&self) -> State {
        self.last_state
    }

    /// `Status(IsFirst)` (`lxx:75-81`).
    pub fn status(&self, is_first: bool) -> State {
        if is_first {
            self.first_state
        } else {
            self.last_state
        }
    }

    /// `SetStatus(S, IsFirst)` (`lxx:63-70`).
    pub fn set_status(&mut self, s: State, is_first: bool) {
        if is_first {
            self.first_state = s;
        } else {
            self.last_state = s;
        }
    }

    /// `GetTypeOfConcavity()` (`lxx:44-47`).
    pub fn get_type_of_concavity(&self) -> TypeOfConcavity {
        self.my_type_of_concavity
    }

    /// `SetTypeOfConcavity(theType)` (`lxx:22-25`).
    pub fn set_type_of_concavity(&mut self, the_type: TypeOfConcavity) {
        self.my_type_of_concavity = the_type;
    }

    /// `IsTangencyExtremity(IsFirst)` (`lxx:93-99`).
    pub fn is_tangency_extremity(&self, is_first: bool) -> bool {
        if is_first {
            self.firstistgt
        } else {
            self.lastistgt
        }
    }

    /// `SetTangencyExtremity(IsTangency, IsFirst)` (`lxx:85-90`).
    pub fn set_tangency_extremity(&mut self, is_tangency: bool, is_first: bool) {
        if is_first {
            self.firstistgt = is_tangency;
        } else {
            self.lastistgt = is_tangency;
        }
    }

    /// `Mode()` (`lxx:158-162`).
    pub fn mode(&self) -> ChamfMode {
        self.my_mode
    }

    /// `ChFiDS_ChamfSpine::SetMode(theMode)`
    /// (`ChFiDS_ChamfSpine.cxx:115-118`): OCCT reaches `myMode` from the
    /// derived class, so the writer lives with the reader.
    pub fn set_mode(&mut self, the_mode: ChamfMode) {
        self.my_mode = the_mode;
    }

    /// `GetTolesp()` (`lxx:166-170`).
    pub fn get_tolesp(&self) -> f64 {
        self.tolesp
    }
}
