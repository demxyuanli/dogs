//! `TopOpeBRepDS_Transition` — how an interference separates the existing
//! geometry in inside and outside.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Transition.hxx:29-125`,
//! `TopOpeBRepDS_Transition.cxx:20-364`.

use crate::abs::{Orientation, ShapeType, State};

/// Source: `TopOpeBRepDS_Transition.hxx:118-124`
#[derive(Debug, Clone)]
pub struct Transition {
    state_before: State,
    state_after: State,
    shape_before: ShapeType,
    shape_after: ShapeType,
    index_before: i32,
    index_after: i32,
}

impl Default for Transition {
    /// `TopOpeBRepDS_Transition::TopOpeBRepDS_Transition()` (`cxx:22-30`).
    fn default() -> Self {
        Self {
            state_before: State::Unknown,
            state_after: State::Unknown,
            shape_before: ShapeType::Face,
            shape_after: ShapeType::Face,
            index_before: 0,
            index_after: 0,
        }
    }
}

impl Transition {
    /// `TopOpeBRepDS_Transition(SB, SA, ONB, ONA)` (`cxx:34-46`).
    pub fn new(
        state_before: State,
        state_after: State,
        shape_before: ShapeType,
        shape_after: ShapeType,
    ) -> Self {
        Self {
            state_before,
            state_after,
            shape_before,
            shape_after,
            index_before: 0,
            index_after: 0,
        }
    }

    /// `TopOpeBRepDS_Transition(const TopAbs_Orientation O)` (`cxx:48-56`):
    /// the shape enums keep their in-class defaults (`TopAbs_FACE`) because the
    /// OCCT member-init list does not list the states; `Set(O)` writes them.
    pub fn from_orientation(o: Orientation) -> Self {
        let mut t = Self::default();
        t.set_orientation(o);
        t
    }

    /// `Set(SB, SA, ONB, ONA)` (`cxx:62-70`).
    pub fn set(
        &mut self,
        state_before: State,
        state_after: State,
        shape_before: ShapeType,
        shape_after: ShapeType,
    ) {
        self.state_before = state_before;
        self.state_after = state_after;
        self.shape_before = shape_before;
        self.shape_after = shape_after;
    }

    /// `StateBefore` (`cxx:74-77`).
    pub fn set_state_before(&mut self, s: State) {
        self.state_before = s;
    }

    /// `StateAfter` (`cxx:80-83`).
    pub fn set_state_after(&mut self, s: State) {
        self.state_after = s;
    }

    /// `ShapeBefore` (`cxx:87-90`).
    pub fn set_shape_before(&mut self, se: ShapeType) {
        self.shape_before = se;
    }

    /// `ShapeAfter` (`cxx:94-97`).
    pub fn set_shape_after(&mut self, se: ShapeType) {
        self.shape_after = se;
    }

    /// `Before(S, ShapeBefore)` (`cxx:101-105`).
    pub fn before(&mut self, s: State, shape_before: ShapeType) {
        self.state_before = s;
        self.shape_before = shape_before;
    }

    /// `After(S, ShapeAfter)` (`cxx:109-113`).
    pub fn after(&mut self, s: State, shape_after: ShapeType) {
        self.state_after = s;
        self.shape_after = shape_after;
    }

    /// `Index(I)` (`cxx:117-120`).
    pub fn set_index(&mut self, i: i32) {
        self.index_before = i;
        self.index_after = i;
    }

    /// `IndexBefore(I)` (`cxx:124-127`).
    pub fn set_index_before(&mut self, i: i32) {
        self.index_before = i;
    }

    /// `IndexAfter(I)` (`cxx:131-134`).
    pub fn set_index_after(&mut self, i: i32) {
        self.index_after = i;
    }

    /// `Before()` (`cxx:138-141`).
    pub fn before_state(&self) -> State {
        self.state_before
    }

    /// `ONBefore()` (`cxx:145-148`).
    pub fn on_before(&self) -> ShapeType {
        self.shape_before
    }

    /// `After()` (`cxx:152-155`).
    pub fn after_state(&self) -> State {
        self.state_after
    }

    /// `ONAfter()` (`cxx:159-162`).
    pub fn on_after(&self) -> ShapeType {
        self.shape_after
    }

    /// `ShapeBefore()` (`cxx:166-169`).
    pub fn shape_before(&self) -> ShapeType {
        self.shape_before
    }

    /// `ShapeAfter()` (`cxx:173-176`).
    pub fn shape_after(&self) -> ShapeType {
        self.shape_after
    }

    /// `Index()` (`cxx:180-187`): throws `Standard_Failure` when the two
    /// indexes differ.
    pub fn index(&self) -> Result<i32, String> {
        if self.index_after != self.index_before {
            return Err("Standard_Failure: Transition::Index() on different shapes".into());
        }
        Ok(self.index_before)
    }

    /// `IndexBefore()` (`cxx:191-194`).
    pub fn index_before(&self) -> i32 {
        self.index_before
    }

    /// `IndexAfter()` (`cxx:198-201`).
    pub fn index_after(&self) -> i32 {
        self.index_after
    }

    /// `Set(const TopAbs_Orientation O)` (`cxx:206-229`).
    pub fn set_orientation(&mut self, o: Orientation) {
        match o {
            Orientation::Forward => {
                self.state_before = State::Out;
                self.state_after = State::In;
            }
            Orientation::Reversed => {
                self.state_before = State::In;
                self.state_after = State::Out;
            }
            Orientation::Internal => {
                self.state_before = State::In;
                self.state_after = State::In;
            }
            Orientation::External => {
                self.state_before = State::Out;
                self.state_after = State::Out;
            }
        }
    }

    /// `Orientation(S, T)` (`cxx:234-265`).
    pub fn orientation(&self, s: State, t: ShapeType) -> Orientation {
        if self.state_before == State::On || self.state_after == State::On {
            return self.orientation_on(s, t);
        }
        if self.state_before == s {
            if self.state_after == s {
                Orientation::Internal
            } else {
                Orientation::Reversed
            }
        } else if self.state_after == s {
            Orientation::Forward
        } else {
            Orientation::External
        }
    }

    /// `OrientationON(S, T)` (`cxx:269-312`). `T` is unused by OCCT (the
    /// parameter is commented out at the definition).
    fn orientation_on(&self, s: State, _t: ShapeType) -> Orientation {
        let mut result = Orientation::Forward;
        if self.state_before == State::On && self.state_after == State::On {
            if s == State::In {
                result = Orientation::Internal;
            } else if s == State::Out {
                result = Orientation::External;
            } else if s == State::On {
                result = Orientation::Internal;
            }
        } else if self.state_before == State::On {
            return if self.state_after == s {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
        } else if self.state_after == State::On {
            return if self.state_before == s {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };
        }
        result
    }

    /// `Complement()` (`cxx:316-357`). The indexes are copied before any
    /// state is written, exactly as in OCCT.
    pub fn complement(&self) -> Transition {
        let mut t = Transition::default();
        t.index_before = self.index_before;
        t.index_after = self.index_after;

        // xpu : 14-01-98
        if self.state_before == State::Unknown && self.state_after == State::Unknown {
            t.set(self.state_after, self.state_before, self.shape_after, self.shape_before);
            return t;
        }

        if self.state_before == State::On || self.state_after == State::On {
            t.set(self.state_after, self.state_before, self.shape_after, self.shape_before);
        } else {
            let o = self.orientation(State::In, ShapeType::Face);
            match o {
                // (OUT,IN) --> (IN,OUT)
                Orientation::Forward => t.set(
                    State::In,
                    State::Out,
                    self.shape_before,
                    self.shape_after,
                ),
                // (IN,OUT) --> (OUT,IN)
                Orientation::Reversed => t.set(
                    State::Out,
                    State::In,
                    self.shape_before,
                    self.shape_after,
                ),
                // (OUT,OUT) --> (IN,IN)
                Orientation::External => t.set(
                    State::In,
                    State::In,
                    self.shape_before,
                    self.shape_after,
                ),
                // (IN,IN) --> (OUT,OUT)
                Orientation::Internal => t.set(
                    State::Out,
                    State::Out,
                    self.shape_before,
                    self.shape_after,
                ),
            }
        }
        t
    }

    /// `IsUnknown()` (`cxx:362-365`).
    pub fn is_unknown(&self) -> bool {
        self.state_before == State::Unknown && self.state_after == State::Unknown
    }
}
