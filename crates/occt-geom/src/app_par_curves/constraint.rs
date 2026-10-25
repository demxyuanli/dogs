//! `AppParCurves_Constraint` and `AppParCurves_ConstraintCouple` — the
//! per-point constraint tag and its (index, tag) pair.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_Constraint.hxx:20-31`
//! and `AppParCurves_ConstraintCouple.hxx:30-56` / `.cxx:16-41` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `AppParCurves_Constraint` is an unscoped C++ enum; the port keeps the
//!   variants and their declaration order.
//! - `AppParCurves_ConstraintCouple` defaults to `myIndex = -1` and
//!   `AppParCurves_NoConstraint` (`cxx:16-22`), including its
//!   `Index()` / `Constraint()` accessors and the two `Set*` mutators.

/// `AppParCurves_Constraint` (`AppParCurves_Constraint.hxx:20-31`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Constraint {
    /// `AppParCurves_NoConstraint`.
    NoConstraint,
    /// `AppParCurves_PassPoint`.
    PassPoint,
    /// `AppParCurves_TangencyPoint`.
    TangencyPoint,
    /// `AppParCurves_CurvaturePoint`.
    CurvaturePoint,
}

/// `AppParCurves_ConstraintCouple` (`AppParCurves_ConstraintCouple.hxx:30-56`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConstraintCouple {
    /// `myIndex` (`hxx:53`).
    index: i32,
    /// `myConstraint` (`hxx:54`).
    constraint: Constraint,
}

impl Default for ConstraintCouple {
    /// `AppParCurves_ConstraintCouple()` (`cxx:16-22`): index `-1`,
    /// `AppParCurves_NoConstraint`.
    fn default() -> Self {
        Self {
            index: -1,
            constraint: Constraint::NoConstraint,
        }
    }
}

impl ConstraintCouple {
    /// `AppParCurves_ConstraintCouple()` (`cxx:16-22`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `AppParCurves_ConstraintCouple(TheIndex, Cons)` (`cxx:24-29`).
    pub fn from_parts(the_index: i32, cons: Constraint) -> Self {
        Self {
            index: the_index,
            constraint: cons,
        }
    }

    /// `Index()` (`cxx:31-34`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `Constraint()` (`cxx:36-39`).
    pub fn constraint(&self) -> Constraint {
        self.constraint
    }

    /// `SetIndex(TheIndex)` (`cxx:41-44`).
    pub fn set_index(&mut self, the_index: i32) {
        self.index = the_index;
    }

    /// `SetConstraint(Cons)` (`cxx:46-49`).
    pub fn set_constraint(&mut self, cons: Constraint) {
        self.constraint = cons;
    }
}
