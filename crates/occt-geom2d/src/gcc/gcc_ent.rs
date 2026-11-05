//! `GccEnt` package (TKG3d `GccEnt`): the position enum every GCC solver
//! qualifies its arguments with, plus the two qualified 2D primitives.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccEnt_Position` | `GccEnt_Position.hxx:16-38` |
//! | `GccEnt_QualifiedCirc` | `GccEnt_QualifiedCirc.cxx:27-53` |
//! | `GccEnt_QualifiedLin` | `GccEnt_QualifiedLin.cxx:27-47` |

use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::lin2d::GpLin2d;

/// `enum GccEnt_Position` (`GccEnt_Position.hxx:16-38`).
///
/// * `Unqualified` - the position of the solution is undefined with respect to
///   the argument,
/// * `Enclosing` - the solution encompasses the argument,
/// * `Enclosed` - the solution is encompassed by the argument,
/// * `Outside` - the solution and the argument are external to one another,
/// * `NoQualifier` - the value returned when consulting the qualifier of an
///   argument defined as `Unqualified`.
///
/// The interior of a line is the left-hand side of the line in relation to its
/// orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GccEntPosition {
    /// `GccEnt_unqualified`
    #[default]
    Unqualified,
    /// `GccEnt_enclosing`
    Enclosing,
    /// `GccEnt_enclosed`
    Enclosed,
    /// `GccEnt_outside`
    Outside,
    /// `GccEnt_noqualifier`
    NoQualifier,
}

impl GccEntPosition {
    /// `TheQualifier == GccEnt_unqualified`.
    pub fn is_unqualified(self) -> bool {
        self == GccEntPosition::Unqualified
    }

    /// `TheQualifier == GccEnt_enclosing`.
    pub fn is_enclosing(self) -> bool {
        self == GccEntPosition::Enclosing
    }

    /// `TheQualifier == GccEnt_enclosed`.
    pub fn is_enclosed(self) -> bool {
        self == GccEntPosition::Enclosed
    }

    /// `TheQualifier == GccEnt_outside`.
    pub fn is_outside(self) -> bool {
        self == GccEntPosition::Outside
    }
}

/// `GccEnt_QualifiedCirc`: a `gp_Circ2d` plus the position the requested
/// solution must take with respect to it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GccEntQualifiedCirc {
    the_qualified: GpCirc2d,
    the_qualifier: GccEntPosition,
}

impl GccEntQualifiedCirc {
    /// `GccEnt_QualifiedCirc(Qualified, Qualifier)`
    /// (`GccEnt_QualifiedCirc.cxx:27-31`).
    pub fn new(qualified: GpCirc2d, qualifier: GccEntPosition) -> Self {
        Self { the_qualified: qualified, the_qualifier: qualifier }
    }

    /// `Qualified()` (`cxx:33-36`).
    pub fn qualified(&self) -> GpCirc2d {
        self.the_qualified
    }

    /// `Qualifier()` (`cxx:38-41`).
    pub fn qualifier(&self) -> GccEntPosition {
        self.the_qualifier
    }

    /// `IsUnqualified()` (`cxx:43-46`).
    pub fn is_unqualified(&self) -> bool {
        self.the_qualifier == GccEntPosition::Unqualified
    }

    /// `IsEnclosing()` (`cxx:48-51`).
    pub fn is_enclosing(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosing
    }

    /// `IsEnclosed()` (`cxx:53-56`).
    pub fn is_enclosed(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosed
    }

    /// `IsOutside()` (`cxx:58-61`).
    pub fn is_outside(&self) -> bool {
        self.the_qualifier == GccEntPosition::Outside
    }
}

/// `GccEnt_QualifiedLin`: a `gp_Lin2d` plus the position the requested
/// solution must take with respect to it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GccEntQualifiedLin {
    the_qualified: GpLin2d,
    the_qualifier: GccEntPosition,
}

impl GccEntQualifiedLin {
    /// `GccEnt_QualifiedLin(Qualified, Qualifier)`
    /// (`GccEnt_QualifiedLin.cxx:27-31`).
    pub fn new(qualified: GpLin2d, qualifier: GccEntPosition) -> Self {
        Self { the_qualified: qualified, the_qualifier: qualifier }
    }

    /// `Qualified()` (`cxx:33-36`).
    pub fn qualified(&self) -> GpLin2d {
        self.the_qualified
    }

    /// `Qualifier()` (`cxx:38-41`).
    pub fn qualifier(&self) -> GccEntPosition {
        self.the_qualifier
    }

    /// `IsUnqualified()` (`cxx:43-46`).
    pub fn is_unqualified(&self) -> bool {
        self.the_qualifier == GccEntPosition::Unqualified
    }

    /// `IsEnclosed()` (`cxx:48-51`).
    pub fn is_enclosed(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosed
    }

    /// `IsOutside()` (`cxx:53-56`).
    pub fn is_outside(&self) -> bool {
        self.the_qualifier == GccEntPosition::Outside
    }

    /// `IsEnclosing()`: `GccEnt_QualifiedLin` has no such query in
    /// `GccEnt_QualifiedLin.hxx`; the two solvers that need it
    /// (`GccAna_Circ2d2TanRad_1.cxx:63-75`, `_4.cxx:67-71`) reject an enclosing
    /// qualification as `GccEnt_BadQualifier` instead, which is what the
    /// dispatch in [`super::circ2d_2tan_rad`] does.
    pub fn is_enclosing(&self) -> bool {
        self.the_qualifier == GccEntPosition::Enclosing
    }
}
