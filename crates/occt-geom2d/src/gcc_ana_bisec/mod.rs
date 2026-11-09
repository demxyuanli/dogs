//! `GccAna` bisecting algorithms (`TKGeomAlgo/GccAna`), the `*2dBisec` family.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccAna_Pnt2dBisec` | `GccAna_Pnt2dBisec.cxx:28-77` |
//! | `GccAna_LinPnt2dBisec` | `GccAna_LinPnt2dBisec.cxx:32-89` |
//! | `GccAna_CircLin2dBisec` | `GccAna_CircLin2dBisec.cxx:35-160` |
//! | `GccAna_CircPnt2dBisec` | `GccAna_CircPnt2dBisec.cxx:36-196` |
//! | `GccAna_Circ2dBisec` | `GccAna_Circ2dBisec.cxx:36-414` |
//! | `GccAna_Lin2dBisec` | `GccAna_Lin2dBisec.cxx:39-194` |
//!
//! Every class answers `GccInt_Bisec` solutions; OCCT reports failure through
//! `StdFail_NotDone` / `GccAna_NoSolution` / `Standard_OutOfRange` exceptions,
//! mirrored here by [`GccAnaBisecError`].

mod circ_circ;mod circ_lin;
mod circ_pnt;
mod lin_lin;
mod lin_pnt;
mod pnt_pnt;

pub use circ_circ::GccAnaCirc2dBisec;
pub use circ_lin::GccAnaCircLin2dBisec;
pub use circ_pnt::GccAnaCircPnt2dBisec;
pub use lin_lin::GccAnaLin2dBisec;
pub use lin_pnt::GccAnaLinPnt2dBisec;
pub use pnt_pnt::GccAnaPnt2dBisec;

/// Exceptions raised by the `GccAna` bisecting algorithms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GccAnaBisecError {
    /// `StdFail_NotDone`
    NotDone,
    /// `GccAna_NoSolution`
    NoSolution,
    /// `Standard_OutOfRange`
    OutOfRange,
}

impl std::fmt::Display for GccAnaBisecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::NotDone => "StdFail_NotDone",
            Self::NoSolution => "GccAna_NoSolution",
            Self::OutOfRange => "Standard_OutOfRange",
        };
        write!(f, "{s}")
    }
}

impl std::error::Error for GccAnaBisecError {}

/// Shared guard of `GccAna_*Bisec::ThisSolution`: `!WellDone` throws
/// `StdFail_NotDone`, an index outside `1..=NbrSol` throws
/// `Standard_OutOfRange`.
fn check_index(done: bool, nbr_sol: i32, index: i32) -> Result<(), GccAnaBisecError> {
    if !done {
        return Err(GccAnaBisecError::NotDone);
    }
    if index <= 0 || index > nbr_sol {
        return Err(GccAnaBisecError::OutOfRange);
    }
    Ok(())
}
