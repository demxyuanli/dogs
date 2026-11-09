//! `Draft_Modification::Perform` (`Draft_Modification_1.cxx:708-1777`): the
//! five phases that turn the records `InternalAdd` / `Propagate` filled into
//! the rebuilt surfaces, curves, pcurves, points and parameters.
//!
//! | OCCT | lines | file |
//! |---|---|---|
//! | `Perform` driver (throw, `myComp`, phase order) | `_1.cxx:708-725` | this module |
//! | "Calculate eventual faces" | `_1.cxx:724-824` | [`faces`] |
//! | "Calculate new edges" | `_1.cxx:826-1445` | [`edges`], [`edges_intersection`], [`edges_spline`] |
//! | "Calculate new vertices" + `Choose` | `_1.cxx:1447-1638`, `:2047-2177` | [`vertices`] |
//! | "small loop of validation/protection" | `_1.cxx:1640-1774` | [`validation`] |
//! | `Parameter` / `SmartParameter` | `_1.cxx:2180-2379` | [`parameters`] |
//!
//! The phase bodies are the `.cxx` control flow; this module only reproduces
//! the driver: the `Standard_ConstructionError` for an already-pending error,
//! the one-shot `myComp` latch, the `Propagate()` guard and the fact that a
//! phase that records `errStat` / `badShape` simply returns - `Perform` itself
//! is `void` and the caller reads `IsDone()` / `Error()`.

mod edges;
mod edges_intersection;
mod edges_spline;
mod faces;
mod parameters;
mod validation;
mod vertices;

use super::DraftModification;

impl DraftModification {
    /// `Draft_Modification::Perform()` (`Draft_Modification_1.cxx:708-1777`).
    ///
    /// `Err` is the reference's `Standard_ConstructionError` thrown while an
    /// error from `Add` is still pending (`cxx:710-713`); a phase failure is
    /// recorded in `err_stat` / `bad_shape` exactly as in OCCT, so `IsDone()`
    /// and `Error()` answer for it.
    pub fn perform(&mut self) -> Result<(), String> {
        // `cxx:710-713`
        if self.bad_shape.is_some() {
            return Err("Standard_ConstructionError: Draft_Modification::Perform".to_string());
        }

        // `cxx:715-720`
        if !self.my_comp {
            self.my_comp = true;
            if !self.propagate() {
                return Ok(());
            }

            // `cxx:724-1637`: every phase stops the whole run on failure.
            if self.perform_faces().is_err() {
                return Ok(());
            }
            if self.perform_edges().is_err() {
                return Ok(());
            }
            if self.perform_vertices().is_err() {
                return Ok(());
            }
            if self.validate_protection().is_err() {
                return Ok(());
            }
        }
        Ok(())
    }
}
