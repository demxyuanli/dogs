//! `Approx_*` — the general approximation interfaces.
//!
//! Ported: `Approx_SweepFunction` (`Approx_SweepFunction.hxx/.cxx`),
//! `Approx_Status`, `Approx_ParametrizationType`.
//! `Approx_SweepApproximation` (the driver) is not ported yet.

mod sweep_function;
mod types;

pub use sweep_function::{ApproxSweepFunction, SweepD1, SweepD2};
pub use types::{ApproxParametrizationType, ApproxStatus};
