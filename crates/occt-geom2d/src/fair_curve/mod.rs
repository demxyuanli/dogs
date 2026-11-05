//! `FairCurve` package (`src/ModelingAlgorithms/TKGeomAlgo/FairCurve`): the
//! energy-minimising fair-curve solvers used by the fillet corner builder
//! (`ChFi3d_Builder_CnCrn.cxx:585-676`, `CalculBatten`).
//!
//! Landed: `FairCurve_AnalysisCode`, `FairCurve_BattenLaw`,
//! `FairCurve_DistributionOfEnergy`, `FairCurve_DistributionOfSagging`,
//! `FairCurve_DistributionOfTension`, `FairCurve_Energy`,
//! `FairCurve_EnergyOfBatten`, `FairCurve_Newton`, `FairCurve_Batten`.
//!
//! UNPORTED (no Rust counterpart yet, do not stub):
//! - `FairCurve_DistributionOfJerk` (`FairCurve_DistributionOfJerk.cxx:1-345`)
//! - `FairCurve_EnergyOfMVC` (`FairCurve_EnergyOfMVC.cxx:1-186`)
//! - `FairCurve_MinimalVariation` (`FairCurve_MinimalVariation.cxx:1-598`)
//! - `FairCurve_Batten::Dump` (`FairCurve_Batten.cxx:633-691`), a debug stream
//!   helper with no Rust analogue.

pub mod analysis_code;
pub mod batten;
pub mod batten_law;
pub mod distribution_of_energy;
pub mod distribution_of_sagging;
pub mod distribution_of_tension;
pub mod energy;
pub mod energy_of_batten;
pub mod newton;

pub use analysis_code::FairCurveAnalysisCode;
pub use batten::Batten;
pub use batten_law::BattenLaw;
pub use distribution_of_energy::SharedPoles;
pub use distribution_of_sagging::DistributionOfSagging;
pub use distribution_of_tension::DistributionOfTension;
pub use energy::Energy;
pub use energy_of_batten::EnergyOfBatten;
pub use newton::FairCurveNewton;
