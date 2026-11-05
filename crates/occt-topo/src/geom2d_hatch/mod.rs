//! `Geom2dHatch` package: 2D hatching of a face region by curve elements.
//!
//! Sources (OCCT 8.0.0): `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/`
//! - `Geom2dHatch_Element` / `Geom2dHatch_Elements` -> [`element`]
//! - `Geom2dHatch_Hatching` -> [`hatching`]
//! - `Geom2dHatch_Intersector` -> [`intersector`]
//! - `Geom2dHatch_FClass2dOfClassifier` / `Geom2dHatch_Classifier` -> [`fclass2d`]
//! - `Geom2dHatch_Hatcher` -> [`hatcher`] + [`hatcher_domains`]
//!
//! The package is the 2D hatching engine `ChFi3d_Builder::SplitKPart`
//! (`ChFi3d_Builder_SpKP.cxx`) drives through `Geom2dHatch_Hatcher`. It builds
//! on [`crate::hatch_gen`] (`HatchGen_*`), on
//! [`crate::top_trans_curve_transition`] for the global transitions and on
//! `occt_geom2d::geom2d_int` for the curve/curve intersections.
//!
//! `Geom2dHatch_Element` indexes its curve by an `Arc<dyn Curve2d>`: the port
//! has no `Geom2dAdaptor_Curve` class, `Curve2d` already carries the `D0` /
//! `D1` / `FirstParameter` / `LastParameter` surface (see [`element`]).

mod element;
mod fclass2d;
mod hatcher;
mod hatcher_domains;
mod hatching;
mod intersector;

pub use element::{Element, Elements};
pub use fclass2d::{Classifier, FClass2dOfClassifier};
pub use hatcher::Hatcher;
pub use hatching::Hatching;
pub use intersector::Intersector;
