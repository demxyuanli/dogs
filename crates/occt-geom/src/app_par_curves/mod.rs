//! `AppParCurves` — the pole containers and basis builders shared by the
//! parallel-curves approximation used by `GeomFill_NSections` and
//! `GeomFill_AppSurf`.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves` (OCCT 8.0.0).
//!
//! Layout follows the OCCT classes:
//! - [`multi_point`] — `AppParCurves_MultiPoint`,
//! - [`multi_curve`] — `AppParCurves_MultiCurve`,
//! - [`multi_bspline_curve`] — `AppParCurves_MultiBSpCurve`,
//! - [`least_square`] — `AppParCurves_LeastSquare`,
//! - [`resol_constraint`] — `AppParCurves_ResolConstraint`,
//! - [`function`] — `AppParCurves_Function`, the parametrisation objective,
//! - [`constraint`] — `AppParCurves_Constraint`,
//!   `AppParCurves_ConstraintCouple`,
//! - [`bernstein`] — `AppParCurves::BernsteinMatrix` / `Bernstein` /
//!   `SecondDerivativeBernstein`,
//! - [`spline_function`] — `AppParCurves::SplineFunction`,
//! - [`gradient`] — `AppParCurves_Gradient`,
//! - [`bsp_gradient`] — `AppParCurves_BSpGradient`.
//!
//! Still UNPORTED from the same package:
//! - `AppParCurves_Projection` (`AppParCurves_Projection.gxx`).

pub mod bernstein;
pub mod bsp_function;
pub mod bsp_gradient;
pub mod constraint;
pub mod function;
pub mod gradient;
pub mod least_square;
pub mod multi_bspline_curve;
pub mod multi_curve;
pub mod multi_point;
pub mod resol_constraint;
pub mod spline_function;

pub use bernstein::AppParCurves;
pub use bsp_function::BspFunction;
pub use bsp_gradient::BspGradient;
pub use constraint::{Constraint, ConstraintCouple};
pub use function::ParFunction;
pub use gradient::Gradient;
pub use least_square::LeastSquare;
pub use multi_bspline_curve::MultiBSpCurve;
pub use multi_curve::MultiCurve;
pub use multi_point::MultiPoint;
pub use resol_constraint::ResolConstraint;
