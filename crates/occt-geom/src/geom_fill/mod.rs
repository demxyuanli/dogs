//! `GeomFill` — the transfinite fillings (`GeomFill_Filling`), the section
//! profiler, and the surface builders bounded by 2/3/4 curves.
//!
//! Source: `TKGeomAlgo/GeomFill` (OCCT 8.0.0).
//!
//! Layout follows the OCCT classes:
//! - [`filling`] — `GeomFill_Filling` (the abstract pole grid),
//! - [`stretch`] — `GeomFill_Stretch`,
//! - [`coons`] — `GeomFill_Coons`,
//! - [`curved`] — `GeomFill_Curved`,
//! - [`bspline_curves`] — `GeomFill_BSplineCurves`,
//! - [`bezier_curves`] — `GeomFill_BezierCurves`,
//! - [`constrained_filling`] — `GeomFill_ConstrainedFilling` and its
//!   boundaries / tangency fields,
//! - [`profiler`] — `GeomFill_Profiler`,
//! - [`generator`] — `GeomFill_Generator`,
//! - [`param_trim`] — the `Geom_TrimmedCurve` handle `BRepFill` builds.

pub mod app_surf;
pub mod bezier_curves;
pub mod bound_with_surf;
pub mod boundary;
pub mod bspline_curves;
pub mod circular_blend_func;
pub mod coons;
pub mod coons_alg_patch;
pub mod constrained_filling;
pub mod constant_bi_normal;
pub mod corner_state;
pub mod corrected_frenet;
pub mod curve_and_trihedron;
pub mod darboux;
pub mod discrete_trihedron;
pub mod draft_trihedron;
pub mod evolved_section;
pub mod curved;
pub mod filling;
pub mod fixed;
pub mod frenet;
pub mod function_draft;
pub mod function_guide;
pub mod generator;
pub mod gordon;
pub mod guide_trihedron_ac;
pub mod guide_trihedron_plan;
pub mod line;
pub mod location_guide;
pub mod location_draft;
pub mod location_law;
pub mod loc_function;
pub mod n_sections;
pub mod param_trim;
pub mod pipe;
pub mod pipe_error;
pub mod plane_curve_inter;
pub mod plan_func;
pub mod polynomial_convertor;
pub mod profiler;
pub mod quasi_angular_convertor;
pub mod section_generator;
pub mod section_convert;
pub mod section_law;
pub mod section_placement;
pub mod snglr_func;
pub mod stretch;
pub mod sweep;
pub mod sweep_function;
pub mod sweep_section_generator;
pub mod tensor;
pub mod tgt_field;
pub mod trihedron;
pub mod trihedron_with_guide;
pub mod trihedron_law;
pub mod uniform_section;
pub mod utilities;

pub use app_surf::{AppBlendSectionGenerator, AppSurf, AppSweep, SectionValues};
pub use bezier_curves::BezierCurves;
pub use bound_with_surf::BoundWithSurf;
pub use boundary::{Boundary, DegeneratedBound, SimpleBound};
pub use bspline_curves::{BSplineCurves, FillingStyle};
pub use circular_blend_func::CircularBlendFunc;
pub use coons::Coons;
pub use coons_alg_patch::CoonsAlgPatch;
pub use constrained_filling::ConstrainedFilling;
pub use constant_bi_normal::ConstantBiNormal;
pub use corner_state::CornerState;
pub use corrected_frenet::CorrectedFrenet;
pub use curve_and_trihedron::CurveAndTrihedron;
pub use darboux::Darboux;
pub use discrete_trihedron::DiscreteTrihedron;
pub use draft_trihedron::DraftTrihedron;
pub use evolved_section::EvolvedSection;
pub use curved::Curved;
pub use filling::Filling;
pub use fixed::Fixed;
pub use frenet::Frenet;
pub use function_draft::FunctionDraft;
pub use function_guide::FunctionGuide;
pub use generator::generator_surface;
pub use gordon::{
    ApproximationMode, BuildReport, BuildStage, Gordon, NetworkSurface, ResultStatus,
};
pub use guide_trihedron_ac::GuideTrihedronAC;
pub use guide_trihedron_plan::GuideTrihedronPlan;
pub use trihedron_with_guide::{TrihedronWithGuide, TrihedronWithGuideState};
pub use line::Line;
pub use location_guide::LocationGuide;
pub use location_draft::LocationDraft;
pub use location_law::{LocationD1, LocationD2, LocationLaw};
pub use loc_function::LocFunction;
pub use n_sections::NSections;
pub use param_trim::ParameterTrim;
pub use pipe::Pipe;
pub use pipe_error::GeomFillPipeError;
pub use plan_func::PlanFunc;
pub use polynomial_convertor::PolynomialConvertor;
pub use profiler::Profiler;
pub use quasi_angular_convertor::QuasiAngularConvertor;
pub use section_generator::SectionGenerator;
pub use section_law::SectionLaw;
pub use section_placement::{SectionGeometry, SectionPlacement};
pub use snglr_func::SnglrFunc;
pub use sweep_section_generator::SweepSectionGenerator;
pub use stretch::Stretch;
pub use sweep::{ApproxStyle, Sweep};
pub use tensor::Tensor;
pub use tgt_field::{TgtField, TgtOnCoons};
pub use trihedron::GeomFillTrihedron;
pub use trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};
pub use uniform_section::UniformSection;
pub use utilities::GeomFill;

