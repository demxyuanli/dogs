//! `GeomLProp` — local properties of curves and surfaces.
//!
//! Source: `ModelingData/TKGeomBase/GeomLProp` and `ModelingData/TKGeomBase/LProp`
//! (OCCT 8.0.0).
//!
//! Ported: [`cl_props`] — `GeomLProp_CLProps` (the 3D instantiation of
//! `GeomLProp_CLPropsBase`, `GeomLProp_CLProps.hxx:47-219`) with the
//! `LProp_CurveUtils::DirectAccess` policy (`LProp_CurveUtils.hxx:53-93`).
//! The 2D instantiation lives in `occt-geom2d::cl_props2d`.

pub mod cl_props;
pub mod sl_props;

pub use cl_props::{ClProps, LPropStatus};
pub use sl_props::SlProps;
