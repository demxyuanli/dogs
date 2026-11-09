// `ShapeFix_Wire::FixEdgeCurves`' `FixSameParameter` arm
// (`ShapeFix_Wire.cxx:953-1005`, including the `CheckPCurveRange` re-projection
// it owns at `cxx:958-993`) is gated by `NeedFix(myFixSameParameterMode)`. That
// field is fed from the *separate* key `FixEdgeSameParameterMode`
// (`ShapeProcess_OperLibrary.cxx:857`), which is `0` on both read paths
// (`DE_ShapeFixParameters.hxx:75` `FixMode::NotFix`, `DEIGES_Parameters.cxx:64-72`,
// `DESTEP_Parameters.cxx:90-93`, `XSAlgo_ShapeProcessor.cxx:681`; the legacy
// maps `IGESData.cxx:152` and `STEPControl_Controller.cxx:247` agree), so
// `NeedFix(0)` is false. The gate now lives in `WireFixModes`
// (`fix_same_parameter_mode`) and the presets carry that value.
//
// The pass that covers the same `ShapeFix_Edge::FixSameParameter` work is the
// whole-shape `ShapeFix_Shape.cxx:259` -> `ShapeFix::SameParameter(Result,
// false, 0.0)` (`ShapeFix.cxx:70-258`), gated by the distinct key
// `FixSameParameterMode` (`ShapeProcess_OperLibrary.cxx:813`, DE default
// `FixOrNot = -1`), which does run.

// ---------------------------------------------------------------------------
// ShapeFix_Wire::FixNotchedEdges / ShapeAnalysis_Wire::CheckNotchedEdges
// ---------------------------------------------------------------------------

mod statics;
mod check_pcurves;
mod perform;

pub use check_pcurves::*;
pub use perform::*;
