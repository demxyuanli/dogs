//! `IntImpParGen_Intersector` (`IntImpParGen_Intersector.gxx:42-822`): the
//! generic intersector between an implicit curve (through an `ImpTool`) and a
//! parametric curve (`ParCurve`), driven by `math_FunctionAllRoots`.
//!
//! The OCCT class is a template over `ImpTool` / `ParCurve` / `ParTool` /
//! `MyImpParTool` / `ProjectOnPCurveTool` and inherits from
//! `IntRes2d_Intersection`. The only instantiation used by `Geom2dInt_GInter`
//! is `IntCurve_IntImpConicParConic`, i.e.
//! `ImpTool = IntCurve_IConicTool`, `ParCurve = IntCurve_PConic`,
//! `ParTool = IntCurve_PConicTool`,
//! `ProjectOnPCurveTool = IntCurve_ProjectOnPConicTool`. This port is
//! monomorphised onto exactly that instantiation instead of introducing traits,
//! which is what the C++ template expansion does anyway.

mod domain_scan;
mod imp_par_tool;
mod param_search;
mod perform;
mod structure;

pub use imp_par_tool::MyImpParTool;
pub use structure::IntImpParGenIntersector;
