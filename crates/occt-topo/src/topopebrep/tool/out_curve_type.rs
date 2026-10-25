//! `TopOpeBRepTool_OutCurveType` — the kind of 3d curve a `MakeCurves` call
//! produces.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_OutCurveType.hxx:20-25`.

/// `enum TopOpeBRepTool_OutCurveType` (`hxx:20-25`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutCurveType {
    /// `TopOpeBRepTool_BSPLINE1` (`hxx:22`): a degree-1 BSpline through the
    /// walking points.
    #[default]
    BSpline1,
    /// `TopOpeBRepTool_APPROX` (`hxx:23`): an approximation curve on the
    /// walking points.
    Approx,
    /// `TopOpeBRepTool_INTERPOL` (`hxx:24`): an interpolating curve.
    ///
    /// `TopOpeBRepTool_CurveTool::MakeCurves` never reads this value
    /// (`CurveTool.cxx:313-321` only tests `CompC3D`), and no `switch` in the
    /// package handles it.
    Interpol,
}
