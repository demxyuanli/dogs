//! `IntCurve` package (TKGeomAlgo, `ModelingAlgorithms/TKGeomAlgo/IntCurve/`).
//!
//! This module holds the parts of `IntCurve` that are expressible against `gp`
//! and `ElCLib` alone, i.e. everything that does not need a `Geom2d` curve tool
//! and does not need the `math_FunctionAllRoots` root finder. Nothing here is
//! wired into a caller yet; it is a zero-wiring tier.
//!
//! ## Ported
//!
//! * [`pconic`] — `IntCurve_PConic` (`IntCurve_PConic.hxx`, `_0.cxx`, `.lxx`)
//! * [`iconic_tool`] — `IntCurve_IConicTool` (`IntCurve_IConicTool.hxx`, `.cxx`)
//! * [`pconic_tool`] — `IntCurve_PConicTool` (`IntCurve_PConicTool.hxx`, `.cxx`)
//! * [`project_on_pconic_tool`] — `IntCurve_ProjectOnPConicTool`
//!   (`IntCurve_ProjectOnPConicTool.hxx`, `.cxx`)
//!
//! The four 2D `ElCLib` overload sets they need (`Value`, `D1`, `D2`, `D3`,
//! `DN`, `Parameter` over `gp_Ax2d` / `gp_Ax22d`) live in
//! [`crate::elib::clib2d`].
//!
//! ## UNPORTED
//!
//! ### Blocked on the `math` root finders
//!
//! `IntImpParGen_Intersector.gxx:28-29` includes `math_FunctionSample.hxx` and
//! `math_FunctionAllRoots.hxx`; `IntCurve_PConic.hxx:55-56` refers to
//! `FunctionAllRoots` for the meaning of `EpsX`. So the whole
//! implicit-vs-parametric family is gated on TKMath code that is not in this
//! tree:
//!
//! * `math_FunctionWithDerivative.hxx` (2 107 B), `math_FunctionSample.hxx`
//!   (1 739 B), `math_FunctionSample.cxx` (1 348 B)
//! * `math_FunctionAllRoots.hxx` / `.lxx` / `.cxx` (4 097 + 1 895 + 7 478 B)
//! * `math_FunctionRoots.hxx` / `.lxx` / `.cxx` (3 368 + 1 484 + 38 347 B)
//! * `math_FunctionSetRoot.hxx` / `.cxx` (7 730 + 42 067 B) and
//!   `math_FunctionSetWithDerivatives.hxx` (2 298 B), which `math_FunctionRoots`
//!   uses
//!
//! Blocked on those:
//!
//! * `IntImpParGen_Intersector.gxx` (31 348 B) — the generic
//!   implicit/parametric intersector
//! * `IntCurve_IntImpConicParConic.hxx` (4 177 B) / `_0.cxx` (1 791 B)
//! * `IntCurve_MyImpParToolOfIntImpConicParConic.hxx` / `_0.cxx`
//! * `IntImpParGen_Tool.hxx` / `.cxx` (1 597 + 5 379 B) and
//!   `IntImpParGen.hxx` / `.cxx` (3 199 + 7 003 B). These two files define
//!   `NormalizeOnDomain` / `Determine_Position` / `Determine_Transition` twice,
//!   with *different* bodies: `IntImpParGen.cxx:28-252` compares curvatures
//!   against `TOLERANCE_ANGULAIRE` (1e-8) and normalizes the period with a
//!   two-condition `while`, whereas `IntImpParGen_Tool.cxx:24-203` compares
//!   against `gp::Resolution()` and uses a one-condition `while`. Both must be
//!   ported separately, not merged.
//! * `IntCurve_IntConicConic.hxx` / `.lxx` / `.cxx` (15 464 + 6 982 + 37 693 B)
//!   / `IntCurve_IntConicConic_1.cxx` (106 021 B). Eleven of the `Perform`
//!   overloads call `IntCurve_IntImpConicParConic::Perform` (see
//!   `IntCurve_IntConicConic.cxx:172`, `:207`, `:213`, `:291`, `:325`, ...);
//!   `_1.cxx:807` (circle/circle), `:1381` (line/line), `:2236` (line/circle)
//!   and `:2861` (line/ellipse) instead work directly against
//!   `IntAna2d_AnaIntersection`, which *is* ported ([`crate::intana2d`]). The
//!   `_1.cxx:41-806` block is a `TOLERANCE_ANGULAIRE` macro/helper prologue.
//! * `IntCurve_IntConicConic_Tool.hxx` / `.cxx` (4 188 + 7 034 B)
//! * `IntCurve_UserIntConicCurveGen.gxx` (30 379 B), `IntCurve_IntConicCurveGen`
//!   (`.gxx` 4 159 B + `.lxx`)
//!
//! ### Blocked on a `Geom2d` curve tool
//!
//! * `IntCurve_IntCurveCurveGen.gxx` (31 189 B) / `.lxx` (6 785 B) — the generic
//!   curve/curve intersector that `Geom2dInt_GInter` instantiates. Its
//!   `TheCurveTool` is `Geom2dInt_Geom2dCurveTool`, which needs a `Geom2d`
//!   curve, so it cannot live in `occt-core`.
//! * `IntCurve_Polygon2dGen.gxx` (11 507 B) / `.lxx` (2 643 B)
//! * `IntCurve_IntPolyPolyGen.gxx` (63 103 B)
//! * `IntCurve_DistBetweenPCurvesGen.gxx` (3 088 B)
//! * `IntCurve_ExactIntersectionPoint.gxx` (9 167 B)

pub mod iconic_tool;
pub mod pconic;
pub mod pconic_tool;
pub mod project_on_pconic_tool;

pub use iconic_tool::IntCurveIConicTool;
pub use pconic::IntCurvePConic;
