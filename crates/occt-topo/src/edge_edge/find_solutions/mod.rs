//! `IntTools_EdgeEdge` general-curve branch — the parameter-box recursion.
//!
//! Faithful port of the OCCT members and file-local helpers that
//! `IntTools_EdgeEdge::Perform` uses for every pair that is not line/line
//! (`IntTools_EdgeEdge.cxx`):
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | `Prepare` (res coeff / resolution / PTol / type swap) | `:89-181` | [`EdgeEdge::prepare`] (`edge_edge.rs`) |
//! | `FindSolutions(ranges1, ranges2, bSplit2)` | `:290-349` | [`EdgeEdge::find_solutions`] |
//! | `FindSolutions(R1, B1, R2, B2, …)` (recursion) | `:353-549` | [`EdgeEdge::find_solutions_range`] |
//! | `FindParameters` | `:553-671` | [`find_parameters`] |
//! | `MergeSolutions` | `:675-779` | [`EdgeEdge::merge_solutions_boxes`] |
//! | `AddSolution` | `:780-825` | [`EdgeEdge::add_solution_occt`] |
//! | `FindBestSolution` | `:826-901` | [`EdgeEdge::find_best_solution_ranges`] |
//! | `IsIntersection` | `:1060-1146` | [`EdgeEdge::is_intersection`] |
//! | `CheckCoincidence` | `:1150-1206` | [`EdgeEdge::check_coincidence`] |
//! | `FindDistPC` / `DistPC` | `:1210-1362` | [`find_dist_pc`] / [`dist_pc`] |
//! | `SplitRangeOnSegments` | `:1366-1406` | [`split_range_on_segments`] |
//! | `BndBuildBox` | `:1410-1419` | [`bnd_build_box`] |
//! | `PointBoxDistance` | `:1423-1452` | [`point_box_distance`] |
//! | `TypeToInteger` | `:1456-1482` | [`curve_kind_to_integer`] |
//! | `ResolutionCoeff` | `:1486-1559` | [`resolution_coeff`] |
//! | `Resolution` | `:1561-1607` | [`resolution_of`] |
//! | `CurveDeflection` | `:1611-1638` | [`curve_deflection`] |
//! | `IsClosed` | `:1642-1659` | [`is_closed`] |
//!
//! Result mapping (the port's [`CommonPrt`] keeps less than
//! `IntTools_CommonPrt`): a `TopAbs_EDGE` solution becomes a
//! [`CommonPartType::Edge`] entry whose `range` is the OCCT `Range1` (already
//! expressed on the *caller's* first edge, see [`EdgeEdge::add_solution_occt`]),
//! and a `TopAbs_VERTEX` solution becomes one [`PntOn2Faces`] whose `uv1`/`uv2`
//! are the `FindBestSolution` parameters on the caller's first/second edge.
//! `IntTools_CommonPrt::{Range2, Edge1, Edge2, BoundingPoints}` are not stored
//! by that reduced struct (pre-existing data-model reduction; no consumer reads
//! them today).
//!
//! The B-rep adaptor range (`myCurve1`/`myCurve2`) is the *edge* range, so the
//! port calls `D0`/`D1` on the edge curve directly with edge parameters —
//! `IntTools_EdgeEdge.cxx:94-95,164-165` initializes `BRepAdaptor_Curve` on the
//! edge, whose `FirstParameter`/`LastParameter` are the edge's.

mod boxes;
mod curve_kind;
mod dist_pc;
mod edge_edge_members;
mod parameters;
mod resolution;

pub(crate) use self::curve_kind::{curve_kind_to_integer, geom_curve_kind, GeomCurveKind};
pub(in crate::edge_edge) use self::resolution::{curve_deflection, resolution_coeff, resolution_of};
