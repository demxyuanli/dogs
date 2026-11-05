//! Port of `BRepFill_Evolved` (`BRepFill_Evolved.cxx:1-3443`,
//! `BRepFill_Evolved.hxx:1-200`).
//!
//! The class sweeps a profile wire along a planar spine face (or a spine
//! wire) and produces the same three families of elementary surfaces the
//! reference builds: tubes around the spine, revolutions around vertical
//! axes and prisms for the planar / vertical profile pieces.
//!
//! The port keeps the reference's decomposition one module per OCCT function
//! so the control flow can be audited against the `.cxx` line by line:
//!
//! - [`side`]: the static predicates (`Side`, `IsVertical`, `IsPlanar`,
//!   `Bubble`, `PosOnFace`, `DoubleOrNotInFace`, `Compare`, `IsInversed`,
//!   `ConcaveSide`, `Relative`, `SimpleExpression`).
//! - [`shape_map`]: `NCollection_DataMap<TopoDS_Shape, ...>` in insertion
//!   order, the iteration order the whole class relies on.
//! - [`prepare_profile`]: `PrepareProfile` (`cxx:1553-1711`) and
//!   `CutEdgeProf` (`cxx:3122-3266`).
//!
//! - [`prepare_spine`]: `PrepareSpine` (`cxx:1713-1830`) and `CutEdge`
//!   (`cxx:3276-3381`).
//! - [`trim_edge`]: `TrimEdge` (`cxx:2770-2866`).
//! - [`trim_face`]: `TrimFace` (`cxx:2659-2711`) and `AddDegeneratedEdge`
//!   (`cxx:3000-3120`).
//! - [`intervals`]: `ComputeIntervals` (`cxx:2866-2946`) and `VertexFromNode`
//!   (`cxx:3391-3443`).
//! - [`put_profil_at`]: `PutProfilAt` (`cxx:2714-2768`).
//! - [`elementary`]: `ElementaryPerform` (`cxx:640-1264`), `MakePipe`
//!   (`cxx:2303-2400`) and `MakeRevol` (`cxx:2266-2300`).
//! - [`planar`]: `PlanarPerform` (`cxx:1267-1412`).
//! - [`vertical`]: `VerticalPerform` (`cxx:1418-1520`).
//! - [`glue_add`]: `Add` (`cxx:1844-1957`).
//! - [`transfert`]: `Transfert` (`cxx:1966-2053`).
//! - [`top_bottom`]: `AddTopAndBottom` (`cxx:2078-2225`) and `MakeSolid`
//!   (`cxx:2229-2262`).
//! - [`continuity`]: `ContinuityOnOffsetEdge` (`cxx:2452-2543`).
//! - [`perform`]: `Perform` / `PrivatePerform` / `SetWork` / `FindLocation` /
//!   `TransformInitWork` (`cxx:309-331`, `:334-547`, `:614-618`, `:2404-2437`,
//!   `:2441-2445`).
//! - [`evolved`]: the `BRepFill_Evolved` data holders (`BRepFill_Evolved.hxx:32-200`).
//!
//! No part of `BRepFill_Evolved.cxx` is left unported. Two prerequisites keep
//! the documented deviations recorded at their call sites: the port stores no
//! per-representation location (`BRep_GCurve::Location`), and `BRepTool::Locked`
//! is not carried, so the locked-shape throws of `BRep_Builder` are absent.

pub mod continuity;
pub mod elementary;
pub mod evolved;
pub mod glue_add;
pub mod intervals;
pub mod make_pipe;
pub mod make_revol;
pub mod perform;
pub mod planar;
pub mod prepare_profile;
pub mod prepare_spine;
pub mod put_profil_at;
pub mod shape_map;
pub mod side;
pub mod top_bottom;
pub mod transfert;
pub mod trim_edge;
pub mod trim_face;
pub mod vertical;

pub use evolved::BRepFillEvolved;
pub use intervals::{compute_intervals, vertex_from_node, Intervals, MapNodeVertex};
pub use prepare_profile::{cut_edge_prof, oy2d_line, prepare_profile};
pub use prepare_spine::{edge_pcurve, prepare_spine};
pub use put_profil_at::put_profil_at;
pub use shape_map::{ShapeListMap, ShapeMap};
pub use trim_edge::trim_edge;
pub use trim_face::{add_degenerated_edge, trim_face};
pub use side::{
    altitud, bubble, compare, complement, concave_side, distance_to_oz, double_or_not_in_face,
    edge_vertices_oriented, is_inversed, is_planar, is_vertical, pos_on_face, relative, side,
    simple_expression, yoz_plane, FILL_CONFUSION,
};
