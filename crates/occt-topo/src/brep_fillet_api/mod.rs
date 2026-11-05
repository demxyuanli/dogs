//! `BRepFilletAPI` — the four front ends of `TKFillet`: fillets and chamfers on
//! a planar face ([`MakeFillet2d`]) and on a solid ([`MakeFillet`],
//! [`MakeChamfer`]), over the shared contour / simulation interface
//! [`LocalOperation`].
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepFilletAPI/` (OCCT 8.0.0):
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepFilletAPI_MakeFillet2d.hxx` | 303 | [`MakeFillet2d`] |
//! | `BRepFilletAPI_MakeFillet2d.cxx` | 171 | [`MakeFillet2d`] |
//! | `BRepFilletAPI_MakeFillet2d.lxx` | 71 | [`MakeFillet2d`] |
//! | `BRepFilletAPI_LocalOperation.hxx` | 101 | [`LocalOperation`] |
//! | `BRepFilletAPI_MakeFillet.hxx` | 387 | [`MakeFillet`] |
//! | `BRepFilletAPI_MakeFillet.cxx` | 545 | [`MakeFillet`] |
//! | `BRepFilletAPI_MakeChamfer.hxx` | 299 | [`MakeChamfer`] |
//! | `BRepFilletAPI_MakeChamfer.cxx` | 371 | [`MakeChamfer`] |
//!
//! [`MakeFillet2d`] is a front end over [`crate::chfi2d::ChFi2dBuilder`]
//! (`myMakeChFi2d`); [`MakeFillet`] and [`MakeChamfer`] are front ends over
//! [`crate::chfi3d::builder::fil_builder::ChFi3dFilBuilder`] and
//! [`crate::chfi3d::builder::ch_builder::ChFi3dChBuilder`]. In all four every
//! parameter and query member is a one-line forwarding, and only `Build` adds
//! the `BRepBuilderAPI_MakeShape` state (`myDone` / `myShape` / `myMap`).
//!
//! Translation notes:
//! - `BRepBuilderAPI_MakeShape::Done()` / `NotDone()` are `done = true` /
//!   `done = false` on the front end's flag; `IsDone()` answers the same flag
//!   (`BRepBuilderAPI_MakeShape.cxx:44-52`),
//! - `Modified(S)` casts `S` with `TopoDS::Edge(S)` (`MakeFillet2d.cxx:146`),
//!   which raises `Standard_TypeMismatch` on a non-edge, and `DescendantEdge`
//!   dereferences `history.Find(E)` (`MakeFillet2d.lxx:78-81`), which raises
//!   `Standard_NoSuchObject` on an unbound edge; the port answers an empty list
//!   for both arms instead of throwing,
//! - `NewEdges(I)` (`MakeFillet2d.cxx:158-171`) reads `FilletEdges()(I)` /
//!   `ChamferEdges()(I - NbFillet())`, and `NCollection_Sequence::Value` raises
//!   `Standard_OutOfRange` outside `1..Length()`; the port answers an empty
//!   list there,
//! - `Descriptor` / `BasisEdge` return a `const TopoDS_Edge&` and OCCT answers
//!   the *input* edge when it is not a descendant (`MakeFillet2d.lxx` /
//!   `MakeFillet2d.cxx:119-122` forwards to `ChFi2d_Builder::BasisEdge`, which
//!   does the same),
//! - the `double&` out parameters of `GetDist` / `Dists` / `GetDistAngle`
//!   (`MakeChamfer.cxx:57-60`, `:84-90`, `:114-117`) became returned values in
//!   [`ChFi3dChBuilder`](crate::chfi3d::builder::ch_builder::ChFi3dChBuilder),
//! - `Edge(I, J)` answers a clone rather than a reference: the port's spine
//!   lives behind a lock guard (see [`LocalOperation::edge`]).

mod local_operation;
mod make_chamfer;
mod make_fillet;
mod make_fillet2d;

pub use local_operation::LocalOperation;
pub use make_chamfer::MakeChamfer;
pub use make_fillet::MakeFillet;
pub use make_fillet2d::MakeFillet2d;
