//! `BRepFeat` package ports (`ModelingAlgorithms/TKFeat/BRepFeat`).
//!
//! This module holds the faithful, branch-by-branch translations of the
//! `BRepFeat` classes. It is **not** [`crate::brepfeat`]: that older module is
//! the self-invented `draft` / `groove` / `neck` / `rib` / `boss_thru_all` API
//! recorded by audit A12, and it stays as-is until those entry points are
//! rewired onto the classes below.
//!
//! Every module mirrors one OCCT class file and keeps the reference's branch
//! structure; the `cxx` / `hxx` / `lxx` line numbers are the physical lines of
//! `D:/source/OCCT-src`.
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | OCCT class | Rust location | Reference | Status |
//! |---|---|---|---|
//! | `BRepFeat_StatusError` | `status_error.rs` | `hxx:21-51` | The enum as [`StatusError`] plus `Print` (`cxx:718-808`) as [`status_error::print`]. |
//! | `BRepFeat_PerfSelection` | `perf_selection.rs` | `hxx:30-36` | The enum as [`PerfSelection`]. |
//! | `BRepFeat` | `statics.rs` | `hxx:87-119`, `cxx:59-808` | `SampleEdges`, `Barycenter`, `ParametricBarycenter`, `ParametricMinMax`, `IsInside`, `IsInOut`, `FaceUntil`, `Tool`; the file-local `IsIn` / `PutInBoundsU` / `PutInBoundsV` are ported in place. |
//!
//! | `BRepFeat_Form` | `form/` | `hxx:1-180`, `cxx:1-1580`, `lxx:1-84` | [`Form::BRepFeatForm`]: the two pure virtuals as [`Form::FormCurves`], `GlobalPerform` (`form/perform.rs`), the gluing route (`form/gluing.rs`), the topological route (`form/topo.rs`), the part selection (`form/selection.rs`), `TransformShapeFU` (`form/transform.rs`), the two `UpdateDescendants` and the static `Descendants` (`form/descendants.rs`). |
//! | `BRepFeat_Builder` | `builder/` | `hxx:1-131`, `cxx:1-827` | [`Builder::BRepFeatBuilder`]: both `Init` / both `SetOperation` / `PartsOfTool` / `KeepParts` / `KeepPart` / `PerformResult` / `RebuildFaces` / `RebuildEdge` / `CheckSolidImages` / `FillRemoved` / the `FillIn3DParts` override / `CheckArgsForOpenSolid -> false`. |
//! | `BRepFeat_Gluer` | `gluer.rs` | `hxx:1-103`, `cxx:1-62`, `lxx:1-69` | `BRepFeatGluer`: the two ctors / `Init` (argument pair inverted, `lxx:25-27`) / both `Bind` / `OpeType` / `BasisShape` / `GluedShape` / `Build` / `IsDeleted` / `Modified`. |
//! | `BRepFeat_SplitShape` | `split_shape.rs` | `hxx:1-121`, `cxx:1-90`, `lxx:1-85` | `BRepFeatSplitShape`: both ctors / `Add(edges)` / `Init` / `SetCheckInterior` / the four `Add` overloads / `DirectLeft` / `Left` / `Right` / `Build` / `IsDeleted` / `Modified`. |
//! | `BRepFeat_MakePrism` | `prism/` | `hxx:1-150`, `cxx:1-1256`, `lxx:1-30` | `MakePrism`: `Init` / `Add` / the eight `Perform` overloads / `Curves` / `BarycCurve`. |
//! | `BRepFeat_MakeRevol` | `revol/` | `hxx:1-105`, `cxx:1-991`, `lxx:1-30` | `MakeRevol`: `Init` / `Add` / `Perform(Angle)` / `Perform(Until)` / `Perform(From, Until)` / `PerformThruAll` / `PerformUntilAngle` / `Curves` / `BarycCurve`. |
//! | `BRepFeat_MakeDPrism` | `dprism/` | `hxx:1-154`, `cxx:1-1315` | `MakeDPrism`: `Init` / `Add` / the seven `Perform` overloads / `BossEdges` / `TopEdges` / `LatEdges` / `Curves` / `BarycCurve`. |
//! | `BRepFeat_MakePipe` | `pipe/` | `hxx:1-122`, `cxx:1-397`, `lxx:19-32` | `MakePipe`: the two ctors / `Init` / `Add` / `Perform` / `Perform(Until)` / `Perform(From, Until)` / `Curves` / `BarycCurve`, plus the file-local `MajMap`. |
//! | `BRepFeat_RibSlot` | `rib_slot/` | `hxx:1-222`, `cxx:1-2671`, `lxx:19-24` | `RibSlot`: the ctor / `LFPerform` / `IsDeleted` / `Modified` / `Generated` / `UpdateDescendants` / `FirstShape` / `LastShape` / `FacesForDraft` / `NewEdges` / `TgtEdges` / `CurrentStatusError`, and the geometry helpers `CheckPoint` (`check_point.rs`), `EdgeExtention` (`edge_extension.rs`), `ChoiceOfFaces` + `HeightMax` (`choice_of_faces.rs`), `PtOnEdgeVertex` (`pt_on_edge_vertex.rs`), `ExtremeFaces` (`extreme_faces.rs`), `SlidingProfile` (`sliding_profile.rs`), `NoSlidingProfile` (`no_sliding_profile.rs`), the file-local `UpdateDescendants` (`descendants.rs`). |
//! | `BRepFeat_MakeLinearForm` | `linear_form/` | `hxx:1-136`, `cxx:1-1366`, `lxx:19-37` | `MakeLinearForm`: `Init` (`init.rs`) / `Perform` (`perform.rs`) / `Propagate` (`propagate.rs`) / the file-local `MajMap` + `SetGluedFaces` (`statics.rs`), on top of `RibSlot`. |
//! | `BRepFeat_MakeRevolutionForm` | `revol_form/` | `hxx:1-125`, `cxx:1-2030`, `lxx:19-43` | `MakeRevolutionForm`: `Init` (`init.rs` + the two profile branches in `init_profile.rs`) / `Perform` (`perform.rs`) / `Propagate` (`propagate.rs`) / `MajMap` + `SetGluedFaces` (`statics.rs`), on top of `RibSlot`. |
//! | `BRepFeat_MakeCylindricalHole` | `cyl_hole/` | `hxx:1-116`, `cxx:1-810`, `lxx:1-54` | `MakeCylindricalHole`: both `Init` / `Status` / `Perform(Radius)` / `PerformThruNext` / `PerformUntilEnd` / `Perform(Radius, PFrom, PTo, Cont)` / `PerformBlind` / `Build` / `Validate`, plus `Baryc` / `BoxParameters` / `GetOffset` / `CreateCyl` (`statics.rs`), on top of `BRepFeat_Builder`. |
//!
//! The Boolean API with history lives next to the features, in `algo_api.rs`
//! (`BRepAlgoAPI_BuilderAlgo` / `_BooleanOperation` / `Cut`).
//!
//! `BRepFeat_Status` (the hole's placement result) is `status.rs`; it is a
//! different enum from `BRepFeat_StatusError` (`status_error.rs`).

pub mod algo_api;
pub mod builder;
pub mod cyl_hole;
pub mod dprism;
pub mod form;
pub mod gluer;
pub mod linear_form;
pub mod perf_selection;
pub mod pipe;
pub mod prism;
pub mod revol_form;
pub mod revol;
pub mod rib_slot;
pub mod split_shape;
pub mod statics;
pub mod status;
pub mod status_error;

pub use algo_api::BooleanOperation;
pub use builder::BRepFeatBuilder;
pub use cyl_hole::MakeCylindricalHole;
pub use dprism::MakeDPrism;
pub use form::{BRepFeatForm, FormCurves};
pub use gluer::FeatureGluer;
pub use linear_form::MakeLinearForm;
pub use perf_selection::PerfSelection;
pub use pipe::MakePipe;
pub use prism::MakePrism;
pub use revol_form::MakeRevolutionForm;
pub use revol::MakeRevol;
pub use rib_slot::RibSlot;
pub use split_shape::BRepFeatSplitShape;
pub use statics::{
    barycenter, face_until, is_in_out, is_inside, parametric_barycenter, parametric_min_max,
    sample_edges, tool,
};
pub use status::Status;
pub use status_error::StatusError;
