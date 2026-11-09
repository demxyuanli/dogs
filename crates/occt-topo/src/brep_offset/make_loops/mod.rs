//! `BRepOffset_MakeLoops` (`BRepOffset_MakeLoops.hxx` 63 lines / `.cxx` 661
//! lines):
//! rebuilds the wires (loops) of the faces produced by an offset, recording
//! every new face and cut edge in a `BRepAlgo_Image`.
//!
//! Ported here so far (control flow taken from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_MakeLoops.cxx`):
//!
//! - `BRepOffset_MakeLoops()` (`:38`) - [`BRepOffsetMakeLoops::new`]
//! - `Build(LF, AsDes, Image, theImageVV, theRange)` (`:42-179`) - `build.rs`
//! - file-local `IsBetweenCorks` (`:184-213`) and
//!   `BuildOnContext(LContext, Analyse, AsDes, Image, InSide, theRange)`
//!   (`:217-451`) - `build_on_context.rs`
//! - `BuildFaces(LF, AsDes, Image, theRange)` (`:455-661`) - `build_faces.rs`
//!
//! The vertex-substitution tail that all three entry points share verbatim
//! (`:125-178`, `:405-450`, `:613-660`) lives once in `substitute.rs`; the three
//! call sites differ only in the shape list they walk.
//!
//! Dependencies the three methods drive are already ported:
//! `BRepAlgo_Loop` (`brep_algo_loop/`, `Init`/`AddEdge`/`AddConstEdge`/
//! `Perform`/`WiresToFaces`/`NewFaces`/`NewEdges`/`VerticesForSubstitute`),
//! `BRepAlgo_AsDes` (`brep_algo/as_des.rs`, `Descendant`/`HasDescendant`/
//! `Ascendant`/`HasAscendant`/`Replace`) and `BRepAlgo_Image`
//! (`brep_algo_image.rs`, `Bind`/`Add`/`Image`/`LastImage`/`HasImage`).
//!
//! Deviations, all "OCCT raises, the port skips":
//!
//! * `Loops.NewEdges(E)` (`:117`, `:397`, `:606`) is
//!   `NCollection_DataMap::operator()`: unbound `E` throws
//!   `Standard_NoSuchObject`. [`BRepAlgoLoop::new_edges`] answers `None` and the
//!   port leaves the image untouched.
//! * `Image.Bind(F, NF)` (`:113`, `:384`, `:592`) throws
//!   `Standard_ConstructionError` when `F` is already imaged; the port ignores
//!   the failure (`BRepAlgoImage::bind_list` answers `Err`).
//! * `Image.Add(E, LoopNE)` (`:119`, `:396`, `:608`) throws when `E` has no
//!   image; the guarded arms above never reach it unbound.
//! * the reference dereferences `myVerVerMap(V)` unguarded (`:161`, `:433`,
//!   `:642`); the port reads the map through `get`, so a dangling binding is
//!   skipped rather than faulting.
//! * `Message_ProgressScope` (`:52-53`, `:234`, `:469`) is not ported.

use std::collections::HashMap;

use crate::shape::TopoShape;
use crate::topo_tools::{is_same, ShapeMapKey};

mod build;
mod build_faces;
mod build_on_context;
mod substitute;

pub use substitute::substitute_vertices;

/// `class BRepOffset_MakeLoops` (`hxx:35-62`).
#[derive(Default)]
pub struct BRepOffsetMakeLoops {
    /// `myVerVerMap` (`hxx:60`):
    /// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>` of the vertices
    /// `BRepAlgo_Loop` wants substituted, shared by the three entry points.
    my_ver_ver_map: HashMap<ShapeMapKey, TopoShape>,
}

impl BRepOffsetMakeLoops {
    /// `BRepOffset_MakeLoops()` (`cxx:38`).
    pub fn new() -> Self {
        Self::default()
    }

    /// The current `myVerVerMap`, for callers that inspect the substitution
    /// table after a build.
    pub fn vertices_for_substitute(&self) -> &HashMap<ShapeMapKey, TopoShape> {
        &self.my_ver_ver_map
    }

    /// The `found = S.IsSame(it2.Value())` scan of `IsBetweenCorks`
    /// (`:198-208`): `true` when every ascendant of `E` also occurs in
    /// `l_context`.
    pub(crate) fn is_between_corks(
        &self,
        e: &TopoShape,
        as_des: &crate::brep_algo::AsDes,
        l_context: &[TopoShape],
    ) -> bool {
        // `if (!AsDes->HasAscendant(E)) return true`.
        if !as_des.has_ascendant(e) {
            return true;
        }
        as_des
            .ascendant(e)
            .iter()
            .all(|s| l_context.iter().any(|c| is_same(s, c)))
    }
}
