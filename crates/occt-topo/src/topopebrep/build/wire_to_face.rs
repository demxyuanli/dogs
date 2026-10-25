//! `TopOpeBRepBuild_WireToFace` — rebuild the faces of a set of wires.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_WireToFace.cxx`
//! (OCCT 8.0.0, 60 lines): `Init` (`:28-32`), `AddWire` (`:36-39`) and
//! `MakeFaces` (`:43-61`). `MakeFaces` feeds the wires into a `WireEdgeSet`,
//! builds a `FaceBuilder` on it with `ForceClass`, then asks a throwaway
//! `Builder` to run `MakeFaces`.

use crate::shape::TopoShape;

use super::builder::Builder;
use super::face_builder::FaceBuilder;
use super::wire_edge_set::WireEdgeSet;

/// `class TopOpeBRepBuild_WireToFace` (`hxx`).
pub struct WireToFace {
    /// `myLW` (`hxx:62`).
    my_lw: Vec<TopoShape>,
}

impl Default for WireToFace {
    fn default() -> Self {
        Self::new()
    }
}

impl WireToFace {
    /// `TopOpeBRepBuild_WireToFace()` (`cxx:23`).
    pub fn new() -> Self {
        Self { my_lw: Vec::new() }
    }

    /// `Init()` (`cxx:28-32`): `myLW.Clear()`.
    pub fn init(&mut self) {
        self.my_lw.clear();
    }

    /// `AddWire(W)` (`cxx:36-39`): `myLW.Append(W)`.
    pub fn add_wire(&mut self, w: &TopoShape) {
        self.my_lw.push(w.clone());
    }

    /// `MakeFaces(F, LF)` (`cxx:43-61`).
    ///
    /// `F` is the face the wires belong to; `LF` is cleared then filled with
    /// the faces `MakeFaces` produces.
    pub fn make_faces(&self, f: &TopoShape, builder: &mut Builder, lf: &mut Vec<TopoShape>) {
        lf.clear();

        // cxx:47-52
        let mut wes = WireEdgeSet::new(f);
        for w in &self.my_lw {
            wes.add_shape(w);
        }

        // cxx:54-57
        let force_class = true;
        let mut fb = FaceBuilder::with_wire_edge_set(&mut wes, f, force_class);

        // cxx:58-59: `TopOpeBRepDS_BuildTool BT(TopOpeBRepTool_APPROX)` then a
        // `Builder` on it. The port's `BuildTool` carries no approximation mode
        // (see `ds/build_tool`), so the caller passes `BuildTool::new()`.

        // cxx:60
        builder.make_faces(f, &mut fb, lf);
    }
}
