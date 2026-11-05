//! `BRepOffsetAPI_MakePipe` — the pipe built by sweeping a profile along a G1
//! spine.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakePipe.cxx:30-133`
//! and `.hxx:47-105` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(Spine, Profile)` | `cxx:30-36` | [`BRepOffsetApiMakePipe::new`] |
//! | `(Spine, Profile, aMode, ForceApproxC1)` | `cxx:50-57` | [`BRepOffsetApiMakePipe::with_mode`] |
//! | `Pipe()` | `cxx:61-64` | [`BRepOffsetApiMakePipe::pipe`] |
//! | `Build(Range)` | `cxx:68-83` | [`BRepOffsetApiMakePipe::build`] |
//! | `FirstShape()` | `cxx:86-89` | [`BRepOffsetApiMakePipe::first_shape`] |
//! | `LastShape()` | `cxx:93-96` | [`BRepOffsetApiMakePipe::last_shape`] |
//! | `Generated(S)` | `cxx:100-104` | [`BRepOffsetApiMakePipe::generated`] |
//! | `Generated(SSpine, SProfile)` | `cxx:111-127` | [`BRepOffsetApiMakePipe::generated_from`] |
//! | `ErrorOnSurface()` | `cxx:130-133` | [`BRepOffsetApiMakePipe::error_on_surface`] |
//!
//! `BRepPrimAPI_MakeSweep` has no `.cxx` and contributes only the two virtual
//! queries; the port models it with [`MakeSweep`].
//!
//! # Translation notes
//!
//! * Both constructors run `Build()` (`cxx:33`, `cxx:55`), so a construction
//!   never returns without a done/not-done verdict.
//! * The two-argument constructor relies on the `BRepFill_Pipe` defaults of
//!   `BRepFill_Pipe.hxx:79-83`: `GeomFill_IsCorrectedFrenet`, `ForceApproxC1 =
//!   false`, `GeneratePartCase = false`.
//! * `NotDone()` / `Done()` are the base-class flags, kept here as `my_done`
//!   because the port mirrors each wrapper with its own state triple.
//! * `Build` uses `TopExp::MapShapes` with `TopTools_ShapeMapHasher`
//!   (`cxx:71-72`), i.e. the dedup key is `TShape` + `Location`;
//!   [`crate::topo_tools::all_subshapes`] uses the same key and skips nothing
//!   but, unlike `TopExp::MapShapes`, it seeds the map with the argument
//!   itself. `TopExp::MapShapes` also maps the shape itself, so the two counts
//!   agree. A null shape maps to `0` in both, which keeps the `Extent() == 1`
//!   test answering `Done()` for an empty pipe exactly as OCCT does.
//! * `Generated(SSpine, SProfile)` returns a null shape for anything that is
//!   neither an edge nor a vertex; the port answers `None` for that `bid`
//!   (`cxx:124-125`).

use occt_geom::geom_fill::GeomFillTrihedron;

use crate::abs::ShapeType;
use crate::brep_fill_pipe::Pipe;
use crate::brep_prim_api_make_sweep::MakeSweep;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::all_subshapes;

/// `BRepOffsetAPI_MakePipe` (`BRepOffsetAPI_MakePipe.hxx:47-105`).
pub struct BRepOffsetApiMakePipe {
    /// `myPipe`.
    my_pipe: Pipe,
    /// `myShape`.
    my_shape: Option<TopoShape>,
    /// `myGenerated`.
    my_generated: Vec<TopoShape>,
    /// `myDone` (the `BRepBuilderAPI_MakeShape` flag).
    my_done: bool,
}

/// The `TopoDS::Edge` / `TopoDS::Vertex` down-cast failure of
/// `cxx:114-121` (`Standard_TypeMismatch`).
fn type_mismatch() -> String {
    "Standard_TypeMismatch: BRepOffsetAPI_MakePipe::Generated".to_string()
}

impl BRepOffsetApiMakePipe {
    /// `BRepOffsetAPI_MakePipe(Spine, Profile)` (`cxx:30-36`): the two-argument
    /// form of `BRepFill_Pipe` plus `Build()`.
    pub fn new(spine: &Wire, profile: &TopoShape) -> Result<Self, String> {
        Self::with_mode(spine, profile, GeomFillTrihedron::CorrectedFrenet, false)
    }

    /// `BRepOffsetAPI_MakePipe(Spine, Profile, aMode, ForceApproxC1)`
    /// (`cxx:50-57`).
    pub fn with_mode(
        spine: &Wire,
        profile: &TopoShape,
        a_mode: GeomFillTrihedron,
        force_approx_c1: bool,
    ) -> Result<Self, String> {
        // `myPipe(Spine, Profile, aMode, ForceApproxC1)` (`cxx:53`); the
        // `GeneratePartCase` default is false (`BRepFill_Pipe.hxx:82`).
        let my_pipe = Pipe::with_spine(spine, profile, a_mode, force_approx_c1, false)?;
        let mut this = Self {
            my_pipe,
            my_shape: None,
            my_generated: Vec::new(),
            my_done: false,
        };
        // `cxx:55`.
        this.build();
        Ok(this)
    }

    /// `Pipe()` (`cxx:61-64`).
    pub fn pipe(&self) -> &Pipe {
        &self.my_pipe
    }

    /// `Build(theRange)` (`cxx:68-83`): copy the pipe's shape and answer not
    /// done when it came out empty.
    pub fn build(&mut self) {
        // `cxx:70`.
        self.my_shape = self.my_pipe.shape().cloned();
        // `cxx:71-72`: `TopExp::MapShapes(myShape, theMap)`.
        let extent = match &self.my_shape {
            Some(shape) if !shape.is_null() => all_subshapes(shape).len(),
            _ => 0,
        };
        // `cxx:73-81`.
        self.my_done = extent != 1;
    }

    /// `IsDone()` — the `BRepBuilderAPI_MakeShape` flag set by `Build`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` — the `myShape` written by `Build`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `FirstShape()` (`cxx:86-89`).
    pub fn first_shape(&self) -> Option<TopoShape> {
        self.my_pipe.first_shape().cloned()
    }

    /// `LastShape()` (`cxx:93-96`).
    pub fn last_shape(&self) -> Option<TopoShape> {
        self.my_pipe.last_shape().cloned()
    }

    /// `Generated(S)` (`cxx:100-104`): refill `myGenerated` from the pipe and
    /// answer it.
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.my_generated = self.my_pipe.generated(s);
        &self.my_generated
    }

    /// `Generated(SSpine, SProfile)` (`cxx:111-127`): an edge profile answers
    /// the swept face, a vertex profile the swept edge, anything else `None`.
    pub fn generated_from(
        &self,
        s_spine: &TopoShape,
        s_profile: &TopoShape,
    ) -> Result<Option<TopoShape>, String> {
        match s_profile.shape_type() {
            // `cxx:114-117`: `TopoDS::Edge(SSpine)` / `TopoDS::Edge(SProfile)`.
            ShapeType::Edge => {
                let spine = Edge::wrap(s_spine.clone()).ok_or_else(type_mismatch)?;
                let profile = Edge::wrap(s_profile.clone()).ok_or_else(type_mismatch)?;
                Ok(Some(TopoShape::from(
                    self.my_pipe.face(&spine, &profile)?,
                )))
            }
            // `cxx:118-121`: `TopoDS::Edge(SSpine)` / `TopoDS::Vertex(SProfile)`.
            ShapeType::Vertex => {
                let spine = Edge::wrap(s_spine.clone()).ok_or_else(type_mismatch)?;
                let profile = Vertex::wrap(s_profile.clone()).ok_or_else(type_mismatch)?;
                Ok(Some(TopoShape::from(
                    self.my_pipe.edge(&spine, &profile)?,
                )))
            }
            // `cxx:124-125`: `bid` is left null ("POP pour NT").
            _ => Ok(None),
        }
    }

    /// `ErrorOnSurface()` (`cxx:130-133`).
    pub fn error_on_surface(&self) -> f64 {
        self.my_pipe.error_on_surface()
    }
}

impl MakeSweep for BRepOffsetApiMakePipe {
    fn first_shape(&self) -> Option<TopoShape> {
        BRepOffsetApiMakePipe::first_shape(self)
    }

    fn last_shape(&self) -> Option<TopoShape> {
        BRepOffsetApiMakePipe::last_shape(self)
    }
}
