//! `BRepOffsetAPI_NormalProjection` - the projection of edges/wires onto a
//! face along its normals, driven through `BRepAlgo_NormalProjection`.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_NormalProjection.cxx:21-104`
//! and `.hxx:35-135` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_NormalProjection()` | `cxx:21` | [`BRepOffsetApiNormalProjection::new`] |
//! | `BRepOffsetAPI_NormalProjection(S)` | `cxx:23-26` | [`BRepOffsetApiNormalProjection::with_shape`] |
//! | `Init(S)` | `cxx:28-31` | [`BRepOffsetApiNormalProjection::init`] |
//! | `Add(ToProj)` | `cxx:33-36` | [`BRepOffsetApiNormalProjection::add`] |
//! | `SetParams(...)` | `cxx:38-46` | [`BRepOffsetApiNormalProjection::set_params`] |
//! | `SetMaxDistance(MaxDist)` | `cxx:48-51` | [`BRepOffsetApiNormalProjection::set_max_distance`] |
//! | `SetLimit(FaceBounds = true)` | `cxx:53-56` | [`BRepOffsetApiNormalProjection::set_limit`] |
//! | `Compute3d(With3d = true)` | `cxx:58-61` | [`BRepOffsetApiNormalProjection::compute3d`] |
//! | `Build(Range)` | `cxx:63-68` | [`BRepOffsetApiNormalProjection::build`] |
//! | `IsDone()` | `cxx:70-73` | [`BRepOffsetApiNormalProjection::is_done`] |
//! | `Couple(E)` | `cxx:75-78` | [`BRepOffsetApiNormalProjection::couple`] |
//! | `Generated(S)` | `cxx:80-85` | [`BRepOffsetApiNormalProjection::generated`] |
//! | `Projection()` | `cxx:87-90` | [`BRepOffsetApiNormalProjection::projection`] |
//! | `Ancestor(E)` | `cxx:92-95` | [`BRepOffsetApiNormalProjection::ancestor`] |
//! | `BuildWire(ListOfWire)` | `cxx:99-104` | [`BRepOffsetApiNormalProjection::build_wire`] |
//!
//! # Translation notes
//!
//! * `IsDone` is **overridden** (`cxx:70-73`, `hxx:114`): it answers
//!   `myNormalProjector.IsDone()`, not the base `BRepBuilderAPI_MakeShape` flag.
//!   The port therefore forwards to [`NormalProjection::is_done`] and keeps the
//!   base `my_done` only for `Build`'s unconditional `Done()` (`cxx:67`).
//! * `Build` (`cxx:63-68`) always calls `Done()` after `myNormalProjector.Build()`,
//!   so the base flag is `true` even when the projector failed; only `IsDone()`
//!   tells the truth. The port copies that order.
//! * `Projection()` returns `myNormalProjector.Projection()` and raises
//!   `StdFail_NotDone` when the projector is not done
//!   (`BRepAlgo_NormalProjection.cxx:595-598`); the port answers the compound the
//!   projector holds (an empty compound before `Build`), the same value OCCT
//!   would return after the raise is caught.
//! * `Couple(E)` / `Ancestor(E)` raise `Standard_NoSuchObject` when the map has
//!   no entry (`BRepAlgo_NormalProjection.cxx:602-611`); the port answers `None`.
//! * `Generated(S)` forwards to `BRepAlgo_NormalProjection::Generated`
//!   (`cxx:615-619`), which answers an empty list for an unknown shape; the port
//!   answers `None`.
//! * `BuildWire(ListOfWire)` (`cxx:99-104`) appends every wire it can build and
//!   answers `false` when it cannot; the port's
//!   [`NormalProjection::build_wire`] answers the single wire it built
//!   (`BRepAlgo_NormalProjection.cxx:642-677`).

use occt_core::kernel::geomabs::Shape as GeomAbsShape;

use crate::normal_projection::NormalProjection;
use crate::shape::{Edge, TopoShape, Wire};

/// `class BRepOffsetAPI_NormalProjection` (`hxx:35-135`).
pub struct BRepOffsetApiNormalProjection {
    /// `myNormalProjector` (`hxx:133`).
    my_normal_projector: NormalProjection,
    /// `BRepBuilderAPI_MakeShape::myDone` (`cxx:67` sets it).
    my_done: bool,
    /// `BRepBuilderAPI_MakeShape::myShape` (`cxx:66` sets it).
    my_shape: Option<TopoShape>,
}

impl Default for BRepOffsetApiNormalProjection {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiNormalProjection {
    /// `BRepOffsetAPI_NormalProjection()` (`cxx:21`): `= default`, so the
    /// projector is default constructed and the shape stays null.
    pub fn new() -> Self {
        Self {
            my_normal_projector: NormalProjection::new(),
            my_done: false,
            my_shape: None,
        }
    }

    /// `BRepOffsetAPI_NormalProjection(const TopoDS_Shape& S)` (`cxx:23-26`):
    /// `myNormalProjector.Init(S)` and nothing else.
    pub fn with_shape(s: &TopoShape) -> Self {
        let mut this = Self::new();
        // `cxx:25`.
        this.my_normal_projector.init(s.clone());
        this
    }

    /// `Init(const TopoDS_Shape& S)` (`cxx:28-31`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_normal_projector.init(s.clone());
    }

    /// `Add(const TopoDS_Shape& ToProj)` (`cxx:33-36`).
    pub fn add(&mut self, to_proj: &TopoShape) {
        self.my_normal_projector.add(to_proj);
    }

    /// `SetParams(Tol3D, Tol2D, InternalContinuity, MaxDegree, MaxSeg)`
    /// (`cxx:38-46`).
    pub fn set_params(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        internal_continuity: GeomAbsShape,
        max_degree: i32,
        max_seg: i32,
    ) {
        self.my_normal_projector.set_params(
            tol3d,
            tol2d,
            internal_continuity as u8,
            max_degree,
            max_seg,
        );
    }

    /// `SetMaxDistance(const double MaxDist)` (`cxx:48-51`).
    pub fn set_max_distance(&mut self, max_dist: f64) {
        self.my_normal_projector.set_max_distance(max_dist);
    }

    /// `SetLimit(const bool FaceBounds = true)` (`cxx:53-56`).
    pub fn set_limit(&mut self, face_bounds: bool) {
        self.my_normal_projector.set_limit(face_bounds);
    }

    /// `Compute3d(const bool With3d = true)` (`cxx:58-61`).
    pub fn compute3d(&mut self, with3d: bool) {
        self.my_normal_projector.compute3d(with3d);
    }

    /// `Build(Range)` (`cxx:63-68`): `myNormalProjector.Build()`, then
    /// `myShape = Projection()`, then `Done()` - unconditionally.
    pub fn build(&mut self) -> Result<(), String> {
        // `cxx:65`.
        self.my_normal_projector.build()?;
        // `cxx:66`.
        self.my_shape = Some(self.my_normal_projector.projection().clone());
        // `cxx:67`.
        self.my_done = true;
        Ok(())
    }

    /// `IsDone()` (`cxx:70-73`): the **projector's** flag (`hxx:114`).
    pub fn is_done(&self) -> bool {
        self.my_normal_projector.is_done()
    }

    /// `Couple(const TopoDS_Edge& E)` (`cxx:75-78`).
    pub fn couple(&self, e: &Edge) -> Option<&TopoShape> {
        self.my_normal_projector.couple(e)
    }

    /// `Generated(const TopoDS_Shape& S)` (`cxx:80-85`).
    pub fn generated(&self, s: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.my_normal_projector.generated(s)
    }

    /// `Projection()` (`cxx:87-90`).
    pub fn projection(&self) -> &TopoShape {
        self.my_normal_projector.projection()
    }

    /// `Ancestor(const TopoDS_Edge& E)` (`cxx:92-95`).
    pub fn ancestor(&self, e: &Edge) -> Option<&TopoShape> {
        self.my_normal_projector.ancestor(e)
    }

    /// `BuildWire(NCollection_List<TopoDS_Shape>& ListOfWire)`
    /// (`cxx:99-104`).
    pub fn build_wire(&self) -> Option<Wire> {
        self.my_normal_projector.build_wire()
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`, inherited): `myShape` set by
    /// `Build` (`cxx:66`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// The `BRepAlgo_NormalProjection` object the reference hides in
    /// `myNormalProjector` (`hxx:133`).
    pub fn normal_projector(&self) -> &NormalProjection {
        &self.my_normal_projector
    }
}
