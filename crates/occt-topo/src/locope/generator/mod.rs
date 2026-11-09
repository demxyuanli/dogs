//! `LocOpe_Generator` (`LocOpe_Generator.hxx:31-67`, `.cxx:71-1230`,
//! `.lxx:20-64`): rebuilds a base shape after a local operation generated faces
//! from some of its sub-shapes.
//!
//! `Perform` (`cxx:71-1218`) walks the base shape against a
//! [`GeneratedShape`]:
//!
//! 1. the faces on the left of the generating edges become the seeds
//!    (`cxx:90-131`);
//! 2. a generated face that is analytically the same plane as a base face is
//!    merged with it, and the connecting edges go to `toRemove`
//!    (`cxx:144-265`);
//! 3. every seam/twin edge of the merged faces is paired in `theEEMap`, which
//!    decides per vertex whether an edge fuses, and records the cases that must
//!    not fuse in `DontFuse` (`cxx:267-439`);
//! 4. the fused edges are rebuilt with the vertices of their twins
//!    (`cxx:454-522`);
//! 5. the merged faces are rebuilt from the surviving edges, with p-curves
//!    recomputed on non-planar surfaces (`cxx:526-1090`);
//! 6. `LocOpe_BuildShape` rebuilds shells and solids from the loose face list,
//!    and a single solid gets the orientation of its reference face back
//!    (`cxx:1130-1218`).
//!
//! `myGen` (`hxx:64`) is declared by the reference but never written by
//! `Perform`, so it is not carried here.

mod perform;
mod statics;

use crate::locope::generated_shape::GeneratedShape;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, Face, TopoShape};

/// `LocOpe_Generator` (`hxx:31-67`).
#[derive(Default)]
pub struct Generator {
    /// `myShape` (`hxx:62`).
    my_shape: Option<TopoShape>,
    /// `myDone` (`hxx:65`).
    my_done: bool,
    /// `myRes` (`hxx:66`).
    my_res: Option<TopoShape>,
    /// `myModShapes` (`hxx:67`): the descendant list of every face of the base
    /// shape *and* of every generated face (`cxx:518-520`, `cxx:1084-1088`,
    /// `cxx:1148-1186`).
    my_mod_shapes: ShapeValueMap<Vec<TopoShape>>,
}

impl Generator {
    /// `LocOpe_Generator()` (`.lxx:22-24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_Generator(const TopoDS_Shape&)` (`.lxx:28-32`).
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            my_shape: Some(s.clone()),
            ..Self::default()
        }
    }

    /// `Init(const TopoDS_Shape&)` (`.lxx:36-40`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_shape = Some(s.clone());
        self.my_done = false;
    }

    /// `Perform(const occ::handle<LocOpe_GeneratedShape>&)` (`cxx:71-1218`).
    ///
    /// `Err` carries `Standard_NullObject` (`cxx:74-77`).
    pub fn perform(&mut self, g: &dyn GeneratedShape) -> Result<(), &'static str> {
        let Some(shape) = self.my_shape.clone() else {
            return Err("Standard_NullObject");
        };
        // `cxx:78-82`.
        self.my_done = false;
        self.my_res = None;
        self.my_mod_shapes.clear();

        let mut stage = perform::Stage::new(shape, g);
        stage.run(&mut self.my_mod_shapes);
        self.my_res = stage.finish(&mut self.my_mod_shapes);
        self.my_done = true;
        Ok(())
    }

    /// `IsDone()` (`.lxx:44-47`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`.lxx:52-55`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `ResultingShape()` (`.lxx:59-64`); `Err` is `StdFail_NotDone`.
    pub fn resulting_shape(&self) -> Result<&TopoShape, &'static str> {
        if !self.my_done {
            return Err("StdFail_NotDone");
        }
        self.my_res.as_ref().ok_or("StdFail_NotDone")
    }

    /// `DescendantFace(const TopoDS_Face&)` (`cxx:1222-1230`).
    ///
    /// `Ok(None)` is the reference's unbound `myModShapes(F)`, which raises
    /// `Standard_NoSuchObject` there; `Err` is `StdFail_NotDone` (`cxx:1226`).
    /// A null entry is a face that disappeared (`hxx:55-58`).
    pub fn descendant_face(&self, f: &Face) -> Result<Option<&Vec<TopoShape>>, &'static str> {
        if !self.my_done {
            return Err("StdFail_NotDone");
        }
        Ok(self.my_mod_shapes.find(&f.0))
    }
}

/// `TopExp::Vertices(const TopoDS_Edge&, Vf, Vl)` on a raw shape: the two
/// vertices of the edge, first then last.
pub(super) fn edge_vertex_shapes(e: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let Some(edge) = Edge::wrap(e.clone()) else {
        return (None, None);
    };
    let (vf, vl) = crate::topo_tools::edge_vertices(&edge);
    (vf.map(|v| v.0), vl.map(|v| v.0))
}

/// `TopoDS_Shape::EmptyCopied()` on an edge: `BRep_TEdge::EmptyCopy`
/// (`BRep_TEdge.cxx:105-133`) keeps the geometry (curve, range, tolerance,
/// flags, every `CurveOnSurface` representation) and drops the children, and
/// `TopoDS_Shape::EmptyCopied` (`TopoDS_Shape.cxx`) keeps the location and
/// resets the orientation to FORWARD.
///
/// The `LocOpe_SplitShape` port carries the same helper without the location;
/// this one keeps it, since `cxx:462-466` copies an edge that may live under a
/// located face.
pub(super) fn empty_copied_edge(e: &Edge) -> Edge {
    let registry = crate::tgeometry::GeometryRegistry::global();
    let mut copy = Edge::new();
    if let Some(geom) = registry.edge_geom(&e.0) {
        registry.set_edge(&copy.0, geom);
    }
    copy.0.set_location(e.0.location());
    copy.0.set_orientation(crate::abs::Orientation::Forward);
    copy
}
