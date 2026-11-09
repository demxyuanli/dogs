//! `Draft_EdgeInfo` (`Draft_EdgeInfo.cxx:26-173`): the new 3D curve and the two
//! pcurves of an edge the draft modification rebuilds.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Draft_EdgeInfo()` | `cxx:26-31` | [`DraftEdgeInfo::new`] |
//! | `Draft_EdgeInfo(HasNewGeometry)` | `cxx:35-41` | [`DraftEdgeInfo::with_new_geometry`] |
//! | `Add(F)` | `cxx:45-57` | [`DraftEdgeInfo::add`] |
//! | `RootFace(F)` | `cxx:61-64` | [`DraftEdgeInfo::root_face`] |
//! | `Tangent(P)` | `cxx:68-72` | [`DraftEdgeInfo::tangent`] |
//! | `IsTangent(P)` | `cxx:76-80` | [`DraftEdgeInfo::is_tangent`] |
//! | `NewGeometry()` / `SetNewGeometry` | `cxx:84-93` | [`DraftEdgeInfo::new_geometry`] / [`DraftEdgeInfo::set_new_geometry`] |
//! | `Geometry()` / `ChangeGeometry()` | `cxx:97-121` | [`DraftEdgeInfo::geometry`] / [`DraftEdgeInfo::change_geometry`] |
//! | `FirstFace()` / `SecondFace()` | `cxx:105-115` | [`DraftEdgeInfo::first_face`] / [`DraftEdgeInfo::second_face`] |
//! | `FirstPC()` / `ChangeFirstPC()` | `cxx:127-147` | [`DraftEdgeInfo::first_pc`] / [`DraftEdgeInfo::change_first_pc`] |
//! | `SecondPC()` / `ChangeSecondPC()` | `cxx:135-156` | [`DraftEdgeInfo::second_pc`] / [`DraftEdgeInfo::change_second_pc`] |
//! | `RootFace()` | `cxx:159-163` | [`DraftEdgeInfo::root_face_value`] |
//! | `Tolerance(tol)` / `Tolerance()` | `cxx:166-173` | [`DraftEdgeInfo::set_tolerance`] / [`DraftEdgeInfo::tolerance`] |
//!
//! # Translation notes
//!
//! * The empty `TopoDS_Face` members (`hxx:75-77`) are `Option<Face>`;
//!   `IsNull()` is `is_none()` and `IsSame(F)` is [`crate::topo_tools::is_same`].
//! * `Add` accumulates the face tolerance with `std::max` (`cxx:56`), and it
//!   records a second face only when that face is not the first one
//!   (`cxx:49-53`); there is no `IsSame` test between the second face and a
//!   later third face, so the second slot keeps the first face it saw.
//! * `myTgt` / `myPt` are the "the new edge is tangent at P" memo filled by
//!   `Propagate`; `myPt` is `gp_Pnt(0, 0, 0)` until `Tangent` runs, which is
//!   what the default-constructed `gp_Pnt` gives in OCCT as well.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::shape::Face;
use crate::topo_tools::is_same;

/// `class Draft_EdgeInfo` (`Draft_EdgeInfo.hxx:29-83`).
#[derive(Clone)]
pub struct DraftEdgeInfo {
    /// `myNewGeom` (`hxx:73`).
    my_new_geom: bool,
    /// `myGeom` (`hxx:74`).
    my_geom: Option<Arc<dyn Curve>>,
    /// `myFirstF` (`hxx:75`).
    my_first_f: Option<Face>,
    /// `mySeconF` (`hxx:76`).
    my_secon_f: Option<Face>,
    /// `myFirstPC` (`hxx:77`).
    my_first_pc: Option<Arc<dyn Curve2d>>,
    /// `mySeconPC` (`hxx:78`).
    my_secon_pc: Option<Arc<dyn Curve2d>>,
    /// `myRootFace` (`hxx:79`).
    my_root_face: Option<Face>,
    /// `myTgt` (`hxx:80`).
    my_tgt: bool,
    /// `myPt` (`hxx:81`).
    my_pt: GpPnt,
    /// `myTol` (`hxx:82`).
    my_tol: f64,
}

impl Default for DraftEdgeInfo {
    fn default() -> Self {
        Self::new()
    }
}

impl DraftEdgeInfo {
    /// `Draft_EdgeInfo()` (`cxx:26-31`).
    pub fn new() -> Self {
        Self::with_new_geometry(false)
    }

    /// `Draft_EdgeInfo(HasNewGeometry)` (`cxx:35-41`).
    pub fn with_new_geometry(has_new_geometry: bool) -> Self {
        Self {
            my_new_geom: has_new_geometry,
            my_geom: None,
            my_first_f: None,
            my_secon_f: None,
            my_first_pc: None,
            my_secon_pc: None,
            my_root_face: None,
            my_tgt: false,
            my_pt: GpPnt::new(0.0, 0.0, 0.0),
            my_tol: 0.0,
        }
    }

    /// `Add(F)` (`cxx:45-57`): files the face in the first free slot and grows
    /// the tolerance.
    pub fn add(&mut self, f: &Face) {
        // `cxx:47-54`
        if self.my_first_f.is_none() {
            self.my_first_f = Some(f.clone());
        } else if !is_same(&self.my_first_f.as_ref().unwrap().0, &f.0) && self.my_secon_f.is_none() {
            self.my_secon_f = Some(f.clone());
        }
        // `cxx:56`
        self.my_tol = self.my_tol.max(BRepTool::face_tolerance(f));
    }

    /// `RootFace(F)` (`cxx:61-64`).
    pub fn root_face(&mut self, f: &Face) {
        self.my_root_face = Some(f.clone());
    }

    /// `Tangent(P)` (`cxx:68-72`).
    pub fn tangent(&mut self, p: &GpPnt) {
        self.my_tgt = true;
        self.my_pt = *p;
    }

    /// `IsTangent(P)` (`cxx:76-80`): answers `myTgt` and copies `myPt` into the
    /// caller's point.
    pub fn is_tangent(&self, p: &mut GpPnt) -> bool {
        *p = self.my_pt;
        self.my_tgt
    }

    /// `NewGeometry()` (`cxx:84-87`).
    pub fn new_geometry(&self) -> bool {
        self.my_new_geom
    }

    /// `SetNewGeometry(NewGeom)` (`cxx:91-94`).
    pub fn set_new_geometry(&mut self, new_geom: bool) {
        self.my_new_geom = new_geom;
    }

    /// `Geometry()` (`cxx:97-100`).
    pub fn geometry(&self) -> Option<&Arc<dyn Curve>> {
        self.my_geom.as_ref()
    }

    /// `ChangeGeometry()` (`cxx:117-121`).
    pub fn change_geometry(&mut self) -> &mut Option<Arc<dyn Curve>> {
        &mut self.my_geom
    }

    /// `FirstFace()` (`cxx:105-109`).
    pub fn first_face(&self) -> Option<&Face> {
        self.my_first_f.as_ref()
    }

    /// `SecondFace()` (`cxx:112-116`).
    pub fn second_face(&self) -> Option<&Face> {
        self.my_secon_f.as_ref()
    }

    /// `FirstPC()` (`cxx:127-131`).
    pub fn first_pc(&self) -> Option<&Arc<dyn Curve2d>> {
        self.my_first_pc.as_ref()
    }

    /// `SecondPC()` (`cxx:135-139`).
    pub fn second_pc(&self) -> Option<&Arc<dyn Curve2d>> {
        self.my_secon_pc.as_ref()
    }

    /// `ChangeFirstPC()` (`cxx:144-148`).
    pub fn change_first_pc(&mut self) -> &mut Option<Arc<dyn Curve2d>> {
        &mut self.my_first_pc
    }

    /// `ChangeSecondPC()` (`cxx:152-156`).
    pub fn change_second_pc(&mut self) -> &mut Option<Arc<dyn Curve2d>> {
        &mut self.my_secon_pc
    }

    /// `RootFace()` (`cxx:159-163`).
    pub fn root_face_value(&self) -> Option<&Face> {
        self.my_root_face.as_ref()
    }

    /// `Tolerance(tol)` (`cxx:166-169`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.my_tol = tol;
    }

    /// `Tolerance()` (`cxx:171-173`).
    pub fn tolerance(&self) -> f64 {
        self.my_tol
    }
}
