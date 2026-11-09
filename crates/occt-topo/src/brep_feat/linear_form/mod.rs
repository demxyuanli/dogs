//! `BRepFeat_MakeLinearForm` (`hxx:1-136`, `.cxx:1-1366`, `.lxx:19-37`).
//!
//! The rib / groove whose profile is a wire swept along a translation between
//! two vectors. The class is a `BRepFeat_RibSlot`: the shared base keeps the
//! maps, the glued faces and `LFPerform`, and this module adds the geometry
//! preparation (`Init`), the prism (`Perform`), the propagation of the profile
//! over the basis faces (`Propagate`) and the two file-local helpers
//! (`MajMap`, `SetGluedFaces`).
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | members (`hxx:125-134`) | [`MakeLinearForm`] |
//! | both ctors (`lxx:19-31`) | [`MakeLinearForm::new`] / [`MakeLinearForm::with_shapes`] |
//! | `Init` (`cxx:84-857`) | `init.rs` |
//! | `Add` (`cxx:858-898`) | [`MakeLinearForm::add`] |
//! | `Perform` (`cxx:905-1028`) | `perform.rs` |
//! | `Propagate` (`cxx:1035-1288`) | `propagate.rs` |
//! | `MajMap` (`cxx:1292-1334`) | [`maj_map`] |
//! | `SetGluedFaces` (`cxx:1336-1366`) | [`set_glued_faces`] |
//!
//! Deviations, all forced by the language or by a documented gap of the
//! underlying port:
//! - `TransformShapeFU(int)` is declared by the header (`hxx:114`) but has no
//!   definition anywhere in OCCT; the port does not declare it either.
//! - `myCrv` (`hxx:126`) is never read by the reference; the port keeps the
//!   field so the reference state stays visible.
//! - `myPln` is a `Geom_Plane` handle that `Init` transforms in place
//!   (`cxx:195`); the port keeps the `gp_Pln` value, so `Perform` sees the same
//!   plane the reference sees.
//! - `BRepAlgoAPI_Section::Approximation(true)` (`cxx:1060-1061`,
//!   `cxx:1212-1213`) has no effect in the port: `IntTools_FaceFace::SetParameters`
//!   is never called by the BOP filler, which is the documented gap of
//!   [`crate::brep_algo_api_section`].
//! - `BRep_Tool::IsClosed(W)` on the built wire and the `BB.UpdateVertex`
//!   tolerances are reproduced through [`crate::brep_tool`] / [`crate::builder`].

mod init;
mod perform;
mod propagate;
mod statics;

pub use statics::{maj_map, set_glued_faces};

use std::sync::Arc;

use occt_core::gp::{GpPln, GpVec};
use occt_geom::curve::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::rib_slot::RibSlot;
use crate::brep_feat::status_error::StatusError;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_MakeLinearForm` (`hxx:53-135`).
#[derive(Clone)]
pub struct MakeLinearForm {
    /// The `BRepFeat_RibSlot` base (`hxx:53`).
    base: RibSlot,
    /// `occ::handle<Geom_Curve> myCrv` (`hxx:126`), never read (`cxx:792-821`
    /// are the same in `BRepFeat_MakeRevolutionForm`, where it is used).
    #[allow(dead_code)]
    my_crv: Option<Arc<dyn Curve>>,
    /// `gp_Vec myDir` (`hxx:127`).
    my_dir: GpVec,
    /// `gp_Vec myDir1` (`hxx:128`).
    my_dir1: GpVec,
    /// `occ::handle<Geom_Plane> myPln` (`hxx:129`).
    my_pln: Option<GpPln>,
    /// `double myBnd` (`hxx:130`).
    my_bnd: f64,
    /// `mySlface` (`hxx:131-132`): the sliding faces with the section edges
    /// found on them.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `NCollection_List<TopoDS_Shape> myListOfEdges` (`hxx:133`).
    my_list_of_edges: Vec<TopoShape>,
    /// `double myTol` (`hxx:134`).
    my_tol: f64,
}

impl Default for MakeLinearForm {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeLinearForm {
    /// `BRepFeat_MakeLinearForm()` (`lxx:19-21`): `myBnd(0.0)`, `myTol(0.0)`.
    pub fn new() -> Self {
        Self {
            base: RibSlot::new(),
            my_crv: None,
            my_dir: GpVec::zero(),
            my_dir1: GpVec::zero(),
            my_pln: None,
            my_bnd: 0.0,
            my_slface: ShapeValueMap::new(),
            my_list_of_edges: Vec::new(),
            my_tol: 0.0,
        }
    }

    /// `BRepFeat_MakeLinearForm(Sbase, W, P, Direction, Direction1, Fuse,
    /// Modify)` (`lxx:25-31`): `Init` with those arguments.
    pub fn with_shapes(
        sbase: &TopoShape,
        w: &Wire,
        p: &GpPln,
        direction: &GpVec,
        direction1: &GpVec,
        fuse: i32,
        modify: bool,
    ) -> Result<Self, String> {
        let mut me = Self::new();
        me.init(sbase, w, p, direction, direction1, fuse, modify)?;
        Ok(me)
    }

    /// `Add(E, OnFace)` (`cxx:858-898`): binds `E` to `F` in `mySlface` when no
    /// sliding face was collected, raising the reference's
    /// `Standard_ConstructionError` when `F` is not a face of the basis.
    pub fn add(&mut self, e: &Edge, on_face: &Face) -> Result<(), String> {
        if !self.my_slface.is_empty() {
            return Ok(());
        }
        // `cxx:866-877`: `F` has to belong to the basis.
        let in_basis = self
            .base
            .my_sbase
            .as_ref()
            .map(|s| {
                explorer_shapes(s, ShapeType::Face)
                    .iter()
                    .any(|f| is_same(f, &on_face.0))
            })
            .unwrap_or(false);
        if !in_basis {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:880-897`.
        if !self.my_slface.is_bound(&on_face.0) {
            self.my_slface.bind(&on_face.0, Vec::new());
        }
        let already = self
            .my_slface
            .find(&on_face.0)
            .map(|list| list.iter().any(|s| is_same(s, &e.0)))
            .unwrap_or(false);
        if !already {
            if let Some(list) = self.my_slface.change(&on_face.0) {
                list.push(e.0.clone());
            }
        }
        Ok(())
    }

    /// `IsDone()` of `BRepBuilderAPI_MakeShape`.
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// `Shape()` of `BRepBuilderAPI_MakeShape`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.base.shape()
    }

    /// `IsDeleted(S)` (`BRepFeat_RibSlot.cxx:263-267`).
    pub fn is_deleted(&self, s: &TopoShape) -> bool {
        self.base.is_deleted(s)
    }

    /// `Modified(S)` (`BRepFeat_RibSlot.cxx:270-293`).
    pub fn modified(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.base.modified(s)
    }

    /// `Generated(S)` (`BRepFeat_RibSlot.cxx:297-356`).
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.base.generated(s)
    }

    /// `FirstShape()` (`BRepFeat_RibSlot.cxx:390-396`).
    pub fn first_shape(&self) -> &[TopoShape] {
        self.base.first_shape()
    }

    /// `LastShape()` (`BRepFeat_RibSlot.cxx:401-407`).
    pub fn last_shape(&self) -> &[TopoShape] {
        self.base.last_shape()
    }

    /// `FacesForDraft()` (`BRepFeat_RibSlot.cxx:412-415`).
    pub fn faces_for_draft(&self) -> &[TopoShape] {
        self.base.faces_for_draft()
    }

    /// `NewEdges()` (`BRepFeat_RibSlot.cxx:419-422`).
    pub fn new_edges(&self) -> &[TopoShape] {
        self.base.new_edges()
    }

    /// `TgtEdges()` (`BRepFeat_RibSlot.cxx:426-429`).
    pub fn tgt_edges(&self) -> &[TopoShape] {
        self.base.tgt_edges()
    }

    /// `CurrentStatusError()` (`BRepFeat_RibSlot.cxx:433-436`).
    pub fn current_status_error(&self) -> StatusError {
        self.base.current_status_error()
    }

    /// The shared base, for the callers that read the `RibSlot` fields.
    pub fn base(&self) -> &RibSlot {
        &self.base
    }
}
