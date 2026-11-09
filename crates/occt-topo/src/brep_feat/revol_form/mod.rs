//! `BRepFeat_MakeRevolutionForm` (`hxx:53-125`, `.cxx:1-2030`, `.lxx:19-43`).
//!
//! The rib / groove whose profile is a wire revolved around `myAxe` between
//! `myAngle1` and `myAngle2`. Like `BRepFeat_MakeLinearForm` the class is a
//! `BRepFeat_RibSlot`, so the shared base keeps the maps, the glued faces and
//! `LFPerform`; this module adds the geometry preparation (`Init`), the
//! revolution and its trimming by two parallel planes (`Perform`), the
//! propagation over the basis faces (`Propagate`) and the two file-local
//! helpers (`MajMap`, `SetGluedFaces`).
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | members (`hxx:112-124`) | [`MakeRevolutionForm`] |
//! | both ctors (`lxx:19-38`) | [`MakeRevolutionForm::new`] / [`MakeRevolutionForm::with_shapes`] |
//! | `Init` (`cxx:101-1121`) | `init.rs` + `init_profile.rs` |
//! | `Add` (`cxx:1124-1164`) | [`MakeRevolutionForm::add`] |
//! | `Perform` (`cxx:1168-1438`) | `perform.rs` |
//! | `Propagate` (`cxx:1444-1935`) | `propagate.rs` |
//! | `MajMap` (`cxx:1937-1981`) | [`maj_map`] |
//! | `SetGluedFaces` (`cxx:1983-2031`) | [`set_glued_faces`] |
//!
//! Deviations, all forced by the language or by a documented gap of the
//! underlying port:
//! - the reference's `TopoDS_Shape&` / `TopoDS_Face&` out-parameters are
//!   `Option`s here, and a raise (`Standard_DomainError`,
//!   `StdFail_NotDone`) comes back as the `Err` of the method.
//! - `myPln` is a `Geom_Plane` handle in the reference; the port keeps the
//!   `gp_Pln` value, which is all `Init` and `Perform` read.
//! - `BRepAlgoAPI_Section::Approximation(true)` (`cxx:1469`, `cxx:1743`,
//!   `cxx:1878`) has no effect in the port: `IntTools_FaceFace::SetParameters`
//!   is never called by the BOP filler, the documented gap of
//!   [`crate::brep_algo_api_section`].
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` / `_FEATRIB()` traces are a
//!   DRAW diagnostic stream and are not reproduced.

mod init;
mod init_profile;
mod perform;
mod propagate;
mod statics;

pub use statics::{maj_map, set_glued_faces};

use occt_core::gp::{GpAx1, GpPln};

use crate::abs::ShapeType;
use crate::brep_feat::rib_slot::RibSlot;
use crate::brep_feat::status_error::StatusError;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_MakeRevolutionForm` (`hxx:53-129`).
#[derive(Clone)]
pub struct MakeRevolutionForm {
    /// The `BRepFeat_RibSlot` base (`hxx:53`).
    base: RibSlot,
    /// `gp_Ax1 myAxe` (`hxx:113`).
    my_axe: GpAx1,
    /// `double myHeight1` (`hxx:114`).
    my_height1: f64,
    /// `double myHeight2` (`hxx:115`).
    my_height2: f64,
    /// `bool mySliding` (`hxx:116`).
    #[allow(dead_code)]
    my_sliding: bool,
    /// `occ::handle<Geom_Plane> myPln` (`hxx:117`).
    my_pln: Option<GpPln>,
    /// `double myBnd` (`hxx:118`).
    my_bnd: f64,
    /// `mySlface` (`hxx:119-120`): the sliding faces with the section edges
    /// found on them.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `NCollection_List<TopoDS_Shape> myListOfEdges` (`hxx:121`).
    my_list_of_edges: Vec<TopoShape>,
    /// `double myTol` (`hxx:122`).
    my_tol: f64,
    /// `double myAngle1` (`hxx:123`).
    my_angle1: f64,
    /// `double myAngle2` (`hxx:124`).
    my_angle2: f64,
}

/// `RealLast()` (`Standard_Real.hxx:136-139`), the port's convention is
/// `f64::MAX` (see `brep_feat::statics`).
const REAL_LAST: f64 = f64::MAX;

impl Default for MakeRevolutionForm {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeRevolutionForm {
    /// `BRepFeat_MakeRevolutionForm()` (`lxx:19-27`): `myHeight1` and
    /// `myHeight2` at `0.0`, `myAngle1` and `myAngle2` at `RealLast()`.
    pub fn new() -> Self {
        Self {
            base: RibSlot::new(),
            my_axe: GpAx1::default(),
            my_height1: 0.0,
            my_height2: 0.0,
            my_sliding: false,
            my_pln: None,
            my_bnd: 0.0,
            my_slface: ShapeValueMap::new(),
            my_list_of_edges: Vec::new(),
            my_tol: 0.0,
            my_angle1: REAL_LAST,
            my_angle2: REAL_LAST,
        }
    }

    /// `BRepFeat_MakeRevolutionForm(Sbase, W, P, Axis, H1, H2, Fuse, Sliding)`
    /// (`lxx:31-38`): `Init` with those arguments.
    pub fn with_shapes(
        sbase: &TopoShape,
        w: &Wire,
        p: &GpPln,
        axis: &GpAx1,
        h1: f64,
        h2: f64,
        fuse: i32,
        modify: bool,
    ) -> Result<Self, String> {
        let mut me = Self::new();
        me.init(sbase, w, p, axis, h1, h2, fuse, modify)?;
        Ok(me)
    }

    /// `Add(E, OnFace)` (`cxx:1124-1164`): binds `E` to `F` in `mySlface` when
    /// no sliding face was collected, raising the reference's
    /// `Standard_ConstructionError` when `F` is not a face of the basis.
    pub fn add(&mut self, e: &Edge, on_face: &Face) -> Result<(), String> {
        if !self.my_slface.is_empty() {
            return Ok(());
        }
        // `cxx:1131-1145`.
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

        // `cxx:1147-1163`.
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
