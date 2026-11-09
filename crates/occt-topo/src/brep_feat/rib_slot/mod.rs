//! `BRepFeat_RibSlot` -- the shared base of the rib / slot features.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.hxx:1-222`,
//! `.cxx:1-2671`, `.lxx:19-24`.
//!
//! The class is `BRepBuilderAPI_MakeShape` plus the two descendant maps
//! (`myMap` for the basis, `myLFMap` for the profile), the glued faces
//! `myGluedF` and the three edge / face lists the caller reads back. `LFPerform`
//! is the topological reconstruction: it tries the gluer first, and falls back
//! to `BRepFeat_Builder` when the glued faces do not stay inside their
//! counterparts (`cxx:82-259`).
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | `hxx:196-218` (members) | [`RibSlot`] |
//! | ctor (`lxx:19-24`) | [`RibSlot::new`] |
//! | `LFPerform` (`cxx:82-259`) | [`RibSlot::lf_perform`] |
//! | `IsDeleted` / `Modified` / `Generated` (`cxx:263-356`) | the three accessors |
//! | `UpdateDescendants(Gluer)` (`cxx:361-387`) | [`RibSlot::update_descendants_gluer`] |
//! | `FirstShape` / `LastShape` / `FacesForDraft` / `NewEdges` / `TgtEdges` / `CurrentStatusError` (`cxx:390-438`) | the accessors |
//! | `CheckPoint` / `Normal` / `IntPar` (`cxx:443-580`) | `check_point.rs` |
//! | `EdgeExtention` (`cxx:582-643`) | `edge_extension.rs` |
//! | `ChoiceOfFaces` / `HeightMax` (`cxx:645-745`) | `choice_of_faces.rs` |
//! | `ExtremeFaces` (`cxx:747-1325`) | `extreme_faces.rs` |
//! | `PtOnEdgeVertex` (`cxx:1327-1416`) | `pt_on_edge_vertex.rs` |
//! | `SlidingProfile` (`cxx:1418-1791`) | `sliding_profile.rs` |
//! | `NoSlidingProfile` (`cxx:1793-2670`) | `no_sliding_profile.rs` |
//! | `UpdateDescendants(BooleanOperation, Shape, SkipFace)` (`cxx:2672-...`) | `descendants.rs` |
//!
//! Deviations, all forced by the language:
//! - `BRepBuilderAPI_MakeShape`'s `myShape` / `myDone` / `myGenerated` live in
//!   this struct instead of a base object, as they do for the `LocOpe` ports.
//! - the reference returns function-local `static` lists from `Modified` and
//!   `Generated`; the port answers a scratch list of the object
//!   (`my_scratch`), which keeps `myGenerated` free for `FirstShape` /
//!   `LastShape` exactly as the reference does.
//! - `myMap(F)` raises `Standard_NoSuchObject` on an unbound shape; `IsDeleted`
//!   answers `false` there, like [`crate::brep_feat::form::BRepFeatForm`].
//! - `cxx:110-128` walks the faces of `mySUntil` against the faces of `mySbase`
//!   and writes nothing; it has no observable effect and is not reproduced.
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` / `_FEATRIB()` traces are a DRAW
//!   diagnostic stream and are not reproduced.

mod check_point;
mod choice_of_faces;
mod descendants;
mod edge_extension;
mod extreme_faces;
mod make_edges;
mod no_sliding_profile;
mod pt_on_edge_vertex;
mod sliding_profile;

pub use check_point::{check_point, int_par, normal};
pub use choice_of_faces::{choice_of_faces, height_max};
pub use edge_extension::edge_extension;
pub use extreme_faces::ExtremeFaces;
pub use no_sliding_profile::NoSlidingProfileArgs;
pub use pt_on_edge_vertex::{pt_on_edge_vertex, PtOnEdgeVertex};
pub use sliding_profile::{RibProfileError, SlidingProfileArgs};

use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_class3d::SolidClassifier;
use crate::brep_feat::builder::BRepFeatBuilder;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::brep_feat::statics::is_inside;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::locope::shape_maps::{ShapeDataMap, ShapeValueMap};
use crate::locope::{FindEdges, Gluer, LocOpeOperation};
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_RibSlot` (`hxx:90-192`).
#[derive(Clone)]
pub struct RibSlot {
    /// `myFirstPnt` (`hxx:196`).
    pub(crate) my_first_pnt: GpPnt,
    /// `myLastPnt` (`hxx:197`).
    pub(crate) my_last_pnt: GpPnt,
    /// `myFuse` (`hxx:198`).
    pub(crate) my_fuse: bool,
    /// `mySliding` (`hxx:199`).
    #[allow(dead_code)]
    pub(crate) my_sliding: bool,
    /// `myMap` (`hxx:200-201`).
    pub(crate) my_map: ShapeValueMap<Vec<TopoShape>>,
    /// `myLFMap` (`hxx:202-203`).
    pub(crate) my_lf_map: ShapeValueMap<Vec<TopoShape>>,
    /// `myFShape` (`hxx:204`).
    pub(crate) my_f_shape: Option<TopoShape>,
    /// `myLShape` (`hxx:205`).
    pub(crate) my_l_shape: Option<TopoShape>,
    /// `myPerfSelection` (`hxx:206`).
    pub(crate) my_perf_selection: PerfSelection,
    /// `myWire` (`hxx:207`).
    #[allow(dead_code)]
    pub(crate) my_wire: Option<Wire>,
    /// `mySbase` (`hxx:208`).
    pub(crate) my_sbase: Option<TopoShape>,
    /// `mySkface` (`hxx:209`).
    pub(crate) my_skface: Option<Face>,
    /// `myPbase` (`hxx:210`).
    pub(crate) my_pbase: Option<Face>,
    /// `myGShape` (`hxx:211`).
    pub(crate) my_gshape: Option<TopoShape>,
    /// `mySUntil` (`hxx:212`).
    #[allow(dead_code)]
    pub(crate) my_suntil: Option<TopoShape>,
    /// `myGluedF` (`hxx:213`).
    pub(crate) my_glued_f: ShapeDataMap,
    /// `myNewEdges` (`hxx:214`).
    pub(crate) my_new_edges: Vec<TopoShape>,
    /// `myTgtEdges` (`hxx:215`).
    pub(crate) my_tgt_edges: Vec<TopoShape>,
    /// `myFacesForDraft` (`hxx:216`).
    pub(crate) my_faces_for_draft: Vec<TopoShape>,
    /// `myStatusError` (`hxx:217`).
    pub(crate) my_status_error: StatusError,

    /// `BRepBuilderAPI_MakeShape::myShape`.
    pub(crate) my_shape: Option<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myDone`.
    my_done: bool,
    /// `BRepBuilderAPI_MakeShape::myGenerated`.
    pub(crate) my_generated: Vec<TopoShape>,
    /// The result buffer of `Modified` / `Generated`, standing for the
    /// function-local `static NCollection_List` of the reference.
    my_scratch: Vec<TopoShape>,
}

impl Default for RibSlot {
    fn default() -> Self {
        Self::new()
    }
}

impl RibSlot {
    /// `BRepFeat_RibSlot()` (`lxx:19-24`): `myFuse(false)`, `mySliding(false)`,
    /// `myStatusError(BRepFeat_OK)`.
    pub fn new() -> Self {
        Self {
            my_first_pnt: GpPnt::default(),
            my_last_pnt: GpPnt::default(),
            my_fuse: false,
            my_sliding: false,
            my_map: ShapeValueMap::new(),
            my_lf_map: ShapeValueMap::new(),
            my_f_shape: None,
            my_l_shape: None,
            my_perf_selection: PerfSelection::NoSelection,
            my_wire: None,
            my_sbase: None,
            my_skface: None,
            my_pbase: None,
            my_gshape: None,
            my_suntil: None,
            my_glued_f: ShapeDataMap::new(),
            my_new_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_faces_for_draft: Vec::new(),
            my_status_error: StatusError::Ok,
            my_shape: None,
            my_done: false,
            my_generated: Vec::new(),
            my_scratch: Vec::new(),
        }
    }

    /// `Done()` of `BRepBuilderAPI_MakeShape`.
    pub(crate) fn done(&mut self) {
        self.my_done = true;
    }

    /// `NotDone()` of `BRepBuilderAPI_MakeShape`.
    pub(crate) fn not_done(&mut self) {
        self.my_done = false;
    }

    /// `IsDone()` of `BRepBuilderAPI_MakeShape`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` of `BRepBuilderAPI_MakeShape`.
    pub fn shape(&self) -> Option<&TopoShape> {
        if !self.my_done {
            return None;
        }
        self.my_shape.as_ref()
    }

    /// `LFPerform()` (`cxx:82-259`): the topological reconstruction of the rib.
    ///
    /// The gluing route runs only when every glued face of `myGShape` stays
    /// inside the basis face it was bound to and the gluer's operation agrees
    /// with `myFuse`; otherwise `BRepFeat_Builder` selects the parts of the tool
    /// around `myFirstPnt` / `myLastPnt`.
    pub fn lf_perform(&mut self) -> Result<(), String> {
        // `cxx:89-101`.
        if self.my_sbase.is_none()
            || self.my_pbase.is_none()
            || self.my_skface.is_none()
            || self.my_gshape.is_none()
            || self.my_lf_map.is_empty()
        {
            self.my_status_error = StatusError::NotInitialized;
            self.not_done();
            return Ok(());
        }

        // `cxx:104-108`.
        let mut the_ope = 2;
        if !self.my_glued_f.is_empty() {
            the_ope = 1;
        }

        // `cxx:134-187`: the gluing route.
        let mut the_glue = Gluer::new();
        if the_ope == 1 {
            let sbase = self.my_sbase.clone().ok_or_else(|| "Standard_ConstructionError".to_string())?;
            let gshape = self.my_gshape.clone().ok_or_else(|| "Standard_ConstructionError".to_string())?;
            let mut collage = true;
            let mut the_fe = FindEdges::new();
            the_glue.init(&sbase, &gshape);
            let glued: Vec<(TopoShape, TopoShape)> =
                self.my_glued_f.iter().cloned().collect();
            for (key, value) in glued {
                let glface = Face::wrap(key).ok_or_else(|| "Standard_DomainError".to_string())?;
                let fac = Face::wrap(value).ok_or_else(|| "Standard_DomainError".to_string())?;
                // `cxx:150-155`: the glued face has to belong to `myGShape`.
                let found = explorer_shapes(&gshape, ShapeType::Face)
                    .iter()
                    .any(|f| is_same(f, &glface.0));
                if found {
                    collage = is_inside(&glface, &fac);
                    if !collage {
                        the_ope = 2;
                        break;
                    }
                    the_glue.bind_face(&glface, &fac)?;
                    the_fe.set(&glface.0, &fac.0);
                    the_fe.init_iterator();
                    while the_fe.more() {
                        let (from, to) = (the_fe.edge_from().cloned(), the_fe.edge_to().cloned());
                        if let (Some(from), Some(to)) = (from, to) {
                            the_glue.bind_edge(&from, &to)?;
                        }
                        the_fe.next();
                    }
                }
            }
            // `cxx:177-185`.
            let ope = the_glue.ope_type();
            if ope == LocOpeOperation::Invalid
                || (self.my_fuse && ope != LocOpeOperation::Fuse)
                || (!self.my_fuse && ope != LocOpeOperation::Cut)
                || !collage
            {
                the_ope = 2;
            }
        }

        // `cxx:190-211`: gluing is always applicable once the binds agree.
        if the_ope == 1 {
            the_glue.perform()?;
            if the_glue.is_done() {
                self.update_descendants_gluer(&the_glue);
                self.my_new_edges = the_glue.edges().to_vec();
                self.my_tgt_edges = the_glue.tgt_edges().to_vec();
                self.done();
                let resulting = the_glue
                    .resulting_shape()
                    .map_err(|e| e.to_string())?
                    .clone();
                self.my_shape = Some(resulting.clone());
                // `BRepLib::SameParameter(myShape, 1.e-7, Standard_True)`.
                same_parameter_shape(&resulting, 1.0e-7, true);
            } else {
                the_ope = 2;
            }
        }

        // `cxx:213-258`: without gluing, `BRepFeat_Builder` picks the parts.
        if the_ope == 2 {
            let sbase = self.my_sbase.clone().ok_or_else(|| "Standard_ConstructionError".to_string())?;
            let gshape = self.my_gshape.clone().ok_or_else(|| "Standard_ConstructionError".to_string())?;
            let b_flag = self.my_perf_selection != PerfSelection::NoSelection;
            let mut the_builder = BRepFeatBuilder::new();
            the_builder.init_with_tool(&sbase, &gshape);
            the_builder.set_operation_with_flag(if self.my_fuse { 1 } else { 0 }, b_flag);
            the_builder.perform()?;
            if b_flag {
                let parts_of_tool = the_builder.parts_of_tool();
                if !parts_of_tool.is_empty() && self.my_perf_selection != PerfSelection::NoSelection
                {
                    // `cxx:233`: `BRep_Tool::Tolerance(myPbase) * 2`.
                    let toler = self
                        .my_pbase
                        .as_ref()
                        .map(|p| BRepTool::face_tolerance(p) * 2.0)
                        .unwrap_or(0.0);
                    let mut oussa = SolidClassifier::new();
                    for s in parts_of_tool {
                        oussa.load(s.clone());
                        oussa.perform(&self.my_first_pnt, toler);
                        let sp1 = oussa.state();
                        oussa.perform(&self.my_last_pnt, toler);
                        let sp2 = oussa.state();
                        if sp1 != FaceState::Out && sp2 != FaceState::Out {
                            the_builder.keep_part(&s);
                        }
                    }
                }
                // `cxx:250`: `PerformResult` is void in the reference.
                let _ = the_builder.perform_result();
                self.my_shape = Some(the_builder.shape().clone());
            } else {
                self.my_shape = Some(the_builder.shape().clone());
            }
            self.done();
        }
        Ok(())
    }

    /// `IsDeleted(F)` (`cxx:263-267`): `myMap(F).IsEmpty()`.
    pub fn is_deleted(&self, f: &TopoShape) -> bool {
        match self.my_map.find(f) {
            Some(list) => list.is_empty(),
            None => false,
        }
    }

    /// `Modified(F)` (`cxx:270-293`): the descendants of `F` that are not `F`.
    pub fn modified(&mut self, f: &TopoShape) -> &[TopoShape] {
        if let Some(list) = self.my_map.find(f) {
            let list = list.clone();
            self.my_scratch.clear();
            for sh in list {
                if !is_same(&sh, f) {
                    self.my_scratch.push(sh);
                }
            }
            return &self.my_scratch;
        }
        self.my_scratch.clear();
        &self.my_scratch
    }

    /// `Generated(S)` (`cxx:297-356`).
    ///
    /// A `FACE` key answers `myGenerated` untouched; otherwise the `myLFMap`
    /// entry of `S` is expanded through `myMap`, and `myGenerated` is cleared as
    /// the reference does (`cxx:307`, `cxx:330`).
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        if s.shape_type() == ShapeType::Face {
            return &self.my_generated;
        }
        self.my_generated.clear();
        if self.my_lf_map.is_empty() || !self.my_lf_map.is_bound(s) {
            if let Some(list) = self.my_map.find(s) {
                let list = list.clone();
                self.my_scratch.clear();
                for sh in list {
                    if !is_same(&sh, s) {
                        self.my_scratch.push(sh);
                    }
                }
                return &self.my_scratch;
            }
            self.my_scratch.clear();
            return &self.my_scratch;
        }

        // `cxx:329-351`.
        let lf = self.my_lf_map.find(s).cloned().unwrap_or_default();
        self.my_scratch.clear();
        for v in lf {
            if let Some(list) = self.my_map.find(&v).cloned() {
                for sh in list {
                    if !is_same(&sh, s) {
                        self.my_scratch.push(sh);
                    }
                }
            }
        }
        &self.my_scratch
    }

    /// `UpdateDescendants(const LocOpe_Gluer&)` (`cxx:361-387`).
    pub fn update_descendants_gluer(&mut self, g: &Gluer) {
        let entries: Vec<(TopoShape, Vec<TopoShape>)> = self.my_map.iter().cloned().collect();
        for (orig, value) in entries {
            let mut newdsc = crate::locope::shape_maps::ShapeSet::new();
            for v in value {
                // `TopoDS::Face(it.Value())`.
                if let Some(fdsc) = Face::wrap(v) {
                    // The gluer is done here (`cxx:193`), so `DescendantFaces`
                    // does not raise.
                    for d in g.descendant_faces(&fdsc).unwrap_or(&[]) {
                        newdsc.add(d);
                    }
                }
            }
            if let Some(slot) = self.my_map.change(&orig) {
                slot.clear();
                slot.extend_from_slice(newdsc.iter());
            }
        }
    }

    /// `FirstShape()` (`cxx:390-396`).
    pub fn first_shape(&self) -> &[TopoShape] {
        match &self.my_f_shape {
            Some(s) => self.my_map.find(s).map(|v| v.as_slice()).unwrap_or(&[]),
            None => &self.my_generated,
        }
    }

    /// `LastShape()` (`cxx:401-407`).
    pub fn last_shape(&self) -> &[TopoShape] {
        match &self.my_l_shape {
            Some(s) => self.my_map.find(s).map(|v| v.as_slice()).unwrap_or(&[]),
            None => &self.my_generated,
        }
    }

    /// `FacesForDraft()` (`cxx:412-415`).
    pub fn faces_for_draft(&self) -> &[TopoShape] {
        &self.my_faces_for_draft
    }

    /// `NewEdges()` (`cxx:419-422`).
    pub fn new_edges(&self) -> &[TopoShape] {
        &self.my_new_edges
    }

    /// `TgtEdges()` (`cxx:426-429`).
    pub fn tgt_edges(&self) -> &[TopoShape] {
        &self.my_tgt_edges
    }

    /// `CurrentStatusError()` (`cxx:433-436`).
    pub fn current_status_error(&self) -> StatusError {
        self.my_status_error
    }
}
