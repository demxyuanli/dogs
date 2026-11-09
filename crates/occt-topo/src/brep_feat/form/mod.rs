//! `BRepFeat_Form` - the shared engine of the form features.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Form.hxx:1-180`,
//! `BRepFeat_Form.cxx:1-1580`, `BRepFeat_Form.lxx:17-84`.
//!
//! `GlobalPerform` has two mutually exclusive routes:
//!
//! - the **gluer** route (`theOpe == 1`, `cxx:332-638`): every glued face must
//!   reproduce the basis face it was bound to, then `LocOpe_Gluer` rebuilds the
//!   shape locally;
//! - the **topological** route (`theOpe == 2`, `cxx:653-1236`): the `From` /
//!   `Until` shapes are turned into solids with [`crate::brep_feat::statics::tool`],
//!   cut away from the generated shape with a history-carrying
//!   [`crate::brep_feat::algo_api::BooleanOperation`], and the retained parts of
//!   the tool are selected by [`crate::brep_feat::builder::BRepFeatBuilder`].
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | `hxx:66-180` (members, `lxx:17-84` valid flags) | this file |
//! | `GlobalPerform` (`cxx:59-1240`) | `perform.rs` |
//! | gluer route (`cxx:332-638`) | `gluing.rs` |
//! | topological route (`cxx:653-811`) | `topo.rs` |
//! | part selection (`cxx:812-1211`) | `selection.rs` |
//! | `TransformShapeFU` (`cxx:1378-1480`) | `transform.rs` |
//! | two `UpdateDescendants` + static `Descendants` (`cxx:1310-1336`, `cxx:1489-1580`) | `descendants.rs` |

mod descendants;
mod gluing;
mod perform;
mod selection;
mod topo;
mod transform;

pub use descendants::{descendants, update_descendants_bop, update_descendants_gluer};
pub use perform::global_perform;
pub use transform::transform_shape_fu;

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_equal;

/// The two pure virtuals of `BRepFeat_Form` (`hxx:131-133`).
pub trait FormCurves {
    /// `Curves(S)` - the axes of the feature. The reference's derived classes
    /// assign their whole sequence (`BRepFeat_MakePrism.cxx:1010-1013`), so the
    /// port answers the sequence instead of filling one.
    fn curves(&mut self) -> Vec<Option<Arc<dyn Curve>>>;

    /// `BarycCurve()` - the axis through the centre of mass.
    fn baryc_curve(&mut self) -> Option<Arc<dyn Curve>>;
}

/// `BRepFeat_Form` (`hxx:66-180`).
///
/// The protected members of the reference are `pub(crate)` so the derived
/// features (`RibSlot`, `MakePrism`, ...) keep assigning them directly, as the
/// C++ member functions do.
#[derive(Clone)]
pub struct BRepFeatForm {
    /// `myFuse`.
    pub(crate) my_fuse: bool,
    /// `myModify` - read by the derived features, not by `GlobalPerform`.
    #[allow(dead_code)]
    pub(crate) my_modify: bool,
    /// `myMap` - the descendants of every recorded sub-shape.
    pub(crate) my_map: ShapeValueMap<Vec<TopoShape>>,
    /// `myFShape` - the shape whose descendants are the bottom of the form.
    pub(crate) my_f_shape: Option<TopoShape>,
    /// `myLShape` - the shape whose descendants are the top of the form.
    pub(crate) my_l_shape: Option<TopoShape>,
    /// `myNewEdges`.
    pub(crate) my_new_edges: Vec<TopoShape>,
    /// `myTgtEdges`.
    pub(crate) my_tgt_edges: Vec<TopoShape>,
    /// `myPerfSelection`.
    pub(crate) my_perf_selection: PerfSelection,
    /// `myJustGluer`.
    pub(crate) my_just_gluer: bool,
    /// `myJustFeat`.
    pub(crate) my_just_feat: bool,
    /// `mySbase`.
    pub(crate) my_sbase: Option<TopoShape>,
    /// `mySkface` - `SketchFace()`, read by the derived features.
    #[allow(dead_code)]
    pub(crate) my_skface: Option<Face>,
    /// `myGShape`.
    pub(crate) my_gshape: Option<TopoShape>,
    /// `mySFrom`.
    pub(crate) my_sfrom: Option<TopoShape>,
    /// `mySUntil`.
    pub(crate) my_suntil: Option<TopoShape>,
    /// `myGluedF`.
    pub(crate) my_glued_f: ShapeValueMap<TopoShape>,

    /// `mySbOK`.
    my_sb_ok: bool,
    /// `mySkOK`.
    my_sk_ok: bool,
    /// `myGSOK`.
    my_gs_ok: bool,
    /// `mySFOK`.
    my_sf_ok: bool,
    /// `mySUOK`.
    my_su_ok: bool,
    /// `myGFOK`.
    my_gf_ok: bool,
    /// `myPSOK`.
    my_ps_ok: bool,
    /// `myStatusError`.
    my_status_error: StatusError,

    /// `BRepBuilderAPI_MakeShape::myShape`.
    pub(crate) my_shape: Option<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myDone`.
    my_done: bool,
    /// `BRepBuilderAPI_MakeShape::myGenerated` - scratch list of `Modified` /
    /// `Generated`.
    my_generated: Vec<TopoShape>,
}

impl Default for BRepFeatForm {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFeatForm {
    /// `BRepFeat_Form()` (`lxx:19-35`).
    pub fn new() -> Self {
        Self {
            my_fuse: false,
            my_modify: false,
            my_map: ShapeValueMap::new(),
            my_f_shape: None,
            my_l_shape: None,
            my_new_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_perf_selection: PerfSelection::NoSelection,
            my_just_gluer: false,
            my_just_feat: false,
            my_sbase: None,
            my_skface: None,
            my_gshape: None,
            my_sfrom: None,
            my_suntil: None,
            my_glued_f: ShapeValueMap::new(),
            my_sb_ok: false,
            my_sk_ok: false,
            my_gs_ok: false,
            my_sf_ok: false,
            my_su_ok: false,
            my_gf_ok: false,
            my_ps_ok: false,
            my_status_error: StatusError::NotInitialized,
            my_shape: None,
            my_done: false,
            my_generated: Vec::new(),
        }
    }

    /// `BasisShapeValid()` (`lxx:39-42`).
    pub fn basis_shape_valid(&mut self) {
        self.my_sb_ok = true;
    }

    /// `GeneratedShapeValid()` (`lxx:54-57`).
    pub fn generated_shape_valid(&mut self) {
        self.my_gs_ok = true;
    }

    /// `ShapeFromValid()` (`lxx:61-64`).
    pub fn shape_from_valid(&mut self) {
        self.my_sf_ok = true;
    }

    /// `ShapeUntilValid()` (`lxx:68-71`).
    pub fn shape_until_valid(&mut self) {
        self.my_su_ok = true;
    }

    /// `GluedFacesValid()` (`lxx:75-78`).
    pub fn glued_faces_valid(&mut self) {
        self.my_gf_ok = true;
    }

    /// `SketchFaceValid()` (`lxx:82-85`).
    pub fn sketch_face_valid(&mut self) {
        self.my_sk_ok = true;
    }

    /// `PerfSelectionValid()` (`lxx:45-48`).
    pub fn perf_selection_valid(&mut self) {
        self.my_ps_ok = true;
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

    /// `CurrentStatusError()` (`cxx:1483-1486`).
    pub fn current_status_error(&self) -> StatusError {
        self.my_status_error
    }

    /// `IsDeleted(S)` (`cxx:1243-1251`).
    pub fn is_deleted(&self, s: &TopoShape) -> bool {
        match self.my_map.find(s) {
            Some(list) => list.is_empty(),
            None => false,
        }
    }

    /// `Modified(F)` (`cxx:1254-1281`).
    pub fn modified(&mut self, f: &TopoShape) -> &[TopoShape] {
        self.my_generated.clear();
        if !self.my_done {
            return &self.my_generated;
        }
        if let Some(base) = &self.my_sbase {
            if is_equal(base, f) {
                if let Some(sh) = &self.my_shape {
                    self.my_generated.push(sh.clone());
                }
                return &self.my_generated;
            }
        }
        if let Some(list) = self.my_map.find(f) {
            let list = list.clone();
            for sh in list {
                if !crate::topo_tools::is_same(&sh, f) && sh.shape_type() == f.shape_type() {
                    self.my_generated.push(sh);
                }
            }
        }
        &self.my_generated
    }

    /// `Generated(S)` (`cxx:1285-1305`).
    ///
    /// A `FACE` key answers an empty list: the reference filters it out before
    /// looking the shape up (`cxx:1292`).
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.my_generated.clear();
        if !self.my_done {
            return &self.my_generated;
        }
        if s.shape_type() != ShapeType::Face {
            if let Some(list) = self.my_map.find(s) {
                let list = list.clone();
                for sh in list {
                    if !crate::topo_tools::is_same(&sh, s) {
                        self.my_generated.push(sh);
                    }
                }
            }
        }
        &self.my_generated
    }

    /// `FirstShape()` (`cxx:1339-1346`).
    pub fn first_shape(&self) -> &[TopoShape] {
        match &self.my_f_shape {
            Some(s) => self.my_map.find(s).map(|v| v.as_slice()).unwrap_or(&[]),
            None => &self.my_generated,
        }
    }

    /// `LastShape()` (`cxx:1350-1357`).
    pub fn last_shape(&self) -> &[TopoShape] {
        match &self.my_l_shape {
            Some(s) => self.my_map.find(s).map(|v| v.as_slice()).unwrap_or(&[]),
            None => &self.my_generated,
        }
    }

    /// `NewEdges()` (`cxx:1361-1364`).
    pub fn new_edges(&self) -> &[TopoShape] {
        &self.my_new_edges
    }

    /// `TgtEdges()` (`cxx:1368-1371`).
    pub fn tgt_edges(&self) -> &[TopoShape] {
        &self.my_tgt_edges
    }

    /// `GlobalPerform()` (`cxx:59-1240`).
    pub fn global_perform(&mut self, curves: &mut dyn FormCurves) -> Result<(), String> {
        perform::global_perform(self, curves)
    }

    /// `TransformShapeFU(flag)` (`cxx:1378-1480`).
    pub fn transform_shape_fu(&mut self, flag: i32) -> bool {
        transform::transform_shape_fu(self, flag)
    }
}
