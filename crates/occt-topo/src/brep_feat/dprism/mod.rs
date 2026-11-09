//! `BRepFeat_MakeDPrism` -- the draft prism feature.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeDPrism.hxx:1-154`,
//! `.cxx:1-1315`.
//!
//! The class is `BRepFeat_Form` plus one draft prism: the face `myPbase` is
//! swept by a [`crate::locope::DPrism`] along the three-segment profile built
//! from `myAngle`, the first / the last wire of the sweep are recorded in
//! `myMap` by `statics::maj_map`, and the feature either glues, fuses, cuts or
//! answers the bare primitive. Every `Perform` overload first decides how long
//! and in which direction to sweep (`HeightMax` / `SensOfPrism` / `TestCurve`),
//! then hands over to `GlobalPerform` when the limit is a face of `mySbase`, or
//! to a plain Boolean when it is not.
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | ctor / `Init` (`cxx:76-148`) | `MakeDPrism::init` |
//! | `Add(E, F)` (`cxx:152-201`) | `MakeDPrism::add` |
//! | `Perform(Height)` (`cxx:205-276`) | `perform_height.rs` |
//! | `Perform(Until)` (`cxx:283-414`) | `perform_until.rs` |
//! | `Perform(From, Until)` (`cxx:418-617`) | `perform_from_until.rs` |
//! | `PerformUntilEnd` / `PerformFromEnd` / `PerformThruAll` (`cxx:620-875`) | `perform_end.rs` |
//! | `PerformUntilHeight` (`cxx:880-997`) | `perform_height.rs` |
//! | `Curves` / `BarycCurve` (`cxx:1004-1008`, `cxx:1152-1156`) | the `FormCurves` impl |
//! | `BossEdges` / `TopEdges` / `LatEdges` (`cxx:1015-1148`) | `boss_edges.rs` |
//! | file-local `MajMap` / `HeightMax` / `SensOfPrism` / `TestCurve` (`cxx:1159-1313`) | `statics.rs` |
//!
//! Deviations, all forced by the language:
//! - the `Perform` overloads are the distinct methods `perform_height`,
//!   `perform_until`, `perform_from_until`, `perform_until_end`,
//!   `perform_from_end`, `perform_thru_all` and `perform_until_height`.
//! - `Standard_ConstructionError` raises come back as `Err` with that name, and
//!   the `StdFail_NotDone` of `LocOpe_DPrism::Shape()` likewise.
//! - `myStatusError` (`hxx:155`) shadows the protected field of the same name in
//!   `BRepFeat_Form` (`Form.hxx:176`), so the codes assigned by the `Perform`
//!   overloads never reach `CurrentStatusError()`; the shadowing is kept.
//! - the default constructor's `myAngle(RealLast())` is `f64::MAX`.
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` traces are a DRAW diagnostic
//!   stream and are not reproduced.

mod boss_edges;
mod perform_end;
mod perform_from_until;
mod perform_height;
mod perform_until;
mod statics;

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::form::{BRepFeatForm, FormCurves};
use crate::brep_feat::status_error::StatusError;
use crate::locope::dprism::DPrism;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_MakeDPrism` (`hxx:54-154`).
#[derive(Clone)]
pub struct MakeDPrism {
    /// The inherited `BRepFeat_Form` half of the object.
    pub(crate) form: BRepFeatForm,
    /// `myPbase` -- the face that is swept.
    my_pbase: Option<Face>,
    /// `mySlface` -- the sliding edges bound to their faces.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `myAngle` -- the draft angle.
    my_angle: f64,
    /// `myCurves` -- the axes of the feature, filled by `LocOpe_DPrism::Curves`.
    my_curves: Vec<Option<Arc<dyn Curve>>>,
    /// `myBCurve` -- the axis through the centre of mass.
    my_bcurve: Option<Arc<dyn Curve>>,
    /// `myTopEdges` (`hxx:152`).
    my_top_edges: Vec<TopoShape>,
    /// `myLatEdges` (`hxx:153`).
    my_lat_edges: Vec<TopoShape>,
    /// `myStatusError` (`hxx:155`) -- shadows the base field, see the docs.
    #[allow(dead_code)]
    my_status_error: StatusError,
}

impl Default for MakeDPrism {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeDPrism {
    /// `BRepFeat_MakeDPrism()` (`hxx:77-79`): `myAngle(RealLast())` and
    /// `myStatusError(BRepFeat_OK)`, every container empty.
    pub fn new() -> Self {
        Self {
            form: BRepFeatForm::new(),
            my_pbase: None,
            my_slface: ShapeValueMap::new(),
            // `RealLast()`.
            my_angle: f64::MAX,
            my_curves: Vec::new(),
            my_bcurve: None,
            my_top_edges: Vec::new(),
            my_lat_edges: Vec::new(),
            my_status_error: StatusError::Ok,
        }
    }

    /// `BRepFeat_MakeDPrism(Sbase, Pbase, Skface, Angle, Fuse, Modify)`
    /// (`hxx:66-73`): construct, then `Init`.
    pub fn with_shapes(
        sbase: &TopoShape,
        pbase: &Face,
        skface: &Face,
        angle: f64,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init(sbase, pbase, skface, angle, fuse, modify);
        this
    }

    /// `Init(Sbase, Pbase, Skface, Angle, Mode, Modify)` (`cxx:76-148`).
    ///
    /// `Mode` selects the operation: `0` cuts, `1` fuses, `2` answers the bare
    /// primitive (`myJustFeat`); any other value leaves the fuse flags alone.
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        pbase: &Face,
        skface: &Face,
        angle: f64,
        mode: i32,
        modify: bool,
    ) {
        // `cxx:89-94`.
        self.form.my_skface = Some(skface.clone());
        self.form.sketch_face_valid();
        self.form.my_sbase = Some(sbase.clone());
        self.form.basis_shape_valid();
        self.my_pbase = Some(pbase.clone());
        self.my_slface.clear();

        // `cxx:96-112`.
        if mode == 0 {
            self.form.my_fuse = false;
            self.form.my_just_feat = false;
        } else if mode == 1 {
            self.form.my_fuse = true;
            self.form.my_just_feat = false;
        } else if mode == 2 {
            self.form.my_fuse = true;
            self.form.my_just_feat = true;
        }

        // `cxx:113-124`. The reference's `mySkface.Nullify()` is commented out
        // (`cxx:117-119`), so the sketch face set above stays.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        self.form.my_shape = None;
        self.form.my_map.clear();
        self.form.my_f_shape = None;
        self.form.my_l_shape = None;
        self.my_top_edges.clear();
        self.my_lat_edges.clear();

        // `cxx:126-133`: every face of the basis is its own descendant.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.form.my_map.bind(&f, vec![f.clone()]);
        }

        // `cxx:134`.
        self.my_angle = angle;
    }

    /// `Add(E, F)` (`cxx:152-201`): binds the sliding edge `e` to the face `f`.
    ///
    /// `Err("Standard_ConstructionError")` when `f` is not a face of `mySbase`
    /// or `e` is not an edge of `myPbase` (`cxx:166-182`).
    pub fn add(&mut self, e: &TopoShape, f: &Face) -> Result<(), String> {
        // `cxx:159-171`.
        let found_f = match &self.form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face)
                .iter()
                .any(|x| is_same(x, &f.0)),
            None => false,
        };
        if !found_f {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:173-183`.
        let found_e = match &self.my_pbase {
            Some(pbase) => explorer_shapes(&pbase.0, ShapeType::Edge)
                .iter()
                .any(|x| is_same(x, e)),
            None => false,
        };
        if !found_e {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:185-200`: append once.
        if !self.my_slface.is_bound(&f.0) {
            self.my_slface.bind(&f.0, Vec::new());
        }
        if let Some(list) = self.my_slface.change(&f.0) {
            if !list.iter().any(|x| is_same(x, e)) {
                list.push(e.clone());
            }
        }
        Ok(())
    }

    /// The profile `myPbase`.
    pub fn pbase(&self) -> Option<&Face> {
        self.my_pbase.as_ref()
    }

    /// The draft angle `myAngle`.
    pub fn angle(&self) -> f64 {
        self.my_angle
    }

    /// The sliding faces `mySlface`.
    pub fn sliding_faces(&self) -> &ShapeValueMap<Vec<TopoShape>> {
        &self.my_slface
    }

    /// `myPbase` as a shape, the profile of the `Perform` overloads.
    pub(crate) fn pbase_face(&self) -> Result<Face, String> {
        self.my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())
    }

    /// `GlobalPerform()` on the inherited half, with `self` as the curve source.
    ///
    /// The form has to leave `self` for the duration of the call because
    /// `GlobalPerform` asks the same object for `Curves` / `BarycCurve`; the
    /// value is restored on both the `Ok` and the `Err` path.
    pub(crate) fn global_perform(&mut self) -> Result<(), String> {
        let mut form = std::mem::take(&mut self.form);
        let res = crate::brep_feat::form::global_perform(&mut form, self);
        self.form = form;
        res
    }

    /// Stores `theDPrism.Curves()` / `theDPrism.BarycCurve()` (`cxx:248-249`
    /// and the other six call sites).
    pub(crate) fn take_dprism_curves(&mut self, the_dprism: &DPrism) {
        self.my_curves = the_dprism.curves().into_iter().map(Some).collect();
        self.my_bcurve = the_dprism
            .baryc_curve()
            .map(|l| Arc::new(l) as Arc<dyn Curve>);
    }

    /// The give-up prelude of the branches that stop early: `NotDone()` plus the
    /// (shadowed) status code.
    pub(crate) fn set_error(&mut self, code: StatusError) {
        self.form.not_done();
        self.my_status_error = code;
    }

    /// `cxx:254-274`, `cxx:383-400`, `cxx:596-615`, `cxx:792-814`,
    /// `cxx:857-864` and `cxx:975-993`: fuse `cutsh` into the basis or cut it
    /// out of the basis.
    ///
    /// Every one of these sites tests `myFuse == 1` / `myFuse == 0` and keeps a
    /// trailing `else { myShape = cutsh; Done(); }`. `myFuse` is a boolean, so
    /// that trailing branch needs a value that is neither and cannot be reached;
    /// the reference never reads `myJustFeat` in this class.
    pub(crate) fn finish_with_boolean_on(&mut self, cutsh: &TopoShape) -> Result<(), String> {
        let sbase = self
            .form
            .my_sbase
            .clone()
            .ok_or_else(|| "Standard_ConstructionError".to_string())?;
        if self.form.my_fuse {
            let (f, shape) = statics::boolean_result(crate::brep_feat::algo_api::fuse(&sbase, cutsh))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &f, &shape, false)?;
            self.form.done();
        } else {
            let (c, shape) = statics::boolean_result(crate::brep_feat::algo_api::cut(&sbase, cutsh))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &c, &shape, false)?;
            self.form.done();
        }
        Ok(())
    }
}

impl FormCurves for MakeDPrism {
    /// `Curves(S)` (`cxx:1004-1008`).
    fn curves(&mut self) -> Vec<Option<Arc<dyn Curve>>> {
        self.my_curves.clone()
    }

    /// `BarycCurve()` (`cxx:1152-1156`).
    fn baryc_curve(&mut self) -> Option<Arc<dyn Curve>> {
        self.my_bcurve.clone()
    }
}
