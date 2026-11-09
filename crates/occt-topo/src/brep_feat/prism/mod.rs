//! `BRepFeat_MakePrism` - the prism feature.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePrism.hxx:1-150`,
//! `.cxx:1-1256`, `.lxx:1-30`.
//!
//! The class is `BRepFeat_Form` plus one prism: the primitive is a
//! [`crate::locope::Prism`] over `myPbase`, its first/last wire are recorded in
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
//! | ctor / `Init` (`cxx:80-148`) | `MakePrism::init` |
//! | `Add(E, F)` (`cxx:155-204`) | `MakePrism::add` |
//! | `Perform(Length)` (`cxx:210-319`) | `perform_length.rs` |
//! | `Perform(Until)` (`cxx:325-433`) | `perform_until.rs` |
//! | `Perform(From, Until)` (`cxx:440-651`) | `perform_from_until.rs` |
//! | `PerformUntilEnd` / `PerformFromEnd` / `PerformThruAll` (`cxx:659-897`) | `perform_end.rs` |
//! | `PerformUntilHeight` (`cxx:903-1005`) | `perform_length.rs` |
//! | `Curves` / `BarycCurve` (`cxx:1010-1024`) | the `FormCurves` impl |
//! | file-local `MajMap` / `ToFuse` / `HeightMax` / `SensOfPrism` / `TestCurve` (`cxx:1136-1256`) | `statics.rs` |
//!
//! Deviations, all forced by the language:
//! - the eight `Perform` overloads are the distinct methods `perform_length`,
//!   `perform_until`, `perform_from_until`, `perform_until_end`,
//!   `perform_from_end`, `perform_thru_all` and `perform_until_height`.
//! - `Standard_ConstructionError` raises come back as `Err` with that name, and
//!   the `StdFail_NotDone` of `Shape()` on an unfinished `LocOpe_Prism` likewise.
//! - `myStatusError` (`hxx:145`) shadows the protected field of the same name in
//!   `BRepFeat_Form` (`Form.hxx:176`), so the codes assigned by the `Perform`
//!   overloads never reach `CurrentStatusError()`; the shadowing is kept.
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` traces are a DRAW diagnostic
//!   stream and are not reproduced.

mod perform_end;
mod perform_from_until;
mod perform_length;
mod perform_until;
mod statics;

use std::sync::Arc;

use occt_core::gp::{GpDir, GpVec};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::form::{BRepFeatForm, FormCurves};
use crate::brep_feat::status_error::StatusError;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;

/// `BRepFeat_MakePrism` (`hxx:54-147`).
#[derive(Clone)]
pub struct MakePrism {
    /// The inherited `BRepFeat_Form` half of the object.
    pub(crate) form: BRepFeatForm,
    /// `myPbase` - the profile that is swept.
    my_pbase: Option<TopoShape>,
    /// `mySlface` - the sliding edges bound to their faces.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `myDir` - the sweep direction.
    my_dir: GpDir,
    /// `myCurves` - the axes of the feature, filled by `LocOpe_Prism::Curves`.
    my_curves: Vec<Option<Arc<dyn Curve>>>,
    /// `myBCurve` - the axis through the centre of mass.
    my_bcurve: Option<Arc<dyn Curve>>,
    /// `myStatusError` (`hxx:145`) - shadows the base field, see the module docs.
    #[allow(dead_code)]
    my_status_error: StatusError,
}

impl Default for MakePrism {
    fn default() -> Self {
        Self::new()
    }
}

impl MakePrism {
    /// `BRepFeat_MakePrism()` (`hxx:61`): `gp_Dir` default-constructed, every
    /// container empty, `myStatusError = BRepFeat_NotInitialized`.
    pub fn new() -> Self {
        Self {
            form: BRepFeatForm::new(),
            my_pbase: None,
            my_slface: ShapeValueMap::new(),
            my_dir: GpDir::default(),
            my_curves: Vec::new(),
            my_bcurve: None,
            my_status_error: StatusError::NotInitialized,
        }
    }

    /// `BRepFeat_MakePrism(Sbase, Pbase, Skface, Direction, Fuse, Modify)`
    /// (`hxx:78-86`): construct, then `Init`.
    pub fn with_shapes(
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        direction: &GpDir,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init(sbase, pbase, skface, direction, fuse, modify);
        this
    }

    /// `Init(Sbase, Pbase, Skface, Direc, Mode, Modify)` (`cxx:80-148`).
    ///
    /// `Mode` selects the operation: `0` cuts, `1` fuses, `2` answers the bare
    /// primitive (`myJustFeat`); any other value leaves the fuse flags alone.
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        direc: &GpDir,
        mode: i32,
        modify: bool,
    ) {
        // `cxx:88-92`.
        self.form.my_skface = Some(skface.clone());
        self.form.sketch_face_valid();
        self.form.my_sbase = Some(sbase.clone());
        self.form.basis_shape_valid();
        self.my_pbase = Some(pbase.clone());
        self.my_slface.clear();
        self.my_dir = *direc;

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

        // `cxx:113-125`. The reference's `mySkface.Nullify()` is commented out
        // (`cxx:117-119`), so the sketch face set above stays.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        self.form.my_shape = None;
        self.form.my_new_edges.clear();
        self.form.my_tgt_edges.clear();
        self.form.my_map.clear();
        self.form.my_f_shape = None;
        self.form.my_l_shape = None;

        // `cxx:127-133`: every face of the basis is its own descendant.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.form.my_map.bind(&f, vec![f.clone()]);
        }
    }

    /// `Add(E, F)` (`cxx:155-204`): binds the sliding edge `e` to the face `f`.
    ///
    /// `Err("Standard_ConstructionError")` when `f` is not a face of `mySbase`
    /// or `e` is not an edge of `myPbase` (`cxx:169-182`).
    pub fn add(&mut self, e: &TopoShape, f: &Face) -> Result<(), String> {
        // `cxx:160-172`.
        let found_f = match &self.form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face)
                .iter()
                .any(|x| crate::topo_tools::is_same(x, &f.0)),
            None => false,
        };
        if !found_f {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:174-183`.
        let found_e = match &self.my_pbase {
            Some(pbase) => explorer_shapes(pbase, ShapeType::Edge)
                .iter()
                .any(|x| crate::topo_tools::is_same(x, e)),
            None => false,
        };
        if !found_e {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:185-203`: append once.
        if !self.my_slface.is_bound(&f.0) {
            self.my_slface.bind(&f.0, Vec::new());
        }
        if let Some(list) = self.my_slface.change(&f.0) {
            if !list.iter().any(|x| crate::topo_tools::is_same(x, e)) {
                list.push(e.clone());
            }
        }
        Ok(())
    }

    /// The profile `myPbase`.
    pub fn pbase(&self) -> Option<&TopoShape> {
        self.my_pbase.as_ref()
    }

    /// The sweep direction `myDir`.
    pub fn dir(&self) -> &GpDir {
        &self.my_dir
    }

    /// The sliding faces `mySlface`.
    pub fn sliding_faces(&self) -> &ShapeValueMap<Vec<TopoShape>> {
        &self.my_slface
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

    /// Stores `thePrism.Curves()` / `thePrism.BarycCurve()` (`cxx:291-292` and
    /// the other seven call sites).
    pub(crate) fn take_prism_curves(&mut self, the_prism: &crate::locope::Prism) {
        self.my_curves = the_prism.curves().into_iter().map(Some).collect();
        self.my_bcurve = the_prism
            .baryc_curve()
            .map(|l| Arc::new(l) as Arc<dyn Curve>);
    }

    /// The give-up prelude of the branches that stop early: `NotDone()` plus the
    /// (shadowed) status code.
    pub(crate) fn set_error(&mut self, code: StatusError) {
        self.form.not_done();
        self.my_status_error = code;
    }

    /// `scalar * myDir` as a `gp_Vec`.
    pub(crate) fn scaled_dir(&self, scalar: f64) -> GpVec {
        statics::scaled_dir(scalar, &self.my_dir)
    }

    /// `myDir` assignment (`myDir.Reverse()` at `cxx:719`).
    pub(crate) fn set_dir(&mut self, d: GpDir) {
        self.my_dir = d;
    }

    /// `cxx:297-318`, `cxx:411-430`, `cxx:632-650`, `cxx:824-842` and
    /// `cxx:984-1001`: fuse `cutsh` into the basis, cut it out of the basis, or
    /// answer `cutsh` itself for a `myJustFeat` feature.
    pub(crate) fn finish_with_boolean_on(&mut self, cutsh: &TopoShape) -> Result<(), String> {
        let sbase = self
            .form
            .my_sbase
            .clone()
            .ok_or_else(|| "Standard_ConstructionError".to_string())?;
        if self.form.my_fuse && !self.form.my_just_feat {
            let (f, shape) = statics::boolean_result(crate::brep_feat::algo_api::fuse(&sbase, cutsh))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &f, &shape, false)?;
            self.form.done();
        } else if !self.form.my_fuse {
            let (c, shape) = statics::boolean_result(crate::brep_feat::algo_api::cut(&sbase, cutsh))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &c, &shape, false)?;
            self.form.done();
        } else {
            self.form.my_shape = Some(cutsh.clone());
            self.form.done();
        }
        Ok(())
    }
}

impl FormCurves for MakePrism {
    /// `Curves(S)` (`cxx:1010-1013`).
    fn curves(&mut self) -> Vec<Option<Arc<dyn Curve>>> {
        self.my_curves.clone()
    }

    /// `BarycCurve()` (`cxx:1020-1023`).
    fn baryc_curve(&mut self) -> Option<Arc<dyn Curve>> {
        self.my_bcurve.clone()
    }
}
