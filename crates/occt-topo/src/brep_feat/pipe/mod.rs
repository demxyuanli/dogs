//! `BRepFeat_MakePipe` - the pipe feature.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePipe.hxx:1-122`,
//! `.cxx:1-397`, `.lxx:19-32`.
//!
//! The class is `BRepFeat_Form` plus one [`crate::locope::Pipe`]: the primitive
//! is the sweep of `myPbase` along `mySpine`, its first/last wire are recorded
//! in `myMap` by the file-local `MajMap`, and the feature either fuses, cuts or
//! answers the bare primitive. When the sketch face produced glued faces the
//! feature hands over to `GlobalPerform` instead, the same split
//! `BRepFeat_MakePrism` uses.
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | ctor / `Init` (`cxx:49-118`) | `MakePipe::new` / `MakePipe::with_shapes` / `MakePipe::init` |
//! | `Add(E, F)` (`cxx:122-171`) | `MakePipe::add` |
//! | `Perform()` (`cxx:177-228`) | `perform.rs` |
//! | `Perform(Until)` (`cxx:234-270`) | `perform.rs` |
//! | `Perform(From, Until)` (`cxx:276-336`) | `perform.rs` |
//! | `Curves` / `BarycCurve` (`cxx:338-352`) | the `FormCurves` impl |
//! | file-local `MajMap` (`cxx:354-397`) | `statics.rs` |
//!
//! Deviations, all forced by the language:
//! - the three `Perform` overloads are `perform` / `perform_until` /
//!   `perform_from_until`.
//! - `Standard_ConstructionError` raises come back as `Err` with that name; a
//!   null `mySpine` / `myPbase` / null `LocOpe_Pipe::Shape()` - all dereferences
//!   or `StdFail_NotDone` in the reference - answer `Err("StdFail_NotDone")`.
//! - `BRepFeat_MakePipe.hxx:113` declares `mySkface` again, shadowing the field
//!   of the same name in `BRepFeat_Form` (`Form.hxx:162`). The shadowing is kept:
//!   `Init` only assigns the derived field, so the base field stays unset and
//!   `GlobalPerform` never sees the sketch face, exactly as in the reference.
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` traces are a DRAW diagnostic
//!   stream and are not reproduced.

mod perform;
mod statics;

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::form::{BRepFeatForm, FormCurves};
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_MakePipe` (`hxx:54-118`).
#[derive(Clone)]
pub struct MakePipe {
    /// The inherited `BRepFeat_Form` half of the object.
    pub(crate) form: BRepFeatForm,
    /// `myPbase` (`hxx:112`) - the profile that is swept.
    my_pbase: Option<TopoShape>,
    /// `mySkface` (`hxx:113`) - shadows the base field, see the module docs.
    my_skface: Option<Face>,
    /// `mySlface` (`hxx:114`) - the sliding edges bound to their faces.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `mySpine` (`hxx:116`) - the generating wire.
    my_spine: Option<Wire>,
    /// `myCurves` (`hxx:117`) - the axes of the feature.
    my_curves: Vec<Option<Arc<dyn Curve>>>,
    /// `myBCurve` (`hxx:118`) - the axis through the centre of mass.
    my_bcurve: Option<Arc<dyn Curve>>,
}

impl Default for MakePipe {
    fn default() -> Self {
        Self::new()
    }
}

impl MakePipe {
    /// `BRepFeat_MakePipe()` (`.lxx:21`): every container empty.
    pub fn new() -> Self {
        Self {
            form: BRepFeatForm::new(),
            my_pbase: None,
            my_skface: None,
            my_slface: ShapeValueMap::new(),
            my_spine: None,
            my_curves: Vec::new(),
            my_bcurve: None,
        }
    }

    /// `BRepFeat_MakePipe(Sbase, Pbase, Skface, Spine, Fuse, Modify)`
    /// (`.lxx:25-32`): construct, then `Init`.
    pub fn with_shapes(
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        spine: &Wire,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init(sbase, pbase, skface, spine, fuse, modify);
        this
    }

    /// `Init(Sbase, Pbase, Skface, Spine, Mode, Modify)` (`cxx:49-118`).
    ///
    /// `Mode` selects the operation: `0` cuts, `1` fuses, `2` answers the bare
    /// primitive (`myJustFeat`); any other value leaves the fuse flags alone.
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        spine: &Wire,
        mode: i32,
        modify: bool,
    ) {
        // `cxx:61-67`.
        self.form.my_sbase = Some(sbase.clone());
        self.form.basis_shape_valid();
        self.my_skface = Some(skface.clone());
        self.form.sketch_face_valid();
        self.my_pbase = Some(pbase.clone());
        self.my_slface.clear();
        self.my_spine = Some(spine.clone());

        // `cxx:68-88`.
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

        // `cxx:89-97`. `mySkface.Nullify()` is commented out in the reference
        // (`cxx:90-92`), so the sketch face set above stays.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        self.form.my_shape = None;
        self.form.my_map.clear();
        self.form.my_f_shape = None;
        self.form.my_l_shape = None;

        // `cxx:96-103`: every face of the basis is its own descendant.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.form.my_map.bind(&f, vec![f.clone()]);
        }
    }

    /// `Add(E, F)` (`cxx:122-171`): binds the sliding edge `e` to the face `f`.
    ///
    /// `Err("Standard_ConstructionError")` when `f` is not a face of `mySbase`
    /// or `e` is not an edge of `myPbase` (`cxx:128-152`).
    pub fn add(&mut self, e: &TopoShape, f: &Face) -> Result<(), String> {
        // `cxx:128-139`.
        let found_f = match &self.form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face)
                .iter()
                .any(|x| is_same(x, &f.0)),
            None => false,
        };
        if !found_f {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:141-152`.
        let found_e = match &self.my_pbase {
            Some(pbase) => explorer_shapes(pbase, ShapeType::Edge)
                .iter()
                .any(|x| is_same(x, e)),
            None => false,
        };
        if !found_e {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:154-170`: append once.
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
    pub fn pbase(&self) -> Option<&TopoShape> {
        self.my_pbase.as_ref()
    }

    /// The generating wire `mySpine`.
    pub fn spine(&self) -> Option<&Wire> {
        self.my_spine.as_ref()
    }

    /// The sketch face `mySkface` of the derived class (`hxx:113`).
    pub fn skface(&self) -> Option<&Face> {
        self.my_skface.as_ref()
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
}

impl FormCurves for MakePipe {
    /// `Curves(S)` (`cxx:338-341`).
    fn curves(&mut self) -> Vec<Option<Arc<dyn Curve>>> {
        self.my_curves.clone()
    }

    /// `BarycCurve()` (`cxx:347-351`).
    fn baryc_curve(&mut self) -> Option<Arc<dyn Curve>> {
        self.my_bcurve.clone()
    }
}
