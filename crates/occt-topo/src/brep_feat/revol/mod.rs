//! `BRepFeat_MakeRevol` - the revolved feature.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevol.hxx:1-105`,
//! `.cxx:1-991`, `.lxx:1-30`.
//!
//! The class is `BRepFeat_Form` plus one revolution: the primitive is a
//! [`crate::locope::Revol`] over `myPbase` around `myAxis`, its first/last wire
//! are recorded in `myMap` by `statics::maj_map`, and the feature either glues,
//! fuses, cuts or answers the bare primitive.
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | ctor / `Init` (`cxx:71-141`) | `MakeRevol::init` |
//! | `Add(E, F)` (`cxx:147-196`) | `MakeRevol::add` |
//! | `Perform(Angle)` (`cxx:200-336`) | `perform_angle.rs` |
//! | `Perform(Until)` (`cxx:340-467`) | `perform_until.rs` |
//! | `Perform(From, Until)` (`cxx:473-648`) | `perform_from_until.rs` |
//! | `PerformThruAll` (`cxx:655-667`) | `perform_until.rs` |
//! | `PerformUntilAngle` (`cxx:670-792`) | `perform_angle.rs` |
//! | `Curves` / `BarycCurve` (`cxx:799-812`) | the `FormCurves` impl |
//! | file-local `MajMap` / `VerifGluedFaces` / `ToFuse` (`cxx:820-990`) | `statics.rs` |
//!
//! Deviations, all forced by the language:
//! - the five `Perform` overloads are the distinct methods `perform_angle`,
//!   `perform_until`, `perform_from_until`, `perform_thru_all` and
//!   `perform_until_angle`.
//! - the `C == 1` / `C == 0` comparisons on the `bool` `myFuse` of the reference
//!   (`cxx:310`, `:318`, `:448`, `:457`, `:774`, `:778`) are the two branches
//!   of the port; the third `else` of the reference is unreachable for a `bool`
//!   and is not written.
//! - `Standard_ConstructionError` raises come back as `Err` with that name, and
//!   the `StdFail_NotDone` of `Shape()` on an unfinished `LocOpe_Revol` likewise.
//! - `myStatusError` (`hxx:102`) shadows the protected field of the same name in
//!   `BRepFeat_Form` (`Form.hxx:176`); the shadowing is kept.
//! - the `OCCT_DEBUG` `BRepFeat_GettraceFEAT()` traces are a DRAW diagnostic
//!   stream and are not reproduced.

mod perform_angle;
mod perform_from_until;
mod perform_until;
mod statics;

use std::sync::Arc;

use occt_core::gp::GpAx1;
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::form::{BRepFeatForm, FormCurves};
use crate::brep_feat::statics::tool;
use crate::brep_feat::status_error::StatusError;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::{CSIntersector, Revol};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFeat_MakeRevol` (`hxx:54-102`).
#[derive(Clone)]
pub struct MakeRevol {
    /// The inherited `BRepFeat_Form` half of the object.
    pub(crate) form: BRepFeatForm,
    /// `myPbase` - the profile that is revolved.
    pub(crate) my_pbase: Option<TopoShape>,
    /// `myAxis` - the revolution axis.
    my_axis: GpAx1,
    /// `mySlface` - the sliding edges bound to their faces.
    my_slface: ShapeValueMap<Vec<TopoShape>>,
    /// `myCurves` - the axes of the feature, filled by `LocOpe_Revol::Curves`.
    pub(crate) my_curves: Vec<Option<Arc<dyn Curve>>>,
    /// `myBCurve` - the axis through the centre of mass.
    pub(crate) my_bcurve: Option<Arc<dyn Curve>>,
    /// `myStatusError` (`hxx:102`) - shadows the base field, see the module docs.
    #[allow(dead_code)]
    my_status_error: StatusError,
}

impl Default for MakeRevol {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeRevol {
    /// `BRepFeat_MakeRevol()` (`hxx:61`): `gp_Ax1` default-constructed, every
    /// container empty, `myStatusError = BRepFeat_NotInitialized`.
    pub fn new() -> Self {
        Self {
            form: BRepFeatForm::new(),
            my_pbase: None,
            my_axis: GpAx1::default(),
            my_slface: ShapeValueMap::new(),
            my_curves: Vec::new(),
            my_bcurve: None,
            my_status_error: StatusError::NotInitialized,
        }
    }

    /// `BRepFeat_MakeRevol(Sbase, Pbase, Skface, Axis, Fuse, Modify)`
    /// (`hxx:78-86`): construct, then `Init`.
    pub fn with_shapes(
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        axis: &GpAx1,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init(sbase, pbase, skface, axis, fuse, modify);
        this
    }

    /// `Init(Sbase, Pbase, Skface, Axis, Mode, Modify)` (`cxx:71-141`).
    ///
    /// `Mode` selects the operation: `0` cuts, `1` fuses, `2` answers the bare
    /// primitive (`myJustFeat`); any other value leaves the fuse flags alone.
    pub fn init(
        &mut self,
        sbase: &TopoShape,
        pbase: &TopoShape,
        skface: &Face,
        axis: &GpAx1,
        mode: i32,
        modify: bool,
    ) {
        // `cxx:83-91`.
        self.my_axis = *axis;
        self.my_pbase = Some(pbase.clone());
        self.form.my_sbase = Some(sbase.clone());
        self.form.basis_shape_valid();
        self.form.my_skface = Some(skface.clone());
        self.form.sketch_face_valid();
        self.my_slface.clear();

        // `cxx:92-108`.
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

        // `cxx:109-117`. The reference's `mySkface.Nullify()` is commented out
        // (`cxx:113-115`), so the sketch face set above stays.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        self.form.my_shape = None;
        self.form.my_map.clear();
        self.form.my_f_shape = None;
        self.form.my_l_shape = None;

        // `cxx:120-127`: every face of the basis is its own descendant.
        for f in explorer_shapes(sbase, ShapeType::Face) {
            self.form.my_map.bind(&f, vec![f.clone()]);
        }
    }

    /// `Add(E, F)` (`cxx:147-196`): binds the sliding edge `e` to the face `f`.
    ///
    /// `Err("Standard_ConstructionError")` when `f` is not a face of `mySbase`
    /// or `e` is not an edge of `myPbase` (`cxx:155-176`).
    pub fn add(&mut self, e: &TopoShape, f: &Face) -> Result<(), String> {
        // `cxx:155-166`.
        let found_f = match &self.form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face)
                .iter()
                .any(|x| is_same(x, &f.0)),
            None => false,
        };
        if !found_f {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:168-176`.
        let found_e = match &self.my_pbase {
            Some(pbase) => explorer_shapes(pbase, ShapeType::Edge)
                .iter()
                .any(|x| is_same(x, e)),
            None => false,
        };
        if !found_e {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:178-195`: append once.
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

    /// The revolution axis `myAxis`.
    pub fn axis(&self) -> &GpAx1 {
        &self.my_axis
    }

    /// The sliding faces `mySlface`.
    pub fn sliding_faces(&self) -> &ShapeValueMap<Vec<TopoShape>> {
        &self.my_slface
    }

    /// `myAxis.Reverse()` (`cxx:498`).
    pub(crate) fn reverse_axis(&mut self) {
        self.my_axis.reverse();
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

    /// Stores `theRevol.Curves()` / `theRevol.BarycCurve()` (`cxx:331-332`,
    /// `:407-408`, `:412-413`, `:732-733`, `:739-740`).
    pub(crate) fn take_revol_curves(&mut self, the_revol: &Revol) {
        self.my_curves = the_revol.curves().into_iter().map(Some).collect();
        self.my_bcurve = the_revol
            .baryc_curve()
            .map(|l| Arc::new(l) as Arc<dyn Curve>);
    }

    /// The give-up prelude of the branches that stop early: `NotDone()` plus the
    /// (shadowed) status code.
    pub(crate) fn set_error(&mut self, code: StatusError) {
        self.form.not_done();
        self.my_status_error = code;
    }

    /// `cxx:411-467` (and `cxx:739-789`): cut the revolution by the tool solid
    /// of `Until` and answer the solid that still holds a face of `myPbase`.
    ///
    /// `Ok(None)` is the reference's "no intersection" path, where nothing is
    /// assigned and the shape is left undone.
    pub(crate) fn cut_pick_solid(
        &mut self,
        vrairevol: &TopoShape,
        until: &TopoShape,
    ) -> Result<Option<TopoShape>, String> {
        // `cxx:414-419`.
        let mut asi = CSIntersector::with_shape(until);
        let scur: Vec<Option<Arc<dyn Curve>>> = vec![self.my_bcurve.clone()];
        let _ = asi.perform_curves(&scur);
        if !(asi.is_done() && asi.nb_points(1).map(|n| n >= 1).unwrap_or(false)) {
            return Ok(None);
        }

        // `cxx:421-431`.
        let point = match asi.point(1, 1) {
            Ok(p) => p,
            Err(_) => return Ok(None),
        };
        let or = point.orientation();
        let f_until = point.face().clone();
        let mut tools: Vec<TopoShape> = Vec::new();
        if let Some(s) = tool(until, &f_until, or) {
            tools.push(s.0);
        }
        let comp = TopoBuilder::new().make_compound_of(&tools);

        // `cxx:431-437`.
        let mut tr_p = crate::brep_feat::algo_api::cut(vrairevol, &comp.0);
        if tr_p.build().is_err() {
            return Ok(None);
        }
        let Some(cutsh) = tr_p.shape().cloned() else {
            return Ok(None);
        };

        // `cxx:438-447`: the first solid that keeps a face of `myPbase`.
        let mut picked: Option<TopoShape> = None;
        for solid in explorer_shapes(&cutsh, ShapeType::Solid) {
            let keeps = match &self.my_pbase {
                Some(pbase) => explorer_shapes(&solid, ShapeType::Face)
                    .iter()
                    .any(|f| is_same(f, pbase)),
                None => false,
            };
            if keeps {
                picked = Some(solid);
                break;
            }
        }
        Ok(picked)
    }

    /// `cxx:448-466` / `cxx:774-789`: fuse the picked revolution into the basis,
    /// cut it out of the basis, or keep it as the answer of a `myJustFeat`
    /// feature.
    ///
    /// The reference's third `else` (`myFuse == 2`) is unreachable for the
    /// `bool` member and is not written; `myJustFeat` is not consulted on this
    /// path either, exactly as in `Perform(Angle)`.
    pub(crate) fn finish_fuse_or_cut(&mut self, vrairevol: &TopoShape) -> Result<(), String> {
        let sbase = self
            .form
            .my_sbase
            .clone()
            .ok_or_else(|| "Standard_ConstructionError".to_string())?;
        if self.form.my_fuse {
            let (f, shape) = statics::boolean_result(crate::brep_feat::algo_api::fuse(
                &sbase, vrairevol,
            ))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &f, &shape, false)?;
            self.form.done();
        } else {
            let (c, shape) =
                statics::boolean_result(crate::brep_feat::algo_api::cut(&sbase, vrairevol))?;
            self.form.my_shape = Some(shape.clone());
            update_descendants_bop(&mut self.form, &c, &shape, false)?;
            self.form.done();
        }
        Ok(())
    }
}

impl FormCurves for MakeRevol {
    /// `Curves(scur)` (`cxx:799-802`).
    fn curves(&mut self) -> Vec<Option<Arc<dyn Curve>>> {
        self.my_curves.clone()
    }

    /// `BarycCurve()` (`cxx:808-811`).
    fn baryc_curve(&mut self) -> Option<Arc<dyn Curve>> {
        self.my_bcurve.clone()
    }
}
