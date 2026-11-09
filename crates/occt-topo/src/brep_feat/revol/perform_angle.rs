//! `BRepFeat_MakeRevol::Perform(Angle)` (`cxx:200-336`) and
//! `PerformUntilAngle` (`cxx:670-792`).

use std::f64::consts::PI;

use occt_core::precision;

use crate::abs::ShapeType;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::locope::Revol;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, faces_of, is_same};

use super::statics;
use super::MakeRevol;

impl MakeRevol {
    /// `Perform(Angle)` (`cxx:200-336`).
    ///
    /// A complete turn drops the sketch face (`cxx:236-238`); the `angledec`
    /// block that would shorten the sweep is commented out in the reference
    /// (`cxx:219-234`), so the revolution is always the plain `Angle`.
    ///
    /// When no face was glued the answer is a Boolean against the generated
    /// shape (`cxx:308-328`); otherwise `GlobalPerform` runs the topological
    /// selection (`cxx:330-335`).
    pub fn perform_angle(&mut self, angle: f64) -> Result<(), String> {
        // `cxx:207-214`.
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();

        // `cxx:215`.
        let revol_comp = 2.0 * PI - angle.abs() <= precision::ANGULAR;

        // `cxx:219-234`: the `angledec` computation is commented out in the
        // reference, so `angledec` stays `0.` and only `mySkface.Nullify()`
        // remains (`cxx:236-238`).
        if revol_comp {
            self.form.my_skface = None;
        }

        // `cxx:239-243`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let mut the_revol = Revol::new();
        the_revol.perform(&pbase, self.axis(), angle);
        let vrairevol = the_revol
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:245-246`.
        statics::maj_map(
            &pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // `cxx:248`.
        self.form.my_gshape = Some(vrairevol);
        self.form.generated_shape_valid();

        // `cxx:249-256`: the first shape must hold exactly one face.
        let first = the_revol
            .first_shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let base_faces = faces_of(&first);
        let the_base = base_faces
            .first()
            .cloned()
            .ok_or_else(|| "Standard_OutOfRange".to_string())?;
        if base_faces.len() > 1 {
            self.set_error(StatusError::InvFirstShape);
            return Ok(());
        }

        // `cxx:258-301`: the search of the last-shape face and its coplanar
        // partner in the basis. Neither `FFace` nor `found` is read after the
        // block; the traversal is kept branch-for-branch.
        if self.form.my_skface.is_some() || !self.my_slface.is_empty() {
            let mut fface: Option<Face> = None;
            if let Some(lshape) = self.form.my_l_shape.clone() {
                if lshape.shape_type() == ShapeType::Wire {
                    if let Some(gshape) = self.form.my_gshape.clone() {
                        'faces: for f in explorer_shapes(&gshape, ShapeType::Face) {
                            for w in explorer_shapes(&f, ShapeType::Wire) {
                                if is_same(&w, &lshape) {
                                    fface = Face::wrap(f.clone());
                                    break 'faces;
                                }
                            }
                        }
                    }
                }
            }
            if let Some(sbase) = self.form.my_sbase.clone() {
                for ff in faces_of(&sbase) {
                    if !statics::to_fuse(Some(&ff), fface.as_ref()) {
                        continue;
                    }
                    // `cxx:295-299`. `ToFuse` is false on a null `FFace`.
                    if let Some(fc) = fface.as_ref() {
                        if !is_same(&fc.0, &pbase) && crate::brep_feat::statics::is_inside(&ff, fc)
                        {
                            break;
                        }
                    }
                }
            }
        }

        // `cxx:302`.
        self.form.glued_faces_valid();

        // `cxx:303-307`.
        if let Some(sk) = self.form.my_skface.clone() {
            statics::verif_glued_faces(self, &sk, &the_base.0, &the_revol)?;
        }

        // `cxx:308-335`.
        if self.form.my_glued_f.is_empty() {
            let sbase = self
                .form
                .my_sbase
                .clone()
                .ok_or_else(|| "Standard_ConstructionError".to_string())?;
            let gshape = self
                .form
                .my_gshape
                .clone()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;
            // `cxx:310-327`. The reference's third `else` (`myFuse == 2`) is
            // unreachable for the `bool` member and is not written.
            if self.form.my_fuse {
                let (f, shape) = statics::boolean_result(crate::brep_feat::algo_api::fuse(
                    &sbase, &gshape,
                ))?;
                self.form.my_shape = Some(shape.clone());
                update_descendants_bop(&mut self.form, &f, &shape, false)?;
                self.form.done();
            } else {
                let (c, shape) =
                    statics::boolean_result(crate::brep_feat::algo_api::cut(&sbase, &gshape))?;
                self.form.my_shape = Some(shape.clone());
                update_descendants_bop(&mut self.form, &c, &shape, false)?;
                self.form.done();
            }
        } else {
            // `cxx:330-335`.
            self.take_revol_curves(&the_revol);
            self.global_perform()?;
        }
        Ok(())
    }

    /// `PerformUntilAngle(Until, Angle)` (`cxx:670-792`).
    ///
    /// The reference does not return after the `Until.IsNull()` /
    /// `Angle == 0` redirections (`cxx:678-685`), so both fall through into the
    /// rest of the body; that is reproduced here.
    pub fn perform_until_angle(
        &mut self,
        until: &TopoShape,
        angle: f64,
    ) -> Result<(), String> {
        // `cxx:678-685`.
        if until.is_null() {
            self.perform_angle(angle)?;
        }
        if angle == 0.0 {
            self.perform_until(until)?;
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:690-694`.
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(until, &sk.0) {
                return self.perform_angle(angle);
            }
        }

        // `cxx:695-703`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:706-710`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let mut the_revol = Revol::new();
        the_revol.perform(&pbase, self.axis(), angle);
        let vrairevol = the_revol
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:710`.
        statics::maj_map(
            &pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if trf {
            // `cxx:712-736`.
            self.form.my_gshape = Some(vrairevol);
            self.form.generated_shape_valid();

            let first = the_revol
                .first_shape()
                .cloned()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;
            let base_faces = faces_of(&first);
            base_faces
                .first()
                .cloned()
                .ok_or_else(|| "Standard_OutOfRange".to_string())?;
            if base_faces.len() > 1 {
                self.set_error(StatusError::InvFirstShape);
                return Ok(());
            }

            // `cxx:727`.
            self.form.glued_faces_valid();
            // `cxx:728`: `VerifGluedFaces` is commented out in the reference.
            self.take_revol_curves(&the_revol);
            self.global_perform()?;
        } else {
            // `cxx:738-790`.
            self.take_revol_curves(&the_revol);
            if let Some(picked) = self.cut_pick_solid(&vrairevol, until)? {
                self.finish_fuse_or_cut(&picked)?;
            }
        }
        Ok(())
    }
}
