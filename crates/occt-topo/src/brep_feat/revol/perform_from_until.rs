//! `BRepFeat_MakeRevol::Perform(From, Until)` (`cxx:474-648`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::clib::in_period;
use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::bop_builder2::BoolOp2;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::{parametric_barycenter, tool};
use crate::brep_feat::status_error::StatusError;
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, Revol};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, faces_of, is_same};

use super::statics;
use super::MakeRevol;

impl MakeRevol {
    /// `Perform(From, Until)` (`cxx:474-648`).
    pub fn perform_from_until(
        &mut self,
        from: &TopoShape,
        until: &TopoShape,
    ) -> Result<(), String> {
        // `cxx:480-483`.
        if from.is_null() || until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:485-507`: the sketch face as a limit means "glue from here", and
        // the reversed axis is the reference's `myAxis.Reverse()` (`cxx:498`).
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(from, &sk.0) {
                self.form.my_just_gluer = true;
                self.perform_until(until)?;
                if self.form.my_just_gluer {
                    return Ok(());
                }
            } else if is_same(until, &sk.0) {
                self.form.my_just_gluer = true;
                self.reverse_axis();
                self.perform_until(from)?;
                if self.form.my_just_gluer {
                    return Ok(());
                }
            }
        }

        // `cxx:508-527`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        if faces_of(from).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        self.form.my_sfrom = Some(from.clone());
        let trff = self.form.transform_shape_fu(0);
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trfu = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        if trfu != trff {
            self.set_error(StatusError::IncTypes);
            return Ok(());
        }

        // `cxx:535-540`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let mut the_revol = Revol::new();
        the_revol.perform(&pbase, self.axis(), 2.0 * PI);
        let vrairevol = the_revol
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        statics::maj_map(
            &pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trff {
            // `cxx:542-552`.
            self.form.my_gshape = Some(vrairevol);
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            self.take_revol_curves(&the_revol);
            self.global_perform()?;
        } else {
            // `cxx:553-647`.
            self.take_revol_curves(&the_revol);
            let scur: Vec<Option<Arc<dyn Curve>>> = vec![self.my_bcurve.clone()];
            let su = self.form.my_suntil.clone();
            let sf = self.form.my_sfrom.clone();
            let (Some(su), Some(sf)) = (su, sf) else {
                return Ok(());
            };
            let mut asi1 = CSIntersector::with_shape(&su);
            let mut asi2 = CSIntersector::with_shape(&sf);
            let _ = asi1.perform_curves(&scur);
            let _ = asi2.perform_curves(&scur);

            // `cxx:567-578`.
            let nb1 = if asi1.is_done() {
                asi1.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb1 < 1 {
                self.set_error(StatusError::NoIntersectU);
                return Ok(());
            }
            let or_u = asi1
                .point(1, 1)
                .ok()
                .map(|p| p.orientation())
                .unwrap_or(Orientation::Forward);
            let f_until = asi1.point(1, 1).ok().map(|p| p.face().clone());
            let pr_u = asi1.point(1, 1).ok().map(|p| p.parameter()).unwrap_or(0.0);

            // `cxx:579-596`.
            let nb2 = if asi2.is_done() {
                asi2.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb2 < 1 {
                self.set_error(StatusError::NoIntersectF);
                return Ok(());
            }
            let pr1 = in_period(
                asi2.point(1, 1).ok().map(|p| p.parameter()).unwrap_or(0.0),
                pr_u - 2.0 * PI,
                pr_u,
            );
            let pr2 = in_period(
                asi2.point(1, nb2).ok().map(|p| p.parameter()).unwrap_or(0.0),
                pr_u - 2.0 * PI,
                pr_u,
            );
            let or_f = or_u.reversed();
            let f_from = asi2.point(1, 1).ok().map(|p| p.face().clone());
            let pr_f = if pr1 > pr2 { pr1 } else { pr2 };

            // `cxx:597-600`.
            if !(pr_u > pr_f) {
                self.set_error(StatusError::IncParameter);
                return Ok(());
            }

            // `cxx:601-617`.
            let mut tools: Vec<TopoShape> = Vec::new();
            if let (Some(sf), Some(f)) = (&self.form.my_sfrom, f_from) {
                if let Some(s) = tool(sf, &f, or_f) {
                    tools.push(s.0);
                }
            }
            if let (Some(su), Some(f)) = (&self.form.my_suntil, f_until) {
                if let Some(s) = tool(su, &f, or_u) {
                    tools.push(s.0);
                }
            }
            let comp = TopoBuilder::new().make_compound_of(&tools);
            let mut tr_p = BooleanOperation::new(BoolOp2::Cut);
            tr_p.set_arguments(&[vrairevol]);
            tr_p.set_tools(&[comp.0]);
            tr_p.build()?;
            let cutsh = tr_p
                .shape()
                .cloned()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;

            // `cxx:618-627`.
            let solids = explorer_shapes(&cutsh, ShapeType::Solid);
            let mut vrairevol2 = solids
                .first()
                .cloned()
                .ok_or_else(|| "Standard_OutOfRange".to_string())?;
            for s in &solids {
                let Some(bcurve) = self.my_bcurve.clone() else {
                    break;
                };
                let pr_cur = parametric_barycenter(s, &bcurve);
                if pr_f <= pr_cur && pr_u >= pr_cur {
                    vrairevol2 = s.clone();
                    break;
                }
            }

            // `cxx:628-646`.
            if self.form.my_fuse && !self.form.my_just_feat {
                let sbase = self
                    .form
                    .my_sbase
                    .clone()
                    .ok_or_else(|| "Standard_ConstructionError".to_string())?;
                let (f, shape) =
                    statics::boolean_result(crate::brep_feat::algo_api::fuse(&sbase, &vrairevol2))?;
                self.form.my_shape = Some(shape.clone());
                crate::brep_feat::form::update_descendants_bop(
                    &mut self.form,
                    &f,
                    &shape,
                    false,
                )?;
                self.form.done();
            } else if !self.form.my_fuse && !self.form.my_just_feat {
                let sbase = self
                    .form
                    .my_sbase
                    .clone()
                    .ok_or_else(|| "Standard_ConstructionError".to_string())?;
                let (c, shape) =
                    statics::boolean_result(crate::brep_feat::algo_api::cut(&sbase, &vrairevol2))?;
                self.form.my_shape = Some(shape.clone());
                crate::brep_feat::form::update_descendants_bop(
                    &mut self.form,
                    &c,
                    &shape,
                    false,
                )?;
                self.form.done();
            } else {
                self.form.my_shape = Some(vrairevol2);
                self.form.done();
            }
        }
        Ok(())
    }
}
