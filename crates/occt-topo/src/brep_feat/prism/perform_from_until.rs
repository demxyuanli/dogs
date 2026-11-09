//! `BRepFeat_MakePrism::Perform(From, Until)` (`cxx:440-651`).
//!
//! The prism is swept far enough to cover both limits, with an optional
//! pre-translation when the two limits lie on opposite sides (`tran < 0`). When
//! both supports are proper faces the regular `GlobalPerform` runs; otherwise
//! both limits are turned into tool solids and cut from the primitive.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::Orientation;
use crate::bop_builder2::BoolOp2;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::brep_feat::statics::tool;
use crate::locope::{CSIntersector, Prism};
use crate::shape::TopoShape;
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakePrism;

impl MakePrism {
    /// `Perform(From, Until)` (`cxx:440-651`).
    pub fn perform_from_until(
        &mut self,
        from: &TopoShape,
        until: &TopoShape,
    ) -> Result<(), String> {
        // `cxx:448-451`.
        if from.is_null() || until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:453-474`: a limit that is the sketch face means "glue from here",
        // which is the single-limit overload.
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(from, &sk.0) {
                self.form.my_just_gluer = true;
                self.perform_until(until)?;
                if self.form.my_just_gluer {
                    return Ok(());
                }
            } else if is_same(until, &sk.0) {
                self.form.my_just_gluer = true;
                self.perform_until(from)?;
                if self.form.my_just_gluer {
                    return Ok(());
                }
            }
        }

        // `cxx:476-485`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        if faces_of(from).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:487-499`.
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

        // `cxx:501-518`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let c = statics::test_curve(&pbase, self.dir());
        let (sens, tran) = if is_same(from, until) {
            (1, -1)
        } else {
            let sens = match &c {
                Some(c) => statics::sens_of_prism(c, until),
                None => 1,
            };
            let sens_from = match (&c, &self.form.my_sfrom) {
                (Some(c), Some(sf)) => statics::sens_of_prism(c, sf),
                _ => 1,
            };
            (sens, sens * sens_from)
        };

        // `cxx:518-527`.
        let the_prism = if tran < 0 {
            let v = self.scaled_dir(3.0 * sens as f64 * height);
            let vtra = self.scaled_dir(-3.0 * height * sens as f64 / 2.0);
            Prism::with_base_and_translation(&pbase, &v, &vtra)
        } else {
            let v = self.scaled_dir(2.0 * sens as f64 * height);
            Prism::with_base(&pbase, &v)
        };
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        if !trff {
            // `cxx:529-537`.
            statics::maj_map(
                &pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            self.take_prism_curves(&the_prism);
            self.global_perform()?;
        } else {
            // `cxx:540-650`.
            statics::maj_map(
                &pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur: Vec<Option<Arc<dyn Curve>>> = match &c {
                Some(c) => vec![Some(c.clone())],
                None => vec![None],
            };
            let mut asi1 = match &self.form.my_suntil {
                Some(su) => CSIntersector::with_shape(su),
                None => return Ok(()),
            };
            let mut asi2 = match &self.form.my_sfrom {
                Some(sf) => CSIntersector::with_shape(sf),
                None => return Ok(()),
            };
            let _ = asi1.perform_curves(&scur);
            let _ = asi2.perform_curves(&scur);

            // `cxx:550-570`.
            let nb1 = if asi1.is_done() {
                asi1.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb1 < 1 {
                self.set_error(StatusError::NoIntersectU);
                return Ok(());
            }
            let mut or_u = if self.form.my_fuse {
                asi1.point(1, 1).ok().map(|p| p.orientation())
            } else {
                asi1.point(1, nb1).ok().map(|p| p.orientation())
            }
            .unwrap_or(Orientation::Forward);
            if sens == -1 {
                or_u = or_u.reversed();
            }
            let f_until = asi1.point(1, 1).ok().map(|p| p.face().clone());
            let par_u = asi1.point(1, 1).ok().map(|p| p.parameter()).unwrap_or(0.0);

            // `cxx:576-592`.
            let nb2 = if asi2.is_done() {
                asi2.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb2 < 1 {
                self.set_error(StatusError::NoIntersectF);
                return Ok(());
            }
            let mut or_f = asi2
                .point(1, 1)
                .map(|p| p.orientation())
                .unwrap_or(Orientation::Forward);
            if sens == 1 {
                or_f = or_f.reversed();
            }
            let f_from = asi2.point(1, 1).ok().map(|p| p.face().clone());
            let par_f = asi2.point(1, 1).ok().map(|p| p.parameter()).unwrap_or(0.0);

            // `cxx:594-598`: the orientations follow the limit closer to the
            // origin of the sweep when the prism was transferred.
            if tran > 0 && par_u.abs() < par_f.abs() {
                std::mem::swap(&mut or_u, &mut or_f);
            }

            // `cxx:601-621`.
            let mut tools: Vec<TopoShape> = Vec::new();
            let tool_u = match (&self.form.my_suntil, f_until) {
                (Some(su), Some(f)) => tool(su, &f, or_u),
                _ => None,
            };
            match tool_u {
                Some(s) => tools.push(s.0),
                None => {
                    self.set_error(StatusError::NullToolU);
                    return Ok(());
                }
            }
            let tool_f = match (&self.form.my_sfrom, f_from) {
                (Some(sf), Some(f)) => tool(sf, &f, or_f),
                _ => None,
            };
            match tool_f {
                Some(s) => tools.push(s.0),
                None => {
                    self.set_error(StatusError::NullToolF);
                    return Ok(());
                }
            }

            // `cxx:623-631`: `BRepAlgoAPI_Cut trP; SetArguments / SetTools / Build`.
            let mut tr_p = BooleanOperation::new(BoolOp2::Cut);
            tr_p.set_arguments(&[vraiprism.clone()]);
            tr_p.set_tools(&tools);
            let (tr_p, tr_shape) = statics::boolean_result(tr_p)?;
            update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;

            // `cxx:632-650`.
            self.finish_with_boolean_on(&tr_shape)?;
        }
        Ok(())
    }
}
