//! `BRepFeat_MakeDPrism::Perform(Height)` (`cxx:205-276`) and
//! `PerformUntilHeight(Until, Height)` (`cxx:880-997`).
//!
//! Both sweep the profile by a known height: `Perform(Height)` divides it by
//! `cos(myAngle)` so the draft side reaches `Height` along the spine normal,
//! `PerformUntilHeight` multiplies it by the direction `SensOfPrism` finds.
//! When the top of the sweep is a single face they go through `GlobalPerform`
//! (gluer) or a plain Boolean on the primitive.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::brep_feat::statics::tool;
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, DPrism};
use crate::shape::TopoShape;
use crate::topo_tools::faces_of;

use super::statics;
use super::MakeDPrism;

impl MakeDPrism {
    /// `Perform(Height)` (`cxx:205-276`): a draft prism of `Height / cos(angle)`
    /// along the spine, then the glue-or-Boolean decision (`cxx:254-275`).
    pub fn perform_height(&mut self, height: f64) -> Result<(), String> {
        // `cxx:212-218`.
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = None;
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();

        // `cxx:220-223`.
        let pbase = self.pbase_face()?;
        let theheight = height / self.angle().cos();
        let the_dprism = DPrism::with_height(&pbase, theheight, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:225-228`.
        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = Some(vraiprism.clone());
        self.form.generated_shape_valid();

        // `cxx:230-238`: the bottom of the sweep has to be a single face.
        if !statics::first_shape_has_one_face(&the_dprism)? {
            self.set_error(StatusError::InvFirstShape);
            return Ok(());
        }

        // `cxx:242-252`: the gluing route.
        self.form.glued_faces_valid();
        if !self.form.my_glued_f.is_empty() {
            self.form.my_just_gluer = true;
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        }

        // `cxx:254-275`: without gluing, a plain Boolean on the primitive.
        if !self.form.my_just_gluer {
            self.finish_with_boolean_on(&vraiprism)?;
        }
        Ok(())
    }

    /// `PerformUntilHeight(Until, Height)` (`cxx:880-997`).
    ///
    /// The reference calls `Perform(Height)` on a null `Until` and
    /// `Perform(Until)` on a zero height *without returning*, so the selected
    /// overload runs and the body continues; that control flow is kept.
    pub fn perform_until_height(&mut self, until: &TopoShape, height: f64) -> Result<(), String> {
        // `cxx:888-895`.
        if until.is_null() {
            self.perform_height(height)?;
        }
        if height == 0.0 {
            self.perform_until(until)?;
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:901-908`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:909-913`.
        let pbase = self.pbase_face()?;
        let c = statics::test_curve(&pbase);
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            // `TestCurve` answers a null handle on a non-planar profile; the
            // direction is undefined there, and the reference would raise.
            None => 1,
        };
        let the_dprism = DPrism::with_height(&pbase, sens as f64 * height, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:916-946` / `cxx:937-945`.
        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trf {
            // `cxx:916-935`.
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            if !statics::first_shape_has_one_face(&the_dprism)? {
                self.set_error(StatusError::InvFirstShape);
                return Ok(());
            }
            self.form.glued_faces_valid();
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        } else {
            // `cxx:937-996`.
            let scur: Vec<Option<Arc<dyn Curve>>> = match &c {
                Some(c) => vec![Some(c.clone())],
                None => vec![None],
            };
            let mut asi = CSIntersector::with_shape(until);
            let _ = asi.perform_curves(&scur);
            let nb = if asi.is_done() {
                asi.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb >= 1 {
                // `cxx:955-964`.
                let or = if self.form.my_fuse {
                    asi.point(1, 1).ok().map(|p| p.orientation())
                } else {
                    asi.point(1, nb).ok().map(|p| p.orientation())
                }
                .unwrap_or(Orientation::Forward);
                // `cxx:966-977`: the tool solid inside a one-element compound.
                let mut tools: Vec<TopoShape> = Vec::new();
                if let Some(f_until) = asi.point(1, 1).ok().map(|p| p.face().clone()) {
                    if let Some(s) = tool(until, &f_until, or) {
                        tools.push(s.0);
                    }
                }
                let comp = TopoBuilder::new().make_compound_of(&tools);
                // `cxx:978-979`.
                let (tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
                update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;
                // `cxx:979-996`.
                self.finish_with_boolean_on(&tr_shape)?;
            }
        }
        Ok(())
    }
}
