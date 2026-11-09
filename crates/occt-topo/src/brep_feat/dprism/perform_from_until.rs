//! `BRepFeat_MakeDPrism::Perform(From, Until)` (`cxx:418-617`).
//!
//! The draft prism is swept by the bounding-box extent in both directions
//! (`Height1 == Height2 == Height`, so it covers `From` and `Until`). When both
//! supports are proper faces the regular `GlobalPerform` runs; otherwise both
//! limits are turned into tool solids and cut from the primitive.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::brep_feat::statics::tool;
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, DPrism};
use crate::shape::TopoShape;
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakeDPrism;

impl MakeDPrism {
    /// `Perform(From, Until)` (`cxx:418-617`).
    pub fn perform_from_until(&mut self, from: &TopoShape, until: &TopoShape) -> Result<(), String> {
        // `cxx:426-429`.
        if from.is_null() || until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:431-450`: a limit that is the sketch face means "glue from here",
        // which is the single-limit overload; the call returns only when the
        // gluer answered (`cxx:439`, `cxx:448`).
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

        // `cxx:452-457`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionFU;
        self.form.perf_selection_valid();

        // `cxx:459-466`.
        if faces_of(from).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:468-478`.
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

        // `cxx:479-492`. `HeightMax` is handed `myPbase` in place of the sketch
        // face, as the reference does (`cxx:490`).
        let pbase = self.pbase_face()?;
        let c = statics::test_curve(&pbase);
        let sens = if is_same(from, until) {
            1
        } else {
            match &c {
                Some(c) => statics::sens_of_prism(c, until),
                // `TestCurve` answers a null handle on a non-planar profile;
                // the direction is undefined there, and the reference raises.
                None => 1,
            }
        };
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.my_pbase.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let the_dprism = DPrism::with_two_heights(
            &pbase,
            sens as f64 * height,
            sens as f64 * height,
            self.angle(),
        );
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:495` and `cxx:513` both start with the same `MajMap`.
        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trff {
            // `cxx:496-511`: both limits are faces of the basis.
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        } else {
            // `cxx:512-616`: both supports were trimmed, so each limit becomes a
            // tool solid.
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

            // `cxx:533-553`.
            let nb1 = if asi1.is_done() {
                asi1.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb1 < 1 {
                self.set_error(StatusError::NoIntersectU);
                return Ok(());
            }
            let or_u = if self.form.my_fuse {
                asi1.point(1, 1).ok().map(|p| p.orientation())
            } else {
                asi1.point(1, nb1).ok().map(|p| p.orientation())
            }
            .unwrap_or(Orientation::Forward);
            let f_until = asi1.point(1, 1).ok().map(|p| p.face().clone());

            // `cxx:555-566`.
            let nb2 = if asi2.is_done() {
                asi2.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb2 < 1 {
                self.set_error(StatusError::NoIntersectF);
                return Ok(());
            }
            let or_f = asi2
                .point(1, 1)
                .map(|p| p.orientation().reversed())
                .unwrap_or(Orientation::Forward);
            let f_from = asi2.point(1, 1).ok().map(|p| p.face().clone());

            // `cxx:568-591`.
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

            // `cxx:593-594`.
            let comp = TopoBuilder::new().make_compound_of(&tools);
            // `cxx:596-615`: no `UpdateDescendants` on the intermediate cut in
            // this overload.
            let (_tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
            self.finish_with_boolean_on(&tr_shape)?;
        }
        Ok(())
    }
}
