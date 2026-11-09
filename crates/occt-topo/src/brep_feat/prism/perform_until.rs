//! `BRepFeat_MakePrism::Perform(Until)` (`cxx:325-433`).
//!
//! The prism is swept twice the bounding-box height along the direction found by
//! `SensOfPrism`, so it always crosses `Until`. When `TransformShapeFU` did not
//! trim an infinite support (`Trf == false`) the limit is a proper face and the
//! regular `GlobalPerform` runs; otherwise the primitive is cut by the tool
//! solid built from `Until` and the kept solid is fused / cut with the basis.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::tool;
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, Prism};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, faces_of};

use super::statics;
use super::MakePrism;

impl MakePrism {
    /// `Perform(Until)` (`cxx:325-433`).
    pub fn perform_until(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:337-343`.
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:344-350`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:351-358`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let c = statics::test_curve(&pbase, self.dir());
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            None => 1,
        };
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let v = self.scaled_dir(2.0 * sens as f64 * height);
        let the_prism = Prism::with_base(&pbase, &v);
        let vraiprism = the_prism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        if !trf {
            // `cxx:360-375`: the limit is a face of the basis.
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
            // `cxx:376-431`: the infinite support was trimmed, so the trim is a
            // Boolean against the tool solid built from `Until`.
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
            let mut asi = CSIntersector::with_shape(until);
            let _ = asi.perform_curves(&scur);
            let nb = if asi.is_done() {
                asi.nb_points(1).unwrap_or(0)
            } else {
                0
            };
            if nb >= 1 {
                // `cxx:384-398`.
                let mut or = if self.form.my_fuse {
                    asi.point(1, 1).ok().map(|p| p.orientation())
                } else {
                    asi.point(1, nb).ok().map(|p| p.orientation())
                }
                .unwrap_or(Orientation::Forward);
                if sens == -1 {
                    or = or.reversed();
                }
                // `cxx:398-403`: the tool solid inside a one-element compound.
                let mut tools: Vec<TopoShape> = Vec::new();
                if let Some(f_until) = asi.point(1, 1).ok().map(|p| p.face().clone()) {
                    if let Some(s) = tool(until, &f_until, or) {
                        tools.push(s.0);
                    }
                }
                let comp = TopoBuilder::new().make_compound_of(&tools);
                // `cxx:404-408`.
                let (tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
                update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;
                // `cxx:409-410`: `ex.Current()` raises on an empty explorer, so
                // the missing first solid is `Standard_OutOfRange`.
                let cutsh = explorer_shapes(&tr_shape, ShapeType::Solid)
                    .into_iter()
                    .next()
                    .ok_or_else(|| "Standard_OutOfRange".to_string())?;
                // `cxx:411-430`.
                self.finish_with_boolean_on(&cutsh)?;
            }
        }
        Ok(())
    }
}
