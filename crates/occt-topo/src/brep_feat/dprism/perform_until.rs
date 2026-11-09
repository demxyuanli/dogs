//! `BRepFeat_MakeDPrism::Perform(Until)` (`cxx:283-414`).
//!
//! The draft prism is swept by the bounding-box extent along the direction
//! `SensOfPrism` finds, so it always crosses `Until`. When `TransformShapeFU`
//! did not trim an infinite support (`Trf == false`) the limit is a proper face
//! and the regular `GlobalPerform` runs; otherwise the primitive is cut by the
//! tool solid built from `Until` and the kept solid is fused / cut with the
//! basis. The reference then re-parameterises every new edge that fails
//! `BRepAlgo::IsValid` (`cxx:403-414`).

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_check;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::form::update_descendants_bop;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::tool;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::{CSIntersector, DPrism};
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, faces_of};

use super::statics;
use super::MakeDPrism;

impl MakeDPrism {
    /// `Perform(Until)` (`cxx:283-414`).
    pub fn perform_until(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:290-299`.
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }

        // `cxx:301-308`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:309-313`.
        let pbase = self.pbase_face()?;
        let c = statics::test_curve(&pbase);
        let sens = match &c {
            Some(c) => statics::sens_of_prism(c, until),
            // `TestCurve` answers a null handle on a non-planar profile; the
            // direction is undefined there, and the reference would raise.
            None => 1,
        };
        let height = statics::height_max(
            self.form.my_sbase.as_ref(),
            self.form.my_skface.as_ref(),
            self.form.my_sfrom.as_ref(),
            self.form.my_suntil.as_ref(),
        );
        let the_dprism = DPrism::with_height(&pbase, sens as f64 * height, self.angle());
        let vraiprism = the_dprism
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;

        // `cxx:318` and `cxx:338` both start with the same `MajMap`.
        statics::maj_map(
            &pbase.0,
            &the_dprism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        if !trf {
            // `cxx:319-336`: the limit is a face of the basis.
            self.form.my_gshape = Some(vraiprism);
            self.form.generated_shape_valid();
            if !statics::first_shape_has_one_face(&the_dprism)? {
                self.set_error(crate::brep_feat::status_error::StatusError::InvFirstShape);
                return Ok(());
            }
            self.form.glued_faces_valid();
            self.take_dprism_curves(&the_dprism);
            self.global_perform()?;
        } else {
            // `cxx:337-401`: the infinite support was trimmed, so the trim is a
            // Boolean against the tool solid built from `Until`.
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
                // `cxx:355-366`.
                let or = if self.form.my_fuse {
                    asi.point(1, 1).ok().map(|p| p.orientation())
                } else {
                    asi.point(1, nb).ok().map(|p| p.orientation())
                }
                .unwrap_or(Orientation::Forward);
                // `cxx:368-376`: the tool solid inside a one-element compound.
                let mut tools: Vec<TopoShape> = Vec::new();
                if let Some(f_until) = asi.point(1, 1).ok().map(|p| p.face().clone()) {
                    if let Some(s) = tool(until, &f_until, or) {
                        tools.push(s.0);
                    }
                }
                let comp = TopoBuilder::new().make_compound_of(&tools);
                // `cxx:378-379`.
                let (tr_p, tr_shape) = statics::boolean_result(cut(&vraiprism, &comp.0))?;
                update_descendants_bop(&mut self.form, &tr_p, &tr_shape, false)?;
                // `cxx:381-382`: `ex.Current()` raises on an empty explorer, so
                // the missing first solid is `Standard_OutOfRange`.
                let cutsh = explorer_shapes(&tr_shape, ShapeType::Solid)
                    .into_iter()
                    .next()
                    .ok_or_else(|| "Standard_OutOfRange".to_string())?;
                // `cxx:383-400`.
                self.finish_with_boolean_on(&cutsh)?;
            }
        }

        // `cxx:403-414`.
        self.fix_new_edges()?;
        Ok(())
    }

    /// `cxx:403-414`: every edge of `myNewEdges` that fails `BRepAlgo::IsValid`
    /// loses its same-range / same-parameter flags and is re-parameterised.
    fn fix_new_edges(&mut self) -> Result<(), String> {
        for e in self.form.my_new_edges.clone() {
            // `TopoDS::Edge(ited.Value())` raises on a non-edge entry.
            let edge = Edge::wrap(e).ok_or_else(|| "Standard_DomainError".to_string())?;
            if !brep_check::is_valid(&edge.0) {
                let reg = GeometryRegistry::global();
                reg.set_same_range(&edge.0, false);
                reg.set_same_parameter(&edge.0, false);
                same_parameter_edge_inplace(&edge, BRepTool::edge_tolerance(&edge));
            }
        }
        Ok(())
    }
}
