//! `BRepFeat_MakeRevol::Perform(Until)` (`cxx:342-467`) and
//! `PerformThruAll` (`cxx:655-667`).

use std::f64::consts::PI;

use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::locope::Revol;
use crate::shape::TopoShape;
use crate::topo_tools::{faces_of, is_same};

use super::statics;
use super::MakeRevol;

impl MakeRevol {
    /// `Perform(Until)` (`cxx:342-467`).
    ///
    /// The sweep is a nearly complete turn (`2*pi - 3 deg`) unless `Until` is
    /// the sketch face, where the full turn is glued. When `TransformShapeFU`
    /// trimmed an infinite support the trimmed revolution is cut by the tool
    /// solid of `Until` and the solid that keeps a face of `myPbase` is used.
    pub fn perform_until(&mut self, until: &TopoShape) -> Result<(), String> {
        // `cxx:348-357`.
        let mut angle = 0.0;
        let mut tour_complet = false;
        if until.is_null() {
            return Err("Standard_ConstructionError".to_string());
        }
        if faces_of(until).is_empty() {
            return Err("Standard_ConstructionError".to_string());
        }
        // `cxx:359-362`.
        if let Some(sk) = self.form.my_skface.clone() {
            if is_same(until, &sk.0) {
                angle = 2.0 * PI;
                tour_complet = true;
            }
        }

        // `cxx:363-371`.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = PerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = None;
        self.form.shape_from_valid();
        self.form.my_suntil = Some(until.clone());
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // `cxx:373-386`.
        let pbase = self
            .my_pbase
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let mut the_revol = Revol::new();
        if !tour_complet {
            angle = 2.0 * PI - 3.0 * PI / 180.0;
        }
        the_revol.perform(&pbase, self.axis(), angle);
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

        if !trf {
            // `cxx:388-410`.
            self.form.my_gshape = Some(vrairevol);
            self.form.generated_shape_valid();
            let first = the_revol
                .first_shape()
                .cloned()
                .ok_or_else(|| "StdFail_NotDone".to_string())?;
            let base_faces = faces_of(&first);
            if base_faces.is_empty() {
                return Err("Standard_OutOfRange".to_string());
            }
            if base_faces.len() > 1 {
                self.set_error(StatusError::InvFirstShape);
                return Ok(());
            }
            self.form.glued_faces_valid();
            self.take_revol_curves(&the_revol);
            self.global_perform()?;
        } else {
            // `cxx:411-466`.
            self.take_revol_curves(&the_revol);
            if let Some(v) = self.cut_pick_solid(&vrairevol, until)? {
                self.finish_fuse_or_cut(&v)?;
            }
        }
        Ok(())
    }

    /// `PerformThruAll()` (`cxx:655-667`): a full revolution.
    pub fn perform_thru_all(&mut self) -> Result<(), String> {
        self.perform_angle(2.0 * PI)
    }
}
