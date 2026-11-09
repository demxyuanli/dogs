//! `BRepFeat_MakeLinearForm::Perform` (`BRepFeat_MakeLinearForm.cxx:905-1028`).
//!
//! The prism of `myPbase` along `Length * myDir`, then the descendant maps
//! (`MajMap`), the glued faces (`SetGluedFaces`) and the topological
//! reconstruction of the base class (`LFPerform`).

use occt_core::gp::GpVec;
use occt_core::precision;

use crate::abs::ShapeType;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::locope::LinearForm;
use crate::topo_tools::explorer_shapes;

use super::{maj_map, set_glued_faces, MakeLinearForm};

impl MakeLinearForm {
    /// `Perform()` (`cxx:905-1028`).
    ///
    /// `Err` is the reference's raised `Standard_NullObject` when `myDir` is
    /// null (`gp_Dir(gp_Vec)` at `cxx:937`).
    pub fn perform(&mut self) -> Result<(), String> {
        // `cxx:912-921`.
        if self.base.my_sbase.is_none()
            || self.base.my_skface.is_none()
            || self.base.my_pbase.is_none()
        {
            self.base.my_status_error = StatusError::NotInitialized;
            self.base.not_done();
            return Ok(());
        }

        // `cxx:924-937`.
        let nulldir = GpVec::zero();
        let length = self.my_dir.magnitude() + self.my_dir1.magnitude();

        self.base.my_glued_f.clear();

        if self.base.my_suntil.is_some() {
            self.base.my_perf_selection = PerfSelection::SelectionU;
        } else {
            self.base.my_perf_selection = PerfSelection::NoSelection;
        }

        let dir = occt_core::gp::GpDir::from_vec(&self.my_dir)
            .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        let v = GpVec::from_xyz(dir.xyz()).multiplied_scalar(length);

        // `cxx:941-951`.
        let pbase = self
            .base
            .my_pbase
            .clone()
            .ok_or_else(|| "Standard_NullObject: myPbase".to_string())?;
        let mut the_form = LinearForm::new();
        if self
            .my_dir1
            .is_equal_tol(&nulldir, precision::CONFUSION, precision::CONFUSION)
        {
            the_form.perform(
                &pbase.0,
                &v,
                &self.base.my_first_pnt,
                &self.base.my_last_pnt,
            );
        } else {
            the_form.perform_with_translation(
                &pbase.0,
                &v,
                &self.my_dir1,
                &self.base.my_first_pnt,
                &self.base.my_last_pnt,
            );
        }

        // `cxx:953-956`.
        let vrai_form = the_form
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: LocOpe_LinearForm".to_string())?;
        // The reference appends `FirstShape()` / `LastShape()` as they are,
        // including a null shape; the port only keeps a real shape because its
        // list holds `TopoShape` values.
        if let Some(f) = the_form.first_shape() {
            self.base.my_faces_for_draft.push(f.clone());
        }
        if let Some(l) = the_form.last_shape() {
            self.base.my_faces_for_draft.push(l.clone());
        }
        maj_map(
            &pbase.0,
            &the_form,
            &mut self.base.my_map,
            &mut self.base.my_f_shape,
            &mut self.base.my_l_shape,
        );

        // `cxx:957-972`: every edge of the profile has to have an image.
        for e in explorer_shapes(&pbase.0, ShapeType::Edge) {
            if !self.base.my_map.is_bound(&e) {
                self.base.my_status_error = StatusError::IncSlidFace;
                self.base.not_done();
                return Ok(());
            }
        }

        // `cxx:974-975`.
        self.base.my_gshape = Some(vrai_form);
        set_glued_faces(&self.my_slface, &the_form, &mut self.base.my_glued_f);
        let _ = nulldir;

        // `cxx:977-989`.
        if !self.base.my_glued_f.is_empty() && self.base.my_suntil.is_some() {
            self.base.my_status_error = StatusError::InvShape;
            self.base.not_done();
            return Ok(());
        }

        // `cxx:991`.
        self.base.lf_perform()?;
        Ok(())
    }
}
