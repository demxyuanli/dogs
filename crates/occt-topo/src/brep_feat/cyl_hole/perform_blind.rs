//! `BRepFeat_MakeCylindricalHole::PerformBlind`.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.cxx:474-626`.
//!
//! A blind hole is a cylinder that starts at the first hit of the axis and
//! reaches `Length` into the material. The reference checks the length a priori
//! (`cxx:519-528`) and again through the bounding box (`cxx:541-547`), then
//! keeps the tool part closest to `First`.

use occt_core::elib::clib::line_parameter_ax1;
use occt_core::gp::{GpAx2, GpPnt};

use crate::abs::Orientation;
use crate::brep_feat::status::Status;
use crate::locope::curve_shape_intersector::CurveShapeIntersector;
use crate::locope::pnt_face::PntFace;
use crate::shape::TopoShape;

use super::statics::{baryc, box_parameters, cylinder};
use super::MakeCylindricalHole;

impl MakeCylindricalHole {
    /// `PerformBlind(Radius, Length, Cont)` (`cxx:474-626`).
    pub fn perform_blind(
        &mut self,
        radius: f64,
        length: f64,
        with_control: bool,
    ) -> Result<(), String> {
        let object = self.checked_object()?;
        if length <= 0.0 {
            return Err("Standard_ConstructionError".to_string());
        }

        self.my_is_blind = true;
        self.my_validate = with_control;
        self.my_status = Status::NoError;

        let asi = CurveShapeIntersector::with_axis(&self.my_axis, &object);
        if !asi.is_done() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let mut first = PntFace::default();
        let mut ind_from = 0usize;
        let mut ind_to = 0usize;
        let mut ok = false;
        if let Some((mut the_or, f, t)) = asi.localize_after_param(0.0) {
            ind_from = f;
            ind_to = t;
            ok = true;
            if the_or == Orientation::Reversed {
                // `cxx:505-508`: reset, it is possible to find the next.
                match asi.localize_before_index(ind_from) {
                    Some((o, f, t)) => {
                        the_or = o;
                        ind_from = f;
                        ind_to = t;
                    }
                    None => ok = false,
                }
            }
            ok = ok && the_or == Orientation::Forward;
        }
        if !ok {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        // `cxx:516-528`: check the length of the hole a priori.
        let Some((_or_next, if_next, _it_next)) = asi.localize_after_index(ind_to) else {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        };
        let Some(pnt_next) = asi.point(if_next) else {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        };
        if pnt_next.parameter() <= length {
            self.my_status = Status::HoleTooLong;
            return Ok(());
        }

        // `cxx:533-537` also fills `theList` with the faces between the hits;
        // the reference never reads it, so the port omits the dead fill.

        if let Some(pnt_first) = asi.point(ind_from) {
            first = pnt_first.clone();
        }
        let first_par = first.parameter();

        // `cxx:539-556`: the finite cylinder that stands in for the blind tool.
        let Some((pmin, _pmax)) = box_parameters(&object, &self.my_axis) else {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        };
        if pmin > length {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }
        let height = 3.0 * (length - pmin) / 2.0;
        let orig = self
            .my_axis
            .location()
            .xyz()
            .add(&self.my_axis.direction().xyz().multiply_scalar(
                (3.0 * pmin - length) / 2.0,
            ));
        let pnt = GpPnt::new(orig.x(), orig.y(), orig.z());
        let a2 = GpAx2::from_axis(pnt, *self.my_axis.direction());
        let (solid, top, _bottom) = cylinder(&a2, radius, height)?;

        self.my_top_face = Some(top);
        self.my_bot_face = None;

        // `cxx:567-570`.
        self.add_tool_and_perform(&solid.0)?;

        let parts = self.base.parts_of_tool();
        if parts.is_empty() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }
        if parts.len() >= 2 {
            // `cxx:587-624`: keep the smallest parameter above `First`.
            let mut tokeep: Option<TopoShape> = None;
            let mut parmin = f64::MAX;
            for part in &parts {
                let bar = baryc(part);
                let parbar = line_parameter_ax1(&self.my_axis, &bar);
                if parbar >= first_par && parbar <= parmin {
                    parmin = parbar;
                    tokeep = Some(part.clone());
                }
            }
            if tokeep.is_none() {
                // `cxx:603-615`: preserve the closest interval.
                let mut dmin = f64::MAX;
                for part in &parts {
                    let bar = baryc(part);
                    let parbar = line_parameter_ax1(&self.my_axis, &bar);
                    if (first_par - parbar).abs() < dmin {
                        dmin = (first_par - parbar).abs();
                        tokeep = Some(part.clone());
                    }
                }
            }
            if let Some(tokeep) = tokeep {
                for part in &parts {
                    if tokeep.same_tshape(part) {
                        self.base.keep_part(part);
                        break;
                    }
                }
            }
        }
        Ok(())
    }
}
