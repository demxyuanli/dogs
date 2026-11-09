//! `BRepFeat_MakeCylindricalHole::Perform*` without the blind variant.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.cxx`:
//!
//! | Reference | Rust |
//! |---|---|
//! | `Perform(Radius)` (`cxx:56-100`) | [`MakeCylindricalHole::perform`] |
//! | `PerformThruNext` (`cxx:106-253`) | [`MakeCylindricalHole::perform_thru_next`] |
//! | `PerformUntilEnd` (`cxx:257-353`) | [`MakeCylindricalHole::perform_until_end`] |
//! | `Perform(Radius, PFrom, PTo, Cont)` (`cxx:357-470`) | [`MakeCylindricalHole::perform_range`] |
//!
//! Every entry point follows the same five steps: locate the hits along the
//! axis, build the cylinder between the two selected hits, run the tool through
//! the inherited `BOPAlgo_BOP::Perform`, count the solids of the result, and
//! keep the parts the reference keeps when there is more than one.

use occt_core::elib::clib::line_parameter_ax1;
use occt_core::gp::{GpAx2, GpPnt};

use crate::abs::Orientation;
use crate::brep_feat::status::Status;
use crate::builder::TopoBuilder;
use crate::locope::curve_shape_intersector::CurveShapeIntersector;
use crate::locope::pnt_face::PntFace;
use crate::shape::TopoShape;

use super::statics::{baryc, box_parameters, create_cyl, cylinder};
use super::MakeCylindricalHole;

impl MakeCylindricalHole {
    /// `Perform(Radius)` (`cxx:56-100`): cut with an infinite cylinder, done
    /// through a finite one that spans twice the shape's extent.
    pub fn perform(&mut self, radius: f64) -> Result<(), String> {
        let object = self.checked_object()?;

        self.my_is_blind = false;
        self.my_status = Status::NoError;

        let asi = CurveShapeIntersector::with_axis(&self.my_axis, &object);
        if !asi.is_done() || asi.nb_points() == 0 {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        // `cxx:76-81`: "It is not possible to use infinite cylinder for
        // topological operations."
        let Some((pmin, pmax)) = box_parameters(&object, &self.my_axis) else {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        };
        let height = 2.0 * (pmax - pmin);
        let orig = self
            .my_axis
            .location()
            .xyz()
            .add(&self.my_axis.direction().xyz().multiply_scalar(
                (3.0 * pmin - pmax) / 2.0,
            ));
        let pnt = GpPnt::new(orig.x(), orig.y(), orig.z());
        let a2 = GpAx2::from_axis(pnt, *self.my_axis.direction());
        let (solid, top, bottom) = cylinder(&a2, radius, height)?;

        self.my_top_face = Some(top);
        self.my_bot_face = Some(bottom);
        self.my_validate = false;
        self.add_tool_and_perform(&solid.0)
    }

    /// `PerformThruNext(Radius, Cont)` (`cxx:106-253`): the first hole the axis
    /// meets, in the direction of the axis.
    pub fn perform_thru_next(&mut self, radius: f64, with_control: bool) -> Result<(), String> {
        let object = self.checked_object()?;

        self.my_is_blind = false;
        self.my_validate = with_control;
        self.my_status = Status::NoError;

        let asi = CurveShapeIntersector::with_axis(&self.my_axis, &object);
        if !asi.is_done() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let mut first = PntFace::default();
        let mut last = PntFace::default();
        let mut ok = false;
        if let Some((the_or, ind_from, ind_to)) = asi.localize_after_param(0.0) {
            ok = true;
            if the_or == Orientation::Forward {
                first = point_of(&asi, ind_from)?;
                match asi.localize_after_index(ind_to) {
                    Some((o, _f, t)) => {
                        if o != Orientation::Reversed {
                            ok = false;
                        } else {
                            last = point_of(&asi, t)?;
                        }
                    }
                    None => ok = false,
                }
            } else {
                // `TopAbs_REVERSED`.
                last = point_of(&asi, ind_to)?;
                match asi.localize_before_index(ind_from) {
                    Some((o, f, _t)) => {
                        if o != Orientation::Forward {
                            ok = false;
                        } else {
                            first = point_of(&asi, f)?;
                        }
                    }
                    None => ok = false,
                }
            }
        }
        if !ok {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let (shell, top, bottom) = create_cyl(&first, &last, radius, &self.my_axis)?;
        let builder = TopoBuilder::new();
        let solid = builder.make_solid(&[shell]);
        self.my_top_face = Some(top);
        self.my_bot_face = Some(bottom);
        self.add_tool_and_perform(&solid.0)?;

        let parts = self.base.parts_of_tool();
        if parts.is_empty() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }
        if parts.len() >= 2 {
            // `cxx:198-252`: keep the smallest parameter along the axis.
            let first_par = first.parameter();
            let last_par = last.parameter();
            let mut tokeep: Option<TopoShape> = None;
            let mut parmin = last_par;
            for part in &parts {
                let bar = baryc(part);
                let parbar = line_parameter_ax1(&self.my_axis, &bar);
                if parbar >= first_par && parbar <= last_par && parbar <= parmin {
                    parmin = parbar;
                    tokeep = Some(part.clone());
                }
            }
            if tokeep.is_none() {
                // `cxx:217-241`: "preserve the closest interval". The
                // reference nests the `parbar > Last` arm in the `else` of
                // `First - parbar < dmin`, so that arm only runs when the part
                // lies before `First`; the port keeps that nesting as written.
                let mut dmin = f64::MAX;
                for part in &parts {
                    let bar = baryc(part);
                    let parbar = line_parameter_ax1(&self.my_axis, &bar);
                    if parbar < first_par {
                        if first_par - parbar < dmin {
                            dmin = first_par - parbar;
                            tokeep = Some(part.clone());
                        } else {
                            // parbar > Last
                            if parbar - last_par < dmin {
                                dmin = parbar - last_par;
                                tokeep = Some(part.clone());
                            }
                        }
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

    /// `PerformUntilEnd(Radius, Cont)` (`cxx:257-353`): every hole after the
    /// origin of the axis.
    pub fn perform_until_end(&mut self, radius: f64, with_control: bool) -> Result<(), String> {
        let object = self.checked_object()?;

        self.my_is_blind = false;
        self.my_validate = with_control;
        self.my_status = Status::NoError;

        let asi = CurveShapeIntersector::with_axis(&self.my_axis, &object);
        if !asi.is_done() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let mut first = PntFace::default();
        let mut last = PntFace::default();
        let mut ok = false;
        if let Some((mut the_or, mut ind_from, _ind_to)) = asi.localize_after_param(0.0) {
            ok = true;
            if the_or == Orientation::Reversed {
                // `cxx:283-286` "on reset; it is possible to search for the next".
                match asi.localize_before_index(ind_from) {
                    Some((o, f, _t)) => {
                        the_or = o;
                        ind_from = f;
                    }
                    None => ok = false,
                }
            }
            if ok && the_or == Orientation::Forward {
                first = point_of(&asi, ind_from)?;
                match asi.localize_before_index(asi.nb_points() + 1) {
                    Some((o, _f, t)) => {
                        if o != Orientation::Reversed {
                            ok = false;
                        } else {
                            last = point_of(&asi, t)?;
                        }
                    }
                    None => ok = false,
                }
            }
        }
        if !ok {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let (shell, top, bottom) = create_cyl(&first, &last, radius, &self.my_axis)?;
        let builder = TopoBuilder::new();
        let solid = builder.make_solid(&[shell]);
        self.my_top_face = Some(top);
        self.my_bot_face = Some(bottom);
        self.add_tool_and_perform(&solid.0)?;

        let parts = self.base.parts_of_tool();
        if parts.is_empty() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }
        if parts.len() >= 2 {
            // `cxx:338-351`: keep everything above `First`.
            let first_par = first.parameter();
            for part in &parts {
                let bar = baryc(part);
                let parbar = line_parameter_ax1(&self.my_axis, &bar);
                if parbar > first_par {
                    self.base.keep_part(part);
                }
            }
        }
        Ok(())
    }

    /// `Perform(Radius, PFrom, PTo, Cont)` (`cxx:357-470`): every hole between
    /// the two parameters.
    pub fn perform_range(
        &mut self,
        radius: f64,
        p_from: f64,
        p_to: f64,
        with_control: bool,
    ) -> Result<(), String> {
        let object = self.checked_object()?;

        self.my_is_blind = false;
        self.my_validate = with_control;
        self.my_status = Status::NoError;

        let asi = CurveShapeIntersector::with_axis(&self.my_axis, &object);
        if !asi.is_done() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let (the_p_from, the_p_to) = if p_from < p_to {
            (p_from, p_to)
        } else {
            (p_to, p_from)
        };

        let mut first = PntFace::default();
        let mut last = PntFace::default();
        let mut ok = false;
        if let Some((mut the_or, mut ind_from, _ind_to)) = asi.localize_after_param(the_p_from) {
            ok = true;
            if the_or == Orientation::Reversed {
                // `cxx:400-403`: reset, it is possible to find the next.
                match asi.localize_before_index(ind_from) {
                    Some((o, f, _t)) => {
                        the_or = o;
                        ind_from = f;
                    }
                    None => ok = false,
                }
            }
            if ok && the_or == Orientation::Forward {
                first = point_of(&asi, ind_from)?;
                match asi.localize_before_param(the_p_to) {
                    Some((o, _f, mut t)) => {
                        the_or = o;
                        if the_or == Orientation::Forward {
                            match asi.localize_after_index(t) {
                                Some((o2, _f2, t2)) => {
                                    the_or = o2;
                                    t = t2;
                                }
                                None => ok = false,
                            }
                        }
                        if ok && the_or == Orientation::Reversed {
                            last = point_of(&asi, t)?;
                        }
                    }
                    None => ok = false,
                }
            }
        }
        if !ok {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }

        let (shell, top, bottom) = create_cyl(&first, &last, radius, &self.my_axis)?;
        let builder = TopoBuilder::new();
        let solid = builder.make_solid(&[shell]);
        self.my_top_face = Some(top);
        self.my_bot_face = Some(bottom);
        self.add_tool_and_perform(&solid.0)?;

        let parts = self.base.parts_of_tool();
        if parts.is_empty() {
            self.my_status = Status::InvalidPlacement;
            return Ok(());
        }
        if parts.len() >= 2 {
            // `cxx:457-469`: keep the parts between `First` and `Last`.
            let first_par = first.parameter();
            let last_par = last.parameter();
            for part in &parts {
                let bar = baryc(part);
                let parbar = line_parameter_ax1(&self.my_axis, &bar);
                if parbar >= first_par && parbar <= last_par {
                    self.base.keep_part(part);
                }
            }
        }
        Ok(())
    }
}

/// `theASI.Point(Ind)`; a missing index is the reference's out-of-range
/// `Point` (`Standard_OutOfRange`).
fn point_of(asi: &CurveShapeIntersector, index: usize) -> Result<PntFace, String> {
    asi.point(index)
        .cloned()
        .ok_or_else(|| "Standard_OutOfRange: LocOpe_CurveShapeIntersector::Point".to_string())
}
