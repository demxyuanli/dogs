//! `BRepFill_PipeShell::Prepare` / `Place` / `ResetLoc` and the statics they
//! need (`BRepFill_PipeShell.cxx:130-148`, `:966-1258`).

use std::sync::Arc;

use occt_core::gp::GpTrsf;
use occt_geom::geom_fill::{
    CorrectedFrenet, CurveAndTrihedron, GeomFillPipeError, GeomFillTrihedron,
    LocationLaw as GeomFillLocationLaw, TrihedronLaw,
};

use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_edge3d_law::Edge3DLaw;
use crate::brep_fill_location_law::LocationLaw;
use crate::brep_fill_n_sections::NSections;
use crate::brep_fill_pipe_shell::PipeShell;
use crate::brep_fill_section::Section;
use crate::brep_fill_section_placement::SectionPlacement;
use crate::brep_fill_shape_law::ShapeLaw;
use crate::brep_fill_sweep::TransitionStyle;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, TopoShape, Wire};
use crate::topo_tools::shape_map_key;

impl PipeShell {
    /// `Prepare()` (`cxx:966-1234`).
    pub fn prepare(&mut self) -> Result<(), String> {
        self.w_seq_mut().clear();
        self.edge_new_edges_mut().clear();

        if !self.is_ready() {
            return Err("StdFail_NotDone: BRepFill_PipeShell::Prepare".to_string());
        }
        if self.location().is_some() && self.section().is_some() {
            return Ok(()); // It is ready
        }

        // `cxx:979-990`: a punctual section is only allowed at the two ends.
        let len = self.seq().len();
        if len >= 3 {
            for i in 1..=(len - 2) {
                if Self::section_is_degenerated(&self.seq()[i]) {
                    return Err("Standard_Failure: Wrong usage of punctual sections".to_string());
                }
            }
        }
        // `cxx:994-1010`: with at most two sections not all of them may be
        // punctual either.
        if len <= 2 {
            let mut wdeg = true;
            for sec in self.seq() {
                wdeg = wdeg && Self::section_is_degenerated(sec);
            }
            if wdeg {
                return Err("Standard_Failure: Wrong usage of punctual sections".to_string());
            }
        }

        // `cxx:1012-1028`: the location law, CorrectedFrenet unless `Set` ran.
        if self.location().is_none() {
            match self.trihedron() {
                GeomFillTrihedron::CorrectedFrenet => {
                    let tlaw: Arc<dyn TrihedronLaw> = Arc::new(CorrectedFrenet::new());
                    let loc: Arc<dyn GeomFillLocationLaw> = Arc::new(CurveAndTrihedron::new(tlaw));
                    let law = Edge3DLaw::new(self.spine_ref(), &loc);
                    *self.location_mut() = Some(Arc::new(law));
                }
                // `cxx:1022-1025`: "Not planned!".
                _ => {
                    return Err(
                        "Standard_ConstructionError: BRepFill_PipeShell::Prepare".to_string()
                    )
                }
            }
        }

        // `cxx:1030`: Transition Management.
        let mode = self.transition();
        let ang_min = self.ang_min();
        match self
            .location_mut()
            .as_mut()
            .and_then(|loc| Arc::get_mut(loc))
        {
            Some(loc) => perform_transition(mode, loc, ang_min)?,
            // The handle is shared only with a `BRepFill_Sweep`, which exists
            // only after `Build`, and `Build` leaves `mySection` set, so the
            // shared case cannot reach here. Report it rather than skipping
            // the transform.
            None => {
                return Err("BRepFill_PipeShell::Prepare: the location law is shared".to_string())
            }
        }

        // `cxx:1032-1062`: the section law, single-section case.
        if self.seq().len() == 1 {
            let (the_sect, _trsf, _p1) = self.place(&self.seq()[0])?;
            let a_local_shape = the_sect.clone();
            let shape_law = if self.seq()[0].is_law() {
                let law = self
                    .law_ref()
                    .cloned()
                    .ok_or_else(|| "BRepFill_PipeShell::Prepare: no homothetic law".to_string())?;
                ShapeLaw::with_wire_and_law(&a_local_shape, law, true)?
            } else {
                ShapeLaw::with_wire(&a_local_shape, true)?
            };
            *self.section_mut() = Some(Arc::new(shape_law));

            self.w_seq_mut().push(the_sect.clone());
            self.ind_of_sec_mut().push(1);

            // `cxx:1052-1061`: every edge of the working section maps to
            // itself, since nothing split it.
            for an_edge in ShapeIterator::of_shape(&the_sect.0) {
                self.edge_new_edges_mut()
                    .insert(shape_map_key(&an_edge), vec![an_edge]);
            }
        } else {
            // `cxx:1063-1205`: the multi-section branch.
            self.prepare_multi_section()?;
        }

        // `cxx:1207-1228`: the contact retiming, guide-only.
        if matches!(
            self.trihedron(),
            GeomFillTrihedron::GuidePlanWithContact | GeomFillTrihedron::GuideACWithContact
        ) {
            self.retime_contact()?;
        }

        // `cxx:1229-1233`.
        let status = self
            .location()
            .map(|loc| loc.base().get_status())
            .unwrap_or(GeomFillPipeError::PipeOk);
        self.set_status(status);
        let is_done = self
            .section()
            .map(|sec| sec.base().is_done())
            .unwrap_or(false);
        if !is_done {
            self.set_status(GeomFillPipeError::PipeNotOk);
        }
        Ok(())
    }

    /// The `mySeq.Length() > 1` arm of `Prepare` (`cxx:1063-1205`): place every
    /// section on the path, close the loop when the path is closed, order the
    /// sections by abscissa, make them compatible, and build the
    /// `BRepFill_NSections` law from the compatible wires.
    fn prepare_multi_section(&mut self) -> Result<(), String> {
        let loc = self
            .location()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Prepare: no location law".to_string())?;

        // `cxx:1068-1072`.
        let nb_l = loc.base().nb_law();
        let (_v1_bounds, v2) = loc.base().curvilinear_bounds(nb_l)?;
        let v1 = 0.0;

        let mut param: Vec<f64> = Vec::new();
        let mut ind_sec: Vec<usize> = Vec::new();
        let mut transformations: Vec<GpTrsf> = Vec::new();
        let mut w_seq: Vec<Wire> = Vec::new();
        let mut ideb = 0usize;
        let mut ifin = 0usize;

        // `cxx:1075-1089`.
        for iseq in 1..=self.seq().len() {
            ind_sec.push(iseq);
            let (the_sect, a_trsf, p) = self.place(&self.seq()[iseq - 1])?;
            param.push(p);
            w_seq.push(the_sect);
            transformations.push(a_trsf);
            if p == v1 {
                ideb = iseq;
            }
            if p == v2 {
                ifin = iseq;
            }
        }

        // `cxx:1092-1139`: looping sections.
        if loc.base().is_closed() {
            if ideb > 0 {
                // Place the initial section at the final position.
                param.push(v2);
                w_seq.push(w_seq[ideb - 1].clone());
            } else if ifin > 0 {
                // Place the final section at the initial position.
                param.push(v1);
                w_seq.push(w_seq[ifin - 1].clone());
            } else {
                // Find the extreme sections and impose them at `V1` and `V2`.
                let mut pmin = f64::MAX;
                let mut pmax = -f64::MAX;
                let mut wmin: Option<Wire> = None;
                let mut wmax: Option<Wire> = None;
                for iseq in 1..=w_seq.len() {
                    if param[iseq - 1] < pmin {
                        pmin = param[iseq - 1];
                        wmin = Some(w_seq[iseq - 1].clone());
                    }
                    if param[iseq - 1] > pmax {
                        pmax = param[iseq - 1];
                        wmax = Some(w_seq[iseq - 1].clone());
                    }
                }
                let dmin = (pmin - v1).abs();
                let dmax = (pmax - v2).abs();
                let (Some(wmin), Some(wmax)) = (wmin, wmax) else {
                    return Err(
                        "BRepFill_PipeShell::Prepare: the closed path has no section".to_string()
                    );
                };
                let wres = compute_section(&wmin, &wmax, dmin, dmax)?;
                // `cxx:1131-1139`: the same medium section at both bounds.
                param.push(v1);
                w_seq.push(wres.clone());
                ind_sec.push(w_seq.len());
                param.push(v2);
                w_seq.push(wres);
                ind_sec.push(w_seq.len());
            }
        }

        // `cxx:1142-1160`: parse the sections by increasing parameter. OCCT's
        // `NCollection_Sequence::Exchange` is a swap, and the outer loop reruns
        // until nothing moves.
        let mut play_again = true;
        while play_again {
            play_again = false;
            for iseq in 1..=w_seq.len() {
                for jseq in (iseq + 1)..=w_seq.len() {
                    if param[iseq - 1] > param[jseq - 1] {
                        param.swap(iseq - 1, jseq - 1);
                        w_seq.swap(iseq - 1, jseq - 1);
                        ind_sec.swap(iseq - 1, jseq - 1);
                        play_again = true;
                    }
                }
            }
        }

        // `cxx:1161-1171`: the real index of every input section.
        for ii in 1..=self.seq().len() {
            for jj in 1..=ind_sec.len() {
                if ind_sec[jj - 1] == ii {
                    self.ind_of_sec_mut().push(jj);
                    break;
                }
            }
        }

        // `cxx:1173-1202`: the compatible working sections.
        let sections: Vec<TopoShape> = w_seq.iter().map(|w| w.0.clone()).collect();
        let mut georges = CompatibleWires::with_sections(&sections);
        georges.set_percent(0.1);
        georges.perform(false)?;
        if !georges.is_done() {
            return Err("Standard_ConstructionError: PipeShell : uncompatible wires".to_string());
        }
        let working_sections = georges.shape().to_vec();

        // For each sub-edge of each section, save its splits.
        for wire in &sections {
            for shape in ShapeIterator::of_shape(wire) {
                let Some(edge) = Edge::wrap(shape) else {
                    continue;
                };
                let new_edges = georges.generated_shapes(&edge.0).to_vec();
                self.edge_new_edges_mut()
                    .insert(shape_map_key(&edge.0), new_edges);
            }
        }

        self.w_seq_mut().extend(w_seq);

        // `cxx:1203`: `mySection = new BRepFill_NSections(WorkingSections,
        // Transformations, Param, V1, V2)`.
        let law = NSections::with_transformations(
            working_sections,
            transformations,
            param,
            v1,
            v2,
            true,
        )?;
        *self.section_mut() = Some(Arc::new(law));
        Ok(())
    }

    /// `Place(Sec, W, Trsf, Param)` (`cxx:1242-1258`): the section moved into
    /// the local frame of the path, its transform and its abscissa.
    pub fn place(&self, sec: &Section) -> Result<(Wire, GpTrsf, f64), String> {
        let loc = self
            .location()
            .ok_or_else(|| "BRepFill_PipeShell::Place: no location law".to_string())?;
        let the_wire = sec
            .wire()
            .ok_or_else(|| "BRepFill_PipeShell::Place: the section has no wire".to_string())?;
        let the_vertex = sec.vertex().filter(|v| !v.0.is_null());

        let place = match the_vertex {
            Some(v) => SectionPlacement::with_vertex(
                loc.base(),
                &the_wire.0,
                &v.0,
                sec.with_contact(),
                sec.with_correction(),
            )?,
            None => SectionPlacement::new(
                loc.base(),
                &the_wire.0,
                sec.with_contact(),
                sec.with_correction(),
            )?,
        };

        let a_trsf = place.transformation().clone();
        // `W = TopoDS::Wire(BRepBuilderAPI_Transform(TmpWire, aTrsf, true))`.
        let moved = crate::transform::transformed(&the_wire.0, &a_trsf);
        let w = Wire::wrap(moved)
            .ok_or_else(|| "BRepFill_PipeShell::Place: the section is not a wire".to_string())?;
        let param = place.abscissa_on_path();
        Ok((w, a_trsf, param))
    }

    /// `ResetLoc()` (`cxx:1264-1275`): drop the rotation the guide laws keep.
    ///
    /// Only the two contact trihedra carry a `GeomFill_LocationGuide`; for
    /// every other trihedron OCCT runs the empty body.
    pub fn reset_loc(&mut self) {
        if !matches!(
            self.my_trihedron,
            GeomFillTrihedron::GuidePlanWithContact | GeomFillTrihedron::GuideACWithContact
        ) {
            return;
        }
        let Some(location) = self.my_location.as_ref() else {
            return;
        };
        for isec in 1..=location.base().nb_law() {
            if let Some(law) = location.base().law(isec) {
                if let Some(guide) = law.as_location_guide() {
                    guide.erase_rotation(); // remove the rotation
                }
            }
        }
    }

    /// The `myTrihedron` contact arm of `Prepare` (`cxx:1207-1228`): rescale
    /// the section law's domain onto every path edge and force the guide laws
    /// to rotate, chaining the previous edge's last angle into the next edge's
    /// first one.
    fn retime_contact(&mut self) -> Result<(), String> {
        let location = self
            .location()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Prepare: no location law".to_string())?;
        let section = self
            .section()
            .cloned()
            .ok_or_else(|| "BRepFill_PipeShell::Prepare: no section law".to_string())?;

        // `Sec = mySection->ConcatenedLaw()`.
        let sec = section
            .concatened_law()
            .ok_or_else(|| "StdFail_NotDone: ConcatenedLaw".to_string())?;

        // `myLocation->CurvilinearBounds(myLocation->NbLaw(), f, Length)`.
        let (_unused_f, length) = location.base().curvilinear_bounds(location.base().nb_law())?;
        // `Sec->GetDomain(fs, l)`.
        let (fs, l) = sec.get_domain();
        // `Delta = (l - fs) / Length`.
        let delta = (l - fs) / length;

        let mut old_angle = 0.0;
        for ipath in 1..=location.base().nb_law() {
            // `myLocation->CurvilinearBounds(ipath, f, l)`.
            let (f, l) = location.base().curvilinear_bounds(ipath)?;
            let law = location.base().law(ipath).ok_or_else(|| {
                "Standard_NoSuchObject: BRepFill_PipeShell::Prepare: Law(ipath)".to_string()
            })?;
            let loc = law.as_location_guide().ok_or_else(|| {
                "Standard_TypeMismatch: BRepFill_PipeShell::Prepare: GeomFill_LocationGuide"
                    .to_string()
            })?;
            // `Loc->Set(Sec, true, fs + f * Delta, fs + l * Delta, old_angle, angle)`.
            let angle = loc.set(
                Some(Arc::clone(&sec)),
                true,
                fs + f * delta,
                fs + l * delta,
                old_angle,
            )?;
            old_angle = angle;
        }
        Ok(())
    }

    /// `cxx:980-990`: whether every edge of the section is degenerated.
    fn section_is_degenerated(sec: &Section) -> bool {
        let Some(wire) = sec.wire() else {
            return true;
        };
        let mut wdeg = true;
        for shape in ShapeIterator::of_shape(&wire.0) {
            let degen = Edge::wrap(shape)
                .map(|e| BRepTool::is_degenerated(&e))
                .unwrap_or(false);
            wdeg = wdeg && degen;
        }
        wdeg
    }
}

/// `PerformTransition(Mode, Loc, angmin)` (`cxx:130-148`).
pub fn perform_transition(
    mode: TransitionStyle,
    loc: &mut dyn LocationLaw,
    ang_min: f64,
) -> Result<(), String> {
    loc.base_mut().delete_transform()?;
    if mode == TransitionStyle::Modified {
        loc.base_mut().transform_in_g0_law()?;
    } else {
        loc.base_mut().transform_in_compatible_law(ang_min)?;
    }
    Ok(())
}

/// `ComputeSection(W1, W2, p1, p2, Wres)` (`BRepFill_PipeShell.cxx:96-124`): the
/// intermediate section of a closed spine, built by making the two wires
/// compatible and reading the `BRepFill_NSections` law at `US`.
pub fn compute_section(w1: &Wire, w2: &Wire, p1: f64, p2: f64) -> Result<Wire, String> {
    // `cxx:103-106`.
    let params = [0.0, 1.0];
    let sections = [w1.0.clone(), w2.0.clone()];

    // `cxx:111-118`.
    let mut cw = CompatibleWires::with_sections(&sections);
    cw.set_percent(0.1);
    cw.perform(false)?;
    if !cw.is_done() {
        return Err("StdFail_NotDone: Uncompatible wires".to_string());
    }

    // `cxx:119-121`: `new BRepFill_NSections(CW.Shape(), EmptyTrsfs, SR, 0., 1.)`
    // then `SL->D0(US, Wres)` with `US = p1 / (p1 + p2)`.
    let law =
        NSections::with_transformations(cw.shape().to_vec(), Vec::new(), params.to_vec(), 0.0, 1.0, true)?;
    let us = p1 / (p1 + p2);
    let shape = law
        .d0_shape(us)
        .ok_or_else(|| "BRepFill_PipeShell::ComputeSection: no section".to_string())?;
    Wire::wrap(shape).ok_or_else(|| "BRepFill_PipeShell::ComputeSection: not a wire".to_string())
}
