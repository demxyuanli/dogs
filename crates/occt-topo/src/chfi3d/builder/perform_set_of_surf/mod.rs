//! `ChFi3d_Builder::PerformSetOfSurfOnElSpine`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:2207-3003`
//! (OCCT 8.0.0).
//!
//! The driver walks one elementary spine of the guideline and, at every
//! section, asks the subclass for the band between the two current support
//! faces, then splices the computed `SurfData` into the stripe:
//!
//! 1. the guideline is extended by `prab = 0.01` at both ends when it is not
//!    periodic, and the `SurfData` strictly between `Guide.Previous()` and
//!    `Guide.Next()` are removed (`cxx:2233-2257`),
//! 2. the starting solution comes either from `StartSol(Stripe, HGuide, ...)`
//!    (`cxx:2290`, when there is no neighbouring `SurfData` at all) or from the
//!    two `StartSol(Spine, HS, ...)` calls that walk the last known contact
//!    onto the next faces (`cxx:2326-2360`),
//! 3. the `while (!fini)` loop (`cxx:2416-2988`) advances one section at a
//!    time, forward first and then, when `fini && complete`, backward from the
//!    start again,
//! 4. the extensions at both extremities are handled by `ChFi3d_BuildPlane`
//!    (`cxx:2514-2523`) and marked by `intf` / `intl` = 5,
//! 5. the slice is spliced with `InsertAfter` / `InsertBefore` and cleaned by
//!    `ChFi3d_Purge` (`cxx:2872-2927`, [`super::stripe_surgery`]).
//!
//! The four "contact" cases the section computation splits into (both faces on
//! an obstacle, one face on an obstacle, no obstacle) live in
//! [`obstacles`]; the splice lives in [`splice`].
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - the `ChFiDS_Spine` / `ChFiDS_Stripe` handles are the crate's shared-mutable
//!   [`SpineHandle`] / [`StripeHandle`], and `Guide` is `&mut Arc<ElSpine>`
//!   because the function rewrites its parameter range,
//! - `NCollection_Sequence<handle(ChFiDS_SurfData)> SeqSD` is
//!   `Vec<Arc<SurfData>>`; `SD` (the handle appended at `cxx:2564-2565`) is the
//!   first element, and the `SimulSurf` arms work through `Arc::make_mut`,
//! - `done` is the `ChFi3d_Builder` member (`hxx:762`) the corner drivers write;
//!   this driver only reads it before the splice and sets it on the failure
//!   arm (`cxx:2847-2869`), exactly as OCCT,
//! - the `throw Standard_Failure` arms are `Err`,
//! - `OffsetHGuide` is the element of `Spine->ChangeOffsetElSpines()` in OCCT,
//!   so its range writes are visible to the spine; the port keeps the list
//!   index and [`Walk::flush_offset`] puts the updated `Arc<ElSpine>` back
//!   after every write (`cxx:2215-2230`),
//! - `ChFi3d_Builder::PerformSetOfSurfOnElSpine` has no ported caller yet
//!   (`PerformOneCorner` / `PerformTwoCornerbyInter` / `PerformMoreSurfdata`
//!   are still UNPORTED, see `builder/mod.rs`), hence the `dead_code` allow.

#![allow(dead_code)]

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_math::MathVector;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{ChamfMode, ElSpine, SpineHandle, State, StripeHandle, SurfData};
use crate::shape::TopoShape;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::call_perform_surf::initialize_topol_tool;
use super::stripe_surgery::remove_sd;
use super::surf_ops::ChFi3dSurfOps;
use super::ChFi3dBuilder;

mod obstacles;
mod splice;
mod start;
mod step;

/// What one turn of the `while (!fini)` loop asks the driver to do.
pub(super) enum StepOutcome {
    /// `fini` stayed `false`: run another section.
    Continue,
    /// `fini` is `true`: leave the loop and restore the guideline range
    /// (`cxx:2989-2999`).
    Done,
    /// The `return;` of the periodic same-arc test (`cxx:2962-2965`): leave the
    /// function **without** restoring the guideline range, as OCCT does.
    EarlyReturn,
}

/// The state the walk carries between sections (`cxx:2252-2276`).
///
/// OCCT keeps these as locals of the one big function; the port keeps them in
/// one struct because the section computation, the three obstacle arms and the
/// splice are translated in their own files and all of them read / write the
/// same set of variables.
pub(super) struct Walk<'a> {
    /// `ChFi3d_Builder` itself, for `myDS->ChangeDS().AddShape`, `tolesp`,
    /// `fleche` and the five `StartSol` / `CallPerformSurf` drivers.
    pub(super) builder: &'a mut ChFi3dBuilder,
    /// The protected virtuals the section computation dispatches to.
    pub(super) ops: &'a mut dyn ChFi3dSurfOps,
    /// `Stripe`
    pub(super) stripe: &'a StripeHandle,
    /// `Spine = Stripe->ChangeSpine()`
    pub(super) spine: SpineHandle,
    /// `ChFiDS_ElSpine& Guide = *HGuide`
    pub(super) guide: &'a mut Arc<ElSpine>,
    /// `OffsetHGuide` (`cxx:2215-2230`).
    pub(super) offset: Option<Arc<ElSpine>>,
    /// Position of `OffsetHGuide` inside `Spine->ChangeOffsetElSpines()`.
    /// OCCT holds the very element of that list, so every
    /// `OffsetHGuide->FirstParameter(...)` writes through to the spine; the
    /// port keeps the index and [`Walk::flush_offset`] puts the updated
    /// `ElSpine` back into the list.
    pub(super) offset_index: Option<usize>,
    /// `It1` / `It2`
    pub(super) it1: &'a Rc<RefCell<TopolTool>>,
    pub(super) it2: &'a Rc<RefCell<TopolTool>>,
    /// `Simul`
    pub(super) simul: bool,

    /// `wf` / `wl`: the domain of the guideline, possibly extended.
    pub(super) wf: f64,
    pub(super) wl: f64,
    /// `wfsav` / `wlsav`: the unextended domain, saved at `cxx:2234`.
    pub(super) wf_sav: f64,
    pub(super) wl_sav: f64,
    /// `locfleche = (wl - wf) * fleche`, computed before the extension.
    pub(super) locfleche: f64,
    /// `nbed`, rewritten to the edge span of the section at `cxx:2390-2395`.
    pub(super) nbed: i32,
    /// `lastedlastp = Spine->LastParameter(nbed)` (`cxx:2250`).
    pub(super) lastedlastp: f64,

    /// `ref` / `raf`: `Guide.Previous()` / `Guide.Next()`.
    pub(super) ref_sd: Option<Arc<SurfData>>,
    pub(super) raf: Option<Arc<SurfData>>,

    /// `HS1` / `HS2`
    pub(super) hs1: SurfaceAdaptor,
    pub(super) hs2: SurfaceAdaptor,
    /// `HS3` / `HS4`: the second-choice faces `StartSol` hands back (`cxx:2261`).
    pub(super) hs3: Option<SurfaceAdaptor>,
    pub(super) hs4: Option<SurfaceAdaptor>,
    /// `HSref1` / `HSref2`: the reference faces, non-null handles in OCCT
    /// (`cxx:2263-2264`).
    pub(super) hsref1: SurfaceAdaptor,
    pub(super) hsref2: SurfaceAdaptor,
    /// `HC1` / `HC2`: the contact arcs, null handles in OCCT (`cxx:2266`).
    pub(super) hc1: Option<Curve2dAdaptor>,
    pub(super) hc2: Option<Curve2dAdaptor>,
    /// `HCref1` / `HCref2`, non-null handles in OCCT (`cxx:2267-2268`).
    pub(super) hcref1: Curve2dAdaptor,
    pub(super) hcref2: Curve2dAdaptor,
    /// `decroch1` / `decroch2` (`cxx:2269`). The two `StartSol` overload the
    /// driver calls take `decroch` by value (`cxx:1179`), so these stay `false`.
    pub(super) decroch1: bool,
    pub(super) decroch2: bool,
    /// `RecP1` / `RecS1` / `RecRst1` and the `2` family (`cxx:2268-2269`).
    pub(super) rec_p1: bool,
    pub(super) rec_s1: bool,
    pub(super) rec_rst1: bool,
    pub(super) obstacleon1: bool,
    pub(super) rec_p2: bool,
    pub(super) rec_s2: bool,
    pub(super) rec_rst2: bool,
    pub(super) obstacleon2: bool,

    /// `pp1` / `pp2`: the two support parameters of the current section, and
    /// `pp3` / `pp4`: the same on the second-choice faces.
    pub(super) pp1: GpPnt2d,
    pub(super) pp2: GpPnt2d,
    pub(super) pp3: GpPnt2d,
    pub(super) pp4: GpPnt2d,
    /// `w1` / `w2`: the guideline parameters of the two contacts.
    pub(super) w1: f64,
    pub(super) w2: f64,
    /// `Soldep` (`cxx:2272`), `SoldepCS` / `SoldepCC` (`cxx:2273-2274`).
    pub(super) soldep: MathVector,
    pub(super) soldep_cs: MathVector,
    pub(super) soldep_cc: MathVector,

    /// `forward`, `Inside`, `First`, `Last` (`cxx:2277-2281`).
    pub(super) forward: bool,
    pub(super) inside: bool,
    pub(super) first: f64,
    pub(super) last: f64,
    /// `Ok1` / `Ok2` (`cxx:2282`), and `Vref` (`cxx:2284`).
    pub(super) ok1: bool,
    pub(super) ok2: bool,
    pub(super) vref: Option<TopoShape>,

    /// `Firstsov` (`cxx:2412`), `intf` / `intl` (`cxx:2413`), `complete`
    /// (`cxx:2387`) and `MaxStep` (`cxx:2411`).
    pub(super) firstsov: f64,
    pub(super) intf: i32,
    pub(super) intl: i32,
    pub(super) complete: bool,
    pub(super) max_step: f64,
}

impl ChFi3dBuilder {
    /// `PerformSetOfSurfOnElSpine(HGuide, Stripe, It1, It2, Simul)`
    /// (`cxx:2207-3003`).
    ///
    /// Walks the guideline `HGuide` and fills `Stripe`'s `SurfData` sequence.
    /// `Err` carries the three `throw Standard_Failure` arms of the function and
    /// the failure modes the port refuses to dereference (a stripe without a
    /// spine, a `StartSol` that found no solution, an empty sequence where OCCT
    /// would read `First()`).
    ///
    /// The three `ChFi3d_Builder_NotImp.cxx` edge-contact overloads answer
    /// `Err("... Not Implemented")` until `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`
    /// are ported, so the two obstacle arms propagate that `Err`.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_set_of_surf_on_el_spine(
        &mut self,
        ops: &mut dyn ChFi3dSurfOps,
        hguide: &mut Arc<ElSpine>,
        stripe: &StripeHandle,
        it1: &Rc<RefCell<TopolTool>>,
        it2: &Rc<RefCell<TopolTool>>,
        simul: bool,
    ) -> Result<(), String> {
        // `occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine();`
        // (`cxx:2214`).
        let spine = {
            let guard = stripe
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned stripe".to_string())?;
            guard.spine().cloned()
        };
        let Some(spine) = spine else {
            return Err("PerformSetOfSurfOnElSpine: the stripe has no spine".into());
        };

        // `if (Spine->Mode() == ChFiDS_ConstThroatWithPenetrationChamfer) { ...
        // OffsetHGuide = ILES_offset.Value(); }` (`cxx:2217-2230`).
        //
        // The two `NCollection_List` iterators advance in step; a shorter
        // offset list would raise `Standard_OutOfRange` on `ILES_offset.Value()`
        // and the port answers "no match" instead.
        let offset = {
            let guard = spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
            if guard.mode() != ChamfMode::ConstThroatWithPenetrationChamfer {
                None
            } else {
                let mut found = None;
                for (i, hel) in guard.el_spines().iter().enumerate() {
                    if Arc::ptr_eq(hel, hguide) {
                        found = guard.offset_el_spines().get(i).map(|o| (i, o.clone()));
                    }
                }
                found
            }
        };
        let offset_index = offset.as_ref().map(|(i, _)| *i);
        let mut offset = offset.map(|(_, o)| o);

        // `double wf = Guide.FirstParameter(); double wl =
        // Guide.LastParameter(); double locfleche = (wl - wf) * fleche;
        // double wfsav = wf, wlsav = wl;` (`cxx:2233-2236`).
        let wf = hguide.first_parameter();
        let wl = hguide.last_parameter();
        let locfleche = (wl - wf) * self.fleche;
        let wf_sav = wf;
        let wl_sav = wl;

        // `if (!Guide.IsPeriodic()) { const double prab = 0.01; ...
        // }` (`cxx:2237-2248`): the ElSpine is artificially extended to help
        // rsnld.
        if !hguide.is_periodic() {
            let prab = 0.01;
            Arc::make_mut(hguide).set_first_parameter(wf - prab * (wl - wf));
            Arc::make_mut(hguide).set_last_parameter(wl + prab * (wl - wf));
            if let (Some(current), Some(index)) = (offset.clone(), offset_index) {
                let mut updated = current;
                Arc::make_mut(&mut updated).set_first_parameter(wf - prab * (wl - wf));
                Arc::make_mut(&mut updated).set_last_parameter(wl + prab * (wl - wf));
                // `OffsetHGuide` is the list element itself in OCCT.
                if let Ok(mut guard) = spine.write() {
                    if let Some(slot) = guard.offset_el_spines_mut().get_mut(index) {
                        *slot = updated.clone();
                    }
                }
                offset = Some(updated);
            }
        }

        // `int ii, nbed = Spine->NbEdges(); double lastedlastp =
        // Spine->LastParameter(nbed);` (`cxx:2249-2250`).
        let (nbed, lastedlastp) = {
            let guard = spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
            let nbed = guard.nb_edges() as i32;
            let lastedlastp = guard
                .last_parameter_of(nbed)
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: empty spine".to_string())?;
            (nbed, lastedlastp)
        };

        // `occ::handle<ChFiDS_SurfData> ref = Guide.Previous(); ...
        // occ::handle<ChFiDS_SurfData> raf = Guide.Next();
        // RemoveSD(Stripe, ref, raf);` (`cxx:2255-2257`).
        let ref_sd = hguide.previous().cloned();
        let raf = hguide.next().cloned();
        remove_sd(stripe, ref_sd.as_ref(), raf.as_ref());

        // `HS1` / `HS2` are `new BRepAdaptor_Surface()` (`cxx:2259-2260`): a
        // null face until `StartSol` / `ConexFaces` initializes them.
        let mut walk = Walk {
            builder: self,
            ops,
            stripe,
            spine,
            guide: hguide,
            offset,
            offset_index,
            it1,
            it2,
            simul,
            wf,
            wl,
            wf_sav,
            wl_sav,
            locfleche,
            nbed,
            lastedlastp,
            ref_sd,
            raf,
            hs1: SurfaceAdaptor::new(),
            hs2: SurfaceAdaptor::new(),
            hs3: None,
            hs4: None,
            hsref1: SurfaceAdaptor::new(),
            hsref2: SurfaceAdaptor::new(),
            hc1: None,
            hc2: None,
            hcref1: Curve2dAdaptor::new(),
            hcref2: Curve2dAdaptor::new(),
            decroch1: false,
            decroch2: false,
            rec_p1: false,
            rec_s1: false,
            rec_rst1: false,
            obstacleon1: false,
            rec_p2: false,
            rec_s2: false,
            rec_rst2: false,
            obstacleon2: false,
            pp1: GpPnt2d::new(0.0, 0.0),
            pp2: GpPnt2d::new(0.0, 0.0),
            pp3: GpPnt2d::new(0.0, 0.0),
            pp4: GpPnt2d::new(0.0, 0.0),
            w1: 0.0,
            w2: 0.0,
            // `math_Vector Soldep(1, 4); SoldepCS(1, 3); SoldepCC(1, 2)`
            // (`cxx:2272-2274`).
            soldep: MathVector::new(1, 4),
            soldep_cs: MathVector::new(1, 3),
            soldep_cc: MathVector::new(1, 2),
            forward: true,
            inside: false,
            first: wf,
            last: wl,
            ok1: true,
            ok2: true,
            vref: None,
            firstsov: 0.0,
            intf: 0,
            intl: 0,
            complete: false,
            max_step: 0.0,
        };

        walk.start()?;
        loop {
            match walk.step()? {
                StepOutcome::Continue => {}
                StepOutcome::Done => break,
                StepOutcome::EarlyReturn => return Ok(()),
            }
        }

        // `if (!Guide.IsPeriodic()) { Guide.FirstParameter(wfsav);
        // Guide.LastParameter(wlsav); ... }` (`cxx:2989-2999`).
        if !walk.guide.is_periodic() {
            let (wf_sav, wl_sav) = (walk.wf_sav, walk.wl_sav);
            walk.set_guide_first(wf_sav);
            walk.set_guide_last(wl_sav);
        }
        Ok(())
    }
}

impl<'a> Walk<'a> {
    /// `Guide.FirstParameter(P)` together with the `OffsetHGuide` twin
    /// (`cxx:2367-2386`, `:2436-2441`, `:2533-2548`).
    fn set_guide_first(&mut self, p: f64) {
        Arc::make_mut(self.guide).set_first_parameter(p);
        if let Some(offset) = self.offset.as_mut() {
            Arc::make_mut(offset).set_first_parameter(p);
        }
        self.flush_offset();
    }

    /// `Guide.LastParameter(P)` together with the `OffsetHGuide` twin.
    fn set_guide_last(&mut self, p: f64) {
        Arc::make_mut(self.guide).set_last_parameter(p);
        if let Some(offset) = self.offset.as_mut() {
            Arc::make_mut(offset).set_last_parameter(p);
        }
        self.flush_offset();
    }

    /// Writes `OffsetHGuide` back into `Spine->ChangeOffsetElSpines()`
    /// (`cxx:2215-2230`). In OCCT `OffsetHGuide` *is* the element of that
    /// list, so the `FirstParameter` / `LastParameter` writes of this function
    /// are visible to the spine afterwards; `Arc::make_mut` splits the shared
    /// `Arc` instead, hence the explicit write-back.
    fn flush_offset(&mut self) {
        let (Some(offset), Some(index)) = (self.offset.as_ref(), self.offset_index) else {
            return;
        };
        let offset = offset.clone();
        if let Ok(mut guard) = self.spine.write() {
            if let Some(slot) = guard.offset_el_spines_mut().get_mut(index) {
                *slot = offset;
            }
        }
    }

    /// `Spine->SetErrorStatus(S)` (`cxx:919-922`).
    fn set_spine_error(&self, status: crate::chfi_ds::ErrorStatus) -> Result<(), String> {
        let mut guard = self
            .spine
            .write()
            .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
        guard.set_error_status(status);
        Ok(())
    }

    /// `SD->ChangeIndexOfS1(DStr.AddShape(HS1->Face()))` and the three
    /// companions (`cxx:2662-2664`, `:2683-2686`, `:2734-2748`,
    /// `:2878-2886`): the face / edge indices the build stage resolves the
    /// `SurfData` through.
    fn set_indices(&mut self, sd: &mut SurfData) -> Result<(), String> {
        let face1 = self
            .hs1
            .face()
            .cloned()
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: S1 without face".to_string())?;
        let index1 = self.builder.my_ds.change_ds().add_shape(&face1.0);
        sd.set_index_of_s1(index1);
        if self.obstacleon1 {
            let edge1 = self
                .hc1
                .as_ref()
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: obstacle on S1 without arc".to_string())?
                .edge()
                .clone();
            let index_c1 = self.builder.my_ds.change_ds().add_shape(&edge1.0);
            sd.set_index_of_c1(index_c1);
        }
        let face2 = self
            .hs2
            .face()
            .cloned()
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: S2 without face".to_string())?;
        let index2 = self.builder.my_ds.change_ds().add_shape(&face2.0);
        sd.set_index_of_s2(index2);
        if self.obstacleon2 {
            let edge2 = self
                .hc2
                .as_ref()
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: obstacle on S2 without arc".to_string())?
                .edge()
                .clone();
            let index_c2 = self.builder.my_ds.change_ds().add_shape(&edge2.0);
            sd.set_index_of_c2(index_c2);
        }
        Ok(())
    }
}
