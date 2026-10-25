//! `ChFi3d_Builder::PerformSetOfSurf` / `PerformSetOfKGen`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3298-3900`
//! (OCCT 8.0.0):
//!
//! - `PerformSetOfSurf` (`cxx:3882-3900`), the entry point: the solid index of
//!   the spine, `PerformSetOfKPart` when the spine is not split yet,
//!   `PerformSetOfKGen`, then `ChFi3d_MakeExtremities` in the non-simultaneous
//!   run,
//! - `PerformSetOfKGen` (`cxx:3298-3881`): the walk of every `ElSpine` of the
//!   spine through [`ChFi3dBuilder::perform_set_of_surf_on_el_spine`], then, in
//!   the non-simultaneous run, the repair of the twisted sections
//!   (`cxx:3318-3648`, [`twist`]), the update of the point tolerances
//!   (`cxx:3650-3754`, [`tolerances`]) and the update of the edge / new-face
//!   links of `myEVIMap` (`cxx:3755-3878`, [`tolerances::update_edge_face_links`]).
//!
//! Translation notes:
//! - `occ::handle<BRepTopAdaptor_TopolTool>` are the `Rc<RefCell<TopolTool>>`
//!   the crate passes down to the walk, created empty exactly like
//!   `cxx:3300-3301` (`new` without `Initialize`),
//! - the `NCollection_List<handle(ChFiDS_ElSpine)>` walk of `cxx:3304-3308`
//!   copies the list out of the spine, walks the copies and writes each
//!   `ElSpine` back into the list, because the port's
//!   [`ChFi3dBuilder::perform_set_of_surf_on_el_spine`] takes the guideline by
//!   `&mut Arc` and reaches `Arc::make_mut`,
//! - `Standard_Failure` (`cxx:3380`, `:3424`, `:3496`) is the `Err`,
//! - `DStr` accesses go through `myDS->ChangeDS()` at each use, exactly like
//!   the OCCT reference which holds a `TopOpeBRepDS_DataStructure&`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{ElSpine, StripeHandle};
use crate::chfi3d::builder_0_helpers::ch_fi3d_solid_index;

use super::surf_ops::ChFi3dSurfOps;
use super::ChFi3dBuilder;

mod tolerances;
mod twist;

impl ChFi3dBuilder {
    /// `PerformSetOfSurf(Stripe, Simul)` (`cxx:3882-3900`).
    pub fn perform_set_of_surf(
        &mut self,
        ops: &mut dyn ChFi3dSurfOps,
        stripe: &StripeHandle,
        simul: bool,
    ) -> Result<(), String> {
        // `TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();`
        // `const Handle(ChFiDS_Spine)& sp = Stripe->Spine();`
        // `int SI = ChFi3d_SolidIndex(sp, DStr, myESoMap, myEShMap);`
        // `Stripe->SetSolidIndex(SI);` (`cxx:3885-3888`).
        let sp = stripe
            .read()
            .map_err(|_| "PerformSetOfSurf: poisoned stripe".to_string())?
            .spine()
            .cloned()
            .ok_or_else(|| "PerformSetOfSurf: the stripe has no spine".to_string())?;
        let si = ch_fi3d_solid_index(
            &sp,
            &mut self.my_ds,
            &self.my_eso_map,
            &self.my_esh_map,
        )?;
        stripe
            .write()
            .map_err(|_| "PerformSetOfSurf: poisoned stripe".to_string())?
            .set_solid_index(si);

        // `if (!sp->SplitDone()) { PerformSetOfKPart(Stripe, Simul); }`
        // (`cxx:3889-3892`).
        let split_done = sp
            .read()
            .map_err(|_| "PerformSetOfSurf: poisoned spine".to_string())?
            .split_done();
        if !split_done {
            self.perform_set_of_k_part(ops, stripe, simul)?;
        }

        // `PerformSetOfKGen(Stripe, Simul);` (`cxx:3894`).
        self.perform_set_of_k_gen(ops, stripe, simul)?;

        // `if (!Simul) { ChFi3d_MakeExtremities(Stripe, DStr, myEFMap,
        // tolapp3d, tol2d); }` (`cxx:3896-3899`).
        if !simul {
            let (tol3d, tol2d) = (self.tolapp3d, self.tol2d);
            self.make_extremities(stripe, tol3d, tol2d)?;
        }
        Ok(())
    }

    /// `PerformSetOfKGen(Stripe, Simul)` (`cxx:3298-3881`).
    pub fn perform_set_of_k_gen(
        &mut self,
        ops: &mut dyn ChFi3dSurfOps,
        stripe: &StripeHandle,
        simul: bool,
    ) -> Result<(), String> {
        // `Handle(BRepTopAdaptor_TopolTool) It1 = new ...; It2 = new ...;`
        // (`cxx:3300-3301`).
        let it1: Rc<RefCell<TopolTool>> = Rc::new(RefCell::new(TopolTool::new()));
        let it2: Rc<RefCell<TopolTool>> = Rc::new(RefCell::new(TopolTool::new()));

        // `Handle(ChFiDS_Spine)& Spine = Stripe->ChangeSpine();`
        // `List(Handle(ChFiDS_ElSpine))& ll = Spine->ChangeElSpines();`
        // (`cxx:3302-3303`).
        let spine = stripe
            .read()
            .map_err(|_| "PerformSetOfKGen: poisoned stripe".to_string())?
            .spine()
            .cloned()
            .ok_or_else(|| "PerformSetOfKGen: the stripe has no spine".to_string())?;

        // `for (; ILES.More(); ILES.Next()) {
        //    PerformSetOfSurfOnElSpine(ILES.Value(), Stripe, It1, It2, Simul); }`
        // (`cxx:3305-3308`).
        let guides: Vec<Arc<ElSpine>> = spine
            .read()
            .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
            .el_spines()
            .to_vec();
        for (i, guide) in guides.iter().enumerate() {
            let mut hguide = Arc::clone(guide);
            self.perform_set_of_surf_on_el_spine(ops, &mut hguide, stripe, &it1, &it2, simul)?;
            // The list element is the handle `PerformSetOfSurfOnElSpine` wrote
            // back through `Change*` in OCCT, so the walked copy replaces it.
            spine
                .write()
                .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
                .el_spines_mut()[i] = hguide;
        }

        // `if (!Simul) { ... }` (`cxx:3309-3880`).
        if !simul {
            let periodic = spine
                .read()
                .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
                .is_periodic();

            // `Handle(NCollection_HSequence<...>)& HData =
            // Stripe->ChangeSetOfSurfData();` (`cxx:3312-3314`).
            let mut seq: Vec<Arc<crate::chfi_ds::SurfData>> = stripe
                .read()
                .map_err(|_| "PerformSetOfKGen: poisoned stripe".to_string())?
                .set_of_surf_data()
                .to_vec();

            // `It is attempted to reprocess the squares that bore.`
            // (`cxx:3322-3648`).
            twist::reprocess_twists(self, &spine, &mut seq)?;

            // `int len = SeqSurf.Length(); int last = len; if (periodic)
            // last++;` (`cxx:3315-3321`, re-read after the repair loop wrote
            // the sequence back).
            let len = seq.len() as i32;
            let last = if periodic { len + 1 } else { len };

            // `The tolerance of points is updated.` (`cxx:3650-3754`).
            tolerances::update_point_tolerances(self, &mut seq, len, last)?;

            // `The connections edge/new faces are updated.`
            // (`cxx:3755-3878`).
            tolerances::update_edge_face_links(self, &spine, &seq, last, periodic)?;

            // The sequence is the stripe's own `HData` in OCCT; the port holds
            // a copy, so it is written back once at the end.
            stripe
                .write()
                .map_err(|_| "PerformSetOfKGen: poisoned stripe".to_string())?
                .change_set_of_surf_data()
                .clear();
            stripe
                .write()
                .map_err(|_| "PerformSetOfKGen: poisoned stripe".to_string())?
                .change_set_of_surf_data()
                .extend(seq);
        }
        Ok(())
    }
}
