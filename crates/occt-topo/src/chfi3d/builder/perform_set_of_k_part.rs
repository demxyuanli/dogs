//! `ChFi3d_Builder::PerformSetOfKPart`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3004-3281`
//! (OCCT 8.0.0).
//!
//! The first of the two walks `PerformSetOfSurf` runs on a stripe
//! (`ChFi3d_Builder_2.cxx:3889-3894`): it scans the edges of the spine, asks
//! [`ch_fi3d_k_particular`] which of them are a plane / (plane | cylinder |
//! cone) corner along a line or a circle, builds those corners with
//! `ChFiKPart_ComputeData::Compute`, splits them with `SplitKPart`, and turns
//! the resulting `SurfData` sequence into the `ElSpine` list `PerformSetOfKGen`
//! walks. Every edge that is not a particular corner is left to the general
//! walk.
//!
//! Divergences that are not heuristics:
//! - `occ::handle<ChFiDS_Stripe>&` / `occ::handle<ChFiDS_Spine>&` are the
//!   shared-mutable [`StripeHandle`] / [`SpineHandle`] of the crate,
//! - `NCollection_Sequence<handle(ChFiDS_SurfData)>` is `Vec<Arc<SurfData>>`,
//!   so `LSD.ChangeValue(j)` is an `Arc` clone,
//! - `gp_Pnt&` / `gp_Vec&` out-parameters of `TgtKP` are a returned tuple,
//! - OCCT writes `CurrentHE->ChangeNext()` *after* `AppendElSpine` because the
//!   list stores the same handle; the port performs the same write before the
//!   move, where `Arc::make_mut` still has a unique owner and therefore lands in
//!   the appended element all the same,
//! - `Standard_Failure` (`cxx:3380`, `:3424`) and the null adaptor
//!   dereferences (`cxx:3057`, `:3060`) are `Err`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpPnt, GpVec};

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_k_particular, ch_fi3d_perform_el_spine,
};
use crate::chfi3d::builder_2_helpers::tgt_kp;
use crate::chfi3d::builder_helpers::in_period;
use crate::chfi3d::chfi_k_part::ch_fi_kpart_compute_data;
use crate::chfi3d::side::next_side;
use crate::chfi_ds::{ChamfMode, ElSpine, SpineHandle, StripeHandle, SurfData};

use super::call_perform_surf::initialize_topol_tool;
use super::split_k_part::split_k_part;
use super::surf_ops::ChFi3dSurfOps;
use super::ChFi3dBuilder;

/// `Spine->Period()`, with the OCCT `Standard_Failure` on a non-periodic spine
/// (`ChFiDS_Spine.cxx:1500-1510`) as the `Err`.
fn spine_period(spine: &SpineHandle) -> Result<f64, String> {
    spine
        .read()
        .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
        .period()
        .ok_or_else(|| "PerformSetOfKPart: ChFiDS_Spine::Period() on a non-periodic spine".into())
}

/// `Spine->IsPeriodic()` (`ChFiDS_Spine.lxx:64-67`).
fn spine_is_periodic(spine: &SpineHandle) -> Result<bool, String> {
    Ok(spine
        .read()
        .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
        .is_periodic())
}

/// `Spine->AppendElSpine(HES)` / `Spine->AppendOffsetElSpine(HES)`.
fn append_el_spine(spine: &SpineHandle, es: Arc<ElSpine>, offset: bool) -> Result<(), String> {
    let mut guard = spine
        .write()
        .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?;
    if offset {
        guard.append_offset_el_spine(es);
    } else {
        guard.append_el_spine(es);
    }
    Ok(())
}

impl ChFi3dBuilder {
    /// `PerformSetOfKPart(Stripe, Simul)` (`cxx:3004-3281`).
    ///
    /// `ops` is the [`ChFi3dSurfOps`] the `SimulKPart` dispatch goes to; `Err`
    /// covers the OCCT throws, an unbound spine / face pair, and the two
    /// UNPORTED call-outs ([`ch_fi_kpart_compute_data`] answers "no data"
    /// instead of failing, [`split_k_part`] answers `Err`).
    pub fn perform_set_of_k_part(
        &mut self,
        ops: &mut dyn ChFi3dSurfOps,
        stripe: &StripeHandle,
        simul: bool,
    ) -> Result<(), String> {
        // `occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine();`
        // (`cxx:3007`).
        let spine = stripe
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned stripe".to_string())?
            .spine()
            .cloned()
            .ok_or_else(|| "PerformSetOfKPart: the stripe has no spine".to_string())?;

        // `Stripe->Reset(); HData = new ...; SeqSurf = HData->ChangeSequence();`
        // (`cxx:3012-3017`). `Reset` already nullifies `myHdata`, and the port's
        // `Stripe::reset` clears the vector, so the extra `clear` is the empty
        // sequence `new NCollection_HSequence` produces.
        {
            let mut guard = stripe
                .write()
                .map_err(|_| "PerformSetOfKPart: poisoned stripe".to_string())?;
            guard.reset();
            guard.change_set_of_surf_data().clear();
        }

        // `StripeOrientations(Spine, RefOr1, RefOr2, RefChoix);
        // Stripe->OrientationOnFace1(RefOr1); ... Stripe->Choix(RefChoix);`
        // (`cxx:3019-3022`).
        let (ref_or1, ref_or2, ref_choix) = self.stripe_orientations(&spine)?;
        {
            let mut guard = stripe
                .write()
                .map_err(|_| "PerformSetOfKPart: poisoned stripe".to_string())?;
            guard.set_orientation_on_face1(ref_or1);
            guard.set_orientation_on_face2(ref_or2);
            guard.set_choix(ref_choix);
        }

        // `occ::handle<BRepTopAdaptor_TopolTool> It1 = new
        // BRepTopAdaptor_TopolTool();` (`cxx:3024-3025`).
        let it1: Rc<RefCell<TopolTool>> = Rc::new(RefCell::new(TopolTool::new()));
        let it2: Rc<RefCell<TopolTool>> = Rc::new(RefCell::new(TopolTool::new()));

        // `double WFirst, WLast = 0.; gp_Vec TFirst, TLast, TEndPeriodic;
        // gp_Pnt PFirst, PLast, PEndPeriodic;` (`cxx:3027-3029`). OCCT leaves
        // `PLast` / `TLast` uninitialized until the first `TgtKP` writes them;
        // the port keeps `Option` so an early read is an `Err` rather than a
        // garbage point.
        let mut w_first: f64;
        let mut w_last: f64 = 0.0;
        let mut p_first_t_first: (GpPnt, GpVec);
        let mut p_last_t_last: Option<(GpPnt, GpVec)>;
        let mut p_end_periodic_t_end_periodic: (GpPnt, GpVec);
        // `bool intf = false, intl = false;` (`cxx:3030`). `intf` is assigned
        // before its only read, so its OCCT initializer is dead; `intl` is read
        // after the loop and keeps the `false`.
        let mut intf: bool;
        let mut intl = false;

        // `Spine->D1(Spine->FirstParameter(), PFirst, TFirst);
        // CurrentHE->FirstParameter(Spine->FirstParameter());
        // CurrentHE->SetFirstPointAndTgt(PFirst, TFirst); ...` (`cxx:3036-3039`).
        let first_par = spine
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .first_parameter();
        p_first_t_first = spine
            .write()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .d1(first_par)
            .ok_or_else(|| {
                "PerformSetOfKPart: ChFiDS_Spine::D1 on the first parameter".to_string()
            })?;
        let mut current_he = Arc::new(ElSpine::new());
        let mut current_offset_he = Arc::new(ElSpine::new());
        {
            let he = Arc::make_mut(&mut current_he);
            he.set_first_parameter(first_par);
            he.set_first_point_and_tgt(p_first_t_first.0.clone(), p_first_t_first.1.clone());
            let ohe = Arc::make_mut(&mut current_offset_he);
            ohe.set_first_parameter(first_par);
            ohe.set_first_point_and_tgt(p_first_t_first.0.clone(), p_first_t_first.1.clone());
        }

        // `bool YaKPart = false; int iedgelastkpart = 0;` (`cxx:3042-3043`).
        let mut ya_kpart = false;
        let mut iedge_last_kpart = 0usize;

        // `double WStartPeriodic = 0.; double WEndPeriodic =
        // Spine->LastParameter(Spine->NbEdges()); Spine->D1(WEndPeriodic,
        // PEndPeriodic, TEndPeriodic);` (`cxx:3045-3047`).
        let mut w_start_periodic = 0.0f64;
        let nb_edges = spine
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .nb_edges();
        if nb_edges == 0 {
            return Err("PerformSetOfKPart: ChFiDS_Spine::Edges on an empty spine".into());
        }
        let mut w_end_periodic = spine
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .last_parameter_of(nb_edges as i32)
            .ok_or_else(|| {
                "PerformSetOfKPart: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
            })?;
        p_end_periodic_t_end_periodic = spine
            .write()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .d1(w_end_periodic)
            .ok_or_else(|| {
                "PerformSetOfKPart: ChFiDS_Spine::D1 on the last parameter".to_string()
            })?;

        let tolesp = self.tolesp;

        // `for (int iedge = 1; iedge <= Spine->NbEdges(); iedge++)`
        // (`cxx:3050`).
        for iedge in 1..=nb_edges {
            // `ConexFaces(Spine, iedge, HS1, HS2);` (`cxx:3053`).
            let (hs1, hs2) = self
                .cone_x_faces(&spine, iedge)
                .ok_or_else(|| format!("PerformSetOfKPart: ConexFaces failed on edge {iedge}"))?;

            // `if (ChFi3d_KParticular(Spine, iedge, *HS1, *HS2))` (`cxx:3055`).
            if !ch_fi3d_k_particular(&spine, iedge as i32, &hs1, &hs2) {
                continue;
            }

            // `intf = ((iedge == 1) && !Spine->IsPeriodic()); intl = ((iedge ==
            // Spine->NbEdges()) && !Spine->IsPeriodic());` (`cxx:3057-3058`).
            let is_periodic = spine_is_periodic(&spine)?;
            intf = (iedge == 1) && !is_periodic;
            intl = (iedge == nb_edges) && !is_periodic;

            // `Or1 = HS1->Face().Orientation(); Or2 = HS2->Face().Orientation();
            // ChFi3d::NextSide(Or1, Or2, RefOr1, RefOr2, RefChoix);`
            // (`cxx:3059-3062`).
            let mut or1 = hs1
                .face()
                .map(|f| f.0.orientation())
                .ok_or_else(|| "PerformSetOfKPart: HS1->Face() on a null adaptor".to_string())?;
            let mut or2 = hs2
                .face()
                .map(|f| f.0.orientation())
                .ok_or_else(|| "PerformSetOfKPart: HS2->Face() on a null adaptor".to_string())?;
            next_side(&mut or1, &mut or2, ref_or1, ref_or2, ref_choix);

            // `It1->Initialize(HSon1); It2->Initialize(HSon2);` (`cxx:3066-3067`).
            initialize_topol_tool(&it1, &hs1)?;
            initialize_topol_tool(&it2, &hs2)?;

            // `occ::handle<ChFiDS_SurfData> SD = new ChFiDS_SurfData();
            // NCollection_Sequence<...> LSD;` (`cxx:3069-3070`).
            let mut sd = SurfData::new();
            let mut lsd: Vec<Arc<SurfData>> = Vec::new();

            // `if (!ChFiKPart_ComputeData::Compute(DStr, SD, HS1, HS2, Or1, Or2,
            // Spine, iedge)) { } else if (!SplitKPart(SD, LSD, Spine, iedge,
            // HS1, It1, HS2, It2, intf, intl)) { LSD.Clear(); } else {
            // iedgelastkpart = iedge; }` (`cxx:3071-3081`).
            let computed = ch_fi_kpart_compute_data(
                self.my_ds.change_ds(),
                &mut sd,
                &hs1,
                &hs2,
                or1,
                or2,
                &spine,
                iedge as i32,
            )?;
            if !computed {
                // The empty arm of OCCT: no data, nothing is appended.
            } else {
                let sd = Arc::new(sd);
                if !split_k_part(
                    &sd,
                    &mut lsd,
                    &spine,
                    iedge as i32,
                    &hs1,
                    &it1,
                    &hs2,
                    &it2,
                    &mut intf,
                    &mut intl,
                )? {
                    lsd.clear();
                } else {
                    iedge_last_kpart = iedge;
                }
            }

            // `if (Spine->IsPeriodic()) { ... debug provisory for SD that arrive
            // in disorder ... }` (`cxx:3082-3121`).
            if is_periodic {
                let nbsd = lsd.len();
                let period = spine_period(&spine)?;
                let mut wfp = w_start_periodic;
                let mut wlp = w_end_periodic;
                // `if (!YaKPart && nbsd > 0) { ... }` (`cxx:3087-3102`).
                if !ya_kpart && nbsd > 0 {
                    let first_sd = &lsd[0];
                    let wwf = in_period(first_sd.first_spine_param(), wfp, wlp, tolesp);
                    let mut wwl = in_period(first_sd.last_spine_param(), wfp, wlp, tolesp);
                    if wwl <= wwf + tolesp {
                        wwl += period;
                    }
                    wfp = wwf;
                    wlp = wfp + period;
                    let _ = wwl;
                }
                // `for (int j = 1; j < nbsd; j++) for (int k = j + 1; k <= nbsd;
                // k++) { ... }` (`cxx:3104-3121`).
                let mut j = 0usize;
                while j + 1 < nbsd {
                    let jwf = in_period(lsd[j].first_spine_param(), wfp, wlp, tolesp);
                    let mut k = j + 1;
                    while k < nbsd {
                        let kwf = in_period(lsd[k].first_spine_param(), wfp, wlp, tolesp);
                        if kwf < jwf {
                            lsd.swap(j, k);
                        }
                        k += 1;
                    }
                    j += 1;
                }
            }

            // `NCollection_List<int> li; for (int j = 1; j <= LSD.Length(); j++)
            // { ... }` (`cxx:3122-3217`).
            let mut li: Vec<i32> = Vec::new();
            for j in 0..lsd.len() {
                // `occ::handle<ChFiDS_SurfData>& curSD = LSD.ChangeValue(j);`
                // (`cxx:3125`).
                let cur_sd = lsd[j].clone();

                // `if (Simul) { SimulKPart(curSD); }` (`cxx:3126-3129`).
                if simul {
                    let mut sd_mut = (*cur_sd).clone();
                    ops.simul_k_part(&mut sd_mut)?;
                }

                // `SeqSurf.Append(curSD);` (`cxx:3130`).
                stripe
                    .write()
                    .map_err(|_| "PerformSetOfKPart: poisoned stripe".to_string())?
                    .change_set_of_surf_data()
                    .push(cur_sd.clone());

                // `if (!Simul) { li.Append(curSD->Surf()); }` (`cxx:3132-3135`).
                if !simul {
                    li.push(cur_sd.surf());
                }

                // `WFirst = LSD.Value(j)->FirstSpineParam(); WLast =
                // LSD.Value(j)->LastSpineParam(); if (Spine->IsPeriodic()) {
                // WFirst = ChFi3d_InPeriod(...); ... }` (`cxx:3136-3143`).
                w_first = cur_sd.first_spine_param();
                w_last = cur_sd.last_spine_param();
                if is_periodic {
                    w_first = in_period(w_first, w_start_periodic, w_end_periodic, tolesp);
                    w_last = in_period(w_last, w_start_periodic, w_end_periodic, tolesp);
                    if w_last <= w_first + tolesp {
                        w_last += spine_period(&spine)?;
                    }
                }

                // `TgtKP(LSD.Value(j), Spine, iedge, true, PFirst, TFirst);
                // TgtKP(LSD.Value(j), Spine, iedge, false, PLast, TLast);`
                // (`cxx:3145-3146`).
                p_first_t_first = tgt_kp(&cur_sd, &spine, iedge as i32, true)
                    .ok_or_else(|| "PerformSetOfKPart: TgtKP on the first section".to_string())?;
                p_last_t_last = Some(
                    tgt_kp(&cur_sd, &spine, iedge as i32, false).ok_or_else(|| {
                        "PerformSetOfKPart: TgtKP on the last section".to_string()
                    })?,
                );
                let (p_first, t_first) = p_first_t_first.clone();
                let (p_last, t_last) = p_last_t_last.clone().ok_or_else(|| {
                    "PerformSetOfKPart: TgtKP on the last section".to_string()
                })?;

                // `if (!YaKPart) { ... } else { ... }` (`cxx:3149-3214`).
                if !ya_kpart {
                    if is_periodic {
                        // `WStartPeriodic = WFirst; WEndPeriodic = WStartPeriodic
                        // + Spine->Period(); WLast = ElCLib::InPeriod(WLast,
                        // WStartPeriodic, WEndPeriodic); ...; Spine->
                        // SetFirstParameter(WStartPeriodic); ...`
                        // (`cxx:3150-3168`).
                        let period = spine_period(&spine)?;
                        w_start_periodic = w_first;
                        w_end_periodic = w_start_periodic + period;
                        // `WLast = ElCLib::InPeriod(WLast, WStartPeriodic,
                        // WEndPeriodic);` (`cxx:3155`).
                        w_last = clib::in_period(w_last, w_start_periodic, w_end_periodic);
                        if w_last <= w_first + tolesp {
                            w_last += period;
                        }
                        p_end_periodic_t_end_periodic = (p_first.clone(), t_first.clone());
                        {
                            let mut guard = spine.write().map_err(|_| {
                                "PerformSetOfKPart: poisoned spine".to_string()
                            })?;
                            guard.set_first_parameter(w_start_periodic);
                            guard.set_last_parameter(w_end_periodic);
                        }
                    } else if !intf || (iedge > 1) {
                        // `Spine->SetFirstTgt(std::min(0., WFirst));
                        // CurrentHE->LastParameter(WFirst); ...; Spine->
                        // AppendElSpine(CurrentHE); ...` (`cxx:3169-3181`).
                        spine
                            .write()
                            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                            .set_first_tgt(0.0f64.min(w_first))
                            .ok_or_else(|| {
                                "PerformSetOfKPart: SetFirstTgt on a periodic spine".to_string()
                            })?;
                        {
                            let he = Arc::make_mut(&mut current_he);
                            he.set_last_parameter(w_first);
                            he.set_last_point_and_tgt(p_first.clone(), t_first.clone());
                            // `CurrentHE->ChangeNext() = LSD.Value(j);`
                            // (`cxx:3174`) — written before the move, see the
                            // module note.
                            *he.next_mut() = Some(lsd[j].clone());
                        }
                        append_el_spine(&spine, current_he.clone(), false)?;
                        current_he = Arc::new(ElSpine::new());

                        {
                            let ohe = Arc::make_mut(&mut current_offset_he);
                            ohe.set_last_parameter(w_first);
                            ohe.set_last_point_and_tgt(p_first.clone(), t_first.clone());
                            *ohe.next_mut() = Some(lsd[j].clone());
                        }
                        append_el_spine(&spine, current_offset_he.clone(), true)?;
                        current_offset_he = Arc::new(ElSpine::new());
                    }
                    // `CurrentHE->FirstParameter(WLast);
                    // CurrentHE->SetFirstPointAndTgt(PLast, TLast);
                    // CurrentHE->ChangePrevious() = LSD.Value(j); ...`
                    // (`cxx:3182-3188`).
                    {
                        let he = Arc::make_mut(&mut current_he);
                        he.set_first_parameter(w_last);
                        he.set_first_point_and_tgt(p_last.clone(), t_last.clone());
                        *he.previous_mut() = Some(lsd[j].clone());
                    }
                    {
                        let ohe = Arc::make_mut(&mut current_offset_he);
                        ohe.set_first_parameter(w_last);
                        ohe.set_first_point_and_tgt(p_last.clone(), t_last.clone());
                        *ohe.previous_mut() = Some(lsd[j].clone());
                    }
                    ya_kpart = true;
                } else {
                    // `if (WFirst - CurrentHE->FirstParameter() > tolesp) { ...
                    // }` (`cxx:3192-3207`).
                    if w_first - current_he.first_parameter() > tolesp {
                        {
                            let he = Arc::make_mut(&mut current_he);
                            he.set_last_parameter(w_first);
                            he.set_last_point_and_tgt(p_first.clone(), t_first.clone());
                            *he.next_mut() = Some(lsd[j].clone());
                        }
                        append_el_spine(&spine, current_he.clone(), false)?;
                        current_he = Arc::new(ElSpine::new());

                        {
                            let ohe = Arc::make_mut(&mut current_offset_he);
                            ohe.set_last_parameter(w_first);
                            ohe.set_last_point_and_tgt(p_first.clone(), t_first.clone());
                            *ohe.next_mut() = Some(lsd[j].clone());
                        }
                        append_el_spine(&spine, current_offset_he.clone(), true)?;
                        current_offset_he = Arc::new(ElSpine::new());
                    }
                    // `CurrentHE->FirstParameter(WLast); ...` (`cxx:3208-3213`).
                    {
                        let he = Arc::make_mut(&mut current_he);
                        he.set_first_parameter(w_last);
                        he.set_first_point_and_tgt(p_last.clone(), t_last.clone());
                        *he.previous_mut() = Some(lsd[j].clone());
                    }
                    {
                        let ohe = Arc::make_mut(&mut current_offset_he);
                        ohe.set_first_parameter(w_last);
                        ohe.set_first_point_and_tgt(p_last.clone(), t_last.clone());
                        *ohe.previous_mut() = Some(lsd[j].clone());
                    }
                }
            }

            // `if (!li.IsEmpty()) { myEVIMap.Bind(Spine->Edges(iedge), li); }`
            // (`cxx:3215-3218`).
            if !li.is_empty() {
                let edge = spine
                    .read()
                    .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                    .edges(iedge)
                    .cloned()
                    .ok_or_else(|| format!("PerformSetOfKPart: Edges({iedge})"))?;
                // `NCollection_DataMap::Bind` replaces an existing binding.
                let slot = self.change_evi(&edge);
                slot.clear();
                slot.extend(li);
            }
        }

        // `if (!intl || (iedgelastkpart < Spine->NbEdges())) { ... }`
        // (`cxx:3222-3265`).
        if !intl || (iedge_last_kpart < nb_edges) {
            if spine_is_periodic(&spine)? {
                // `if (WEndPeriodic - WLast > tolesp) { ... }` (`cxx:3228-3248`).
                if w_end_periodic - w_last > tolesp {
                    {
                        let he = Arc::make_mut(&mut current_he);
                        he.set_last_parameter(w_end_periodic);
                        he.set_last_point_and_tgt(
                            p_end_periodic_t_end_periodic.0.clone(),
                            p_end_periodic_t_end_periodic.1.clone(),
                        );
                        // `if (!YaKPart) { CurrentHE->SetPeriodic(true); }`
                        if !ya_kpart {
                            he.set_periodic(true);
                        }
                    }
                    append_el_spine(&spine, current_he.clone(), false)?;
                    {
                        let ohe = Arc::make_mut(&mut current_offset_he);
                        ohe.set_last_parameter(w_end_periodic);
                        ohe.set_last_point_and_tgt(
                            p_end_periodic_t_end_periodic.0.clone(),
                            p_end_periodic_t_end_periodic.1.clone(),
                        );
                        if !ya_kpart {
                            ohe.set_periodic(true);
                        }
                    }
                    append_el_spine(&spine, current_offset_he.clone(), true)?;
                }
            } else {
                // `Spine->D1(Spine->LastParameter(), PLast, TLast);
                // Spine->SetLastTgt(std::max(Spine->LastParameter(Spine->
                // NbEdges()), WLast)); ...` (`cxx:3250-3263`).
                let (last_par, d1) = {
                    let mut guard = spine
                        .write()
                        .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?;
                    let last_par = guard.last_parameter().ok_or_else(|| {
                        "PerformSetOfKPart: ChFiDS_Spine::LastParameter".to_string()
                    })?;
                    let d1 = guard.d1(last_par);
                    (last_par, d1)
                };
                p_last_t_last = d1;
                let (p_last, t_last) = p_last_t_last.clone().ok_or_else(|| {
                    "PerformSetOfKPart: ChFiDS_Spine::D1 on the last parameter".to_string()
                })?;
                {
                    let mut guard = spine
                        .write()
                        .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?;
                    let w = guard
                        .last_parameter_of(nb_edges as i32)
                        .ok_or_else(|| {
                            "PerformSetOfKPart: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
                        })?
                        .max(w_last);
                    guard.set_last_tgt(w).ok_or_else(|| {
                        "PerformSetOfKPart: SetLastTgt on a periodic spine".to_string()
                    })?;
                }
                if last_par - w_last > tolesp {
                    {
                        let he = Arc::make_mut(&mut current_he);
                        he.set_last_parameter(last_par);
                        he.set_last_point_and_tgt(p_last.clone(), t_last.clone());
                    }
                    append_el_spine(&spine, current_he.clone(), false)?;
                    {
                        let ohe = Arc::make_mut(&mut current_offset_he);
                        ohe.set_last_parameter(last_par);
                        ohe.set_last_point_and_tgt(p_last, t_last);
                    }
                    append_el_spine(&spine, current_offset_he.clone(), true)?;
                }
            }
        }

        // `NCollection_List<handle(ChFiDS_ElSpine)>& ll =
        // Spine->ChangeElSpines(); for (...) ChFi3d_PerformElSpine(...);`
        // (`cxx:3267-3270`). The list is walked by index so that the write lock
        // is not held across `ChFi3d_PerformElSpine`, which reads the spine.
        let continuity = self.my_conti;
        let el_spine_count = spine
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .el_spines()
            .len();
        for i in 0..el_spine_count {
            let mut es = spine
                .read()
                .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                .el_spines()[i]
                .clone();
            ch_fi3d_perform_el_spine(&mut es, &spine, continuity, tolesp, false)?;
            spine
                .write()
                .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                .el_spines_mut()[i] = es;
        }

        // `if (Spine->Mode() == ChFiDS_ConstThroatWithPenetrationChamfer) { ...
        // ChFi3d_PerformElSpine(ILES.ChangeValue(), Spine, myConti, tolesp, true);
        // }` (`cxx:3271-3278`).
        let mode = spine
            .read()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .mode();
        if mode == ChamfMode::ConstThroatWithPenetrationChamfer {
            let offset_count = spine
                .read()
                .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                .offset_el_spines()
                .len();
            for i in 0..offset_count {
                let mut es = spine
                    .read()
                    .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                    .offset_el_spines()[i]
                    .clone();
                ch_fi3d_perform_el_spine(&mut es, &spine, continuity, tolesp, true)?;
                spine
                    .write()
                    .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
                    .offset_el_spines_mut()[i] = es;
            }
        }

        // `Spine->SplitDone(true);` (`cxx:3279`).
        spine
            .write()
            .map_err(|_| "PerformSetOfKPart: poisoned spine".to_string())?
            .set_split_done(true);

        Ok(())
    }
}
