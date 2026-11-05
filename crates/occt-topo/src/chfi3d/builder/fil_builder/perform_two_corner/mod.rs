//! `ChFi3d_FilBuilder::PerformTwoCorner` - the corner shared by two fillet
//! stripes of a vertex.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder_C2.cxx:143-1153`
//! (OCCT 8.0.0). The body is split by OCCT line ranges:
//!
//! - this file - `cxx:143-514`: `ToricRotule` (`cxx:87-122`) and `RemoveSD`
//!   (`cxx:124-139`) of the same file, the two stripe ends and their spine
//!   edges, the `ang1` tangency test of `cxx:230-236`, the `ChFi3d_IsInFront`
//!   call of `cxx:239-257` with its `nbsurf1` / `nbsurf2` retry loop
//!   (`cxx:259-322`) and the `RemoveSD` trims of `cxx:329-345`, the two spine
//!   statuses that pick `c1biseau` / `c1rotule` (`cxx:360-368`), the `pivot`
//!   arc with its `sameparam` test (`cxx:376-405`) and the `FF1` / `FF2` /
//!   `myEFMap(pivot)` checks (`cxx:418-451`), then the bevel and knee
//!   `PerformTwoCornerbyInter` calls (`cxx:459-514`),
//! - [`fill`] - `cxx:516-1152`: the corner `SurfData` built by hand when the
//!   two stripes are neither bevel nor knee.
//!
//! Translation notes:
//! - OCCT's `done` is [`ChFi3dBuilder::done`](super::super::ChFi3dBuilder):
//!   every `return` of the body sets it the way the C++ does,
//! - `st1` / `st2` are `ChFiDS_Stripe` handles, so `st2 == st1` of `cxx:184` is
//!   [`std::sync::Arc::ptr_eq`],
//! - `sd1` / `sd2` are **copies** of the stripes' `ChFiDS_SurfData` (OCCT
//!   holds references into the sequences); [`CornerCtx::flush`] writes them
//!   back, which is also what keeps the two writes coherent when the two
//!   stripes are the same handle,
//! - `UIntPC1` / `UIntPC2` are in/out in OCCT and are never reset between the
//!   `ChFi3d_IsInFront` calls of the retry loop; [`IsInFrontState::p1`] /
//!   `p2` carry them, and each call is seeded with the previous answer,
//! - a null handle OCCT dereferences (a missing stripe, an out-of-range
//!   `SurfData`, a common point without an arc, a null `DStr.Shape`) is the
//!   `Err` of this port, and so is `LProp_NotDefined` from
//!   `BRepLProp_CLProps::Tangent` (`cxx:220-221`).

mod fill;

use std::sync::Arc;

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_index_of_surf_data, ch_fi3d_is_in_front, IsInFrontState,
};
use crate::chfi_ds::{State, Stripe, StripeHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::lprop::CLProps;

use super::super::perform_more_three_corner::spine_edge_at;
use super::super::ChFi3dBuilder;
use fill::CornerCtx;

/// `Corner->SetOfSurfData()->Value(index)` with the 1-based `index`.
fn surf_data(stripe: &StripeHandle, index: i32) -> Result<SurfData, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
    index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .and_then(|i| guard.set_of_surf_data().get(i))
        .map(|sd| sd.as_ref().clone())
        .ok_or_else(|| format!("PerformTwoCorner: SurfData {index} is out of range"))
}

/// `Corner->SetOfSurfData()->Length()`.
fn surf_data_len(stripe: &StripeHandle) -> Result<i32, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
    i32::try_from(guard.set_of_surf_data().len())
        .map_err(|_| "PerformTwoCorner: the SurfData sequence is too long".to_string())
}

/// `Corner->Spine()->Status(isfirst)` (`cxx:366-367`).
fn spine_status(stripe: &StripeHandle, is_first: bool) -> Result<State, String> {
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "PerformTwoCorner: the stripe has no spine".to_string())?
    };
    let guard = spine
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned spine".to_string())?;
    Ok(guard.status(is_first))
}

/// `RemoveSD(Stripe, num1, num2)` (`cxx:124-139`).
///
/// OCCT's empty-sequence return and its `num1 == num2` single removal are
/// both the inclusive range `[num1, num2]`.
fn remove_sd(stripe: &mut Stripe, num1: i32, num2: i32) {
    let seq = stripe.change_set_of_surf_data();
    if seq.is_empty() || num1 < 1 || num2 < num1 {
        return;
    }
    let start = (num1 - 1) as usize;
    let end = (num2 as usize).min(seq.len());
    if start >= end {
        return;
    }
    seq.drain(start..end);
}

/// One `ChFi3d_IsInFront` call of `cxx:239-257` / `cxx:300-317`.
///
/// `p1` / `p2` are in/out in OCCT: the retry loop of `cxx:297-321` seeds each
/// call with what the previous one answered.
#[allow(clippy::too_many_arguments)]
fn is_in_front(
    b: &ChFi3dBuilder,
    st1: &StripeHandle,
    st2: &StripeHandle,
    isd1: i32,
    isd2: i32,
    sens1: i32,
    sens2: i32,
    vtx: &TopoShape,
    p1: f64,
    p2: f64,
) -> Result<IsInFrontState, String> {
    let cd1 = st1
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
    let cd2 = st2
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
    let mut st = IsInFrontState::new(p1, p2);
    ch_fi3d_is_in_front(
        &b.my_ds, &cd1, &cd2, isd1, isd2, sens1, sens2, &mut st, Some(vtx), true, false,
    )?;
    Ok(st)
}

/// `PerformTwoCorner(Index)` (`ChFi3d_FilBuilder_C2.cxx:143-1153`).
///
/// `shape` is the `myShape` of `ChFi3d_FilBuilder` (`hxx:377`), the only
/// subclass member the body reads (through the two `CompleteData` calls of
/// [`fill`]).
pub(super) fn perform_two_corner(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    index: usize,
) -> Result<(), String> {
    // `done = false; const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);`
    // (`cxx:150-151`).
    base.done = false;
    let vtx = base
        .my_v_data_map
        .find_key(index)
        .cloned()
        .ok_or_else(|| format!("PerformTwoCorner: myVDataMap has no key {index}"))?;
    // `It.Initialize(myVDataMap(Index)); ... st1 = It.Value();`
    // (`cxx:154-176`).
    let stripes: Vec<StripeHandle> = base
        .my_v_data_map
        .get_from_index(index)
        .map(|s| s.to_vec())
        .ok_or_else(|| format!("PerformTwoCorner: myVDataMap({index}) is out of range"))?;
    let st1 = stripes
        .first()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: the vertex has no stripe".to_string())?;
    let mut sens1 = 0;
    // `Isd1 = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);` (`cxx:177`).
    let mut isd1 = ch_fi3d_index_of_surf_data(&vtx, &st1, &mut sens1)?;

    // `It.Next(); st2 = It.Value(); if (st2 == st1) { Sens2 = -1; Isd2 =
    // st2->SetOfSurfData()->Length(); } else { Isd2 =
    // ChFi3d_IndexOfSurfData(Vtx, st2, Sens2); }` (`cxx:181-192`).
    let st2 = stripes
        .get(1)
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: the vertex has a single stripe".to_string())?;
    let (mut sens2, mut isd2) = if Arc::ptr_eq(&st1, &st2) {
        (-1, surf_data_len(&st2)?)
    } else {
        let mut s = 0;
        let i = ch_fi3d_index_of_surf_data(&vtx, &st2, &mut s)?;
        (s, i)
    };

    // `if (Sens1 == 1) E1 = st1->Spine()->Edges(1); else E1 =
    // st1->Spine()->Edges(st1->Spine()->NbEdges());` and the same for `E2`
    // (`cxx:196-212`).
    let e1 = spine_edge_at(&st1, sens1 == 1)?;
    let e2 = spine_edge_at(&st2, sens2 == 1)?;
    // `BRepAdaptor_Curve BCurv1(E1); ... parE1 = BRep_Tool::Parameter(Vtx,
    // E1);` (`cxx:214-217`).
    let bc1 = CurveAdaptor::of_edge(&e1);
    let bc2 = CurveAdaptor::of_edge(&e2);
    let par_e1 = crate::brep_tool::BRepTool::parameter_of_vertex(&e1.0, &vtx);
    let par_e2 = crate::brep_tool::BRepTool::parameter_of_vertex(&e2.0, &vtx);
    // `BRepLProp_CLProps CL1(BCurv1, parE1, 1, 1.e-4); CL1.Tangent(dir1);
    // ... if (Sens1 == -1) dir1.Reverse(); if (Sens2 == -1)
    // dir2.Reverse();` (`cxx:218-229`).
    let mut cl1 = CLProps::new(&bc1, par_e1, 1, 1.0e-4)?;
    let mut cl2 = CLProps::new(&bc2, par_e2, 1, 1.0e-4)?;
    let mut dir1 = cl1.tangent(&bc1, par_e1, &mut 0)?;
    let mut dir2 = cl2.tangent(&bc2, par_e2, &mut 0)?;
    if sens1 == -1 {
        dir1 = dir1.reversed();
    }
    if sens2 == -1 {
        dir2 = dir2.reversed();
    }
    // `ang1 = std::abs(dir1.Angle(dir2)); if (ang1 < M_PI / 180.) {
    // PerformMoreThreeCorner(Index, 2); done = true; return; }`
    // (`cxx:230-236`).
    let ang1 = dir1.angle(&dir2).abs();
    if ang1 < std::f64::consts::PI / 180.0 {
        base.perform_more_three_corner(index, 2)?;
        base.done = true;
        return Ok(());
    }

    // `OkinterCC = ChFi3d_IsInFront(DStr, st1, st2, Isd1, Isd2, Sens1, Sens2,
    // UIntPC1, UIntPC2, FaCo, SameSide, IFaCo1, IFaCo2, Okvisavis, Vtx,
    // true);` (`cxx:239-257`). `UIntPC1` / `UIntPC2` are uninitialized
    // out-parameters in OCCT; the port seeds them with `0.0`.
    let mut st = is_in_front(base, &st1, &st2, isd1, isd2, sens1, sens2, &vtx, 0.0, 0.0)?;
    let mut okinter_cc = st.ok;
    if !st.visavis {
        // `nbsurf1 = st1->SetOfSurfData()->Length(); ... deb1 = 1; ...`
        // (`cxx:261-292`).
        let nbsurf1 = surf_data_len(&st1)?;
        let nbsurf2 = surf_data_len(&st2)?;
        let mut deb1 = 1;
        let mut deb2 = 1;
        let mut fin1 = 1;
        let mut fin2 = 1;
        if nbsurf1 != 1 {
            if sens1 == 1 {
                deb1 = 1;
                fin1 = 2;
            } else {
                deb1 = nbsurf1 - 1;
                fin1 = nbsurf1;
            }
        }
        if nbsurf2 != 1 {
            if sens2 == 1 {
                deb2 = 1;
                fin2 = 2;
            } else {
                deb2 = nbsurf2 - 1;
                fin2 = nbsurf2;
            }
        }
        // `for (i1 = deb1; i1 <= fin1 && !trouve; i1++) { Isd1 = i1; for (i2 =
        // deb2; i2 <= fin2 && !trouve; i2++) { Isd2 = i2; OkinterCC =
        // ChFi3d_IsInFront(...); trouve = Okvisavis; } }` (`cxx:295-321`).
        let mut trouve = false;
        'outer: for i1 in deb1..=fin1 {
            isd1 = i1;
            for i2 in deb2..=fin2 {
                isd2 = i2;
                st = is_in_front(
                    base, &st1, &st2, isd1, isd2, sens1, sens2, &vtx, st.p1, st.p2,
                )?;
                okinter_cc = st.ok;
                if st.visavis {
                    trouve = true;
                    break 'outer;
                }
            }
        }
        // `if (!trouve) { PerformMoreThreeCorner(Index, 2); done = true;
        // return; }` (`cxx:322-327`).
        if !trouve {
            base.perform_more_three_corner(index, 2)?;
            base.done = true;
            return Ok(());
        }
        // `if (Sens1 == 1 && Isd1 != 1) RemoveSD(st1, 1, 1); if (Sens1 != 1 &&
        // Isd1 != nbsurf1) RemoveSD(st1, fin1, fin1); ...` (`cxx:329-345`).
        if sens1 == 1 && isd1 != 1 {
            remove_sd_of(&st1, 1, 1)?;
        }
        if sens1 != 1 && isd1 != nbsurf1 {
            remove_sd_of(&st1, fin1, fin1)?;
        }
        if sens2 == 1 && isd2 != 1 {
            remove_sd_of(&st2, 1, 1)?;
        }
        if sens2 != 1 && isd2 != nbsurf2 {
            remove_sd_of(&st2, fin2, fin2)?;
        }
        // `Isd1 = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1); Isd2 =
        // ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);` (`cxx:346-347`).
        isd1 = ch_fi3d_index_of_surf_data(&vtx, &st1, &mut sens1)?;
        isd2 = ch_fi3d_index_of_surf_data(&vtx, &st2, &mut sens2)?;
    }

    // `int IFaArc1 = 3 - IFaCo1, IFaArc2 = 3 - IFaCo2;` (`cxx:350`).
    let i_fa_co1 = st.jf1;
    let i_fa_co2 = st.jf2;
    let i_fa_arc1 = 3 - i_fa_co1;
    let i_fa_arc2 = 3 - i_fa_co2;
    // `sd1 = SeqFil1.ChangeValue(Isd1); ... sd2 = SeqFil2.ChangeValue(Isd2);`
    // (`cxx:355-357`).
    let sd1 = surf_data(&st1, isd1)?;
    let sd2 = surf_data(&st2, isd2)?;
    // `TopAbs_Orientation OFaCo = FaCo.Orientation();` (`cxx:358`).
    let fa_co = st
        .face
        .clone()
        .ok_or_else(|| "PerformTwoCorner: ChFi3d_IsInFront found no common face".to_string())?;
    let o_fa_co = fa_co.0.orientation();

    // `isfirst1 = (Sens1 == 1); isfirst2 = (Sens2 == 1); Stat1 =
    // st1->Spine()->Status(isfirst1); Stat2 = ...; bool c1biseau = (Stat1 ==
    // ChFiDS_AllSame); bool c1rotule = (Stat1 == ChFiDS_OnSame && Stat2 ==
    // ChFiDS_OnSame);` (`cxx:362-368`).
    let is_first1 = sens1 == 1;
    let is_first2 = sens2 == 1;
    let stat1 = spine_status(&st1, is_first1)?;
    let stat2 = spine_status(&st2, is_first2)?;
    let c1biseau = stat1 == State::AllSame;
    let c1rotule = stat1 == State::OnSame && stat2 == State::OnSame;

    // `ChFiDS_CommonPoint& CP1 = sd1->ChangeVertex(isfirst1, IFaArc1);
    // ChFiDS_CommonPoint& CP2 = sd2->ChangeVertex(isfirst2, IFaArc2);`
    // (`cxx:372-373`).
    let cp1 = sd1.vertex(is_first1, i_fa_arc1).clone();
    let cp2 = sd2.vertex(is_first2, i_fa_arc2).clone();

    // `TopoDS_Edge pivot; bool yapiv = false; if (CP1.IsOnArc()) { pivot =
    // CP1.Arc(); } else { PerformMoreThreeCorner(Index, 2); done = true;
    // return; }` (`cxx:378-388`).
    let pivot: TopoShape = match cp1.arc() {
        Some(arc) => arc.clone(),
        None => {
            base.perform_more_three_corner(index, 2)?;
            base.done = true;
            return Ok(());
        }
    };
    // `if (CP1.IsOnArc() && CP2.IsOnArc()) yapiv =
    // (pivot.IsSame(CP2.Arc()));` (`cxx:390-393`).
    let yapiv = match cp2.arc() {
        Some(arc) => is_same(&pivot, arc),
        None => false,
    };
    // `occ::handle<BRepAdaptor_Curve> Hpivot; bool sameparam = false; double
    // parCP1 = 0., parCP2 = 0.; if (yapiv) { Hpivot = new
    // BRepAdaptor_Curve(pivot); parCP1 = CP1.ParameterOnArc(); parCP2 =
    // CP2.ParameterOnArc(); gp_Pnt tst1 = Hpivot->Value(parCP1); gp_Pnt tst2
    // = Hpivot->Value(parCP2); sameparam = tst1.Distance(tst2) <= tolapp3d;
    // }` (`cxx:394-406`).
    let mut par_cp1 = 0.0;
    let mut par_cp2 = 0.0;
    let mut same_param = false;
    let hpivot = if yapiv {
        let hp = CurveAdaptor::of_edge(&Edge(pivot.clone()));
        par_cp1 = cp1
            .parameter_on_arc()
            .ok_or_else(|| "PerformTwoCorner: CP1 has no parameter on arc".to_string())?;
        par_cp2 = cp2
            .parameter_on_arc()
            .ok_or_else(|| "PerformTwoCorner: CP2 has no parameter on arc".to_string())?;
        let tst1 = hp
            .d0(par_cp1)
            .ok_or_else(|| "PerformTwoCorner: the pivot has no value".to_string())?;
        let tst2 = hp
            .d0(par_cp2)
            .ok_or_else(|| "PerformTwoCorner: the pivot has no value".to_string())?;
        same_param = tst1.distance(&tst2) <= base.tolapp3d;
        Some(hp)
    } else {
        None
    };

    // `FF1 = TopoDS::Face(DStr.Shape(sd1->Index(IFaArc1))); FF2 = ...; if
    // (FF1.IsNull() || FF2.IsNull()) { PerformMoreThreeCorner(Index, 2);
    // done = true; return; }` (`cxx:419-427`).
    let ff1 = base
        .my_ds
        .ds()
        .shape(sd1.index(i_fa_arc1), true)
        .and_then(Face::wrap);
    let ff2 = base
        .my_ds
        .ds()
        .shape(sd2.index(i_fa_arc2), true)
        .and_then(Face::wrap);
    let (Some(ff1), Some(ff2)) = (ff1, ff2) else {
        base.perform_more_three_corner(index, 2)?;
        base.done = true;
        return Ok(());
    };

    // `if (yapiv) { for (Kt.Initialize(myEFMap(pivot)); ...) { F =
    // TopoDS::Face(Kt.Value()); if (!ok1 && FF1.IsSame(F)) ok1 = true; if
    // (!ok2 && FF2.IsSame(F)) ok2 = true; } if (!ok1 || !ok2) {
    // PerformMoreThreeCorner(Index, 2); done = true; return; } }`
    // (`cxx:430-451`).
    if yapiv {
        let faces = base.my_ef_map.get_from_key(&pivot);
        let mut ok1 = false;
        let mut ok2 = false;
        for f in faces {
            if !ok1 && is_same(f, &ff1.0) {
                ok1 = true;
            }
            if !ok2 && is_same(f, &ff2.0) {
                ok2 = true;
            }
        }
        if !ok1 || !ok2 {
            base.perform_more_three_corner(index, 2)?;
            base.done = true;
            return Ok(());
        }
    }

    // `if (c1biseau) { done = PerformTwoCornerbyInter(Index); if (!done) {
    // PerformMoreThreeCorner(Index, 2); done = true; return; } }`
    // (`cxx:467-481`).
    if c1biseau {
        base.perform_two_corner_by_inter(index)?;
        if !base.done {
            base.perform_more_three_corner(index, 2)?;
            base.done = true;
            return Ok(());
        }
    }

    let mut ctx = CornerCtx {
        st1,
        st2,
        sd1,
        sd2,
        isd1,
        isd2,
        is_first1,
        is_first2,
        i_fa_co1,
        i_fa_co2,
        i_fa_arc1,
        i_fa_arc2,
        o_fa_co,
        fa_co,
        ff1,
        ff2,
        pivot,
        hpivot,
        yapiv,
        same_param,
        par_cp1,
        par_cp2,
        okinter_cc,
        stat1,
        stat2,
        same_side: st.sameside,
        u_int_pc1: st.p1,
        u_int_pc2: st.p2,
        cp1,
        cp2,
        reset_cp1: false,
        reset_cp2: false,
        pctrans: crate::abs::Orientation::Forward,
        fa_piv: None,
        sens1,
        sens2,
        vtx: vtx.clone(),
    };

    // `else if (c1rotule) { // save. cp11 = sd1->Vertex(isfirst1, 1); ... done
    // = PerformTwoCornerbyInter(Index); if (!done) { restore ... done =
    // false; } }` (`cxx:482-514`).
    if c1rotule {
        let cp11 = ctx.sd1.vertex(is_first1, 1).clone();
        let cp12 = ctx.sd1.vertex(is_first1, 2).clone();
        let cp21 = ctx.sd2.vertex(is_first2, 1).clone();
        let cp22 = ctx.sd2.vertex(is_first2, 2).clone();
        let intf11 = ctx.sd1.interference_on_s1().clone();
        let intf12 = ctx.sd1.interference_on_s2().clone();
        let intf21 = ctx.sd2.interference_on_s1().clone();
        let intf22 = ctx.sd2.interference_on_s2().clone();
        base.perform_two_corner_by_inter(index)?;
        if !base.done {
            // The two `SurfData` copies must first catch up with what
            // `PerformTwoCornerbyInter` wrote into the stripes, so that the
            // restore below is the four-field restore of `cxx:503-511` and not
            // a wholesale rollback.
            ctx.sd1 = surf_data(&ctx.st1, ctx.isd1)?;
            ctx.sd2 = surf_data(&ctx.st2, ctx.isd2)?;
            ctx.sd1.vertex_mut(is_first1, 1).clone_from(&cp11);
            ctx.sd1.vertex_mut(is_first1, 2).clone_from(&cp12);
            ctx.sd2.vertex_mut(is_first2, 1).clone_from(&cp21);
            ctx.sd2.vertex_mut(is_first2, 2).clone_from(&cp22);
            ctx.sd1.set_interference_on_s1(intf11);
            ctx.sd1.set_interference_on_s2(intf12);
            ctx.sd2.set_interference_on_s1(intf21);
            ctx.sd2.set_interference_on_s2(intf22);
            ctx.flush()?;
            base.done = false;
        }
    }

    if !c1biseau && !base.done {
        // `cxx:516-1152`.
        return fill::two_corner_fill(base, shape, &mut ctx);
    }
    Ok(())
}

/// `Arc::make_mut` on a slot of a `SurfData` sequence, used by
/// [`CornerCtx::flush`].
pub(super) fn index_slot(seq: &mut [Arc<SurfData>], index: i32) -> Option<&mut Arc<SurfData>> {
    index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .and_then(move |i| seq.get_mut(i))
}

/// `RemoveSD(stripe, num1, num2)` on a handle, the `Err` being a poisoned
/// stripe lock.
fn remove_sd_of(stripe: &StripeHandle, num1: i32, num2: i32) -> Result<(), String> {
    let mut guard = stripe
        .write()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
    remove_sd(&mut guard, num1, num2);
    Ok(())
}
