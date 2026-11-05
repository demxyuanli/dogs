//! `PerformIntersectionAtEnd`'s shrink block and its tail —
//! `ChFi3d_Builder_C1.cxx:3534-3572` and `:3611-3766`.
//!
//! - [`shrink`] is the `if (checkShrink && IsShrink(Ps, p1, p2, ...))` block:
//!   the p-curve of `Cc` on the previous `SurfData`'s surface clamped to
//!   `prevSDParam`, its `SurfaceCurveInterference` on `IsurfPrev`, and the
//!   `midP2d` / `midIpoint` bookkeeping of the shrunk end,
//! - [`tail`] is the block after `store_face`: `stripe->InDS(isfirst)` in the
//!   plain case, else the `ChFi3d_ComputeArete` rebuild of the curve at `Fd`'s
//!   end together with its two `SurfaceCurveInterference`s, the
//!   `stripe->InDS(isfirst, 2)` / `SetLineIndex(0)` / `CommonPoint` carry and
//!   the `Geom2dInt_GInter` self-intersection tolerance of `midIpoint`.
//!
//! Translation notes:
//! - `SDprev` is `SeqFil(num1)` (`cxx:2023-2028`), so its edits go through the
//!   sequence slot `num1 - 1` and `st.sdprev` is re-read afterwards; the
//!   `SurfData` is behind an `Arc`, and mutating through the sequence while a
//!   second handle (`st.sdprev`, `fd_of`) exists would otherwise leave that
//!   handle on the pre-edit copy,
//! - `Crv.Tolerance(std::max(Crv.Tolerance(), aTolreached))` (`cxx:3685`)
//!   writes the local `TopOpeBRepDS_Curve Crv` copy, not the DS curve
//!   `DStr.AddCurve` stored earlier (`cxx:3670-3671`), so the DS tolerance is
//!   unchanged. The port keeps that behaviour: no `change_curve` call,
//! - `Precision::Parametric(tolreached)` is `tolreached * 0.01`
//!   (`Precision.hxx:328`), spelled inline because
//!   [`occt_core::precision::parametric`] is `p / t` and is unused.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_compute_pcurv, ch_fi3d_fil_curve_in_ds,
    ch_fi3d_fil_point_in_ds, ComputeArete,
};
use crate::chfi3d::builder_c1::is_shrink;
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::ChFi3dBuilder;
use super::ds_store::{coord, set_coord};
use super::face_loop::{add_point, fi_param_on_surf, FaceLoop};
use super::{fd_of, EndIntersection};

/// `SeqFil(i)` for a 1-based index, as a mutable slot.
fn seq_slot(st: &mut EndIntersection, one_based: i64) -> Option<&mut Arc<crate::chfi_ds::SurfData>> {
    let i = usize::try_from(one_based).ok()?.checked_sub(1)?;
    st.seq.get_mut(i)
}

/// `cxx:3534-3572`: the shrink block. Answers whether it ran.
pub(super) fn shrink(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &mut FaceLoop,
    nb: i32,
) -> Result<bool, String> {
    let ps = match lp.ps.clone() {
        Some(ps) if st.check_shrink => ps,
        _ => return Ok(false),
    };
    // `IsShrink(Ps, p1, p2, checkShrParam, isUShrink,
    //  Precision::Parametric(tolreached))` (`cxx:3534-3536`).
    let tol = 0.01 * lp.tol_reached;
    if !is_shrink(
        ps.as_ref(),
        lp.p1,
        lp.p2,
        st.check_shr_param,
        st.is_u_shrink,
        tol,
    ) {
        return Ok(false);
    }
    // `shrink[nb - 1] = 1;` (`cxx:3538`).
    lp.shrink[(nb - 1) as usize] = 1;
    let sdprev = st
        .sdprev
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: SDprev is null".to_string())?;
    // `if (!isShrink) { const ChFiDS_FaceInterference& Fi =
    //  SDprev->InterferenceOnS1(); gp_Pnt2d UV =
    //  Fi.PCurveOnSurf()->Value(Fi.Parameter(isfirst));
    //  prevSDParam = isUShrink ? UV.X() : UV.Y(); }` (`cxx:3540-3547`).
    if !st.is_shrink {
        let uv = fi_param_on_surf(sdprev.interference_on_s1(), st.is_first)?;
        st.prev_sd_param = coord(&uv, if st.is_u_shrink { 1 } else { 2 });
    }
    // `gp_Pnt2d UV1 = p2d1, UV2 = p2d2; UV1.SetCoord(isUShrink ? 1 : 2,
    //  prevSDParam); UV2.SetCoord(...);` (`cxx:3548-3550`).
    let mut uv1 = lp.p2d1;
    let mut uv2 = lp.p2d2;
    let k = if st.is_u_shrink { 1 } else { 2 };
    set_coord(&mut uv1, k, st.prev_sd_param);
    set_coord(&mut uv2, k, st.prev_sd_param);
    let cc = lp
        .cc
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Cc is null".to_string())?;
    // `DStr.Surface(SDprev->Surf()).Surface()`.
    let surf = {
        let ds = b.my_ds.change_ds();
        ds.surface(sdprev.surf()).surface().cloned()
    }
    .ok_or_else(|| "IntersectionAtEnd: SDprev->Surf() has no surface".to_string())?;
    // `ChFi3d_ComputePCurv(Cc, UV1, UV2, Ps, surf, p1, p2, tolapp3d,
    //  aTolreached)` (`cxx:3551-3559`).
    let mut ps_new = Some(ps);
    let a_tol = ch_fi3d_compute_pcurv(
        cc.as_ref(),
        &uv1,
        &uv2,
        &mut ps_new,
        surf.as_ref(),
        lp.p1,
        lp.p2,
        b.tolapp3d,
        false,
    )?;
    lp.ps = ps_new;
    // `TopOpeBRepDS_Curve& TCurv = DStr.ChangeCurve(indcurve[nb - 1]);
    //  TCurv.Tolerance(std::max(TCurv.Tolerance(), aTolreached));`
    // (`cxx:3560-3561`).
    let indcurv = lp.ind_curve[(nb - 1) as usize];
    {
        let ds = b.my_ds.change_ds();
        let curv = ds.change_curve(indcurv);
        let t = curv.tolerance();
        curv.set_tolerance(t.max(a_tol));
    }
    // `InterfPS[nb - 1] = ChFi3d_FilCurveInDS(indcurve[nb - 1], IsurfPrev, Ps,
    //  orcourbe); DStr.ChangeSurfaceInterferences(IsurfPrev).Append(...)`
    // (`cxx:3563-3564`).
    let interf = ch_fi3d_fil_curve_in_ds(indcurv, st.isurf_prev, lp.ps.clone(), st.or_courbe);
    {
        let ds = b.my_ds.change_ds();
        ds.change_surface_interferences(st.isurf_prev).push(interf.clone());
    }
    lp.interf_ps[(nb - 1) as usize] = Some(interf);
    // `if (isOnSame2) { midP2d = p2d2; midIpoint = indpoint2; } else if
    //  (!isShrink) { midP2d = p2d1; midIpoint = indpoint1; } isShrink =
    //  true;` (`cxx:3566-3572`).
    if st.is_on_same2 {
        st.mid_p2d = lp.p2d2;
        st.mid_ipoint = lp.ind_point2;
    } else if !st.is_shrink {
        st.mid_p2d = lp.p2d1;
        st.mid_ipoint = lp.ind_point1;
    }
    st.is_shrink = true;
    Ok(true)
}

/// `cxx:3611-3766`: the `if (!isShrink) stripe->InDS(isfirst); else ...` tail.
pub(super) fn tail(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &mut FaceLoop,
) -> Result<(), String> {
    let is_first = st.is_first;
    if !st.is_shrink {
        // `stripe->InDS(isfirst);` (`cxx:3614`).
        if let Ok(mut stripe) = st.stripe.write() {
            stripe.in_ds(is_first, 1);
        }
        return Ok(());
    }
    // `gp_Pnt2d UV, UV1 = midP2d, UV2 = midP2d; if (isOnSame1) { UV = UV2 =
    //  Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst)); } else { UV = UV1 =
    //  Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst)); }`
    // (`cxx:3622-3630`).
    let mut uv1 = st.mid_p2d;
    let mut uv2 = st.mid_p2d;
    if st.is_on_same1 {
        uv2 = fi_param_on_surf(&st.fi1, !is_first)?;
    } else {
        uv1 = fi_param_on_surf(&st.fi2, !is_first)?;
    }
    // `aSurf = DStr.Surface(Fd->Surf()).Surface();` (`cxx:3633`).
    let a_surf = {
        let ds = b.my_ds.change_ds();
        ds.surface(fd_of(st)?.surf()).surface().cloned()
    }
    .ok_or_else(|| "IntersectionAtEnd: Fd->Surf() has no surface".to_string())?;
    let (cv1, cv2) = super::end_points(st)?;
    // `ChFi3d_ComputeArete(CV1, UV1, CV2, UV2, aSurf, C3d, Ps, p1, p2,
    //  tolapp3d, tol2d, aTolreached, 0)` (`cxx:3636-3648`).
    let mut arete = ComputeArete {
        pardeb: lp.p1,
        parfin: lp.p2,
        ..Default::default()
    };
    ch_fi3d_compute_arete(
        &cv1,
        &uv1,
        &cv2,
        &uv2,
        &a_surf,
        &mut arete,
        b.tolapp3d,
        b.tol2d,
        0,
    )?;
    let c3d = arete
        .c3d
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: ChFi3d_ComputeArete returned no C3d".to_string())?;
    lp.ps = arete.p_curve.clone();
    lp.p1 = arete.pardeb;
    lp.p2 = arete.parfin;
    let p1 = lp.p1;
    let p2 = lp.p2;
    let a_tol = arete.tolreached;
    // `indpoint1 = indpoint2 = midIpoint;` then the `C3d->Value(...)` point of
    // the side `isOnSame1` picks (`cxx:3650-3667`).
    lp.ind_point1 = st.mid_ipoint;
    lp.ind_point2 = st.mid_ipoint;
    if st.is_on_same1 {
        lp.ind_point2 = add_point(b, c3d.d0(p2), a_tol);
    } else {
        lp.ind_point1 = add_point(b, c3d.d0(p1), a_tol);
    }
    // `TopOpeBRepDS_Curve Crv(C3d, aTolreached); int Icurv =
    //  DStr.AddCurve(Crv);` (`cxx:3669-3671`).
    let icurv = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&c3d)),
        a_tol,
        false,
    ));
    // `Interfp1` / `Interfp2` and `Interfc = ChFi3d_FilCurveInDS(Icurv, Isurf,
    //  Ps, orcourbe)` (`cxx:3672-3677`).
    let interfp1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, icurv, lp.ind_point1, p1, false);
    let interfp2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, icurv, lp.ind_point2, p2, false);
    let interfc = ch_fi3d_fil_curve_in_ds(icurv, st.isurf, lp.ps.clone(), st.or_courbe);
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(icurv);
        list.push(interfp1);
        list.push(interfp2);
        ds.change_surface_interferences(st.isurf).push(interfc);
    }
    // `aSurf = DStr.Surface(SDprev->Surf()).Surface(); UV1.SetCoord(...);
    //  UV2.SetCoord(...); ChFi3d_ComputePCurv(C3d, UV1, UV2, Pc, aSurf, p1, p2,
    //  tolapp3d, aTolreached);` (`cxx:3680-3683`).
    let sdprev = st
        .sdprev
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: SDprev is null".to_string())?;
    let prev_surf = {
        let ds = b.my_ds.change_ds();
        ds.surface(sdprev.surf()).surface().cloned()
    }
    .ok_or_else(|| "IntersectionAtEnd: SDprev->Surf() has no surface".to_string())?;
    let k = if st.is_u_shrink { 1 } else { 2 };
    set_coord(&mut uv1, k, st.prev_sd_param);
    set_coord(&mut uv2, k, st.prev_sd_param);
    let mut pc_new = lp.pc.clone();
    let a_tol2 = ch_fi3d_compute_pcurv(
        c3d.as_ref(),
        &uv1,
        &uv2,
        &mut pc_new,
        prev_surf.as_ref(),
        p1,
        p2,
        b.tolapp3d,
        false,
    )?;
    lp.pc = pc_new;
    // `Crv.Tolerance(std::max(Crv.Tolerance(), aTolreached));` (`cxx:3685`) —
    // a write to the local `Crv` copy, not to the DS curve; see the module
    // note.
    let _ = a_tol2;
    // `Interfc = ChFi3d_FilCurveInDS(Icurv, IsurfPrev, Pc,
    //  TopAbs::Reverse(orcourbe)); DStr.ChangeSurfaceInterferences(
    //  IsurfPrev).Append(Interfc);` (`cxx:3686-3687`).
    let interfc2 = ch_fi3d_fil_curve_in_ds(
        icurv,
        st.isurf_prev,
        lp.pc.clone(),
        st.or_courbe.reversed(),
    );
    b.my_ds
        .change_ds()
        .change_surface_interferences(st.isurf_prev)
        .push(interfc2);
    // `stripe->InDS(isfirst, 2);` (`cxx:3694`).
    if let Ok(mut stripe) = st.stripe.write() {
        stripe.in_ds(is_first, 2);
    }
    // `Fd->ChangeInterference(isOnSame1 ? 2 : 1).SetLineIndex(0);`
    // (`cxx:3695`).
    let side = if st.is_on_same1 { 2 } else { 1 };
    let num = st.num as i64;
    if let Some(slot) = seq_slot(st, num) {
        Arc::make_mut(slot)
            .interference_mut(side)
            .set_line_index(0);
    }
    // `ChFiDS_CommonPoint& CPprev1 = SDprev->ChangeVertex(isfirst, side);
    //  ... CPlast1 = Fd->ChangeVertex(isfirst, side); CPlast2 =
    //  Fd->ChangeVertex(!isfirst, side); if (CPprev1.IsOnArc()) { CPlast1 =
    //  CPprev1; CPprev1.Reset(); CPprev1.SetPoint(CPlast1.Point());
    //  CPlast2.Reset(); CPlast2.SetPoint(CPlast1.Point()); }`
    // (`cxx:3696-3708`).
    let cp_prev1 = sdprev.vertex(is_first, side).clone();
    if cp_prev1.is_on_arc() {
        let point = *cp_prev1.point();
        if let Some(slot) = seq_slot(st, num) {
            let fd = Arc::make_mut(slot);
            *fd.vertex_mut(is_first, side) = cp_prev1.clone();
            let last = fd.vertex_mut(!is_first, side);
            last.reset();
            last.set_point(&point);
        }
        if let Some(slot) = seq_slot(st, st.num1) {
            let prev = Arc::make_mut(slot);
            let cp = prev.vertex_mut(is_first, side);
            cp.reset();
            cp.set_point(&point);
        }
        st.sdprev = seq_slot(st, st.num1).cloned();
    }
    // `for (nb = 0; nb < nbface; nb++) if (isOnSame1 ? shrink[nb + 1] :
    //  !shrink[nb]) break;` (`cxx:3715-3719`).
    let nbface = st.nbface as usize;
    let mut nb_found = nbface;
    for i in 0..nbface {
        let hit = if st.is_on_same1 {
            lp.shrink.get(i + 1).copied().unwrap_or(0) != 0
        } else {
            lp.shrink.get(i).copied().unwrap_or(0) == 0
        };
        if hit {
            nb_found = i;
            break;
        }
    }
    // `Cend = DStr.Curve(indcurve[nb]).Curve(); PCend = InterfPS[nb]->PCurve();`
    // (`cxx:3720-3722`).
    let ind_end = *lp
        .ind_curve
        .get(nb_found)
        .ok_or_else(|| "IntersectionAtEnd: indcurve[nb] is out of range".to_string())?;
    let cend = {
        let ds = b.my_ds.change_ds();
        ds.curve(ind_end).curve().cloned()
    }
    .ok_or_else(|| "IntersectionAtEnd: DStr.Curve(indcurve[nb]) is null".to_string())?;
    let pcend = lp
        .interf_ps
        .get(nb_found)
        .and_then(|i| i.as_ref())
        .and_then(|i| i.read().ok().and_then(|g| g.p_curve().cloned()))
        .ok_or_else(|| "IntersectionAtEnd: InterfPS[nb] has no p-curve".to_string())?;
    // `TopOpeBRepDS_Point& Pds = DStr.ChangePoint(midIpoint); const gp_Pnt&
    //  Pvert = Pds.Point(); double tol = Pds.Tolerance();`
    // (`cxx:3724-3727`).
    let (pvert, mut tol): (GpPnt, f64) = {
        let ds = b.my_ds.change_ds();
        let pds = ds
            .change_point(st.mid_ipoint)
            .ok_or_else(|| "IntersectionAtEnd: midIpoint is not in the DS".to_string())?;
        (*pds.change_point(), pds.tolerance())
    };
    // `Geom2dAdaptor_Curve PC1(Ps), PC2(PCend); Geom2dInt_GInter Intersector(
    //  PC1, PC2, PConfusion, PConfusion);` (`cxx:3729-3731`).
    let ps = lp
        .ps
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Ps is null".to_string())?;
    let mut intersector = Geom2dIntGInter::new();
    intersector.perform_curves(ps.as_ref(), pcend.as_ref(), PCONFUSION, PCONFUSION);
    // `if (!Intersector.IsDone()) return;` (`cxx:3732-3735`).
    if !intersector.is_done() {
        return Ok(());
    }
    // `for (nb = 1; nb <= Intersector.NbPoints(); nb++) { ... tol =
    //  max(tol, Pvert.Distance(C3d->Value(ip.ParamOnFirst())));
    //  tol = max(tol, Pvert.Distance(Cend->Value(ip.ParamOnSecond()))); }`
    // (`cxx:3736-3742`).
    for i in 1..=intersector.nb_points() {
        let ip = intersector.point(i);
        tol = tol.max(pvert.distance(&c3d.d0(ip.param_on_first())));
        tol = tol.max(pvert.distance(&cend.d0(ip.param_on_second())));
    }
    // `for (nb = 1; nb <= Intersector.NbSegments(); nb++) { ... }`
    // (`cxx:3743-3761`).
    for i in 1..=intersector.nb_segments() {
        let seg = intersector.segment(i);
        let mut acc = |ip: &occt_core::intres2d::IntRes2dIntersectionPoint| {
            tol = tol.max(pvert.distance(&c3d.d0(ip.param_on_first())));
            tol = tol.max(pvert.distance(&cend.d0(ip.param_on_second())));
        };
        if seg.has_first_point() {
            acc(seg.first_point());
        }
        if seg.has_last_point() {
            acc(seg.last_point());
        }
    }
    // `Pds.Tolerance(tol);` (`cxx:3762`).
    let ds = b.my_ds.change_ds();
    if let Some(pds) = ds.change_point(st.mid_ipoint) {
        pds.set_tolerance(tol);
    }
    Ok(())
}
