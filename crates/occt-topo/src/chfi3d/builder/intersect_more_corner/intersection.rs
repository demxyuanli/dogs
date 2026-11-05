//! `IntersectMoreCorner`'s intersection step — `ChFi3d_Builder_C1.cxx:5004-5123`.
//!
//! `ChFi3d_BoundSurf` bounds the fillet surface, `ChFi3d_BoundFac` the face,
//! `ChFi3d_ComputeCurves` intersects them into `Cc` with its two p-curves, and
//! the sewing edge of `Fv` is tested against `Cc` with `Extrema_ExtCC`: a hit
//! inside `Precision::SquareConfusion()` splits `Cc` at the sewing parameter
//! (`intcouture`).
//!
//! Divergences that are not heuristics: `ChFi3d_Recale` runs unconditionally
//! here (the `PerformOneCorner` body guards it with `if (onsame)`), and the
//! p-curve values are read with OCCT's `LineIndex() != 0` choice like
//! `PerformOneCorner`.

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom::extrema_cc::curve_curve_extrema_all;
use occt_geom::Curve;

use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_boite, ch_fi3d_bound_fac, ch_fi3d_bound_surf, ch_fi3d_compute_curves, ComputeCurves,
};
use crate::chfi3d::builder_c1::ch_fi3d_recale;
use crate::chfi3d::builder_helpers::couture;
use crate::shape::Edge;
use crate::topopebrep::ds::{Curve as DsCurve, Point};
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

use super::super::ChFi3dBuilder;
use super::MoreCorner;

/// `IntersectMoreCorner`'s `if (inters)` block (`cxx:5004-5123`).
pub(super) fn intersect(b: &mut ChFi3dBuilder, mc: &mut MoreCorner) -> Result<(), String> {
    let is_first = mc.is_first;
    let fv = mc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;

    // `HGs = ChFi3d_BoundSurf(DStr, Fd, 1, 2);` (`cxx:5006`).
    let hgs = ch_fi3d_bound_surf(b.my_ds.change_ds(), &mc.fd, 1, 2)?;
    // `const ChFiDS_FaceInterference& Fi1 = Fd->InterferenceOnS1();` and `Fi2`
    // (`cxx:5007-5008`).
    let fi1 = mc.fd.interference_on_s1().clone();
    let fi2 = mc.fd.interference_on_s2().clone();
    let cv1 = mc.fd.vertex(is_first, 1).clone();
    let cv2 = mc.fd.vertex(is_first, 2).clone();

    // `if (IFopArc == 1) pfac1 = p2dbout; else { Hc1 = CurveOnSurface(CV1.Arc(),
    //  Fv, ...); pfac1 = Hc1->Value(CV1.ParameterOnArc()); }` (`cxx:5012-5024`).
    let pfac1 = if mc.i_fop_arc == 1 {
        mc.p2dbout
    } else {
        corner_p2d(&fv, &cv1)?
    };
    // `cxx:5025-5037`.
    let pfac2 = if mc.i_fop_arc == 2 {
        mc.p2dbout
    } else {
        corner_p2d(&fv, &cv2)?
    };

    // `pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));`, with the
    // `!isfirst` arm when `LineIndex() == 0` (`cxx:5038-5053`).
    let pfil1 = surf_pcurve_point(&fi1, is_first)?;
    let pfil2 = surf_pcurve_point(&fi2, is_first)?;

    // `ChFi3d_Recale(Bs, pfac1, pfac2, (IFadArc == 1));` (`cxx:5054`).
    let mut pfac1 = pfac1;
    let mut pfac2 = pfac2;
    {
        let hbs = mc
            .hbs
            .as_ref()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
        ch_fi3d_recale(hbs, &mut pfac1, &mut pfac2, mc.i_fad_arc == 1);
    }

    // `Pardeb(1 .. 4)` / `Parfin(1 .. 4)` (`cxx:5055-5062`).
    let pardeb = [pfil1.x(), pfil1.y(), pfac1.x(), pfac1.y()];
    let parfin = [pfil2.x(), pfil2.y(), pfac2.x(), pfac2.y()];

    // `ChFi3d_Boite(pfac1, pfac2, ...); ChFi3d_BoundFac(Bs, ...);`
    // (`cxx:5064-5065`).
    let boite = ch_fi3d_boite(&pfac1, &pfac2);
    let hbs = ch_fi3d_bound_fac(
        mc.hbs
            .as_ref()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?,
        boite.u_min,
        boite.u_max,
        boite.v_min,
        boite.v_max,
    );
    mc.hbs = Some(hbs);

    // `if (!ChFi3d_ComputeCurves(HGs, HBs, Pardeb, Parfin, Cc, Ps, Pc, tolapp3d,
    //  tol2d, tolreached)) throw Standard_Failure("OneCorner : failed calculation
    //  intersection");` (`cxx:5067-5070`).
    let mut out = ComputeCurves::default();
    if !ch_fi3d_compute_curves(
        &hgs,
        mc.hbs
            .as_ref()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?,
        pardeb,
        parfin,
        &mut out,
        b.tolapp3d,
        b.tol2d,
        false,
    ) {
        return Err("OneCorner : failed calculation intersection".to_string());
    }
    let cc: Arc<dyn Curve> = out
        .c3d
        .clone()
        .ok_or_else(|| "OneCorner : failed calculation intersection".to_string())?;
    mc.ps = out.pc1.clone();
    mc.pc = out.pc2.clone();
    mc.tolreached = out.tolreached;

    // `Udeb = Cc->FirstParameter(); Ufin = Cc->LastParameter();`
    // (`cxx:5072-5073`).
    mc.udeb = cc.first_parameter();
    mc.ufin = cc.last_parameter();
    mc.cc = Some(Arc::clone(&cc));

    // `ChFi3d_Couture(Fv, couture, edgecouture); if (couture && !BRep_Tool::
    //  Degenerated(edgecouture))` (`cxx:5077-5079`).
    mc.edgecouture = couture(&fv);
    mc.couture = mc.edgecouture.is_some();
    let Some(edgecouture) = mc.edgecouture.clone() else {
        return Ok(());
    };
    if BRepTool::is_degenerated(&edgecouture) {
        return Ok(());
    }

    // `C = BRep_Tool::Curve(edgecouture, Ubid, Vbid); Ctrim = new
    //  Geom_TrimmedCurve(C, Ubid, Vbid); GeomAdaptor_Curve cur1(Ctrim->
    //  BasisCurve()); GeomAdaptor_Curve cur2(Cc);` (`cxx:5083-5086`).
    let c1 = BRepTool::edge_curve(&edgecouture)
        .ok_or_else(|| "IntersectMoreCorner: sewing edge has no curve".to_string())?;
    let cur1 = GeomAdaptorCurve::of_curve(c1);
    let cur2 = GeomAdaptorCurve::of_curve(Arc::clone(&cc));
    let all = curve_curve_extrema_all(&cur1, &cur2);
    // `if (extCC.IsDone() && extCC.NbExt() != 0)` (`cxx:5088`).
    if all.is_empty() {
        return Ok(());
    }

    // `cxx:5090-5099`: the smallest square distance.
    let mut imin = 0usize;
    let mut distmin2 = f64::MAX; // `RealLast()`
    for (i, e) in all.iter().enumerate() {
        let d2 = e.distance * e.distance;
        if d2 < distmin2 {
            distmin2 = d2;
            imin = i;
        }
    }
    // `if (dist2min <= Precision::SquareConfusion())` (`cxx:5100`).
    if distmin2 > CONFUSION * CONFUSION {
        return Ok(());
    }

    // `cxx:5102-5107`.
    let e = &all[imin];
    mc.par1 = e.u1;
    mc.par2 = e.u2;
    let tol = 1.0e-4;
    if (mc.par2 - mc.udeb).abs() <= tol || (mc.ufin - mc.par2).abs() <= tol {
        return Ok(());
    }

    // `TopOpeBRepDS_Point tpoint(P1, Tol); indpt = DStr.AddPoint(tpoint);
    //  intcouture = true;` (`cxx:5109-5112`).
    mc.indpt = b
        .my_ds
        .change_ds()
        .add_point(&Point::with_tolerance(e.p1, tol));
    mc.intcouture = true;
    // `curv1 = new Geom_TrimmedCurve(Cc, Udeb, par2); curv2 = ...` — see the
    // module note on the `GeomAdaptor_Curve` stand-in.
    let curv1: Arc<dyn Curve> =
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cc), mc.udeb, mc.par2));
    let curv2: Arc<dyn Curve> =
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cc), mc.par2, mc.ufin));
    // `TopOpeBRepDS_Curve tcurv1(curv1, tolreached); ... Icurv1 =
    //  DStr.AddCurve(tcurv1); Icurv2 = DStr.AddCurve(tcurv2);`
    // (`cxx:5115-5118`).
    let tolreached = mc.tolreached;
    let ds = b.my_ds.change_ds();
    mc.icurv1 = ds.add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&curv1)),
        tolreached,
        false,
    ));
    mc.icurv2 = ds.add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&curv2)),
        tolreached,
        false,
    ));
    mc.curv1 = Some(curv1);
    mc.curv2 = Some(curv2);
    Ok(())
}

/// `Hc = BRep_Tool::CurveOnSurface(CV.Arc(), Fv, Ubid, Ubid); Hc->Value(
/// CV.ParameterOnArc());` (`cxx:5018-5023` / `:5031-5036`).
fn corner_p2d(
    fv: &crate::shape::Face,
    cv: &crate::chfi_ds::CommonPoint,
) -> Result<occt_core::gp::GpPnt2d, String> {
    let arc = cv
        .arc()
        .cloned()
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let (hc, _, _) = curve_on_surface_range(&Edge(arc), fv)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let par = cv
        .parameter_on_arc()
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    Ok(hc.value(par))
}

/// `Fi.PCurveOnSurf()->Value(Fi.Parameter(isfirst))`, with `!isfirst` when
/// `LineIndex() == 0` (`cxx:5040-5052`).
fn surf_pcurve_point(
    fi: &crate::chfi_ds::FaceInterference,
    is_first: bool,
) -> Result<occt_core::gp::GpPnt2d, String> {
    let pc = fi
        .p_curve_on_surf()
        .ok_or_else(|| "IntersectMoreCorner: Fi has no p-curve on surf".to_string())?;
    let par = if fi.line_index() != 0 {
        fi.parameter(is_first)
    } else {
        fi.parameter(!is_first)
    };
    Ok(pc.value(par))
}
