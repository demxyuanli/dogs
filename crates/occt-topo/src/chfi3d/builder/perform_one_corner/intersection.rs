//! `PerformOneCorner`'s intersection step — `ChFi3d_Builder_C1.cxx:933-1098`.
//!
//! `ChFi3d_BoundSurf` bounds the fillet surface, `ChFi3d_BoundFac` the face,
//! `ChFi3d_ComputeCurves` intersects them into `Cc` with its two p-curves, and
//! the sewing edge of `Fv` is tested against `Cc` with `Extrema_ExtCC`: a hit
//! inside `Precision::SquareConfusion()` splits `Cc` at the sewing parameter
//! (`intcouture`).

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
use super::OneCorner;

/// `PerformOneCorner`'s `if (inters)` block (`cxx:938-1062`).
pub(super) fn intersect(b: &mut ChFi3dBuilder, oc: &mut OneCorner) -> Result<(), String> {
    let is_first = oc.is_first;
    let fv = oc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;

    // `HGs = ChFi3d_BoundSurf(DStr, Fd, 1, 2);` (`cxx:940`).
    let hgs_adaptor = {
        let ds = b.my_ds.change_ds();
        ch_fi3d_bound_surf(ds, &oc.fd, 1, 2)?
    };
    // `const ChFiDS_FaceInterference& Fi1 = Fd->InterferenceOnS1();` and `Fi2`
    // (`cxx:941-942`).
    let fi1 = oc.fd.interference_on_s1().clone();
    let fi2 = oc.fd.interference_on_s2().clone();
    let cv1 = oc.fd.vertex(is_first, 1).clone();
    let cv2 = oc.fd.vertex(is_first, 2).clone();

    // `gp_Pnt2d pfil1, pfac1, pfil2, pfac2;` (`cxx:944`).
    let pfac1 = if oc.onsame && oc.i_fop_arc == 1 {
        // `pfac1 = p2dbout;` (`cxx:948-949`).
        oc.p2dbout
    } else {
        // `Hc1 = BRep_Tool::CurveOnSurface(CV1.Arc(), Fv, Ubid, Ubid); ...
        //  pfac1 = Hc1->Value(CV1.ParameterOnArc());` (`cxx:951-956`).
        let arc = cv1
            .arc()
            .cloned()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let (hc1, _, _) = curve_on_surface_range(&Edge(arc), &fv)
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let par = cv1
            .parameter_on_arc()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        hc1.value(par)
    };
    let pfac2 = if oc.onsame && oc.i_fop_arc == 2 {
        // `pfac2 = p2dbout;` (`cxx:959-960`).
        oc.p2dbout
    } else {
        // `cxx:962-970`.
        let arc = cv2
            .arc()
            .cloned()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let (hc2, _, _) = curve_on_surface_range(&Edge(arc), &fv)
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let par = cv2
            .parameter_on_arc()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        hc2.value(par)
    };

    // `pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));`, with the
    // `!isfirst` arm when `LineIndex() == 0` (`cxx:972-987`).
    let pfil1 = {
        let pc = fi1
            .p_curve_on_surf()
            .ok_or_else(|| "PerformOneCorner: Fi1 has no p-curve on surf".to_string())?;
        let par = if fi1.line_index() != 0 {
            fi1.parameter(is_first)
        } else {
            fi1.parameter(!is_first)
        };
        pc.value(par)
    };
    let pfil2 = {
        let pc = fi2
            .p_curve_on_surf()
            .ok_or_else(|| "PerformOneCorner: Fi2 has no p-curve on surf".to_string())?;
        let par = if fi2.line_index() != 0 {
            fi2.parameter(is_first)
        } else {
            fi2.parameter(!is_first)
        };
        pc.value(par)
    };

    // `if (onsame) { ChFi3d_Recale(Bs, pfac1, pfac2, (IFadArc == 1)); }`
    // (`cxx:989-992`).
    let mut pfac1 = pfac1;
    let mut pfac2 = pfac2;
    if oc.onsame {
        let hbs = oc
            .hbs
            .as_ref()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
        ch_fi3d_recale(hbs, &mut pfac1, &mut pfac2, oc.i_fad_arc == 1);
    }

    // `Pardeb(1 .. 4)` / `Parfin(1 .. 4)` (`cxx:994-1000`).
    let pardeb = [pfil1.x(), pfil1.y(), pfac1.x(), pfac1.y()];
    let parfin = [pfil2.x(), pfil2.y(), pfac2.x(), pfac2.y()];

    // `ChFi3d_Boite(pfac1, pfac2, uu1, uu2, vv1, vv2); ChFi3d_BoundFac(Bs, ...);`
    // (`cxx:1002-1005`).
    let boite = ch_fi3d_boite(&pfac1, &pfac2);
    let hbs = ch_fi3d_bound_fac(
        oc.hbs
            .as_ref()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?,
        boite.u_min,
        boite.u_max,
        boite.v_min,
        boite.v_max,
    );
    oc.hbs = Some(hbs);

    // `if (!ChFi3d_ComputeCurves(HGs, HBs, Pardeb, Parfin, Cc, Ps, Pc, tolapp3d,
    //  tol2d, tolreached)) throw Standard_Failure("OneCorner : echec calcul
    //  intersection");` (`cxx:1007-1010`).
    let mut out = ComputeCurves::default();
    if !ch_fi3d_compute_curves(
        &hgs_adaptor,
        oc.hbs
            .as_ref()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?,
        pardeb,
        parfin,
        &mut out,
        b.tolapp3d,
        b.tol2d,
        false,
    ) {
        return Err("OneCorner : echec calcul intersection".to_string());
    }
    let cc: Arc<dyn Curve> = out
        .c3d
        .clone()
        .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
    oc.ps = out.pc1.clone();
    oc.pc = out.pc2.clone();
    oc.tolreached = out.tolreached;

    // `Udeb = Cc->FirstParameter(); Ufin = Cc->LastParameter();` (`cxx:1012-1013`).
    oc.udeb = cc.first_parameter();
    oc.ufin = cc.last_parameter();
    oc.cc = Some(Arc::clone(&cc));

    // `ChFi3d_Couture(Fv, couture, edgecouture); if (couture &&
    //  !BRep_Tool::Degenerated(edgecouture))` (`cxx:1017-1019`).
    oc.edgecouture = couture(&fv);
    oc.couture = oc.edgecouture.is_some();
    let Some(edgecouture) = oc.edgecouture.clone() else {
        return Ok(());
    };
    if BRepTool::is_degenerated(&edgecouture) {
        return Ok(());
    }

    // `C = BRep_Tool::Curve(edgecouture, Ubid, Vbid); Ctrim = new
    //  Geom_TrimmedCurve(C, Ubid, Vbid); GeomAdaptor_Curve cur1(Ctrim->
    //  BasisCurve()); GeomAdaptor_Curve cur2(Cc);` (`cxx:1022-1026`). The
    //  basis curve keeps its own range, so `cur1` is the plain curve adaptor.
    let c1 = BRepTool::edge_curve(&edgecouture)
        .ok_or_else(|| "PerformOneCorner: sewing edge has no curve".to_string())?;
    let cur1 = GeomAdaptorCurve::of_curve(c1);
    let cur2 = GeomAdaptorCurve::of_curve(Arc::clone(&cc));
    let all = curve_curve_extrema_all(&cur1, &cur2);
    // `if (extCC.IsDone() && extCC.NbExt() != 0)` (`cxx:1028`).
    if all.is_empty() {
        return Ok(());
    }

    // `cxx:1029-1038`: the smallest square distance.
    let mut imin = 0usize;
    let mut distmin2 = f64::MAX; // `RealLast()`
    for (i, e) in all.iter().enumerate() {
        let d2 = e.distance * e.distance;
        if d2 < distmin2 {
            distmin2 = d2;
            imin = i;
        }
    }
    // `if (distmin2 <= Precision::SquareConfusion())` (`cxx:1039`).
    if distmin2 > CONFUSION * CONFUSION {
        return Ok(());
    }

    // `cxx:1041-1047`.
    let e = &all[imin];
    oc.par1 = e.u1;
    oc.par2 = e.u2;
    let tol = 1.0e-4;
    if (oc.par2 - oc.udeb).abs() <= tol || (oc.ufin - oc.par2).abs() <= tol {
        return Ok(());
    }

    // `cxx:1048-1058`.
    oc.indpt = b
        .my_ds
        .change_ds()
        .add_point(&Point::with_tolerance(e.p1, tol));
    oc.intcouture = true;
    // `curv1 = new Geom_TrimmedCurve(Cc, Udeb, par2);` / `curv2 = ...` — see
    // the module note on the `GeomAdaptor_Curve` stand-in.
    let curv1: Arc<dyn Curve> =
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cc), oc.udeb, oc.par2));
    let curv2: Arc<dyn Curve> =
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cc), oc.par2, oc.ufin));
    let tolreached = oc.tolreached;
    let ds = b.my_ds.change_ds();
    oc.icurv1 = ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&curv1)), tolreached, false));
    oc.icurv2 = ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&curv2)), tolreached, false));
    oc.curv1 = Some(curv1);
    oc.curv2 = Some(curv2);
    Ok(())
}
