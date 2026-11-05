//! `IntersectMoreCorner`'s tail — `ChFi3d_Builder_C1.cxx:5243-5445`.
//!
//! - `cxx:5243` — `ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2)`;
//! - `cxx:5245-5423` — the "missing end of curve" `IZob` between `Bop` and
//!   `Bs`: the `Arcspine` vertex interference, the `Arcprolbis` search, the
//!   `ChFi3d_ComputeCurves(Bop, Bs, ...)` intersection, and the `IZob`
//!   interferences on `Fv` / `Fop`;
//! - `cxx:5424-5440` — the final `ChFi3d_EnlargeBox` calls and the point
//!   tolerances.
//!
//! Divergences that are not heuristics:
//! - `ChFi3d_Couture(Fop, couture, edgecouture)` (`cxx:5316`) is called and its
//!   results are never read afterwards in this body; the port keeps the call
//!   and drops the values.
//! - `ComputeCurve2d` is ported (`builder_c1/compute_curve2d.rs`) but this body
//!   does not reach the split write, so no `c2d1` / `c2d2` is taken here.
//! - `Geom_TrimmedCurve`-valued locals are not built here; `zob3d` keeps its own
//!   parameters as in `PerformOneCorner`'s `VARIANT 1`.

use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_boite, ch_fi3d_bound_fac, ch_fi3d_compute_curves, ch_fi3d_enlarge_box_edge,
    ch_fi3d_enlarge_box_pcurve, ch_fi3d_enlarge_box_surf_data, ch_fi3d_fil_curve_in_ds,
    ch_fi3d_fil_point_in_ds, ch_fi3d_fil_vertex_in_ds, ch_fi3d_set_point_tolerance, ComputeCurves,
};
use crate::chfi3d::builder_c1::ch_fi3d_recale;
use crate::chfi3d::builder_helpers::couture;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::ds::Curve as DsCurve;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;
use super::MoreCorner;

/// `ChFi3d_BoundFac(B, uu1, uu2, vv1, vv2)` on a face adaptor, answering the
/// new adaptor (`cxx:5372`, `:5374`).
fn bound_fac(
    adaptor: Option<SurfaceAdaptor>,
    p1: &occt_core::gp::GpPnt2d,
    p2: &occt_core::gp::GpPnt2d,
) -> Option<SurfaceAdaptor> {
    let b_adaptor = adaptor?;
    let boite = ch_fi3d_boite(p1, p2);
    Some(ch_fi3d_bound_fac(
        &b_adaptor,
        boite.u_min,
        boite.u_max,
        boite.v_min,
        boite.v_max,
    ))
}

/// `cxx:5245-5423` — the missing end of curve between `Bop` and `Bs`.
fn missing_curve(b: &mut ChFi3dBuilder, mc: &mut MoreCorner) -> Result<(), String> {
    let is_first = mc.is_first;
    let fv = mc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let fop = mc
        .fop
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let fopbis = mc
        .fopbis
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let arcspine = mc.arcspine.clone();
    let arcprol = mc.arcprol.clone();

    // `int IArcspine = DStr.AddShape(Arcspine); int IVtx = DStr.AddShape(Vtx);
    //  TopAbs_Orientation OVtx2 = OVtx = TopAbs_FORWARD;` (`cxx:5253-5256`).
    let i_arcspine = b.my_ds.change_ds().add_shape(&arcspine);
    let ivtx = b.my_ds.change_ds().add_shape(&mc.vtx);
    let mut ovtx = Orientation::Forward;
    let mut ovtx2 = Orientation::Forward;
    // `for (ex.Init(Arcspine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ...)`
    // (`cxx:5257-5264`).
    for v in explorer_shapes(&arcspine.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &mc.vtx) {
            ovtx = v.orientation();
            break;
        }
    }
    ovtx = ovtx.reversed();
    // `parVtx = BRep_Tool::Parameter(Vtx, Arcspine); interfv =
    //  ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx, parVtx);`
    // (`cxx:5266-5269`).
    let par_vtx = BRepTool::parameter_of_vertex(&arcspine, &mc.vtx);
    let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, ivtx, par_vtx);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&arcspine)
        .push(interfv);

    // `cxx:5272-5283`: the suite of `Arcprol` in the other face.
    let mut arcprolbis = None;
    for item in b.my_ve_map.get_from_key(&mc.vtx) {
        let same = is_same(item, &arcprol)
            || is_same(item, &arcspine)
            || mc.arcpiv.as_ref().map(|a| is_same(item, a)).unwrap_or(false);
        if !same {
            arcprolbis = Some(item.clone());
            break;
        }
    }
    // `if (Arcprolbis.IsNull()) throw Standard_ConstructionError(
    //  "IntersectMoreCorner: continuation edge (Arcprolbis) is not found");`
    // (`cxx:5289-5293`).
    let Some(arcprolbis) = arcprolbis else {
        return Err(
            "IntersectMoreCorner: continuation edge (Arcprolbis) is not found".to_string(),
        );
    };
    mc.arcprolbis = Some(arcprolbis.clone());

    // `cxx:5296-5311`: the two vertex orientations.
    for v in explorer_shapes(&arcprolbis.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &mc.vtx) {
            ovtx2 = v.orientation();
            break;
        }
    }
    for v in explorer_shapes(&arcprol.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &mc.vtx) {
            ovtx = v.orientation();
            break;
        }
    }

    // `ChFi3d_Couture(Fop, couture, edgecouture);` (`cxx:5316`) — the results are
    // not read again in this body.
    let _ = couture(&fop);

    // `const ChFiDS_FaceInterference& Fiop = Fd->Interference(IFopArc);`
    // (`cxx:5319`).
    let fiop = mc.fd.interference(mc.i_fop_arc).clone();
    // `parVtx = BRep_Tool::Parameter(Vtx, Arcprolbis);` (`cxx:5322`).
    let mut par_vtx = BRepTool::parameter_of_vertex(&arcprolbis, &mc.vtx);
    // `if (Fop.IsSame(Fopbis)) OArcprolbis = OArcprolop; else { for (ex.Init(
    //  Fop, TopAbs_EDGE); ...) if (Arcprolbis.IsSame(ex.Current()))
    //  OArcprolbis = ex.Current().Orientation(); }` (`cxx:5326-5340`).
    mc.o_arcprolbis = if is_same(&fop.0, &fopbis.0) {
        mc.o_arcprol_op
    } else {
        let mut ori = Orientation::Forward;
        for e in crate::chfi3d::builder::perform_one_corner::setup::face_edges(&fop) {
            if is_same(&arcprolbis, &e) {
                ori = e.orientation();
                break;
            }
        }
        ori
    };

    // `Hc = BRep_Tool::CurveOnSurface(Arcprolbis, Fop, Ubid, Ubid); pop1 =
    //  Hc->Value(parVtx); pop2 = Fiop.PCurveOnFace()->Value(Fiop
    //  .Parameter(isfirst));` (`cxx:5343-5349`).
    let (hc, _, _) = curve_on_surface_range(&Edge(arcprolbis.clone()), &fop)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let pop1 = hc.value(par_vtx);
    let pop2 = fiop
        .p_curve_on_face()
        .ok_or_else(|| "IntersectMoreCorner: Fiop has no p-curve on face".to_string())?
        .value(fiop.parameter(is_first));

    // `Hc = BRep_Tool::CurveOnSurface(Arcprol, Fv, Ubid, Ubid); parVtx =
    //  BRep_Tool::Parameter(Vtx, Arcprol); pv1 = Hc->Value(parVtx); pv2 =
    //  p2dbout; ChFi3d_Recale(Bs, pv1, pv2, true);` (`cxx:5350-5360`).
    let (hc, _, _) = curve_on_surface_range(&Edge(arcprol.clone()), &fv)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    par_vtx = BRepTool::parameter_of_vertex(&arcprol, &mc.vtx);
    let mut pv1 = hc.value(par_vtx);
    let mut pv2 = mc.p2dbout;
    {
        let hbs = mc
            .hbs
            .as_ref()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
        ch_fi3d_recale(hbs, &mut pv1, &mut pv2, true);
    }

    // `cxx:5361-5374`: both bounding boxes.
    let pardeb = [pop1.x(), pop1.y(), pv1.x(), pv1.y()];
    let parfin = [pop2.x(), pop2.y(), pv2.x(), pv2.y()];
    mc.hbs = bound_fac(mc.hbs.take(), &pv1, &pv2);
    mc.hbop = bound_fac(mc.hbop.take(), &pop1, &pop2);

    // `if (!ChFi3d_ComputeCurves(HBop, HBs, Pardeb, Parfin, zob3d, zob2dop,
    //  zob2dv, tolapp3d, tol2d, tolreached)) throw Standard_Failure("OneCorner :
    //  echec calcul intersection");` (`cxx:5379-5391`).
    let mut out = ComputeCurves::default();
    let hs = mc
        .hbs
        .as_ref()
        .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
    let hop = mc
        .hbop
        .as_ref()
        .ok_or_else(|| "IntersectMoreCorner: Bop has no surface".to_string())?;
    if !ch_fi3d_compute_curves(
        hop,
        hs,
        pardeb,
        parfin,
        &mut out,
        b.tolapp3d,
        b.tol2d,
        false,
    ) {
        return Err("OneCorner : echec calcul intersection".to_string());
    }
    let zob3d = out
        .c3d
        .clone()
        .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
    // `zob2dop` is the p-curve on `HBop` (`Fop`), `zob2dv` the one on `HBs`.
    let zob2dop = out.pc1.clone();
    let zob2dv = out.pc2.clone();
    mc.tolreached = out.tolreached;

    // `Udeb = zob3d->FirstParameter(); Ufin = zob3d->LastParameter();
    //  TopOpeBRepDS_Curve Zob(zob3d, tolreached); int IZob =
    //  DStr.AddCurve(Zob);` (`cxx:5393-5396`).
    mc.udeb = zob3d.first_parameter();
    mc.ufin = zob3d.last_parameter();
    let izob = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&zob3d)),
        mc.tolreached,
        false,
    ));

    // `cxx:5400-5421`.
    let ishape = b.my_ds.change_ds().add_shape(&fv.0);
    let iop = b.my_ds.change_ds().add_shape(&fop.0);
    // `Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv)); InterFv =
    //  ChFi3d_FilCurveInDS(IZob, IShape, zob2dv, Et);`
    let et = Orientation::compose(ovtx, mc.o_arcprol_v).reversed();
    let inter_fv = ch_fi3d_fil_curve_in_ds(izob, ishape, zob2dv, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .push(inter_fv);
    // `Et = TopAbs::Reverse(TopAbs::Compose(OVtx2, OArcprolbis)); Interfop =
    //  ChFi3d_FilCurveInDS(IZob, Iop, zob2dop, Et);`
    let et = Orientation::compose(ovtx2, mc.o_arcprolbis).reversed();
    let inter_fop = ch_fi3d_fil_curve_in_ds(izob, iop, zob2dop, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fop.0)
        .push(inter_fop);
    // `interfprol = ChFi3d_FilVertexInDS(TopAbs_FORWARD, IZob, IVtx, Udeb);`
    let i1 = ch_fi3d_fil_vertex_in_ds(Orientation::Forward, izob, ivtx, mc.udeb);
    b.my_ds
        .change_ds()
        .change_curve_interferences(izob)
        .push(i1);
    // `int icc = stripe->IndexPoint(isfirst, IFopArc); interfprol =
    //  ChFi3d_FilPointInDS(TopAbs_REVERSED, IZob, icc, Ufin);`
    let icc = {
        let stripe = mc
            .stripe
            .read()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        stripe.index_point(is_first, mc.i_fop_arc)
    };
    let i2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, izob, icc, mc.ufin, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(izob)
        .push(i2);
    Ok(())
}

/// `cxx:5243-5445`.
pub(super) fn extend(b: &mut ChFi3dBuilder, mc: &mut MoreCorner) -> Result<(), String> {
    let is_first = mc.is_first;
    // `ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2);` (`cxx:5243`).
    {
        let hbs = mc
            .hbs
            .as_ref()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
        let surf = hbs
            .surface()
            .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
        let pc = mc
            .pc
            .as_ref()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        ch_fi3d_enlarge_box_pcurve(surf, pc, mc.udeb, mc.ufin, &mut mc.box1, &mut mc.box2);
    }

    // `if (inters)` (`cxx:5245`) — the body throws when `inters` is false, so
    // this arm is the only one reached from `intersect_more_corner`.
    missing_curve(b, mc)?;

    // `ChFi3d_EnlargeBox(DStr, stripe, Fd, box1, box2, isfirst);` (`cxx:5424`).
    {
        let stripe = mc
            .stripe
            .read()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        ch_fi3d_enlarge_box_surf_data(
            b.my_ds.change_ds(),
            Some(&stripe),
            &mc.fd,
            &mut mc.box1,
            &mut mc.box2,
            is_first,
        )?;
    }

    let cv1 = mc.fd.vertex(is_first, 1).clone();
    let cv2 = mc.fd.vertex(is_first, 2).clone();
    // `if (CV1.IsOnArc()) ChFi3d_EnlargeBox(CV1.Arc(), myEFMap(CV1.Arc()),
    //  CV1.ParameterOnArc(), box1);` and `CV2` / `box2` (`cxx:5425-5432`).
    if let (Some(arc), Some(par)) = (cv1.arc(), cv1.parameter_on_arc()) {
        let lf: Vec<TopoShape> = b.my_ef_map.get_from_key(arc).to_vec();
        ch_fi3d_enlarge_box_edge(arc, &lf, par, &mut mc.box1);
    }
    if let (Some(arc), Some(par)) = (cv2.arc(), cv2.parameter_on_arc()) {
        let lf: Vec<TopoShape> = b.my_ef_map.get_from_key(arc).to_vec();
        ch_fi3d_enlarge_box_edge(arc, &lf, par, &mut mc.box2);
    }

    // `if (!CV1.IsVertex()) ChFi3d_SetPointTolerance(DStr, box1, stripe->
    //  IndexPoint(isfirst, 1));` and `CV2` / `box2` (`cxx:5433-5440`).
    if !cv1.is_vertex() {
        let ip = {
            let stripe = mc
                .stripe
                .read()
                .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, 1)
        };
        ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &mc.box1, ip);
    }
    if !cv2.is_vertex() {
        let ip = {
            let stripe = mc
                .stripe
                .read()
                .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, 2)
        };
        ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &mc.box2, ip);
    }
    Ok(())
}
