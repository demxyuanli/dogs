//! `PerformOneCorner`'s tail — `ChFi3d_Builder_C1.cxx:1394-1647`.
//!
//! - `cxx:1394` — `ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2)`;
//! - `cxx:1396-1630` — the `onsame && inters` block: the `VARIANT 1` missing end
//!   of curve (`IZob`) between `Bop` and `Bs`, its optional splitting on the
//!   sewing edge of `Fop`, and the `IShape` / `Iop` / `Icurv1` / `Icurv2`
//!   interferences;
//! - `cxx:1631-1646` — the final `ChFi3d_EnlargeBox` calls and the point
//!   tolerances.
//!
//! Divergences that are not heuristics:
//! - `BRepLib_MakeEdge Bedge(zob3d); BRepExtrema_ExtCC(edgecouture, edg)`
//!   (`cxx:1505-1508`) is [`curve_curve_extrema_all`] between `BRepAdaptor_Curve`
//!   of `edgecouture` and `GeomAdaptor_Curve` of `zob3d` — the two adaptors
//!   `BRepExtrema_ExtCC` builds internally for the two edges. The made edge is
//!   never used again (`cxx:1508` only).
//! - `Geom_TrimmedCurve(zob3d, Udeb, par2)` (`cxx:1522-1523`) is
//!   [`GeomAdaptorCurve::with_range`], as in the rest of `perform_one_corner`.
//! - `ComputeCurve2d` is ported (`builder_c1/compute_curve2d.rs`), so `c2d1`
//!   / `c2d2` are the real projected p-curves (`cxx:1537-1541`, `:1569-1575`).
//! - `VARIANT 2` (`cxx:1399`, `:1441-1443`, `:1600-1602`, `:1621-1627`) is not
//!   compiled: OCCT `#define`s `VARIANT1` (`cxx:1403`), and `V2` is only read
//!   by `VARIANT 2`. The `V2` search of `cxx:1424-1435` is kept because it is
//!   also where `OVtx` is re-read from `Arcprol`.

use std::sync::Arc;

use occt_geom::extrema_cc::curve_curve_extrema_all;
use occt_geom::Curve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_boite, ch_fi3d_bound_fac, ch_fi3d_compute_curves, ch_fi3d_enlarge_box_edge,
    ch_fi3d_enlarge_box_pcurve, ch_fi3d_enlarge_box_surf_data, ch_fi3d_fil_curve_in_ds,
    ch_fi3d_fil_point_in_ds, ch_fi3d_fil_vertex_in_ds, ch_fi3d_set_point_tolerance, ComputeCurves,
};
use crate::chfi3d::builder_c1::{ch_fi3d_recale, compute_curve2d};
use crate::chfi3d::builder_helpers::couture;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::ds::{Curve as DsCurve, Point};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, GeomAdaptorCurve};

use super::super::ChFi3dBuilder;
use super::OneCorner;

/// `ChFi3d_BoundFac(B, uu1, uu2, vv1, vv2)` on a face adaptor, answering the
/// new adaptor (`cxx:1471`, `:1474`). The OCCT call mutates the adaptor in
/// place; the port's helper returns the bounded adaptor.
fn bound_fac(
    adaptor: Option<crate::topopebrep::tool::adaptor::SurfaceAdaptor>,
    p1: &occt_core::gp::GpPnt2d,
    p2: &occt_core::gp::GpPnt2d,
) -> Option<crate::topopebrep::tool::adaptor::SurfaceAdaptor> {
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

/// `cxx:1396-1630` — the `onsame && inters` "missing end of curve" block.
fn missing_curve(b: &mut ChFi3dBuilder, oc: &mut OneCorner) -> Result<(), String> {
    let is_first = oc.is_first;
    let fv = oc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let fop = oc
        .fop
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let arcspine = oc
        .arcspine
        .clone()
        .ok_or_else(|| "PerformOneCorner: the spine has no edge".to_string())?;
    let arcprol = oc
        .arcprol
        .clone()
        .ok_or_else(|| "OneCorner : edge a prolonger non trouve".to_string())?;

    // `int IArcspine = DStr.AddShape(Arcspine); int IVtx = DStr.AddShape(Vtx);`
    // (`cxx:1405-1406`).
    let i_arcspine = b.my_ds.change_ds().add_shape(&arcspine);
    let ivtx = b.my_ds.change_ds().add_shape(&oc.vtx);
    // `TopAbs_Orientation OVtx = TopAbs_FORWARD; for (ex.Init(Arcspine
    //  .Oriented(FORWARD), VERTEX); ...)` (`cxx:1407-1415`).
    let mut ovtx = Orientation::Forward;
    for v in explorer_shapes(&arcspine.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &oc.vtx) {
            ovtx = v.orientation();
            break;
        }
    }
    ovtx = ovtx.reversed();
    // `parVtx = BRep_Tool::Parameter(Vtx, Arcspine);` (`cxx:1417`).
    let mut par_vtx = BRepTool::parameter_of_vertex(&arcspine, &oc.vtx);
    let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, ivtx, par_vtx);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&arcspine)
        .push(interfv);

    // `cxx:1424-1435`: the `V2` / `OVtx` walk of `Arcprol`. `V2` is only read
    // by `VARIANT 2`, but `OVtx` is re-read here and used below.
    for v in explorer_shapes(&arcprol.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &oc.vtx) {
            ovtx = v.orientation();
        }
    }

    // `parVtx = BRep_Tool::Parameter(Vtx, Arcprol);` (`cxx:1439`, VARIANT 1).
    par_vtx = BRepTool::parameter_of_vertex(&arcprol, &oc.vtx);

    // `const ChFiDS_FaceInterference& Fiop = Fd->Interference(IFopArc);`
    // (`cxx:1444`).
    let fiop = oc.fd.interference(oc.i_fop_arc);
    // `Hc = BRep_Tool::CurveOnSurface(Arcprol, Fop, Ubid, Ubid); pop1 =
    //  Hc->Value(parVtx); pop2 = Fiop.PCurveOnFace()->Value(Fiop
    //  .Parameter(isfirst));` (`cxx:1447-1451`).
    let (hc, _, _) = curve_on_surface_range(&Edge(arcprol.clone()), &fop)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let pop1 = hc.value(par_vtx);
    let pop2 = fiop
        .p_curve_on_face()
        .ok_or_else(|| "OneCorner : Fiop has no p-curve".to_string())?
        .value(fiop.parameter(is_first));
    // `Hc = BRep_Tool::CurveOnSurface(Arcprol, Fv, Ubid, Ubid); pv1 =
    //  Hc->Value(parVtx); pv2 = p2dbout; ChFi3d_Recale(Bs, pv1, pv2, true);`
    // (`cxx:1452-1459`).
    let (hc, _, _) = curve_on_surface_range(&Edge(arcprol.clone()), &fv)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let mut pv1 = hc.value(par_vtx);
    let mut pv2 = oc.p2dbout;
    {
        let hbs = oc
            .hbs
            .as_ref()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
        ch_fi3d_recale(hbs, &mut pv1, &mut pv2, true);
    }

    // `cxx:1460-1474`: both bounding boxes.
    let pardeb = [pop1.x(), pop1.y(), pv1.x(), pv1.y()];
    let parfin = [pop2.x(), pop2.y(), pv2.x(), pv2.y()];
    oc.hbs = bound_fac(oc.hbs.take(), &pv1, &pv2);
    oc.hbop = bound_fac(oc.hbop.take(), &pop1, &pop2);

    // `ChFi3d_ComputeCurves(HBop, HBs, Pardeb, Parfin, zob3d, zob2dop, zob2dv,
    //  tolapp3d, tol2d, tolreached)` (`cxx:1477-1491`).
    let mut out = ComputeCurves::default();
    let hs = oc
        .hbs
        .as_ref()
        .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
    let hop = oc
        .hbop
        .as_ref()
        .ok_or_else(|| "PerformOneCorner: Bop has no surface".to_string())?;
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
    // `zob2dop` is the p-curve on the first surface (`HBop`), `zob2dv` the one
    // on `HBs`.
    let zob2dop = out.pc1.clone();
    let zob2dv = out.pc2.clone();
    oc.tolreached = out.tolreached;
    // `Udeb = zob3d->FirstParameter(); Ufin = zob3d->LastParameter();
    //  TopOpeBRepDS_Curve Zob(zob3d, tolreached); int IZob = DStr.AddCurve(Zob);`
    // (`cxx:1493-1496`).
    oc.udeb = zob3d.first_parameter();
    oc.ufin = zob3d.last_parameter();
    let izob = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&zob3d)),
        oc.tolreached,
        false,
    ));

    // `ChFi3d_Couture(Fop, couture, edgecouture); if (couture && !BRep_Tool::
    //  Degenerated(edgecouture))` (`cxx:1502-1504`).
    let edgecouture = couture(&fop);
    if let Some(ec) = edgecouture.as_ref() {
        if !BRepTool::is_degenerated(ec) {
            // `cxx:1505-1508`: `BRepLib_MakeEdge Bedge(zob3d); edg = Bedge.Edge();
            //  BRepExtrema_ExtCC extCC(edgecouture, edg);`.
            let e1 = CurveAdaptor::of_edge(ec);
            let e2 = GeomAdaptorCurve::of_curve(Arc::clone(&zob3d));
            let all = curve_curve_extrema_all(&e1, &e2);
            for e in all {
                // `for (int i = 1; i <= extCC.NbExt() && !intcouture; i++)`
                // (`cxx:1511`).
                if oc.intcouture {
                    break;
                }
                // `if (extCC.SquareDistance(i) <= 1.e-8)` (`cxx:1513`).
                if e.distance * e.distance > 1e-8 {
                    continue;
                }
                oc.par1 = e.u1;
                oc.par2 = e.u2;
                // `TopOpeBRepDS_Point tpoint(P1, 1.e-4); indpt = DStr.AddPoint(...)`
                // (`cxx:1517-1519`).
                oc.indpt = b
                    .my_ds
                    .change_ds()
                    .add_point(&Point::with_tolerance(e.p1, 1e-4));
                oc.intcouture = true;
                // `curv1 = new Geom_TrimmedCurve(zob3d, Udeb, par2); curv2 =
                //  new Geom_TrimmedCurve(zob3d, par2, Ufin);` (`cxx:1522-1523`).
                let curv1: Arc<dyn Curve> =
                    Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&zob3d), oc.udeb, oc.par2));
                let curv2: Arc<dyn Curve> =
                    Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&zob3d), oc.par2, oc.ufin));
                let tolreached = oc.tolreached;
                let ds = b.my_ds.change_ds();
                oc.icurv1 = ds.add_curve(&mut DsCurve::with_curve(
                    Some(Arc::clone(&curv1)),
                    tolreached,
                    false,
                ));
                oc.icurv2 = ds.add_curve(&mut DsCurve::with_curve(
                    Some(Arc::clone(&curv2)),
                    tolreached,
                    false,
                ));
                oc.curv1 = Some(curv1);
                oc.curv2 = Some(curv2);
            }
        }
    }

    // `IShape` / `Iop` are `DStr.AddShape` of `Fv` / `Fop`
    // (`cxx:1066`, `:1566`, `:1592`).
    let ishape = b.my_ds.change_ds().add_shape(&fv.0);
    let iop = b.my_ds.change_ds().add_shape(&fop.0);

    if oc.intcouture {
        // `Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));`
        // (`cxx:1534`).
        let et = Orientation::compose(ovtx, oc.o_arcprol_v).reversed();
        let curv1 = oc
            .curv1
            .clone()
            .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
        let curv2 = oc
            .curv2
            .clone()
            .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
        // `ComputeCurve2d(curv1, Fop, c2d1); InterFv = ChFi3d_FilCurveInDS(Icurv1,
        //  IShape, c2d1, Et);` (`cxx:1535-1539`).
        let c2d1 = compute_curve2d(&curv1, &fop);
        let inter = ch_fi3d_fil_curve_in_ds(oc.icurv1, ishape, c2d1, et);
        b.my_ds.change_ds().change_shape_interferences(&fv.0).push(inter);
        // `cxx:1540-1542`.
        let c2d2 = compute_curve2d(&curv2, &fop);
        let inter = ch_fi3d_fil_curve_in_ds(oc.icurv2, ishape, c2d2, et);
        b.my_ds.change_ds().change_shape_interferences(&fv.0).push(inter);

        // `cxx:1544-1563`: limitation of the sewing edge.
        let ec = edgecouture
            .as_ref()
            .ok_or_else(|| "OneCorner : no sewing edge".to_string())?;
        let iarc = b.my_ds.change_ds().add_shape(&ec.0);
        let vdeb = first_vertex(&ec.0)
            .ok_or_else(|| "OneCorner : sewing edge has no first vertex".to_string())?;
        let vfin = last_vertex(&ec.0)
            .ok_or_else(|| "OneCorner : sewing edge has no last vertex".to_string())?;
        let pard = BRepTool::parameter_of_vertex(&ec.0, &vdeb);
        let parf = BRepTool::parameter_of_vertex(&ec.0, &vfin);
        // `if (std::abs(par1 - pard) < std::abs(parf - par1)) ori = REVERSED;
        //  else ori = FORWARD;` (`cxx:1555-1562`).
        let ori = if (oc.par1 - pard).abs() < (parf - oc.par1).abs() {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        let inter = ch_fi3d_fil_point_in_ds(ori, iarc, oc.indpt, oc.par1, false);
        b.my_ds
            .change_ds()
            .change_shape_interferences(&ec.0)
            .push(inter);

        // `cxx:1565-1576`: both pieces on `Iop`, still with the p-curves taken
        // from `Fop`.
        let et = Orientation::compose(ovtx, oc.o_arcprol_op).reversed();
        let c2d1 = compute_curve2d(&curv1, &fop);
        let inter = ch_fi3d_fil_curve_in_ds(oc.icurv1, iop, c2d1, et);
        b.my_ds.change_ds().change_shape_interferences(&fop.0).push(inter);
        let c2d2 = compute_curve2d(&curv2, &fop);
        let inter = ch_fi3d_fil_curve_in_ds(oc.icurv2, iop, c2d2, et);
        b.my_ds.change_ds().change_shape_interferences(&fop.0).push(inter);

        // `cxx:1577-1584`.
        let i1 = ch_fi3d_fil_vertex_in_ds(Orientation::Forward, oc.icurv1, ivtx, oc.udeb);
        b.my_ds
            .change_ds()
            .change_curve_interferences(oc.icurv1)
            .push(i1);
        let i2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, oc.icurv1, oc.indpt, oc.par2, false);
        b.my_ds
            .change_ds()
            .change_curve_interferences(oc.icurv1)
            .push(i2);
        let icc = {
            let stripe = oc
                .stripe
                .read()
                .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, oc.i_fop_arc)
        };
        let i3 = ch_fi3d_fil_point_in_ds(Orientation::Forward, oc.icurv2, oc.indpt, oc.par2, false);
        b.my_ds
            .change_ds()
            .change_curve_interferences(oc.icurv2)
            .push(i3);
        let i4 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, oc.icurv2, icc, oc.ufin, false);
        b.my_ds
            .change_ds()
            .change_curve_interferences(oc.icurv2)
            .push(i4);
    } else {
        // `Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));
        //  InterFv = ChFi3d_FilCurveInDS(IZob, IShape, zob2dv, Et);`
        // (`cxx:1588-1591`).
        let et = Orientation::compose(ovtx, oc.o_arcprol_v).reversed();
        let inter = ch_fi3d_fil_curve_in_ds(izob, ishape, zob2dv, et);
        b.my_ds.change_ds().change_shape_interferences(&fv.0).push(inter);
        // `Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolop)); Iop =
        //  DStr.AddShape(Fop); Interfop = ChFi3d_FilCurveInDS(IZob, Iop,
        //  zob2dop, Et);` (`cxx:1592-1597`).
        let et = Orientation::compose(ovtx, oc.o_arcprol_op).reversed();
        let inter = ch_fi3d_fil_curve_in_ds(izob, iop, zob2dop, et);
        b.my_ds.change_ds().change_shape_interferences(&fop.0).push(inter);
        // `interfprol = ChFi3d_FilVertexInDS(FORWARD, IZob, IVtx, Udeb);`
        // (`cxx:1599`).
        let i1 = ch_fi3d_fil_vertex_in_ds(Orientation::Forward, izob, ivtx, oc.udeb);
        b.my_ds
            .change_ds()
            .change_curve_interferences(izob)
            .push(i1);
        // `int icc = stripe->IndexPoint(isfirst, IFopArc); interfprol =
        //  ChFi3d_FilPointInDS(REVERSED, IZob, icc, Ufin);` (`cxx:1607-1609`).
        let icc = {
            let stripe = oc
                .stripe
                .read()
                .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, oc.i_fop_arc)
        };
        let i2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, izob, icc, oc.ufin, false);
        b.my_ds
            .change_ds()
            .change_curve_interferences(izob)
            .push(i2);
        // `if (IFopArc == 1) box1.Add(zob3d->Value(Ufin)); else box2.Add(...)`
        // (`cxx:1611-1619`).
        let p = zob3d.d0(oc.ufin);
        if oc.i_fop_arc == 1 {
            oc.box1.add_point(&p);
        } else {
            oc.box2.add_point(&p);
        }
    }
    Ok(())
}

/// `cxx:1394-1647`.
pub(super) fn extend(b: &mut ChFi3dBuilder, oc: &mut OneCorner) -> Result<(), String> {
    let is_first = oc.is_first;
    // `ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2);` (`cxx:1394`).
    {
        let hbs = oc
            .hbs
            .as_ref()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
        let surf = hbs
            .surface()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
        let pc = oc
            .pc
            .as_ref()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        ch_fi3d_enlarge_box_pcurve(
            surf,
            pc,
            oc.udeb,
            oc.ufin,
            &mut oc.box1,
            &mut oc.box2,
        );
    }

    // `if (onsame && inters)` (`cxx:1396`).
    if oc.onsame && oc.inters {
        missing_curve(b, oc)?;
    }

    // `ChFi3d_EnlargeBox(DStr, stripe, Fd, box1, box2, isfirst);` (`cxx:1631`).
    {
        let stripe = oc
            .stripe
            .read()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        ch_fi3d_enlarge_box_surf_data(
            b.my_ds.change_ds(),
            Some(&stripe),
            &oc.fd,
            &mut oc.box1,
            &mut oc.box2,
            is_first,
        )?;
    }

    let cv1 = oc.fd.vertex(is_first, 1).clone();
    let cv2 = oc.fd.vertex(is_first, 2).clone();
    // `if (CV1.IsOnArc()) ChFi3d_EnlargeBox(CV1.Arc(), myEFMap(CV1.Arc()),
    //  CV1.ParameterOnArc(), box1);` and `CV2` / `box2` (`cxx:1632-1638`).
    if let (Some(arc), Some(par)) = (cv1.arc(), cv1.parameter_on_arc()) {
        let lf: Vec<TopoShape> = b.my_ef_map.get_from_key(arc).to_vec();
        ch_fi3d_enlarge_box_edge(arc, &lf, par, &mut oc.box1);
    }
    if let (Some(arc), Some(par)) = (cv2.arc(), cv2.parameter_on_arc()) {
        let lf: Vec<TopoShape> = b.my_ef_map.get_from_key(arc).to_vec();
        ch_fi3d_enlarge_box_edge(arc, &lf, par, &mut oc.box2);
    }

    // `if (!CV1.IsVertex()) ChFi3d_SetPointTolerance(DStr, box1,
    //  stripe->IndexPoint(isfirst, 1));` and `CV2` / `box2` (`cxx:1639-1646`).
    if !cv1.is_vertex() {
        let ip = {
            let stripe = oc
                .stripe
                .read()
                .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, 1)
        };
        ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &oc.box1, ip);
    }
    if !cv2.is_vertex() {
        let ip = {
            let stripe = oc
                .stripe
                .read()
                .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
            stripe.index_point(is_first, 2)
        };
        ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &oc.box2, ip);
    }
    Ok(())
}
