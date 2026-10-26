//! `PerformOneCorner`'s DS write stage — `ChFi3d_Builder_C1.cxx:1066-1393`.
//!
//! Ported blocks, in OCCT order:
//! - `cxx:1066-1092` — `IShape = DStr.AddShape(Fv)` and the orientation `Et` of
//!   the arc of `CV1` (or `CV2` when `IFadArc != 1`) inside `Fv`;
//! - `cxx:1100-1101` — `SetIndexPoint(ChFi3d_IndexPointInDS(CV1 / CV2, DStr))`;
//! - `cxx:1103-1330` — the `!intcouture` write: `Tc(Cc, tolreached)` in the DS,
//!   the OCC119 overlap guard (`cxx:1116-1204`), the pseudo-seam block
//!   (`cxx:1207-1327`) and the stripe bookkeeping (`cxx:1328-1329`);
//! - `cxx:1332-1392` — the `intcouture` write: both pieces of `Cc` on `Fv` and
//!   on `Isurf`, the limitation of the sewing edge and the four
//!   `CurvePointInterference`s of `Icurv1` / `Icurv2`.
//!
//! Divergences that are not heuristics:
//! - `Geom2dAdaptor_Curve(C, First, Last)` becomes
//!   [`Geom2dIntGInter::domain_over_range`], the same bounded domain
//!   `ComputeDomain` derives from the narrowed adaptor (as
//!   `chfi3d::builder_0_helpers::stripe_edge_inter` does).
//! - `occ::down_cast<Geom_TrimmedCurve>(DStr.Curve(GI).Curve())`
//!   (`cxx:1189-1193`) has no downcast in this port's `Curve` trait, so the
//!   second OCC119 loop uses the stored curve's own `first_parameter()` /
//!   `last_parameter()`; for the `GeomAdaptor_Curve::with_range` curves this
//!   stage writes, those are the OCCT basis limits. A null curve is skipped
//!   like the failed downcast.
//! - `occ::down_cast<Geom2d_TrimmedCurve>(Cd) -> BasisCurve()` (`cxx:1309-1313`)
//!   has no 2D downcast either. `Geom2dTrimmedCurve` keeps the basis
//!   parameterisation in this port, so `par` and `fd` stay in one
//!   parameterisation and the projection is taken on the stored p-curve.
//! - `Geom_TrimmedCurve(Ct, prm1, prm2)` (`cxx:1250`) is
//!   [`GeomAdaptorCurve::with_range`] in the port, for the reason stated in
//!   `perform_one_corner::mod`.
//! - `ComputeCurve2d` is UNPORTED (`builder_c1/compute_curve2d.rs`), so the
//!   `intcouture` write receives `None` for `c2d1` / `c2d2`
//!   (`cxx:1338-1342`, `:1537-1541`).

use std::sync::Arc;

use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api::project_point_on_curve_2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ch_fi3d_index_point_in_ds,
};
use crate::chfi3d::builder_c1::compute_curve2d;
use crate::chfi3d::builder_helpers::is_pseudo_seam;
use crate::shape::{Edge, Face};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::ds::Curve as DsCurve;
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

use super::super::ChFi3dBuilder;
use super::OneCorner;

/// `IsEqual(A, B)` of `Standard_Real.hxx:148-151`.
fn is_equal_real(a: f64, b: f64) -> bool {
    (a - b).abs() < REAL_SMALL
}

/// `cxx:1066-1092`: `IShape` and `Et`.
fn face_index_and_orientation(
    oc: &OneCorner,
    b: &mut ChFi3dBuilder,
    fv: &Face,
) -> Result<(i32, Orientation), String> {
    // `int IShape = DStr.AddShape(Fv);` (`cxx:1066`).
    let ishape = b.my_ds.change_ds().add_shape(&fv.0);
    // `TopAbs_Orientation Et = TopAbs_FORWARD;` (`cxx:1067`).
    let mut et = Orientation::Forward;
    let on_s = if oc.i_fad_arc == 1 { 1 } else { 2 };
    let cv = oc.fd.vertex(oc.is_first, on_s);
    let Some(target) = cv.arc() else {
        // `CV1.Arc()` is not compared when the common point is off arc; the
        // loop body never matches and `Et` stays FORWARD.
        return Ok((ishape, et));
    };
    let trans = cv
        .transition_on_arc()
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let forward = fv.0.oriented(Orientation::Forward);
    for e in explorer_shapes(&forward, ShapeType::Edge) {
        if is_same(&e, target) {
            let composed = Orientation::compose(e.orientation(), trans);
            // `cxx:1076-1078` reverses, `cxx:1088-1090` does not.
            et = if oc.i_fad_arc == 1 {
                composed.reversed()
            } else {
                composed
            };
            break;
        }
    }
    Ok((ishape, et))
}

/// `cxx:1116-1204` — the OCC119 "fillets have too big radiuses" guard.
fn check_cork_overlap(
    oc: &OneCorner,
    b: &mut ChFi3dBuilder,
    ishape: i32,
    fv: &Face,
    pc: &Arc<dyn Curve2d>,
) -> Result<(), String> {
    // `Geom2dAdaptor_Curve aCorkPCurve(Pc, Udeb, Ufin);` (`cxx:1118`).
    let cork = pc.as_ref();
    let tol_domain = b.tol2d.max(PCONFUSION);
    let d_cork = Geom2dIntGInter::domain_over_range(cork, oc.udeb, oc.ufin, tol_domain);
    let mut an_intersector = Geom2dIntGInter::new();

    // `cxx:1122-1176`: every part of every stripe, on whichever of its two
    // faces `IShape` is.
    let mut others: Vec<(Arc<dyn Curve2d>, f64, f64)> = Vec::new();
    for handle in b.my_list_stripe.iter() {
        let stripe = handle
            .read()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        for data in stripe.set_of_surf_data() {
            let fi = if ishape == data.index_of_s1() {
                data.interference_on_s1()
            } else if ishape == data.index_of_s2() {
                data.interference_on_s2()
            } else {
                // `cxx:1159-1162`: no common surface.
                continue;
            };
            let Some(p) = fi.p_curve_on_face() else {
                // `cxx:1138-1141` / `:1150-1153`.
                continue;
            };
            let (first, last) = (fi.first_parameter(), fi.last_parameter());
            // `cxx:1164-1168`: degenerates.
            if is_equal_real(last, first) {
                continue;
            }
            others.push((Arc::clone(p), first, last));
        }
    }
    for (p, first, last) in others {
        let other = p.as_ref();
        let d_other = Geom2dIntGInter::domain_over_range(other, first, last, tol_domain);
        an_intersector.perform(cork, &d_cork, other, &d_other, b.tol2d, PCONFUSION);
        if an_intersector.nb_segments() > 0 || an_intersector.nb_points() > 0 {
            return Err("OneCorner : fillets have too big radiuses".to_string());
        }
    }

    // `cxx:1177-1204`: the interferences already stored on `IShape`.
    let stored = b
        .my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .clone();
    for handle in stored {
        let other = handle
            .read()
            .map_err(|_| "PerformOneCorner: poisoned interference".to_string())?;
        // `occ::down_cast<TopOpeBRepDS_SurfaceCurveInterference>` (`cxx:1180`).
        if other.as_surface_curve().is_none() {
            continue;
        }
        // `occ::down_cast<Geom_TrimmedCurve>(DStr.Curve(GI).Curve())`
        // (`cxx:1189-1193`) — see the module note.
        let Some(cur) = b.my_ds.change_ds().curve(other.geometry()).curve() else {
            continue;
        };
        let (first, last) = (cur.first_parameter(), cur.last_parameter());
        let Some(p) = other.p_curve() else {
            continue;
        };
        let op = p.as_ref();
        let d_other = Geom2dIntGInter::domain_over_range(op, first, last, tol_domain);
        an_intersector.perform(cork, &d_cork, op, &d_other, b.tol2d, PCONFUSION);
        if an_intersector.nb_segments() > 0 || an_intersector.nb_points() > 0 {
            return Err("OneCorner : fillets have too big radiuses".to_string());
        }
    }
    Ok(())
}

/// `cxx:1207-1327` — the pseudo seam piece of the straightened curve.
fn write_pseudo_seam(
    oc: &mut OneCorner,
    b: &mut ChFi3dBuilder,
    fv: &Face,
    pc: &Arc<dyn Curve2d>,
) -> Result<(), String> {
    let is_first = oc.is_first;
    let cv1 = oc.fd.vertex(is_first, 1).clone();
    let cv2 = oc.fd.vertex(is_first, 2).clone();
    for i in 0..2 {
        // `CV[0] = CV1; CV[1] = CV2;` (`cxx:1207-1209`).
        let cv = if i == 0 { &cv1 } else { &cv2 };
        if !cv.is_on_arc() {
            continue;
        }
        let arc = cv
            .arc()
            .cloned()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        // `if (CV[i].IsOnArc() && ChFi3d_IsPseudoSeam(CV[i].Arc(), Fv))`
        // (`cxx:1211`).
        if !is_pseudo_seam(&Edge(arc.clone()), fv) {
            continue;
        }
        // `Hc = BRep_Tool::CurveOnSurface(CV[i].Arc(), Fv, first, last);`
        // (`cxx:1218`).
        let (hc, first_pc, last_pc) = curve_on_surface_range(&Edge(arc.clone()), fv)
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let par_arc = cv
            .parameter_on_arc()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        // `pfac1 = Hc->Value(CV[i].ParameterOnArc());` (`cxx:1223`).
        let pfac1 = hc.value(par_arc);
        // `PcF = Pc->Value(Udeb); PcL = Pc->Value(Ufin);` (`cxx:1224-1225`).
        let pcf = pc.value(oc.udeb);
        let pcl = pc.value(oc.ufin);
        // `onfirst = pfac1.Distance(PcF) < pfac1.Distance(PcL);` (`cxx:1226`).
        let onfirst = pfac1.distance(&pcf) < pfac1.distance(&pcl);
        // `cxx:1227-1235`.
        let der_pc = if onfirst {
            pc.d1(oc.udeb).1
        } else {
            pc.d1(oc.ufin).1.reversed()
        };
        let der_hc = hc.d1(par_arc).1;
        // `cxx:1237-1247`. Note `prm2 = last` reads the *p-curve* limit here:
        // `BRep_Tool::CurveOnSurface` is what set `first` / `last` last.
        let (prm1, prm2, first_to_par) = if der_hc.dot(&der_pc) > 0.0 {
            (par_arc, last_pc, false)
        } else {
            (first_pc, par_arc, true)
        };
        // `Ct = BRep_Tool::Curve(CV[i].Arc(), first, last);
        //  Ct = new Geom_TrimmedCurve(Ct, prm1, prm2);` (`cxx:1249-1250`).
        let ct = BRepTool::edge_curve(&Edge(arc.clone()))
            .ok_or_else(|| "Failed to get curve of edge".to_string())?;
        let ct: Arc<dyn Curve> = Arc::new(GeomAdaptorCurve::with_range(ct, prm1, prm2));
        // `double toled = BRep_Tool::Tolerance(CV[i].Arc());` (`cxx:1251`).
        let toled = BRepTool::edge_tolerance(&Edge(arc.clone()));
        let indcurv = b
            .my_ds
            .change_ds()
            .add_curve(&mut DsCurve::with_curve(Some(ct), toled, false));
        // `indpoint = isfirst ? IndexFirstPointOnS1() : IndexLastPointOnS1()`
        // (`cxx:1257`).
        let indpoint = {
            let stripe = oc
                .stripe
                .read()
                .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
            if is_first {
                stripe.index_first_point_on_s1()
            } else {
                stripe.index_last_point_on_s1()
            }
        };
        let indvertex = b.my_ds.change_ds().add_shape(&oc.vtx);
        // `cxx:1259-1268`.
        let (interfp1, interfp2) = if first_to_par {
            (
                ch_fi3d_fil_point_in_ds(Orientation::Forward, indcurv, indvertex, prm1, true),
                ch_fi3d_fil_point_in_ds(Orientation::Reversed, indcurv, indpoint, prm2, false),
            )
        } else {
            (
                ch_fi3d_fil_point_in_ds(Orientation::Forward, indcurv, indpoint, prm1, false),
                ch_fi3d_fil_point_in_ds(Orientation::Reversed, indcurv, indvertex, prm2, true),
            )
        };
        {
            let ds = b.my_ds.change_ds();
            ds.change_curve_interferences(indcurv).push(interfp1);
            ds.change_curve_interferences(indcurv).push(interfp2);
        }
        // `indface = DStr.AddShape(Fv); Interfc = ChFi3d_FilCurveInDS(indcurv,
        //  indface, Hc, CV[i].Arc().Orientation());` (`cxx:1270-1272`).
        let indface = b.my_ds.change_ds().add_shape(&fv.0);
        let interfc = ch_fi3d_fil_curve_in_ds(
            indcurv,
            indface,
            Some(Arc::clone(&hc)),
            arc.orientation(),
        );
        b.my_ds
            .change_ds()
            .change_shape_interferences(&fv.0)
            .push(interfc);
        // `TopoDS_Edge aLocalEdge = CV[i].Arc(); aLocalEdge.Reverse();
        //  HcR = BRep_Tool::CurveOnSurface(aLocalEdge, Fv, first, last);`
        // (`cxx:1275-1278`).
        let mut local_edge = arc.clone();
        local_edge.reverse();
        let (hcr, _, _) = curve_on_surface_range(&Edge(local_edge.clone()), fv)
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let interfc = ch_fi3d_fil_curve_in_ds(
            indcurv,
            indface,
            Some(hcr),
            local_edge.orientation(),
        );
        b.my_ds
            .change_ds()
            .change_shape_interferences(&fv.0)
            .push(interfc);

        // `cxx:1283-1301`: the degenerated edge of `Fv` that holds `Vtx`.
        let mut degen = None;
        for e in explorer_shapes(&fv.0, ShapeType::Edge) {
            if !BRepTool::is_degenerated(&Edge(e.clone())) {
                continue;
            }
            let (vf, vl) = edge_vertices(&Edge(e.clone()));
            let at_vtx = vf.map(|v| is_same(&v.0, &oc.vtx)).unwrap_or(false)
                || vl.map(|v| is_same(&v.0, &oc.vtx)).unwrap_or(false);
            if at_vtx {
                degen = Some(e);
                break;
            }
        }
        // `cxx:1302-1324`.
        if let Some(edeg) = degen {
            let (cd, fd, _ld) = curve_on_surface_range(&Edge(edeg.clone()), fv)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            // `P2d = (FirstToPar) ? Hc->Value(first) : Hc->Value(last);`
            // (`cxx:1315`).
            let p2d = if first_to_par {
                hc.value(first_pc)
            } else {
                hc.value(last_pc)
            };
            // `Geom2dAPI_ProjectPointOnCurve Projector(P2d, Cd);
            //  par = Projector.LowerDistanceParameter();` (`cxx:1316-1317`).
            let par = project_point_on_curve_2d(cd.as_ref(), &p2d, CONFUSION)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?
                .parameter;
            let ideg = b.my_ds.change_ds().add_shape(&edeg);
            // `ori = (par < fd) ? TopAbs_FORWARD : TopAbs_REVERSED;`
            // (`cxx:1319`).
            let ori = if par < fd {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            let interfp1 = ch_fi3d_fil_point_in_ds(ori, ideg, indvertex, par, true);
            b.my_ds
                .change_ds()
                .change_shape_interferences(&edeg)
                .push(interfp1);
        }
    }
    Ok(())
}

/// `cxx:1332-1392` — the split (`intcouture`) DS write.
fn write_split(
    oc: &mut OneCorner,
    b: &mut ChFi3dBuilder,
    fv: &Face,
    ishape: i32,
    et: Orientation,
) -> Result<(), String> {
    let is_first = oc.is_first;
    let curv1 = oc
        .curv1
        .clone()
        .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
    let curv2 = oc
        .curv2
        .clone()
        .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;

    // `ComputeCurve2d(curv1, Fv, c2d1); InterFv = ChFi3d_FilCurveInDS(Icurv1,
    //  IShape, c2d1, Et); DStr.ChangeShapeInterferences(IShape).Append(...)`
    // (`cxx:1338-1342`).
    let c2d1 = compute_curve2d(&curv1, fv);
    let inter_fv = ch_fi3d_fil_curve_in_ds(oc.icurv1, ishape, c2d1, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .push(inter_fv);
    // `cxx:1343-1345`.
    let c2d2 = compute_curve2d(&curv2, fv);
    let inter_fv = ch_fi3d_fil_curve_in_ds(oc.icurv2, ishape, c2d2, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .push(inter_fv);

    // `if (Fd->Orientation() == Fv.Orientation()) Et = TopAbs::Reverse(Et);`
    // (`cxx:1347-1350`).
    let mut et = et;
    if oc.fd.orientation() == fv.0.orientation() {
        et = et.reversed();
    }
    // `c2d1 = new Geom2d_TrimmedCurve(Ps, Udeb, par2); InterFv =
    //  ChFi3d_FilCurveInDS(Icurv1, Isurf, c2d1, Et);
    //  DStr.ChangeSurfaceInterferences(Isurf).Append(InterFv);`
    // (`cxx:1351-1353`).
    let isurf = oc.fd.surf();
    let ps = oc
        .ps
        .clone()
        .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
    let c2d1 = Arc::new(Geom2dTrimmedCurve::new(Arc::clone(&ps), oc.udeb, oc.par2))
        as Arc<dyn Curve2d>;
    let inter = ch_fi3d_fil_curve_in_ds(oc.icurv1, isurf, Some(c2d1), et);
    b.my_ds
        .change_ds()
        .change_surface_interferences(isurf)
        .push(inter);
    // `c2d2 = new Geom2d_TrimmedCurve(Ps, par2, Ufin); ...` (`cxx:1354-1356`).
    let c2d2 = Arc::new(Geom2dTrimmedCurve::new(ps, oc.par2, oc.ufin)) as Arc<dyn Curve2d>;
    let inter = ch_fi3d_fil_curve_in_ds(oc.icurv2, isurf, Some(c2d2), et);
    b.my_ds
        .change_ds()
        .change_surface_interferences(isurf)
        .push(inter);

    // `cxx:1358-1377`: limitation of the sewing edge.
    let edgecouture = oc
        .edgecouture
        .clone()
        .ok_or_else(|| "OneCorner : no sewing edge".to_string())?;
    let iarc = b.my_ds.change_ds().add_shape(&edgecouture.0);
    let vdeb = first_vertex(&edgecouture.0)
        .ok_or_else(|| "OneCorner : sewing edge has no first vertex".to_string())?;
    let vfin = last_vertex(&edgecouture.0)
        .ok_or_else(|| "OneCorner : sewing edge has no last vertex".to_string())?;
    // `pard = BRep_Tool::Parameter(Vdeb, edgecouture);` and `parf`
    // (`cxx:1367-1368`).
    let pard = BRepTool::parameter_of_vertex(&edgecouture.0, &vdeb);
    let parf = BRepTool::parameter_of_vertex(&edgecouture.0, &vfin);
    // `cxx:1369-1375`.
    let ori = if (oc.par1 - pard).abs() < (parf - oc.par1).abs() {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let inter = ch_fi3d_fil_point_in_ds(ori, iarc, oc.indpt, oc.par1, false);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&edgecouture.0)
        .push(inter);

    // `cxx:1379-1391`: `stripe->InDS(isfirst);` then the two indexes, then the
    // four `CurvePointInterference`s.
    let (ind1, ind2) = {
        let mut stripe = oc
            .stripe
            .write()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        stripe.in_ds(is_first, 1);
        (stripe.index_point(is_first, 1), stripe.index_point(is_first, 2))
    };
    let i1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, oc.icurv1, ind1, oc.udeb, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(oc.icurv1)
        .push(i1);
    let i2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, oc.icurv1, oc.indpt, oc.par2, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(oc.icurv1)
        .push(i2);
    let i3 = ch_fi3d_fil_point_in_ds(Orientation::Forward, oc.icurv2, oc.indpt, oc.par2, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(oc.icurv2)
        .push(i3);
    let i4 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, oc.icurv2, ind2, oc.ufin, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(oc.icurv2)
        .push(i4);
    Ok(())
}

/// `cxx:1066-1393`.
pub(super) fn write(b: &mut ChFi3dBuilder, oc: &mut OneCorner) -> Result<(), String> {
    let is_first = oc.is_first;
    let fv = oc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let cv1 = oc.fd.vertex(is_first, 1).clone();
    let cv2 = oc.fd.vertex(is_first, 2).clone();

    // `cxx:1066-1092`.
    let (ishape, et) = face_index_and_orientation(oc, b, &fv)?;
    // `stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);`
    // and `CV2` on 2 (`cxx:1100-1101`).
    let ind_cv1 = ch_fi3d_index_point_in_ds(&cv1, b.my_ds.change_ds());
    let ind_cv2 = ch_fi3d_index_point_in_ds(&cv2, b.my_ds.change_ds());
    {
        let mut stripe = oc
            .stripe
            .write()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        stripe.set_index_point(ind_cv1, is_first, 1);
        stripe.set_index_point(ind_cv2, is_first, 2);
    }

    if !oc.intcouture {
        // `TopOpeBRepDS_Curve Tc(Cc, tolreached); ICurve = DStr.AddCurve(Tc);`
        // (`cxx:1110-1112`).
        let cc = oc
            .cc
            .clone()
            .ok_or_else(|| "OneCorner : echec calcul intersection".to_string())?;
        let icurve = b
            .my_ds
            .change_ds()
            .add_curve(&mut DsCurve::with_curve(Some(cc), oc.tolreached, false));
        // `Interfc = ChFi3d_FilCurveInDS(ICurve, IShape, Pc, Et);`
        // (`cxx:1113-1114`).
        let pc = oc
            .pc
            .clone()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let interfc = ch_fi3d_fil_curve_in_ds(icurve, ishape, Some(Arc::clone(&pc)), et);

        // `cxx:1116-1204` — OCC119.
        check_cork_overlap(oc, b, ishape, &fv, &pc)?;
        // `DStr.ChangeShapeInterferences(IShape).Append(Interfc);`
        // (`cxx:1205`).
        b.my_ds
            .change_ds()
            .change_shape_interferences(&fv.0)
            .push(interfc);

        // `cxx:1207-1327`.
        write_pseudo_seam(oc, b, &fv, &pc)?;

        // `stripe->ChangePCurve(isfirst) = Ps; stripe->SetCurve(ICurve,
        //  isfirst); stripe->SetParameters(isfirst, Udeb, Ufin);`
        // (`cxx:1328-1329`).
        let mut stripe = oc
            .stripe
            .write()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        *stripe.change_p_curve(is_first) = oc.ps.clone();
        stripe.set_curve(icurve, is_first);
        stripe.set_parameters(is_first, oc.udeb, oc.ufin);
    } else {
        // `cxx:1332-1392`.
        write_split(oc, b, &fv, ishape, et)?;
    }
    Ok(())
}
