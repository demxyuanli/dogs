//! `PerformIntersectionAtEnd`'s edge limitations —
//! `ChFi3d_Builder_C1.cxx:3262-3532`.
//!
//! The block is entered for every face that is not the last one, and for the
//! last one as well when the end is a single-intersection corner
//! (`nb != nbface || oneintersection1 || oneintersection2`). Its two arms:
//!
//! - `!extend && !(oneintersection1 || oneintersection2)` (`cxx:3320-3327`):
//!   the plain limitation of `Edge[nb]` by the parameter `paredge2` already
//!   computed by [`super::edge_cut`],
//! - the `else` (`cxx:3328-3531`): a new curve `Ct` on the face is built from
//!   the closest point of `Edge[nb]`'s curve (or of `edgesau` in the
//!   single-intersection case) to `Vtx`, then stored with its two p-curves —
//!   on `Sfacemoins1` and, in the single-intersection case, on `Sface` or
//!   `facesau`.
//!
//! Translation notes:
//! - `Geom_TrimmedCurve(csau, par1, par2)` is
//!   [`GeomAdaptorCurve::with_range`]: the port's trimmed 3D curve reports the
//!   basis parameters, so `FirstParameter()` / `D0(FirstParameter())` are
//!   `par1` / `csau(par1)` exactly as in OCCT, and
//!   `ChFi3d_EvalTolReached` samples `Ct` and its p-curves at the same
//!   parameters,
//! - `GeomAdaptor_Curve cad` is [`GeomAdaptorCurve`]; the three-argument
//!   `Extrema_ExtPC(P, C, Tol)` constructor uses the adaptor's whole range
//!   (`Extrema_GGExtPC.hxx:99-106`) and its third argument is the
//!   function-value tolerance `theTolF` (`hxx:84-92`), which OCCT passes
//!   `tolpt` for, so the port keeps `Point(1)` — the first solution of the
//!   interval sweep, not the closest one,
//! - `GeomLib::ExtendCurveToPoint(handle<Geom_BoundedCurve>&, ...)`
//!   (`cxx:3352`, `:3363`) is [`extend_curve_to_point`], which needs the
//!   B-spline OCCT's `GeomConvert_CompCurveToBSplineCurve` builds inside. The
//!   extension only runs when the curve really descends from
//!   `Geom_BoundedCurve` ([`is_bounded_curve`]), because OCCT's `down_cast`
//!   answers a null handle for `Geom_Line` / `Geom_Parabola` /
//!   `Geom_Hyperbola` and then skips the whole block,
//! - `ComputeCurve2d(Ct, faceprol[0], C2dint1)` (`cxx:3434`, `:3504`) is
//!   [`crate::chfi3d::builder_c1::compute_curve2d`], wired in T-336 once
//!   `BRepAlgo_NormalProjection` was ported; the `extend` arms still leave
//!   `C2dint1` / `C2dint2` null when the projection produces no edge.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_core::convert::ParameterisationType;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::extrema_pc::extrema_ext_pc_range;
use occt_geom::geom_lib_extend_curve::extend_curve_to_point;
use occt_geom::Curve as GeomCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_eval_tol_reached, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds,
};
use crate::chfi3d::builder_c1::{compute_curve2d, contain_v};
use crate::chfi_ds::CommonPoint;
use crate::shape::{Edge, Vertex};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::Curve as DsCurve;
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

use super::super::ChFi3dBuilder;
use super::corner_faces::face_at;
use super::face_loop::{project_2d, FaceLoop};
use super::EndIntersection;

/// `occ::down_cast<Geom_BoundedCurve>(C1)` for the curve types the port
/// stores: `Geom_Line`, `Geom_Parabola` and `Geom_Hyperbola` are the three
/// adaptor types that do not descend from `Geom_BoundedCurve`.
fn is_bounded_curve(c: &dyn GeomCurve) -> bool {
    c.gp_line().is_none() && c.gp_parabola().is_none() && c.gp_hyperbola().is_none()
}

/// The `ori` choice of `cxx:3264-3318`, which throws
/// `Standard_Failure("IntersectionAtEnd : pb orientation")` when neither
/// vertex of `Edge[nb]` belongs to `F1` or `F2`.
fn end_orientation(
    st: &EndIntersection,
    e_nb: &Edge,
    nface: usize,
) -> Result<Orientation, String> {
    let (v1, v2) = edge_vertices(e_nb);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return Err("IntersectionAtEnd : pb orientation".to_string());
    };
    // `if (nface == 3)` (`cxx:3264-3309`).
    if nface == 3 {
        let on_f1 = |v: &Vertex| st.f1.as_ref().is_some_and(|f| contain_v(f, &v.0));
        let on_f2 = |v: &Vertex| st.f2.as_ref().is_some_and(|f| contain_v(f, &v.0));
        let mut ori = if on_f1(&v1) || on_f2(&v1) {
            Orientation::Forward
        } else if on_f1(&v2) || on_f2(&v2) {
            Orientation::Reversed
        } else {
            return Err("IntersectionAtEnd : pb orientation".to_string());
        };
        let vtx = BRepTool::vertex_point(&Vertex(st.vtx.clone()));
        // `if (containV(F1, V1) && containV(F1, V2))` (`cxx:3282-3295`).
        if on_f1(&v1) && on_f1(&v2) {
            let d1 = BRepTool::vertex_point(&v1).distance(&vtx);
            let d2 = BRepTool::vertex_point(&v2).distance(&vtx);
            ori = if d1 < d2 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
        }
        // `if (containV(F2, V1) && containV(F2, V2))` (`cxx:3296-3308`).
        if on_f2(&v1) && on_f2(&v2) {
            let d1 = BRepTool::vertex_point(&v1).distance(&vtx);
            let d2 = BRepTool::vertex_point(&v2).distance(&vtx);
            ori = if d1 < d2 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
        }
        return Ok(ori);
    }
    // `if (TopExp::FirstVertex(Edge[nb]).IsSame(Vtx)) ori = FORWARD; else ori
    //  = REVERSED;` (`cxx:3311-3318`).
    Ok(if is_same(&v1.0, &st.vtx) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}

/// `cxx:3262-3532`.
pub(super) fn limit(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &mut FaceLoop,
    cv1: &CommonPoint,
    cv2: &CommonPoint,
    e_nb: &Edge,
    nb: i32,
) -> Result<(), String> {
    // `if (nb != nbface || oneintersection1 || oneintersection2)`
    // (`cxx:3262`).
    let one = st.one_intersection1 || st.one_intersection2;
    if nb == st.nbface && !one {
        return Ok(());
    }
    // `ori` (`cxx:3237`).
    lp.ori = end_orientation(st, e_nb, st.nface)?;

    // `if (!extend && !(oneintersection1 || oneintersection2))`
    // (`cxx:3320-3327`).
    if !lp.extend && !one {
        let iarc2 = b.my_ds.change_ds().add_shape(&e_nb.0);
        // `ChFi3d_FilPointInDS(ori, Iarc2, indpoint2, paredge2)` — the
        // default `IsVertex = false`.
        lp.interf_edge[nb as usize] = Some(ch_fi3d_fil_point_in_ds(
            lp.ori,
            iarc2,
            lp.ind_point2,
            lp.paredge2,
            false,
        ));
        return Ok(());
    }

    // `if (!(oneintersection1 || oneintersection2)) proledge[nb] = true;`
    // (`cxx:3328-3331`).
    if !one {
        lp.prol_edge[nb as usize] = 1;
    }
    // `pext = BRep_Tool::Pnt(Vtx);` (`cxx:3334`).
    let mut pext = BRepTool::vertex_point(&Vertex(st.vtx.clone()));
    // `GeomAdaptor_Curve cad; handle<Geom_Curve> csau;` (`cxx:3335-3336`).
    let csau: Arc<dyn GeomCurve>;
    if !one {
        // `cad.Load(cint); csau = cint;` (`cxx:3339-3341`) — `cint` is the
        // `GeomInt_IntSS` line of the `extend` arm of `cxx:2868-2886`.
        csau = lp
            .cint
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: cint is null".to_string())?;
    } else {
        // `csau = BRep_Tool::Curve(edgesau, ubid, vbid);` (`cxx:3346`).
        let edgesau = st
            .edge_sau
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: edgesau is null".to_string())?;
        let mut cur = BRepTool::edge_curve(&edgesau)
            .ok_or_else(|| "IntersectionAtEnd: edgesau has no curve".to_string())?;
        let bounded = is_bounded_curve(cur.as_ref());
        // `if (oneintersection1 && extend) { Pl = C1->Value(C1->
        //  LastParameter()); sens = Pl.Distance(pext) < tolpt;
        //  GeomLib::ExtendCurveToPoint(C1, CV1.Point(), 1, sens); csau = C1; }`
        // (`cxx:3347-3370`).
        if st.one_intersection1 && lp.extend && bounded {
            let pl = cur.d0(cur.last_parameter());
            let sens = pl.distance(&pext) < st.tolpt;
            cur = extend_with_point(cur, cv1.point(), sens)?;
        } else if st.one_intersection2 && lp.extend && bounded {
            let pl = cur.d0(cur.last_parameter());
            let sens = pl.distance(&pext) < st.tolpt;
            cur = extend_with_point(cur, cv2.point(), sens)?;
        }
        csau = cur;
    }
    // `cad.Load(csau); Extrema_ExtPC ext(pext, cad, tolpt);
    //  par1 = ext.Point(1).Parameter(); ParVtx = par1;` (`cxx:3371-3377`).
    let cad = GeomAdaptorCurve::of_curve(Arc::clone(&csau));
    let par1_ext = ext_point_parameter(&cad, &pext)?;
    // `cxx:3378-3398`.
    let (par2_ext, ind) = if one {
        if st.one_intersection2 {
            pext = *cv2.point();
            (ext_point_parameter(&cad, &pext)?, lp.ind_point2)
        } else {
            pext = *cv1.point();
            (ext_point_parameter(&cad, &pext)?, lp.ind_point1)
        }
    } else {
        (lp.paredge2, lp.ind_point2)
    };
    // `if (par1 > par2) { indp1 = ind; indp2 = DStr.AddShape(Vtx); vtx2 =
    //  true; par = par1; par1 = par2; par2 = par; } else { indp1 =
    //  DStr.AddShape(Vtx); indp2 = ind; vtx1 = true; }` (`cxx:3399-3415`).
    let ind_vtx = b.my_ds.change_ds().add_shape(&st.vtx);
    let (mut par1, mut par2, indp1, indp2, vtx1, vtx2) = if par1_ext > par2_ext {
        (par2_ext, par1_ext, ind, ind_vtx, false, true)
    } else {
        (par1_ext, par2_ext, ind_vtx, ind, true, false)
    };
    // `Ct = new Geom_TrimmedCurve(csau, par1, par2);` (`cxx:3416`).
    let mut ct: Arc<dyn GeomCurve> =
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&csau), par1, par2));
    // `Cc->D0(Cc->FirstParameter(), P1); Cc->D0(Cc->LastParameter(), P2);
    //  Ct->D0(Ct->FirstParameter(), P3); Ct->D0(Ct->LastParameter(), P4);
    //  if (P2.Distance(P3) < tolpt || P1.Distance(P4) < tolpt) orient = orien;
    //  else orient = TopAbs::Reverse(orien);` (`cxx:3417-3427`).
    let cc = lp
        .cc
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Cc is null".to_string())?;
    let pa = cc.d0(cc.first_parameter());
    let pb = cc.d0(cc.last_parameter());
    let pc_ = ct.d0(par1);
    let pd = ct.d0(par2);
    let mut orient = if pb.distance(&pc_) < st.tolpt || pa.distance(&pd) < st.tolpt {
        lp.orien
    } else {
        lp.orien.reversed()
    };

    // `int indice;` (`cxx:3237` region) — the index the two `Interfc` are
    // appended to.
    let mut indice;
    if one {
        // `indice = DStr.AddShape(Face[0]); if (extend) {
        //  DStr.SetNewSurface(Face[0], Sfacemoins1); ComputeCurve2d(Ct,
        //  faceprol[0], C2dint1); } else { TopoDS_Edge aLocalEdge = edgesau;
        //  if (edgesau.Orientation() != orient) aLocalEdge.Reverse();
        //  C2dint1 = CurveOnSurface(aLocalEdge, Face[0], ubid, vbid); }`
        // (`cxx:3429-3445`).
        let face0 = face_at(st, 0, 3430)?;
        indice = b.my_ds.change_ds().add_shape(&face0.0);
        if lp.extend {
            let surf = lp
                .sface_moins1
                .clone()
                .ok_or_else(|| "IntersectionAtEnd: Sfacemoins1 is null".to_string())?;
            b.my_ds.change_ds().set_new_surface(&face0.0, surf);
            lp.c2d_int1 = lp
                .face_prol
                .first()
                .cloned()
                .flatten()
                .and_then(|f| compute_curve2d(&ct, &f));
        } else {
            let edgesau = st
                .edge_sau
                .clone()
                .ok_or_else(|| "IntersectionAtEnd: edgesau is null".to_string())?;
            let local = local_edge(&edgesau, orient);
            lp.c2d_int1 = curve_on_surface_range(&local, &face0).map(|(c, _, _)| c);
        }
    } else {
        // `indice = DStr.AddShape(Face[nb - 1]);
        //  DStr.SetNewSurface(Face[nb - 1], Sfacemoins1);` (`cxx:3447-3450`).
        let face_nb1 = face_at(st, (nb - 1) as usize, 3447)?;
        indice = b.my_ds.change_ds().add_shape(&face_nb1.0);
        let surf = lp
            .sface_moins1
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: Sfacemoins1 is null".to_string())?;
        b.my_ds.change_ds().set_new_surface(&face_nb1.0, surf);
    }

    // `if (cad.IsPeriodic() && !C2dint1.IsNull()) { P2d =
    //  BRep_Tool::Parameters(Vtx, Face[0]); Geom2dAPI_ProjectPointOnCurve
    //  Projector(P2d, C2dint1); par = Projector.LowerDistanceParameter();
    //  shift = par - ParVtx; if (|shift| > Confusion) { par1 += shift;
    //  par2 += shift; } }` (`cxx:3451-3464`).
    if csau.is_periodic() {
        if let Some(c2d1) = lp.c2d_int1.clone() {
            let face0 = face_at(st, 0, 3453)?;
            let p2d = BRepTool::parameters_of_vertex_on_face(&st.vtx, &face0)?;
            let par = project_2d(&c2d1, &p2d)?;
            let shift = par - par1_ext;
            if shift.abs() > CONFUSION {
                par1 += shift;
                par2 += shift;
            }
        }
    }
    // `Ct = new Geom_TrimmedCurve(csau, par1, par2);` (`cxx:3467`).
    ct = Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&csau), par1, par2));
    // `if (oneintersection1 || oneintersection2) tolex = 10 *
    //  BRep_Tool::Tolerance(edgesau);` (`cxx:3468-3471`).
    if one {
        let edgesau = st
            .edge_sau
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: edgesau is null".to_string())?;
        lp.tolex = 10.0 * BRepTool::edge_tolerance(&edgesau);
    }
    // `if (extend) { H1 = new GeomAdaptor_Surface(Sfacemoins1); if
    //  (Sface.IsNull()) tolex = max(tolex, ChFi3d_EvalTolReached(H1, C2dint1,
    //  H1, C2dint1, Ct)); else { H2 = new GeomAdaptor_Surface(Sface); tolex =
    //  max(tolex, ChFi3d_EvalTolReached(H1, C2dint1, H2, C2dint2, Ct)); } }`
    // (`cxx:3472-3485`).
    if lp.extend {
        if let Some(c2d1) = lp.c2d_int1.clone() {
            let h1 = lp
                .sface_moins1
                .clone()
                .ok_or_else(|| "IntersectionAtEnd: Sfacemoins1 is null".to_string())?;
            let tol = match lp.sface.clone() {
                None => ch_fi3d_eval_tol_reached(
                    h1.as_ref(),
                    c2d1.as_ref(),
                    h1.as_ref(),
                    c2d1.as_ref(),
                    ct.as_ref(),
                ),
                Some(h2) => {
                    let c2d2 = lp
                        .c2d_int2
                        .clone()
                        .ok_or_else(|| "IntersectionAtEnd: C2dint2 is null".to_string())?;
                    ch_fi3d_eval_tol_reached(
                        h1.as_ref(),
                        c2d1.as_ref(),
                        h2.as_ref(),
                        c2d2.as_ref(),
                        ct.as_ref(),
                    )
                }
            };
            lp.tolex = lp.tolex.max(tol);
        }
    }

    // `TopOpeBRepDS_Curve tcurv(Ct, tolex); int indcurv =
    //  DStr.AddCurve(tcurv);` (`cxx:3486-3488`).
    let indcurv = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&ct)),
        lp.tolex,
        false,
    ));
    // `Interfp1 = ChFi3d_FilPointInDS(FORWARD, indcurv, indp1, par1, vtx1);
    //  Interfp2 = ... (REVERSED, indp2, par2, vtx2);` and the two appends
    // (`cxx:3489-3493`).
    let interfp1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, indcurv, indp1, par1, vtx1);
    let interfp2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, indcurv, indp2, par2, vtx2);
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(indcurv);
        list.push(interfp1);
        list.push(interfp2);
    }
    // `Interfc = ChFi3d_FilCurveInDS(indcurv, indice, C2dint1, orient);
    //  DStr.ChangeShapeInterferences(indice).Append(Interfc);`
    // (`cxx:3494-3497`).
    let interfc = ch_fi3d_fil_curve_in_ds(indcurv, indice, lp.c2d_int1.clone(), orient);
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_shape_interferences_at(indice);
        list.push(interfc);
    }

    if one {
        // `indice = DStr.AddShape(facesau); if (facesau.Orientation() ==
        //  Face[0].Orientation()) orient = TopAbs::Reverse(orient); if
        //  (extend) ComputeCurve2d(Ct, faceprol[1], C2dint2); else {
        //  aLocalEdge = edgesau; if (edgesau.Orientation() != orient)
        //  aLocalEdge.Reverse(); C2dint2 = CurveOnSurface(aLocalEdge,
        //  facesau, ubid, vbid); }` (`cxx:3498-3515`).
        let facesau = st
            .face_sau
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: facesau is null".to_string())?;
        indice = b.my_ds.change_ds().add_shape(&facesau.0);
        let face0 = face_at(st, 0, 3499)?;
        if facesau.0.orientation() == face0.0.orientation() {
            orient = orient.reversed();
        }
        if lp.extend {
            lp.c2d_int2 = lp
                .face_prol
                .get(1)
                .cloned()
                .flatten()
                .and_then(|f| compute_curve2d(&ct, &f));
        } else {
            let edgesau = st
                .edge_sau
                .clone()
                .ok_or_else(|| "IntersectionAtEnd: edgesau is null".to_string())?;
            let local = local_edge(&edgesau, orient);
            lp.c2d_int2 = curve_on_surface_range(&local, &facesau).map(|(c, _, _)| c);
        }
        // `if (!bordlibre) { Interfc = ChFi3d_FilCurveInDS(indcurv, indice,
        //  C2dint2, orient); DStr.ChangeShapeInterferences(indice).Append(
        //  Interfc); }` (`cxx:3526-3529`).
        if !st.bord_libre {
            let interfc = ch_fi3d_fil_curve_in_ds(indcurv, indice, lp.c2d_int2.clone(), orient);
            let ds = b.my_ds.change_ds();
            let list = ds.change_shape_interferences_at(indice);
            list.push(interfc);
        }
    } else {
        // `indice = DStr.AddShape(Face[nb]); DStr.SetNewSurface(Face[nb],
        //  Sface); if (Face[nb].Orientation() == Face[nb - 1].Orientation())
        //  orient = TopAbs::Reverse(orient);` (`cxx:3517-3524`).
        let face_nb = face_at(st, nb as usize, 3518)?;
        indice = b.my_ds.change_ds().add_shape(&face_nb.0);
        let sface = lp
            .sface
            .clone()
            .ok_or_else(|| "IntersectionAtEnd: Sface is null".to_string())?;
        b.my_ds.change_ds().set_new_surface(&face_nb.0, sface);
        let face_nb1 = face_at(st, (nb - 1) as usize, 3521)?;
        if face_nb.0.orientation() == face_nb1.0.orientation() {
            orient = orient.reversed();
        }
        // `if (!bordlibre) { ... }` (`cxx:3526-3529`).
        if !st.bord_libre {
            let interfc = ch_fi3d_fil_curve_in_ds(indcurv, indice, lp.c2d_int2.clone(), orient);
            let ds = b.my_ds.change_ds();
            let list = ds.change_shape_interferences_at(indice);
            list.push(interfc);
        }
    }
    Ok(())
}

/// `TopoDS_Edge aLocalEdge = edgesau; if (edgesau.Orientation() != orient)
/// aLocalEdge.Reverse();`.
fn local_edge(edgesau: &Edge, orient: Orientation) -> Edge {
    if edgesau.0.orientation() != orient {
        let mut e = edgesau.0.clone();
        e.reverse();
        Edge(e)
    } else {
        edgesau.clone()
    }
}

/// `ExtendCurveToPoint(C1, P, 1, sens)` on the B-spline OCCT's `Concat`
/// builds, keeping the port's `Arc<dyn Curve>` view of the result.
fn extend_with_point(
    cur: Arc<dyn GeomCurve>,
    point: &GpPnt,
    after: bool,
) -> Result<Arc<dyn GeomCurve>, String> {
    let mut bs = curve_to_bspline_curve(cur.as_ref(), ParameterisationType::QuasiAngular)
        .map_err(|_| "IntersectionAtEnd: ExtendCurveToPoint on a non-bounded curve".to_string())?;
    extend_curve_to_point(&mut bs, point, 1, after)?;
    Ok(Arc::new(bs))
}

/// `Extrema_ExtPC(P, C, Tol)` then `Point(1).Parameter()`.
fn ext_point_parameter(cad: &GeomAdaptorCurve, p: &GpPnt) -> Result<f64, String> {
    let sols = extrema_ext_pc_range(cad, p, cad.first_parameter(), cad.last_parameter());
    sols.first()
        .map(|s| s.u)
        .ok_or_else(|| "IntersectionAtEnd: Extrema_ExtPC found no solution".to_string())
}
