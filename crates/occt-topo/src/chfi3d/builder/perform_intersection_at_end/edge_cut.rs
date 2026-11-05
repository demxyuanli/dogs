//! `PerformIntersectionAtEnd`'s end-edge cutting of the fillet surface —
//! `ChFi3d_Builder_C1.cxx:2767-2952`, `:2953-2992`, `:2998-3052`.
//!
//! [`cut`] is the `if (nb != nbface || isOnSame2)` block: the curve of `E2`
//! against the fillet surface `HGs`
//! (`IntCurveSurface_HInter::Perform`), with the two OCCT fallbacks — the
//! `GeomLib::ExtendSurfByLength` extension of the fillet surface
//! (`cxx:2798-2841`) and the `GeomInt_IntSS` surface/surface lines
//! (`cxx:2842-2889`) — then the nearest intersection point to `Vtx` and its
//! p-curves (`cxx:2890-2951`).
//!
//! [`last_edge_pcurve`] is the `else` arm (`cxx:2953-2992`) that answers
//! `Hc2` and `pfac2` when `nb == nbface` and the end is not `OnSame`.
//!
//! [`same_end_update`] is the `isOnSame2` close (`cxx:2998-3052`): the
//! `IntersUpdateOnSame` on `Fd->ChangeInterferenceOnS2()`, the `paredge2`
//! re-projection when the intersection at the end failed, the point in the DS
//! and `CV2.SetArc`.
//!
//! Translation notes:
//! - `Geom_TrimmedCurve(C, Ubid, Vbid)` only feeds `IsPeriodic`, `Period` and
//!   `BasisCurve()->First/LastParameter`, all of which the basis curve answers
//!   as well, so the trim is not materialised; the ctor's `U1 > U2` rejection
//!   is kept as an error.
//! - `Cad` aliases `HC` (`cxx:2796`), so `Cad.Load(cint)` re-points the adaptor
//!   the following `inters.Perform(HC, HGs)` uses; the loop re-runs
//!   `perform_curve_surface` on `cint` with its own range.
//! - `Gs.Load(S1)` (`cxx:2810`) is written and never read in OCCT either: both
//!   `Perform` calls use `HGs`, which points at the surface
//!   `ExtendSurfByLength` mutated in place. The port rebuilds `st.hgs` on the
//!   extended surface over the same parameter box, which is the state
//!   `HGs` has after the in-place mutation.
//! - `faceprol` slots are `BRep_Builder::MakeFace(faceprol[i], S, L, Tol)`.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::extrema_surf::point_surface_extrema;
use occt_geom::geom_lib_extend_surf::extend_surf_by_length;
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_c1::{ch_fi3d_extend_surface, inter_update_on_same};
use crate::chfi_ds::{CommonPoint, FaceInterference};
use crate::geom_int::IntSS;
use crate::intcurvesurface::{perform_curve_surface, HInterResult};
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::edge_vertices;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;
use super::corner_faces::{edge_at, face_at};
use super::face_loop::{add_point, fi_param_on_surf, make_prol_face, project_2d, FaceLoop};
use super::{fd_of, EndIntersection, EndIntersectionExit};

/// The `HGs` parameter box as the `uv_range` of `perform_curve_surface`.
fn uv_of(a: &SurfaceAdaptor) -> (f64, f64, f64, f64) {
    (
        a.first_u_parameter(),
        a.first_v_parameter(),
        a.last_u_parameter(),
        a.last_v_parameter(),
    )
}

/// `GeomAdaptor_Surface::Load(S)` + `Perform(HC, HGs)`: the curve of `e2` (or
/// `cint`) against the fillet surface over the adaptor's box.
fn perform_on_curve(
    curve: &dyn occt_geom::Curve,
    cu_range: (f64, f64),
    hgs: &SurfaceAdaptor,
) -> Option<HInterResult> {
    let surface = hgs.surface()?;
    perform_curve_surface(curve, surface.as_ref(), cu_range, uv_of(hgs)).ok()
}

/// `TopExp::Vertices(E, v1, v2)` then `Vtx.IsSame(v1) ? FORWARD : REVERSED`
/// (`cxx:2711-2713`, `:3048-3050`).
pub(super) fn arc_orientation(e: &Edge, vtx: &crate::shape::TopoShape) -> Orientation {
    let (v1, _v2) = edge_vertices(e);
    if v1.map(|v| crate::topo_tools::is_same(&v.0, vtx)).unwrap_or(false) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

/// `cxx:2767-2952`. `Ok(Some(exit))` is the `PerformMoreThreeCorner(Index, 1)`
/// of `cxx:2901-2902`.
#[allow(clippy::too_many_arguments)]
pub(super) fn cut(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &mut FaceLoop,
    e2: &Edge,
    f: &Face,
    sface_moins1: &mut Arc<dyn Surface>,
    nb: i32,
) -> Result<Option<EndIntersectionExit>, String> {
    let is_first = st.is_first;
    // `C = BRep_Tool::Curve(E2, Ubid, Vbid); Ctrim = new Geom_TrimmedCurve(C,
    //  Ubid, Vbid); Utrim = Ctrim->BasisCurve()->FirstParameter(); Vtrim = ...
    //  (`cxx:2771-2776`).
    let c = BRepTool::edge_curve(e2)
        .ok_or_else(|| "IntersectionAtEnd: E2 has no curve".to_string())?;
    let (ubid0, vbid0) = BRepTool::edge_parameters(e2);
    if ubid0 > vbid0 {
        // `Geom_TrimmedCurve` refuses `U1 > U2`.
        return Err("IntersectionAtEnd: E2 curve parameters are invalid".to_string());
    }
    let utrim = c.first_parameter();
    let vtrim = c.last_parameter();
    // `cxx:2777-2794`.
    let (ubid, vbid) = if c.is_periodic() {
        if ubid0 > c.period() {
            (0.5 * (utrim + vtrim), vtrim)
        } else {
            (utrim, 0.5 * (utrim + vtrim))
        }
    } else {
        (utrim, vtrim)
    };
    // `GeomAdaptor_Curve HC(C, Ubid, Vbid); inters.Perform(HC, HGs);`
    // (`cxx:2795-2797`).
    let mut inters = perform_on_curve(c.as_ref(), (ubid, vbid), &st.hgs);

    // `if (!prolface[nn] && (!inters.IsDone() || inters.NbPoints() == 0))`
    // (`cxx:2798`).
    let no_point = |i: &Option<HInterResult>| i.as_ref().map(|r| r.nb_points() == 0).unwrap_or(true);
    if lp.prol_face[super::NN] == 0 && no_point(&inters) {
        // `S1 = down_cast<Geom_BoundedSurface>(DStr.Surface(Fd->Surf()).Surface())`
        // (`cxx:2801-2802`). Every `Geom_Surface` the port implements descends
        // from `Geom_BoundedSurface`, so the down cast succeeds whenever the
        // DS surface is not null.
        let isurf = fd_of(st)?.surf();
        let s1 = {
            let ds = b.my_ds.change_ds();
            ds.surface(isurf).surface().cloned()
        };
        if let Some(s1) = s1 {
            // `length = 0.5 * std::max(Fi1Length, Fi2Length);
            //  GeomLib::ExtendSurfByLength(S1, length, 1, false, !isfirst);`
            // (`cxx:2805-2806`).
            let length = 0.5 * st.fi1_length.max(st.fi2_length);
            let extended = extend_surf_by_length(&s1, length, 1, false, !is_first)
                .unwrap_or_else(|| Arc::clone(&s1));
            // `prolface[nn] = 1;` (`cxx:2807`).
            lp.prol_face[super::NN] = 1;
            // `if (!stripe->IsInDS(!isfirst))` (`cxx:2808`).
            let in_ds = st
                .stripe
                .read()
                .map(|s| s.is_in_ds(!is_first))
                .unwrap_or(0);
            if in_ds == 0 {
                // `Gs.Load(S1); inters.Perform(HC, HGs);` (`cxx:2810-2811`) —
                // `HGs` sees the extended surface in place.
                let uv = uv_of(&st.hgs);
                st.hgs = SurfaceAdaptor::of_surface_with_range(
                    Arc::clone(&extended),
                    uv.0,
                    uv.2,
                    uv.1,
                    uv.3,
                );
                inters = perform_on_curve(c.as_ref(), (ubid, vbid), &st.hgs);
                // `if (inters.IsDone() && inters.NbPoints() != 0)`
                // (`cxx:2812`).
                if !no_point(&inters) {
                    // `Fd->ChangeSurf(DStr.AddSurface(TopOpeBRepDS_Surface(S1,
                    //  DStr.ChangeSurface(Isurf).Tolerance())));`
                    // (`cxx:2814-2815`).
                    let new_index = {
                        let ds = b.my_ds.change_ds();
                        let tol = ds.change_surface(isurf).tolerance();
                        ds.add_surface(&crate::topopebrep::ds::Surface::with_tolerance(
                            Arc::clone(&extended),
                            tol,
                        ))
                    };
                    if let Some(slot) = st.num.checked_sub(1).and_then(|i| st.seq.get_mut(i)) {
                        Arc::make_mut(slot).set_surf(new_index);
                    }
                    // `if (myEVIMap.IsBound(EdgeSpine)) { ... remove Isurf ...
                    //  Append(Fd->Surf()); } else { IndexList = [Fd->Surf()];
                    //  myEVIMap.Bind(EdgeSpine, IndexList); }` (`cxx:2817-2835`).
                    let list = b.change_evi(&st.edge_spine);
                    if let Some(pos) = list.iter().position(|&x| x == isurf) {
                        list.remove(pos);
                    }
                    list.push(new_index);
                    // `Isurf = Fd->Surf();` (`cxx:2837`).
                    st.isurf = new_index;
                }
            }
        }
    }

    // `if (!inters.IsDone() || inters.NbPoints() == 0)` (`cxx:2842`).
    if no_point(&inters) {
        // `cd = down_cast<Geom_BSplineCurve>(C); cd1 =
        //  down_cast<Geom_BezierCurve>(C); if (!cd.IsNull() || !cd1.IsNull())`
        // (`cxx:2844-2846`).
        let is_bspline = c.bspline_knots().is_some();
        let is_bezier = c.bezier_poles().is_some();
        if is_bspline || is_bezier {
            // `Sface = BRep_Tool::Surface(Face[nb]);` (`cxx:2849`).
            let face_nb = face_at(st, nb as usize, 2849)?;
            let mut sface = BRepTool::face_surface(&face_nb)
                .ok_or_else(|| "IntersectionAtEnd: Face[nb] has no surface".to_string())?;
            // `ChFi3d_ExtendSurface(Sface, prolface[nb]); tol =
            //  BRep_Tool::Tolerance(F); BRE.MakeFace(faceprol[nb], Sface,
            //  Face[nb].Location(), tol);` (`cxx:2850-2852`).
            sface = ch_fi3d_extend_surface(&sface, &mut lp.prol_face[nb as usize]);
            lp.sface = Some(Arc::clone(&sface));
            let tol = BRepTool::face_tolerance(f);
            lp.face_prol[nb as usize] = Some(make_prol_face(Arc::clone(&sface), tol));
            // `cxx:2853-2866`.
            if nb != 0 && lp.prol_face[(nb - 1) as usize] == 0 {
                *sface_moins1 =
                    ch_fi3d_extend_surface(sface_moins1, &mut lp.prol_face[(nb - 1) as usize]);
                if lp.prol_face[(nb - 1) as usize] != 0 {
                    let tol = BRepTool::face_tolerance(f);
                    lp.face_prol[(nb - 1) as usize] =
                        Some(make_prol_face(Arc::clone(sface_moins1), tol));
                }
            } else {
                // `int prol = 0; ChFi3d_ExtendSurface(Sfacemoins1, prol);` — the
                // extended surface is discarded (`cxx:2864-2865`).
                let mut prol = 0;
                let _ = ch_fi3d_extend_surface(sface_moins1, &mut prol);
            }
            // `GeomInt_IntSS InterSS(Sfacemoins1, Sface, 1.e-7, true, true,
            //  true);` (`cxx:2867`).
            let mut inter_ss = IntSS::new();
            inter_ss.perform(
                Arc::clone(sface_moins1),
                Arc::clone(&sface),
                1.0e-7,
                true,
                true,
                true,
            );
            // `cxx:2868-2887`.
            if inter_ss.is_done() {
                let mut trouve = false;
                let mut i = 1;
                while i <= inter_ss.nb_lines() && !trouve {
                    lp.extend = true;
                    let cint = inter_ss.line(i).map(|l| Arc::clone(&l.curve));
                    lp.cint = cint.clone();
                    lp.c2d_int1 = inter_ss.line_on_s1(i).cloned();
                    lp.c2d_int2 = inter_ss.line_on_s2(i).cloned();
                    if let Some(cint) = cint {
                        let range = (cint.first_parameter(), cint.last_parameter());
                        inters = perform_on_curve(cint.as_ref(), range, &st.hgs);
                        trouve = !no_point(&inters);
                    }
                    // `tolex = InterSS.TolReached3d();` (`cxx:2885`).
                    lp.tolex = inter_ss.tol_reached_3d();
                    i += 1;
                }
            }
        }
    }

    // `if (inters.IsDone())` (`cxx:2890`).
    let Some(result) = inters else {
        // `throw Standard_Failure("IntersectionAtEnd: pb intersection Face
        //  cb");` (`cxx:2950`).
        return Err("IntersectionAtEnd: pb intersection Face cb".to_string());
    };
    let nbp = result.nb_points();
    if nbp == 0 {
        // `if (nb == 0 || nb == nbface) intersOnSameFailed = true; else {
        //  PerformMoreThreeCorner(Index, 1); return; }` (`cxx:2893-2903`).
        if nb == 0 || nb == st.nbface {
            lp.inters_on_same_failed = true;
        } else {
            return Ok(Some(EndIntersectionExit::MoreThreeCorner(1)));
        }
    } else {
        // `gp_Pnt P = BRep_Tool::Pnt(Vtx); distmin = P.Distance(Point(1).Pnt());
        //  nbp = 1; for (i = 2; i <= NbPoints(); i++) ...` (`cxx:2907-2918`).
        let p = BRepTool::vertex_point(&Vertex(st.vtx.clone()));
        let mut distmin = p.distance(&result.point(0).pnt());
        let mut best = 1usize;
        for i in 1..nbp {
            let dist = p.distance(&result.point(i).pnt());
            if dist < distmin {
                distmin = dist;
                best = i + 1;
            }
        }
        // `gp_Pnt2d pt2d(Point(nbp).U(), Point(nbp).V()); pfil2 = pt2d;
        //  paredge2 = Point(nbp).W();` (`cxx:2919-2921`).
        let pt = result.point(best - 1);
        lp.pfil2 = GpPnt2d::new(pt.u(), pt.v());
        lp.paredge2 = pt.param();
        if !lp.extend {
            // `cfacemoins1 = CurveOnSurface(E2, F, u2, v2); ... cface =
            //  CurveOnSurface(E2, Face[nb], u2, v2); cfacemoins1->D0(paredge2,
            //  pfac2); cface->D0(paredge2, pint);` (`cxx:2924-2935`).
            let (cface_moins1, _, _) = curve_on_surface_range(e2, f)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            let face_nb = face_at(st, nb as usize, 2929)?;
            let (cface, _, _) = curve_on_surface_range(e2, &face_nb)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            lp.pfac2 = cface_moins1.value(lp.paredge2);
            lp.pint = cface.value(lp.paredge2);
        } else if lp.c2d_int1.is_none() || lp.c2d_int2.is_none() {
            // `throw Standard_ConstructionError("Failed to get p-curve of
            //  edge");` (`cxx:2939`).
            return Err("Failed to get p-curve of edge".to_string());
        } else {
            // `C2dint1->D0(paredge2, pfac2); C2dint2->D0(paredge2, pint);`
            // (`cxx:2943-2944`).
            let c2d1 = lp.c2d_int1.as_ref().expect("checked above");
            let c2d2 = lp.c2d_int2.as_ref().expect("checked above");
            lp.pfac2 = GpPnt2d::new(0.0, 0.0);
            let v1 = c2d1.value(lp.paredge2);
            let v2 = c2d2.value(lp.paredge2);
            lp.pfac2 = v1;
            lp.pint = v2;
        }
    }
    Ok(None)
}

/// `cxx:2953-2992`: the `nb == nbface && !isOnSame2` arm.
pub(super) fn last_edge_pcurve(
    st: &EndIntersection,
    lp: &mut FaceLoop,
    cv2: &mut CommonPoint,
    fi2: &FaceInterference,
    e2: &Edge,
    sface_moins1: &mut Arc<dyn Surface>,
) -> Result<(), String> {
    // `Hc2 = BRep_Tool::CurveOnSurface(E2, Face[nbface - 1], Ubid, Ubid);`
    // (`cxx:2955`).
    let face_nbm1 = face_at(st, (st.nbface - 1) as usize, 2955)?;
    lp.hc2 = curve_on_surface_range(e2, &face_nbm1).map(|(c, _, _)| c);
    if lp.hc2.is_none() {
        // `ChFi3d_ExtendSurface(Sfacemoins1, prolface[0]); if (prolface[0]) {
        //  extend = true; BRE.MakeFace(faceprol[nb - 1], Sfacemoins1,
        //  F.Location(), tol); Extrema_ExtPS ext(CV2.Point(), Asurf, tol, tol,
        //  Extrema_ExtFlag_MIN); ... }` (`cxx:2960-2977`).
        *sface_moins1 = ch_fi3d_extend_surface(sface_moins1, &mut lp.prol_face[0]);
        if lp.prol_face[0] != 0 {
            lp.extend = true;
            let face_f = face_at(st, (st.nbface - 1) as usize, 2965)?;
            let tol = BRepTool::face_tolerance(&face_f);
            lp.face_prol[(st.nbface - 1) as usize] =
                Some(make_prol_face(Arc::clone(sface_moins1), tol));
            let p = *cv2.point();
            if let Some(pair) = point_surface_extrema(sface_moins1.as_ref(), &p) {
                lp.pfac2 = GpPnt2d::new(pair.u1, pair.v1.unwrap_or(0.0));
            }
        }
    } else {
        // `pfac2 = Hc2->Value(CV2.ParameterOnArc());` (`cxx:2981`).
        let hc2 = lp.hc2.as_ref().expect("checked above");
        let par = cv2
            .parameter_on_arc()
            .ok_or_else(|| "IntersectionAtEnd: CV2 is not on an arc".to_string())?;
        lp.pfac2 = hc2.value(par);
    }
    // `paredge2 = CV2.ParameterOnArc();` (`cxx:2983`).
    lp.paredge2 = cv2
        .parameter_on_arc()
        .ok_or_else(|| "IntersectionAtEnd: CV2 is not on an arc".to_string())?;
    // `pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(...));` (`cxx:2984-2991`).
    lp.pfil2 = fi_param_on_surf(fi2, st.is_first)?;
    Ok(())
}

/// `cxx:2998-3052`: the `isOnSame2` close.
#[allow(clippy::too_many_arguments)]
pub(super) fn same_end_update(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &mut FaceLoop,
    cv2: &mut CommonPoint,
    fi2: &mut FaceInterference,
    e2: &Edge,
    nb: i32,
) -> Result<(), String> {
    // `Bs.Initialize(prolface[nb - 1] ? faceprol[nb - 1] : Face[nb - 1]);`
    // (`cxx:3001-3008`).
    let bs = lp.face_adaptor(st, (nb - 1) as usize, true)?;
    let hgs = st
        .hgs
        .surface()
        .cloned()
        .ok_or_else(|| "PerformIntersectionAtEnd: HGs has no surface".to_string())?;
    let hbs = bs
        .surface()
        .cloned()
        .ok_or_else(|| "PerformIntersectionAtEnd: Bs has no surface".to_string())?;
    // `const handle<Geom_Curve>& c3df = DStr.Curve(Fi2.LineIndex()).Curve();`
    // (`cxx:3009`).
    let c3d_fi = {
        let ds = b.my_ds.change_ds();
        let line = fi2.line_index();
        if line <= 0 {
            return Err(
                "PerformIntersectionAtEnd: Fi2.LineIndex() is not in the DS".to_string(),
            );
        }
        ds.curve(line).curve().cloned().ok_or_else(|| {
            "PerformIntersectionAtEnd: DStr.Curve(Fi2.LineIndex()) is null".to_string()
        })?
    };
    // `ChFiDS_FaceInterference& Fi = Fd->ChangeInterferenceOnS2(); ...
    //  IntersUpdateOnSame(HGs, HBs, c3df, F2, F, Edge[nb], Vtx, isfirst,
    //  10 * tolapp3d, Fi, CV2, pfac2, Ufi)` (`cxx:3010-3027`).
    let f2 = st
        .f2
        .clone()
        .ok_or_else(|| "PerformIntersectionAtEnd: F2 is null".to_string())?;
    let f = face_at(st, (nb - 1) as usize, 3016)?;
    let e_nb = edge_at(st, nb as usize, 3017)?;
    let mut ufi = st.fi1.parameter(st.is_first);
    let ok = inter_update_on_same(
        &hgs,
        &hbs,
        &c3d_fi,
        &f2,
        &f,
        &e_nb,
        &st.vtx,
        st.is_first,
        10.0 * b.tolapp3d,
        fi2,
        cv2,
        &mut lp.pfac2,
        &mut ufi,
    )?;
    if !ok {
        // `throw Standard_Failure("IntersectionAtEnd: pb intersection Face -
        //  Fi");` (`cxx:3026`).
        return Err("IntersectionAtEnd: pb intersection Face - Fi".to_string());
    }
    // `if (intersOnSameFailed) { proj.Init(pfac2, extend ? C2dint2 :
    //  CurveOnSurface(E2, Face[nbface - 1], Ubid, Ubid)); paredge2 =
    //  proj.LowerDistanceParameter(); }` (`cxx:3029-3042`).
    if lp.inters_on_same_failed {
        let curve = if lp.extend {
            lp.c2d_int2.clone()
        } else {
            let face_nbm1 = face_at(st, (st.nbface - 1) as usize, 3039)?;
            curve_on_surface_range(e2, &face_nbm1).map(|(c, _, _)| c)
        };
        let curve = curve.ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        lp.paredge2 = project_2d(&curve, &lp.pfac2)?;
    }
    // `TopOpeBRepDS_Point tpoint(CV2.Point(), tolapp3d); indpoint2 =
    //  DStr.AddPoint(tpoint); stripe->SetIndexPoint(indpoint2, isfirst, 2);`
    // (`cxx:3044-3046`).
    let p = *cv2.point();
    lp.ind_point2 = add_point(b, p, b.tolapp3d);
    if let Ok(mut stripe) = st.stripe.write() {
        stripe.set_index_point(lp.ind_point2, st.is_first, 2);
    }
    // `TopExp::Vertices(Edge[nbface], vert1, vert2); arcOri = Vtx.IsSame(vert1)
    //  ? FORWARD : REVERSED; CV2.SetArc(tolapp3d, Edge[nbface], paredge2,
    //  arcOri);` (`cxx:3048-3051`).
    let e_end = edge_at(st, st.nbface as usize, 3049)?;
    let ori = arc_orientation(&e_end, &st.vtx);
    cv2.set_arc(b.tolapp3d, &e_end.0, lp.paredge2, ori);
    Ok(())
}
