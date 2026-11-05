//! `PerformIntersectionAtEnd`'s DS storage of the per-face intersection —
//! `ChFi3d_Builder_C1.cxx:3146-3247` and `:3576-3610`.
//!
//! Ported blocks, in OCCT order:
//! - [`store`] — `cxx:3146-3247`: the `p1` / `p2` range of `Cc`, the two
//!   evaluated tolerances `to1` / `to2`, the two end points in the DS
//!   (`IndexPoint` of the stripe for the first and last face, a fresh
//!   `DStr.AddPoint` in between), the curve `Cc` itself with its two
//!   `CurvePointInterference`s, the orientation `orien` of the p-curves and
//!   the two edge limitations of `Edge[0]` / `Edge[nbface]`,
//! - [`store_face`] — `cxx:3576-3585`: the `InterfPC` / `InterfPS` of the
//!   iteration and the `indpoint1 = indpoint2` carry,
//! - [`finish`] — `cxx:3587-3610`: the `CV1` / `CV2` resets, then the
//!   appends of `InterfPC` / `InterfPS` / `Interfedge` onto the faces and the
//!   edges, plus `Interfedge[nbface]` onto `Edge[nbface]`.
//!
//! Translation notes:
//! - `Interfedge` / `InterfPC` / `InterfPS` are `occ::handle` arrays whose
//!   slots may stay null (`Interfedge[nb]` on the `extend` arm of
//!   `cxx:3330-3531`, `InterfPS[nb - 1]` under `shrink[nb - 1]`); the port
//!   keeps them as `Option` and `finish` appends only the slots OCCT has
//!   filled, since a `null` handle cannot be stored in a
//!   `Vec<InterferenceHandle>`.
//! - `TopOpeBRepDS_Point& tpt = DStr.ChangePoint(indpoint1);
//!   tpt.Tolerance(std::max(tpt.Tolerance(), to1));` is
//!   [`DataStructure::change_point`] with the same in-place maximum.
//! - `occ::down_cast` of the OCC119-style handles does not exist here: the
//!   port stores the built `InterferenceHandle`s themselves, exactly the
//!   objects `ChFi3d_FilCurveInDS` / `ChFi3d_FilPointInDS` answer.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds};
use crate::topopebrep::ds::{Curve as DsCurve, Point};

use super::super::ChFi3dBuilder;
use super::corner_faces::{edge_at, face_at};
use super::face_loop::FaceLoop;
use super::{end_points, store_end_points, EndIntersection};

/// `cxx:3146-3247`: the tolerances, the two points and the curve in the DS.
pub(super) fn store(
    b: &mut ChFi3dBuilder,
    st: &EndIntersection,
    lp: &mut FaceLoop,
    cv1: &crate::chfi_ds::CommonPoint,
    cv2: &crate::chfi_ds::CommonPoint,
    nb: i32,
) -> Result<(), String> {
    let cc = lp
        .cc
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Cc is null".to_string())?;
    // `p1 = Cc->FirstParameter(); p2 = Cc->LastParameter();`
    // (`cxx:3149-3150`).
    lp.p1 = cc.first_parameter();
    lp.p2 = cc.last_parameter();
    let (p1, p2) = (lp.p1, lp.p2);

    // `HGs->D0(Pardeb(1), Pardeb(2), P1); ... HBs->D0(Parfin(3), Parfin(4),
    //  P4);` (`cxx:3155-3158`).
    let hgs = st
        .hgs
        .surface()
        .cloned()
        .ok_or_else(|| "IntersectionAtEnd: HGs has no surface".to_string())?;
    let hbs = lp
        .hbs
        .as_ref()
        .and_then(|a| a.surface().cloned())
        .ok_or_else(|| "IntersectionAtEnd: HBs has no surface".to_string())?;
    let pa = hgs.d0(lp.pardeb[0], lp.pardeb[1]);
    let pb = hgs.d0(lp.parfin[0], lp.parfin[1]);
    let pc_ = hbs.d0(lp.pardeb[2], lp.pardeb[3]);
    let pd = hbs.d0(lp.parfin[2], lp.parfin[3]);
    // `Pc->D0(p1, p2d1); Pc->D0(p2, p2d2); HBs->D0(p2d1..., P7);
    //  HBs->D0(p2d2..., P8);` (`cxx:3159-3162`).
    let pc_curve = lp
        .pc
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Pc is null".to_string())?;
    let uv7 = pc_curve.d0(p1);
    let uv8 = pc_curve.d0(p2);
    let p7 = hbs.d0(uv7.x(), uv7.y());
    let p8 = hbs.d0(uv8.x(), uv8.y());
    // `Ps->D0(p1, p2d1); Ps->D0(p2, p2d2); HGs->D0(p2d1..., P5);
    //  HGs->D0(p2d2..., P6);` (`cxx:3163-3166`) — this is where the shared
    // `p2d1` / `p2d2` stop holding the `Pc` values.
    let ps_curve = lp
        .ps
        .clone()
        .ok_or_else(|| "IntersectionAtEnd: Ps is null".to_string())?;
    lp.p2d1 = ps_curve.d0(p1);
    lp.p2d2 = ps_curve.d0(p2);
    let p5 = hgs.d0(lp.p2d1.x(), lp.p2d1.y());
    let p6 = hgs.d0(lp.p2d2.x(), lp.p2d2.y());

    // `to1 = std::max(P1.Distance(P5) + P3.Distance(P7), tolreached); to2 =
    //  std::max(P2.Distance(P6) + P4.Distance(P8), tolreached);`
    // (`cxx:3167-3168`).
    let to1 = (pa.distance(&p5) + pc_.distance(&p7)).max(lp.tol_reached);
    let to2 = (pb.distance(&p6) + pd.distance(&p8)).max(lp.tol_reached);

    // `bool Isvtx1 = false; bool Isvtx2 = false;` (`cxx:3174-3175`).
    let mut is_vtx1 = false;
    let mut is_vtx2 = false;

    // `if (nb == 1) { indpoint1 = stripe->IndexPoint(isfirst, 1); if
    //  (!CV1.IsVertex()) { tpt = DStr.ChangePoint(indpoint1);
    //  tpt.Tolerance(max(tpt.Tolerance(), to1)); } else Isvtx1 = true; }`
    // (`cxx:3177-3189`).
    if nb == 1 {
        lp.ind_point1 = st
            .stripe
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?
            .index_point(st.is_first, 1);
        if !cv1.is_vertex() {
            let ds = b.my_ds.change_ds();
            if let Some(tpt) = ds.change_point(lp.ind_point1) {
                let t = tpt.tolerance();
                tpt.set_tolerance(t.max(to1));
            }
        } else {
            is_vtx1 = true;
        }
    }
    // `if (nb == nbface) { ... indpoint2 = stripe->IndexPoint(isfirst, 2);
    //  ... } else { point = Cc->Value(Cc->LastParameter());
    //  TopOpeBRepDS_Point tpoint(point, to2); indpoint2 =
    //  DStr.AddPoint(tpoint); }` (`cxx:3190-3208`).
    if nb == st.nbface {
        lp.ind_point2 = st
            .stripe
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?
            .index_point(st.is_first, 2);
        if !cv2.is_vertex() {
            let ds = b.my_ds.change_ds();
            if let Some(tpt) = ds.change_point(lp.ind_point2) {
                let t = tpt.tolerance();
                tpt.set_tolerance(t.max(to2));
            }
        } else {
            is_vtx2 = true;
        }
    } else {
        let point = cc.d0(p2);
        lp.ind_point2 = b
            .my_ds
            .change_ds()
            .add_point(&Point::with_tolerance(point, to2));
    }
    // `if (nb != 1) { tpt = DStr.ChangePoint(indpoint1);
    //  tpt.Tolerance(max(tpt.Tolerance(), to1)); }` (`cxx:3210-3214`).
    if nb != 1 {
        let ds = b.my_ds.change_ds();
        if let Some(tpt) = ds.change_point(lp.ind_point1) {
            let t = tpt.tolerance();
            tpt.set_tolerance(t.max(to1));
        }
    }

    // `TopOpeBRepDS_Curve tcurv3d(Cc, tolreached);
    //  indcurve[nb - 1] = DStr.AddCurve(tcurv3d);` (`cxx:3215-3216`).
    let indcurv = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&cc)),
        lp.tol_reached,
        false,
    ));
    lp.ind_curve[(nb - 1) as usize] = indcurv;

    // `Interfp1 = ChFi3d_FilPointInDS(FORWARD, indcurve[nb - 1], indpoint1,
    //  Cc->FirstParameter(), Isvtx1); Interfp2 = ... (REVERSED, ...,
    //  Cc->LastParameter(), Isvtx2);` and the two appends (`cxx:3218-3231`).
    let interfp1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, indcurv, lp.ind_point1, p1, is_vtx1);
    let interfp2 =
        ch_fi3d_fil_point_in_ds(Orientation::Reversed, indcurv, lp.ind_point2, p2, is_vtx2);
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(indcurv);
        list.push(interfp1);
        list.push(interfp2);
    }

    // `TopAbs_Orientation ori = TopAbs_FORWARD; orface = Face[nb -
    //  1].Orientation(); if (orface == orsurfdata) orien =
    //  TopAbs::Reverse(orcourbe); else orien = orcourbe;`
    // (`cxx:3236-3245`).
    lp.ori = Orientation::Forward;
    lp.orface = face_at(st, (nb - 1) as usize, 3238)?.0.orientation();
    lp.orien = if lp.orface == st.or_surf_data {
        st.or_courbe.reversed()
    } else {
        st.or_courbe
    };

    // `if (nb == 1) { int Iarc1 = DStr.AddShape(Edge[0]); Interfedge[0] =
    //  ChFi3d_FilPointInDS(CV1.TransitionOnArc(), Iarc1, indpoint1, paredge1,
    //  Isvtx1); }` (`cxx:3248-3255`).
    if nb == 1 {
        let e0 = edge_at(st, 0, 3250)?;
        let iarc1 = b.my_ds.change_ds().add_shape(&e0.0);
        let tra = cv1
            .transition_on_arc()
            .ok_or_else(|| "ChFiDS_CommonPoint::TransitionOnArc".to_string())?;
        lp.interf_edge[0] = Some(ch_fi3d_fil_point_in_ds(
            tra,
            iarc1,
            lp.ind_point1,
            lp.paredge1,
            is_vtx1,
        ));
    }
    // `if (nb == nbface) { int Iarc2 = DStr.AddShape(Edge[nb]);
    //  Interfedge[nb] = ChFi3d_FilPointInDS(CV2.TransitionOnArc(), Iarc2,
    //  indpoint2, paredge2, Isvtx2); }` (`cxx:3257-3264`).
    if nb == st.nbface {
        let enb = edge_at(st, nb as usize, 3258)?;
        let iarc2 = b.my_ds.change_ds().add_shape(&enb.0);
        let tra = cv2
            .transition_on_arc()
            .ok_or_else(|| "ChFiDS_CommonPoint::TransitionOnArc".to_string())?;
        lp.interf_edge[nb as usize] = Some(ch_fi3d_fil_point_in_ds(
            tra,
            iarc2,
            lp.ind_point2,
            lp.paredge2,
            is_vtx2,
        ));
    }
    Ok(())
}/// `cxx:3576-3585`: the two face p-curves of the iteration.
pub(super) fn store_face(
    b: &mut ChFi3dBuilder,
    st: &EndIntersection,
    lp: &mut FaceLoop,
    nb: i32,
) -> Result<(), String> {
    // `indice = DStr.AddShape(Face[nb - 1]); InterfPC[nb - 1] =
    //  ChFi3d_FilCurveInDS(indcurve[nb - 1], indice, Pc, orien);`
    // (`cxx:3576-3578`).
    let face = face_at(st, (nb - 1) as usize, 3576)?;
    let indice = b.my_ds.change_ds().add_shape(&face.0);
    let indcurv = lp.ind_curve[(nb - 1) as usize];
    lp.interf_pc[(nb - 1) as usize] = Some(ch_fi3d_fil_curve_in_ds(
        indcurv,
        indice,
        lp.pc.clone(),
        lp.orien,
    ));
    // `if (!shrink[nb - 1]) { InterfPS[nb - 1] = ChFi3d_FilCurveInDS(
    //  indcurve[nb - 1], Isurf, Ps, orcourbe); }` (`cxx:3579-3582`).
    if lp.shrink[(nb - 1) as usize] == 0 {
        lp.interf_ps[(nb - 1) as usize] = Some(ch_fi3d_fil_curve_in_ds(
            indcurv,
            st.isurf,
            lp.ps.clone(),
            st.or_courbe,
        ));
    }
    // `indpoint1 = indpoint2;` (`cxx:3585`).
    lp.ind_point1 = lp.ind_point2;
    Ok(())
}

/// `cxx:3587-3610`: the tail of the loop over the faces.
pub(super) fn finish(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
    lp: &FaceLoop,
) -> Result<(), String> {
    // `if (isOnSame1) CV1.Reset(); if (isOnSame2) CV2.Reset();`
    // (`cxx:3587-3592`).
    if st.is_on_same1 || st.is_on_same2 {
        let (mut cv1, mut cv2) = end_points(st)?;
        if st.is_on_same1 {
            cv1.reset();
        }
        if st.is_on_same2 {
            cv2.reset();
        }
        store_end_points(st, &cv1, &cv2);
    }
    let nbface = st.nbface;
    for nb in 1..=nbface {
        // `int indice = DStr.AddShape(Face[nb - 1]);
        //  DStr.ChangeShapeInterferences(indice).Append(InterfPC[nb - 1]);`
        // (`cxx:3599-3600`).
        let face = face_at(st, (nb - 1) as usize, 3599)?;
        b.my_ds.change_ds().add_shape(&face.0);
        if let Some(interf) = lp.interf_pc[(nb - 1) as usize].clone() {
            b.my_ds
                .change_ds()
                .change_shape_interferences(&face.0)
                .push(interf);
        }
        // `if (!shrink[nb - 1]) DStr.ChangeSurfaceInterferences(Isurf).Append(
        //  InterfPS[nb - 1]);` (`cxx:3601-3604`).
        if lp.shrink[(nb - 1) as usize] == 0 {
            if let Some(interf) = lp.interf_ps[(nb - 1) as usize].clone() {
                b.my_ds
                    .change_ds()
                    .change_surface_interferences(st.isurf)
                    .push(interf);
            }
        }
        // `if (!proledge[nb - 1]) DStr.ChangeShapeInterferences(
        //  Edge[nb - 1]).Append(Interfedge[nb - 1]);` (`cxx:3605-3608`).
        if lp.prol_edge[(nb - 1) as usize] == 0 {
            let edge = edge_at(st, (nb - 1) as usize, 3606)?;
            if let Some(interf) = lp.interf_edge[(nb - 1) as usize].clone() {
                b.my_ds
                    .change_ds()
                    .change_shape_interferences(&edge.0)
                    .push(interf);
            }
        }
    }
    // `DStr.ChangeShapeInterferences(Edge[nbface]).Append(
    //  Interfedge[nbface]);` (`cxx:3610`).
    let e_end = edge_at(st, nbface as usize, 3610)?;
    if let Some(interf) = lp.interf_edge[nbface as usize].clone() {
        b.my_ds
            .change_ds()
            .change_shape_interferences(&e_end.0)
            .push(interf);
    }
    Ok(())
}

/// `gp_Pnt2d::SetCoord(isUShrink ? 1 : 2, prevSDParam)` — the 1-based
/// `SetCoord` of `gp_Pnt2d` (`gp_Pnt2d.hxx`), used by the shrink block.
pub(super) fn set_coord(p: &mut GpPnt2d, index: i32, value: f64) {
    if index == 1 {
        p.set_x(value);
    } else {
        p.set_y(value);
    }
}

/// `gp_Pnt2d::Coord(isUShrink ? 1 : 2)`.
pub(super) fn coord(p: &GpPnt2d, index: i32) -> f64 {
    if index == 1 {
        p.x()
    } else {
        p.y()
    }
}