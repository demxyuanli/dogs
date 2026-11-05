//! `IntersectMoreCorner`'s DS write stage — `ChFi3d_Builder_C1.cxx:5129-5241`.
//!
//! Ported blocks, in OCCT order:
//! - `cxx:5129-5156` — `IShape = DStr.AddShape(Fv)` and the orientation `Et` of
//!   the arc of `CV1` (or `CV2` when `IFadArc != 1`) inside `Fv`;
//! - `cxx:5163-5164` — `SetIndexPoint(ChFi3d_IndexPointInDS(CV1 / CV2, DStr))`;
//! - `cxx:5166-5180` — the `!intcouture` write: `Tc(Cc, tolreached)` in the DS,
//!   one interference on `Fv`, then the stripe bookkeeping. Unlike
//!   `PerformOneCorner` this arm has **no** OCC119 guard and **no** pseudo-seam
//!   block;
//! - `cxx:5181-5241` — the `intcouture` write: both pieces of `Cc` on `Fv` and
//!   on `Isurf`, the limitation of the sewing edge and the four
//!   `CurvePointInterference`s of `Icurv1` / `Icurv2`.
//!
//! Divergences that are not heuristics: the sewing-edge limitation chooses
//! `TopAbs_FORWARD` when `par1` is nearer `pard` here (`cxx:5217-5224`), the
//! opposite of `PerformOneCorner`'s `write_split` (`cxx:1555-1562`); the source
//! is followed as written.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ch_fi3d_index_point_in_ds,
};
use crate::chfi3d::builder_c1::compute_curve2d;
use crate::shape::Face;
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::ChFi3dBuilder;
use super::MoreCorner;

/// `cxx:5129-5156`: `IShape` and `Et`.
fn face_index_and_orientation(
    mc: &MoreCorner,
    b: &mut ChFi3dBuilder,
    fv: &Face,
) -> Result<(i32, Orientation), String> {
    // `int IShape = DStr.AddShape(Fv);` (`cxx:5129`).
    let ishape = b.my_ds.change_ds().add_shape(&fv.0);
    // `TopAbs_Orientation Et = TopAbs_FORWARD;` (`cxx:5130`).
    let mut et = Orientation::Forward;
    let cv = mc.fd.vertex(mc.is_first, mc.i_fad_arc);
    let Some(target) = cv.arc() else {
        // The loop body never matches when the common point is off arc.
        return Ok((ishape, et));
    };
    let trans = cv
        .transition_on_arc()
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let forward = fv.0.oriented(Orientation::Forward);
    for e in explorer_shapes(&forward, ShapeType::Edge) {
        if is_same(&e, target) {
            let composed = Orientation::compose(e.orientation(), trans);
            // `cxx:5138` reverses, `cxx:5150` does not.
            et = if mc.i_fad_arc == 1 {
                composed.reversed()
            } else {
                composed
            };
            break;
        }
    }
    Ok((ishape, et))
}

/// `cxx:5129-5241`.
pub(super) fn write(b: &mut ChFi3dBuilder, mc: &mut MoreCorner) -> Result<(), String> {
    let is_first = mc.is_first;
    let fv = mc
        .fv
        .clone()
        .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
    let cv1 = mc.fd.vertex(is_first, 1).clone();
    let cv2 = mc.fd.vertex(is_first, 2).clone();

    // `cxx:5129-5156`.
    let (ishape, et) = face_index_and_orientation(mc, b, &fv)?;
    // `stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);`
    // and `CV2` on 2 (`cxx:5163-5164`).
    let ind_cv1 = ch_fi3d_index_point_in_ds(&cv1, b.my_ds.change_ds());
    let ind_cv2 = ch_fi3d_index_point_in_ds(&cv2, b.my_ds.change_ds());
    {
        let mut stripe = mc
            .stripe
            .write()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        stripe.set_index_point(ind_cv1, is_first, 1);
        stripe.set_index_point(ind_cv2, is_first, 2);
    }

    if !mc.intcouture {
        // `TopOpeBRepDS_Curve Tc(Cc, tolreached); ICurve = DStr.AddCurve(Tc);`
        // (`cxx:5172-5173`).
        let cc = mc
            .cc
            .clone()
            .ok_or_else(|| "OneCorner : failed calculation intersection".to_string())?;
        let icurve = b
            .my_ds
            .change_ds()
            .add_curve(&mut DsCurve::with_curve(Some(cc), mc.tolreached, false));
        // `Interfc = ChFi3d_FilCurveInDS(ICurve, IShape, Pc, Et); DStr
        //  .ChangeShapeInterferences(IShape).Append(Interfc);` (`cxx:5174-5176`).
        let pc = mc
            .pc
            .clone()
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        let interfc = ch_fi3d_fil_curve_in_ds(icurve, ishape, Some(pc), et);
        b.my_ds
            .change_ds()
            .change_shape_interferences(&fv.0)
            .push(interfc);

        // `stripe->ChangePCurve(isfirst) = Ps; stripe->SetCurve(ICurve,
        //  isfirst); stripe->SetParameters(isfirst, Udeb, Ufin);`
        // (`cxx:5177-5179`).
        let mut stripe = mc
            .stripe
            .write()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        *stripe.change_p_curve(is_first) = mc.ps.clone();
        stripe.set_curve(icurve, is_first);
        stripe.set_parameters(is_first, mc.udeb, mc.ufin);
        return Ok(());
    }

    // `cxx:5181-5241`.
    let curv1 = mc
        .curv1
        .clone()
        .ok_or_else(|| "OneCorner : failed calculation intersection".to_string())?;
    let curv2 = mc
        .curv2
        .clone()
        .ok_or_else(|| "OneCorner : failed calculation intersection".to_string())?;

    // `ComputeCurve2d(curv1, Fv, c2d1); InterFv = ChFi3d_FilCurveInDS(Icurv1,
    //  IShape, c2d1, Et); DStr.ChangeShapeInterferences(IShape).Append(...)`
    // (`cxx:5188-5194`).
    let c2d1 = compute_curve2d(&curv1, &fv);
    let inter_fv = ch_fi3d_fil_curve_in_ds(mc.icurv1, ishape, c2d1, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .push(inter_fv);
    let c2d2 = compute_curve2d(&curv2, &fv);
    let inter_fv = ch_fi3d_fil_curve_in_ds(mc.icurv2, ishape, c2d2, et);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&fv.0)
        .push(inter_fv);

    // `if (Fd->Orientation() == Fv.Orientation()) Et = TopAbs::Reverse(Et);`
    // (`cxx:5196-5199`).
    let mut et = et;
    if mc.fd.orientation() == fv.0.orientation() {
        et = et.reversed();
    }
    // `c2d1 = new Geom2d_TrimmedCurve(Ps, Udeb, par2); InterFv =
    //  ChFi3d_FilCurveInDS(Icurv1, Isurf, c2d1, Et); DStr
    //  .ChangeSurfaceInterferences(Isurf).Append(InterFv);` (`cxx:5200-5202`).
    let isurf = mc.fd.surf();
    let ps = mc
        .ps
        .clone()
        .ok_or_else(|| "OneCorner : failed calculation intersection".to_string())?;
    let c2d1 = Arc::new(Geom2dTrimmedCurve::new(Arc::clone(&ps), mc.udeb, mc.par2))
        as Arc<dyn Curve2d>;
    let inter = ch_fi3d_fil_curve_in_ds(mc.icurv1, isurf, Some(c2d1), et);
    b.my_ds
        .change_ds()
        .change_surface_interferences(isurf)
        .push(inter);
    // `c2d2 = new Geom2d_TrimmedCurve(Ps, par2, Ufin); ...` (`cxx:5203-5205`).
    let c2d2 = Arc::new(Geom2dTrimmedCurve::new(ps, mc.par2, mc.ufin)) as Arc<dyn Curve2d>;
    let inter = ch_fi3d_fil_curve_in_ds(mc.icurv2, isurf, Some(c2d2), et);
    b.my_ds
        .change_ds()
        .change_surface_interferences(isurf)
        .push(inter);

    // `cxx:5207-5226`: limitation of the sewing edge.
    let edgecouture = mc
        .edgecouture
        .clone()
        .ok_or_else(|| "IntersectMoreCorner: no sewing edge".to_string())?;
    let iarc = b.my_ds.change_ds().add_shape(&edgecouture.0);
    let vdeb = first_vertex(&edgecouture.0)
        .ok_or_else(|| "IntersectMoreCorner: sewing edge has no first vertex".to_string())?;
    let vfin = last_vertex(&edgecouture.0)
        .ok_or_else(|| "IntersectMoreCorner: sewing edge has no last vertex".to_string())?;
    let pard = BRepTool::parameter_of_vertex(&edgecouture.0, &vdeb);
    let parf = BRepTool::parameter_of_vertex(&edgecouture.0, &vfin);
    // `if (std::abs(par1 - pard) < std::abs(parf - par1)) ori = TopAbs_FORWARD;
    //  else ori = TopAbs_REVERSED;` (`cxx:5217-5224`) — the opposite choice of
    // `PerformOneCorner`'s `write_split`.
    let ori = if (mc.par1 - pard).abs() < (parf - mc.par1).abs() {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let inter = ch_fi3d_fil_point_in_ds(ori, iarc, mc.indpt, mc.par1, false);
    b.my_ds
        .change_ds()
        .change_shape_interferences(&edgecouture.0)
        .push(inter);

    // `cxx:5228-5240`: `stripe->InDS(isfirst);` then the two indexes and the
    // four `CurvePointInterference`s.
    let (ind1, ind2) = {
        let mut stripe = mc
            .stripe
            .write()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        stripe.in_ds(is_first, 1);
        (stripe.index_point(is_first, 1), stripe.index_point(is_first, 2))
    };
    let i1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, mc.icurv1, ind1, mc.udeb, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(mc.icurv1)
        .push(i1);
    let i2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, mc.icurv1, mc.indpt, mc.par2, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(mc.icurv1)
        .push(i2);
    let i3 = ch_fi3d_fil_point_in_ds(Orientation::Forward, mc.icurv2, mc.indpt, mc.par2, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(mc.icurv2)
        .push(i3);
    let i4 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, mc.icurv2, ind2, mc.ufin, false);
    b.my_ds
        .change_ds()
        .change_curve_interferences(mc.icurv2)
        .push(i4);
    Ok(())
}
