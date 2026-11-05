//! `PerformTwoCornerSameExt` - the intersection of two stripes that share
//! their common points - and the `nconges == 2` stage that reaches it.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:810-951`
//! for the static function and `:2642-2705` for the block of
//! `PerformMoreThreeCorner` that calls it (OCCT 8.0.0).
//!
//! The stage is only reached with two fillets on the same vertex
//! (`nconges == 2`) that are not tangent (`!deuxconges`). It intersects the two
//! fillet surfaces and keeps the result only when its two ends fall on the
//! common points the two stripes already share; otherwise the caller falls
//! through to the `GeomPlate` section (`cxx:2707` on, ported in
//! `plate_stripe_constraints` through `plate_tails`).
//!
//! Translation notes:
//! - OCCT hands `Com11` / `Com12` / `Com21` out as `ChFiDS_CommonPoint&`
//!   references into the stripe's `SurfData` (`cxx:835-838`); they are only
//!   read here, so the port works on copies,
//! - `PCurveOnSurf()->Value(U)` is [`Curve2d::value`] and
//!   `PCurveOnFace()->D0(U, p2d)` is [`Curve2d::value`] as well,
//! - `ChFi3d_ComputeCurves` keeps the *untrimmed* intersection curve plus its
//!   `(Pardeb, Parfin)` range (`builder_0_helpers::compute_curves`), so the
//!   `cint->FirstParameter()` / `LastParameter()` of `cxx:882-883` and
//!   `cxx:938-939` are [`ComputeCurves::pardeb`] / [`ComputeCurves::parfin`],
//! - the call at `cxx:880` passes ten arguments, so `wholeCurv` keeps its
//!   `true` default (`ChFi3d_Builder_0.hxx:396-406`): the intersection curve is
//!   not trimmed to the two known extremities,
//! - `TopOpeBRepDS_Curve tcurv3d(cint, tol)` is [`DsCurve::with_curve`] with
//!   `IsWalk = false` (the two-argument OCCT constructor,
//!   `TopOpeBRepDS_Curve.cxx:42-52`),
//! - `DStr.Shape(I).Orientation()` uses the default `FindKeep = true`
//!   (`TopOpeBRepDS_DataStructure.hxx:254`), and a non-face shape there is
//!   OCCT's `TopoDS::Face` throw, answered as `Err`,
//! - `trouve` in the stage (`cxx:2651`, `:2657`) is a dead store: the next read
//!   of that local is at `cxx:2777`, after it is reset, so the port drops it.

use std::sync::{Arc, RwLockReadGuard};

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_surf, ch_fi3d_compute_curves, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds,
    ch_fi3d_index_point_in_ds, ComputeCurves,
};
use crate::chfi_ds::{CommonPoint, FaceInterference, Stripe, StripeHandle, SurfData};
use crate::shape::Face;
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::calcul::calcul_p2d_on_surf;
use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;

/// `1.e-4`: the coincidence tolerance of `cxx:858`, `:868`, `:884-886`,
/// `:925` and `:933`.
const SAME_POINT: f64 = 1.0e-4;

/// `1.e-7`: the tolerance of the stage's `intersection` test (`cxx:2686-2688`).
const INTERSECTION_TOL: f64 = 1.0e-7;

/// `ChFi3d_ComputeCurves(..., 1.e-4, 1.e-5, tol)` (`cxx:880`).
const TOL3D: f64 = 1.0e-4;
/// The `tol2d` argument of that call (`cxx:880`).
const TOL2D: f64 = 1.0e-5;

/// The read half of [`StripeHandle`], with OCCT's null-handle failure as the
/// error arm.
fn read_stripe(handle: &StripeHandle) -> Result<RwLockReadGuard<'_, Stripe>, String> {
    handle
        .read()
        .map_err(|_| "PerformTwoCornerSameExt: poisoned stripe".to_string())
}

/// `stripe->SetOfSurfData()->Value(index)` (`cxx:835-846`), `index` 1-based.
fn surf_data(stripe: &StripeHandle, index: i32) -> Result<Arc<SurfData>, String> {
    if index < 1 {
        return Err(format!(
            "PerformTwoCornerSameExt: SurfData {index} is out of range"
        ));
    }
    let guard = read_stripe(stripe)?;
    guard
        .set_of_surf_data()
        .get(index as usize - 1)
        .cloned()
        .ok_or_else(|| format!("PerformTwoCornerSameExt: SurfData {index} is out of range"))
}

/// `Fi.PCurveOnSurf()->Value(Fi.Parameter(isfirst))` (`cxx:853-866`).
fn surf_value(fi: &FaceInterference, is_first: bool) -> Result<GpPnt2d, String> {
    let pc = fi
        .p_curve_on_surf()
        .ok_or_else(|| "PerformTwoCornerSameExt: the interference has no p-curve".to_string())?;
    Ok(pc.value(fi.parameter(is_first)))
}

/// Everything `PerformTwoCornerSameExt` reads out of the two stripes before it
/// starts writing to the DS.
struct CornerData {
    com11: CommonPoint,
    com12: CommonPoint,
    com21: CommonPoint,
    indic1: i32,
    indic2: i32,
    /// `Pardeb(1..4)` (`cxx:869-874`).
    pardeb: [f64; 4],
    /// `Parfin(1..4)` (`cxx:872-877`).
    parfin: [f64; 4],
    /// `Fi11.Transition()` (`cxx:912`).
    transition: Orientation,
    /// `Fi11.FirstParameter()` / `Fi11.LastParameter()` (`cxx:918-923`).
    fi11_first: f64,
    fi11_last: f64,
    /// `Fi11.PCurveOnFace()` (`cxx:918`, `:922`).
    pc_face: Option<Arc<dyn Curve2d>>,
    /// `Fd1->Orientation()` (`cxx:909`).
    orsurf1: Orientation,
    /// `Fd2->Orientation()` (`cxx:914`).
    orsurf2: Orientation,
    /// `Fd1->IndexOfS1()` (`cxx:910`, `:920`).
    index_of_s1: i32,
}

/// `cxx:834-878`: the common points, the surface indices, the four UV
/// extremities and the orientation / p-curve data, together with the two
/// bounded surfaces of `ChFi3d_BoundSurf(DStr, Fd, 1, 2)` (`cxx:877-878`).
fn read_corners(
    ds: &DataStructure,
    stripe1: &StripeHandle,
    index1: i32,
    sens1: i32,
    stripe2: &StripeHandle,
    index2: i32,
    sens2: i32,
) -> Result<(CornerData, SurfaceAdaptor, SurfaceAdaptor), String> {
    let fd1 = surf_data(stripe1, index1)?;
    let fd2 = surf_data(stripe2, index2)?;

    // `isfirst = sens1 == 1;` then `ChangeVertex(isfirst, 1)` / `(isfirst, 2)`
    // (`cxx:834-837`).
    let is_first1 = sens1 == 1;
    let com11 = fd1.vertex(is_first1, 1).clone();
    let com12 = fd1.vertex(is_first1, 2).clone();
    // `isfirst = sens2 == 1; ChangeVertex(isfirst, 1)` (`cxx:838-839`).
    let is_first2 = sens2 == 1;
    let com21 = fd2.vertex(is_first2, 1).clone();

    // `indic1 = ...Surf(); indic2 = ...Surf();` (`cxx:843-844`) and the four
    // `InterferenceOnS1()` / `InterferenceOnS2()` references (`cxx:848-851`).
    let fi11 = fd1.interference(1);
    let fi12 = fd1.interference(2);
    let fi21 = fd2.interference(1);
    let fi22 = fd2.interference(2);

    // `pfi11 = Fi11.PCurveOnSurf()->Value(Fi11.Parameter(isfirst));` and
    // `pfi12` (`cxx:853-855`).
    let pfi11 = surf_value(fi11, is_first1)?;
    let pfi12 = surf_value(fi12, is_first1)?;
    // `cxx:858-866`: the ends are swapped when the two stripes meet crosswise.
    let (pfi21, pfi22) = if com11.point().distance(com21.point()) < SAME_POINT {
        (surf_value(fi21, is_first2)?, surf_value(fi22, is_first2)?)
    } else {
        (surf_value(fi22, is_first2)?, surf_value(fi21, is_first2)?)
    };

    let hs1 = ch_fi3d_bound_surf(ds, &fd1, 1, 2)?;
    let hs2 = ch_fi3d_bound_surf(ds, &fd2, 1, 2)?;

    let corners = CornerData {
        com11,
        com12,
        com21,
        indic1: fd1.surf(),
        indic2: fd2.surf(),
        // `Pardeb(1) = pfi11.X(); ... Parfin(4) = pfi22.Y();` (`cxx:868-877`).
        pardeb: [pfi11.x(), pfi11.y(), pfi21.x(), pfi21.y()],
        parfin: [pfi12.x(), pfi12.y(), pfi22.x(), pfi22.y()],
        transition: fi11.transition(),
        fi11_first: fi11.first_parameter(),
        fi11_last: fi11.last_parameter(),
        pc_face: fi11.p_curve_on_face().cloned(),
        orsurf1: fd1.orientation(),
        orsurf2: fd2.orientation(),
        index_of_s1: fd1.index_of_s1(),
    };
    Ok((corners, hs1, hs2))
}

/// `PerformTwoCornerSameExt(DStr, stripe1, index1, sens1, stripe2, index2,
/// sens2, trouve)` (`cxx:810-951`).
///
/// Answers OCCT's `trouve`: `false` leaves the DS untouched, as in the source.
#[allow(clippy::too_many_arguments)]
pub(super) fn perform_two_corner_same_ext(
    ds: &mut DataStructure,
    stripe1: &StripeHandle,
    index1: i32,
    sens1: i32,
    stripe2: &StripeHandle,
    index2: i32,
    sens2: i32,
) -> Result<bool, String> {
    let (c, hs1, hs2) = read_corners(ds, stripe1, index1, sens1, stripe2, index2, sens2)?;

    // `if (ChFi3d_ComputeCurves(HS1, HS2, Pardeb, Parfin, cint, C2dint1,
    // C2dint2, 1.e-4, 1.e-5, tol))` (`cxx:880`). The ten-argument call keeps
    // `wholeCurv = true` (`ChFi3d_Builder_0.hxx:406`).
    let mut computed = ComputeCurves::default();
    if !ch_fi3d_compute_curves(
        &hs1,
        &hs2,
        c.pardeb,
        c.parfin,
        &mut computed,
        TOL3D,
        TOL2D,
        true,
    ) {
        return Ok(false);
    }
    let cint: Arc<dyn Curve> = computed.c3d.clone().ok_or_else(|| {
        "PerformTwoCornerSameExt: ChFi3d_ComputeCurves answers true without a curve".to_string()
    })?;

    // `cint->D0(cint->FirstParameter(), P1); cint->D0(cint->LastParameter(),
    // P2);` (`cxx:882-883`).
    let p1 = cint.d0(computed.pardeb);
    let p2 = cint.d0(computed.parfin);

    // `trouve = ((Com11.Point().Distance(P1) < 1.e-4 || Com11.Point()
    // .Distance(P2) < 1.e-4) && (Com12 ...));` (`cxx:884-885`).
    let near = |com: &CommonPoint, p: &GpPnt| com.point().distance(p) < SAME_POINT;
    let trouve = (near(&c.com11, &p1) || near(&c.com11, &p2))
        && (near(&c.com12, &p1) || near(&c.com12, &p2));
    if !trouve {
        return Ok(false);
    }

    // `isfirst = sens1 == 1; stripe1->InDS(isfirst);` (`cxx:890-891`).
    let is_first1 = sens1 == 1;
    {
        let mut guard = stripe1
            .write()
            .map_err(|_| "PerformTwoCornerSameExt: poisoned stripe".to_string())?;
        guard.in_ds(is_first1, 1);
    }
    // `indpoint1 = ChFi3d_IndexPointInDS(Com11, DStr); indpoint2 =
    // ChFi3d_IndexPointInDS(Com12, DStr);` (`cxx:892-893`).
    let indpoint1 = ch_fi3d_index_point_in_ds(&c.com11, ds);
    let indpoint2 = ch_fi3d_index_point_in_ds(&c.com12, ds);
    // `stripe1->SetIndexPoint(indpoint1, isfirst, 1); ... (indpoint2, 2);`
    // (`cxx:894-895`).
    {
        let mut guard = stripe1
            .write()
            .map_err(|_| "PerformTwoCornerSameExt: poisoned stripe".to_string())?;
        guard.set_index_point(indpoint1, is_first1, 1);
        guard.set_index_point(indpoint2, is_first1, 2);
    }

    // `isfirst = sens2 == 1; stripe2->InDS(isfirst);` (`cxx:896-897`).
    let is_first2 = sens2 == 1;
    {
        let mut guard = stripe2
            .write()
            .map_err(|_| "PerformTwoCornerSameExt: poisoned stripe".to_string())?;
        guard.in_ds(is_first2, 1);
    }
    // `cxx:898-907`: the two point indices follow the `Com11` / `Com21`
    // coincidence test of `cxx:858`.
    {
        let mut guard = stripe2
            .write()
            .map_err(|_| "PerformTwoCornerSameExt: poisoned stripe".to_string())?;
        if c.com11.point().distance(c.com21.point()) < SAME_POINT {
            guard.set_index_point(indpoint1, is_first2, 1);
            guard.set_index_point(indpoint2, is_first2, 2);
        } else {
            guard.set_index_point(indpoint2, is_first2, 1);
            guard.set_index_point(indpoint1, is_first2, 2);
        }
    }

    // `orsurf1 = Fd1->Orientation(); trafil1 = DStr.Shape(Fd1->IndexOfS1())
    // .Orientation(); trafil1 = TopAbs::Compose(trafil1, Fd1->Orientation());
    // trafil1 = TopAbs::Compose(TopAbs::Reverse(Fi11.Transition()), trafil1);`
    // (`cxx:909-913`).
    let orsurf1 = c.orsurf1;
    let face_shape = ds.shape(c.index_of_s1, true).ok_or_else(|| {
        format!(
            "PerformTwoCornerSameExt: DStr.Shape({}) is null",
            c.index_of_s1
        )
    })?;
    let mut trafil1 = face_shape.orientation();
    trafil1 = Orientation::compose(trafil1, c.orsurf1);
    trafil1 = Orientation::compose(c.transition.reversed(), trafil1);

    // `orsurf2 = Fd2->Orientation();` (`cxx:914`).
    let orsurf2 = c.orsurf2;

    // `TopOpeBRepDS_Curve tcurv3d(cint, tol); indcurve =
    // DStr.AddCurve(tcurv3d);` (`cxx:915-916`).
    let tolreached = computed.tolreached;
    let indcurve = ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&cint)), tolreached, false));

    // `Fi11.PCurveOnFace()->D0(Fi11.LastParameter(), p2d); Stemp =
    // BRep_Tool::Surface(TopoDS::Face(DStr.Shape(Fd1->IndexOfS1())));
    // Stemp->D0(p2d.X(), p2d.Y(), P4); Fi11.PCurveOnFace()->D0(
    // Fi11.FirstParameter(), p2d); Stemp->D0(p2d.X(), p2d.Y(), P3);`
    // (`cxx:918-923`).
    let pc_face = c.pc_face.clone().ok_or_else(|| {
        "PerformTwoCornerSameExt: Fi11 has no p-curve on the face".to_string()
    })?;
    let face = Face::wrap(face_shape).ok_or_else(|| {
        format!(
            "PerformTwoCornerSameExt: DStr.Shape({}) is not a face",
            c.index_of_s1
        )
    })?;
    let stemp = BRepTool::face_surface(&face)
        .ok_or_else(|| "PerformTwoCornerSameExt: BRep_Tool::Surface is null".to_string())?;
    let uv_last = pc_face.value(c.fi11_last);
    let p4 = stemp.d0(uv_last.x(), uv_last.y());
    let uv_first = pc_face.value(c.fi11_first);
    let p3 = stemp.d0(uv_first.x(), uv_first.y());

    // `if (P1.Distance(P4) < 1.e-4 || P2.Distance(P3) < 1.e-4) orpcurve =
    // trafil1; else orpcurve = TopAbs::Reverse(trafil1);` (`cxx:925-931`).
    let orpcurve = if p1.distance(&p4) < SAME_POINT || p2.distance(&p3) < SAME_POINT {
        trafil1
    } else {
        trafil1.reversed()
    };

    // `if (Com11.Point().Distance(P1) > 1.e-4) { ind = indpoint1; indpoint1 =
    // indpoint2; indpoint2 = ind; }` (`cxx:933-938`).
    let (indpoint1, indpoint2) = if c.com11.point().distance(&p1) > SAME_POINT {
        (indpoint2, indpoint1)
    } else {
        (indpoint1, indpoint2)
    };

    // `Interfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indcurve, indpoint1,
    // cint->FirstParameter()); Interfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED,
    // indcurve, indpoint2, cint->LastParameter());` (`cxx:937-939`).
    let interfp1 = ch_fi3d_fil_point_in_ds(Orientation::Forward, indcurve, indpoint1, computed.pardeb, false);
    let interfp2 = ch_fi3d_fil_point_in_ds(Orientation::Reversed, indcurve, indpoint2, computed.parfin, false);
    ds.change_curve_interferences(indcurve).push(interfp1);
    ds.change_curve_interferences(indcurve).push(interfp2);

    // `Interfc = ChFi3d_FilCurveInDS(indcurve, indic1, C2dint1, orpcurve);
    // DStr.ChangeSurfaceInterferences(indic1).Append(Interfc);`
    // (`cxx:942-944`).
    let interfc = ch_fi3d_fil_curve_in_ds(indcurve, c.indic1, computed.pc1.clone(), orpcurve);
    ds.change_surface_interferences(c.indic1).push(interfc);

    // `if (orsurf1 == orsurf2) orpcurve = TopAbs::Reverse(orpcurve);
    // Interfc = ChFi3d_FilCurveInDS(indcurve, indic2, C2dint2, orpcurve);`
    // (`cxx:945-948`).
    let orpcurve = if orsurf1 == orsurf2 {
        orpcurve.reversed()
    } else {
        orpcurve
    };
    let interfc = ch_fi3d_fil_curve_in_ds(indcurve, c.indic2, computed.pc2.clone(), orpcurve);
    ds.change_surface_interferences(c.indic2).push(interfc);

    Ok(true)
}

impl MoreThreeCorner {
    /// `cxx:2642-2705`: the two-fillet / non-tangent tail of
    /// `PerformMoreThreeCorner`.
    ///
    /// Answers `true` when `PerformTwoCornerSameExt` succeeded, which is the
    /// `return` of `cxx:2698-2701`. A `false` answer leaves the caller to the
    /// `GeomPlate` section (`cxx:2707` on, `plate_stripe_constraints`).
    pub(super) fn two_corner_same_ext(
        &self,
        ds: &mut DataStructure,
        nconges: usize,
    ) -> Result<bool, String> {
        // `if (nconges == 2 && !deuxconges)` (`cxx:2647`).
        if nconges != 2 || self.deuxconges {
            return Ok(false);
        }
        let nedge = self.nedge as usize;

        // `trouve = false; for (ic = 0; ic < nedge && !trouve; ic++) { if
        // (!sharp.Value(ic)) { ic1 = ic; trouve = true; } }` (`cxx:2651-2658`).
        // The `trouve` store has no reader before `cxx:2777` resets it.
        let mut ic1 = 0usize;
        for ic in 0..nedge {
            if !self.sharp[ic] {
                ic1 = ic;
                break;
            }
        }
        // `for (ic = 0; ic < nedge; ic++) { if (!sharp.Value(ic) && ic != ic1)
        // ic2 = ic; }` (`cxx:2660-2666`).
        let mut ic2 = 0usize;
        for ic in 0..nedge {
            if !self.sharp[ic] && ic != ic1 {
                ic2 = ic;
            }
        }

        // `jfp = 3 - jf.Value(ic1); Indices(nedge, ic1, icplus, icmoins);
        // Calcul_P2dOnSurf(CD.Value(ic1), jfp, i.Value(ic1, icmoins),
        // p.Value(ic1, icmoins), p2); indice = SurfIndex(CD, ic1,
        // i.Value(ic1, icmoins), ChFiSURFACE); DStr.Surface(indice)
        // .Surface()->D0(p2.X(), p2.Y(), P1);` (`cxx:2668-2671`).
        let (icplus1, icmoins1) = indices(nedge, ic1);
        let p1 = self.surf_point(ds, ic1, 3 - self.jf[ic1], icmoins1)?;
        let p2 = self.surf_point(ds, ic1, self.jf[ic1], icplus1)?;
        // `cxx:2674-2684`: the same pair on the second stripe.
        let (icplus2, icmoins2) = indices(nedge, ic2);
        let p3 = self.surf_point(ds, ic2, 3 - self.jf[ic2], icmoins2)?;
        let p4 = self.surf_point(ds, ic2, self.jf[ic2], icplus2)?;

        // `intersection = (P1.Distance(P4) <= 1.e-7 || P1.Distance(P3) <=
        // 1.e-7) && (P2.Distance(P4) <= 1.e-7 || P2.Distance(P3) <= 1.e-7);`
        // (`cxx:2686-2688`).
        let intersection = (p1.distance(&p4) <= INTERSECTION_TOL
            || p1.distance(&p3) <= INTERSECTION_TOL)
            && (p2.distance(&p4) <= INTERSECTION_TOL || p2.distance(&p3) <= INTERSECTION_TOL);
        if !intersection {
            return Ok(false);
        }

        // `PerformTwoCornerSameExt(DStr, CD.Value(ic1), Index.Value(ic1),
        // sens.Value(ic1), CD.Value(ic2), Index.Value(ic2), sens.Value(ic2),
        // introuve); if (introuve) return;` (`cxx:2690-2701`).
        let stripe1 = self.stripe_at(ic1)?.clone();
        let stripe2 = self.stripe_at(ic2)?.clone();
        let index1 = self.index[ic1];
        let sens1 = self.sens[ic1];
        let index2 = self.index[ic2];
        let sens2 = self.sens[ic2];
        perform_two_corner_same_ext(ds, &stripe1, index1, sens1, &stripe2, index2, sens2)
    }

    /// `Calcul_P2dOnSurf` then `DStr.Surface(SurfIndex(...))
    /// .Surface()->D0(p2.X(), p2.Y(), P)` (`cxx:2669-2671`, `:2674-2684`): the
    /// 3d point of the `(Jfposit, I, P)` cell of `slot`.
    fn surf_point(
        &self,
        ds: &DataStructure,
        slot: usize,
        jf_posit: i32,
        cell: usize,
    ) -> Result<GpPnt, String> {
        let indice = self.i[slot][cell];
        let param = self.p[slot][cell];
        let p2 = {
            let stripe = self.stripe_at(slot)?;
            let guard = read_stripe(stripe)?;
            calcul_p2d_on_surf(&guard, jf_posit, indice, param)?
        };
        let surface_index = surf_index(&self.cd, slot, indice as usize, SurfType::ChFiSurface)?;
        let surface = ds
            .surface(surface_index)
            .surface()
            .cloned()
            .ok_or_else(|| {
                format!("PerformTwoCornerSameExt: DStr.Surface({surface_index}) is null")
            })?;
        Ok(surface.d0(p2.x(), p2.y()))
    }
}
