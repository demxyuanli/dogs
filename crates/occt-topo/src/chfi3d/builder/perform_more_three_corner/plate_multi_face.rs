//! The multi-face projected curves and the tangent free borders of the
//! `GeomPlate` section.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx`
//! (OCCT 8.0.0):
//!
//! - `cxx:3262-3407`: for every `moresurf` slot, one `GeomPlate` constraint per
//!   projected edge `Eproj(nb)` across `Fproj(nb)`, with the point interference
//!   of the shared `Ecom(nb)` arc and the curve interference on the forward
//!   face of `Fproj(nb)`,
//! - `cxx:3408-3488`: the `droit` case, where the two free borders are tangent:
//!   one trimmed copy of the living edge, reversed when the free end is `V1`,
//!   as a G0 constraint carrying the `V1` shape on one end.
//!
//! Translation notes:
//! - `Eproj(1..nbface)` / `Ecom(1..nbface)` / `parcom(1..nbface)` are the
//!   1-based sequences of the state, read `get(nb - 1)`; `Fproj(1..nbface)` is
//!   the 0-based [`MoreThreeCorner::fproj`], read `[nb - 1]`,
//! - `BRep_Tool::Curve(E, u1, u2)` is [`BRepTool::edge_curve`] plus
//!   [`BRepTool::edge_parameters`], the same reading as `curve_hermite.rs`,
//! - the `TopExp_Explorer` walk of `cxx:3351-3358` is [`cumulated_children`]
//!   over the forward face, the same reading as
//!   [`super::plate_helpers::orientation_arete_vive_consecutive`],
//! - `tpt1.Tolerance(tpt1.Tolerance() + error)` then
//!   `tpt2.Tolerance(tpt1.Tolerance() + error)` (`cxx:3279-3280`) reads the
//!   *updated* `tpt1` tolerance for `tpt2`; the port keeps that order,
//! - `proj2d` is the `Geom2d_TrimmedCurve` over the range
//!   `BRep_Tool::CurveOnSurface` answered, and `cproj` the `Geom_TrimmedCurve`
//!   over the edge's own range (`cxx:3306-3311` overwrite `up1` / `up2`),
//! - `ctrim->Reverse()` / `ctrim2d->Reverse()` (`cxx:3447-3461`) are
//!   [`Curve::reverse`] / [`Curve2d::reverse`] before the handles are shared,
//! - `GeomAdaptor_Surface` is the stored `BRep_Tool::Surface` handle itself, the
//!   convention of `plate_stripe_constraints.rs`.

use std::sync::Arc;

use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_plate::{BoundaryCurve, BuildPlateSurface, CurveConstraint};
use occt_geom::{Curve, GeomTrimmedCurve};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds};
use crate::iterator::cumulated_children;
use crate::shape::{Face, Vertex};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure, Point as DsPoint};

use super::surf_index::indices;
use super::MoreThreeCorner;

/// `GeomPlate_CurveConstraint(..., 10, ...)` (`cxx:3340`, `:3476`).
const NBCURVPNT: i32 = 10;
/// `1.e-4`, the tolerance `TopOpeBRepDS_Curve tcurv3d(ctrim, 1.e-4)` carries
/// (`cxx:3479`).
const CURVE_TOL: f64 = 1.0e-4;
/// `1.e-3`, the `P1` / `tpt1.Point()` coincidence of `cxx:3317` and `:3321`.
const POINT_SAME: f64 = 1.0e-3;

impl MoreThreeCorner {
    /// `cxx:3262-3407`: the DS storage of the curves projected on several
    /// faces.
    pub(super) fn plate_multi_face_curves(
        &mut self,
        ds: &mut DataStructure,
        psurf: &mut BuildPlateSurface,
        n3d: &mut usize,
        tolapp3d: f64,
        angular: f64,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;
        let nbface = self.nbface as usize;

        // `V1` as a vertex: `BRep_Tool::Pnt(V1)` (`cxx:3290`).
        let v1_point = {
            let v1 = Vertex::wrap(self.vtx.clone())
                .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;
            BRepTool::vertex_point(&v1)
        };

        // `for (ic = 0; ic < nedge; ic++)` (`cxx:3263-3406`).
        for ic in 0..nedge {
            if !self.moresurf[ic] {
                continue;
            }
            let (icplus, _icmoins) = indices(nedge, ic);

            // `TopOpeBRepDS_Point& tpt1 = DStr.ChangePoint(indpoint(ic, 1));
            //  TopOpeBRepDS_Point& tpt2 = DStr.ChangePoint(indpoint(icplus, 0));
            //  tpt1.Tolerance(tpt1.Tolerance() + error);
            //  tpt2.Tolerance(tpt1.Tolerance() + error);` (`cxx:3277-3280`).
            let error = self.hermite_error;
            let t1 = self.indpoint[ic][1];
            let t2 = self.indpoint[icplus][0];
            let tpt1_tol = ds
                .change_point(t1)
                .map(|tpt| {
                    let tol = tpt.tolerance() + error;
                    tpt.set_tolerance(tol);
                    tol
                })
                .ok_or_else(|| format!("PerformMoreThreeCorner: DS point {t1} is out of range"))?;
            ds.change_point(t2)
                .map(|tpt| tpt.set_tolerance(tpt1_tol + error))
                .ok_or_else(|| format!("PerformMoreThreeCorner: DS point {t2} is out of range"))?;

            // `ind = 0; // must be initialized because of possible use, see
            //  L2249` (`cxx:3268`) - a function-level local of the `ic` block,
            // read by `indpoint1 = ind` (`cxx:3321`) and `indice = ind`
            // (`cxx:3388`).
            let mut ind = 0i32;

            // `for (nb = 1; nb <= nbface; nb++)` (`cxx:3281-3405`).
            for nb in 1..=nbface {
                // `cxx:3283-3292`: the `Ecom(nb)` ends and its orientation.
                let ecom = self
                    .ecom
                    .get(nb - 1)
                    .cloned()
                    .ok_or_else(|| format!("PerformMoreThreeCorner: Ecom({nb}) is missing"))?;
                let (vf, vl) = edge_vertices(&ecom);
                let vf = vf.ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Ecom({nb}) has no first vertex")
                })?;
                let vl = vl.ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Ecom({nb}) has no last vertex")
                })?;
                let pf = BRepTool::vertex_point(&vf);
                let pl = BRepTool::vertex_point(&vl);
                // `para = parcom.Value(nb);
                //  Pcom = BRep_Tool::Curve(Ecom(nb), up1, up2)->Value(para);`
                // (`cxx:3288-3289`).
                let para = *self
                    .parcom
                    .get(nb - 1)
                    .ok_or_else(|| format!("PerformMoreThreeCorner: parcom({nb}) is missing"))?;
                let ecom_curve = BRepTool::edge_curve(&ecom).ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Ecom({nb}) has no 3d curve")
                })?;
                let pcom = ecom_curve.d0(para);
                let orvt = if pf.distance(&v1_point) < pl.distance(&v1_point) {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };

                if let Some(eproj) = self.eproj.get(nb - 1).cloned().flatten() {
                    // `n3d++;` (`cxx:3295`).
                    *n3d += 1;

                    // `proj = BRep_Tool::CurveOnSurface(Eproj(nb),
                    //  Face(Fproj(nb)), up1, up2); proj2d = new
                    //  Geom2d_TrimmedCurve(proj, up1, up2);` (`cxx:3297-3305`).
                    let fproj_shape = self
                        .fproj
                        .get(nb - 1)
                        .cloned()
                        .ok_or_else(|| format!("PerformMoreThreeCorner: Fproj({nb}) is missing"))?;
                    let fproj = Face::wrap(fproj_shape.clone()).ok_or_else(|| {
                        format!("PerformMoreThreeCorner: Fproj({nb}) is not a face")
                    })?;
                    let (proj, up1, up2) = curve_on_surface_range(&eproj, &fproj)
                        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
                    let proj2d: Arc<dyn Curve2d> =
                        Arc::new(Geom2dTrimmedCurve::new(proj, up1, up2));

                    // `projc = BRep_Tool::Curve(Eproj(nb), up1, up2);
                    //  cproj = new Geom_TrimmedCurve(projc, up1, up2);
                    //  pardeb = cproj->FirstParameter(); parfin =
                    //  cproj->LastParameter(); P1 = cproj->Value(pardeb);
                    //  P2 = cproj->Value(parfin);` (`cxx:3306-3314`).
                    let (c_first, c_last) = BRepTool::edge_parameters(&eproj);
                    let projc = BRepTool::edge_curve(&eproj)
                        .ok_or_else(|| "Failed to get 3D curve of edge".to_string())?;
                    let cproj = GeomTrimmedCurve::new(projc, c_first, c_last);
                    let pardeb = cproj.first_parameter();
                    let parfin = cproj.last_parameter();
                    let p1 = cproj.d0(pardeb);
                    let p2 = cproj.d0(parfin);

                    // `cxx:3316-3331`: the two point indices. `indpoint1` keeps
                    // `ind` when `P1` is not the DS point; `indpoint2` becomes
                    // the new point and is what `ind` carries afterwards.
                    let tpt1_point = *ds
                        .point(t1)
                        .ok_or_else(|| {
                            format!("PerformMoreThreeCorner: DS point {t1} is out of range")
                        })?
                        .point();
                    let tpt2_point = *ds
                        .point(t2)
                        .ok_or_else(|| {
                            format!("PerformMoreThreeCorner: DS point {t2} is out of range")
                        })?
                        .point();
                    let indpoint1 = if p1.distance(&tpt1_point) < POINT_SAME {
                        t1
                    } else {
                        ind
                    };
                    let indpoint2 = if p2.distance(&tpt2_point) < POINT_SAME {
                        t2
                    } else {
                        let added = ds.add_point(&DsPoint::with_tolerance(p2, error));
                        ind = added;
                        added
                    };

                    // `Asurf = ...; CurvOnS(Acurv, Asurf); HCons = new
                    //  Adaptor3d_CurveOnSurface(CurvOnS); Order.SetValue(n3d, 1);
                    //  Cont = new GeomPlate_CurveConstraint(HCons, Order(n3d),
                    //  10, tolapp3d, angular, 0.1); PSurf.Add(Cont);`
                    // (`cxx:3332-3342`).
                    let fproj_surface = BRepTool::face_surface(&fproj).ok_or_else(|| {
                        "PerformMoreThreeCorner: BRep_Tool::Surface is null".to_string()
                    })?;
                    let cons =
                        CurveOnSurface::new(Arc::clone(&proj2d), Arc::clone(&fproj_surface));
                    self.order[*n3d] = 1;
                    let cont = CurveConstraint::from_boundary(
                        BoundaryCurve::OnSurface(Arc::new(cons)),
                        1,
                        NBCURVPNT,
                        tolapp3d,
                        angular,
                        0.1,
                    )?;
                    psurf.add_curve_constraint(cont);

                    // `TopOpeBRepDS_Curve tcurv3d(cproj, error);
                    //  indcurve3d.SetValue(n3d, DStr.AddCurve(tcurv3d));
                    //  Interfp1 = ChFi3d_FilPointInDS(FORWARD, ..., indpoint1,
                    //  pardeb); Interfp2 = ...REVERSED..., indpoint2, parfin);
                    //  DStr.ChangeCurveInterferences(...).Append(...)` twice
                    // (`cxx:3343-3348`).
                    let cproj: Arc<dyn Curve> = Arc::new(cproj);
                    let indcurve = ds.add_curve(&mut DsCurve::with_curve(
                        Some(Arc::clone(&cproj)),
                        error,
                        false,
                    ));
                    self.indcurve3d[*n3d] = indcurve;
                    let interfp1 = ch_fi3d_fil_point_in_ds(
                        Orientation::Forward,
                        indcurve,
                        indpoint1,
                        pardeb,
                        false,
                    );
                    let interfp2 = ch_fi3d_fil_point_in_ds(
                        Orientation::Reversed,
                        indcurve,
                        indpoint2,
                        parfin,
                        false,
                    );
                    ds.change_curve_interferences(indcurve).push(interfp1);
                    ds.change_curve_interferences(indcurve).push(interfp2);

                    // `num = DStr.AddShape(Fproj(nb));` (`cxx:3349`).
                    let num = ds.add_shape(&fproj_shape);

                    // `cxx:3350-3358`: the orientation of `Ecom(nb)` in the
                    // forward `Fproj(nb)`.
                    let mut oredge = Orientation::Forward;
                    for child in cumulated_children(&fproj_shape.oriented(Orientation::Forward)) {
                        if child.shape_type() == crate::abs::ShapeType::Edge
                            && is_same(&child, &ecom.0)
                        {
                            oredge = child.orientation();
                            break;
                        }
                    }

                    // `cxx:3360-3384`: `orinterf` from whether the trimmed
                    // start is `Pcom`.
                    let orinterf = if p1.distance(&pcom) > CURVE_TOL {
                        if orvt.is_forward() {
                            oredge
                        } else {
                            oredge.reversed()
                        }
                    } else if orvt.is_forward() {
                        oredge.reversed()
                    } else {
                        oredge
                    };

                    // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d), num,
                    //  proj2d, orinterf);
                    //  DStr.ChangeShapeInterferences(num).Append(Interfc);`
                    // (`cxx:3385-3386`).
                    let interfc = ch_fi3d_fil_curve_in_ds(indcurve, num, Some(proj2d), orinterf);
                    ds.change_shape_interferences_at(num).push(interfc);
                }

                // `indice = ind; if (nb != nbface) { ... }` (`cxx:3388-3402`).
                let mut indice = ind;
                if nb != nbface {
                    if self.eproj.get(nb - 1).cloned().flatten().is_none() {
                        indice = self.indpoint[ic][1];
                    }
                    if self.eproj.get(nb).cloned().flatten().is_none() {
                        indice = self.indpoint[icplus][0];
                    }
                    self.indice[*n3d] = indice;

                    // `int Iarc1 = DStr.AddShape(TopoDS::Edge(Ecom(nb)));
                    //  Interfp1 = ChFi3d_FilPointInDS(orvt, Iarc1, indice,
                    //  parcom.Value(nb));
                    //  DStr.ChangeShapeInterferences(Iarc1).Append(Interfp1);`
                    // (`cxx:3399-3401`).
                    let iarc1 = ds.add_shape(&ecom.0);
                    let interfp1 = ch_fi3d_fil_point_in_ds(orvt, iarc1, indice, para, false);
                    ds.change_shape_interferences_at(iarc1).push(interfp1);
                }
            }
        }

        Ok(())
    }

    /// `cxx:3408-3488`: the two tangent free borders.
    pub(super) fn plate_droit(
        &mut self,
        ds: &mut DataStructure,
        psurf: &mut BuildPlateSurface,
        n3d: &mut usize,
        tolapp3d: f64,
        angular: f64,
    ) -> Result<(), String> {
        // `if (droit)` (`cxx:3408`).
        if !self.droit {
            return Ok(());
        }
        let nedge = self.nedge as usize;

        // `for (ic = 0; ic < nedge; ic++)` (`cxx:3410-3487`).
        for ic in 0..nedge {
            let (icplus, icmoins) = indices(nedge, ic);
            let ecur = self.edge_at(ic)?;
            // `if (ecur.IsSame(edgelibre1) || ecur.IsSame(edgelibre2))`
            // (`cxx:3419`).
            if !self.is_free_border(&ecur) {
                continue;
            }
            *n3d += 1;

            // `curve2d = BRep_Tool::CurveOnSurface(Evive(ic),
            //  Face(Fvive(ic, icplus)), ufirst, ulast);
            //  curve = BRep_Tool::Curve(Evive(ic), ufirst, ulast);`
            // (`cxx:3422-3426`).
            let face = self.fvive_face(ic, icplus)?.ok_or_else(|| {
                format!("PerformMoreThreeCorner: Fvive({ic}, {icplus}) is null")
            })?;
            let (curve2d, _pfirst, _plast) = curve_on_surface_range(&ecur, &face)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            let (ufirst, ulast) = BRepTool::edge_parameters(&ecur);
            let curve = BRepTool::edge_curve(&ecur)
                .ok_or_else(|| "Failed to get 3D curve of edge".to_string())?;
            let p_icmoins = self.p[ic][icmoins];

            // `cxx:3427-3442`: the trimmed border, then the ends. `V1` sits
            // either first (`ctrim` keeps `ufirst -> p(ic, icmoins)`) or last.
            let (first, last) = edge_vertices(&ecur);
            let first_is_v1 = first.as_ref().is_some_and(|v| is_same(&v.0, &self.vtx));
            let last_is_v1 = last.as_ref().is_some_and(|v| is_same(&v.0, &self.vtx));

            let (mut ctrim, mut ctrim2d, mut indpoint1, mut indpoint2, mut isvt1, mut isvt2);
            if first_is_v1 {
                ctrim = GeomTrimmedCurve::new(Arc::clone(&curve), ufirst, p_icmoins);
                ctrim2d = Geom2dTrimmedCurve::new(Arc::clone(&curve2d), ufirst, p_icmoins);
                indpoint1 = ds.add_shape(&self.vtx);
                indpoint2 = self.indpoint[ic][1];
                isvt1 = true;
                isvt2 = false;
            } else {
                ctrim = GeomTrimmedCurve::new(Arc::clone(&curve), p_icmoins, ulast);
                ctrim2d = Geom2dTrimmedCurve::new(Arc::clone(&curve2d), p_icmoins, ulast);
                indpoint2 = ds.add_shape(&self.vtx);
                indpoint1 = self.indpoint[ic][1];
                isvt1 = false;
                isvt2 = true;
            }

            // `cxx:3443-3464`: `libre(ic)` keeps the first-vertex arm, the
            // other case keeps the last-vertex one; both flip the border and
            // move the `V1` end to index 2 (the `!libre` arm to index 1).
            if self.libre[ic] {
                if first_is_v1 {
                    ctrim.reverse();
                    ctrim2d.reverse();
                    indpoint2 = ds.add_shape(&self.vtx);
                    isvt2 = true;
                    isvt1 = false;
                    indpoint1 = self.indpoint[ic][1];
                }
            } else if last_is_v1 {
                ctrim.reverse();
                ctrim2d.reverse();
                indpoint1 = ds.add_shape(&self.vtx);
                isvt1 = true;
                isvt2 = false;
                indpoint2 = self.indpoint[ic][1];
            }

            // `ufirst = ctrim->FirstParameter();
            //  ulast = ctrim->LastParameter();` (`cxx:3466-3467`).
            let ufirst = ctrim.first_parameter();
            let ulast = ctrim.last_parameter();

            // `Asurf = ...; CurvOnS(Acurv, Asurf); HCons = new
            //  Adaptor3d_CurveOnSurface(CurvOnS); Order.SetValue(n3d, 0);
            //  Cont = new GeomPlate_CurveConstraint(HCons, Order(n3d), 10,
            //  tolapp3d, angular, 0.1); PSurf.Add(Cont);` (`cxx:3469-3477`).
            let face_surface = BRepTool::face_surface(&face)
                .ok_or_else(|| "PerformMoreThreeCorner: BRep_Tool::Surface is null".to_string())?;
            let ctrim2d: Arc<dyn Curve2d> = Arc::new(ctrim2d);
            let cons = CurveOnSurface::new(Arc::clone(&ctrim2d), Arc::clone(&face_surface));
            self.order[*n3d] = 0;
            let cont = CurveConstraint::from_boundary(
                BoundaryCurve::OnSurface(Arc::new(cons)),
                0,
                NBCURVPNT,
                tolapp3d,
                angular,
                0.1,
            )?;
            psurf.add_curve_constraint(cont);

            // `TopOpeBRepDS_Curve tcurv3d(ctrim, 1.e-4);
            //  indcurve3d.SetValue(n3d, DStr.AddCurve(tcurv3d));
            //  Interfp1 = ChFi3d_FilPointInDS(FORWARD, ..., indpoint1, ufirst,
            //  isvt1); Interfp2 = ...REVERSED..., indpoint2, ulast, isvt2);
            //  DStr.ChangeCurveInterferences(...).Append(...)` twice
            // (`cxx:3478-3485`).
            let ctrim: Arc<dyn Curve> = Arc::new(ctrim);
            let indcurve =
                ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&ctrim)), CURVE_TOL, false));
            self.indcurve3d[*n3d] = indcurve;
            let interfp1 = ch_fi3d_fil_point_in_ds(
                Orientation::Forward,
                indcurve,
                indpoint1,
                ufirst,
                isvt1,
            );
            let interfp2 = ch_fi3d_fil_point_in_ds(
                Orientation::Reversed,
                indcurve,
                indpoint2,
                ulast,
                isvt2,
            );
            ds.change_curve_interferences(indcurve).push(interfp1);
            ds.change_curve_interferences(indcurve).push(interfp2);
        }

        Ok(())
    }
}
