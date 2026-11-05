//! The first two stages of the `GeomPlate` section of
//! `PerformMoreThreeCorner`: the p-curve constraints of the non-sharp stripes
//! and the DS point indices of the living edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx`
//! (OCCT 8.0.0):
//!
//! - `cxx:2716-2852`: for every non-sharp slot, the straight 2D border between
//!   its two common points, the `GeomPlate_CurveConstraint` added to `PSurf`,
//!   the 3D curve built from it, its two DS points (shared with the previous
//!   slot when the coordinates coincide within `1.e-4`) and the write back into
//!   the stripe's surf data,
//! - `cxx:2854-2905`: for every sharp slot, the DS point indices of the living
//!   edge at the corner, taken from the neighbouring stripes' common points
//!   when they coincide with the edge point within `max(1.e-4, tolcp)`.
//!
//! Translation notes:
//! - `SurfIndex(CD, ic, i(ic, icplus), ChFiSURFACE)` is [`surf_index`] with
//!   [`SurfType::ChFiSurface`],
//! - `DStr.Surface(indice).Surface()` is read through the DS surface table; a
//!   null handle is the `Err` of the stage, which is what the `GeomAdaptor_Surface`
//!   construction would raise on,
//! - `GeomLib::BuildCurve3d(tolapp, CurvOnS, First, Last, Curv3d, maxapp,
//!   avedev)` keeps every default of `GeomLib.hxx:85-94` (`GeomAbs_C1`,
//!   `MaxDegree = 14`, `MaxSegment = 30`),
//! - `CD.Value(ic)->ChangePCurve(isfirst) = pcurve` is
//!   [`crate::chfi_ds::Stripe::change_p_curve`],
//! - the OCCT `point1.Distance(tpt.Point()) < 1.e-4` scans keep the *loop*
//!   index `ii` and then read `indpoint.Value(ii - 1, ...)` (`cxx:2792`,
//!   `:2812`), so the port does the same and answers an error when that entry
//!   would be out of range.

use std::sync::Arc;

use occt_core::gp::GpDir2d;
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_lib::geom_lib_build_curve3d;
use occt_geom::geom_plate::{BoundaryCurve, BuildPlateSurface, CurveConstraint};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::CommonPoint;
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure, Point as DsPoint};

use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;

/// `double tolapp = 1.e-4` (`cxx:1313`).
const TOLAPP: f64 = 1.0e-4;

/// `int nbcurvpnt = 10` (`cxx:2713`).
const NBCURVPNT: i32 = 10;

/// `GeomLib::BuildCurve3d`'s default `Continuity` (`GeomLib.hxx:92`).
const BUILD_CONTINUITY: Shape = Shape::C1;
/// `GeomLib::BuildCurve3d`'s default `MaxDegree` (`GeomLib.hxx:93`).
const BUILD_MAX_DEGREE: i32 = 14;
/// `GeomLib::BuildCurve3d`'s default `MaxSegment` (`GeomLib.hxx:94`).
const BUILD_MAX_SEGMENT: i32 = 30;

/// `1.e-4`, the point-coincidence tolerance of `cxx:2782`, `:2803`, `:2877`
/// and `:2889`.
const SAME_POINT: f64 = 1.0e-4;

/// `ChFiDS_CommonPoint cp; cp.SetPoint(P);` — the OCCT two-step form, kept
/// because `ChFiDS_CommonPoint` has no point-only constructor.
fn point_only_cp(p: &occt_core::gp::GpPnt) -> CommonPoint {
    let mut cp = CommonPoint::new();
    cp.set_point(p);
    cp
}

impl MoreThreeCorner {
    /// `cxx:2716-2852`: the p-curve constraints of the non-sharp stripes.
    ///
    /// `n3d` is the function-scope curve counter (`cxx:1303`), shared with the
    /// stages that follow. `tolapp3d` / `angular` are the builder's
    /// `myTolApp3d` / `myAngular`.
    pub(super) fn plate_stripe_constraints(
        &mut self,
        ds: &mut DataStructure,
        psurf: &mut BuildPlateSurface,
        n3d: &mut usize,
        tolapp3d: f64,
        angular: f64,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;

        // `for (ic = 0; ic < nedge; ic++)` (`cxx:2719-2852`).
        for ic in 0..nedge {
            if self.sharp[ic] {
                continue;
            }
            // `n3d++; Indices(nedge, ic, icplus, icmoins);
            //  jfp = 3 - jf.Value(ic);` (`cxx:2723-2725`).
            *n3d += 1;
            let (icplus, icmoins) = indices(nedge, ic);
            self.jfp = 3 - self.jf[ic];

            // `Calcul_P2dOnSurf(CD.Value(ic), jfp, i(ic, icmoins),
            //  p(ic, icmoins), p2d1); Calcul_P2dOnSurf(CD.Value(ic), jf(ic),
            //  i(ic, icplus), p(ic, icplus), p2d2);` (`cxx:2726-2728`).
            let stripe = self.stripe_at(ic)?.clone();
            let p2d1 = {
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                super::calcul::calcul_p2d_on_surf(
                    &guard,
                    self.jfp,
                    self.i[ic][icmoins],
                    self.p[ic][icmoins],
                )?
            };
            let p2d2 = {
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                super::calcul::calcul_p2d_on_surf(
                    &guard,
                    self.jf[ic],
                    self.i[ic][icplus],
                    self.p[ic][icplus],
                )?
            };

            // `indice = SurfIndex(CD, ic, i(ic, icplus), ChFiSURFACE);
            //  Asurf = new GeomAdaptor_Surface(DStr.Surface(indice).Surface());`
            // (`cxx:2730-2732`).
            let indice =
                surf_index(&self.cd, ic, self.i[ic][icplus] as usize, SurfType::ChFiSurface)?;
            let asurf = ds
                .surface(indice)
                .surface()
                .cloned()
                .ok_or_else(|| format!("PerformMoreThreeCorner: DStr.Surface({indice}) is null"))?;

            // `xdir = p2d2.X() - p2d1.X(); ydir = ...; l0 = sqrt(...);
            //  dir(xdir, ydir); l = new Geom2d_Line(p2d1, dir);
            //  pcurve = new Geom2d_TrimmedCurve(l, 0, l0);` (`cxx:2734-2739`).
            let xdir = p2d2.x() - p2d1.x();
            let ydir = p2d2.y() - p2d1.y();
            let l0 = (xdir * xdir + ydir * ydir).sqrt();
            let dir = GpDir2d::new(xdir, ydir)
                .map_err(|_| "PerformMoreThreeCorner: gp_Dir2d on a null border".to_string())?;
            let line = Geom2dLine::from_pnt_dir(p2d1, dir);
            let pcurve: Arc<dyn Curve2d> =
                Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), 0.0, l0));

            // `Order.SetValue(ic, constr)` with `constr = 1` (`cxx:2714`,
            // `:2743`).
            let order_ic = 1;
            self.order[ic] = order_ic;

            // `Cont = new GeomPlate_CurveConstraint(HCons, Order(ic), nbcurvpnt,
            //  tolapp3d, angular, 0.1); PSurf.Add(Cont);` (`cxx:2745-2748`).
            let cons = CurveOnSurface::new(Arc::clone(&pcurve), Arc::clone(&asurf));
            let cont = CurveConstraint::from_boundary(
                BoundaryCurve::OnSurface(Arc::new(cons.clone())),
                order_ic,
                NBCURVPNT,
                tolapp3d,
                angular,
                0.1,
            )?;
            psurf.add_curve_constraint(cont);

            // `isfirst = (sens.Value(ic) == 1); GeomLib::BuildCurve3d(...)`
            // (`cxx:2750-2757`).
            let isfirst = self.sens[ic] == 1;
            let built = geom_lib_build_curve3d(
                TOLAPP,
                &cons,
                &asurf,
                &pcurve,
                cons.first_parameter(),
                cons.last_parameter(),
                BUILD_CONTINUITY,
                BUILD_MAX_DEGREE,
                BUILD_MAX_SEGMENT,
            );
            let maxapp = built.max_deviation;
            let curv3d = built.curve.ok_or_else(|| {
                "PerformMoreThreeCorner: GeomLib::BuildCurve3d is null".to_string()
            })?;

            // `TopOpeBRepDS_Curve tcurv3d(Curv3d, maxapp);
            //  indcurve3d.SetValue(n3d, DStr.AddCurve(tcurv3d));`
            // (`cxx:2758-2759`).
            let icurv =
                ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&curv3d)), maxapp, false));
            self.indcurve3d[*n3d] = icurv;

            // `point1 = CurvOnS.Value(FirstParameter()); point2 =
            //  CurvOnS.Value(LastParameter());` (`cxx:2760-2762`).
            let point1 = cons.d0(cons.first_parameter());
            let point2 = cons.d0(cons.last_parameter());
            self.errapp[ic] = maxapp;

            // `cxx:2766-2818`: the two DS points, shared with the previous
            // non-sharp slot when they coincide.
            let (ip0, ip1) = if ic == 0 {
                let i0 = ds.add_point(&DsPoint::with_tolerance(point1, maxapp));
                let i1 = ds.add_point(&DsPoint::with_tolerance(point2, maxapp));
                (i0, i1)
            } else {
                // `trouve = false; for (ii = 0; ii < ic && !trouve; ii++)
                //  if (!sharp(ii)) { tpt = DStr.ChangePoint(indpoint(ii, 1));
                //  if (point1.Distance(tpt.Point()) < 1.e-4) trouve = true; }`
                // (`cxx:2776-2785`).
                let mut found_at: Option<usize> = None;
                for ii in 0..ic {
                    if self.sharp[ii] {
                        continue;
                    }
                    let idx = self.indpoint[ii][1];
                    let tpt = ds.point(idx).ok_or_else(|| {
                        format!("PerformMoreThreeCorner: DS point {idx} is out of range")
                    })?;
                    if point1.distance(tpt.point()) < SAME_POINT {
                        found_at = Some(ii);
                        break;
                    }
                }
                let first = match found_at {
                    // `indpoint.SetValue(ic, 0, indpoint.Value(ii - 1, 1));`
                    // (`cxx:2791-2793`).
                    Some(ii) => {
                        let prev = ii.checked_sub(1).ok_or_else(|| {
                            "PerformMoreThreeCorner: indpoint(ii - 1, 1) is out of range"
                                .to_string()
                        })?;
                        self.indpoint[prev][1]
                    }
                    None => ds.add_point(&DsPoint::with_tolerance(point1, maxapp)),
                };

                // `trouve = false; for (ii = 0; ii < ic && !trouve; ii++)
                //  if (!sharp(ii)) { tpt = DStr.ChangePoint(indpoint(ii, 0));
                //  if (point2.Distance(tpt.Point()) < 1.e-4) trouve = true; }`
                // (`cxx:2798-2808`).
                let mut found_at: Option<usize> = None;
                for ii in 0..ic {
                    if self.sharp[ii] {
                        continue;
                    }
                    let idx = self.indpoint[ii][0];
                    let tpt = ds.point(idx).ok_or_else(|| {
                        format!("PerformMoreThreeCorner: DS point {idx} is out of range")
                    })?;
                    if point2.distance(tpt.point()) < SAME_POINT {
                        found_at = Some(ii);
                        break;
                    }
                }
                let second = match found_at {
                    // `indpoint.SetValue(ic, 1, indpoint.Value(ii - 1, 0));`
                    // (`cxx:2811-2813`).
                    Some(ii) => {
                        let prev = ii.checked_sub(1).ok_or_else(|| {
                            "PerformMoreThreeCorner: indpoint(ii - 1, 0) is out of range"
                                .to_string()
                        })?;
                        self.indpoint[prev][0]
                    }
                    None => ds.add_point(&DsPoint::with_tolerance(point2, maxapp)),
                };
                (first, second)
            };
            self.indpoint[ic][0] = ip0;
            self.indpoint[ic][1] = ip1;

            // `cxx:2820-2830`: the stripe update.
            let isurf1 = 3 - self.jf[ic];
            let isurf2 = self.jf[ic];
            let indcurve = self.indcurve3d[*n3d];
            {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                if isurf1 == 2 {
                    guard.set_orientation(Orientation::Reversed, isfirst);
                }
                guard.set_curve(indcurve, isfirst);
                guard.set_index_point(self.indpoint[ic][0], isfirst, isurf1);
                guard.set_index_point(self.indpoint[ic][1], isfirst, isurf2);
                guard.set_parameters(isfirst, pcurve.first_parameter(), pcurve.last_parameter());
            }

            // `cp1.SetPoint(point1); cp2.SetPoint(point2);
            //  ...->Value(i(ic, icmoins))->ChangeVertex(isfirst, isurf1) = cp1;
            //  ...->ChangeVertex(isfirst, isurf2) = cp2;`
            // (`cxx:2831-2839`).
            let idx_sd = self.i[ic][icmoins];
            let cp1 = point_only_cp(&point1);
            let cp2 = point_only_cp(&point2);
            let p_icmoins = self.p[ic][icmoins];
            let p_icplus = self.p[ic][icplus];
            {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                let seq = guard.change_set_of_surf_data();
                let pos = usize::try_from(idx_sd)
                    .ok()
                    .and_then(|i| i.checked_sub(1))
                    .filter(|i| *i < seq.len())
                    .ok_or_else(|| {
                        format!("PerformMoreThreeCorner: surf data index {idx_sd} is out of range")
                    })?;
                let sd = Arc::make_mut(&mut seq[pos]);
                *sd.vertex_mut(isfirst, isurf1) = cp1;
                *sd.vertex_mut(isfirst, isurf2) = cp2;
                // `ChangeInterference(isurf1).SetParameter(p(ic, icmoins),
                //  isfirst); ChangeInterference(isurf2).SetParameter(
                //  p(ic, icplus), isfirst);` (`cxx:2840-2848`).
                sd.interference_mut(isurf1)
                    .set_parameter(p_icmoins, isfirst);
                sd.interference_mut(isurf2)
                    .set_parameter(p_icplus, isfirst);
            }

            // `CD.Value(ic)->ChangePCurve(isfirst) = pcurve;` (`cxx:2849`).
            {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                *guard.change_p_curve(isfirst) = Some(pcurve);
            }
        }

        Ok(())
    }

    /// `cxx:2854-2905`: the DS point indices of the living edges.
    ///
    /// `tolcp` is the largest common-point tolerance the corner accepted
    /// (`cxx:2201-2315`).
    pub(super) fn plate_living_point_indices(
        &mut self,
        ds: &mut DataStructure,
        tolcp: f64,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;

        // `for (ic = 0; ic < nedge; ic++)` (`cxx:2857-2905`).
        for ic in 0..nedge {
            if !self.sharp[ic] {
                continue;
            }
            let (icplus, icmoins) = indices(nedge, ic);

            // `BRepAdaptor_Curve C(TopoDS::Edge(Evive.Value(ic)));
            //  PE = C.Value(p(ic, icplus));
            //  TopOpeBRepDS_Point TPE(PE, BRep_Tool::Tolerance(...));`
            // (`cxx:2859-2861`).
            let edge = self.edge_at(ic)?;
            let curve = crate::topopebrep::tool::adaptor::CurveAdaptor::of_edge(&edge);
            let pe = curve
                .d0(self.p[ic][icplus])
                .ok_or_else(|| "PerformMoreThreeCorner: Evive has no 3d curve".to_string())?;
            let tpe = DsPoint::with_tolerance(pe, BRepTool::edge_tolerance(&edge));

            // `cxx:2864-2868`: with two fillets on the vertex, `V1` is the edge
            // point.
            if self.deuxconges {
                let vtx = self.vtx.clone();
                let ivtx = ds.add_shape(&vtx);
                self.ivtx = ivtx;
                self.indpoint[ic][0] = ivtx;
                self.indpoint[ic][1] = ivtx;
            }

            // `if (!sharp(icplus)) { isfirst = (sens(icplus) == 1);
            //  jfp = 3 - jf(icplus); cp = CD(icplus)->...->ChangeVertex(isfirst,
            //  jfp); if (cp.Point().Distance(PE) <= max(1.e-4, tolcp)) {
            //  indpoint(ic, 0) = indpoint(ic, 1) = indpoint(icplus, 0); } }`
            // (`cxx:2869-2881`).
            if !self.sharp[icplus] {
                let isfirst = self.sens[icplus] == 1;
                let jfp = 3 - self.jf[icplus];
                let cp = self.common_point_at(icplus, isfirst, jfp)?;
                if cp.point().distance(&pe) <= SAME_POINT.max(tolcp) {
                    self.indpoint[ic][0] = self.indpoint[icplus][0];
                    self.indpoint[ic][1] = self.indpoint[icplus][0];
                }
            }

            // `if (!sharp(icmoins)) { isfirst = (sens(icmoins) == 1);
            //  cp = CD(icmoins)->...->ChangeVertex(isfirst, jf(icmoins));
            //  if (cp.Point().Distance(PE) <= max(1.e-4, tolcp)) {
            //  if (indpoint(ic, 0) == 0) { indpoint(ic, 0) =
            //  indpoint(ic, 1) = indpoint(icmoins, 1); } } }` (`cxx:2882-2896`).
            if !self.sharp[icmoins] {
                let isfirst = self.sens[icmoins] == 1;
                let cp = self.common_point_at(icmoins, isfirst, self.jf[icmoins])?;
                if cp.point().distance(&pe) <= SAME_POINT.max(tolcp) && self.indpoint[ic][0] == 0 {
                    self.indpoint[ic][0] = self.indpoint[icmoins][1];
                    self.indpoint[ic][1] = self.indpoint[icmoins][1];
                }
            }

            // `if (indpoint(ic, 0) == 0) { indpoint(ic, 0) =
            //  DStr.AddPoint(TPE); indpoint(ic, 1) = indpoint(ic, 0); }`
            // (`cxx:2898-2902`).
            if self.indpoint[ic][0] == 0 {
                let ip = ds.add_point(&tpe);
                self.indpoint[ic][0] = ip;
                self.indpoint[ic][1] = ip;
            }
        }

        Ok(())
    }
}
