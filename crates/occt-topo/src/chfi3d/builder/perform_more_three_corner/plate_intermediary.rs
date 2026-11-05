//! The intermediary connecting curves of the `GeomPlate` section: the curves
//! that bridge two stripes when the corner has no intersection.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:2907-3260`
//! (OCCT 8.0.0).
//!
//! For every slot that is neither `oksea` nor `moresurf` nor `libre`, the stage
//! takes the two 2D borders meeting at the corner (`BRep_Tool::CurveOnSurface`
//! on a sharp slot, `Calcul_C2dOnFace` on a stripe), decides how to join them -
//! a straight line, one of the `proj2d2` / `proj2d1` projected curves, or a
//! `FairCurve_Batten` - adds the resulting `GeomPlate_CurveConstraint`, stores
//! the 3D curve in the DS with its two point interferences and, finally, the
//! curve interference on the forward face.
//!
//! Translation notes:
//! - `pcurve`, `curveint`, `maxapp1` and `raccordbatten` are per-slot locals in
//!   OCCT; the port declares them inside the loop body for the same reason,
//! - `Geom2dAdaptor_Curve(pcurve)` + `Adaptor3d_CurveOnSurface(Acurv, Asurf)` is
//!   [`CurveOnSurface::new`]: the adaptor only wraps the pcurve and the surface,
//! - `Order.SetValue(n3d, ...)` is the 1-based `Order` array of the state, the
//!   same array the non-sharp stage writes at `cxx:2744`,
//! - `BRep_Tool::CurveOnSurface(E, F, u1, u2)` is
//!   [`curve_on_surface_range`]; a missing p-curve is the
//!   `Standard_ConstructionError("Failed to get p-curve of edge")` of
//!   `cxx:2951-2954` / `:2971-2974`,
//! - `tpt1.Tolerance(tpt1.Tolerance() + maxapp1)` is
//!   [`DataStructure::change_point`] followed by `set_tolerance`,
//! - `couture` reframing (`cxx:2990-3024`) mutates only the U of the two UV
//!   points; the port keeps the same two `set_x` calls.

use std::sync::Arc;

use occt_geom::approx_same_parameter::{u_resolution, v_resolution};use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_lib::geom_lib_build_curve3d;
use occt_geom::geom_plate::{BoundaryCurve, BuildPlateSurface, CurveConstraint};
use occt_geom::{Curve, GeomTrimmedCurve};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure, Point as DsPoint};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::plate_helpers::{
    calcul_batten, calcul_droite, orientation_arete_vive_consecutive, orientation_ic_non_vive,
    orientation_icplus_non_vive,
};
use super::surf_index::indices;
use super::MoreThreeCorner;

/// `double tolapp = 1.e-4` (`cxx:1313`).
const TOLAPP: f64 = 1.0e-4;
/// `GeomPlate_CurveConstraint(..., 10, ...)` (`cxx:3170`).
const NBCURVPNT: i32 = 10;
/// The `1.e-6` the projection arms store in `maxapp1` (`cxx:3072`, `:3084`).
const PROJ_TOL: f64 = 1.0e-6;

/// `GeomLib::BuildCurve3d`'s defaults (`GeomLib.hxx:92-94`).
const BUILD_CONTINUITY: occt_core::kernel::geomabs::Shape = occt_core::kernel::geomabs::Shape::C1;
const BUILD_MAX_DEGREE: i32 = 14;
const BUILD_MAX_SEGMENT: i32 = 30;

/// `Calcul_C2dOnFace(CD, jfposit, indice, c2d)` (`cxx:187-201`).
fn calcul_c2d_on_face(
    stripe: &crate::chfi_ds::Stripe,
    jfposit: i32,
    indice: i32,
) -> Result<Option<Arc<dyn Curve2d>>, String> {
    if indice < 1 {
        return Err(format!(
            "PerformMoreThreeCorner: surf data index {indice} is out of range"
        ));
    }
    let data = stripe
        .set_of_surf_data()
        .get(indice as usize - 1)
        .ok_or_else(|| {
            format!("PerformMoreThreeCorner: surf data index {indice} is out of range")
        })?;
    let fi = if jfposit == 1 {
        data.interference_on_s1()
    } else {
        data.interference_on_s2()
    };
    Ok(fi.p_curve_on_face().cloned())
}

impl MoreThreeCorner {
    /// `cxx:2907-3260`: the intermediary connecting curves.
    pub(super) fn plate_intermediary_curves(
        &mut self,
        ds: &mut DataStructure,
        psurf: &mut BuildPlateSurface,
        n3d: &mut usize,
        nconges: usize,
        tolapp3d: f64,
        angular: f64,
    ) -> Result<(), String> {
        // `bool raccordbatten; if (!inters) { ... }` (`cxx:2910-2912`).
        if self.inters {
            return Ok(());
        }
        let nedge = self.nedge as usize;

        // `for (ic = 0; ic < nedge; ic++)` (`cxx:2914-3259`).
        for ic in 0..nedge {
            // `if (!oksea(ic) && !moresurf(ic) && !libre(ic))` (`cxx:2917`).
            if self.oksea[ic] || self.moresurf[ic] || self.libre[ic] {
                continue;
            }
            let (icplus, icmoins) = indices(nedge, ic);

            // `raccordbatten = false; if (!regul(ic)) { raccordbatten = true;
            //  if (regul(icplus)) raccordbatten = false; }` (`cxx:2920-2927`).
            let mut raccordbatten = !self.regul[ic] && !self.regul[icplus];

            *n3d += 1;

            // `cxx:2936-2953`: the first 2D border.
            let c1_opt: Option<Arc<dyn Curve2d>> = if self.sharp[ic] {
                let edge = self.edge_at(ic)?;
                let face = self.fvive_face(ic, icplus)?.ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Fvive({ic}, {icplus}) is null")
                })?;
                curve_on_surface_range(&edge, &face).map(|(c, _u1, _u2)| c)
            } else {
                let stripe = self.stripe_at(ic)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                calcul_c2d_on_face(&guard, self.jf[ic], self.i[ic][icplus])?
            };
            // `if (curv2d1.IsNull()) throw Standard_ConstructionError(
            //  "Failed to get p-curve of edge");` (`cxx:2950-2953`).
            let c1 = c1_opt.ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            let mut p2d1 = c1.d0(self.p[ic][icplus]);

            // `cxx:2956-2970`: the second 2D border.
            let curv2d2: Arc<dyn Curve2d> = if self.sharp[icplus] {
                let edge = self.edge_at(icplus)?;
                let face = self.fvive_face(ic, icplus)?.ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Fvive({ic}, {icplus}) is null")
                })?;
                curve_on_surface_range(&edge, &face)
                    .map(|(c, _u1, _u2)| c)
                    .ok_or_else(|| "Failed to get p-curve of edge".to_string())?
            } else {
                self.jfp = 3 - self.jf[icplus];
                let stripe = self.stripe_at(icplus)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                calcul_c2d_on_face(&guard, self.jfp, self.i[icplus][ic])?
                    .ok_or_else(|| "Failed to get p-curve of edge".to_string())?
            };
            let mut p2d2 = curv2d2.d0(self.p[icplus][ic]);

            // `Asurf = new GeomAdaptor_Surface(BRep_Tool::Surface(Fvive));
            //  tolu = Asurf->UResolution(1.e-3); tolv = VResolution(1.e-3);
            //  ratio = max / min` (`cxx:2976-2988`).
            let face = self.fvive_face(ic, icplus)?.ok_or_else(|| {
                format!("PerformMoreThreeCorner: Fvive({ic}, {icplus}) is null")
            })?;
            let face_surface = BRepTool::face_surface(&face)
                .ok_or_else(|| "PerformMoreThreeCorner: BRep_Tool::Surface is null".to_string())?;
            let _asurf = SurfaceAdaptor::of_surface(Arc::clone(&face_surface));
            let tolu = u_resolution(face_surface.as_ref(), 1.0e-3);
            let tolv = v_resolution(face_surface.as_ref(), 1.0e-3);
            let ratio = if tolu > tolv { tolu / tolv } else { tolv / tolu };

            // `cxx:2990-3024`: the sewing-edge reframing.
            if self.couture {
                let pi1 = (0.0..=std::f64::consts::PI).contains(&p2d1.x());
                let pi2 = (0.0..=std::f64::consts::PI).contains(&p2d2.x());
                if let Some(edge_couture) = self.edge_couture.as_ref() {
                    let edge_ic = self.edge_at(ic)?;
                    if is_same(&edge_ic.0, &edge_couture.0) {
                        let mut xx = p2d1.x();
                        if pi2 && !pi1 {
                            xx -= 2.0 * std::f64::consts::PI;
                        }
                        if !pi2 && pi1 {
                            xx += 2.0 * std::f64::consts::PI;
                        }
                        p2d1.set_x(xx);
                    }
                    let edge_icplus = self.edge_at(icplus)?;
                    if is_same(&edge_icplus.0, &edge_couture.0) {
                        let mut xx = p2d2.x();
                        if pi2 && !pi1 {
                            xx += 2.0 * std::f64::consts::PI;
                        }
                        if !pi2 && pi1 {
                            xx -= 2.0 * std::f64::consts::PI;
                        }
                        p2d2.set_x(xx);
                    }
                }
            }

            // `xdir = p2d2.X() - p2d1.X(); ydir = ...;
            //  l0 = sqrt(...)` (`cxx:3025-3028`).
            let xdir = p2d2.x() - p2d1.x();
            let ydir = p2d2.y() - p2d1.y();
            let l0 = (xdir * xdir + ydir * ydir).sqrt();

            // `if (l0 < 1.e-7 || ponctuel(ic))` (`cxx:3029-3048`).
            if l0 < 1.0e-7 || self.ponctuel[ic] {
                *n3d -= 1;
                self.ponctuel[ic] = true;
                if !self.deuxconges {
                    if self.sharp[icplus] && self.indpoint[icplus][0] == 0 {
                        self.indpoint[icplus][0] = self.indpoint[ic][1];
                        self.indpoint[icplus][1] = self.indpoint[ic][1];
                    }
                    if self.sharp[ic] && self.indpoint[ic][0] == 0 {
                        self.indpoint[ic][0] = self.indpoint[icmoins][1];
                        self.indpoint[ic][1] = self.indpoint[icmoins][1];
                    }
                }
                continue;
            }

            // `cxx:3050-3088`: the choice between line, projection and batten.
            let mut pcurve_opt: Option<Arc<dyn Curve2d>> = None;
            let mut curveint: Option<Arc<dyn Curve>> = None;
            let mut maxapp1 = 0.0f64;

            if ratio > 10.0 && nconges == 1 {
                raccordbatten = true;
            }
            if ratio > 10.0 && raccordbatten {
                pcurve_opt = Some(calcul_droite(&p2d1, xdir, ydir)?);
                raccordbatten = false;
            } else if !raccordbatten {
                // `cxx:3060-3086`: the projected curves are returned.
                if self.regul[ic] {
                    match (self.proj2d2[ic].clone(), self.cproj2[ic].clone()) {
                        (Some(pc), Some(ci)) => {
                            pcurve_opt = Some(pc);
                            curveint = Some(ci);
                            maxapp1 = PROJ_TOL;
                        }
                        _ => raccordbatten = true,
                    }
                } else {
                    let slot = ic + 1;
                    let pc = self.proj2d1.get(slot).cloned().flatten();
                    let ci = self.cproj1.get(slot).cloned().flatten();
                    match (pc, ci) {
                        (Some(pc), Some(ci)) => {
                            pcurve_opt = Some(pc);
                            curveint = Some(ci);
                            maxapp1 = PROJ_TOL;
                        }
                        _ => raccordbatten = true,
                    }
                }
            }

            // `bool contraint1 = true, contraint2 = true;` (`cxx:3088`).
            let mut contraint1 = true;
            let mut contraint2 = true;

            // `if (raccordbatten) { ... CalculBatten(...); }` (`cxx:3089-3138`).
            if raccordbatten {
                // `inverseic` / `inverseicplus` (`cxx:3095-3110`).
                let inverseic = if self.sharp[ic] {
                    let edge = self.edge_at(ic)?;
                    let first = edge_vertices(&edge).0;
                    first.as_ref().is_some_and(|v| is_same(&v.0, &self.vtx))
                } else {
                    self.sens[ic] == 1
                };
                let inverseicplus = if self.sharp[icplus] {
                    let edge = self.edge_at(icplus)?;
                    let first = edge_vertices(&edge).0;
                    first.as_ref().is_some_and(|v| is_same(&v.0, &self.vtx))
                } else {
                    self.sens[icplus] == 1
                };
                // `cxx:3111-3119`: a free border never constrains the batten.
                let edge_ic = self.edge_at(ic)?;
                if self.is_free_border(&edge_ic) {
                    contraint1 = false;
                }
                let edge_icplus = self.edge_at(icplus)?;
                if self.is_free_border(&edge_icplus) {
                    contraint2 = false;
                }

                // `CalculBatten(Asurf, Face(Fvive(ic, icplus)), xdir, ydir,
                //  p2d1, p2d2, contraint1, contraint2, curv2d1, curv2d2,
                //  p(ic, icplus), p(icplus, ic), inverseic, inverseicplus,
                //  pcurve)` (`cxx:3120-3138`).
                pcurve_opt = Some(calcul_batten(
                    &face_surface,
                    &face,
                    xdir,
                    ydir,
                    &p2d1,
                    &p2d2,
                    contraint1,
                    contraint2,
                    &c1,
                    &curv2d2,
                    self.p[ic][icplus],
                    self.p[icplus][ic],
                    inverseic,
                    inverseicplus,
                )?);
            }

            let pcurve = pcurve_opt.ok_or_else(|| {
                "PerformMoreThreeCorner: the intermediary p-curve is null".to_string()
            })?;

            // `Acurv = new Geom2dAdaptor_Curve(pcurve);
            //  CurvOnS(Acurv, Asurf); HCons = new Adaptor3d_CurveOnSurface(
            //  CurvOnS);` (`cxx:3140-3143`).
            let cons = CurveOnSurface::new(Arc::clone(&pcurve), Arc::clone(&face_surface));

            // `Order.SetValue(n3d, 0);` then the five G1 promotions
            // (`cxx:3147-3168`).
            let mut order = 0;
            if !self.sharp[ic] && !self.sharp[icplus] {
                order = 1;
            }
            if !contraint1 && !self.sharp[icplus] {
                order = 1;
            }
            if !contraint2 && !self.sharp[ic] {
                order = 1;
            }
            if self.tangentregul[ic] || self.tangentregul[icplus] {
                order = 1;
            }
            if self.is_g1[ic] {
                order = 1;
            }
            self.order[*n3d] = order;

            // `Cont = new GeomPlate_CurveConstraint(HCons, Order(n3d), 10,
            //  tolapp3d, angular, 0.1); PSurf.Add(Cont);` (`cxx:3169-3172`).
            let cont = CurveConstraint::from_boundary(
                BoundaryCurve::OnSurface(Arc::new(cons.clone())),
                order,
                NBCURVPNT,
                tolapp3d,
                angular,
                0.1,
            )?;
            psurf.add_curve_constraint(cont);

            // `if (curveint.IsNull()) { GeomLib::BuildCurve3d(tolapp, CurvOnS,
            //  CurvOnS.FirstParameter(), CurvOnS.LastParameter(), Curv3d,
            //  maxapp1, avedev); pardeb = ...; parfin = ...;
            //  curveint = new Geom_TrimmedCurve(Curv3d, pardeb, parfin); }`
            // (`cxx:3174-3186`).
            if curveint.is_none() {
                let built = geom_lib_build_curve3d(
                    TOLAPP,
                    &cons,
                    &face_surface,
                    &pcurve,
                    cons.first_parameter(),
                    cons.last_parameter(),
                    BUILD_CONTINUITY,
                    BUILD_MAX_DEGREE,
                    BUILD_MAX_SEGMENT,
                );
                maxapp1 = built.max_deviation;
                let curv3d = built.curve.ok_or_else(|| {
                    "PerformMoreThreeCorner: GeomLib::BuildCurve3d is null".to_string()
                })?;
                let pardeb = cons.first_parameter();
                let parfin = cons.last_parameter();
                curveint = Some(Arc::new(GeomTrimmedCurve::new(curv3d, pardeb, parfin)));
            }
            let needle = curveint
                .ok_or_else(|| "PerformMoreThreeCorner: curveint is null".to_string())?;

            // `TopOpeBRepDS_Curve tcurv3d(curveint, maxapp1);
            //  indcurve3d.SetValue(n3d, DStr.AddCurve(tcurv3d));
            //  pardeb = curveint->FirstParameter();
            //  parfin = curveint->LastParameter();` (`cxx:3188-3192`).
            let indcurve = ds.add_curve(&mut DsCurve::with_curve(Some(Arc::clone(&needle)), maxapp1, false));
            self.indcurve3d[*n3d] = indcurve;
            let pardeb = needle.first_parameter();
            let parfin = needle.last_parameter();

            // `if (sharp(icplus) && indpoint(icplus, 0) == 0) { point2 =
            //  curveint->Value(parfin); tpoint2(point2, maxapp);
            //  indpoint(icplus, 0) = DStr.AddPoint(tpoint2);
            //  indpoint(icplus, 1) = indpoint(icplus, 0); }` (`cxx:3193-3201`).
            if self.sharp[icplus] && self.indpoint[icplus][0] == 0 {
                let point2 = needle.d0(parfin);
                // `TopOpeBRepDS_Point tpoint2(point2, maxapp)` uses the
                // `maxapp` of `cxx:1313` (`0.`), not `maxapp1`.
                let tpoint2 = DsPoint::with_tolerance(point2, 0.0);
                let idx = ds.add_point(&tpoint2);
                self.indpoint[icplus][0] = idx;
                self.indpoint[icplus][1] = idx;
            }

            // `bool IsVt1 = false, IsVt2 = false; if (deuxconges) {
            //  IsVt1 = sharp(ic); IsVt2 = sharp(icplus); }` (`cxx:3202-3207`).
            let (is_vt1, is_vt2) = if self.deuxconges {
                (self.sharp[ic], self.sharp[icplus])
            } else {
                (false, false)
            };

            // `Interfp1 = ChFi3d_FilPointInDS(FORWARD, indcurve3d(n3d),
            //  indpoint(ic, 1), pardeb, IsVt1); Interfp2 = ...REVERSED...,
            //  indpoint(icplus, 0), parfin, IsVt2);
            //  DStr.ChangeCurveInterferences(indcurve3d(n3d)).Append(...)` twice
            // (`cxx:3208-3221`).
            let interfp1 = ch_fi3d_fil_point_in_ds(
                Orientation::Forward,
                indcurve,
                self.indpoint[ic][1],
                pardeb,
                is_vt1,
            );
            let interfp2 = ch_fi3d_fil_point_in_ds(
                Orientation::Reversed,
                indcurve,
                self.indpoint[icplus][0],
                parfin,
                is_vt2,
            );
            ds.change_curve_interferences(indcurve).push(interfp1);
            ds.change_curve_interferences(indcurve).push(interfp2);

            // `if (!IsVt1) { tpt1 = DStr.ChangePoint(indpoint(ic, 1));
            //  tpt1.Tolerance(tpt1.Tolerance() + maxapp1); }` (`cxx:3222-3227`)
            // and its `icplus, 0` mirror (`cxx:3228-3232`).
            if !is_vt1 {
                let idx = self.indpoint[ic][1];
                if let Some(tpt1) = ds.change_point(idx) {
                    let tol = tpt1.tolerance();
                    tpt1.set_tolerance(tol + maxapp1);
                }
            }
            if !is_vt2 {
                let idx = self.indpoint[icplus][0];
                if let Some(tpt2) = ds.change_point(idx) {
                    let tol = tpt2.tolerance();
                    tpt2.set_tolerance(tol + maxapp1);
                }
            }

            // `cxx:3234-3253`: the orientation of the connecting curve.
            let orinterf = if !self.sharp[ic] {
                let stripe = self.stripe_at(ic)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                orientation_ic_non_vive(
                    &guard,
                    self.jf[ic],
                    self.i[ic][icplus],
                    self.sens[ic],
                )?
            } else if !self.sharp[icplus] {
                let stripe = self.stripe_at(icplus)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                orientation_icplus_non_vive(
                    &guard,
                    self.jf[icplus],
                    self.i[icplus][ic],
                    self.sens[icplus],
                )?
            } else {
                let face = self.fvive_face(ic, icplus)?.ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Fvive({ic}, {icplus}) is null")
                })?;
                let edge = self.edge_at(ic)?;
                orientation_arete_vive_consecutive(&face.0, &edge, &self.vtx)?
            };

            // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d),
            //  numfa(ic, icplus), pcurve, orinterf);
            //  DStr.ChangeShapeInterferences(numfa(ic, icplus)).Append(Interfc);`
            // (`cxx:3254-3257`).
            let numfa = self.numfa[ic][icplus];
            let interfc = ch_fi3d_fil_curve_in_ds(indcurve, numfa, Some(pcurve), orinterf);
            ds.change_shape_interferences_at(numfa).push(interfc);
        }

        Ok(())
    }
}
