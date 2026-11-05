//! The tails of the `GeomPlate` section: p-curve orientations, the connection
//! curves, the regular-strip storage, and the partial-result branch.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx`
//! (OCCT 8.0.0):
//!
//! - `cxx:3628-3674`: `trafil1` / `trafil2` / `orpcurve` / `orien`, the
//!   orientation every plate p-curve interference is stored with,
//! - `cxx:3676-3730`: the free-border connection curves (`!droit`), one
//!   `ChFi3d_ComputeArete` per `libre(ic)` slot,
//! - `cxx:3732-3776`: the stripe curves against `Isurf`, with the `errapp` /
//!   `apperror` tolerance accumulation and the regular-strip records,
//! - `cxx:3778-3836`: the living edges (`sharp`) and the connection curves
//!   (`!ponctuel && !libre && !moresurf`),
//! - `cxx:3838-3872`: the multi-face projected curves against `Isurf`,
//! - `cxx:3874-3890`: the tangent free borders (`droit`),
//! - `cxx:3891-3926`: the `else` of `PSurf.IsDone()`, where `done` is cleared
//!   and `hasresult` kept so that `Compute` builds a partial result.
//!
//! Translation notes:
//! - the `!droit` loop runs `ic <= nedge` (`cxx:3679`) while `libre` is only
//!   written for `ic < nedge` (`cxx:1793`), so the extra iteration never fires;
//!   the port loops `0..nedge` and keeps the bound note here,
//! - `PSurf.Curves2d()->Value(n3d)` is `pc_arr[n3d - 1]`;
//!   [`BuildPlateSurface::curves2d`] answers the post-`Perform` ordering,
//! - `DStr.ChangeCurve(i).Tolerance(t + apperror)` and the matching
//!   `ChangePoint(i).Tolerance(...)` are read-modify-writes on the DS tables,
//! - `regular` is the function-scope `ChFiDS_Regul` of `cxx:1316`: the port
//!   keeps one `Regul` across the loops and appends a copy at each site, so a
//!   field left over from an earlier block would carry exactly as in OCCT,
//! - `ChFi3d_ComputeArete(..., IFlag = 0)` is
//!   [`ch_fi3d_compute_arete`] with an empty incoming [`ComputeArete`].

use std::sync::Arc;

use occt_geom::geom_plate::BuildPlateSurface;

use crate::abs::Orientation;
use crate::chfi_ds::{CommonPoint, Regul};
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ComputeArete,
};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::plate_helpers::parametre_plate;
use super::surf_index::indices;
use super::MoreThreeCorner;

/// `if (Order.Value(n3d) == 1)` (`cxx:3825`, `:3861`) — the regular strips are
/// only recorded for the G1 constraints.
const ORDER_G1: i32 = 1;

impl MoreThreeCorner {
    /// `cxx:3628-3926`. Answers `true` for the `PSurf.IsDone()` branch,
    /// `false` for the partial-result one.
    pub(super) fn plate_tails(
        &mut self,
        ds: &mut DataStructure,
        regul: &mut Vec<Regul>,
        psurf: &BuildPlateSurface,
    ) -> Result<bool, String> {
        let nedge = self.nedge as usize;
        let isurf = self.isurf;
        let apperror = self.apperror;
        let surf = self
            .plate_surface
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: the plate surface is missing".to_string())?;
        let pc_arr = psurf.curves2d()?;

        // --- `cxx:3628-3674`: the p-curve orientation of every plate curve.
        // `Ishape1 = Fd->IndexOfS1(); Ishape2 = Fd->IndexOfS2();
        //  Fi1 = Fd->InterferenceOnS1(); Fi2 = ...OnS2();
        //  trafil1 = trafil2 = FORWARD;` (`cxx:3632-3637`).
        let fd = {
            let (_icplus0, icmoins0) = indices(nedge, 0);
            let stripe = self.stripe_at(0)?.clone();
            let guard = stripe
                .read()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            let idx = self.i[0][icmoins0];
            let pos = usize::try_from(idx)
                .ok()
                .and_then(|i| i.checked_sub(1))
                .ok_or_else(|| {
                    format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
                })?;
            Arc::clone(guard.set_of_surf_data().get(pos).ok_or_else(|| {
                format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
            })?)
        };
        let ishape1 = fd.index_of_s1();
        let ishape2 = fd.index_of_s2();
        let fi1 = fd.interference_on_s1().clone();
        let fi2 = fd.interference_on_s2().clone();
        let mut trafil1 = Orientation::Forward;
        let mut trafil2 = Orientation::Forward;

        // `if (Ishape1 != 0) { if (Ishape1 > 0) trafil1 =
        //  DStr.Shape(Ishape1).Orientation(); trafil1 = Compose(trafil1,
        //  Fd->Orientation()); trafil1 = Compose(Reverse(Fi1.Transition()),
        //  trafil1); trafil2 = Reverse(trafil1); }` (`cxx:3638-3648`).
        if ishape1 != 0 {
            if ishape1 > 0 {
                trafil1 = self.shape_orientation_of(ds, ishape1)?;
            }
            trafil1 = Orientation::compose(trafil1, fd.orientation());
            trafil1 = Orientation::compose(fi1.transition().reversed(), trafil1);
            trafil2 = trafil1.reversed();
        } else {
            // `cxx:3649-3656`: the `Ishape2` mirror. The `Ishape2 > 0` arm is
            // unconditional in OCCT, so a `trafil2` of `FORWARD` reaches the
            // composition when `Ishape2 <= 0`.
            if ishape2 > 0 {
                trafil2 = self.shape_orientation_of(ds, ishape2)?;
            }
            trafil2 = Orientation::compose(trafil2, fd.orientation());
            trafil2 = Orientation::compose(fi2.transition().reversed(), trafil2);
            trafil1 = trafil2.reversed();
        }
        // OCCT keeps computing `trafil2` for both arms (`cxx:3647`, `:3655`)
        // and then only reads `trafil1`; the port keeps the writes and marks
        // the value unused rather than dropping a branch of the source.
        let _ = trafil2;

        // `if (isfirst) { orpcurve = Reverse(trafil1); orpcurve =
        //  Compose(orpcurve, CD.Value(0)->FirstPCurveOrientation()); } else {
        //  orpcurve = trafil1; orpcurve = Compose(orpcurve,
        //  CD.Value(0)->LastPCurveOrientation()); }` (`cxx:3658-3665`).
        let first_p_curve_orientation = {
            let stripe = self.stripe_at(0)?.clone();
            let guard = stripe
                .read()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            if self.is_first_plate {
                (true, guard.first_p_curve_orientation())
            } else {
                (false, guard.last_p_curve_orientation())
            }
        };
        let mut orpcurve = if first_p_curve_orientation.0 {
            trafil1.reversed()
        } else {
            trafil1
        };
        orpcurve = Orientation::compose(orpcurve, first_p_curve_orientation.1);

        // `if (orsurfdata == orplate) orien = Reverse(orpcurve);
        //  else orien = orpcurve;` (`cxx:3667-3673`).
        let orien = if self.or_surfdata == self.or_plate {
            orpcurve.reversed()
        } else {
            orpcurve
        };

        // --- `cxx:3676-3730`: the connection curves of the free borders.
        if !self.droit {
            for ic in 0..nedge {
                if !self.libre[ic] {
                    continue;
                }
                let (icplus, _icmoins) = indices(nedge, ic);

                // `par1 = p(ic, icplus); par2 = p(icplus, ic);
                //  BCurv1.D0(par1, ptic); BCurv2.D0(par2, pticplus);`
                // (`cxx:3693-3695`).
                let edge1 = self.edge_at(ic)?;
                let edge2 = self.edge_at(icplus)?;
                let par1 = self.p[ic][icplus];
                let par2 = self.p[icplus][ic];
                let ptic = CurveAdaptor::of_edge(&edge1)
                    .d0(par1)
                    .ok_or_else(|| "PerformMoreThreeCorner: Evive has no 3d curve".to_string())?;
                let pticplus = CurveAdaptor::of_edge(&edge2)
                    .d0(par2)
                    .ok_or_else(|| "PerformMoreThreeCorner: Evive has no 3d curve".to_string())?;

                // `ParametrePlate(n3d, PSurf, Surf, ptic, apperror, UV1);`
                // (`cxx:3696-3697`).
                let uv1 = parametre_plate(self.n3d_curves, psurf, &surf, &ptic, apperror)?;
                let uv2 = parametre_plate(self.n3d_curves, psurf, &surf, &pticplus, apperror)?;

                // `CP1.SetArc(1.e-3, Edge(Evive(ic)), par1, FORWARD);
                //  CP1.SetPoint(ptic);` and the `CP2` mirror
                // (`cxx:3698-3703`).
                let mut cp1 = CommonPoint::new();
                cp1.set_arc(1.0e-3, &edge1.0, par1, Orientation::Forward);
                cp1.set_point(&ptic);
                let mut cp2 = CommonPoint::new();
                cp2.set_arc(1.0e-3, &edge2.0, par2, Orientation::Forward);
                cp2.set_point(&pticplus);

                // `ChFi3d_ComputeArete(CP1, UV1, CP2, UV2, Surf, C3d, C2d,
                //  param1, param2, to3d = 1.e-3, to2d = 1.e-6, tolreached, 0);`
                // (`cxx:3704-3717`).
                let mut inout = ComputeArete {
                    c3d: None,
                    p_curve: None,
                    pardeb: 0.0,
                    parfin: 0.0,
                    tolreached: 0.0,
                };
                ch_fi3d_compute_arete(&cp1, &uv1, &cp2, &uv2, &surf, &mut inout, 1.0e-3, 1.0e-6, 0)?;
                let c3d = inout
                    .c3d
                    .clone()
                    .ok_or_else(|| "PerformMoreThreeCorner: ComputeArete C3d is null".to_string())?;
                let c2d = inout.p_curve.clone();

                // `TopOpeBRepDS_Curve tcurv3d(C3d, tolreached);
                //  ind1 = indpoint(ic, 0); ind2 = indpoint(icplus, 0);
                //  indcurv = DStr.AddCurve(tcurv3d);` (`cxx:3718-3722`).
                let ind1 = self.indpoint[ic][0];
                let ind2 = self.indpoint[icplus][0];
                let indcurv = ds.add_curve(&mut crate::topopebrep::ds::Curve::with_curve(
                    Some(c3d),
                    inout.tolreached,
                    false,
                ));

                // `Interfp1 = ChFi3d_FilPointInDS(FORWARD, indcurv, ind1,
                //  param1); Interfp2 = ...REVERSED..., ind2, param2);`
                // (`cxx:3723-3726`).
                let interfp1 = ch_fi3d_fil_point_in_ds(
                    Orientation::Forward,
                    indcurv,
                    ind1,
                    inout.pardeb,
                    false,
                );
                let interfp2 = ch_fi3d_fil_point_in_ds(
                    Orientation::Reversed,
                    indcurv,
                    ind2,
                    inout.parfin,
                    false,
                );
                ds.change_curve_interferences(indcurv).push(interfp1);
                ds.change_curve_interferences(indcurv).push(interfp2);

                // `Interfc = ChFi3d_FilCurveInDS(indcurv, Isurf, C2d, orien);
                //  DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc);`
                // (`cxx:3727-3728`).
                let interfc = ch_fi3d_fil_curve_in_ds(indcurv, isurf, c2d, orien);
                ds.change_surface_interferences(isurf).push(interfc);
            }
        }

        // --- `cxx:3732-3776`: the stripe curves.
        // `n3d = 0;` (`cxx:3733`).
        self.n3d_curves = 0;
        let mut regular = Regul::new();
        for ic in 0..nedge {
            if self.sharp[ic] {
                continue;
            }
            self.n3d_curves += 1;
            let n3d = self.n3d_curves;
            let (icplus, icmoins) = indices(nedge, ic);
            let _ = icplus;

            // `isfirst = (sens(ic) == 1); apperror = Mapp.CriterionError() *
            //  coef; pardeb = CD.Value(ic)->PCurve(isfirst)->FirstParameter();
            //  parfin = ...->LastParameter();` (`cxx:3742-3746`).
            let isfirst = self.sens[ic] == 1;
            let p_curve = {
                let stripe = self.stripe_at(ic)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                guard
                    .p_curve(isfirst)
                    .cloned()
                    .ok_or_else(|| "PerformMoreThreeCorner: the stripe has no p-curve".to_string())?
            };
            let pardeb = p_curve.first_parameter();
            let parfin = p_curve.last_parameter();
            let indcurve = self.indcurve3d[n3d];

            // `Interfp1 = ChFi3d_FilPointInDS(FORWARD, indcurve3d(n3d),
            //  indpoint(ic, 0), pardeb); Interfp2 = ...REVERSED...,
            //  indpoint(ic, 1), parfin);` (`cxx:3748-3755`).
            let interfp1 = ch_fi3d_fil_point_in_ds(
                Orientation::Forward,
                indcurve,
                self.indpoint[ic][0],
                pardeb,
                false,
            );
            let interfp2 = ch_fi3d_fil_point_in_ds(
                Orientation::Reversed,
                indcurve,
                self.indpoint[ic][1],
                parfin,
                false,
            );
            ds.change_curve_interferences(indcurve).push(interfp1);
            ds.change_curve_interferences(indcurve).push(interfp2);

            // `tcourb.Tolerance(errapp.Value(ic) + apperror);`
            // (`cxx:3758`).
            self.set_curve_tolerance(ds, indcurve, self.errapp[ic] + apperror);
            // `tpt1.Tolerance(tpt1.Tolerance() + apperror);` for both ends
            // (`cxx:3759-3762`) — a missing DS point is the out-of-range raise
            // of `ChangePoint`.
            self.bump_point_tolerance(ds, self.indpoint[ic][0], apperror)?;
            self.bump_point_tolerance(ds, self.indpoint[ic][1], apperror)?;

            // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d), Isurf,
            //  PSurf.Curves2d()->Value(n3d), orien);
            //  DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc);
            //  regular.SetCurve(indcurve3d(n3d)); regular.SetS1(Isurf, false);
            //  indice = CD.Value(ic)->SetOfSurfData()->Value(i(ic, icmoins))
            //  ->Surf(); regular.SetS2(indice, false); myRegul.Append(regular);`
            // (`cxx:3766-3773`).
            let pc = pc_at(&pc_arr, n3d)?;
            let interfc = ch_fi3d_fil_curve_in_ds(indcurve, isurf, Some(pc), orien);
            ds.change_surface_interferences(isurf).push(interfc);
            regular.set_curve(indcurve);
            regular.set_s1(isurf, false);
            let indice = {
                let stripe = self.stripe_at(ic)?.clone();
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                let idx = self.i[ic][icmoins];
                let pos = usize::try_from(idx)
                    .ok()
                    .and_then(|i| i.checked_sub(1))
                    .ok_or_else(|| {
                        format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
                    })?;
                guard
                    .set_of_surf_data()
                    .get(pos)
                    .ok_or_else(|| {
                        format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
                    })?
                    .surf()
            };
            regular.set_s2(indice, false);
            regul.push(regular);
        }

        // --- `cxx:3778-3836`: living edges and connection curves.
        for ic in 0..nedge {
            let (icplus, _icmoins) = indices(nedge, ic);
            if self.oksea[ic] {
                continue;
            }

            // `if (sharp(ic) && !deuxconges)`: the living edge limitation
            // (`cxx:3782-3803`).
            if self.sharp[ic] && !self.deuxconges {
                let edge = self.edge_at(ic)?;
                let (vd, vf) = crate::topo_tools::edge_vertices(&edge);
                let vd = vd.ok_or_else(|| {
                    "PerformMoreThreeCorner: Evive has no first vertex".to_string()
                })?;
                let vf = vf.ok_or_else(|| {
                    "PerformMoreThreeCorner: Evive has no last vertex".to_string()
                })?;
                let pf = crate::brep_tool::BRepTool::vertex_point(&vd);
                let pl = crate::brep_tool::BRepTool::vertex_point(&vf);
                let sommet1 = {
                    let v1 = crate::shape::Vertex::wrap(self.vtx.clone()).ok_or_else(|| {
                        "PerformMoreThreeCorner: V1 is not a vertex".to_string()
                    })?;
                    crate::brep_tool::BRepTool::vertex_point(&v1)
                };
                let ori = if pf.distance(&sommet1) < pl.distance(&sommet1) {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                // `int Iarc1 = DStr.AddShape(Edge(Evive(ic)));
                //  Interfp1 = ChFi3d_FilPointInDS(ori, Iarc1,
                //  indpoint(ic, 1), p(ic, icplus));
                //  DStr.ChangeShapeInterferences(Edge).Append(Interfp1);`
                // (`cxx:3801-3803`).
                let iarc1 = ds.add_shape(&edge.0);
                let interfp1 = ch_fi3d_fil_point_in_ds(
                    ori,
                    iarc1,
                    self.indpoint[ic][1],
                    self.p[ic][icplus],
                    false,
                );
                ds.change_shape_interferences_at(iarc1).push(interfp1);
            }

            // `if (!ponctuel(ic) && !libre(ic))` (`cxx:3805-3835`).
            if self.ponctuel[ic] || self.libre[ic] {
                continue;
            }
            if self.moresurf[ic] {
                continue;
            }
            self.n3d_curves += 1;
            let n3d = self.n3d_curves;
            let indcurve = self.indcurve3d[n3d];

            // `tcourb1.Tolerance(tcourb1.Tolerance() + apperror);`
            // (`cxx:3813`).
            self.bump_curve_tolerance(ds, indcurve, apperror);
            // `if (!deuxconges) { tpt11 = ChangePoint(indpoint(ic, 1));
            //  tpt21 = ChangePoint(indpoint(icplus, 0)); ... }`
            // (`cxx:3814-3821`).
            if !self.deuxconges {
                self.bump_point_tolerance(ds, self.indpoint[ic][1], apperror)?;
                self.bump_point_tolerance(ds, self.indpoint[icplus][0], apperror)?;
            }

            // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d), Isurf,
            //  PSurf.Curves2d()->Value(n3d), orien);
            //  DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc);`
            // (`cxx:3822-3824`).
            let pc = pc_at(&pc_arr, n3d)?;
            let interfc = ch_fi3d_fil_curve_in_ds(indcurve, isurf, Some(pc), orien);
            ds.change_surface_interferences(isurf).push(interfc);

            // `if (Order.Value(n3d) == 1) { regular.SetCurve(...);
            //  regular.SetS1(Isurf, false); regular.SetS2(numfa.Value(ic,
            //  icplus)); myRegul.Append(regular); }` (`cxx:3825-3833`).
            if self.order[n3d] == ORDER_G1 {
                regular.set_curve(indcurve);
                regular.set_s1(isurf, false);
                regular.set_s2(self.numfa[ic][icplus], true);
                regul.push(regular);
            }
        }

        // --- `cxx:3838-3872`: the multi-face projected curves.
        let nbface = self.nbface as usize;
        for ic in 0..nedge {
            if !self.moresurf[ic] {
                continue;
            }
            for nb in 1..=nbface {
                if self.eproj.get(nb - 1).cloned().flatten().is_none() {
                    continue;
                }
                self.n3d_curves += 1;
                let n3d = self.n3d_curves;
                let indcurve = self.indcurve3d[n3d];

                // `tcourb1.Tolerance(tcourb1.Tolerance() + apperror);`
                // (`cxx:3850`).
                self.bump_curve_tolerance(ds, indcurve, apperror);
                // `if (Indice.Value(n3d) != 0) tpt11 = ChangePoint(
                //  Indice.Value(n3d)); tpt11.Tolerance(... + apperror);`
                // (`cxx:3851-3855`).
                if self.indice[n3d] != 0 {
                    self.bump_point_tolerance(ds, self.indice[n3d], apperror)?;
                }

                // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d), Isurf,
                //  PSurf.Curves2d()->Value(n3d), orien);
                //  DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc);`
                // (`cxx:3857-3860`).
                let pc = pc_at(&pc_arr, n3d)?;
                let interfc = ch_fi3d_fil_curve_in_ds(indcurve, isurf, Some(pc), orien);
                ds.change_surface_interferences(isurf).push(interfc);

                // `if (Order.Value(n3d) == 1) { ... regular.SetS2(
                //  DStr.AddShape(Face(Fproj.Value(nb)))); myRegul.Append(
                //  regular); }` (`cxx:3861-3868`).
                if self.order[n3d] == ORDER_G1 {
                    let fproj_shape = self.fproj.get(nb - 1).cloned().ok_or_else(|| {
                        format!("PerformMoreThreeCorner: Fproj({nb}) is missing")
                    })?;
                    let iface = ds.add_shape(&fproj_shape);
                    regular.set_curve(indcurve);
                    regular.set_s1(isurf, false);
                    regular.set_s2(iface, true);
                    regul.push(regular);
                }
            }
        }

        // --- `cxx:3874-3890`: the tangent free borders.
        if self.droit {
            for ic in 0..nedge {
                let ecur = self.edge_at(ic)?;
                if !self.is_free_border(&ecur) {
                    continue;
                }
                self.n3d_curves += 1;
                let n3d = self.n3d_curves;
                let indcurve = self.indcurve3d[n3d];
                // `tcourb1.Tolerance(tcourb1.Tolerance() + apperror);`
                // (`cxx:3885`).
                self.bump_curve_tolerance(ds, indcurve, apperror);
                // `Interfc = ChFi3d_FilCurveInDS(indcurve3d(n3d), Isurf,
                //  PSurf.Curves2d()->Value(n3d), orien);`
                // (`cxx:3886-3887`).
                let pc = pc_at(&pc_arr, n3d)?;
                let interfc = ch_fi3d_fil_curve_in_ds(indcurve, isurf, Some(pc), orien);
                ds.change_surface_interferences(isurf).push(interfc);
            }
        }

        // `else { done = false; hasresult = true; ... }` (`cxx:3891-3926`) is
        // answered `false` here; the caller writes the two flags.
        Ok(true)
    }

    /// `cxx:3891-3926`: the partial-result branch of the `else` of
    /// `PSurf.IsDone()`.
    pub(super) fn plate_partial_tail(&mut self, ds: &mut DataStructure) -> Result<(), String> {
        let nedge = self.nedge as usize;
        for ic in 0..nedge {
            let (icplus, _icmoins) = indices(nedge, ic);
            if self.oksea[ic] {
                continue;
            }
            if !(self.sharp[ic] && !self.deuxconges) {
                continue;
            }
            // `Vd = FirstVertex(Evive(ic)); Vf = LastVertex(...);
            //  sommet1 = Pnt(V1); if (Pf.Distance(sommet1) <
            //  Pl.Distance(sommet1)) ori = FORWARD; else ori = REVERSED;`
            // (`cxx:3903-3919`).
            let edge = self.edge_at(ic)?;
            let (vd, vf) = crate::topo_tools::edge_vertices(&edge);
            let vd = vd
                .ok_or_else(|| "PerformMoreThreeCorner: Evive has no first vertex".to_string())?;
            let vf = vf
                .ok_or_else(|| "PerformMoreThreeCorner: Evive has no last vertex".to_string())?;
            let pf = crate::brep_tool::BRepTool::vertex_point(&vd);
            let pl = crate::brep_tool::BRepTool::vertex_point(&vf);
            let sommet1 = {
                let v1 = crate::shape::Vertex::wrap(self.vtx.clone())
                    .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;
                crate::brep_tool::BRepTool::vertex_point(&v1)
            };
            let ori = if pf.distance(&sommet1) < pl.distance(&sommet1) {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            // `int Iarc1 = DStr.AddShape(Edge(Evive(ic)));
            //  Interfp1 = ChFi3d_FilPointInDS(ori, Iarc1, indpoint(ic, 1),
            //  p(ic, icplus)); DStr.ChangeShapeInterferences(Edge).Append(...);`
            // (`cxx:3920-3923`).
            let iarc1 = ds.add_shape(&edge.0);
            let interfp1 = ch_fi3d_fil_point_in_ds(
                ori,
                iarc1,
                self.indpoint[ic][1],
                self.p[ic][icplus],
                false,
            );
            ds.change_shape_interferences_at(iarc1).push(interfp1);
        }
        Ok(())
    }

    /// `DStr.Shape(i).Orientation()` (`cxx:3640`, `:3649`).
    fn shape_orientation_of(&self, ds: &DataStructure, i: i32) -> Result<Orientation, String> {
        ds.shape(i, false)
            .map(|s| s.orientation())
            .ok_or_else(|| format!("PerformMoreThreeCorner: DS shape {i} is not bound"))
    }

    /// `tcourb.Tolerance(errapp(ic) + apperror)` (`cxx:3758`).
    fn set_curve_tolerance(&self, ds: &mut DataStructure, i: i32, t: f64) {
        ds.change_curve(i).set_tolerance(t);
    }

    /// `DStr.ChangeCurve(i).Tolerance(t + apperror)` (`cxx:3813`).
    fn bump_curve_tolerance(&self, ds: &mut DataStructure, i: i32, apperror: f64) {
        let tol = ds.curve(i).tolerance();
        ds.change_curve(i).set_tolerance(tol + apperror);
    }

    /// `DStr.ChangePoint(i).Tolerance(t + apperror)` (`cxx:3760-3762` etc).
    fn bump_point_tolerance(
        &self,
        ds: &mut DataStructure,
        i: i32,
        apperror: f64,
    ) -> Result<(), String> {
        let tpt = ds
            .change_point(i)
            .ok_or_else(|| format!("PerformMoreThreeCorner: DS point {i} is out of range"))?;
        tpt.set_tolerance(tpt.tolerance() + apperror);
        Ok(())
    }
}

/// `PSurf.Curves2d()->Value(n3d)` — 1-based, and a null entry is the
/// `Standard_NullObject` the OCCT call would raise on.
fn pc_at(
    pc_arr: &[Option<Arc<dyn occt_geom2d::Curve2d>>],
    n3d: usize,
) -> Result<Arc<dyn occt_geom2d::Curve2d>, String> {
    pc_arr
        .get(n3d.checked_sub(1).ok_or_else(|| {
            "PerformMoreThreeCorner: PSurf.Curves2d()->Value(0) is out of range".to_string()
        })?)
        .and_then(|c| c.as_ref())
        .cloned()
        .ok_or_else(|| format!("PerformMoreThreeCorner: PSurf.Curves2d()->Value({n3d}) is null"))
}
