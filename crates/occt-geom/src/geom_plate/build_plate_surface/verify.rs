//! `EcartContraintesMil` (`cxx:751-865`), `CourbeJointive`
//! (`cxx:1349-1444`), `VerifSurface` (`cxx:2588-2734`) and `VerifPoints`
//! (`cxx:2737-2782`) of `GeomPlate_BuildPlateSurface`.
//!
//! Translation notes:
//! - Every `NCollection_HArray1<double>(1, N)` becomes a `Vec<f64>` of length
//!   `N`; the OCCT `d->ChangeValue(i)` writes and `d->Value(i)` reads with
//!   `i` in `1..=N` become `d[i - 1]`.
//! - `EcartContraintesMil` writes its three arrays in place and `VerifSurface`
//!   sizes them `1..=NbPts_i-1`, so the port passes slices and drops a write that
//!   would fall outside them (OCCT would raise an out-of-range error there).
//! - `CG2` is default-constructed in OCCT (`cxx:839`, `:2772`); the port seeds
//!   it through [`SurfaceContinuity::new`] with zero thresholds, which
//!   `ComputeAnalysis` overwrites for the G2 study it runs.
//! - `myGeomPlateSurface` is only read here; a null handle is an error, where
//!   OCCT would dereference it.

use std::f64::consts::PI;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;

use crate::Surface;
use crate::geom_lprop::SlProps;
use crate::geom_plate::PlateSurface;
use crate::local_analysis::SurfaceContinuity;

use super::BuildPlateSurface;

impl BuildPlateSurface {
    /// `myGeomPlateSurface` with a null-handle error.
    fn plate_surface_ref(&self) -> Result<&PlateSurface, String> {
        self.geom_plate_surface
            .as_ref()
            .map(|s| s.as_ref())
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no surface".to_string())
    }

    /// The `ProjectedCurve` / `Curve2dOnSurf` / `ProjectPoint` cascade of
    /// `cxx:779-791` and `:807-819`.
    pub(super) fn constraint_parameter(&self, c: usize, u: f64, p3d: &GpPnt) -> Result<(f64, f64), String> {
        let cc = self.lin(c);
        if let Some(projected) = cc.projected_curve() {
            let p2d = projected.d0(u);
            return Ok((p2d.x(), p2d.y()));
        }
        if let Some(curve2d) = cc.curve2d_on_surf()? {
            let p2d = curve2d.d0(u);
            return Ok((p2d.x(), p2d.y()));
        }
        let p2d = self.project_point(p3d)?;
        Ok((p2d.x(), p2d.y()))
    }

    /// `EcartContraintesMil(c, d, an, courb)` (`cxx:751-865`): the gap, the
    /// angle and the curvature gap at the midpoints between consecutive
    /// constraint parameters of curve `c`.
    pub(super) fn ecart_contraintes_mil(
        &mut self,
        c: usize,
        d: &mut [f64],
        an: &mut [f64],
        courb: &mut [f64],
    ) -> Result<(), String> {
        // `cxx:758-767`.
        let nb_pt = if self.par_len(c) < 3 {
            4
        } else {
            self.par_len(c)
        };
        match self.lin(c).order() {
            0 => {
                // `cxx:775-798`.
                for i in 1..nb_pt {
                    let u = (self.par(c)[i - 1] + self.par(c)[i]) / 2.0;
                    let pi = self.lin(c).d0(u);
                    let (pu, pv) = self.constraint_parameter(c, u, &pi)?;
                    let pf = self.plate_surface_ref()?.d0(pu, pv);
                    an.fill(0.0);
                    courb.fill(0.0);
                    if let Some(slot) = d.get_mut(i - 1) {
                        *slot = pf.distance(&pi);
                    }
                }
            }
            1 => {
                // `cxx:801-835`.
                for i in 1..nb_pt {
                    let u = (self.par(c)[i - 1] + self.par(c)[i]) / 2.0;
                    let (pi, v1i, v2i) = self.lin(c).d1(u)?;
                    let (pu, pv) = self.constraint_parameter(c, u, &pi)?;
                    let (pf, v1f, v2f) = self.plate_surface_ref()?.d1(pu, pv);
                    if let Some(slot) = d.get_mut(i - 1) {
                        *slot = pf.distance(&pi);
                    }
                    let v3i = v1i.crossed(&v2i);
                    let v3f = v1f.crossed(&v2f);
                    let angle = v3f.angle(&v3i);
                    let folded = if angle > PI / 2.0 { PI - angle } else { angle };
                    if let Some(slot) = an.get_mut(i - 1) {
                        *slot = folded;
                    }
                    courb.fill(0.0);
                }
            }
            2 => {
                // `cxx:838-863`.
                for i in 1..nb_pt {
                    let u = (self.par(c)[i - 1] + self.par(c)[i]) / 2.0;
                    let pi = self.lin(c).d0(u);
                    let (pu, pv) = self.constraint_parameter(c, u, &pi)?;
                    let plate = self
                        .geom_plate_surface
                        .clone()
                        .ok_or_else(|| "GeomPlate_BuildPlateSurface: no surface".to_string())?;
                    // `GeomLProp_SLProps Prop(Splate, u, v, 2, 0.001)`.
                    let mut prop = SlProps::new(plate.as_ref(), pu, pv, 2, 0.001)?;
                    let mut cg2 =
                        SurfaceContinuity::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    let (lprop, surf) = self.lin_mut(c).l_prop_surf(u)?;
                    let surf = surf.clone();
                    cg2.compute_analysis(
                        &mut prop,
                        plate.as_ref(),
                        lprop,
                        surf.as_ref(),
                        Shape::G2,
                    )?;
                    if let Some(slot) = d.get_mut(i - 1) {
                        *slot = cg2.c0_value()?;
                    }
                    if let Some(slot) = an.get_mut(i - 1) {
                        *slot = cg2.g1_angle()?;
                    }
                    if let Some(slot) = courb.get_mut(i - 1) {
                        *slot = cg2.g2_curvature_gap()?;
                    }
                }
            }
            // `cxx:864`: the `default:` arm does nothing.
            _ => {}
        }
        Ok(())
    }

    /// `CourbeJointive(tolerance)` (`cxx:1349-1444`): chains the curves in
    /// order, swapping each neighbour into place and recording its orientation
    /// in `mySense`; `true` when the chain closes.
    pub(super) fn courbe_jointive(&mut self, tolerance: f64) -> Result<bool, String> {
        let nbf = self.lin_cont.len();
        self.sense = vec![0; nbf];
        let mut result = true;
        let mut j = 1i32;
        // `cxx:1357`.
        while j <= self.nb_bounds - 1 {
            let mut a = 0;
            let mut i = j + 1;
            // `cxx:1361`.
            if i > self.nb_bounds {
                result = false;
                a = 2;
            }
            // `cxx:1366`.
            while a < 1 {
                if i > self.nb_bounds {
                    result = false;
                    a = 2;
                } else {
                    // `cxx:1376-1379`.
                    let jc = (j - 1) as usize;
                    let ic = (i - 1) as usize;
                    let uinit1 = self.lin(jc).first_parameter();
                    let mut ufinal1 = self.lin(jc).last_parameter();
                    let uinit2 = self.lin(ic).first_parameter();
                    let ufinal2 = self.lin(ic).last_parameter();
                    // `cxx:1380-1383`.
                    if self.sense[jc] == 1 {
                        ufinal1 = uinit1;
                    }
                    let p1 = self.lin(jc).d0(ufinal1);
                    let p2 = self.lin(ic).d0(uinit2);
                    if p1.distance(&p2) < tolerance {
                        // `cxx:1388-1399`.
                        self.swap_neighbour(jc + 1, ic)?;
                        a = 2;
                        self.sense[jc + 1] = 0;
                    } else {
                        let p2b = self.lin(ic).d0(ufinal2);
                        if p1.distance(&p2b) < tolerance {
                            // `cxx:1405-1418`.
                            self.swap_neighbour(jc + 1, ic)?;
                            a = 2;
                            self.sense[jc + 1] = 1;
                        }
                    }
                }
                i += 1;
            }
            j += 1;
        }
        // `cxx:1425-1429`.
        let nb = self.nb_bounds;
        let uinit1 = self.lin((nb - 1) as usize).first_parameter();
        let ufinal1 = self.lin((nb - 1) as usize).last_parameter();
        let uinit2 = self.lin(0).first_parameter();
        let p1 = self.lin((nb - 1) as usize).d0(ufinal1);
        let p2 = self.lin(0).d0(uinit2);
        if self.sense[(nb - 1) as usize] == 0 && p1.distance(&p2) < tolerance {
            return Ok(result);
        }
        let p1 = self.lin((nb - 1) as usize).d0(uinit1);
        if self.sense[(nb - 1) as usize] == 1 && p1.distance(&p2) < tolerance {
            Ok(result)
        } else {
            Ok(false)
        }
    }

    /// The `if (i != j + 1)` swap of `cxx:1389-1397` / `:1408-1416`: exchanges
    /// the two curve constraints and their `myInitOrder`
    /// entries.
    fn swap_neighbour(&mut self, a: usize, b: usize) -> Result<(), String> {
        if a == b {
            return Ok(());
        }
        self.lin_cont.swap(a, b);
        let a_order = self.init_order.get(a).copied();
        let b_order = self.init_order.get(b).copied();
        match (a_order, b_order) {
            (Some(ao), Some(bo)) => {
                self.init_order[a] = bo;
                self.init_order[b] = ao;
                Ok(())
            }
            _ => Err("GeomPlate_BuildPlateSurface: no initial order".to_string()),
        }
    }

    /// `VerifSurface(NbBoucle)` (`cxx:2588-2734`): accumulates the three
    /// maximum errors and refines `NbPoints` on the constraints that exceed
    /// their criteria; `false` asks `Perform` for another loop.
    pub(super) fn verif_surface(&mut self, nb_boucle: i32) -> Result<bool, String> {
        let nt_lin_cont = self.lin_cont.len();
        let mut result = true;
        // `cxx:2597`.
        self.g0_error = 0.0;
        self.g1_error = 0.0;
        self.g2_error = 0.0;
        for i in 0..nt_lin_cont {
            if self.lin(i).order() == -1 {
                continue;
            }
            // `cxx:2605-2617`.
            let mut nb_pts_i = self.par_len(i);
            if nb_pts_i < 3 {
                nb_pts_i = 4;
            }
            let mut tdist = vec![0.0; nb_pts_i - 1];
            let mut tang = vec![0.0; nb_pts_i - 1];
            let mut tcourb = vec![0.0; nb_pts_i - 1];
            self.ecart_contraintes_mil(i, &mut tdist, &mut tang, &mut tcourb)?;
            let mut diff_dist_max = 0.0;
            let mut diff_ang_max = 0.0;
            let mut n_diff_dist = 0;
            let mut n_diff_ang = 0;
            for j in 1..nb_pts_i {
                // `cxx:2627-2638`.
                if tdist[j - 1] > self.g0_error {
                    self.g0_error = tdist[j - 1];
                }
                if tang[j - 1] > self.g1_error {
                    self.g1_error = tang[j - 1];
                }
                if tcourb[j - 1] > self.g2_error {
                    self.g2_error = tcourb[j - 1];
                }
                // `cxx:2640-2646`.
                let u = if self.par_len(i) > 3 {
                    (self.par(i)[j - 1] + self.par(i)[j]) / 2.0
                } else {
                    let first = self.lin(i).first_parameter();
                    let last = self.lin(i).last_parameter();
                    first + (last - first) * ((j as f64) - 1.0) / ((nb_pts_i as f64) - 2.0)
                };
                // `cxx:2645-2655`.
                let mut diff_dist = tdist[j - 1] - self.lin(i).g0_criterion(u)?;
                let diff_ang = if self.lin(i).order() > 0 {
                    tang[j - 1] - self.lin(i).g1_criterion(u)?
                } else {
                    0.0
                };
                // `cxx:2652-2676`.
                if diff_dist > 0.0 {
                    diff_dist /= self.lin(i).g0_criterion(u)?;
                    if diff_dist > diff_dist_max {
                        diff_dist_max = diff_dist;
                    }
                    n_diff_dist += 1;
                } else if diff_ang > 0.0 && self.lin(i).order() == 1 {
                    let diff_ang = diff_ang / self.lin(i).g1_criterion(u)?;
                    if diff_ang > diff_ang_max {
                        diff_ang_max = diff_ang;
                    }
                    n_diff_ang += 1;
                }
            }
            // `cxx:2679-2708`.
            if n_diff_dist > 0 {
                let mut coef = if self.lin(i).order() == 0 {
                    0.6 * (diff_dist_max + 7.4).ln()
                } else {
                    (diff_dist_max + 3.3).ln()
                };
                if coef > 3.0 {
                    coef = 3.0;
                }
                if nb_boucle > 1 && diff_dist_max > 2.0 {
                    coef = 1.6;
                }
                let nb_points = self.lin(i).nb_points();
                if (nb_points as f64) >= ((nb_points as f64) * coef).floor() {
                    coef = 2.0;
                }
                self.lin_mut(i)
                    .set_nb_points((nb_points as f64 * coef) as i32);
                result = false;
            } else if n_diff_ang > 0 {
                // `cxx:2710-2719`.
                let mut coef = 1.5;
                let nb_points = self.lin(i).nb_points();
                if (nb_points + 1) as f64 >= ((nb_points as f64) * coef).floor() {
                    coef = 2.0;
                }
                self.lin_mut(i)
                    .set_nb_points((nb_points as f64 * coef) as i32);
                result = false;
            }
        }
        // `cxx:2723-2729`.
        if !result {
            if self.free && nb_boucle == 1 {
                self.prev_plate = self.plate.clone();
            }
            self.plate.init();
        }
        Ok(result)
    }

    /// `VerifPoints(Dist, Ang, Curv)` (`cxx:2737-2782`): the errors of the last
    /// point constraint (OCCT overwrites the three outputs on every iteration).
    pub(super) fn verif_points(&mut self) -> Result<(f64, f64, f64), String> {
        let mut ang = 0.0;
        let mut dist = 0.0;
        let mut curv = 0.0;
        let nt_pnt_cont = self.pnt_cont.len();
        for i in 0..nt_pnt_cont {
            match self.pnt(i).order() {
                0 => {
                    // `cxx:2752-2756`.
                    let p2d = self.pnt(i).pnt2d_on_surf();
                    let pi = self.pnt(i).d0();
                    let pf = self.plate_surface_ref()?.d0(p2d.x(), p2d.y());
                    dist = pf.distance(&pi);
                }
                1 => {
                    // `cxx:2757-2769`.
                    let (pi, v1i, v2i) = self.pnt(i).d1();
                    let p2d = self.pnt(i).pnt2d_on_surf();
                    let (pf, v1f, v2f) = self.plate_surface_ref()?.d1(p2d.x(), p2d.y());
                    dist = pf.distance(&pi);
                    let v3i: GpVec = v1i.crossed(&v2i);
                    let v3f: GpVec = v1f.crossed(&v2f);
                    let angle = v3f.angle(&v3i);
                    ang = if angle > PI / 2.0 { PI - angle } else { angle };
                }
                2 => {
                    // `cxx:2770-2779`.
                    let p2d = self.pnt(i).pnt2d_on_surf();
                    let plate = self
                        .geom_plate_surface
                        .clone()
                        .ok_or_else(|| "GeomPlate_BuildPlateSurface: no surface".to_string())?;
                    let mut prop = SlProps::new(plate.as_ref(), p2d.x(), p2d.y(), 2, 0.001)?;
                    let mut cg2 =
                        SurfaceContinuity::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    let (lprop, surf) = self.pnt_mut(i).l_prop_surf();
                    let surf = surf
                        .cloned()
                        .ok_or_else(|| "GeomPlate_PointConstraint: no surface".to_string())?;
                    cg2.compute_analysis(
                        &mut prop,
                        plate.as_ref(),
                        lprop,
                        surf.as_ref(),
                        Shape::G2,
                    )?;
                    dist = cg2.c0_value()?;
                    ang = cg2.g1_angle()?;
                    curv = cg2.g2_curvature_gap()?;
                }
                // `cxx:2781`: the `default:` arm does nothing.
                _ => {}
            }
        }
        Ok((dist, ang, curv))
    }
}

