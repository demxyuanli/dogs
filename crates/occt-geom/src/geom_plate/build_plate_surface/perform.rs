//! `Perform` (`GeomPlate_BuildPlateSurface.cxx:438-747`), `ComputeSurfInit`
//! (`cxx:1454-1902`) and `CalculNbPtsInit` (`cxx:2366-2404`).
//!
//! Translation notes:
//! - `Perform` is `void` and returns early on every failure arm. The port keeps
//!   those arms as [`PerformOutcome`] values so the caller can tell "no
//!   constraints" / "no initial surface" / "plate failed" from a finished run,
//!   and reserves `Err` for the OCCT exceptions (`SetBounds`,
//!   `Standard_ConstructionError`) and the UNPORTED arms.
//! - `Message_ProgressRange` / `Message_ProgressScope` (`cxx:471`, `:1875`) have
//!   no port; the `UserBreak()` early returns cannot trigger.
//! - `myInitOrder` is a `Vec<i32>` holding the 1-based permutation that
//!   `TrierTab` rewrites in place.
//! - `GeomAdaptor_Surface::UResolution` / `VResolution` are
//!   [`crate::approx_same_parameter::u_resolution`] /
//!   [`crate::approx_same_parameter::v_resolution`].
//! - **UNPORTED**: the `if (!Ok)` arm of `Perform` (`cxx:535-581`) builds a
//!   `GeomPlate_MakeApprox` from `myGeomPlateSurface`, which `Perform` nulled at
//!   `cxx:440`, so OCCT dereferences a null handle there. The port reports that
//!   instead of reproducing the crash.
//! - `Geom_RectangularTrimmedSurface(BAP.Plane(), u1, u2, v1, v2)` is
//!   [`GeomRectangularTrimmedSurface::uv`]; `Adaptor3d_CurveOnSurface` is
//!   [`CurveOnSurface`] over the `ProjLib_CompProjectedCurve` wrapper.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;

use crate::Surface;
use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::curve::Curve;
use crate::curve_on_surface::CurveOnSurface;
use crate::geom_plate::{
    Aij, BuildAveragePlane, MakeApprox, PlateSurface, half_space,
};
use crate::projlib::{CompProjectedCurve, ProjectedCurve2d};
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;

use super::BuildPlateSurface;

/// How far `Perform` got, one variant per early `return` of `cxx:438-747`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerformOutcome {
    /// `cxx:460-465`: both constraint lists are empty.
    NoConstraints,
    /// `cxx:507-510` (and the `ComputeSurfInit` failures that leave
    /// `mySurfInit` null): nothing to build on.
    NoInitSurface,
    /// `cxx:672-677`, `:711-717`: `myPlate.IsDone()` is false.
    PlateFailed,
    /// `cxx:747`: the loop converged.
    Done,
}

impl BuildPlateSurface {
    /// `CalculNbPtsInit()` (`cxx:2366-2404`): distributes `myNbPtsOnCur`
    /// points over the curves by length, with the `1.2` / `1` / `0.7` weights of
    /// the G0 / G1 / G2 orders, and never below three points.
    pub(super) fn calcul_nb_pts_init(&mut self) -> Result<(), String> {
        let nt_lin_cont = self.lin_cont.len();
        let nt_point = self.nb_pts_on_cur * nt_lin_cont as i32;
        // `cxx:2374-2377`.
        let mut len_t = 0.0;
        for i in 0..nt_lin_cont {
            len_t += self.lin(i).length()?;
        }
        // `cxx:2378-2399`.
        for i in 0..nt_lin_cont {
            let len = self.lin(i).length()?;
            let nb_points = match self.lin(i).order() {
                // `cxx:2382-2386`: case 0, G0 * 1.2.
                0 => (1.2 * nt_point as f64 * len / len_t) as i32,
                // `cxx:2387-2390`: case 1, G1 * 1.
                1 => (nt_point as f64 * len / len_t) as i32,
                // `cxx:2391-2394`: case 2, G2 * 0.7.
                2 => (0.7 * nt_point as f64 * len / len_t) as i32,
                // `cxx:2379-2395`: no `default:` arm, the count is untouched.
                _ => self.lin(i).nb_points(),
            };
            self.lin_mut(i).set_nb_points(nb_points);
            if self.lin(i).nb_points() < 3 {
                self.lin_mut(i).set_nb_points(3);
            }
        }
        Ok(())
    }

    /// `Perform()` (`cxx:438-747`): builds the initial surface, projects the
    /// constraints on it and iterates `Discretise` -> `Load*` -> `SolveTI` ->
    /// `VerifSurface` until the criteria hold.
    pub fn perform(&mut self) -> Result<PerformOutcome, String> {
        // `cxx:440`.
        self.geom_plate_surface = None;
        // `cxx:449-452`.
        if self.nb_bounds == 0 {
            self.nb_bounds = self.lin_cont.len() as i32;
        }
        // `cxx:454`.
        self.plate.init();
        let nt_lin_cont = self.lin_cont.len();
        let nt_pnt_cont = self.pnt_cont.len();
        // `cxx:460-465`.
        if nt_lin_cont + nt_pnt_cont == 0 {
            return Ok(PerformOutcome::NoConstraints);
        }
        // `cxx:474-505`: the initial surface.
        if !self.surf_init_is_give {
            self.compute_surf_init()?;
        } else if nt_lin_cont >= 2 {
            // `cxx:483-498`.
            self.init_order = (1..=nt_lin_cont as i32).collect();
            // `cxx:491`: a non-jointed curve set only prints a warning.
            self.courbe_jointive(self.tol3d)?;
            Self::trier_tab(&mut self.init_order);
        } else if nt_lin_cont > 0 {
            // `cxx:499-504`: a patch.
            self.sense = vec![0; nt_lin_cont];
            self.init_order = vec![1; nt_lin_cont];
        }
        // `cxx:507-510`.
        let surf = match self.surf_init.clone() {
            Some(s) => s,
            None => return Ok(PerformOutcome::NoInitSurface),
        };
        // `cxx:512-517`.
        self.tol_u = u_resolution(surf.as_ref(), self.tol3d);
        self.tol_v = v_resolution(surf.as_ref(), self.tol3d);
        // `cxx:522-534`: projection of the curves that have no p-curve yet.
        for i in 0..nt_lin_cont {
            if self.lin(i).curve2d_on_surf()?.is_none() {
                let curve3d = self
                    .lin(i)
                    .curve3d()
                    .ok_or_else(|| "GeomPlate_BuildPlateSurface: no 3D curve".to_string())?;
                match self.project_curve(&curve3d)? {
                    Some(curve2d) => self.lin_mut(i).set_curve2d_on_surf(curve2d),
                    None => {
                        // `cxx:535-581`: the fallback runs `GeomPlate_MakeApprox`
                        // on `myGeomPlateSurface`, which is null here.
                        return Err(
                            "UNPORTED GeomPlate_BuildPlateSurface.cxx:535-581: ProjectCurve \
                             failed and the MakeApprox fallback reads the null myGeomPlateSurface"
                                .to_string(),
                        );
                    }
                }
            }
        }
        // `cxx:584-596`: projection of the points.
        for i in 0..nt_pnt_cont {
            if !self.pnt(i).has_pnt2d_on_surf() {
                let p3d = self.pnt(i).d0();
                let p2d = self.project_point(&p3d)?;
                self.pnt_mut(i).set_pnt2d_on_surf(p2d);
            }
        }
        // `cxx:599-602`.
        if nt_lin_cont != 0 && self.nb_pts_on_cur != 0 {
            self.calcul_nb_pts_init()?;
        }
        // `cxx:608-616`.
        let mut pnt_inter: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
        let mut pnt_g1g1: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
        if nt_lin_cont != 0 {
            self.intersect(&mut pnt_inter, &mut pnt_g1g1)?;
        }
        // `cxx:620`.
        self.free = !self.is_linear;
        let mut nb_boucle = 0;
        loop {
            nb_boucle += 1;
            let mut fini;
            if nt_lin_cont != 0 {
                // `cxx:643-655`: `NPointMax` is computed and never read.
                self.discretise(&pnt_inter, &pnt_g1g1)?;
                // `cxx:655-659`.
                self.load_curve(nb_boucle, 2)?;
                if !self.pnt_cont.is_empty() {
                    self.load_point(2)?;
                }
                // `cxx:664`.
                self.plate
                    .solve_ti(self.degree, self.compute_anisotropie());
                // `cxx:672-677`.
                if !self.plate.is_done() {
                    return Ok(PerformOutcome::PlateFailed);
                }
                // `cxx:679-682`.
                self.set_plate_surface(&surf)?;
                // `cxx:684-692`.
                fini = self.verif_surface(nb_boucle)?;
                if nb_boucle >= self.nb_iter && !fini {
                    fini = true;
                }
                // `cxx:694-697`.
                if nt_pnt_cont != 0 && fini {
                    self.verif_points()?;
                }
            } else {
                // `cxx:700-728`.
                self.load_point(2)?;
                self.plate
                    .solve_ti(self.degree, self.compute_anisotropie());
                if !self.plate.is_done() {
                    return Ok(PerformOutcome::PlateFailed);
                }
                self.set_plate_surface(&surf)?;
                self.verif_points()?;
                fini = true;
            }
            if fini {
                break;
            }
        }
        Ok(PerformOutcome::Done)
    }

    /// `cxx:679-682` / `:720-723`: a new `GeomPlate_Surface` over the initial
    /// surface and the current plate solution, bounded by `Plate::UVBox`.
    fn set_plate_surface(&mut self, surf_init: &Arc<dyn Surface>) -> Result<(), String> {
        let mut plate_surface = PlateSurface::new(surf_init.clone(), self.plate.clone());
        let (umin, umax, vmin, vmax) = self.plate.uv_box();
        plate_surface.set_bounds(umin, umax, vmin, vmax)?;
        self.geom_plate_surface = Some(Arc::new(plate_surface));
        Ok(())
    }

    /// `ComputeSurfInit()` (`cxx:1454-1902`): the initial surface, either the
    /// half-space average plane of the G1 normals (`nopt = 3`) or the plane of
    /// inertia of a point cloud (`nopt = 1` / `2`), then - when the curve
    /// metrics do not match their projections - a first plate solve whose
    /// approximation becomes the new initial surface.
    pub(super) fn compute_surf_init(&mut self) -> Result<(), String> {
        // `cxx:1456-1458`.
        let mut nopt = 2;
        let popt = 2;
        // `cxx:1455`: `Np = 1` indexes the 1-based `Pts` array OCCT fills with
        // `Pts->ChangeValue(++pnum)`; the port pushes onto a `Vec` instead, so
        // the counter has no use here.
        let _np = 1usize;
        let mut is_half_space = true;
        let lin_tol = 0.001;
        let ang_tol = 0.001;
        let nt_lin_cont = self.lin_cont.len();
        let nt_pnt_cont = self.pnt_cont.len();
        // `cxx:1462-1470`.
        if nt_lin_cont != 0 {
            self.init_order = (1..=nt_lin_cont as i32).collect();
        }
        // `cxx:1472`.
        let courbe_joint = nt_lin_cont != 0 && self.courbe_jointive(self.tol3d)?;
        if courbe_joint && self.is_order_g1() {
            // `cxx:1476-1649`: the half-space / max-flux plane.
            nopt = 3;
            let nb_point = 20usize;
            let discr = nb_point / 4;
            let mut pts: Vec<GpPnt> =
                Vec::with_capacity((nb_point + 1) * nt_lin_cont + nt_pnt_cont);
            let mut vecs: Vec<GpVec> = Vec::new();
            let mut aset: Vec<Aij> = Vec::new();
            let mut last_vec = GpVec::zero();
            for i in 0..nt_lin_cont {
                let order = self.lin(i).order();
                let mut new_vecs: Vec<GpVec> = Vec::new();
                // `cxx:1492-1499`.
                let mut uinit = self.lin(i).first_parameter();
                let ufinal = self.lin(i).last_parameter();
                let mut uif = ufinal - uinit;
                if self.sense[i] == 1 {
                    uinit = ufinal;
                    uif = -uif;
                }
                // `cxx:1500-1512`.
                let mut normal = GpVec::zero();
                let to_reverse = if i > 0 && order >= 1 {
                    let (_, v1, v2) = self.lin(i).d1(uinit)?;
                    normal = v1.crossed(&v2);
                    last_vec.is_opposite_ang(&normal, ang_tol)
                } else {
                    false
                };
                // `cxx:1514-1556`: 20 points per curve, every `Discr` of them
                // contributing a normal.
                for j in 0..=nb_point {
                    let inter = (j as f64) * uif / (nb_point as f64);
                    if order < 1 || j % discr != 0 {
                        pts.push(self.lin(i).d0(uinit + inter));
                    } else {
                        let (p, v1, v2) = self.lin(i).d1(uinit + inter)?;
                        pts.push(p);
                        normal = v1.crossed(&v2).normalized();
                        if to_reverse {
                            normal.reverse();
                        }
                        let mut is_new = true;
                        for v in &vecs {
                            if v.is_equal_tol(&normal, lin_tol, ang_tol) {
                                is_new = false;
                                break;
                            }
                        }
                        if is_new {
                            for v in &new_vecs {
                                if v.is_equal_tol(&normal, lin_tol, ang_tol) {
                                    is_new = false;
                                    break;
                                }
                            }
                        }
                        if is_new {
                            new_vecs.push(normal);
                        }
                    }
                }
                // `cxx:1558-1564`.
                if order >= 1 {
                    is_half_space =
                        half_space(&new_vecs, &mut vecs, &mut aset, lin_tol, ang_tol);
                    if !is_half_space {
                        break;
                    }
                    last_vec = normal;
                }
            }
            if is_half_space {
                // `cxx:1571-1610`: the point constraints.
                for i in 0..nt_pnt_cont {
                    let order = self.pnt(i).order();
                    let mut new_vecs: Vec<GpVec> = Vec::new();
                    if order < 1 {
                        pts.push(self.pnt(i).d0());
                    } else {
                        let (p, v1, v2) = self.pnt(i).d1();
                        pts.push(p);
                        let normal = v1.crossed(&v2).normalized();
                        let mut is_new = true;
                        for v in &vecs {
                            if v.is_equal_tol(&normal, lin_tol, ang_tol) {
                                is_new = false;
                                break;
                            }
                        }
                        if is_new {
                            new_vecs.push(normal);
                            is_half_space =
                                half_space(&new_vecs, &mut vecs, &mut aset, lin_tol, ang_tol);
                            if !is_half_space {
                                new_vecs[0].reverse();
                                is_half_space = half_space(
                                    &new_vecs, &mut vecs, &mut aset, lin_tol, ang_tol,
                                );
                            }
                            if !is_half_space {
                                break;
                            }
                        }
                    }
                }
                if is_half_space {
                    // `cxx:1615-1627`: drop the null normals.
                    loop {
                        let mut null_exist = false;
                        let mut k = 0usize;
                        while k < vecs.len() {
                            if vecs[k].square_magnitude() == 0.0 {
                                null_exist = true;
                                vecs.remove(k);
                                break;
                            }
                            k += 1;
                        }
                        if !null_exist {
                            break;
                        }
                    }
                    // `cxx:1628-1639`.
                    let bap = BuildAveragePlane::from_normals(&vecs, pts.clone())?;
                    let (u1, u2, v1, v2) = bap.min_max_box();
                    let du = u2 - u1;
                    let dv = v2 - v1;
                    let trimmed = GeomRectangularTrimmedSurface::uv(
                        Arc::new(bap.plane()?.clone()),
                        u1 - du,
                        u2 + du,
                        v1 - dv,
                        v2 + dv,
                    )?;
                    self.surf_init = Some(Arc::new(trimmed));
                }
            }
            // `cxx:1642-1649`.
            if !is_half_space {
                self.is_linear = false;
                nopt = 2;
            }
        }
        // `cxx:1653-1656`.
        if nt_lin_cont != 0 {
            Self::trier_tab(&mut self.init_order);
        }
        if nopt != 3 {
            // `cxx:1659-1668`.
            if nt_pnt_cont != 0 {
                nopt = 1;
            } else if !courbe_joint || nt_lin_cont != self.nb_bounds.max(0) as usize {
                nopt = 1;
            }
            // `cxx:1671-1689`: the per-curve point count only feeds the debug
            // `Npt`, which OCCT also marks "(void) unused but set for debug".
            let mut len_t = 0.0;
            for i in 0..nt_lin_cont {
                len_t += self.lin(i).length()?;
            }
            let nt_point = 20 * nt_lin_cont;
            for i in 0..nt_lin_cont {
                let mut nb_point = (nt_point as f64 * self.lin(i).length()? / len_t) as i32;
                if nb_point < 10 {
                    nb_point = 10;
                }
                let _ = nb_point;
            }
            // `cxx:1691-1719`: the point cloud of the plane of inertia.
            let mut pts: Vec<GpPnt> = Vec::with_capacity(20 * nt_lin_cont + nt_pnt_cont);
            let nb_point = 20usize;
            for i in 0..nt_lin_cont {
                let mut uinit = self.lin(i).first_parameter();
                let ufinal = self.lin(i).last_parameter();
                let mut uif = ufinal - uinit;
                if self.sense[i] == 1 {
                    uinit = ufinal;
                    uif = -uif;
                }
                for j in 0..nb_point {
                    let inter = (j as f64) * uif / (nb_point as f64);
                    pts.push(self.lin(i).d0(uinit + inter));
                }
            }
            // `cxx:1713-1718`.
            for i in 0..nt_pnt_cont {
                pts.push(self.pnt(i).d0());
            }
            // `cxx:1720-1723`.
            if !courbe_joint {
                self.nb_bounds = 0;
            }
            let bap = BuildAveragePlane::from_points(
                pts,
                nb_point as i32 * self.nb_bounds,
                self.tol3d / 1000.0,
                popt,
                nopt,
            )?;
            // `cxx:1726-1731`.
            if !bap.is_plane() {
                return Ok(());
            }
            // `cxx:1732-1745`.
            let (u1, u2, v1, v2) = bap.min_max_box();
            let du = u2 - u1;
            let dv = v2 - v1;
            let trimmed = GeomRectangularTrimmedSurface::uv(
                Arc::new(bap.plane()?.clone()),
                u1 - du,
                u2 + du,
                v1 - dv,
                v2 + dv,
            )?;
            self.surf_init = Some(Arc::new(trimmed));
        }
        // `cxx:1748-1805`: comparing the metrics of the curves and of their
        // projections on the plane.
        if nt_lin_cont != 0 && self.is_linear {
            let init_plane = self
                .surf_init
                .as_ref()
                .and_then(|s| s.rectangular_trimmed_basis())
                .ok_or_else(|| {
                    "GeomPlate_BuildPlateSurface: the initial surface is not trimmed".to_string()
                })?;
            let mut ratio = 0.0;
            let r1 = 2.0;
            let r2 = 0.6;
            let nb_point = 20usize;
            for i in 0..nt_lin_cont {
                if !self.is_linear {
                    break;
                }
                let first_par = self.lin(i).first_parameter();
                let last_par = self.lin(i).last_parameter();
                let uif = (last_par - first_par) / (nb_point as f64);
                let curve = self
                    .lin(i)
                    .curve3d()
                    .ok_or_else(|| "GeomPlate_BuildPlateSurface: no 3D curve".to_string())?;
                let proj = CompProjectedCurve::with_tolerances(
                    init_plane.clone(),
                    curve.clone(),
                    self.tol3d,
                    self.tol3d,
                );
                let aproj = CurveOnSurface::new(
                    Arc::new(ProjectedCurve2d::new(proj)),
                    init_plane.clone(),
                );
                for j in 1..nb_point {
                    if !self.is_linear {
                        break;
                    }
                    let inter = first_par + (j as f64) * uif;
                    let (_, der_c) = curve.d1(inter);
                    let (_, der_proj) = aproj.d1(inter);
                    let a1 = der_c.magnitude();
                    let a2 = der_proj.magnitude();
                    ratio = if a2 <= 1.0e-20 { 1.0e20 } else { a1 / a2 };
                    if ratio > r1 || ratio < r2 {
                        self.is_linear = false;
                        break;
                    }
                }
            }
            let _ = ratio;
        }
        // `cxx:1806-1902`: the first plate solve when the metrics differ.
        if !self.is_linear {
            self.planar_surf_init = self.surf_init.clone();
            let surf = self
                .surf_init
                .clone()
                .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;
            self.tol_u = u_resolution(surf.as_ref(), self.tol3d);
            self.tol_v = v_resolution(surf.as_ref(), self.tol3d);
            // `cxx:1817-1825`.
            for i in 0..nt_lin_cont {
                if self.lin(i).curve2d_on_surf()?.is_none() {
                    let curve3d = self.lin(i).curve3d().ok_or_else(|| {
                        "GeomPlate_BuildPlateSurface: no 3D curve".to_string()
                    })?;
                    if let Some(curve2d) = self.project_curve(&curve3d)? {
                        self.lin_mut(i).set_curve2d_on_surf(curve2d);
                    }
                }
            }
            // `cxx:1828-1838`.
            for i in 0..nt_pnt_cont {
                let p3d = self.pnt(i).d0();
                if !self.pnt(i).has_pnt2d_on_surf() {
                    let p2d = self.project_point(&p3d)?;
                    self.pnt_mut(i).set_pnt2d_on_surf(p2d);
                }
            }
            // `cxx:1841-1844`.
            if nt_lin_cont != 0 && self.nb_pts_on_cur != 0 {
                self.calcul_nb_pts_init()?;
            }
            // `cxx:1848-1857`.
            let mut pnt_inter: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
            let mut pnt_g1g1: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
            if nt_lin_cont != 0 {
                self.intersect(&mut pnt_inter, &mut pnt_g1g1)?;
            }
            // `cxx:1862-1870`.
            self.discretise(&pnt_inter, &pnt_g1g1)?;
            self.load_curve(0, 0)?;
            if !self.pnt_cont.is_empty() {
                self.load_point(0)?;
            }
            // `cxx:1875-1887`.
            self.plate.solve_ti(2, self.compute_anisotropie());
            if !self.plate.is_done() {
                return Ok(());
            }
            // `cxx:1889-1897`.
            let plate_surface =
                Arc::new(PlateSurface::new(surf.clone(), self.plate.clone()));
            let app = MakeApprox::from_crit_order(
                plate_surface,
                self.tol3d,
                1,
                3,
                15.0 * self.tol3d,
                -1,
                Shape::C0,
                1.3,
            )?;
            self.surf_init = app
                .surface()
                .cloned()
                .map(|s| Arc::new(s) as Arc<dyn Surface>);
            self.surf_init_is_give = true;
            // `cxx:1896`.
            self.plate.init();
            // `cxx:1898-1901`: the p-curves no longer match the new surface.
            for i in 0..nt_lin_cont {
                self.lin_mut(i).clear_curve2d_on_surf();
            }
        }
        Ok(())
    }
}
