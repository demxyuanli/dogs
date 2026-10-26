//! `ProjLib_CompProjectedCurve::ShallowCopy` / `Perform`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx`
//! `ShallowCopy` `cxx:630-656`, `Perform` `cxx:1229-1408` (OCCT 8.0.0).
//!
//! `Perform` walks the continuous parts built by [`super::init::Init`]:
//!
//! - a punctual part becomes a 3D point on the surface (`cxx:1269-1284`),
//! - a `U`-isoparametric part becomes a trimmed `v` line at the constant `u`,
//!   a `V`-isoparametric part a trimmed `u` line at the constant `v`
//!   (`cxx:1289-1327`), and both switch the 2D approximation off,
//! - every other part is approximated by `Approx_CurveOnSurface` over a
//!   shallow copy of the projection whose surface has its parametric domain
//!   widened by `10 * TolU` / `10 * TolV` on the non-periodic sides
//!   (`cxx:1328-1375`).
//!
//! Translation notes:
//! - `HS = HS->UTrim(U1, U2, 0.0)` / `VTrim` (`cxx:1333-1352`) is
//!   `GeomAdaptor_Surface::UTrim`/`VTrim`
//!   (`GeomAdaptor_Surface.cxx:886-905`), i.e. a `GeomAdaptor_Surface` window
//!   with no tightness normalisation: [`GeomRectangularTrimmedSurface::window`].
//! - `HPCur = new Geom2dAdaptor_Curve(PCur2d)` (`cxx:1324`, `:1345`) resolves
//!   to the trimmed curve's basis for every type query
//!   (`Geom2dAdaptor_Curve.cxx:96-120`), which is exactly what
//!   [`Geom2dTrimmedCurve`] forwards, so the trimmed curve itself is used.
//! - `ShallowCopy` shares the handles in OCCT and copies them here; the
//!   `mySequence` samples are read-only, so the behaviour is the same.
//! - `approx2d` / `approx3d` are function-scope in OCCT and are switched off
//!   by the first isoparametric part for the **remaining** parts, as are the
//!   loop-carried `PCur2d` / `PCur3d` handles; both are reproduced verbatim.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt, GpPnt2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::approx::curve_on_surface::ApproxCurveOnSurface;
use crate::curve::Curve;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;

use super::adaptor::ProjectedCurve2d;
use super::CompProjectedCurve;

impl CompProjectedCurve {
    /// `ShallowCopy()` (`cxx:632-654`): a copy sharing the surface, the
    /// curve and the sample sequence, and **not** carrying `myTabInt` (which
    /// the fresh `ProjLib_CompProjectedCurve()` never set).
    pub fn shallow_copy(&self) -> Self {
        Self {
            surface: self.surface.clone(),
            curve: self.curve.clone(),
            nb_curves: self.nb_curves,
            sequence: self.sequence.clone(),
            ui_so: self.ui_so.clone(),
            v_iso: self.v_iso.clone(),
            sngl_pnts: self.sngl_pnts.clone(),
            max_distance: self.max_distance.clone(),
            tab_int: None,
            tol3d: self.tol3d,
            continuity: self.continuity,
            max_degree: self.max_degree,
            max_seg: self.max_seg,
            proj2d: self.proj2d,
            proj3d: self.proj3d,
            max_dist: self.max_dist,
            tol_u: self.tol_u,
            tol_v: self.tol_v,
            result_is_point: Vec::new(),
            result_2d_u_approx_error: Vec::new(),
            result_2d_v_approx_error: Vec::new(),
            result_3d_approx_error: Vec::new(),
            result_3d_point: Vec::new(),
            result_2d_point: Vec::new(),
            result_3d_curve: Vec::new(),
            result_2d_curve: Vec::new(),
        }
    }

    /// `Perform()` (`cxx:1229-1408`).
    pub fn perform(&mut self) -> Result<(), String> {
        if self.nb_curves == 0 {
            return Ok(());
        }

        let mut approx2d = self.proj2d;
        let approx3d = self.proj3d;
        let nb = self.nb_curves as usize;

        // `cxx:1245-1267`: the result arrays.
        if self.proj2d {
            self.result_2d_point = vec![GpPnt2d::new(0.0, 0.0); nb];
            self.result_2d_curve = vec![None; nb];
        }
        if self.proj3d {
            self.result_3d_point = vec![GpPnt::zero(); nb];
            self.result_3d_curve = vec![None; nb];
        }
        self.result_is_point = vec![false; nb];
        self.result_3d_approx_error = vec![0.0; nb];
        self.result_2d_u_approx_error = vec![0.0; nb];
        self.result_2d_v_approx_error = vec![0.0; nb];

        let surface = self
            .surface
            .clone()
            .ok_or_else(|| "ProjLib_CompProjectedCurve::Perform: no surface".to_string())?;
        let s = surface.as_ref();
        let (fu1, fu2) = (s.u_range().0, s.u_range().1);
        let (fv1, fv2) = (s.v_range().0, s.v_range().1);
        // `HS` and the window it carries (`cxx:1241`, `:1333-1352`).
        let mut hsurf: Arc<dyn Surface> = surface.clone();
        let (mut wu1, mut wu2, mut wv1, mut wv2) = (fu1, fu2, fv1, fv2);
        // `PCur2d` / `PCur3d` are declared outside the loop (`cxx:1243-1244`).
        let mut pcur2d: Option<Arc<dyn Curve2d>> = None;
        let mut pcur3d: Option<Arc<dyn Curve>> = None;

        for k in 1..=self.nb_curves {
            let (p2d, is_pnt) = self
                .is_single_pnt(k)
                .ok_or_else(|| "ProjLib_CompProjectedCurve::Perform: no such part".to_string())?;
            if is_pnt {
                // `cxx:1269-1284`.
                let p = s.d0(p2d.x(), p2d.y());
                if self.proj2d {
                    self.result_2d_point[k as usize - 1] = p2d;
                }
                if self.proj3d {
                    self.result_3d_point[k as usize - 1] = p;
                }
                self.result_is_point[k as usize - 1] = true;
                continue;
            }

            let (udeb, ufin) = self
                .bounds(k)
                .ok_or_else(|| "ProjLib_CompProjectedCurve::Perform: no such part".to_string())?;
            let (mut udeb, mut ufin) = (udeb, ufin);

            let (u_iso, is_u_iso) = self
                .is_u_iso(k)
                .ok_or_else(|| "ProjLib_CompProjectedCurve::Perform: no such part".to_string())?;
            let (v_iso, is_v_iso) = self
                .is_v_iso(k)
                .ok_or_else(|| "ProjLib_CompProjectedCurve::Perform: no such part".to_string())?;

            let h_pcur: Arc<dyn Curve2d>;
            if is_u_iso {
                // `cxx:1289-1308`: the part is the `u = UIso` isoparametric
                // curve, parametrised by `v`.
                approx2d = false;
                let pdeb = self.d0(udeb)?;
                let pfin = self.d0(ufin)?;
                udeb = pdeb.y();
                ufin = pfin.y();
                let dir;
                if udeb > ufin {
                    dir = GpDir2d { x: 0.0, y: -1.0 };
                    udeb = -udeb;
                    ufin = -ufin;
                } else {
                    dir = GpDir2d { x: 0.0, y: 1.0 };
                }
                let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_iso, 0.0), dir);
                let trimmed = Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), udeb, ufin));
                pcur2d = Some(trimmed.clone());
                h_pcur = trimmed;
            } else if is_v_iso {
                // `cxx:1309-1327`: the part is the `v = VIso` isoparametric
                // curve, parametrised by `u`.
                approx2d = false;
                let pdeb = self.d0(udeb)?;
                let pfin = self.d0(ufin)?;
                udeb = pdeb.x();
                ufin = pfin.x();
                let dir;
                if udeb > ufin {
                    dir = GpDir2d { x: -1.0, y: 0.0 };
                    udeb = -udeb;
                    ufin = -ufin;
                } else {
                    dir = GpDir2d { x: 1.0, y: 0.0 };
                }
                let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v_iso), dir);
                let trimmed = Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), udeb, ufin));
                pcur2d = Some(trimmed.clone());
                h_pcur = trimmed;
            } else {
                // `cxx:1328-1375`: widen the domain so a projection touching
                // the border is inside it, then approximate against a shallow
                // copy of this projection.
                if !s.is_u_periodic() {
                    let du = 10.0 * self.tol_u;
                    wu1 = fu1 - du;
                    wu2 = fu2 + du;
                    hsurf = Arc::new(GeomRectangularTrimmedSurface::window(
                        hsurf.clone(),
                        wu1,
                        wu2,
                        wv1,
                        wv2,
                    )?);
                }
                if !s.is_v_periodic() {
                    let dv = 10.0 * self.tol_v;
                    wv1 = fv1 - dv;
                    wv2 = fv2 + dv;
                    hsurf = Arc::new(GeomRectangularTrimmedSurface::window(
                        hsurf.clone(),
                        wu1,
                        wu2,
                        wv1,
                        wv2,
                    )?);
                }
                let mut hp = self.shallow_copy();
                hp.load_surface(hsurf.clone());
                h_pcur = Arc::new(ProjectedCurve2d::new(hp));
            }

            if approx2d || approx3d {
                let (only2d, only3d) = if approx2d && approx3d {
                    (false, false)
                } else {
                    (approx2d, approx3d)
                };
                let mut appr =
                    ApproxCurveOnSurface::new(h_pcur, hsurf.clone(), udeb, ufin, self.tol3d);
                appr.perform(
                    self.max_seg,
                    self.max_degree,
                    self.continuity,
                    only3d,
                    only2d,
                )?;
                if approx2d {
                    pcur2d = appr
                        .curve2d()
                        .map(|c| Arc::new(c.clone()) as Arc<dyn Curve2d>);
                    self.result_2d_u_approx_error[k as usize - 1] = appr.max_error2du();
                    self.result_2d_v_approx_error[k as usize - 1] = appr.max_error2dv();
                }
                if approx3d {
                    pcur3d = appr
                        .curve3d()
                        .map(|c| Arc::new(c.clone()) as Arc<dyn Curve>);
                    self.result_3d_approx_error[k as usize - 1] = appr.max_error3d();
                }
            }

            if self.proj2d {
                self.result_2d_curve[k as usize - 1] = pcur2d.clone();
            }
            if self.proj3d {
                self.result_3d_curve[k as usize - 1] = pcur3d.clone();
            }
        }
        Ok(())
    }
}
