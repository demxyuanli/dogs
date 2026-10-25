//! `GeomFill_LocationDraft` — the location law of a drafted sweep: a
//! `GeomFill_DraftTrihedron` frame plus, when a stop surface is given, the two
//! 2d poles of the generatrix / stop-surface intersection.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationDraft.cxx:45-806` and
//! `GeomFill_LocationDraft.hxx:35-184` (OCCT 8.0.0).
//!
//! Layout:
//! - `mod.rs` (this file): the state, the constructor, `SetStopSurf` /
//!   `SetAngle` / `SetCurve` / `Copy` / `SetTrsf` / `Prepare`, the interval and
//!   norm queries, and the `GeomFill_LocationLaw` implementation,
//! - `eval.rs`: `D0` / `D1` / `D2` (`cxx:203-645`),
//! - `intersection.rs`: the generatrix / stop-surface intersection
//!   (`cxx:168-169`).
//!
//! Translation notes:
//! - `myLaw` is the concrete `GeomFill_DraftTrihedron` (`hxx:174`), so
//!   `SetAngle` reaches [`DraftTrihedron::set_angle`] directly.
//! - `myPoles2d` (`hxx:169`) is `(1, 2*myNbPts)`; OCCT index `i` is
//!   `poles2d[i - 1]` here, odd entries `(W, t)` on the curve and even entries
//!   `(U, V)` on the surface (`cxx:186-187`).
//! - `myTrimmed = myCurve->Trim(First, Last, 0)` (`cxx:695`) is
//!   [`GeomTrimmedCurveBasis`], the parameter-preserving `Geom_TrimmedCurve`.
//! - **Partial UNPORTED**: the generatrix / stop-surface intersection
//!   (`IntCurveSurface_HInter::Perform`) is ported for a planar stop surface
//!   in [`intersection`], the exact line / plane arm of
//!   `IntCurveSurface_Inter::PerformConicSurf`. Every other surface type
//!   returns [`intersection::UNPORTED_NON_PLANE`], because the faithful
//!   `IntCurveSurface_HInter` port lives in
//!   `occt_topo::intcurvesurface::perform_curve_surface` and `occt-geom` is the
//!   lower crate.

pub mod eval;
pub mod intersection;

use std::sync::{Arc, Mutex, MutexGuard};

use occt_core::gp::{GpDir, GpMat, GpPnt, GpPnt2d, GpVec};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::draft_trihedron::DraftTrihedron;
use crate::geom_fill::location_draft::intersection::{
    closest_intersection, line_surface_intersections,
};
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::trihedron_law::TrihedronLaw;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

/// `myNbPts` (`cxx:55`): the number of generatrix samples of `Prepare`.
pub const NB_PTS: usize = 41;

/// The mutable part of `GeomFill_LocationDraft` (`hxx:172-183`).
pub struct LdState {
    /// `Trans` (`hxx:173`).
    pub trans: GpMat,
    /// `myLaw` (`hxx:174`).
    pub law: DraftTrihedron,
    /// `mySurf` (`hxx:175`).
    pub surf: Option<Arc<dyn Surface>>,
    /// `myCurve` (`hxx:176`).
    pub curve: Option<Arc<dyn Curve>>,
    /// `myTrimmed` (`hxx:177`).
    pub trimmed: Option<Arc<dyn Curve>>,
    /// `myDir` (`hxx:178`).
    pub dir: GpDir,
    /// `myAngle` (`hxx:179`).
    pub angle: f64,
    /// `myNbPts` (`hxx:180`).
    pub nb_pts: usize,
    /// `Intersec` (`hxx:181`): "intersection avec surface d'arret ?".
    pub intersec: bool,
    /// `WithTrans` (`hxx:182`).
    pub with_trans: bool,
    /// `myPoles2d` (`hxx:169`): OCCT index `i` is `poles2d[i - 1]`.
    pub poles2d: Vec<GpPnt2d>,
}

/// `GeomFill_LocationDraft` (`GeomFill_LocationDraft.hxx:35-184`).
pub struct LocationDraft {
    st: Mutex<LdState>,
}

impl LocationDraft {
    /// `GeomFill_LocationDraft(Direction, Angle)` (`cxx:45-59`).
    pub fn new(direction: &GpDir, angle: f64) -> Self {
        let bi_normal = GpVec::from_xyz(direction.xyz());
        Self {
            st: Mutex::new(LdState {
                trans: GpMat::identity(),
                law: DraftTrihedron::new(&bi_normal, angle),
                surf: None,
                curve: None,
                trimmed: None,
                dir: *direction,
                angle,
                nb_pts: NB_PTS,
                intersec: false,
                with_trans: false,
                poles2d: vec![GpPnt2d::zero(); 2 * NB_PTS],
            }),
        }
    }

    /// The shared state (`hxx:169-183`).
    pub(crate) fn lock(&self) -> MutexGuard<'_, LdState> {
        self.st.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `SetStopSurf(Surf)` (`cxx:115-120`): stores the surface and runs
    /// `Prepare`. The `UNPORTED` intersection is reported here.
    pub fn set_stop_surf(&self, surf: Option<Arc<dyn Surface>>) -> Result<(), String> {
        {
            let mut st = self.lock();
            st.surf = surf;
        }
        self.prepare()
    }

    /// `SetAngle(Angle)` (`cxx:124-129`).
    pub fn set_angle(&self, angle: f64) -> Result<(), String> {
        {
            let mut st = self.lock();
            st.angle = angle;
            st.law.set_angle(angle);
        }
        self.prepare()
    }

    /// `Prepare()` (`cxx:135-195`): the 41 generatrix samples and, per sample,
    /// the poles of the generatrix / stop-surface intersection.
    pub fn prepare(&self) -> Result<(), String> {
        let (surf, curve, nb_pts, angle) = {
            let mut st = self.lock();
            if st.surf.is_none() {
                // `cxx:137-141`: no stop surface, no intersection poles.
                st.intersec = false;
                return Ok(());
            }
            // `cxx:144`: a stop surface is set.
            st.intersec = true;
            (st.surf.clone(), st.curve.clone(), st.nb_pts, st.angle)
        };
        let surf = surf.expect("surf checked above");
        let curve = curve.ok_or_else(|| {
            "Standard_NullObject: GeomFill_LocationDraft::Prepare".to_string()
        })?;

        let f = curve.first_parameter();
        let l = curve.last_parameter();
        let mut intersec = true;
        let mut poles2d = vec![GpPnt2d::zero(); 2 * nb_pts];
        for ii in 1..=nb_pts {
            // `cxx:157-158`: the sample parameter.
            let t = ((nb_pts - ii) as f64 * f + (ii - 1) as f64 * l) / (nb_pts - 1) as f64;
            let p = curve.d0(t);
            let (_, n, b) = {
                let st = self.lock();
                st.law.d0(t)?
            };
            // `cxx:163`: the generatrix.
            let d = b
                .multiplied_scalar(angle.cos())
                .added(&n.multiplied_scalar(angle.sin()));
            // `cxx:168-169`: `IntCurveSurface_HInter::Perform(AC, mySurf)`.
            let points = line_surface_intersections(&p, &d, surf.as_ref())?;
            if let Some((u, v, w)) = closest_intersection(&points) {
                // `cxx:184-187`.
                poles2d[2 * ii - 2] = GpPnt2d::new(w, t);
                poles2d[2 * ii - 1] = GpPnt2d::new(u, v);
            } else {
                // `cxx:189-192`: at least one sample has no intersection.
                intersec = false;
            }
        }
        let mut st = self.lock();
        st.poles2d = poles2d;
        st.intersec = intersec;
        Ok(())
    }

    /// `GetCurve()` (`cxx:198-201`).
    pub fn get_curve_handle(&self) -> Option<Arc<dyn Curve>> {
        self.lock().curve.clone()
    }

    /// `IsIntersec()` (`cxx:795-798`).
    pub fn is_intersec(&self) -> bool {
        self.lock().intersec
    }

    /// `Direction()` (`cxx:803-806`).
    pub fn direction(&self) -> GpDir {
        self.lock().dir
    }

    /// `D0` of the frame only, shared with `eval.rs`.
    pub(crate) fn frame(&self, param: f64, use_trimmed: bool) -> Result<(GpMat, GpVec), String> {
        let st = self.lock();
        let curve = if use_trimmed {
            st.trimmed.clone()
        } else {
            st.curve.clone()
        }
        .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D0".to_string())?;
        let p = curve.d0(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let (t, n, b) = st.law.d0(param)?;
        let mut m = GpMat::from_cols(n.xyz(), b.xyz(), t.xyz());
        if st.with_trans {
            m = m.multiply(&st.trans);
        }
        Ok((m, v))
    }
}

impl LocationLaw for LocationDraft {
    /// `SetCurve(C)` (`cxx:106-112`): `myCurve = myTrimmed = C` and
    /// `myLaw->SetCurve(C)`, then `Prepare()`.
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        let ok = {
            let mut st = self.lock();
            st.curve = Some(curve.clone());
            st.trimmed = Some(curve.clone());
            st.law.set_curve(curve)
        };
        ok && self.prepare().is_ok()
    }

    /// `GetCurve()` (`cxx:198-201`).
    fn get_curve(&self) -> Option<Arc<dyn Curve>> {
        self.lock().curve.clone()
    }

    /// `SetTrsf(Transfo)` (`cxx:81-97`): `WithTrans` is false only when every
    /// entry of `I - Transfo` is `<= 1.e-14` in magnitude.
    fn set_trsf(&self, transfo: &GpMat) {
        let mut st = self.lock();
        st.trans = *transfo;
        let aux = GpMat::identity().subtract(transfo);
        let mut with_trans = false;
        for ii in 1..=3 {
            for jj in 1..=3 {
                if !with_trans && aux.value(ii, jj).abs() > 1.0e-14 {
                    with_trans = true;
                }
            }
        }
        st.with_trans = with_trans;
    }

    /// `Copy()` (`cxx:62-77`): `myLaw->Copy()` is computed and discarded
    /// (`cxx:65-66`), the copy rebuilds its own `DraftTrihedron`.
    fn copy(&self) -> Arc<dyn LocationLaw> {
        let (dir, angle, curve, surf, with_trans, trans) = {
            let st = self.lock();
            (
                st.dir,
                st.angle,
                st.curve.clone(),
                st.surf.clone(),
                st.with_trans,
                st.trans,
            )
        };
        let copy = LocationDraft::new(&dir, angle);
        if let Some(curve) = curve {
            copy.set_curve(curve);
        }
        copy.set_stop_surf(surf).ok();
        if with_trans {
            copy.set_trsf(&trans);
        }
        Arc::new(copy)
    }

    /// `D0(Param, M, V)` (`cxx:203-229`).
    fn d0(&self, param: f64) -> Result<(GpMat, GpVec), String> {
        self.frame(param, true)
    }

    /// `D0(Param, M, V, Poles2d)` (`cxx:233-341`).
    fn d0_with_poles(&self, param: f64) -> Result<(GpMat, GpVec, Vec<GpPnt2d>), String> {
        self.d0_poles(param)
    }

    /// `D1` (`cxx:346-478`).
    fn d1(&self, param: f64) -> Result<crate::geom_fill::location_law::LocationD1, String> {
        self.d1_frame(param)
    }

    /// `D2` (`cxx:483-645`).
    fn d2(&self, param: f64) -> Result<crate::geom_fill::location_law::LocationD2, String> {
        self.d2_frame(param)
    }

    /// `HasFirstRestriction()` (`cxx:647-650`).
    fn has_first_restriction(&self) -> bool {
        false
    }

    /// `HasLastRestriction()` (`cxx:653-658`): `Intersec`.
    fn has_last_restriction(&self) -> bool {
        self.lock().intersec
    }

    /// `TraceNumber()` (`cxx:660-673`).
    fn trace_number(&self) -> i32 {
        if self.lock().intersec {
            1
        } else {
            0
        }
    }

    /// `NbIntervals(S)` (`cxx:675-680`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let st = self.lock();
        st.law.nb_intervals(s)
    }

    /// `Intervals(T, S)` (`cxx:682-688`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let st = self.lock();
        st.law.intervals(s)
    }

    /// `SetInterval(First, Last)` (`cxx:690-695`).
    fn set_interval(&self, first: f64, last: f64) {
        let mut st = self.lock();
        st.law.set_interval(first, last).ok();
        if let Some(curve) = st.curve.clone() {
            st.trimmed = Some(Arc::new(GeomTrimmedCurveBasis::new(curve, first, last)));
        }
    }

    /// `GetInterval` (`cxx:697-703`).
    fn get_interval(&self) -> (f64, f64) {
        let st = self.lock();
        match &st.trimmed {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        }
    }

    /// `GetDomain` (`cxx:705-712`).
    fn get_domain(&self) -> (f64, f64) {
        let st = self.lock();
        match &st.curve {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        }
    }

    /// `Resolution(Index, Tol, TolU, TolV)` (`cxx:714-732`): index 1 asks the
    /// stop surface for its parametric resolutions.
    fn resolution(&self, index: i32, tol: f64) -> Result<(f64, f64), String> {
        if index == 1 {
            let surf = self.lock().surf.clone();
            let surf = surf
                .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::Resolution".to_string())?;
            return surf
                .uv_resolution(tol)
                .ok_or_else(|| {
                    "UNPORTED: Adaptor3d_Surface::UResolution / VResolution for this surface".to_string()
                });
        }
        Ok((tol, tol))
    }

    /// `GetMaximalNorm()` (`cxx:738-741`): the trihedrons are normalized.
    fn get_maximal_norm(&self) -> f64 {
        1.0
    }

    /// `GetAverageLaw(AM, AV)` (`cxx:744-762`).
    fn get_average_law(&self) -> (GpMat, GpVec) {
        let st = self.lock();
        let (v1, v2, v3) = st.law.get_average_law();
        let am = GpMat::from_cols(v1.xyz(), v2.xyz(), v3.xyz());
        let Some(trimmed) = &st.trimmed else {
            return (am, GpVec::zero());
        };
        let mut av = GpVec::zero();
        let delta = (trimmed.last_parameter() - trimmed.first_parameter()) / 10.0;
        let mut u = trimmed.first_parameter();
        for _ in 0..=10 {
            let p = trimmed.d0(u);
            av = av.add(&GpVec::new(p.x(), p.y(), p.z()));
            u += delta;
        }
        (am, av.multiply_scalar(1.0 / 11.0))
    }

    /// `IsTranslation(Error)` (`cxx:766-770`): `myLaw->IsConstant()`.
    fn is_translation(&self, _error: &mut f64) -> bool {
        let st = self.lock();
        st.law.is_constant()
    }

    /// `IsRotation(Error)` (`cxx:773-784`): a circle whose draft trihedron only
    /// depends on the 3d curve.
    fn is_rotation(&self, error: &mut f64) -> bool {
        *error = 0.0;
        let st = self.lock();
        let Some(curve) = &st.curve else {
            return false;
        };
        if extrema_curve_type(curve.as_ref()) == GeomAbsCurveType::Circle {
            return st.law.is_only_by_3d_curve();
        }
        false
    }

    /// `Rotation(Centre)` (`cxx:787-792`): `myCurve->Circle().Location()`.
    fn rotation(&self) -> Result<GpPnt, String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::Rotation".to_string())?;
        curve
            .gp_circ()
            .map(|c| c.location())
            .ok_or_else(|| "Standard_NoSuchObject: GeomFill_LocationDraft::Rotation".to_string())
    }
}
