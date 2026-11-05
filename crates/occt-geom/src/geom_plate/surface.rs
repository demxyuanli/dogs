//! `GeomPlate_Surface` — the `Plate_Plate` solution as a `Geom_Surface`.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_Surface.hxx`,
//! `GeomPlate_Surface.cxx:42-309` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Plate_Plate mySurfinter` is held behind a `Mutex` (see
//!   [`PlateSurface::plate`]). OCCT's `Plate_Plate` memoises `SolEm` through
//!   `const_cast` and is therefore not thread-safe; the port's [`Plate`] keeps
//!   that memoisation in `Cell`s, which are `Send` but not `Sync`, so the
//!   `Surface` bound (`Send + Sync`) is met by serialising access. Every method
//!   is `const` in OCCT, so the lock never nests.
//! - `Bounds()` (`cxx:97-111`) branches on
//!   `mySurfinit->DynamicType() == STANDARD_TYPE(GeomPlate_Surface)`; that test
//!   is [`Surface::as_plate_surface`], and the branch delegates to the init
//!   surface's own bounds, exactly as OCCT does.
//! - `EvalD3` / `EvalDN` throw `Geom_UndefinedDerivative` (`cxx:241-251`); the
//!   port reports them through [`Surface::try_eval_dn`].
//! - `Copy()` (`cxx:254-258`) rebuilds the surface from the same init surface
//!   and plate solution, so the `SetBounds` values are **not** copied (the
//!   constructor resets them to zero, `cxx:43-51`). [`PlateSurface::copy`]
//!   keeps that.

use std::sync::{Arc, Mutex, MutexGuard};

use occt_core::gp::{GpPnt, GpVec, GpTrsf, GpXY};

use crate::plate::Plate;
use crate::surface::Surface;

/// `class GeomPlate_Surface : public Geom_Surface` (`hxx:41-256`).
pub struct PlateSurface {
    /// `mySurfinter` (`hxx:248`).
    plate: Mutex<Plate>,
    /// `mySurfinit` (`hxx:249`).
    surfinit: Arc<dyn Surface>,
    /// `myUmin` (`hxx:250`).
    u_min: f64,
    /// `myUmax` (`hxx:251`).
    u_max: f64,
    /// `myVmin` (`hxx:252`).
    v_min: f64,
    /// `myVmax` (`hxx:253`).
    v_max: f64,
}

impl PlateSurface {
    /// `GeomPlate_Surface(Surfinit, Surfinter)` (`cxx:42-51`).
    pub fn new(surfinit: Arc<dyn Surface>, surfinter: Plate) -> Self {
        Self {
            plate: Mutex::new(surfinter),
            surfinit,
            u_min: 0.0,
            u_max: 0.0,
            v_min: 0.0,
            v_max: 0.0,
        }
    }

    /// The `Plate_Plate` solution (`mySurfinter`). OCCT exposes it only
    /// internally; the plate builder reads it back to re-solve.
    pub fn plate(&self) -> MutexGuard<'_, Plate> {
        self.plate.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `CallSurfinit()` (`cxx:270-273`).
    pub fn call_surfinit(&self) -> Arc<dyn Surface> {
        self.surfinit.clone()
    }

    /// `SetBounds(Umin, Umax, Vmin, Vmax)` (`cxx:277-296`).
    ///
    /// Raises when the order or the span is degenerate (`cxx:283-291`).
    pub fn set_bounds(
        &mut self,
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
    ) -> Result<(), String> {
        if u_min > u_max || v_min > v_max {
            return Err("Bounds haven't the good sense".to_string());
        }
        if u_min == u_max || v_min == v_max {
            return Err("Bounds are equal".to_string());
        }
        self.u_min = u_min;
        self.u_max = u_max;
        self.v_min = v_min;
        self.v_max = v_max;
        Ok(())
    }

    /// `Bounds(U1, U2, V1, V2)` (`cxx:97-111`).
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        if let Some(inner) = self.surfinit.as_plate_surface() {
            return inner.bounds();
        }
        (self.u_min, self.u_max, self.v_min, self.v_max)
    }

    /// `RealBounds(U1, U2, V1, V2)` (`cxx:299-302`): `mySurfinter.UVBox`.
    pub fn real_bounds(&self) -> (f64, f64, f64, f64) {
        self.plate().uv_box()
    }

    /// `Constraints(Seq)` (`cxx:306-309`): `mySurfinter.UVConstraints`.
    pub fn constraints(&self) -> Vec<GpXY> {
        self.plate().uv_constraints()
    }

    /// `Copy()` (`cxx:254-258`): a new `GeomPlate_Surface` over the same init
    /// surface and plate solution, with the `SetBounds` values reset.
    pub fn copy(&self) -> Self {
        Self::new(self.surfinit.clone(), self.plate().clone())
    }
}

impl Surface for PlateSurface {
    /// `UReverse()` (`cxx:55-58`) does nothing in OCCT, so the reversed surface
    /// is a plain copy.
    fn u_reversed(&self) -> Option<Arc<dyn Surface>> {
        Some(Arc::new(self.copy()))
    }

    /// `UReversedParameter(U)` (`cxx:62-65`): `-U`.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        -u
    }

    /// `VReverse()` (`cxx:69-71`) does nothing in OCCT.
    fn v_reversed(&self) -> Option<Arc<dyn Surface>> {
        Some(Arc::new(self.copy()))
    }

    /// `VReversedParameter(V)` (`cxx:75-78`): `-V`.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }

    /// `Bounds` (`cxx:97-111`).
    fn u_range(&self) -> (f64, f64) {
        let (u1, u2, _, _) = self.bounds();
        (u1, u2)
    }

    /// `Bounds` (`cxx:97-111`).
    fn v_range(&self) -> (f64, f64) {
        let (_, _, v1, v2) = self.bounds();
        (v1, v2)
    }

    /// `IsUClosed()` (`cxx:115-119`) returns `false`.
    fn is_u_closed(&self) -> bool {
        false
    }

    /// `IsVClosed()` (`cxx:123-127`) returns `false`.
    fn is_v_closed(&self) -> bool {
        false
    }

    /// `IsUPeriodic()` (`cxx:131-134`) returns `false`.
    fn is_u_periodic(&self) -> bool {
        false
    }

    /// `UPeriod()` (`cxx:138-141`) returns `false`, i.e. zero.
    fn u_period(&self) -> f64 {
        0.0
    }

    /// `IsVPeriodic()` (`cxx:145-148`) returns `false`.
    fn is_v_periodic(&self) -> bool {
        false
    }

    /// `VPeriod()` (`cxx:152-155`) returns `false`, i.e. zero.
    fn v_period(&self) -> f64 {
        0.0
    }

    /// `UIso(U)` (`cxx:159-162`) returns a null handle.
    fn u_iso_curve(&self, _u: f64) -> Option<Arc<dyn crate::curve::Curve>> {
        None
    }

    /// `VIso(V)` (`cxx:166-170`) returns a null handle.
    fn v_iso_curve(&self, _v: f64) -> Option<Arc<dyn crate::curve::Curve>> {
        None
    }

    /// `Continuity()` (`cxx:175-178`) returns a default `GeomAbs_Shape`, i.e.
    /// `GeomAbs_C0`.
    fn continuity(&self) -> u8 {
        0
    }

    /// `IsCNu(N)` (`cxx:182-185`) raises.
    fn is_cnu(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_Surface::IsCNu: negative order".into());
        }
        Err("GeomPlate_Surface::IsCNu".into())
    }

    /// `IsCNv(N)` (`cxx:189-192`) raises.
    fn is_cnv(&self, n: i32) -> Result<bool, String> {
        if n < 0 {
            return Err("Geom_Surface::IsCNv: negative order".into());
        }
        Err("GeomPlate_Surface::IsCNv".into())
    }

    /// `EvalD0` (`cxx:197-204`): `Surfinit->EvalD0(U, V) + Surfinter.Evaluate`.
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        let base = self.surfinit.d0(u, v);
        let delta = self.plate().evaluate(GpXY::new(u, v));
        GpPnt::from_xyz(&base.xyz().added(&delta))
    }

    /// `EvalD1` (`cxx:208-217`).
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1u, d1v) = self.surfinit.d1(u, v);
        let point2d = GpXY::new(u, v);
        let (delta, du, dv) = {
            let plate = self.plate();
            (
                plate.evaluate(point2d),
                plate.evaluate_derivative(point2d, 1, 0),
                plate.evaluate_derivative(point2d, 0, 1),
            )
        };
        (
            GpPnt::from_xyz(&p.xyz().added(&delta)),
            GpVec::from_xyz(&d1u.xyz().added(&du)),
            GpVec::from_xyz(&d1v.xyz().added(&dv)),
        )
    }

    /// `EvalD2` (`cxx:221-236`).
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        let (p, d1u, d1v, d2u, d2v, d2uv) = self.surfinit.d2(u, v);
        let point2d = GpXY::new(u, v);
        let (delta, du, dv, d2u_i, d2v_i, d2uv_i) = {
            let plate = self.plate();
            (
                plate.evaluate(point2d),
                plate.evaluate_derivative(point2d, 1, 0),
                plate.evaluate_derivative(point2d, 0, 1),
                plate.evaluate_derivative(point2d, 2, 0),
                plate.evaluate_derivative(point2d, 0, 2),
                plate.evaluate_derivative(point2d, 1, 1),
            )
        };
        (
            GpPnt::from_xyz(&p.xyz().added(&delta)),
            GpVec::from_xyz(&d1u.xyz().added(&du)),
            GpVec::from_xyz(&d1v.xyz().added(&dv)),
            GpVec::from_xyz(&d2u.xyz().added(&d2u_i)),
            GpVec::from_xyz(&d2v.xyz().added(&d2v_i)),
            GpVec::from_xyz(&d2uv.xyz().added(&d2uv_i)),
        )
    }

    /// `EvalDN` (`cxx:248-251`) throws `Geom_UndefinedDerivative`.
    fn eval_dn(&self, _u: f64, _v: f64, _nu: i32, _nv: i32) -> GpVec {
        GpVec::zero()
    }

    /// `EvalDN` (`cxx:248-251`) throws `Geom_UndefinedDerivative`; `EvalD3`
    /// (`cxx:241-245`) throws the same. Both report the missing derivative so
    /// callers never read the zero above as a real one.
    fn try_eval_dn(&self, _u: f64, _v: f64, _nu: i32, _nv: i32) -> Result<GpVec, String> {
        Err("GeomPlate_Surface::EvalDN".to_string())
    }

    /// `Transform(T)` (`cxx:261-263`) does nothing.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Copy()` (`cxx:254-258`).
    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.copy())
    }

    /// `mySurfinit->DynamicType() == STANDARD_TYPE(GeomPlate_Surface)`
    /// (`cxx:99`).
    fn as_plate_surface(&self) -> Option<&PlateSurface> {
        Some(self)
    }
}
