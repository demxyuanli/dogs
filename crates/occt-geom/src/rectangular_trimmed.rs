//! Rectangular trimmed surface. Source: `Geom_RectangularTrimmedSurface.hxx`

use std::sync::Arc;

use occt_core::gp::{GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpTrsf, GpVec};
use occt_core::precision::Precision;

use crate::curve::Curve;
use crate::geom_adaptor_local::{try_local_d1, try_local_d2};
use crate::offset_surface::GeomOffsetSurface;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurve;

/// Portion of a surface limited by U and/or V parameter bounds.
/// Source: `Geom_RectangularTrimmedSurface.cxx`.
#[derive(Clone)]
pub struct GeomRectangularTrimmedSurface {
    basis: Arc<dyn Surface>,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    u_trimmed: bool,
    v_trimmed: bool,
}

/// Adaptor domain tolerance used by `UTrim`/`VTrim` (`Adaptor3d_CurveOnSurface`
/// EvalFirstLastSurf passes `Precision::PConfusion`).
const ADAPTOR_TOL: f64 = Precision::PCONFUSION;

fn copy_untrimmed(s: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    s.rectangular_trimmed_basis()
        .unwrap_or_else(|| Arc::from(s.clone_dyn()))
}

fn ieee_remainder(x: f64, y: f64) -> f64 {
    if !y.is_finite() || y.abs() < f64::EPSILON {
        return x;
    }
    x - y * (x / y).round()
}

fn closed_when_trimmed(periodic: bool, period: f64, a: f64, b: f64) -> bool {
    if !periodic {
        return false;
    }
    let length = b - a;
    length > Precision::PCONFUSION
        && ieee_remainder(length, period).abs() <= Precision::PCONFUSION
}

impl GeomRectangularTrimmedSurface {
    /// UV trim. Source: `Geom_RectangularTrimmedSurface.cxx:73-111`.
    pub fn uv(s: Arc<dyn Surface>, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        let mut basis = copy_untrimmed(&s);
        if let (Some(ob), Some(d)) = (basis.offset_basis_surface(), basis.offset_distance()) {
            let inner = Arc::new(Self::uv_raw(ob, u1, u2, v1, v2));
            basis = Arc::new(GeomOffsetSurface::new(inner, d));
        }
        Self::uv_raw(basis, u1, u2, v1, v2)
    }

    fn uv_raw(basis: Arc<dyn Surface>, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        Self {
            basis,
            u1,
            u2,
            v1,
            v2,
            u_trimmed: true,
            v_trimmed: true,
        }
    }

    /// Trim in one parametric direction (`UTrim = true` trims U).
    /// Source: `Geom_RectangularTrimmedSurface.cxx:115-152`.
    pub fn one_param(s: Arc<dyn Surface>, param1: f64, param2: f64, u_trim: bool) -> Self {
        let mut basis = copy_untrimmed(&s);
        if let (Some(ob), Some(d)) = (basis.offset_basis_surface(), basis.offset_distance()) {
            let inner = Arc::new(Self::one_param_raw(ob, param1, param2, u_trim));
            basis = Arc::new(GeomOffsetSurface::new(inner, d));
        }
        Self::one_param_raw(basis, param1, param2, u_trim)
    }

    fn one_param_raw(basis: Arc<dyn Surface>, param1: f64, param2: f64, u_trim: bool) -> Self {
        let (bu1, bu2) = basis.u_range();
        let (bv1, bv2) = basis.v_range();
        if u_trim {
            Self {
                basis,
                u1: param1,
                u2: param2,
                v1: bv1,
                v2: bv2,
                u_trimmed: true,
                v_trimmed: false,
            }
        } else {
            Self {
                basis,
                u1: bu1,
                u2: bu2,
                v1: param1,
                v2: param2,
                u_trimmed: false,
                v_trimmed: true,
            }
        }
    }

    /// The basis surface (`Geom_RectangularTrimmedSurface::BasisSurface`).
    pub fn basis_surface(&self) -> &Arc<dyn Surface> {
        &self.basis
    }
}

impl Surface for GeomRectangularTrimmedSurface {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.basis.d0(u, v)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        // `GeomAdaptor_Surface::EvalD1` BSpline LocalD1 when UV is on the
        // restricted adaptor domain end (`cxx:1129-1195`).
        if let Some(bs) = self.basis.osculating_bspline() {
            if let Some(r) = try_local_d1(
                &bs,
                u,
                v,
                self.u1,
                self.u2,
                self.v1,
                self.v2,
                ADAPTOR_TOL,
                ADAPTOR_TOL,
            ) {
                return r;
            }
        }
        self.basis.d1(u, v)
    }

    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        // `GeomAdaptor_Surface::EvalD2` BSpline LocalD2 (`cxx:1273-1292`).
        if let Some(bs) = self.basis.osculating_bspline() {
            if let Some(r) = try_local_d2(
                &bs,
                u,
                v,
                self.u1,
                self.u2,
                self.v1,
                self.v2,
                ADAPTOR_TOL,
                ADAPTOR_TOL,
            ) {
                return r;
            }
        }
        self.basis.d2(u, v)
    }

    fn u_range(&self) -> (f64, f64) {
        (self.u1, self.u2)
    }

    fn v_range(&self) -> (f64, f64) {
        (self.v1, self.v2)
    }

    fn is_u_periodic(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:509-523`.
        if !self.basis.is_u_periodic() {
            return false;
        }
        if !self.u_trimmed {
            return true;
        }
        closed_when_trimmed(true, self.basis.u_period(), self.u1, self.u2)
    }

    fn is_v_periodic(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:534-548`.
        if !self.basis.is_v_periodic() {
            return false;
        }
        if !self.v_trimmed {
            return true;
        }
        closed_when_trimmed(true, self.basis.v_period(), self.v1, self.v2)
    }

    fn is_u_closed(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:559-576`.
        if !self.u_trimmed {
            return self.basis.is_u_closed();
        }
        closed_when_trimmed(
            self.basis.is_u_periodic(),
            self.basis.u_period(),
            self.u1,
            self.u2,
        )
    }

    fn is_v_closed(&self) -> bool {
        // `Geom_RectangularTrimmedSurface.cxx:580-597`.
        if !self.v_trimmed {
            return self.basis.is_v_closed();
        }
        closed_when_trimmed(
            self.basis.is_v_periodic(),
            self.basis.v_period(),
            self.v1,
            self.v2,
        )
    }

    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }

    fn transform(&mut self, t: &GpTrsf) {
        let mut b = self.basis.clone_dyn();
        b.transform(t);
        self.basis = Arc::from(b);
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        Box::new(self.clone())
    }

    fn rectangular_trimmed_basis(&self) -> Option<Arc<dyn Surface>> {
        Some(self.basis.clone())
    }

    // `GeomAdaptor_Surface::load` unwraps a `Geom_RectangularTrimmedSurface` to
    // its basis surface and keeps only the parameter range
    // (`GeomAdaptor_Surface.cxx:423-425`), so the `GetType()`-equivalent queries
    // must report the basis type.
    fn gp_pln(&self) -> Option<GpPln> {
        self.basis.gp_pln()
    }
    fn gp_sphere(&self) -> Option<GpSphere> {
        self.basis.gp_sphere()
    }
    fn gp_cylinder(&self) -> Option<GpCylinder> {
        self.basis.gp_cylinder()
    }
    fn gp_cone(&self) -> Option<GpCone> {
        self.basis.gp_cone()
    }
    fn gp_torus(&self) -> Option<GpTorus> {
        self.basis.gp_torus()
    }

    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        // `Geom_RectangularTrimmedSurface.cxx:444-458`.
        let c = self.basis.u_iso_curve(u)?;
        if self.v_trimmed {
            Some(Arc::new(GeomTrimmedCurve::new(c, self.v1, self.v2)))
        } else {
            Some(c)
        }
    }

    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        // `Geom_RectangularTrimmedSurface.cxx:463-477`.
        let c = self.basis.v_iso_curve(v)?;
        if self.u_trimmed {
            Some(Arc::new(GeomTrimmedCurve::new(c, self.u1, self.u2)))
        } else {
            Some(c)
        }
    }

    fn u_period(&self) -> f64 {
        self.basis.u_period()
    }

    fn v_period(&self) -> f64 {
        self.basis.v_period()
    }

    fn nb_u_poles(&self) -> i32 {
        self.basis.nb_u_poles()
    }
    fn nb_v_poles(&self) -> i32 {
        self.basis.nb_v_poles()
    }
    fn nb_u_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_u_intervals(continuity)
    }
    fn nb_v_intervals(&self, continuity: u8) -> i32 {
        self.basis.nb_v_intervals(continuity)
    }
    fn u_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.u_intervals(continuity)
    }
    fn v_intervals(&self, continuity: u8) -> Vec<f64> {
        self.basis.v_intervals(continuity)
    }
    fn u_degree(&self) -> i32 {
        self.basis.u_degree()
    }
    fn v_degree(&self) -> i32 {
        self.basis.v_degree()
    }
    fn revolution_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.basis.revolution_basis_curve()
    }
    fn extrusion_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.basis.extrusion_basis_curve()
    }

    fn uv_resolution(&self, r3d: f64) -> Option<(f64, f64)> {
        self.basis.uv_resolution(r3d)
    }
}
