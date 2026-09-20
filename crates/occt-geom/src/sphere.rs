use std::sync::Arc;
use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpSphere, GpPnt, GpVec, GpTrsf};
use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::surface::Surface;
use crate::trimmed::GeomTrimmedCurveBasis;

#[derive(Debug, Clone)]
pub struct GeomSphere { pos: GpSphere }

impl GeomSphere {
    pub fn new(s: GpSphere) -> Self { Self { pos: s } }
    /// `Geom_SphericalSurface::Sphere`.
    pub fn sphere(&self) -> &GpSphere { &self.pos }
}

impl Surface for GeomSphere {
    fn d0(&self, u: f64, v: f64) -> GpPnt { slib::sphere_value(&self.pos, u, v) }
    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        surface_eval::sphere_d1(&self.pos, u, v)
    }
    /// `Geom_SphericalSurface::D2` → `ElSLib::SphereD2` (`ElSLib.cxx:975-1037`).
    /// Replaces the trait's central-difference default (audit A15).
    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        surface_eval::sphere_d2(&self.pos, u, v)
    }
    fn u_range(&self) -> (f64, f64) { (0.0, 2.0*std::f64::consts::PI) }
    fn v_range(&self) -> (f64, f64) { (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2) }
    fn is_u_periodic(&self) -> bool { true }
    fn gp_sphere(&self) -> Option<GpSphere> { Some(self.pos) }
    fn continuity(&self) -> u8 { 3 } // T-64: OCCT=GeomAbs_CN(6), blocked by consumer
    fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); }
    fn clone_dyn(&self) -> Box<dyn Surface> { Box::new(self.clone()) }

    /// `Geom_SphericalSurface::UIso` (`Geom_SphericalSurface.cxx:292-297`):
    /// the meridian `Geom_Circle(ElSLib::SphereUIso(pos, radius, U))` trimmed
    /// to `[-PI/2, PI/2]` (OCCT wraps it in a `Geom_TrimmedCurve`).
    /// `ElSLib::SphereUIso` (`ElSLib.cxx:1738-1747`) builds the circle on the
    /// axis `N = cx x dz`, `XDir = cx`, so `ElCLib::CircleValue`'s parameter is
    /// the sphere's `V`; `Geom_TrimmedCurve::EvalD0` (`Geom_TrimmedCurve.cxx:212`)
    /// forwards it unchanged and `FirstParameter`/`LastParameter`
    /// (`cxx:255-267`) are the trimming parameters, so the iso curve's
    /// parameter is the sphere's `V` on `[-PI/2, PI/2]`.
    fn u_iso_curve(&self, u: f64) -> Option<Arc<dyn Curve>> {
        let circ = GeomCircle::new(slib::sphere_u_iso(&self.pos.pos, self.pos.radius, u));
        Some(Arc::new(GeomTrimmedCurveBasis::new(
            Arc::new(circ),
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::FRAC_PI_2,
        )))
    }

    /// `Geom_SphericalSurface::VIso` (`Geom_SphericalSurface.cxx:301-305`):
    /// `Geom_Circle(ElSLib::SphereVIso(pos, radius, V))`, the parallel at
    /// latitude `V`. `ElSLib::SphereVIso` (`ElSLib.cxx:1815-1832`) flips the
    /// circle axis when the radius `radius*cos(V)` is negative (V outside
    /// `[-PI/2, PI/2]`).
    fn v_iso_curve(&self, v: f64) -> Option<Arc<dyn Curve>> {
        Some(Arc::new(GeomCircle::new(slib::sphere_v_iso(
            &self.pos.pos,
            self.pos.radius,
            v,
        ))))
    }
}
