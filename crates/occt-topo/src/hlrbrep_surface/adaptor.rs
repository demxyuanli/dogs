//! `HLRBRep_Surface` seen through the `IntCurveSurface` `SurfaceTool` protocol.
//!
//! `HLRBRep_Surface` is a `BRepAdaptor_Surface` plus the cached `GetType()`
//! (`HLRBRep_Surface.hxx:185-193`). The `IntCurveSurface_Inter` template that
//! `HLRBRep_InterCSurf` instantiates (`HLRBRep_InterCSurf.cxx:47-48`) only ever
//! calls it through `HLRBRep_SurfaceTool`, whose members are these adaptor
//! queries (`HLRBRep_SurfaceTool.hxx` + `HLRBRep_Surface.lxx`). This view is the
//! same surface as a [`occt_geom::Surface`], which is the interface the ported
//! `intcurvesurface` engine takes, so the engine can be instantiated on the HLR
//! surface exactly as OCCT instantiates its template on
//! `HLRBRep_Surface*`:
//!
//! * `u_range` / `v_range` are the **face** window
//!   (`HLRBRep_SurfaceTool::First/LastUParameter`, `lxx:37-61`), not the
//!   underlying `Geom_Surface` bounds. Every `ComputeAppendPoint` acceptance
//!   test (`IntCurveSurface_InterUtils.pxx:1117-1176`) uses this window, so it
//!   must be the restricted one.
//! * the quadric accessors answer the `HLRBRep_SurfaceTool::Plane/Cylinder/
//!   Cone/Sphere/Torus` objects (`lxx:198-222`, `HLRBRep_Surface.cxx:319-335`),
//!   i.e. `HLRBRep_Surface::Plane()`'s tangent plane for a Bezier patch.
//! * `u_period` / `v_period` are the adaptor's periods
//!   (`GeomAdaptor_Surface::UPeriod`), which for a periodic face differ from the
//!   face's `u_range` width.
//!
//! `transform` and `clone_dyn` are the two `Geom_Surface` operations with no
//! `HLRBRep_SurfaceTool` counterpart (`HLRBRep_Surface.hxx:42-198` declares
//! neither). The intersection engine never calls them, so they panic with a
//! message naming this fact, the same convention as
//! `Curve2d for HlrBRepCurve` (`hlrbrep_curve_2d.rs`).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpTrsf, GpVec};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom::surface_d3::SurfD3;

use super::HlrBRepSurface;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;

/// A borrowed `HLRBRep_Surface` as an [`occt_geom::Surface`].
#[derive(Clone, Copy)]
pub struct HlrBRepSurfaceView<'a> {
    hlr: &'a HlrBRepSurface,
    /// `HLRBRep_Surface::Surface()` (`lxx:23-33`): the adaptor's
    /// `Handle(Geom_Surface)`, the object `HLRBRep_SurfaceTool::D0/...` evaluate
    /// (`HLRBRep_SurfaceTool` forwards to `BRepAdaptor_Surface`, which forwards
    /// to `GeomAdaptor_Surface`, which evaluates this same handle).
    geom: &'a dyn Surface,
}

impl<'a> HlrBRepSurfaceView<'a> {
    /// `None` when the HLR surface has no underlying `Geom_Surface`, the state
    /// in which every `HLRBRep_SurfaceTool` evaluation throws
    /// (`BRepAdaptor_Surface::Surface` on an uninitialized adaptor).
    pub fn new(hlr: &'a HlrBRepSurface) -> Option<Self> {
        let geom = hlr.surface()?.as_ref();
        Some(Self { hlr, geom })
    }

    /// The `HLRBRep_Surface` this view exposes.
    pub fn hlr(&self) -> &'a HlrBRepSurface {
        self.hlr
    }
}

impl Surface for HlrBRepSurfaceView<'_> {
    fn d0(&self, u: f64, v: f64) -> GpPnt {
        self.geom.d0(u, v)
    }

    fn d1(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
        self.geom.d1(u, v)
    }

    fn d2(&self, u: f64, v: f64) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        self.geom.d2(u, v)
    }

    fn d3(&self, u: f64, v: f64) -> Option<SurfD3> {
        self.geom.d3(u, v)
    }

    fn eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> GpVec {
        self.geom.eval_dn(u, v, nu, nv)
    }

    fn try_eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        self.geom.try_eval_dn(u, v, nu, nv)
    }

    /// `HLRBRep_SurfaceTool::FirstUParameter` / `LastUParameter` (`lxx:37-47`).
    fn u_range(&self) -> (f64, f64) {
        (
            self.hlr.first_u_parameter(),
            self.hlr.last_u_parameter(),
        )
    }

    /// `HLRBRep_SurfaceTool::FirstVParameter` / `LastVParameter` (`lxx:51-61`).
    fn v_range(&self) -> (f64, f64) {
        (
            self.hlr.first_v_parameter(),
            self.hlr.last_v_parameter(),
        )
    }

    fn is_u_periodic(&self) -> bool {
        self.hlr.is_u_periodic()
    }

    fn is_v_periodic(&self) -> bool {
        self.hlr.is_v_periodic()
    }

    fn u_period(&self) -> f64 {
        self.hlr.u_period()
    }

    fn v_period(&self) -> f64 {
        self.hlr.v_period()
    }

    fn is_u_closed(&self) -> bool {
        self.hlr.is_u_closed()
    }

    fn is_v_closed(&self) -> bool {
        self.hlr.is_v_closed()
    }

    fn continuity(&self) -> u8 {
        self.geom.continuity()
    }

    fn nb_u_intervals(&self, s: u8) -> i32 {
        self.geom.nb_u_intervals(s)
    }

    fn nb_v_intervals(&self, s: u8) -> i32 {
        self.geom.nb_v_intervals(s)
    }

    fn u_intervals(&self, s: u8) -> Vec<f64> {
        self.geom.u_intervals(s)
    }

    fn v_intervals(&self, s: u8) -> Vec<f64> {
        self.geom.v_intervals(s)
    }

    fn is_bspline_surface(&self) -> bool {
        self.hlr.adaptor_type() == GeomAbsSurfaceType::BSplineSurface
    }

    fn is_bezier_surface(&self) -> bool {
        self.hlr.adaptor_type() == GeomAbsSurfaceType::BezierSurface
    }

    fn nb_u_poles(&self) -> i32 {
        self.hlr.nb_u_poles().unwrap_or(2)
    }

    fn nb_v_poles(&self) -> i32 {
        self.hlr.nb_v_poles().unwrap_or(2)
    }

    fn u_degree(&self) -> i32 {
        self.hlr.u_degree().unwrap_or(1)
    }

    fn v_degree(&self) -> i32 {
        self.hlr.v_degree().unwrap_or(1)
    }

    /// `HLRBRep_SurfaceTool::Plane` (`lxx:190-195` = `HLRBRep_Surface::Plane`,
    /// `cxx:319-335`).
    fn gp_pln(&self) -> Option<GpPln> {
        self.hlr.plane()
    }

    fn gp_cylinder(&self) -> Option<GpCylinder> {
        self.hlr.cylinder()
    }

    fn gp_cone(&self) -> Option<GpCone> {
        self.hlr.cone()
    }

    fn gp_sphere(&self) -> Option<GpSphere> {
        self.hlr.sphere()
    }

    fn gp_torus(&self) -> Option<GpTorus> {
        self.hlr.torus()
    }

    fn is_surface_of_revolution(&self) -> bool {
        self.hlr.adaptor_type() == GeomAbsSurfaceType::SurfaceOfRevolution
    }

    fn is_surface_of_linear_extrusion(&self) -> bool {
        self.hlr.adaptor_type() == GeomAbsSurfaceType::SurfaceOfExtrusion
    }

    fn revolution_axis(&self) -> Option<GpAx1> {
        self.hlr.axis()
    }

    fn revolution_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.geom.revolution_basis_curve()
    }

    fn extrusion_basis_curve(&self) -> Option<Arc<dyn Curve>> {
        self.geom.extrusion_basis_curve()
    }

    fn transform(&mut self, _t: &GpTrsf) {
        unreachable!(
            "HlrBRepSurfaceView::transform: HLRBRep_Surface has no Trsf() \
             (HLRBRep_Surface.hxx:42-198), so HLRBRep_SurfaceTool never transforms \
             the surface it evaluates"
        )
    }

    fn clone_dyn(&self) -> Box<dyn Surface> {
        unreachable!(
            "HlrBRepSurfaceView::clone_dyn: the view borrows an HLRBRep_Surface and \
             HLRBRep_Surface.hxx:42-198 declares no Copy(); the intersection engine \
             never clones its SurfaceTool argument"
        )
    }
}
