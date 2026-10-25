//! `BRepBlend_Extremity` — the first or last point of a blending line, on a
//! surface or on a curve.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_Extremity.hxx:42-155`,
//! `BRepBlend_Extremity.lxx:16-87`, `BRepBlend_Extremity.cxx:26-142`
//! (OCCT 8.0.0).
//!
//! The four OCCT constructors are the four parameter shapes OCCT overloads on
//! (surface, surface+vertex, curve); each has a matching `SetValue` overload
//! with the same written members, so the port exposes one pair of methods per
//! shape, like [`crate::blend::BlendPoint`].
//!
//! Divergences that are not heuristics:
//! - `Standard_DomainError` from `Tangent()` / `Vertex()` (`lxx:31-35`,
//!   `:60-64`) becomes `Err`,
//! - `occ::handle<Adaptor3d_HVertex>` becomes `Option<HVertex>`: the null case
//!   is observable through [`Extremity::vertex`] instead of raising,
//! - `occ::handle<Adaptor2d_Curve2d>` becomes `Arc<dyn Curve2d>`, and
//!   `NCollection_Sequence<BRepBlend_PointOnRst>` becomes `Vec<PointOnRst>`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom2d::curve::Curve2d;

use crate::brep_top_adaptor_hvertex::HVertex;
use crate::intsurf_transition::IntSurfTransition;

use super::point_on_rst::PointOnRst;

/// `class BRepBlend_Extremity` (`hxx:42-155`).
#[derive(Debug, Clone)]
pub struct Extremity {
    /// `vtx` (`hxx:143`)
    vtx: Option<HVertex>,
    /// `seqpt` (`hxx:144`)
    seqpt: Vec<PointOnRst>,
    /// `pt` (`hxx:145`)
    pt: GpPnt,
    /// `tang` (`hxx:146`)
    tang: GpVec,
    /// `param` (`hxx:147`)
    param: f64,
    /// `u` (`hxx:148`)
    u: f64,
    /// `v` (`hxx:149`)
    v: f64,
    /// `tol` (`hxx:150`)
    tol: f64,
    /// `isvtx` (`hxx:151`)
    isvtx: bool,
    /// `hastang` (`hxx:152`)
    hastang: bool,
}

impl Default for Extremity {
    fn default() -> Self {
        Self::new()
    }
}

impl Extremity {
    /// `BRepBlend_Extremity()` (`cxx:26-36`).
    pub fn new() -> Self {
        Self {
            vtx: None,
            seqpt: Vec::new(),
            pt: GpPnt::zero(),
            tang: GpVec::zero(),
            param: 0.0,
            u: 0.0,
            v: 0.0,
            tol: 0.0,
            isvtx: false,
            hastang: false,
        }
    }

    /// `BRepBlend_Extremity(P, U, V, Param, Tol)` (`cxx:38-51`): an extremity on
    /// a surface.
    pub fn on_surface(p: GpPnt, u: f64, v: f64, param: f64, tol: f64) -> Self {
        let mut e = Self::new();
        e.set_on_surface(p, u, v, param, tol);
        e
    }

    /// `BRepBlend_Extremity(P, U, V, Param, Tol, Vtx)` (`cxx:53-71`): an
    /// extremity on a surface that matches the vertex `vtx`.
    pub fn on_surface_at_vertex(
        p: GpPnt,
        u: f64,
        v: f64,
        param: f64,
        tol: f64,
        vtx: HVertex,
    ) -> Self {
        let mut e = Self::new();
        e.set_on_surface_at_vertex(p, u, v, param, tol, vtx);
        e
    }

    /// `BRepBlend_Extremity(P, W, Param, Tol)` (`cxx:73-86`): an extremity on a
    /// curve. `W` lands in `u` and `v` is set to `0.`.
    pub fn on_curve(p: GpPnt, w: f64, param: f64, tol: f64) -> Self {
        let mut e = Self::new();
        e.set_on_curve(p, w, param, tol);
        e
    }

    /// `SetValue(P, U, V, Param, Tol)` (`cxx:88-101`).
    pub fn set_on_surface(&mut self, p: GpPnt, u: f64, v: f64, param: f64, tol: f64) {
        self.pt = p;
        self.u = u;
        self.v = v;
        self.param = param;
        self.tol = tol;
        self.isvtx = false;
        self.seqpt.clear();
    }

    /// `SetValue(P, U, V, Param, Tol, Vtx)` (`cxx:103-121`).
    pub fn set_on_surface_at_vertex(
        &mut self,
        p: GpPnt,
        u: f64,
        v: f64,
        param: f64,
        tol: f64,
        vtx: HVertex,
    ) {
        self.pt = p;
        self.u = u;
        self.v = v;
        self.param = param;
        self.tol = tol;
        self.isvtx = true;
        self.vtx = Some(vtx);
        self.seqpt.clear();
    }

    /// `SetValue(P, W, Param, Tol)` (`cxx:123-135`): writes `u` and leaves `v`
    /// untouched, as in OCCT (the `v` of the previous state survives).
    pub fn set_on_curve(&mut self, p: GpPnt, w: f64, param: f64, tol: f64) {
        self.pt = p;
        self.u = w;
        self.param = param;
        self.tol = tol;
        self.isvtx = false;
        self.seqpt.clear();
    }

    /// `Value()` (`lxx:18-21`).
    pub fn value(&self) -> &GpPnt {
        &self.pt
    }

    /// `SetTangent(Tangent)` (`lxx:23-27`).
    pub fn set_tangent(&mut self, tangent: GpVec) {
        self.hastang = true;
        self.tang = tangent;
    }

    /// `HasTangent()` (`lxx:29-32`).
    pub fn has_tangent(&self) -> bool {
        self.hastang
    }

    /// `Tangent()` (`lxx:34-41`): `Standard_DomainError` when no tangent was
    /// stored.
    pub fn tangent(&self) -> Result<&GpVec, String> {
        if !self.hastang {
            return Err("BRepBlend_Extremity::Tangent: Standard_DomainError".to_string());
        }
        Ok(&self.tang)
    }

    /// `Tolerance()` (`lxx:47-50`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `SetVertex(V)` (`cxx:137-141`).
    pub fn set_vertex(&mut self, v: HVertex) {
        self.isvtx = true;
        self.vtx = Some(v);
    }

    /// `AddArc(A, Param, TLine, TArc)` (`cxx:143-150`).
    pub fn add_arc(
        &mut self,
        arc: Arc<dyn Curve2d>,
        param: f64,
        tline: IntSurfTransition,
        tarc: IntSurfTransition,
    ) {
        self.seqpt
            .push(PointOnRst::with_arc(arc, param, tline, tarc));
    }

    /// `Parameters(U, V)` (`lxx:43-46`).
    pub fn parameters(&self) -> (f64, f64) {
        (self.u, self.v)
    }

    /// `IsVertex()` (`lxx:52-55`).
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// `Vertex()` (`lxx:57-65`): `Standard_DomainError` outside a vertex
    /// extremity; the null handle OCCT would return is `Ok(None)` here.
    pub fn vertex(&self) -> Result<Option<&HVertex>, String> {
        if !self.isvtx {
            return Err("BRepBlend_Extremity::Vertex: Standard_DomainError".to_string());
        }
        Ok(self.vtx.as_ref())
    }

    /// `NbPointOnRst()` (`lxx:67-70`).
    pub fn nb_point_on_rst(&self) -> usize {
        self.seqpt.len()
    }

    /// `PointOnRst(Index)` (`lxx:72-75`), `Index` `1`-based. OCCT reads the
    /// sequence out of range; the port answers `None`.
    pub fn point_on_rst(&self, index: usize) -> Option<&PointOnRst> {
        self.seqpt.get(index.checked_sub(1)?)
    }

    /// `Parameter()` (`lxx:77-80`): the curve parameter `W` (or the surface `u`).
    pub fn parameter(&self) -> f64 {
        self.u
    }

    /// `ParameterOnGuide()` (`lxx:82-86`).
    pub fn parameter_on_guide(&self) -> f64 {
        self.param
    }
}
