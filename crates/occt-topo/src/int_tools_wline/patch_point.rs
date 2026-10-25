//! `IntPatch_Point` subset used by the GeomInt line tools, plus the arc
//! identity the point records.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_Point.hxx:33-263`,
//! `.lxx:16-303`, `.cxx:26-261`.
//!
//! Translation notes:
//! - `PatchPoint` is the reduced `IntPatch_Point` the port's line structures
//!   store; the full class with its owned `Adaptor2d_Curve2d` handle is
//!   `crate::intpatch_point::IntPatchPoint`.
//! - `arcS1` / `arcS2` (`hxx:229`, `:237`) are `occ::handle<Adaptor2d_Curve2d>`
//!   and every use compares two of them with `==` / `!=`, which is handle
//!   *pointer* identity (`Standard_Handle.hxx`). The port therefore stores the
//!   address of the `Arc<Curve2dAdaptor>` as an identity token instead of an
//!   owning handle. A null handle maps to `None`, and `None == None` keeps the
//!   "two null arcs are the same" behaviour of `IntPatch_Point::ArcOnS1()`
//!   returning a null handle.
//! - `prm1` / `prm2` (`hxx:232`, `:240`) are kept as `param_on_arc1` /
//!   `param_on_arc2` (`ParameterOnArc1/2`, `lxx:155-161`, `:206-212`).
//! - `vS1` / `vS2` (`hxx:220`, `:231`) are the `Adaptor3d_HVertex` handles
//!   that `SetVertex` stores and that `TopOpeBRep_FacesIntersector` compares
//!   with `IsSame`. They are owned here (`vertex_s1` / `vertex_s2`), which is
//!   why `PatchPoint` is `Clone` and not `Copy`.

use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::int_tools_wline::PntOn2S;

/// `IntSurf_TypeTrans`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TransType {
    In,
    Out,
    Touch,
    #[default]
    Undecided,
}

/// Identity token of an `occ::handle<Adaptor2d_Curve2d>`: the address of the
/// adaptation, or `None` for a null handle.
pub type ArcId = Option<usize>;

/// [`ArcId`] of a non-null arc handle.
pub fn arc_id(a: &Arc<Curve2dAdaptor>) -> ArcId {
    Some(Arc::as_ptr(a) as usize)
}

/// `IntPatch_Point` subset used by the GeomInt LineConstructor / LineTool.
#[derive(Debug, Clone)]
pub struct PatchPoint {
    /// 3D point.
    pub p: GpPnt,
    /// Parameter on the intersection line (`ParameterOnLine`).
    pub param_on_line: f64,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    pub on_dom_s1: bool,
    pub on_dom_s2: bool,
    /// `vtxonS1` / `vtxonS2` (`IntPatch_Point.hxx:220`, `:229`): the point is a
    /// vertex of the original shape, not only a point on a restriction.
    pub vertex_on_s1: bool,
    pub vertex_on_s2: bool,
    /// `vS1` / `vS2` (`IntPatch_Point.hxx:220`, `:231`): the vertex the point
    /// stands on, stored by `SetVertex` (`cxx:38-55`).
    pub vertex_s1: Option<HVertex>,
    pub vertex_s2: Option<HVertex>,
    pub trans1: TransType,
    pub trans2: TransType,
    pub is_multiple: bool,
    /// `tgt` (`IntPatch_Point.hxx:212`), written by `SetValue` (`cxx:31`).
    pub tangent: bool,
    pub tolerance: f64,
    /// `arcS1` (`IntPatch_Point.hxx:229`), see [`ArcId`].
    pub arc_on_s1: ArcId,
    /// `arcS2` (`IntPatch_Point.hxx:237`).
    pub arc_on_s2: ArcId,
    /// `prm1` (`IntPatch_Point.hxx:232`).
    pub param_on_arc1: f64,
    /// `prm2` (`IntPatch_Point.hxx:240`).
    pub param_on_arc2: f64,
}

impl PatchPoint {
    /// Point with UV on both surfaces and no domain flags.
    pub fn new(p: GpPnt, param_on_line: f64, u1: f64, v1: f64, u2: f64, v2: f64) -> Self {
        Self {
            p,
            param_on_line,
            u1,
            v1,
            u2,
            v2,
            on_dom_s1: false,
            on_dom_s2: false,
            vertex_on_s1: false,
            vertex_on_s2: false,
            vertex_s1: None,
            vertex_s2: None,
            trans1: TransType::Undecided,
            trans2: TransType::Undecided,
            is_multiple: false,
            tangent: false,
            tolerance: 0.0,
            arc_on_s1: None,
            arc_on_s2: None,
            param_on_arc1: 0.0,
            param_on_arc2: 0.0,
        }
    }

    /// `ParameterOnLine`.
    pub fn parameter_on_line(&self) -> f64 {
        self.param_on_line
    }

    /// `SetParameter(Para)` (`IntPatch_Point.lxx:63-66`).
    pub fn set_parameter(&mut self, para: f64) {
        self.param_on_line = para;
    }

    /// `ParametersOnS1`.
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        (self.u1, self.v1)
    }

    /// `ParametersOnS2`.
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        (self.u2, self.v2)
    }

    /// `Parameters(U1, V1, U2, V2)` (`IntPatch_Point.lxx:219-225`).
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.v1, self.u2, self.v2)
    }

    /// `SetParameters(U1, V1, U2, V2)` (`IntPatch_Point.lxx:50-61`).
    pub fn set_parameters(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.u1 = u1;
        self.v1 = v1;
        self.u2 = u2;
        self.v2 = v2;
    }

    /// `PntOn2S()` (`IntPatch_Point.lxx:214-217`).
    pub fn pnt_on_2s(&self) -> PntOn2S {
        PntOn2S {
            p: self.p,
            u1: self.u1,
            v1: self.v1,
            u2: self.u2,
            v2: self.v2,
        }
    }

    /// `SetValue(thePOn2S)` (`IntPatch_Point.lxx:40-43`).
    pub fn set_pnt_on_2s(&mut self, p: &PntOn2S) {
        self.p = p.p;
        self.u1 = p.u1;
        self.v1 = p.v1;
        self.u2 = p.u2;
        self.v2 = p.v2;
    }

    /// `SetTolerance(Tol)` (`IntPatch_Point.lxx:45-48`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// `IsOnDomS1()` (`IntPatch_Point.lxx:112-115`).
    pub fn is_on_dom_s1(&self) -> bool {
        self.on_dom_s1
    }

    /// `IsOnDomS2()` (`IntPatch_Point.lxx:163-166`).
    pub fn is_on_dom_s2(&self) -> bool {
        self.on_dom_s2
    }

    /// `SetArc(OnFirst, A, Param, TLine, TArc)` (`IntPatch_Point.cxx:57-85`):
    /// the arc identity, its parameter and the domain flag of one side.
    ///
    /// `t_line` / `t_arc` collapse onto the single `trans1` / `trans2` slot
    /// (`IntPatchPoint::from_patch_point` maps both OCCT transitions onto it).
    pub fn set_arc(&mut self, on_first: bool, arc: ArcId, param: f64, t_line: TransType) {
        if on_first {
            self.on_dom_s1 = true;
            self.arc_on_s1 = arc;
            self.param_on_arc1 = param;
            self.trans1 = t_line;
        } else {
            self.on_dom_s2 = true;
            self.arc_on_s2 = arc;
            self.param_on_arc2 = param;
            self.trans2 = t_line;
        }
    }

    /// `ArcOnS1()` / `ArcOnS2()` identity token (`lxx:131-137`, `:182-188`).
    pub fn arc_id(&self, on_first: bool) -> ArcId {
        if on_first {
            self.arc_on_s1
        } else {
            self.arc_on_s2
        }
    }

    /// `ParameterOnArc1()` / `ParameterOnArc2()` (`lxx:155-161`, `:206-212`).
    pub fn param_on_arc(&self, on_first: bool) -> f64 {
        if on_first {
            self.param_on_arc1
        } else {
            self.param_on_arc2
        }
    }

    /// `SetVertex(OnFirst, V)` (`IntPatch_Point.cxx:38-55`): the domain and
    /// vertex flags only, the arc slot is untouched.
    pub fn set_vertex_flag(&mut self, on_first: bool) {
        if on_first {
            self.on_dom_s1 = true;
            self.vertex_on_s1 = true;
        } else {
            self.on_dom_s2 = true;
            self.vertex_on_s2 = true;
        }
    }

    /// `SetValue(P, Tol, Tangent)` (`IntPatch_Point.cxx:25-35`): the point, its
    /// tolerance and the tangency flag, with every domain / vertex flag cleared.
    pub fn set_value(&mut self, p: GpPnt, tol: f64, tangent: bool) {
        self.on_dom_s1 = false;
        self.on_dom_s2 = false;
        self.vertex_on_s1 = false;
        self.vertex_on_s2 = false;
        self.vertex_s1 = None;
        self.vertex_s2 = None;
        self.is_multiple = false;
        self.tangent = tangent;
        self.p = p;
        self.tolerance = tol;
    }

    /// `SetVertex(OnFirst, V)` (`IntPatch_Point.cxx:38-55`): the domain flag,
    /// the vertex flag and the stored `Adaptor3d_HVertex` handle.
    pub fn set_vertex(&mut self, on_first: bool, v: HVertex) {
        if on_first {
            self.on_dom_s1 = true;
            self.vertex_on_s1 = true;
            self.vertex_s1 = Some(v);
        } else {
            self.on_dom_s2 = true;
            self.vertex_on_s2 = true;
            self.vertex_s2 = Some(v);
        }
    }

    /// `VertexOnS1()` / `VertexOnS2()` (`IntPatch_Point.lxx:145-152`, `:196-203`).
    pub fn vertex_handle(&self, on_first: bool) -> Option<&HVertex> {
        if on_first {
            self.vertex_s1.as_ref()
        } else {
            self.vertex_s2.as_ref()
        }
    }

    /// `ReverseTransition()` (`IntPatch_Point.cxx:87-137`).
    ///
    /// Each side rebuilds its transition from a fresh `IntSurf_Transition` and
    /// only the `In` / `Out` arms swap the type; the `default` arm leaves the
    /// fresh UNDECIDED transition in place, so a `Touch` transition is lost,
    /// exactly as in OCCT. `PatchPoint` collapses `traline` / `tra` onto the
    /// single `trans1` / `trans2` slot (`IntPatchPoint::from_patch_point`), so
    /// the two switches of the source become one write.
    pub fn reverse_transition(&mut self) {
        if self.on_dom_s1 {
            self.trans1 = match self.trans1 {
                TransType::In => TransType::Out,
                TransType::Out => TransType::In,
                _ => TransType::Undecided,
            };
        }
        if self.on_dom_s2 {
            self.trans2 = match self.trans2 {
                TransType::In => TransType::Out,
                TransType::Out => TransType::In,
                _ => TransType::Undecided,
            };
        }
    }

    /// `IsVertexOnS1()` (`IntPatch_Point.lxx:139-142`).
    pub fn is_vertex_on_s1(&self) -> bool {
        self.vertex_on_s1
    }    /// `IsVertexOnS2()` (`IntPatch_Point.lxx:190-193`).
    pub fn is_vertex_on_s2(&self) -> bool {
        self.vertex_on_s2
    }
}
