//! `Contap_Point` - one point of a contour line, with the geometric and the
//! logical information the contour algorithms attach to it.
//!
//! Sources: `Contap_Point.hxx:31-132`, `.cxx:22-44`, `.lxx:19-153`.
//!
//! The five accessors that OCCT guards with `Standard_DomainError`
//! (`lxx:103-107`, `:113-117`, `:122-126`, `:134-139`, `:147-151`) answer
//! `Err` here.

use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::intsurf_transition::IntSurfTransition;

/// `Contap_Point` (`hxx:31-132`): a vertex of a contour line.
///
/// No `Debug`: the watched arc is a `Curve2dAdaptor`, which has none.
#[derive(Clone)]
pub struct ContapPoint {
    /// `pt` (`hxx:119`).
    pt: GpPnt,
    /// `uparam` (`hxx:120`).
    uparam: f64,
    /// `vparam` (`hxx:121`).
    vparam: f64,
    /// `paraline` (`hxx:122`).
    paraline: f64,
    /// `onarc` (`hxx:123`).
    onarc: bool,
    /// `arc` (`hxx:124`): the restriction arc the point lies on.
    arc: Option<Arc<Curve2dAdaptor>>,
    /// `traline` (`hxx:125`).
    traline: IntSurfTransition,
    /// `traarc` (`hxx:126`).
    traarc: IntSurfTransition,
    /// `prmarc` (`hxx:127`).
    prmarc: f64,
    /// `isvtx` (`hxx:128`).
    isvtx: bool,
    /// `vtx` (`hxx:129`).
    vtx: Option<HVertex>,
    /// `ismult` (`hxx:130`).
    ismult: bool,
    /// `myInternal` (`hxx:131`).
    internal: bool,
}

impl Default for ContapPoint {
    /// `Contap_Point()` (`cxx:22-30`): `uparam`/`vparam`/`paraline`/`prmarc`
    /// zero, the three flags false, both transitions default (UNDECIDED).
    /// `pt` is a default `gp_Pnt`, i.e. the origin.
    fn default() -> Self {
        Self::new()
    }
}

impl ContapPoint {
    /// `Contap_Point()` (`cxx:22-30`).
    pub fn new() -> Self {
        Self {
            pt: GpPnt::default(),
            uparam: 0.0,
            vparam: 0.0,
            paraline: 0.0,
            onarc: false,
            arc: None,
            traline: IntSurfTransition::new(),
            traarc: IntSurfTransition::new(),
            prmarc: 0.0,
            isvtx: false,
            vtx: None,
            ismult: false,
            internal: false,
        }
    }

    /// `Contap_Point(Pt, U, V)` (`cxx:33-44`).
    pub fn with_uv(pt: GpPnt, u: f64, v: f64) -> Self {
        Self {
            pt,
            uparam: u,
            vparam: v,
            ..Self::new()
        }
    }

    /// `SetValue(Pt, U, V)` (`lxx:19-28`): resets the three logical flags.
    pub fn set_value(&mut self, pt: GpPnt, u: f64, v: f64) {
        self.pt = pt;
        self.uparam = u;
        self.vparam = v;
        self.onarc = false;
        self.isvtx = false;
        self.ismult = false;
        self.internal = false;
    }

    /// `SetParameter(Para)` (`lxx:30-33`).
    pub fn set_parameter(&mut self, para: f64) {
        self.paraline = para;
    }

    /// `SetVertex(V)` (`lxx:35-41`).
    pub fn set_vertex(&mut self, v: HVertex) {
        self.isvtx = true;
        self.vtx = Some(v);
    }

    /// `SetArc(A, Param, TLine, TArc)` (`lxx:43-56`).
    pub fn set_arc(
        &mut self,
        a: Arc<Curve2dAdaptor>,
        param: f64,
        t_line: IntSurfTransition,
        t_arc: IntSurfTransition,
    ) {
        self.onarc = true;
        self.arc = Some(a);
        self.prmarc = param;
        self.traline = t_line;
        self.traarc = t_arc;
    }

    /// `SetMultiple()` (`lxx:58-61`).
    pub fn set_multiple(&mut self) {
        self.ismult = true;
    }

    /// `SetInternal()` (`lxx:63-66`).
    pub fn set_internal(&mut self) {
        self.internal = true;
    }

    /// `IsMultiple()` (`lxx:68-71`).
    pub fn is_multiple(&self) -> bool {
        self.ismult
    }

    /// `IsInternal()` (`lxx:73-76`).
    pub fn is_internal(&self) -> bool {
        self.internal
    }

    /// `Value()` (`lxx:78-81`).
    pub fn value(&self) -> &GpPnt {
        &self.pt
    }

    /// `ParameterOnLine()` (`lxx:83-86`).
    pub fn parameter_on_line(&self) -> f64 {
        self.paraline
    }

    /// `Parameters(U1, V1)` (`lxx:88-93`).
    pub fn parameters(&self) -> (f64, f64) {
        (self.uparam, self.vparam)
    }

    /// `IsOnArc()` (`lxx:95-98`).
    pub fn is_on_arc(&self) -> bool {
        self.onarc
    }

    /// `Arc()` (`lxx:100-108`): `Standard_DomainError` when `IsOnArc()` is
    /// false.
    pub fn arc(&self) -> Result<&Arc<Curve2dAdaptor>, &'static str> {
        self.arc
            .as_ref()
            .ok_or("Contap_Point::Arc: Standard_DomainError")
    }

    /// `ParameterOnArc()` (`lxx:128-137`).
    pub fn parameter_on_arc(&self) -> Result<f64, &'static str> {
        if !self.onarc {
            return Err("Contap_Point::ParameterOnArc: Standard_DomainError");
        }
        Ok(self.prmarc)
    }

    /// `TransitionOnLine()` (`lxx:110-118`).
    pub fn transition_on_line(&self) -> Result<&IntSurfTransition, &'static str> {
        if !self.onarc {
            return Err("Contap_Point::TransitionOnLine: Standard_DomainError");
        }
        Ok(&self.traline)
    }

    /// `TransitionOnArc()` (`lxx:120-128`).
    pub fn transition_on_arc(&self) -> Result<&IntSurfTransition, &'static str> {
        if !self.onarc {
            return Err("Contap_Point::TransitionOnArc: Standard_DomainError");
        }
        Ok(&self.traarc)
    }

    /// `IsVertex()` (`lxx:142-145`).
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// `Vertex()` (`lxx:147-153`): `Standard_DomainError` when `IsVertex()` is
    /// false.
    pub fn vertex(&self) -> Result<&HVertex, &'static str> {
        self.vtx
            .as_ref()
            .ok_or("Contap_Point::Vertex: Standard_DomainError")
    }
}
