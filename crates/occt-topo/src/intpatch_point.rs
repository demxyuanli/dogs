//! `IntPatch_Point` — one point of an intersection line between two surfaces:
//! its 3D value, its parameters on each surface, and its position relative to
//! the restriction arcs of either domain.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_Point.hxx:33-263`,
//! `.lxx:16-303`, `.cxx:26-261`.
//!
//! Translation notes:
//! - The two `occ::handle<Adaptor3d_HVertex> vS1/vS2` members (`hxx:233`, `:241`)
//!   are `HVertex` values: `TopoDS_Vertex` is a value type in OCCT, so an unset
//!   handle and a `HVertex` over a null vertex behave identically (only
//!   `IsVertexOnS*` gates access to them).
//! - The two `occ::handle<Adaptor2d_Curve2d> arcS1/arcS2` (`hxx:235`, `:243`)
//!   become `Option<Arc<Curve2dAdaptor>>`. `TopOpeBRep` only ever hands this
//!   class `BRepAdaptor_Curve2d`-backed arcs (`TopOpeBRep_VPointInter::ArcOnS1`
//!   down-casts to it at `VPointInter.cxx:63-68`), so the abstract
//!   `Adaptor2d_Curve2d` handle is collapsed onto that one class.
//! - Every accessor gated by `onS1` / `onS2` / `vtxonS1` / `vtxonS2` raises
//!   `Standard_DomainError` in OCCT (`lxx:132-136`, `:146-150`, `:153-157`,
//!   `:160-164`, `:167-171`, `:186-190`, `:200-204`, `:207-211`, `:214-218`,
//!   `:221-225`, `:235-239`); the port returns `Err`.
//! - `Dump()` (`cxx:186-261`) writes to `std::cout` only, there is no state in
//!   it; the port emits the same text through `std::fmt`-free `String` building
//!   in [`IntPatchPoint::dump`].
//!
//! `Dump` is kept because `IntPatch_Point::Dump` is a released (non-debug)
//! member in `cxx:186-261`, unlike the `OCCT_DEBUG`-only dumps elsewhere.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;

use crate::brep_top_adaptor_hvertex::HVertex;
use crate::int_tools_wline::PntOn2S;
use crate::intsurf_transition::{IntSurfTransition, IntSurfTypeTrans};

/// `class IntPatch_Point` (`hxx:33-263`).
#[derive(Clone)]
pub struct IntPatchPoint {
    /// `pt` (`hxx:218`): 3D value + the four surface parameters.
    pt: PntOn2S,
    /// `para` (`hxx:219`): parameter on the intersection line.
    para: f64,
    /// `tol` (`hxx:220`).
    tol: f64,
    /// `tgt` (`hxx:221`): a tangency point between the two surfaces.
    tgt: bool,
    /// `mult` (`hxx:224`): the point is on several intersection lines.
    mult: bool,
    /// `onS1` (`hxx:226`).
    on_s1: bool,
    /// `vtxonS1` (`hxx:227`).
    vtx_on_s1: bool,
    /// `vS1` (`hxx:228`).
    v_s1: HVertex,
    /// `arcS1` (`hxx:229`).
    arc_s1: Option<Arc<Curve2dAdaptor>>,
    /// `traline1` (`hxx:230`).
    tra_line1: IntSurfTransition,
    /// `tra1` (`hxx:231`).
    tra1: IntSurfTransition,
    /// `prm1` (`hxx:232`).
    prm1: f64,
    /// `onS2` (`hxx:234`).
    on_s2: bool,
    /// `vtxonS2` (`hxx:235`).
    vtx_on_s2: bool,
    /// `vS2` (`hxx:236`).
    v_s2: HVertex,
    /// `arcS2` (`hxx:237`).
    arc_s2: Option<Arc<Curve2dAdaptor>>,
    /// `traline2` (`hxx:238`).
    tra_line2: IntSurfTransition,
    /// `tra2` (`hxx:239`).
    tra2: IntSurfTransition,
    /// `prm2` (`hxx:240`).
    prm2: f64,
}

impl Default for IntPatchPoint {
    /// `IntPatch_Point()` (`lxx:17-31`): the init list is exactly
    /// `para/tol/tgt/mult/onS1/vtxonS1/prm1/onS2/vtxonS2/prm2`; `pt` and the
    /// four transitions keep their class defaults.
    fn default() -> Self {
        Self {
            pt: PntOn2S {
                p: GpPnt::new(0.0, 0.0, 0.0),
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            },
            para: 0.0,
            tol: 0.0,
            tgt: false,
            mult: false,
            on_s1: false,
            vtx_on_s1: false,
            v_s1: HVertex::new(&crate::shape::Vertex::new()),
            arc_s1: None,
            tra_line1: IntSurfTransition::new(),
            tra1: IntSurfTransition::new(),
            prm1: 0.0,
            on_s2: false,
            vtx_on_s2: false,
            v_s2: HVertex::new(&crate::shape::Vertex::new()),
            arc_s2: None,
            tra_line2: IntSurfTransition::new(),
            tra2: IntSurfTransition::new(),
            prm2: 0.0,
        }
    }
}

impl IntPatchPoint {
    /// `IntPatch_Point()` (`lxx:17-31`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetValue(Pt, Tol, Tangent)` (`cxx:27-38`): resets the domain flags,
    /// then stores the 3D point and the tolerance.
    pub fn set_value_with_tolerance(&mut self, p: &GpPnt, tol: f64, tangent: bool) {
        self.on_s1 = false;
        self.on_s2 = false;
        self.vtx_on_s1 = false;
        self.vtx_on_s2 = false;
        self.mult = false;
        self.tgt = tangent;
        self.pt.set_point(*p);
        self.tol = tol;
    }

    /// `SetValue(Pt)` (`lxx:33-36`): the 3D point only.
    pub fn set_value(&mut self, p: &GpPnt) {
        self.pt.set_point(*p);
    }

    /// `SetValue(thePOn2S)` (`lxx:40-43`).
    pub fn set_value_pnt_on_2s(&mut self, p: &PntOn2S) {
        self.pt = *p;
    }

    /// `SetTolerance(Tol)` (`lxx:45-48`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// `SetParameters(U1, V1, U2, V2)` (`lxx:50-61`).
    pub fn set_parameters(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.pt.set_parameters(u1, v1, u2, v2);
    }

    /// `SetParameter(Para)` (`lxx:63-66`).
    pub fn set_parameter(&mut self, para: f64) {
        self.para = para;
    }

    /// `SetVertex(OnFirst, V)` (`cxx:40-55`).
    pub fn set_vertex(&mut self, on_first: bool, v: &HVertex) {
        if on_first {
            self.on_s1 = true;
            self.vtx_on_s1 = true;
            self.v_s1 = v.clone();
        } else {
            self.on_s2 = true;
            self.vtx_on_s2 = true;
            self.v_s2 = v.clone();
        }
    }

    /// `SetArc(OnFirst, A, Param, TLine, TArc)` (`cxx:57-85`).
    pub fn set_arc(
        &mut self,
        on_first: bool,
        a: Arc<Curve2dAdaptor>,
        param: f64,
        t_line: &IntSurfTransition,
        t_arc: &IntSurfTransition,
    ) {
        if on_first {
            self.on_s1 = true;
            self.arc_s1 = Some(a);
            self.tra_line1 = *t_line;
            self.tra1 = *t_arc;
            self.prm1 = param;
        } else {
            self.on_s2 = true;
            self.arc_s2 = Some(a);
            self.tra_line2 = *t_line;
            self.tra2 = *t_arc;
            self.prm2 = param;
        }
    }

    /// `SetMultiple(IsMult)` (`lxx:68-71`).
    pub fn set_multiple(&mut self, is_mult: bool) {
        self.mult = is_mult;
    }

    /// `Value()` (`lxx:73-76`).
    pub fn value(&self) -> &GpPnt {
        &self.pt.p
    }

    /// `ParameterOnLine()` (`lxx:82-85`).
    pub fn parameter_on_line(&self) -> f64 {
        self.para
    }

    /// `Tolerance()` (`lxx:87-90`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `IsTangencyPoint()` (`lxx:92-95`).
    pub fn is_tangency_point(&self) -> bool {
        self.tgt
    }

    /// `ParametersOnS1(U1, V1)` (`lxx:97-100`).
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        self.pt.parameters_on_s1()
    }

    /// `ParametersOnS2(U2, V2)` (`lxx:102-105`).
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        self.pt.parameters_on_s2()
    }

    /// `IsMultiple()` (`lxx:107-110`).
    pub fn is_multiple(&self) -> bool {
        self.mult
    }

    /// `IsOnDomS1()` (`lxx:112-115`).
    pub fn is_on_dom_s1(&self) -> bool {
        self.on_s1
    }

    /// `IsVertexOnS1()` (`lxx:117-120`).
    pub fn is_vertex_on_s1(&self) -> bool {
        self.vtx_on_s1
    }

    /// `VertexOnS1()` (`lxx:122-129`).
    pub fn vertex_on_s1(&self) -> Result<&HVertex, String> {
        if !self.vtx_on_s1 {
            return Err("IntPatch_Point::VertexOnS1: Standard_DomainError".to_string());
        }
        Ok(&self.v_s1)
    }

    /// `ArcOnS1()` (`lxx:131-137`).
    pub fn arc_on_s1(&self) -> Result<&Arc<Curve2dAdaptor>, String> {
        if !self.on_s1 {
            return Err("IntPatch_Point::ArcOnS1: Standard_DomainError".to_string());
        }
        self.arc_s1
            .as_ref()
            .ok_or_else(|| "IntPatch_Point::ArcOnS1".to_string())
    }

    /// `TransitionLineArc1()` (`lxx:139-145`).
    pub fn transition_line_arc1(&self) -> Result<&IntSurfTransition, String> {
        if !self.on_s1 {
            return Err("IntPatch_Point::TransitionLineArc1: Standard_DomainError".to_string());
        }
        Ok(&self.tra_line1)
    }

    /// `TransitionOnS1()` (`lxx:147-153`).
    pub fn transition_on_s1(&self) -> Result<&IntSurfTransition, String> {
        if !self.on_s1 {
            return Err("IntPatch_Point::TransitionOnS1: Standard_DomainError".to_string());
        }
        Ok(&self.tra1)
    }

    /// `ParameterOnArc1()` (`lxx:155-161`).
    pub fn parameter_on_arc1(&self) -> Result<f64, String> {
        if !self.on_s1 {
            return Err("IntPatch_Point::ParameterOnArc1: Standard_DomainError".to_string());
        }
        Ok(self.prm1)
    }

    /// `IsOnDomS2()` (`lxx:163-166`).
    pub fn is_on_dom_s2(&self) -> bool {
        self.on_s2
    }

    /// `IsVertexOnS2()` (`lxx:168-171`).
    pub fn is_vertex_on_s2(&self) -> bool {
        self.vtx_on_s2
    }

    /// `VertexOnS2()` (`lxx:173-180`).
    pub fn vertex_on_s2(&self) -> Result<&HVertex, String> {
        if !self.vtx_on_s2 {
            return Err("IntPatch_Point::VertexOnS2: Standard_DomainError".to_string());
        }
        Ok(&self.v_s2)
    }

    /// `ArcOnS2()` (`lxx:182-188`).
    pub fn arc_on_s2(&self) -> Result<&Arc<Curve2dAdaptor>, String> {
        if !self.on_s2 {
            return Err("IntPatch_Point::ArcOnS2: Standard_DomainError".to_string());
        }
        self.arc_s2
            .as_ref()
            .ok_or_else(|| "IntPatch_Point::ArcOnS2".to_string())
    }

    /// `TransitionLineArc2()` (`lxx:190-196`).
    pub fn transition_line_arc2(&self) -> Result<&IntSurfTransition, String> {
        if !self.on_s2 {
            return Err("IntPatch_Point::TransitionLineArc2: Standard_DomainError".to_string());
        }
        Ok(&self.tra_line2)
    }

    /// `TransitionOnS2()` (`lxx:198-204`).
    pub fn transition_on_s2(&self) -> Result<&IntSurfTransition, String> {
        if !self.on_s2 {
            return Err("IntPatch_Point::TransitionOnS2: Standard_DomainError".to_string());
        }
        Ok(&self.tra2)
    }

    /// `ParameterOnArc2()` (`lxx:206-212`).
    pub fn parameter_on_arc2(&self) -> Result<f64, String> {
        if !self.on_s2 {
            return Err("IntPatch_Point::ParameterOnArc2: Standard_DomainError".to_string());
        }
        Ok(self.prm2)
    }

    /// `PntOn2S()` (`lxx:214-217`).
    pub fn pnt_on_2s(&self) -> &PntOn2S {
        &self.pt
    }

    /// `Parameters(U1, V1, U2, V2)` (`lxx:219-225`).
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        self.pt.parameters()
    }

    /// `ReverseTransition()` (`cxx:87-150`).
    ///
    /// Each of the four transitions is rebuilt from a fresh
    /// `IntSurf_Transition T` and only the `In` / `Out` arms swap the type; the
    /// `default` arm leaves the fresh (UNDECIDED) transition in place, so a
    /// `Touch` transition is *lost*, exactly as in OCCT.
    pub fn reverse_transition(&mut self) {
        if self.on_s1 {
            let mut t_line = IntSurfTransition::new();
            match self.tra_line1.transition_type() {
                IntSurfTypeTrans::In => {
                    t_line.set_value_with_type(false, IntSurfTypeTrans::Out)
                }
                IntSurfTypeTrans::Out => {
                    t_line.set_value_with_type(false, IntSurfTypeTrans::In)
                }
                _ => {}
            }
            self.tra_line1 = t_line;
            let mut t_arc = IntSurfTransition::new();
            match self.tra1.transition_type() {
                IntSurfTypeTrans::In => t_arc.set_value_with_type(false, IntSurfTypeTrans::Out),
                IntSurfTypeTrans::Out => t_arc.set_value_with_type(false, IntSurfTypeTrans::In),
                _ => {}
            }
            self.tra1 = t_arc;
        }
        if self.on_s2 {
            let mut t_line = IntSurfTransition::new();
            match self.tra_line2.transition_type() {
                IntSurfTypeTrans::In => {
                    t_line.set_value_with_type(false, IntSurfTypeTrans::Out)
                }
                IntSurfTypeTrans::Out => {
                    t_line.set_value_with_type(false, IntSurfTypeTrans::In)
                }
                _ => {}
            }
            self.tra_line2 = t_line;
            let mut t_arc = IntSurfTransition::new();
            match self.tra2.transition_type() {
                IntSurfTypeTrans::In => t_arc.set_value_with_type(false, IntSurfTypeTrans::Out),
                IntSurfTypeTrans::Out => t_arc.set_value_with_type(false, IntSurfTypeTrans::In),
                _ => {}
            }
            self.tra2 = t_arc;
        }
    }

    /// `Dump()` (`cxx:186-261`).
    ///
    /// OCCT prints handles with `%p`; the port prints `null` / `set` instead,
    /// the only difference in the emitted text.
    pub fn dump(&self) -> String {
        let mut s = String::new();
        let (u1, v1, u2, v2) = self.pt.parameters();
        s.push_str("----------- IntPatch_Point : \n");
        s.push_str(&format!(
            "P({:+10.20},{:+10.20},{:+10.20}) UV1({:+10.20},{:+10.20})  UV2({:+10.20},{:+10.20}) (Para:{:+10.20})\n",
            self.pt.p.x(),
            self.pt.p.y(),
            self.pt.p.z(),
            u1,
            v1,
            u2,
            v2,
            self.para
        ));
        if self.on_s1 {
            s.push_str(&format!(
                "*OnS1*  par={:+10.20} arc1={}",
                self.prm1,
                if self.arc_s1.is_some() { "set" } else { "null" }
            ));
        }
        if self.vtx_on_s1 {
            s.push_str(" *Vtx1*  vtx1=set");
        }
        if self.on_s1 || self.vtx_on_s1 {
            s.push('\n');
        }
        if self.on_s2 {
            s.push_str(&format!(
                "*OnS2*  par={:+10.20} arc2={}",
                self.prm2,
                if self.arc_s2.is_some() { "set" } else { "null" }
            ));
        }
        if self.vtx_on_s2 {
            s.push_str(" *Vtx2*  vtx2=set");
        }
        if self.on_s2 || self.vtx_on_s2 {
            s.push('\n');
        }
        if self.on_s1 || self.on_s2 {
            s.push_str(&format!(
                " tgt:{}   mult:{}",
                if self.tgt { 1 } else { 0 },
                if self.mult { 1 } else { 0 }
            ));
            if self.on_s1 {
                s.push_str("\n-traline1  : ");
                s.push_str(&dump_transition(&self.tra_line1));
                s.push_str("  -tra1  : ");
                s.push_str(&dump_transition(&self.tra1));
            }
            if self.on_s2 {
                s.push_str("\n-traline2  : ");
                s.push_str(&dump_transition(&self.tra_line2));
                s.push_str("  -tra2  : ");
                s.push_str(&dump_transition(&self.tra2));
            }
            s.push('\n');
        }
        s.push_str(&format!("\nToler = {}\n", self.tol));
        s
    }
}

/// `static void DumpTransition(const IntSurf_Transition& T)` (`cxx:154-184`).
fn dump_transition(t: &IntSurfTransition) -> String {
    let mut s = String::new();
    let typetrans = t.transition_type();
    match typetrans {
        IntSurfTypeTrans::In => s.push_str(" In  "),
        IntSurfTypeTrans::Out => s.push_str(" Out "),
        IntSurfTypeTrans::Undecided => s.push_str(" Und."),
        IntSurfTypeTrans::Touch => {
            s.push_str(" Touch ");
            match t.situation() {
                Ok(crate::intsurf_transition::IntSurfSituation::Inside) => s.push_str(" Inside "),
                Ok(crate::intsurf_transition::IntSurfSituation::Outside) => s.push_str(" Outside"),
                _ => s.push_str(" Unknown"),
            }
            if t.is_opposite().unwrap_or(false) {
                s.push_str("     Opposite ");
            } else {
                s.push_str(" Non_Opposite ");
            }
        }
    }
    if typetrans != IntSurfTypeTrans::Undecided {
        if t.is_tangent().unwrap_or(false) {
            s.push_str(" Tangent    ");
        } else {
            s.push_str(" Non_Tangent");
        }
    }
    s
}

impl IntPatchPoint {
    /// The full class behind the port's `PatchPoint` subset
    /// (`int_tools_wline::PatchPoint`, "`IntPatch_Point` subset used by GeomInt
    /// LineConstructor / LineTool").
    ///
    /// `PatchPoint` carries the domain flags and the transitions but not the
    /// restriction-arc handle an `IntPatch_Point` is built with
    /// (`IntPatch_Point.cxx:57-85`); the caller that owns the line supplies it
    /// (an `IntPatch_RLine` hands the same arc to every vertex). With `None`
    /// the point keeps OCCT's "on a domain, no arc" state, whose arc-gated
    /// accessors answer `Err` here.
    pub fn from_patch_point(
        p: &crate::int_tools_wline::PatchPoint,
        arc1: Option<Arc<Curve2dAdaptor>>,
        arc2: Option<Arc<Curve2dAdaptor>>,
    ) -> Self {
        use crate::intsurf_transition::IntSurfTypeTrans;

        let mut r = Self::default();
        r.pt = PntOn2S {
            p: p.p,
            u1: p.u1,
            v1: p.v1,
            u2: p.u2,
            v2: p.v2,
        };
        r.para = p.param_on_line;
        r.tol = p.tolerance;
        r.mult = p.is_multiple;
        let trans = |t: crate::int_tools_wline::TransType| {
            let ty = match t {
                crate::int_tools_wline::TransType::In => IntSurfTypeTrans::In,
                crate::int_tools_wline::TransType::Out => IntSurfTypeTrans::Out,
                crate::int_tools_wline::TransType::Touch => IntSurfTypeTrans::Touch,
                crate::int_tools_wline::TransType::Undecided => IntSurfTypeTrans::Undecided,
            };
            IntSurfTransition::with_type(false, ty)
        };
        if p.on_dom_s1 {
            // The flag effect of `SetArc(OnFirst=true, A, Param, TLine, TArc)`
            // (`IntPatch_Point.cxx:57-85`); the parameter on the arc is not
            // part of `PatchPoint`.
            r.on_s1 = true;
            r.arc_s1 = arc1;
            r.tra_line1 = trans(p.trans1);
            r.tra1 = trans(p.trans1);
        }
        if p.on_dom_s2 {
            r.on_s2 = true;
            r.arc_s2 = arc2;
            r.tra_line2 = trans(p.trans2);
            r.tra2 = trans(p.trans2);
        }
        r
    }
}
