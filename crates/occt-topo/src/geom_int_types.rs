//! GeomInt / IntPatch line types used by LineConstructor and IntSS.

use std::sync::Arc;

use occt_core::gp::{GpCirc, GpElips, GpHypr, GpLin, GpParab};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::intsurf_transition::{IntSurfSituation, IntSurfTypeTrans};
use crate::int_tools_wline::{PatchPoint, PntOn2S, WLine};
use crate::intpatch_line::IntPatchLine;
use super::topol::RestrictionArc;

/// `IntPatch_IType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntPatchIType {
    Lin,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    Analytic,
    Walking,
    Restriction,
}

/// Geometric line (`IntPatch_GLine`).
#[derive(Debug, Clone)]
pub enum GLineKind {
    Lin(GpLin),
    Circ(GpCirc),
    Elips(GpElips),
    Parab(GpParab),
    Hypr(GpHypr),
}

/// `IntPatch_GLine`.
#[derive(Debug, Clone)]
pub struct GLine {
    /// `IntPatch_Line` base state (`hxx:35-137`).
    pub base: IntPatchLine,
    pub kind: GLineKind,
    pub vertices: Vec<PatchPoint>,
    pub has_first_point: bool,
    pub has_last_point: bool,
    /// `myIndFirst` / `myIndLast` (`IntPatch_GLine.cxx:SetFirstPoint`).
    pub first_index: i32,
    pub last_index: i32,
}

impl GLine {
    pub fn new(kind: GLineKind) -> Self {
        Self {
            base: IntPatchLine::default(),
            kind,
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            first_index: 0,
            last_index: 0,
        }
    }

    /// `SetFirstPoint(IndFirst)` (`IntPatch_GLine.cxx`).
    pub fn set_first_point(&mut self, ind_first: i32) {
        self.has_first_point = true;
        self.first_index = ind_first;
    }

    /// `SetLastPoint(IndLast)` (`IntPatch_GLine.cxx`).
    pub fn set_last_point(&mut self, ind_last: i32) {
        self.has_last_point = true;
        self.last_index = ind_last;
    }

    /// `FirstPoint()` (`IntPatch_GLine.cxx`).
    pub fn first_point(&self) -> &PatchPoint {
        self.vertex(self.first_index)
    }

    /// `LastPoint()` (`IntPatch_GLine.cxx`).
    pub fn last_point(&self) -> &PatchPoint {
        self.vertex(self.last_index)
    }

    pub fn ityp(&self) -> IntPatchIType {
        match self.kind {
            GLineKind::Lin(_) => IntPatchIType::Lin,
            GLineKind::Circ(_) => IntPatchIType::Circle,
            GLineKind::Elips(_) => IntPatchIType::Ellipse,
            GLineKind::Parab(_) => IntPatchIType::Parabola,
            GLineKind::Hypr(_) => IntPatchIType::Hyperbola,
        }
    }

    pub fn nb_vertex(&self) -> i32 {
        self.vertices.len() as i32
    }

    pub fn vertex(&self, i1: i32) -> &PatchPoint {
        &self.vertices[(i1 as usize).saturating_sub(1)]
    }

    /// 1-based mutable vertex access (`ChangeValue`).
    pub fn vertex_mut(&mut self, i1: i32) -> &mut PatchPoint {
        let i = (i1 as usize).saturating_sub(1);
        &mut self.vertices[i]
    }

    /// `IntPatch_GLine::AddVertex`.
    pub fn add_vertex(&mut self, p: PatchPoint) {
        self.vertices.push(p);
    }
}

/// Restriction line (`IntPatch_RLine`).
#[derive(Clone)]
pub struct RLine {
    /// `IntPatch_Line` base state (`hxx:35-137`).
    pub base: IntPatchLine,
    /// `curv` (`IntPatch_RLine.hxx:182`): the 2D/3D point sequence of the
    /// line. Only `Add(LineOn2S)` (`IntPatch_RLine.cxx`) fills it, which is
    /// what [`crate::intpatch_line_ctor::wline_to_rline`] does when it encodes
    /// a walking line as a restriction line (`cxx:1357`); `ProcessSegments`
    /// builds an `RLine` from its arc alone and leaves it empty, as OCCT does
    /// with its `curv`-less constructor.
    pub points: Vec<PntOn2S>,
    pub vertices: Vec<PatchPoint>,
    pub has_first_point: bool,
    pub has_last_point: bool,
    /// `myIndFirst` / `myIndLast` (`IntPatch_PointLine.cxx:SetFirstPoint`).
    pub first_index: i32,
    pub last_index: i32,
    pub arc_on_s1: bool,
    pub arc_on_s2: bool,
    pub c2d: Option<Arc<dyn Curve2d>>,
    pub param_f: f64,
    pub param_l: f64,
    pub uv_arc: Option<RestrictionArc>,
    /// `BRepAdaptor_Curve2d::Edge()` of OCCT's `myArcOnS1` / `myArcOnS2`
    /// (`IntPatch_RLine.hxx:myArcOnS1`). OCCT reads that edge back through the
    /// stored adaptor (`TopOpeBRep_LineInter.cxx:534-547`); the port stores the
    /// p-curve instead of the adaptor, so it keeps the edge next to it.
    pub arc_edge: Option<crate::shape::Edge>,
}

impl RLine {
    /// Restriction along `arc` on S1 (`on_first`) or S2.
    pub fn from_arc(arc: RestrictionArc, on_first: bool) -> Self {
        Self {
            base: IntPatchLine::default(),
            points: Vec::new(),
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            first_index: 0,
            last_index: 0,
            arc_on_s1: on_first,
            arc_on_s2: !on_first,
            c2d: Some(arc.to_curve2d()),
            param_f: arc.first,
            param_l: arc.last,
            uv_arc: Some(arc),
            arc_edge: None,
        }
    }

    /// The adaptor `IntPatch_RLine` hands to its points: a
    /// `BRepAdaptor_Curve2d` over the arc's p-curve (`arc_edge` is its
    /// `Edge()`), or `None` while the edge is unknown.
    pub fn arc_adaptor(&self) -> Option<Arc<crate::brep_adaptor_curve2d::Curve2dAdaptor>> {
        let e = self.arc_edge.as_ref()?;
        let c = self.c2d.as_ref()?;
        Some(Arc::new(
            crate::brep_adaptor_curve2d::Curve2dAdaptor::from_stored(
                e,
                c.clone(),
                self.param_f,
                self.param_l,
            ),
        ))
    }

    pub fn nb_vertex(&self) -> i32 {
        self.vertices.len() as i32
    }

    pub fn vertex(&self, i1: i32) -> &PatchPoint {
        &self.vertices[(i1 as usize).saturating_sub(1)]
    }

    /// `IntPatch_RLine::AddVertex`.
    pub fn add_vertex(&mut self, p: PatchPoint) {
        self.vertices.push(p);
    }

    /// `SetFirstPoint(IndFirst)` (`IntPatch_PointLine.cxx`).
    pub fn set_first_point(&mut self, ind_first: i32) {
        self.has_first_point = true;
        self.first_index = ind_first;
    }

    /// `SetLastPoint(IndLast)` (`IntPatch_PointLine.cxx`).
    pub fn set_last_point(&mut self, ind_last: i32) {
        self.has_last_point = true;
        self.last_index = ind_last;
    }

    /// `FirstPoint()` (`IntPatch_PointLine.cxx`).
    pub fn first_point(&self) -> &PatchPoint {
        self.vertex(self.first_index)
    }

    /// `LastPoint()` (`IntPatch_PointLine.cxx`).
    pub fn last_point(&self) -> &PatchPoint {
        self.vertex(self.last_index)
    }

    /// `IntPatch_RLine::Add(LineOn2S)` (`IntPatch_RLine.cxx:120-123`):
    /// `curv = L`.
    pub fn add_points(&mut self, points: Vec<PntOn2S>) {
        self.points = points;
    }

    /// `NbPnts()` (`IntPatch_RLine.cxx`): `curv->NbPoints()`.
    pub fn nb_pnts(&self) -> i32 {
        self.points.len() as i32
    }

    /// 1-based `Point(i)`.
    pub fn point(&self, i1: i32) -> &PntOn2S {
        &self.points[(i1 as usize).saturating_sub(1)]
    }
}

/// Analytic line (`IntPatch_ALine`) as a 3D evaluator plus vertices.
#[derive(Clone)]
pub struct ALine {
    /// `IntPatch_Line` base state (`hxx:35-137`).
    pub base: IntPatchLine,
    pub curve: Arc<dyn Curve>,
    pub vertices: Vec<PatchPoint>,
    pub has_first_point: bool,
    pub has_last_point: bool,
    pub first_index: i32,
    pub last_index: i32,
}

impl ALine {
    pub fn from_curve(curve: Arc<dyn Curve>) -> Self {
        Self {
            base: IntPatchLine::default(),
            curve,
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            first_index: 0,
            last_index: 0,
        }
    }

    pub fn nb_vertex(&self) -> i32 {
        self.vertices.len() as i32
    }

    pub fn vertex(&self, i1: i32) -> &PatchPoint {
        &self.vertices[(i1 as usize).saturating_sub(1)]
    }

    pub fn value(&self, t: f64) -> occt_core::gp::GpPnt {
        self.curve.d0(t)
    }

    /// `IntPatch_ALine::AddVertex`.
    pub fn add_vertex(&mut self, p: PatchPoint) {
        self.vertices.push(p);
    }

    /// `IntPatch_ALine::Replace` (1-based).
    pub fn replace(&mut self, k: i32, p: PatchPoint) {
        let i = (k as usize).saturating_sub(1);
        if i < self.vertices.len() {
            self.vertices[i] = p;
        }
    }

    /// `IntPatch_ALine::SetFirstPoint` (1-based).
    pub fn set_first_point(&mut self, i1: i32) {
        self.has_first_point = true;
        self.first_index = i1;
    }

    /// `IntPatch_ALine::SetLastPoint` (1-based).
    pub fn set_last_point(&mut self, i1: i32) {
        self.has_last_point = true;
        self.last_index = i1;
    }

    /// `IntPatch_ALine::FirstPoint()` (`IntPatch_ALine.lxx`).
    pub fn first_point(&self) -> &PatchPoint {
        self.vertex(self.first_index)
    }

    /// `IntPatch_ALine::LastPoint()` (`IntPatch_ALine.lxx`).
    pub fn last_point(&self) -> &PatchPoint {
        self.vertex(self.last_index)
    }
}

/// `IntPatch_Line` sum type.
#[derive(Clone)]
pub enum GeomIntLine {
    Walking(WLine),
    Geometric(GLine),
    Analytic(ALine),
    Restriction(RLine),
}

impl GeomIntLine {
    pub fn arc_type(&self) -> IntPatchIType {
        match self {
            GeomIntLine::Walking(_) => IntPatchIType::Walking,
            GeomIntLine::Geometric(g) => g.ityp(),
            GeomIntLine::Analytic(_) => IntPatchIType::Analytic,
            GeomIntLine::Restriction(_) => IntPatchIType::Restriction,
        }
    }

    pub fn as_wline(&self) -> Option<&WLine> {
        match self {
            GeomIntLine::Walking(w) => Some(w),
            _ => None,
        }
    }

    pub fn as_gline(&self) -> Option<&GLine> {
        match self {
            GeomIntLine::Geometric(g) => Some(g),
            _ => None,
        }
    }

    pub fn as_gline_mut(&mut self) -> Option<&mut GLine> {
        match self {
            GeomIntLine::Geometric(g) => Some(g),
            _ => None,
        }
    }

    pub fn as_aline(&self) -> Option<&ALine> {
        match self {
            GeomIntLine::Analytic(a) => Some(a),
            _ => None,
        }
    }

    pub fn as_aline_mut(&mut self) -> Option<&mut ALine> {
        match self {
            GeomIntLine::Analytic(a) => Some(a),
            _ => None,
        }
    }

    pub fn as_rline(&self) -> Option<&RLine> {
        match self {
            GeomIntLine::Restriction(r) => Some(r),
            _ => None,
        }
    }

    pub fn as_rline_mut(&mut self) -> Option<&mut RLine> {
        match self {
            GeomIntLine::Restriction(r) => Some(r),
            _ => None,
        }
    }

    /// `NbVertex()` of whichever line kind this is (`IntPatch_Line` family).
    pub fn nb_vertex(&self) -> i32 {
        match self {
            GeomIntLine::Walking(w) => w.nb_vertex(),
            GeomIntLine::Geometric(g) => g.nb_vertex(),
            GeomIntLine::Analytic(a) => a.nb_vertex(),
            GeomIntLine::Restriction(r) => r.nb_vertex(),
        }
    }

    /// `Vertex(I)` (1-based) of whichever line kind this is.
    pub fn vertex(&self, i: i32) -> Option<&PatchPoint> {
        if i < 1 || i > self.nb_vertex() {
            return None;
        }
        match self {
            GeomIntLine::Walking(w) => Some(w.vertex(i)),
            GeomIntLine::Geometric(g) => Some(g.vertex(i)),
            GeomIntLine::Analytic(a) => Some(a.vertex(i)),
            GeomIntLine::Restriction(r) => Some(r.vertex(i)),
        }
    }

    /// The `IntPatch_Line` base state of whichever line kind this is
    /// (`IntPatch_Line.hxx:35-137`).
    pub fn base(&self) -> &IntPatchLine {
        match self {
            GeomIntLine::Walking(w) => &w.base,
            GeomIntLine::Geometric(g) => &g.base,
            GeomIntLine::Analytic(a) => &a.base,
            GeomIntLine::Restriction(r) => &r.base,
        }
    }

    /// Mutable [`Self::base`].
    pub fn base_mut(&mut self) -> &mut IntPatchLine {
        match self {
            GeomIntLine::Walking(w) => &mut w.base,
            GeomIntLine::Geometric(g) => &mut g.base,
            GeomIntLine::Analytic(a) => &mut a.base,
            GeomIntLine::Restriction(r) => &mut r.base,
        }
    }

    /// `IntPatch_Line::IsTangent()` (`lxx:33-37`).
    pub fn is_tangent(&self) -> bool {
        self.base().is_tangent()
    }

    /// `IntPatch_Line::TransitionOnS1()` (`lxx:39-43`).
    pub fn transition_on_s1(&self) -> IntSurfTypeTrans {
        self.base().transition_on_s1()
    }

    /// `IntPatch_Line::TransitionOnS2()` (`lxx:45-49`).
    pub fn transition_on_s2(&self) -> IntSurfTypeTrans {
        self.base().transition_on_s2()
    }

    /// `IntPatch_Line::SituationS1()` (`lxx:51-59`).
    pub fn situation_s1(&self) -> Result<IntSurfSituation, String> {
        self.base().situation_s1()
    }

    /// `IntPatch_Line::SituationS2()` (`lxx:61-69`).
    pub fn situation_s2(&self) -> Result<IntSurfSituation, String> {
        self.base().situation_s2()
    }
}
