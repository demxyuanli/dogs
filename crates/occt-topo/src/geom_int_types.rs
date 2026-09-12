//! GeomInt / IntPatch line types used by LineConstructor and IntSS.

use std::sync::Arc;

use occt_core::gp::{GpCirc, GpElips, GpHypr, GpLin, GpParab};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::int_tools_wline::{PatchPoint, WLine};
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
    pub kind: GLineKind,
    pub vertices: Vec<PatchPoint>,
    pub has_first_point: bool,
    pub has_last_point: bool,
}

impl GLine {
    pub fn new(kind: GLineKind) -> Self {
        Self {
            kind,
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
        }
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

    /// `IntPatch_GLine::AddVertex`.
    pub fn add_vertex(&mut self, p: PatchPoint) {
        self.vertices.push(p);
    }
}

/// Restriction line (`IntPatch_RLine`).
#[derive(Clone)]
pub struct RLine {
    pub vertices: Vec<PatchPoint>,
    pub has_first_point: bool,
    pub has_last_point: bool,
    pub arc_on_s1: bool,
    pub arc_on_s2: bool,
    pub c2d: Option<Arc<dyn Curve2d>>,
    pub param_f: f64,
    pub param_l: f64,
    pub uv_arc: Option<RestrictionArc>,
}

impl RLine {
    /// Restriction along `arc` on S1 (`on_first`) or S2.
    pub fn from_arc(arc: RestrictionArc, on_first: bool) -> Self {
        Self {
            vertices: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            arc_on_s1: on_first,
            arc_on_s2: !on_first,
            c2d: Some(arc.to_curve2d()),
            param_f: arc.first,
            param_l: arc.last,
            uv_arc: Some(arc),
        }
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
}

/// Analytic line (`IntPatch_ALine`) as a 3D evaluator plus vertices.
#[derive(Clone)]
pub struct ALine {
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
}
