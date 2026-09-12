//! `BRepCheck` — shape validity analyzer.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepCheck`. [`Analyzer`] walks the
//! shape (`BRepCheck_Analyzer`) and records a [`Status`] per sub-shape. The
//! geometric controls reuse [`crate::shell_check`] and [`crate::shape_checks`]
//! rather than a parallel heuristic.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_occt_util::{explore, iter_children, shape_key};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_check::shell_is_closed;
use crate::shape_checks;
use crate::topo_tools_full::{edges_of, faces_of, wires_of_face};

/// `BRepCheck_Status`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    NoError,
    No3DCurve,
    Invalid3DCurve,
    NoCurveOnSurface,
    InvalidRange,
    EmptyWire,
    RedundantEdge,
    NoSurface,
    InvalidWire,
    EmptyShell,
    NotClosed,
    NotConnected,
    UnorientableShape,
    CheckFail,
}

/// `BRepCheck_Result` for one sub-shape.
#[derive(Debug, Clone)]
pub struct CheckResult {
    pub shape: TopoShape,
    pub statuses: Vec<Status>,
}

impl CheckResult {
    /// True when every recorded status is [`Status::NoError`].
    pub fn is_valid(&self) -> bool {
        self.statuses.iter().all(|s| *s == Status::NoError)
    }
}

fn vertex_check(v: &Vertex) -> Vec<Status> {
    let _ = BRepTool::vertex_point(v);
    vec![Status::NoError]
}

fn edge_check(e: &Edge) -> Vec<Status> {
    if BRepTool::is_degenerated(e) {
        return vec![Status::NoError];
    }
    match BRepTool::edge_curve(e) {
        None => vec![Status::No3DCurve],
        Some(c) => {
            let (a, b) = BRepTool::edge_parameters(e);
            if !a.is_finite() || !b.is_finite() || b < a {
                vec![Status::InvalidRange]
            } else {
                let _ = c.d0(0.5 * (a + b));
                vec![Status::NoError]
            }
        }
    }
}

fn wire_check(w: &Wire) -> Vec<Status> {
    let edges = crate::topo_tools_full::edges_of_wire(w);
    if edges.is_empty() {
        return vec![Status::EmptyWire];
    }
    vec![Status::NoError]
}

fn face_check(f: &Face) -> Vec<Status> {
    if BRepTool::face_surface(f).is_none() {
        return vec![Status::NoSurface];
    }
    let wires = wires_of_face(f);
    if wires.is_empty() {
        return vec![Status::InvalidWire];
    }
    for w in &wires {
        if wire_check(w).iter().any(|s| *s != Status::NoError) {
            return vec![Status::InvalidWire];
        }
    }
    vec![Status::NoError]
}

fn shell_check_one(sh: &Shell) -> Vec<Status> {
    let faces = faces_of(&sh.0);
    if faces.is_empty() {
        return vec![Status::EmptyShell];
    }
    if !shell_is_closed(sh) {
        return vec![Status::NotClosed];
    }
    match shape_checks::shell_closed_check(sh) {
        shape_checks::CheckLevel::Ok => vec![Status::NoError],
        _ => vec![Status::NotClosed],
    }
}

fn solid_check(s: &Solid) -> Vec<Status> {
    let shells = explore(&s.0, ShapeType::Shell);
    if shells.is_empty() {
        return vec![Status::EmptyShell];
    }
    for sh in &shells {
        if !shell_is_closed(&Shell(sh.clone())) {
            return vec![Status::NotClosed];
        }
    }
    vec![Status::NoError]
}

/// `BRepCheck_Analyzer`.
#[derive(Debug, Clone)]
pub struct Analyzer {
    shape: TopoShape,
    results: HashMap<usize, CheckResult>,
}

impl Analyzer {
    /// `BRepCheck_Analyzer(S, GeomControls=true)`.
    pub fn new(shape: TopoShape) -> Self {
        let mut a = Self {
            shape: shape.clone(),
            results: HashMap::new(),
        };
        a.perform(&shape);
        a
    }

    fn perform(&mut self, s: &TopoShape) {
        self.check_one(s);
        for c in iter_children(s) {
            self.perform(&c);
        }
    }

    fn check_one(&mut self, s: &TopoShape) {
        let k = shape_key(s);
        if self.results.contains_key(&k) {
            return;
        }
        let statuses = match s.shape_type() {
            ShapeType::Vertex => vertex_check(&Vertex(s.clone())),
            ShapeType::Edge => edge_check(&Edge(s.clone())),
            ShapeType::Wire => wire_check(&Wire(s.clone())),
            ShapeType::Face => face_check(&Face(s.clone())),
            ShapeType::Shell => shell_check_one(&Shell(s.clone())),
            ShapeType::Solid | ShapeType::CompSolid => solid_check(&Solid(s.clone())),
            _ => vec![Status::NoError],
        };
        self.results.insert(
            k,
            CheckResult {
                shape: s.clone(),
                statuses,
            },
        );
    }

    /// `IsValid()`.
    pub fn is_valid(&self) -> bool {
        self.results.values().all(CheckResult::is_valid)
    }

    /// Statuses recorded for `s`.
    pub fn result(&self, s: &TopoShape) -> Option<&CheckResult> {
        self.results.get(&shape_key(s))
    }

    /// All recorded results.
    pub fn results(&self) -> impl Iterator<Item = &CheckResult> {
        self.results.values()
    }

    /// The shape that was analyzed.
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }
}

/// `BRepCheck_Analyzer::IsValid` convenience.
pub fn is_valid(shape: &TopoShape) -> bool {
    Analyzer::new(shape.clone()).is_valid()
}

/// Face-only check (`BRepCheck_Face`).
pub fn check_face(f: &Face) -> CheckResult {
    CheckResult {
        shape: f.0.clone(),
        statuses: face_check(f),
    }
}

/// Shell-only check (`BRepCheck_Shell`).
pub fn check_shell(sh: &Shell) -> CheckResult {
    CheckResult {
        shape: sh.0.clone(),
        statuses: shell_check_one(sh),
    }
}

/// Solid-only check (`BRepCheck_Solid`).
pub fn check_solid(s: &Solid) -> CheckResult {
    CheckResult {
        shape: s.0.clone(),
        statuses: solid_check(s),
    }
}

/// Wire-only check (`BRepCheck_Wire`).
pub fn check_wire(w: &Wire) -> CheckResult {
    CheckResult {
        shape: w.0.clone(),
        statuses: wire_check(w),
    }
}

/// Edge-only check (`BRepCheck_Edge`).
pub fn check_edge(e: &Edge) -> CheckResult {
    CheckResult {
        shape: e.0.clone(),
        statuses: edge_check(e),
    }
}

/// Unused: keep a status that maps to redundant edges when a wire repeats an edge.
#[allow(dead_code)]
fn wire_redundant(w: &Wire) -> Option<Status> {
    let mut seen = HashMap::new();
    for e in edges_of(&w.0) {
        *seen.entry(shape_key(&e.0)).or_insert(0) += 1;
    }
    if seen.values().any(|&n| n > 2) {
        Some(Status::RedundantEdge)
    } else {
        None
    }
}
