//! `BRepCheck` — shape validity analyzer.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepCheck`. [`Analyzer`] walks the
//! shape (`BRepCheck_Analyzer`) and records a [`Status`] per sub-shape. The
//! geometric controls reuse [`crate::shell_check`] and [`crate::shape_checks`]
//! rather than a parallel heuristic.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::bop_occt_util::{explore, iter_children, shape_key};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_check::shell_is_closed;
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::topo_tools::{faces_of, wires_of_face, ShapeMapKey};

/// `BRepCheck_Status` (`BRepCheck_Status.hxx:20-59`), declaration order kept.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    NoError,
    InvalidPointOnCurve,
    InvalidPointOnCurveOnSurface,
    InvalidPointOnSurface,
    No3DCurve,
    Multiple3DCurve,
    Invalid3DCurve,
    NoCurveOnSurface,
    InvalidCurveOnSurface,
    InvalidCurveOnClosedSurface,
    InvalidSameRangeFlag,
    InvalidSameParameterFlag,
    InvalidDegeneratedFlag,
    FreeEdge,
    InvalidMultiConnexity,
    InvalidRange,
    EmptyWire,
    RedundantEdge,
    SelfIntersectingWire,
    NoSurface,
    InvalidWire,
    RedundantWire,
    IntersectingWires,
    InvalidImbricationOfWires,
    EmptyShell,
    RedundantFace,
    InvalidImbricationOfShells,
    UnorientableShape,
    NotClosed,
    NotConnected,
    SubshapeNotInShape,
    BadOrientation,
    BadOrientationOfSubshape,
    InvalidPolygonOnTriangulation,
    InvalidToleranceValue,
    EnclosedRegion,
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

/// The `BRep_CurveRepresentation` `BRepCheck_Edge::Minimum` keeps in `myCref`
/// (`BRepCheck_Edge.hxx:68`): the first representation whose curve is non-null.
enum Cref {
    Curve3D(Arc<dyn Curve>),
    CurveOnSurface(Arc<dyn Curve2d>),
}

/// `Geom_Curve::IsPeriodic/Period/FirstParameter/LastParameter` with the
/// `Geom_TrimmedCurve` unwrap of `BRepCheck_Edge::Minimum`
/// (`BRepCheck_Edge.cxx:175-189`): on a trimmed curve the *basis* supplies the
/// range and the period.
///
/// Deviation: OCCT unwraps one `Geom_TrimmedCurve` level (`BasisCurve()`, same
/// `cxx:179-188` as the `DynamicType` test above it), while `Curve::untrimmed_basis`
/// unwraps nested trims the way `GeomAdaptor_Curve::load` (`cxx:252-254`) does,
/// so a doubly-trimmed curve answers the periodicity query from its innermost
/// basis.
fn basis_window_3d(c: &dyn Curve) -> (bool, f64, f64, f64) {
    match c.untrimmed_basis() {
        Some((basis, _, _)) => (
            basis.is_periodic(),
            basis.period(),
            basis.first_parameter(),
            basis.last_parameter(),
        ),
        None => (c.is_periodic(), c.period(), c.first_parameter(), c.last_parameter()),
    }
}

/// 2D counterpart of [`basis_window_3d`] (`BRepCheck_Edge.cxx:219-232`):
/// `Geom2d_TrimmedCurve::BasisCurve()` supplies the range and the period.
fn basis_window_2d(c: &dyn Curve2d) -> (bool, f64, f64, f64) {
    match c.trimmed_basis() {
        Some(basis) => (
            basis.is_periodic(),
            basis.period(),
            basis.first_parameter(),
            basis.last_parameter(),
        ),
        None => (c.is_periodic(), c.period(), c.first_parameter(), c.last_parameter()),
    }
}

/// The first `BRep_CurveOnSurface` representation of an edge and its
/// `BRep_GCurve::First/Last` (`BRepCheck_Edge.cxx:132-141`): a seam edge's
/// forward pcurve, taken in `BRep_TEdge::myCurves` order.
fn first_curve_on_surface(g: &EdgeGeom) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let (key, curves) = g.pcurves.iter().next()?;
    let pc = curves.first()?.clone();
    let (f, l) = g
        .pcurve_ranges
        .get(key)
        .copied()
        .unwrap_or((pc.first_parameter(), pc.last_parameter()));
    Some((pc, f, l))
}

/// `BRepCheck_Edge::Minimum` (`BRepCheck_Edge.cxx:74-258`).
///
/// `myCref` is `BRep_Curve3D` when the edge has a 3D curve, otherwise the first
/// `BRep_CurveOnSurface` (the pcurve-only state that
/// `BRep_Builder::Degenerated(E, true)` installs, `BRep_Builder.cxx:1081-1085`,
/// and that this port models as `EdgeGeomCore::curve3d_null`).
///
/// UNPORTED, part of `Minimum`:
/// * the per-representation `TopLoc_Location` (`cxx:163-171`): this port's edge
///   representations carry no location, so `TransformedParameter` is not
///   applied to the range before the basis comparison; and
/// * `myHCurve` (the `GeomAdaptor_Curve` / `Adaptor3d_CurveOnSurface` kept for
///   `InContext`): `InContext` itself is not ported here.
fn edge_check(e: &Edge) -> Vec<Status> {
    let reg = GeometryRegistry::global();
    let geom = reg.edge_geom(e);
    let mut lst: Vec<Status> = Vec::new();

    // `TE->Degenerated/SameParameter/SameRange` (`cxx:91-93`).
    let degenerated = reg.is_degenerated_edge(e);
    let same_parameter = reg.same_parameter(e);
    let same_range = reg.same_range(e);
    if !same_range && same_parameter {
        // `cxx:94-97`.
        lst.push(Status::InvalidSameParameterFlag);
    }

    // Existence and uniqueness of a 3D representation (`cxx:100-128`).
    // `EdgeGeomCore::curve` is the single `BRep_Curve3D` this model can hold, so
    // `exist` is "the representation is present" (a present-but-null curve
    // counts, as `UpdateCurves` keeps the representation, `BRep_Builder.cxx:85-92`)
    // and `unique` is always true: `Multiple3DCurve` (`cxx:126-128`) is
    // unreachable, not skipped.
    let exists = geom.is_some();
    if !exists {
        lst.push(Status::No3DCurve);
    }

    // `myCref` (`cxx:100-141`); `GeometryRegistry::edge_curve` answers what
    // `BRep_Tool::Curve` answers, i.e. null on a null-curve representation.
    let mut cref = reg.edge_curve(e).map(Cref::Curve3D);
    let mut range = reg.edge_parameters(e);
    if cref.is_none() && !degenerated {
        if let Some((pc, f, l)) = geom.as_ref().and_then(first_curve_on_surface) {
            cref = Some(Cref::CurveOnSurface(pc));
            range = (f, l);
        }
    } else if cref.is_some() && degenerated {
        // `cxx:144-147`.
        lst.push(Status::InvalidDegeneratedFlag);
    }

    if let Some(c) = cref {
        let (first, last) = range;
        if last <= first {
            // `cxx:157-160`.
            lst.push(Status::InvalidRange);
        } else {
            let (is_periodic, period, f, l) = match &c {
                Cref::Curve3D(c3d) => basis_window_3d(c3d.as_ref()),
                Cref::CurveOnSurface(pc) => basis_window_2d(pc.as_ref()),
            };
            if is_periodic && (last - first > period + PCONFUSION) {
                // `cxx:190-194` / `cxx:233-237`.
                lst.push(Status::InvalidRange);
            } else if !is_periodic && (first < f - PCONFUSION || last > l + PCONFUSION) {
                // `cxx:195-199` / `cxx:238-242`.
                lst.push(Status::InvalidRange);
            }
        }
    }

    if lst.is_empty() {
        lst.push(Status::NoError);
    }
    lst
}

fn wire_check(w: &Wire) -> Vec<Status> {
    let edges = crate::topo_tools::edges_of_wire(w);
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
    // `BRepCheck_Shell::Closed` (`BRepCheck_Shell.cxx:404-438`) is the single
    // control flow for shell closure; `shape_checks` delegates here too.
    if !shell_is_closed(sh) {
        return vec![Status::NotClosed];
    }
    vec![Status::NoError]
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
    results: HashMap<ShapeMapKey, CheckResult>,
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
