//! Unified BRepBuilderAPI construction surface — MakeEdge, MakeWire, MakeFace,
//! MakeShell, MakeSolid, MakePolygon and MakeVertex as ergonomic builders.
//! Source: `BRepBuilderAPI_MakeEdge`, `BRepBuilderAPI_MakeWire`,
//! `BRepBuilderAPI_MakeFace`, `BRepBuilderAPI_MakeShell`,
//! `BRepBuilderAPI_MakeSolid`, `BRepBuilderAPI_MakePolygon`,
//! `BRepBuilderAPI_MakeVertex` (TKTopAlgo).
//!
//! The port is organised as a small set of namespace structs, one per OCCT
//! builder class, plus a handful of free convenience functions. Each builder
//! mirrors the constructor overloads of its OCCT counterpart:
//!
//! * [`BRepBuilderEdge`] — `BRepBuilderAPI_MakeEdge` (points, curves, circles,
//!   ellipses, arcs through three points, B-splines).
//! * [`BRepBuilderWire`] — `BRepBuilderAPI_MakeWire` (validated chains of
//!   edges, open polylines, auto-closed wires).
//! * [`BRepBuilderFace`] — `BRepBuilderAPI_MakeFace` (planar faces from wires,
//!   faces from a surface and a trimming wire, bare surfaces).
//! * [`BRepBuilderShell`] — `BRepBuilderAPI_MakeShell` (faces → shell, closed
//!   box shells, closedness queries).
//! * [`BRepBuilderSolid`] — `BRepBuilderAPI_MakeSolid` (shells → solid,
//!   faces → solid, box primitives).
//! * [`BRepBuilderPolygon`] — `BRepBuilderAPI_MakePolygon` (open/closed
//!   polyline wires).
//!
//! Where an equivalent already exists in [`crate::brep_builder_api`] the new
//! builder delegates to it rather than duplicating geometry logic; this module
//! only *adds* the missing constructor surface (validated wire assembly, the
//! ellipse/B-spline edge constructors, box shells, two-vertex edges and the
//! auto-closing wire helper).
//!
//! # Worked example
//!
//! Building a closed square box solid, then recovering its geometry:
//!
//! ```text
//! use occt_core::gp::{GpAx2, GpPnt};
//! use occt_topo::brep_builder_full::{
//!     BRepBuilderEdge, BRepBuilderFace, BRepBuilderPolygon, BRepBuilderSolid,
//!     BRepBuilderWire,
//! };
//!
//! // Edge chain of a unit square.
//! let e1 = BRepBuilderEdge::from_points(GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.));
//! let e2 = BRepBuilderEdge::from_points(GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.));
//! let e3 = BRepBuilderEdge::from_points(GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.));
//! let e4 = BRepBuilderEdge::from_points(GpPnt::new(0., 1., 0.), GpPnt::new(0., 0., 0.));
//! let square = BRepBuilderWire::from_edges(&[e1, e2, e3, e4]);   // closed wire
//! ```
//!
//! Every builder validates its inputs and returns `Result<_, String>` with a
//! descriptive error, so a malformed construction (coincident points, a
//! disjoint wire chain, a degenerate box) fails loudly instead of producing a
//! broken shape. See [`wire_from_edges_auto`] for a best-effort variant that
//! bridges gaps automatically.
//!
//! # Delegation
//!
//! The low-level shape building is done through [`TopoBuilder`], so every
//! shape returned here carries registered geometry in the
//! [`GeometryRegistry`](crate::tgeometry::GeometryRegistry) side-table and can
//! be queried back through [`BRepTool`]. Existing polygon / arc / primitive
//! helpers in [`crate::brep_builder_api`] are reused verbatim — for example
//! [`BRepBuilderEdge::from_arc_3pts`] and [`BRepBuilderPolygon::from_points`]
//! are thin delegations to `brep_builder_api::make_edge_arc` /
//! `make_polygon`.

use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir, GpElips, GpLin, GpPln, GpPnt, GpVec};
use occt_geom::{Curve, GeomBSplineCurve, GeomEllipse, GeomLine, GeomPlane, Surface};

use crate::brep_builder_api;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::topo_tools_full;

/// Tolerance (model units) used when validating that consecutive wire edges
/// share an endpoint, and when deciding whether a chain of edges is closed.
///
/// The value matches the endpoint-matching tolerance used throughout the
/// topology helpers (`topo_tools_full::wire_is_closed`, `sweep` ring
/// reconstruction) so a wire built here also reads as closed there.
const WIRE_TOL: f64 = 1e-9;

/// Edge construction — a port of `BRepBuilderAPI_MakeEdge`.
///
/// Provides the most common edge constructors: a straight segment between two
/// points, an edge over an arbitrary curve on a parameter range, circular and
/// elliptic arcs, an arc through three points, and a B-spline edge through a
/// set of poles. Every constructor registers the edge's curve in the geometry
/// side-table, so [`BRepTool::edge_curve`] and [`BRepTool::edge_vertices`]
/// recover it later exactly like OCCT's `BRep_Tool::Curve`.
///
/// The constructors mirror the OCCT overloads by argument shape:
///
/// | OCCT constructor                     | this port                        |
/// |--------------------------------------|----------------------------------|
/// | `MakeEdge(P1, P2)`                   | [`Self::from_points`]            |
/// | `MakeEdge(curve, first, last)`       | [`Self::from_curve`]             |
/// | `MakeEdge(gp_Circ, first, last)`     | [`Self::from_circle`]            |
/// | `MakeEdge(gp_Elips, first, last)`    | [`Self::from_ellipse`]           |
/// | `MakeEdge(P1, P2, P3)`               | [`Self::from_arc_3pts`]          |
/// | `MakeEdge(Handle(Geom_BSplineCurve))`| [`Self::from_bspline`]           |
///
/// Point-based constructors attach the endpoint vertices as edge children
/// (shared, registered points); curve-parameter constructors leave vertex
/// attachment to the caller, matching OCCT's distinction between edges built
/// from points and edges built from an already-parameterised curve.
#[derive(Debug, Clone, Copy, Default)]
pub struct BRepBuilderEdge;

impl BRepBuilderEdge {
    /// Straight segment edge between two distinct points
    /// (`BRepBuilderAPI_MakeEdge(P1, P2)`).
    ///
    /// The edge's curve is a `GeomLine` parameterised on `[0, |P2−P1|]`, and
    /// the two endpoints are attached as vertex children (shared, registered
    /// points), matching the OCCT behaviour.
    ///
    /// # Errors
    ///
    /// Returns `Err` when the two points coincide — a segment of zero length
    /// has no unique direction.
    pub fn from_points(p1: GpPnt, p2: GpPnt) -> Result<Edge, String> {
        if p1.distance(&p2) <= 1e-12 {
            return Err("BRepBuilderEdge::from_points: coincident points".into());
        }
        Ok(TopoBuilder::new().make_edge_segment(&p1, &p2))
    }

    /// Edge over an arbitrary curve on the parameter range `[first, last]`
    /// (`BRepBuilderAPI_MakeEdge(Handle(Geom_Curve), first, last)`).
    ///
    /// The caller owns the parameterisation; the curve is registered as-is and
    /// no vertex children are attached. Use [`Self::from_points`] or
    /// [`Self::from_arc_3pts`] when endpoint vertices are needed.
    pub fn from_curve(curve: Arc<dyn Curve>, first: f64, last: f64) -> Edge {
        TopoBuilder::new().make_edge(curve, first, last)
    }

    /// Circular arc edge in the plane `ax` with the given radius
    /// (`BRepBuilderAPI_MakeEdge(gp_Circ, first, last)`).
    ///
    /// `ax` is the circle's placement: its Z axis is the circle normal and its
    /// X axis points at the zero-angle point. `first`/`last` are the angular
    /// parameters in radians (a full circle is `0..2π`).
    pub fn from_circle(ax: GpAx2, r: f64, first: f64, last: f64) -> Edge {
        TopoBuilder::new().make_edge_circle(&ax, r, first, last)
    }

    /// Elliptic arc edge with the given major and minor radii
    /// (`BRepBuilderAPI_MakeEdge(gp_Elips, first, last)`).
    ///
    /// The ellipse lies in the plane `ax`: its major axis is `ax`'s X
    /// direction, its minor axis the Y direction, and `first`/`last` are the
    /// angular parameters in radians.
    ///
    /// # Errors
    ///
    /// Returns `Err` when either radius is non-positive or when the major
    /// radius is smaller than the minor radius (a degenerate or
    /// inconsistently-oriented ellipse).
    pub fn from_ellipse(ax: GpAx2, major: f64, minor: f64, first: f64, last: f64) -> Result<Edge, String> {
        if major <= 0.0 || minor <= 0.0 {
            return Err("BRepBuilderEdge::from_ellipse: radii must be positive".into());
        }
        if major < minor {
            return Err("BRepBuilderEdge::from_ellipse: major radius must be >= minor radius".into());
        }
        let e = GeomEllipse::new(GpElips::new(ax, major, minor));
        Ok(TopoBuilder::new().make_edge(Arc::new(e), first, last))
    }

    /// Circular arc through three non-collinear points
    /// (`BRepBuilderAPI_MakeEdge(P1, P2, P3)`).
    ///
    /// The curve is the circle through `p1`, `p2`, `p3` trimmed to the minor
    /// arc that starts at `p1`, passes through `p2` and ends at `p3`. The
    /// registered curve is a `GeomTrimmedCurve` normalised to `[0, 1]`, so the
    /// middle point sits at parameter `0.5`. Endpoint vertices at `p1` and
    /// `p3` are attached.
    ///
    /// # Errors
    ///
    /// Returns `Err` when the three points are collinear (no circle passes
    /// through them) or when `p1` coincides with the circumcenter.
    pub fn from_arc_3pts(p1: GpPnt, p2: GpPnt, p3: GpPnt) -> Result<Edge, String> {
        brep_builder_api::make_edge_arc(&p1, &p2, &p3)
    }

    /// B-spline edge through the given poles
    /// (`BRepBuilderAPI_MakeEdge(Handle(Geom_BSplineCurve))`).
    ///
    /// Builds a non-rational `GeomBSplineCurve` from `poles`, `knots` and
    /// `degree` and registers it on its full natural parameter range. The knot
    /// vector must be clamped (first `degree+1` knots equal, last `degree+1`
    /// knots equal) and satisfy the OCCT count rule
    /// `nb_knots = nb_poles + degree + 1`.
    ///
    /// # Errors
    ///
    /// Returns `Err` when the knot/pole/degree counts are inconsistent (the
    /// underlying `GeomBSplineCurve::new` validation).
    pub fn from_bspline(poles: &[GpPnt], knots: &[f64], degree: usize) -> Result<Edge, String> {
        let curve = GeomBSplineCurve::new(poles.to_vec(), knots.to_vec(), degree)
            .map_err(|e| format!("BRepBuilderEdge::from_bspline: {e}"))?;
        let first = curve.first_parameter();
        let last = curve.last_parameter();
        Ok(TopoBuilder::new().make_edge(Arc::new(curve), first, last))
    }
}

/// Wire construction — a port of `BRepBuilderAPI_MakeWire`.
///
/// A wire is an ordered chain of edges. This builder validates that
/// consecutive edges actually share an endpoint (within [`WIRE_TOL`]) so that
/// a malformed chain is rejected up front rather than silently producing a
/// broken boundary. The struct keeps the running edge list together with the
/// first-start / last-end points so it can both validate incrementally and
/// report closedness without re-walking the chain.
///
/// ```text
/// // A triangle wire: three segments, corner vertices shared.
/// let a = GpPnt::new(0., 0., 0.);
/// let b = GpPnt::new(1., 0., 0.);
/// let c = GpPnt::new(0., 1., 0.);
/// let mut w = BRepBuilderWire::new();
/// w.add_edge(BRepBuilderEdge::from_points(a, b)?)?;  // -> Result<(), String>
/// w.add_edge(BRepBuilderEdge::from_points(b, c)?)?;
/// w.add_edge(BRepBuilderEdge::from_points(c, a)?)?;  // closes the loop
/// let wire = w.build()?;
/// ```
///
/// The incremental [`Self::add_edge`] path is useful when edges arrive one at
/// a time; the batch constructors [`Self::from_edges`] and [`Self::closed`]
/// validate and assemble a whole slice at once. Both paths produce wires whose
/// `closed` flag reflects whether the first start meets the last end.
#[derive(Debug, Clone, Default)]
pub struct BRepBuilderWire {
    /// Edges accumulated so far, in traversal order.
    pub edges: Vec<Edge>,
    /// Start point of the first edge (the wire's beginning).
    first_start: Option<GpPnt>,
    /// End point of the most recently added edge (the wire's current end).
    last_end: Option<GpPnt>,
}

impl BRepBuilderWire {
    /// A new empty wire builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append an edge, validating that its start meets the previous edge's end.
    ///
    /// The first edge is accepted unconditionally. Every subsequent edge must
    /// start within [`WIRE_TOL`] of the previously added edge's end, mirroring
    /// the OCCT `MakeWire::Add` error "the edge is not connected to the wire".
    ///
    /// # Errors
    ///
    /// Returns `Err` when the edge has no registered curve (so its endpoints
    /// cannot be evaluated) or when it does not connect to the chain.
    pub fn add_edge(&mut self, e: Edge) -> Result<(), String> {
        let (start, end) = BRepTool::edge_vertices(&e)
            .ok_or("BRepBuilderWire::add_edge: edge has no registered curve")?;
        if let Some(le) = self.last_end {
            if start.distance(&le) > WIRE_TOL {
                return Err(format!(
                    "BRepBuilderWire::add_edge: consecutive edges are disjoint (gap {})",
                    start.distance(&le)
                ));
            }
        }
        if self.first_start.is_none() {
            self.first_start = Some(start);
        }
        self.last_end = Some(end);
        self.edges.push(e);
        Ok(())
    }

    /// Assemble the accumulated edges into a wire.
    ///
    /// The wire's `closed` flag is set when the first edge's start coincides
    /// with the last edge's end. This is the completion step for the
    /// incremental [`Self::add_edge`] path.
    ///
    /// # Errors
    ///
    /// Returns `Err` when no edges have been added.
    pub fn build(&self) -> Result<Wire, String> {
        if self.edges.is_empty() {
            return Err("BRepBuilderWire::build: no edges".into());
        }
        let wire = TopoBuilder::new().make_wire(&self.edges);
        let closed = match (self.first_start, self.last_end) {
            (Some(fs), Some(le)) => fs.distance(&le) <= WIRE_TOL,
            _ => false,
        };
        wire.set_closed(closed);
        Ok(wire)
    }

    /// Wire from a validated chain of consecutive edges
    /// (`BRepBuilderAPI_MakeWire(edge1, edge2, ...)`).
    ///
    /// Every edge must connect to the next; an open chain is allowed (it
    /// produces an open wire), a chain whose last end meets its first start
    /// produces a closed wire.
    ///
    /// # Errors
    ///
    /// Returns `Err` on an empty slice or when any consecutive pair of edges
    /// is disjoint.
    pub fn from_edges(edges: &[Edge]) -> Result<Wire, String> {
        if edges.is_empty() {
            return Err("BRepBuilderWire::from_edges: no edges".into());
        }
        let mut bld = BRepBuilderWire::new();
        for e in edges {
            bld.add_edge(e.clone())?;
        }
        bld.build()
    }

    /// Open polyline wire through the given points
    /// (`BRepBuilderAPI_MakePolygon`).
    ///
    /// Produces `n − 1` straight segment edges through `n` points, sharing the
    /// corner vertices. The wire is open (no closing edge back to the first
    /// point); see [`BRepBuilderPolygon::from_points`] for the closed form.
    ///
    /// # Errors
    ///
    /// Returns `Err` when fewer than 2 points are given.
    pub fn from_points(pts: &[GpPnt]) -> Result<Wire, String> {
        brep_builder_api::make_wire_from_points(pts)
    }

    /// Closed wire from a chain of edges, auto-adding the closing segment.
    ///
    /// If the chain already closes (last end meets first start), it is
    /// returned as-is. Otherwise a straight segment edge from the last end to
    /// the first start is appended, mirroring the OCCT behaviour of calling
    /// `MakeWire::Add` with an explicit closing edge.
    ///
    /// # Errors
    ///
    /// Returns `Err` on an empty slice, a disjoint chain, or an edge without a
    /// registered curve.
    pub fn closed(edges: &[Edge]) -> Result<Wire, String> {
        if edges.is_empty() {
            return Err("BRepBuilderWire::closed: no edges".into());
        }
        let mut bld = BRepBuilderWire::new();
        for e in edges {
            bld.add_edge(e.clone())?;
        }
        let already_closed = match (bld.first_start, bld.last_end) {
            (Some(fs), Some(le)) => fs.distance(&le) <= WIRE_TOL,
            _ => false,
        };
        if already_closed {
            return bld.build();
        }
        let fs = bld.first_start.ok_or("BRepBuilderWire::closed: no start point")?;
        let le = bld.last_end.ok_or("BRepBuilderWire::closed: no end point")?;
        bld.add_edge(TopoBuilder::new().make_edge_segment(&le, &fs))?;
        bld.build()
    }
}

/// Face construction — a port of `BRepBuilderAPI_MakeFace`.
///
/// Faces are built from a surface plus a set of boundary wires (or no wires
/// for an unbounded natural face). Planar faces can be constructed directly
/// from a wire and a plane, or from a coplanar polygon of points.
///
/// The OCCT `MakeFace` class distinguishes three construction paths, all
/// present here:
///
/// * a plane (or any surface) plus one or more boundary wires — the trimmed
///   face; [`Self::from_wire`], [`Self::from_surface_and_wire`];
/// * a set of planar points that are turned into the boundary automatically —
///   [`Self::from_planar_points`];
/// * a bare surface with no trimming — [`Self::from_surface`] (the
///   natural-restriction face, e.g. a full sphere).
///
/// A face does not validate that its wire lies on the surface (that is a
/// separate `BRepCheck` pass in OCCT); it simply pairs the surface with the
/// wire so `BRepTool::face_surface` and `BRepTool::uv_bounds` can recover the
/// geometry.
#[derive(Debug, Clone, Copy, Default)]
pub struct BRepBuilderFace;

impl BRepBuilderFace {
    /// Planar face bounded by `wire`, lying on `plane`
    /// (`BRepBuilderAPI_MakeFace(gp_Pln, wire)`).
    ///
    /// The face surface is a `GeomPlane` on the given plane and the wire is
    /// attached as its outer boundary. The wire's vertices are *not* required
    /// to lie exactly on the plane — the caller chooses the plane — but a wire
    /// that wanders far off the plane produces a geometrically inconsistent
    /// face.
    ///
    /// # Errors
    ///
    /// Returns `Err` when the wire contains no edges.
    pub fn from_wire(wire: &Wire, plane: &GpPln) -> Result<Face, String> {
        if topo_tools_full::edges_of_wire(wire).is_empty() {
            return Err("BRepBuilderFace::from_wire: wire has no edges".into());
        }
        Ok(TopoBuilder::new().make_face(Arc::new(GeomPlane::new(plane.clone())), &[wire.clone()]))
    }

    /// Planar face from a coplanar polygon of points
    /// (`BRepBuilderAPI_MakeFace(P1, P2, P3, ..., 1e-6)`).
    ///
    /// Computes the plane through the first three non-collinear points,
    /// validates that every point is coplanar within tolerance, and builds the
    /// closed polygon wire and its planar face.
    ///
    /// # Errors
    ///
    /// Returns `Err` when the points are all collinear, are not coplanar, or
    /// form a degenerate (zero-area) polygon.
    pub fn from_planar_points(pts: &[GpPnt]) -> Result<Face, String> {
        brep_builder_api::make_face_from_polygon(pts)
    }

    /// Face over `surface` bounded by `wire` (`BRepBuilderAPI_MakeFace(surface, wire)`).
    ///
    /// The surface is registered as the face's geometry and the wire is
    /// attached as the trimming boundary (e.g. a circle wire on a cylinder
    /// surface for a pipe face).
    ///
    /// # Errors
    ///
    /// Returns `Err` when the wire contains no edges.
    pub fn from_surface_and_wire(surface: Arc<dyn Surface>, wire: &Wire) -> Result<Face, String> {
        if topo_tools_full::edges_of_wire(wire).is_empty() {
            return Err("BRepBuilderFace::from_surface_and_wire: wire has no edges".into());
        }
        Ok(TopoBuilder::new().make_face(surface, &[wire.clone()]))
    }

    /// Unbounded face over `surface` with no trimming wires
    /// (`BRepBuilderAPI_MakeFace(surface)`).
    ///
    /// Equivalent to OCCT's natural-restriction face: the whole parametric
    /// surface is the face.
    pub fn from_surface(surface: Arc<dyn Surface>) -> Face {
        TopoBuilder::new().make_face(surface, &[])
    }
}

/// Shell construction — a port of `BRepBuilderAPI_MakeShell`.
///
/// Wraps a [`Shell`] (a set of faces forming a boundary) together with
/// closedness queries. A closed shell is a manifold boundary in which every
/// edge is used by exactly two faces.
#[derive(Debug, Clone)]
pub struct BRepBuilderShell {
    /// The wrapped shell.
    pub shell: Shell,
}

impl BRepBuilderShell {
    /// Wrap an existing shell.
    pub fn new(shell: Shell) -> Self {
        Self { shell }
    }

    /// Shell containing `faces` (`BRepBuilderAPI_MakeShell(face1, face2, ...)`).
    ///
    /// # Errors
    ///
    /// Returns `Err` when `faces` is empty.
    pub fn from_faces(faces: &[Face]) -> Result<Shell, String> {
        if faces.is_empty() {
            return Err("BRepBuilderShell::from_faces: no faces".into());
        }
        Ok(TopoBuilder::new().make_shell(faces))
    }

    /// Closed box shell with dimensions `w × h × d` (`BRepPrimAPI_MakeBox`).
    ///
    /// Builds the full box primitive and returns its single shell: 6 planar
    /// faces, 12 line edges, 8 vertices, all with registered geometry. The
    /// shell is guaranteed closed (`every edge used by 2 faces`).
    ///
    /// # Errors
    ///
    /// Returns `Err` when any dimension is non-positive.
    pub fn closed_box(w: f64, h: f64, d: f64) -> Result<Shell, String> {
        if w <= 0.0 || h <= 0.0 || d <= 0.0 {
            return Err("BRepBuilderShell::closed_box: dimensions must be positive".into());
        }
        let box_ = crate::primitives::BRepPrimBox::make_box(w, h, d);
        let shell = box_.solid.0.tshape.read().unwrap().children[0].clone();
        Ok(Shell(shell))
    }

    /// Whether the wrapped shell is a closed manifold boundary
    /// (`BRepCheck_Shell::Closed`).
    ///
    /// Delegates to [`crate::shell_check::shell_is_closed`]: true when every
    /// boundary edge is referenced by exactly 2 faces.
    pub fn is_closed(&self) -> bool {
        crate::shell_check::shell_is_closed(&self.shell)
    }
}

/// Solid construction — a port of `BRepBuilderAPI_MakeSolid`.
///
/// Assembles shells into solids and provides the box primitive convenience.
///
/// A solid in the OCCT data model is a set of shells (usually exactly one
/// closed boundary shell). [`Self::from_shell`] wraps a single shell,
/// [`Self::from_faces`] builds the shell implicitly from the given faces, and
/// [`Self::box_corners`] constructs a full axis-aligned box primitive with
/// registered geometry (8 vertices, 12 line edges, 6 planar faces).
///
/// As with OCCT, assembling a solid from faces does *not* check that the
/// faces form a closed manifold — use
/// [`crate::shell_check::shell_is_closed`] on the resulting shell to verify
/// watertightness. The box constructors always produce a closed solid.
#[derive(Debug, Clone, Copy, Default)]
pub struct BRepBuilderSolid;

impl BRepBuilderSolid {
    /// Solid containing a single `shell` (`BRepBuilderAPI_MakeSolid(shell)`).
    ///
    /// # Errors
    ///
    /// Returns `Err` when the shell has no faces.
    pub fn from_shell(shell: &Shell) -> Result<Solid, String> {
        if topo_tools_full::faces_of(&shell.0).is_empty() {
            return Err("BRepBuilderSolid::from_shell: shell has no faces".into());
        }
        Ok(TopoBuilder::new().make_solid(&[shell.clone()]))
    }

    /// Solid from a set of faces: one shell containing them, then the solid
    /// (`BRepBuilderAPI_MakeSolid`).
    ///
    /// The faces should already form a closed manifold boundary if the result
    /// is meant to be a watertight solid.
    ///
    /// # Errors
    ///
    /// Returns `Err` when `faces` is empty.
    pub fn from_faces(faces: &[Face]) -> Result<Solid, String> {
        brep_builder_api::make_solid_from_faces(faces)
    }

    /// Axis-aligned box solid spanning two opposite corner points
    /// (`BRepPrimAPI_MakeBox(p1, p2)`).
    ///
    /// The corners may be given in any order; the box spans the axis-aligned
    /// bounding box between them. Delegates to
    /// [`crate::primitives::BRepPrimBox::make_box_corner`].
    ///
    /// # Errors
    ///
    /// Returns `Err` when the two corners share a coordinate plane (zero
    /// extent in any axis).
    pub fn box_corners(p1: GpPnt, p2: GpPnt) -> Result<Solid, String> {
        let dx = (p2.x() - p1.x()).abs();
        let dy = (p2.y() - p1.y()).abs();
        let dz = (p2.z() - p1.z()).abs();
        if dx <= 1e-12 || dy <= 1e-12 || dz <= 1e-12 {
            return Err("BRepBuilderSolid::box_corners: degenerate box".into());
        }
        Ok(crate::primitives::BRepPrimBox::make_box_corner(&p1, &p2).solid)
    }
}

/// Polygon wire construction — a port of `BRepBuilderAPI_MakePolygon`.
///
/// A polygon is a polyline wire whose edges are straight segments. `closed`
/// controls whether a closing edge is added back to the first point.
///
/// This is the wire-level half of `BRepBuilderAPI_MakePolygon`: the OCCT class
/// also exposes the accumulated vertices and edges, which
/// [`crate::brep_builder_api::make_polygon`] returns as a
/// [`PolygonBuilder`](crate::brep_builder_api::PolygonBuilder). When only the
/// wire is needed, [`Self::from_points`] is the ergonomic entry point; the
/// richer struct is available through the delegated `brep_builder_api` call.
#[derive(Debug, Clone, Copy, Default)]
pub struct BRepBuilderPolygon;

impl BRepBuilderPolygon {
    /// Polyline wire through `pts` (`BRepBuilderAPI_MakePolygon`).
    ///
    /// When `closed` is true the wire has `n` edges (a closing edge returns to
    /// the first point) and is flagged closed; otherwise it has `n − 1` edges
    /// and is open. Corner vertices are shared between adjacent edges.
    ///
    /// # Errors
    ///
    /// Returns `Err` when there are fewer than 2 points for an open polygon or
    /// fewer than 3 for a closed one.
    pub fn from_points(pts: &[GpPnt], closed: bool) -> Result<Wire, String> {
        if closed {
            brep_builder_api::make_polygon(pts).map(|p| p.wire)
        } else {
            brep_builder_api::make_wire_from_points(pts)
        }
    }
}

/// Convenience vertex constructor (`BRepBuilderAPI_MakeVertex`).
///
/// Registers a zero-tolerance vertex at `p` and returns it wrapped as a
/// [`Vertex`]. Equivalently `TopoBuilder::new().make_vertex(p, 0.0)`.
///
/// A vertex built here can later be queried through [`BRepTool::vertex_point`]
/// and can be shared by multiple edges (see [`edge_from_two_vertices`]).
pub fn make_vertex(p: GpPnt) -> Vertex {
    TopoBuilder::new().make_vertex(p, 0.0)
}

/// Straight segment edge between two existing vertices.
///
/// Reuses the two vertices' registered points (via [`BRepTool::vertex_point`])
/// and attaches the *same* vertex shapes as the edge's children, so the edge
/// shares its boundary vertices with any other shape using them — the OCCT
/// `BRepBuilderAPI_MakeEdge` shared-vertex model.
///
/// The edge's curve is a `GeomLine` parameterised on
/// `[0, distance(v1, v2)]`. Because the input vertex `TShape`s are reused
/// (not copied), this is the right constructor for building a wire out of
/// pre-created vertices while keeping the topology connected.
///
/// # Errors
///
/// Returns `Err` when the two vertices are coincident.
pub fn edge_from_two_vertices(v1: &Vertex, v2: &Vertex) -> Result<Edge, String> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    if p1.distance(&p2) <= 1e-12 {
        return Err("edge_from_two_vertices: coincident vertices".into());
    }
    let b = TopoBuilder::new();
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2))
        .map_err(|_| "edge_from_two_vertices: degenerate segment")?;
    let mut e = b.make_edge(
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))),
        0.0,
        p1.distance(&p2),
    );
    b.add(&mut e.0, &v1.0);
    b.add(&mut e.0, &v2.0);
    Ok(e)
}

/// Build a wire from a chain of edges, closing it automatically when needed.
///
/// The behaviour depends on the chain's connectivity (measured on edge curve
/// endpoints within [`WIRE_TOL`]):
///
/// * All consecutive edges connected **and** the chain closed by construction
///   (last end meets first start) → the wire is returned flagged closed.
/// * All consecutive edges connected but open → an open wire is returned.
/// * A gap exists between consecutive edges → the gap is bridged with straight
///   segment edges and, if the chain is still open, a closing edge is appended
///   back to the first start.
///
/// The bridging is *best-effort* (a documented simplification): the inserted
/// segments are straight, so a chain of curved edges whose endpoints do not
/// coincide gains linear bridging edges.
///
/// ```text
/// // Two disjoint segments get bridged into a single connected wire:
/// //   (0,0,0)-(1,0,0)  +  (5,0,0)-(6,0,0)
/// // → (0,0,0)-(1,0,0)-(5,0,0)-(6,0,0), with a (1,0,0)-(5,0,0) bridge and
/// //   a closing (6,0,0)-(0,0,0) edge (best-effort, straight segments).
/// let wire = wire_from_edges_auto(&[e1, e2])?;   // closed
/// ```
///
/// # Errors
///
/// Returns `Err` on an empty slice or when any edge lacks a registered curve.
pub fn wire_from_edges_auto(edges: &[Edge]) -> Result<Wire, String> {
    if edges.is_empty() {
        return Err("wire_from_edges_auto: no edges".into());
    }
    let endpoints: Vec<(GpPnt, GpPnt)> = edges
        .iter()
        .map(|e| {
            BRepTool::edge_vertices(e).ok_or("wire_from_edges_auto: edge has no registered curve")
        })
        .collect::<Result<_, _>>()?;
    let n = endpoints.len();
    let mut gaps: Vec<usize> = Vec::new();
    for i in 0..n.saturating_sub(1) {
        if endpoints[i].1.distance(&endpoints[i + 1].0) > WIRE_TOL {
            gaps.push(i);
        }
    }
    let naturally_closed = n >= 2 && endpoints[n - 1].1.distance(&endpoints[0].0) <= WIRE_TOL;

    let b = TopoBuilder::new();
    if gaps.is_empty() {
        // Already a connected chain: keep it as-is, flagging closedness.
        let wire = b.make_wire(edges);
        wire.set_closed(naturally_closed);
        return Ok(wire);
    }

    // ponytail: best-effort auto-close — bridge each gap with a segment edge
    // and append the closing edge when the chain is open.
    let mut out: Vec<Edge> = Vec::with_capacity(edges.len() + gaps.len() + 1);
    for i in 0..n {
        out.push(edges[i].clone());
        if i + 1 < n && endpoints[i].1.distance(&endpoints[i + 1].0) > WIRE_TOL {
            out.push(b.make_edge_segment(&endpoints[i].1, &endpoints[i + 1].0));
        }
    }
    if !naturally_closed {
        out.push(b.make_edge_segment(&endpoints[n - 1].1, &endpoints[0].0));
    }
    let wire = b.make_wire(&out);
    wire.set_closed(true);
    Ok(wire)
}

/// Planar face from a polygon of points on the given `plane`.
///
/// Builds the closed polygon wire through `pts` and attaches it to a
/// `GeomPlane` face on `plane`. Unlike [`BRepBuilderFace::from_planar_points`],
/// the plane is supplied by the caller rather than inferred from the points.
///
/// # Errors
///
/// Returns `Err` when fewer than 3 points are given (a face needs a closed
/// wire).
pub fn face_from_polygon(pts: &[GpPnt], plane: &GpPln) -> Result<Face, String> {
    if pts.len() < 3 {
        return Err("face_from_polygon: need at least 3 points".into());
    }
    let wire = BRepBuilderPolygon::from_points(pts, true)?;
    Ok(TopoBuilder::new().make_face(Arc::new(GeomPlane::new(plane.clone())), &[wire]))
}

/// Solid by extruding a base polygon along `dir` by `height`
/// (`BRepPrimAPI_MakePrism`).
///
/// The base polygon (≥ 3 points, boundary order) is turned into a planar base
/// face and swept: the result is a closed prism solid with `2 + n` faces.
/// `dir` is a direction; the displacement is `dir` scaled to `height`.
/// Delegates to [`crate::sweep::prism_from_polygon`].
///
/// # Errors
///
/// Returns `Err` when the base has fewer than 3 points, `height` is
/// non-positive, or `dir` is the zero vector.
pub fn solid_from_prism(base: &[GpPnt], dir: &GpVec, height: f64) -> Result<Solid, String> {
    if base.len() < 3 {
        return Err("solid_from_prism: need at least 3 base points".into());
    }
    if height <= 0.0 {
        return Err("solid_from_prism: height must be positive".into());
    }
    if dir.xyz().square_modulus() <= 1e-30 {
        return Err("solid_from_prism: zero sweep direction".into());
    }
    Ok(crate::sweep::prism_from_polygon(base, dir, height).solid)
}

// ============================================================================
// Compatibility notes
// ============================================================================
//
// The builders here map OCCT's `BRepBuilderAPI_*` classes to Rust `Result`
// returning free functions. Differences from OCCT worth knowing:
//
// * **No exception model.** OCCT's `MakeEdge(P1, P2)` raises `StdFail_NotDone`
//   when the points coincide; here the same condition is an `Err(String)`.
//   Every builder validates its preconditions up front and returns a
//   descriptive message instead of panicking (except where noted, e.g. the
//   low-level `BRepPrimBox` dimensions are asserted in `primitives.rs`).
//
// * **Validation is constructive.** `BRepBuilderWire` rejects a chain whose
//   consecutive edges are disjoint, whereas OCCT's `MakeWire` accepts any
//   edges and only reports `WireDone() == false` afterwards. The port
//   deliberately fails fast so a broken boundary is caught at construction
//   time.
//
// * **Vertex sharing.** Edges built from points share their endpoint vertex
//   `TShape`s across the wire (the OCCT shared-TShape model), so a box built
//   from `BRepBuilderWire` has 4 corner vertices, not 8 duplicate ones. This
//   is what makes `shell_euler_characteristic` come out to 2 for a closed box.
//
// * **Planes are unbounded.** A face built with [`BRepBuilderFace::from_wire`]
//   carries a `GeomPlane` whose natural parameter range is infinite; the
//   bounding wire defines the actual trimmed region, exactly as in OCCT.
//
// * **Closedness is a flag, not a guarantee.** [`BRepBuilderWire`] sets the
//   wire's `closed` flag based on endpoint coincidence within [`WIRE_TOL`].
//   The flag is advisory; use `topo_tools_full::wire_is_closed` for a
//   structural re-check on an arbitrary wire.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_measure::edge_length;
    use crate::brep_surface::{classify_surface, face_is_planar, SurfaceKind};
    use crate::fillet_edge::classify_surface_full;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::edges_of_wire;
    use occt_core::gp::{GpAx3, GpCylinder};
    use occt_geom::GeomCylinder;
    use std::f64::consts::PI;

    #[test]
    fn edge_from_points_and_curve() {
        // Straight segment of length 1.
        let e = BRepBuilderEdge::from_points(GpPnt::zero(), GpPnt::new(1.0, 0.0, 0.0)).expect("segment");
        assert!((edge_length(&e, 16) - 1.0).abs() < 1e-9, "segment length");
        let (f, l) = BRepTool::edge_parameters(&e);
        assert!((f - 0.0).abs() < 1e-12 && (l - 1.0).abs() < 1e-12);

        // Circle edge: parameters span first..last, curve sits on radius 2.
        let c = BRepBuilderEdge::from_circle(GpAx2::standard(), 2.0, 0.0, PI);
        let (f, l) = BRepTool::edge_parameters(&c);
        assert!((f - 0.0).abs() < 1e-12 && (l - PI).abs() < 1e-12);
        let curve = BRepTool::edge_curve(&c).expect("circle curve");
        assert!(curve.d0(0.0).distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-9);
        assert!(curve.d0(PI).distance(&GpPnt::new(-2.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn edge_from_ellipse() {
        let e = BRepBuilderEdge::from_ellipse(GpAx2::standard(), 4.0, 2.0, 0.0, 2.0 * PI).expect("ellipse edge");
        let curve = BRepTool::edge_curve(&e).expect("ellipse curve");
        // Major vertex at angle 0, minor vertex at π/2. The OCCT ellipse
        // convention is x = a·cos u, y = −b·sin u, so the π/2 point lies at
        // (0, −minor).
        assert!(curve.d0(0.0).distance(&GpPnt::new(4.0, 0.0, 0.0)) < 1e-9, "major vertex");
        assert!(curve.d0(0.5 * PI).distance(&GpPnt::new(0.0, -2.0, 0.0)) < 1e-9, "minor vertex");
        // Invalid radii rejected.
        assert!(BRepBuilderEdge::from_ellipse(GpAx2::standard(), 2.0, 3.0, 0.0, 1.0).is_err());
        assert!(BRepBuilderEdge::from_ellipse(GpAx2::standard(), 0.0, 2.0, 0.0, 1.0).is_err());
    }

    #[test]
    fn edge_arc_3pts() {
        let p1 = GpPnt::new(1.0, 0.0, 0.0);
        let p2 = GpPnt::new(0.0, 1.0, 0.0);
        let p3 = GpPnt::new(-1.0, 0.0, 0.0);
        let e = BRepBuilderEdge::from_arc_3pts(p1, p2, p3).expect("arc");
        let curve = BRepTool::edge_curve(&e).expect("arc curve");
        // The arc curve passes through the middle point at the half parameter.
        assert!(curve.d0(0.5).distance(&p2) < 1e-9, "arc midpoint at p2");
        assert!(curve.d0(0.0).distance(&p1) < 1e-9, "arc start at p1");
        assert!(curve.d0(1.0).distance(&p3) < 1e-9, "arc end at p3");
        // Endpoint vertex children match p1/p3.
        let (a, z) = crate::topo_tools_full::edge_vertices(&e);
        assert!(crate::topo_tools_full::vertex_position(&a.unwrap()).distance(&p1) < 1e-9);
        assert!(crate::topo_tools_full::vertex_position(&z.unwrap()).distance(&p3) < 1e-9);
    }

    #[test]
    fn edge_bspline() {
        let poles = vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(2., 0., 0.)];
        let knots = vec![0., 0., 0., 1., 1., 1.];
        let e = BRepBuilderEdge::from_bspline(&poles, &knots, 2).expect("bspline edge");
        let curve = BRepTool::edge_curve(&e).expect("bspline curve");
        let (f, l) = BRepTool::edge_parameters(&e);
        assert!(curve.d0(f).distance(&poles[0]) < 1e-6, "start at first pole");
        assert!(curve.d0(l).distance(&poles[2]) < 1e-6, "end at last pole");
        // Bad knot count rejected.
        assert!(BRepBuilderEdge::from_bspline(&poles, &[0., 0., 0., 1., 1.], 2).is_err());
    }

    #[test]
    fn wire_from_edges_validates() {
        let e1 = BRepBuilderEdge::from_points(GpPnt::zero(), GpPnt::new(1., 0., 0.)).unwrap();
        let e2 = BRepBuilderEdge::from_points(GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.)).unwrap();
        assert!(BRepBuilderWire::from_edges(&[e1.clone(), e2]).is_ok(), "connected chain ok");

        let e3 = BRepBuilderEdge::from_points(GpPnt::new(5., 0., 0.), GpPnt::new(6., 0., 0.)).unwrap();
        assert!(BRepBuilderWire::from_edges(&[e1, e3]).is_err(), "disjoint chain rejected");
        assert!(BRepBuilderWire::from_edges(&[]).is_err(), "empty rejected");
    }

    #[test]
    fn wire_from_points_and_closed() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let open = BRepBuilderWire::from_points(&pts).expect("open wire");
        assert_eq!(edges_of_wire(&open).len(), 3, "open: n-1 edges");
        assert!(!open.closed());

        let closed = BRepBuilderPolygon::from_points(&pts, true).expect("closed wire");
        assert_eq!(edges_of_wire(&closed).len(), 4, "closed: n edges");
        assert!(closed.closed());
        assert!(crate::topo_tools_full::wire_is_closed(&closed), "closed chains end-to-end");
    }

    #[test]
    fn wire_incremental_builder() {
        let a = GpPnt::new(0., 0., 0.);
        let b = GpPnt::new(1., 0., 0.);
        let c = GpPnt::new(1., 1., 0.);
        let mut w = BRepBuilderWire::new();
        w.add_edge(BRepBuilderEdge::from_points(a, b).unwrap()).expect("add e1");
        w.add_edge(BRepBuilderEdge::from_points(b, c).unwrap()).expect("add e2");
        // Disjoint third edge is rejected.
        let far = BRepBuilderEdge::from_points(GpPnt::new(9., 0., 0.), GpPnt::new(8., 0., 0.)).unwrap();
        assert!(w.add_edge(far).is_err(), "disjoint edge rejected");
        let wire = w.build().expect("build wire");
        assert_eq!(edges_of_wire(&wire).len(), 2);
        assert!(!wire.closed());
    }

    #[test]
    fn face_from_wire_planar() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let wire = BRepBuilderPolygon::from_points(&pts, true).expect("square wire");
        let plane = GpPln::new(GpAx3::standard());
        let face = BRepBuilderFace::from_wire(&wire, &plane).expect("planar face");
        assert!(face_is_planar(&face));
        let surf = BRepTool::face_surface(&face).expect("surface");
        assert_eq!(classify_surface(surf.as_ref()), SurfaceKind::Plane);
    }

    #[test]
    fn face_from_surface_and_wire() {
        // A cylinder face bounded by two circle wires (bottom + top).
        let b = TopoBuilder::new();
        let ax = GpAx2::standard();
        let r = 2.0;
        let height = 3.0;
        let bottom = b.make_edge_circle(&ax, r, 0.0, 2.0 * PI);
        let mut top_ax = ax;
        top_ax.set_location(GpPnt::new(0.0, 0.0, height));
        let top = b.make_edge_circle(&top_ax, r, 0.0, 2.0 * PI);
        let wire_b = b.make_wire(&[bottom]);
        let wire_t = b.make_wire(&[top]);
        let cyl = GpCylinder::new(GpAx3::standard(), r).expect("cylinder");
        let face = b.make_face(Arc::new(GeomCylinder::new(cyl)), &[wire_b.clone(), wire_t]);
        let surf = BRepTool::face_surface(&face).expect("surface");
        assert_eq!(classify_surface_full(surf.as_ref()), SurfaceKind::Cylinder);

        // Builder path with an explicit wire.
        let f2 = BRepBuilderFace::from_surface_and_wire(
            Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), r).unwrap())),
            &wire_b,
        )
        .expect("surface+wire face");
        assert!(BRepTool::face_surface(&f2).is_some());
    }

    #[test]
    fn shell_closed_box() {
        let shell = BRepBuilderShell::closed_box(1.0, 2.0, 3.0).expect("box shell");
        assert!(crate::shell_check::shell_is_closed(&shell), "closed box shell");
        assert_eq!(crate::topo_tools_full::faces_of(&shell.0).len(), 6);
        assert!(BRepBuilderShell::new(shell).is_closed());
        assert!(BRepBuilderShell::closed_box(0.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn solid_from_faces() {
        let bx = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = crate::topo_tools_full::faces_of(&bx.solid.0);
        let solid = BRepBuilderSolid::from_faces(&faces).expect("solid");
        let shell = Shell(solid.0.tshape.read().unwrap().children[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "solid shell closed");
        assert_eq!(crate::topo_tools_full::vertices_of(&solid.0).len(), 8);
        assert!(BRepBuilderSolid::from_faces(&[]).is_err());
    }

    #[test]
    fn polygon_open_closed() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.),
        ];
        let open = BRepBuilderPolygon::from_points(&pts, false).expect("open polygon");
        assert_eq!(edges_of_wire(&open).len(), pts.len() - 1, "open: n-1 edges");
        assert!(!open.closed());
        let closed = BRepBuilderPolygon::from_points(&pts, true).expect("closed polygon");
        assert_eq!(edges_of_wire(&closed).len(), pts.len(), "closed: n edges");
        assert!(closed.closed());
    }

    #[test]
    fn edge_between_vertices() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0., 0., 0.), 0.0);
        let v2 = b.make_vertex(GpPnt::new(2., 0., 0.), 0.0);
        let e = edge_from_two_vertices(&v1, &v2).expect("segment");
        assert!((edge_length(&e, 16) - 2.0).abs() < 1e-9, "length 2");
        let (a, z) = BRepTool::edge_vertices(&e).expect("endpoints");
        assert!(a.distance(&GpPnt::zero()) < 1e-9);
        assert!(z.distance(&GpPnt::new(2., 0., 0.)) < 1e-9);
        // The two input vertices are reused as the edge's boundary children.
        let (va, vz) = crate::topo_tools_full::edge_vertices(&e);
        assert!(crate::topo_tools_full::is_same(&va.unwrap().0, &v1.0));
        assert!(crate::topo_tools_full::is_same(&vz.unwrap().0, &v2.0));
    }

    #[test]
    fn prism_solid_from_base() {
        let tri = [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)];
        let solid = solid_from_prism(&tri, &GpVec::new(0., 0., 1.), 3.0).expect("prism");
        let shell = Shell(solid.0.tshape.read().unwrap().children[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "triangle prism closed");
        // Boundary counts: 6 vertices, 9 edges, 5 faces.
        assert_eq!(crate::topo_tools_full::vertices_of(&solid.0).len(), 6);
        assert_eq!(crate::topo_tools_full::edges_of(&solid.0).len(), 9);
        assert_eq!(crate::topo_tools_full::faces_of(&solid.0).len(), 5);
        // Analytic prism volume = base area × sweep height, recovered from the
        // solid's own geometry (the mesh-based volume is unreliable on a
        // triangle base: the face grid overhangs the boundary).
        let faces = crate::topo_tools_full::faces_of(&solid.0);
        let zc = |f: &Face| crate::brep_surface::face_centroid(f, 8, 8).map(|p| p.z()).unwrap_or(f64::NAN);
        let base = faces.iter().min_by(|a, b| zc(a).partial_cmp(&zc(b)).unwrap()).unwrap();
        let top = faces.iter().max_by(|a, b| zc(a).partial_cmp(&zc(b)).unwrap()).unwrap();
        let base_plane = crate::brep_surface::face_plane(base).expect("base plane");
        let top_plane = crate::brep_surface::face_plane(top).expect("top plane");
        let height = top_plane.location().z() - base_plane.location().z();
        assert!((height - 3.0).abs() < 1e-9, "sweep height {height}");
        let wire = crate::topo_tools_full::wires_of_face(base).into_iter().next().expect("base wire");
        let ring: Vec<GpPnt> = crate::topo_tools_full::edges_of_wire(&wire)
            .iter()
            .filter_map(|e| crate::topo_tools_full::edge_vertices(e).0)
            .map(|v| crate::topo_tools_full::vertex_position(&v))
            .collect();
        let area = crate::brep_builder_api::polygon_area(&ring);
        assert!((area - 0.5).abs() < 1e-9, "base area {area}");
        let vol = area * height;
        assert!(vol > 0.0, "volume positive");
        assert!((vol - 1.5).abs() < 1e-9, "volume {vol} ~ 1.5");
        // Invalid inputs rejected.
        assert!(solid_from_prism(&tri, &GpVec::new(0., 0., 1.), 0.0).is_err());
        assert!(solid_from_prism(&[GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.)], &GpVec::new(0., 0., 1.), 1.0).is_err());
        assert!(solid_from_prism(&tri, &GpVec::new(0., 0., 0.), 1.0).is_err());
    }

    #[test]
    fn collinear_arc_errors() {
        let e = BRepBuilderEdge::from_arc_3pts(
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(2., 0., 0.),
        );
        assert!(e.is_err(), "collinear points cannot define an arc");
    }
}
