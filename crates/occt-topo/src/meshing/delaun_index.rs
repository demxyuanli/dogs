//! Port of OCCT delaun_index — BRepMesh_{VertexTool, VertexInspector, CircleTool, CircleInspector}.
//!
//! Uniform-grid spatial indexes over 2D UV space. [`VertexTool`] hashes
//! [`DelaunVertex`]es into a `cells_u × cells_v` grid covering a [`Box2d`] and
//! answers "is there already a vertex within `tol` of this UV point?" so that
//! coincident mesh nodes collapse to a single index. [`CircleTool`] does the
//! same for circumcircles ([`DelaunCircle`]), keyed by circle center and
//! queried with "which circumcircles are shot by this point?". The matching
//! `*Inspector` types implement the per-cell predicates (coincidence within
//! tolerance / circle hit by a query point), mirroring OCCT's `NCollection_CellFilter`
//! inspection callbacks.
//!
//! Indices are 0-based (Rust style) rather than OCCT's 1-based arrays.

// Imported so the local stand-ins below are swapped for the real `delaun_types`
// types once that stub is filled (see ponytail note).
#[allow(unused_imports)]
use super::delaun_types::*;

// ---------------------------------------------------------------------------
// Value types (local stand-ins)
// ---------------------------------------------------------------------------
// ponytail: `delaun_types.rs` is still a stub — the value types below are the
// minimal local equivalents the indexes need. Once the Wave-2 Delaun types land
// they will be resolved by the `use super::delaun_types::*;` glob above and the
// local definitions can be deleted.

/// 2D axis-aligned box in UV (parametric) space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Box2d {
    pub min_u: f64,
    pub min_v: f64,
    pub max_u: f64,
    pub max_v: f64,
}

impl Box2d {
    /// Builds a box from its two opposite corners.
    pub const fn new(min_u: f64, min_v: f64, max_u: f64, max_v: f64) -> Self {
        Self { min_u, min_v, max_u, max_v }
    }

    /// Width along the U axis.
    pub fn width(&self) -> f64 {
        self.max_u - self.min_u
    }

    /// Height along the V axis.
    pub fn height(&self) -> f64 {
        self.max_v - self.min_v
    }

    /// True when the box has zero or negative extent on either axis.
    pub fn is_empty(&self) -> bool {
        self.max_u <= self.min_u || self.max_v <= self.min_v
    }
}

/// Movability of a mesh vertex / element. Maps to `BRepMesh_DegreeOfFreedom`;
/// only the four variants used by the Delaunay tools are modeled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Movability {
    /// `BRepMesh_Free`
    Free,
    /// `BRepMesh_Frontier`
    Frontier,
    /// `BRepMesh_Fixed`
    Fixed,
    /// `BRepMesh_Deleted`
    Deleted,
}

/// Mesh vertex in 2D UV space.
///
/// Minimal stand-in for `delaun_types::DelaunVertex` (which adds a 3D-location
/// index); see the ponytail note at the top of this module.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunVertex {
    pub u: f64,
    pub v: f64,
    pub movability: Movability,
}

impl DelaunVertex {
    /// Builds a free vertex at `(u, v)`.
    pub const fn new(u: f64, v: f64) -> Self {
        Self { u, v, movability: Movability::Free }
    }

    /// Builds a vertex with the given movability.
    pub const fn with_movability(u: f64, v: f64, movability: Movability) -> Self {
        Self { u, v, movability }
    }

    /// `(u, v)` coordinates.
    pub fn coord(&self) -> (f64, f64) {
        (self.u, self.v)
    }

    /// True when the vertex has been deleted and may be replaced.
    pub fn is_deleted(&self) -> bool {
        self.movability == Movability::Deleted
    }

    /// Sets the movability of the vertex.
    pub fn set_movability(&mut self, movability: Movability) {
        self.movability = movability;
    }

    /// Squared Euclidean distance to `(u, v)`.
    pub fn square_distance(&self, u: f64, v: f64) -> f64 {
        let du = self.u - u;
        let dv = self.v - v;
        du * du + dv * dv
    }
}

/// Circumcircle of a Delaunay triangle in UV space.
///
/// Minimal stand-in for `delaun_types::DelaunCircle`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunCircle {
    pub u: f64,
    pub v: f64,
    pub radius: f64,
}

impl DelaunCircle {
    /// Builds a circle with the given center and radius.
    pub const fn new(u: f64, v: f64, radius: f64) -> Self {
        Self { u, v, radius }
    }

    /// Center of the circle as `(u, v)`.
    pub fn center(&self) -> (f64, f64) {
        (self.u, self.v)
    }

    /// OCCT marks deleted circles with a negative radius.
    pub fn is_deleted(&self) -> bool {
        self.radius < 0.0
    }
}

// ---------------------------------------------------------------------------
// Cell grid
// ---------------------------------------------------------------------------

/// Uniform `cells_u × cells_v` grid over a [`Box2d`]; each cell holds the
/// indices of the entries hashed into it.
#[derive(Debug, Clone)]
struct CellGrid {
    box2d: Box2d,
    cells_u: usize,
    cells_v: usize,
    cells: Vec<Vec<usize>>,
}

impl CellGrid {
    fn new(box2d: Box2d, cells_u: usize, cells_v: usize) -> Self {
        let cells_u = cells_u.max(1);
        let cells_v = cells_v.max(1);
        Self { box2d, cells_u, cells_v, cells: vec![Vec::new(); cells_u * cells_v] }
    }

    #[inline]
    fn flat(&self, ci: usize, cj: usize) -> usize {
        cj * self.cells_u + ci
    }

    /// Grid cell containing the point `(u, v)`; out-of-box points clamp to the
    /// boundary cells.
    fn cell_of(&self, u: f64, v: f64) -> (usize, usize) {
        let w = self.box2d.width();
        let h = self.box2d.height();
        let fu = if w > 0.0 { (u - self.box2d.min_u) / w } else { 0.0 };
        let fv = if h > 0.0 { (v - self.box2d.min_v) / h } else { 0.0 };
        let nu = if fu.is_finite() { fu } else { 0.0 };
        let nv = if fv.is_finite() { fv } else { 0.0 };
        let ci = (nu * self.cells_u as f64).floor().clamp(0.0, (self.cells_u - 1) as f64) as usize;
        let cj = (nv * self.cells_v as f64).floor().clamp(0.0, (self.cells_v - 1) as f64) as usize;
        (ci, cj)
    }

    /// Inclusive cell range covering the box `[u0, u1] × [v0, v1]`.
    fn cell_range(&self, u0: f64, v0: f64, u1: f64, v1: f64) -> ((usize, usize), (usize, usize)) {
        let a = self.cell_of(u0, v0);
        let b = self.cell_of(u1, v1);
        ((a.0.min(b.0), a.1.min(b.1)), (a.0.max(b.0), a.1.max(b.1)))
    }

    /// Hashes `index` into the cell containing `(u, v)`.
    fn push(&mut self, u: f64, v: f64, index: usize) {
        let (ci, cj) = self.cell_of(u, v);
        let cell = self.flat(ci, cj);
        self.cells[cell].push(index);
    }

    /// Removes `index` from the cell containing `(u, v)`.
    fn remove(&mut self, u: f64, v: f64, index: usize) {
        let (ci, cj) = self.cell_of(u, v);
        let cell = self.flat(ci, cj);
        self.cells[cell].retain(|&i| i != index);
    }

    /// Indices hashed into the given cell.
    fn entries(&self, ci: usize, cj: usize) -> &[usize] {
        &self.cells[self.flat(ci, cj)]
    }
}

// ---------------------------------------------------------------------------
// VertexInspector / VertexTool
// ---------------------------------------------------------------------------

/// Cell-query predicate for vertex coincidence within tolerance.
///
/// Port of `BRepMesh_VertexInspector`: given a reference point it reports
/// whether a stored vertex lies within the (radial or per-axis) tolerance and
/// keeps track of the closest match.
#[derive(Debug, Clone)]
pub struct VertexInspector {
    tolerance_sq: f64,
    tolerance_u_sq: f64,
    tolerance_v_sq: f64,
    use_axis: bool,
    target_u: f64,
    target_v: f64,
    best_index: Option<usize>,
    best_sq_dist: f64,
}

impl VertexInspector {
    /// Builds an inspector with the given radial tolerance.
    pub fn new(tolerance: f64) -> Self {
        let mut inspector = Self {
            tolerance_sq: 0.0,
            tolerance_u_sq: 0.0,
            tolerance_v_sq: 0.0,
            use_axis: false,
            target_u: 0.0,
            target_v: 0.0,
            best_index: None,
            best_sq_dist: f64::MAX,
        };
        inspector.set_tolerance(tolerance);
        inspector
    }

    /// Sets a single radial tolerance used for both dimensions.
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tolerance_sq = tolerance * tolerance;
        self.use_axis = false;
    }

    /// Sets separate tolerances for the U and V dimensions.
    pub fn set_tolerance_xy(&mut self, tol_u: f64, tol_v: f64) {
        self.tolerance_u_sq = tol_u * tol_u;
        self.tolerance_v_sq = tol_v * tol_v;
        self.use_axis = true;
    }

    /// Sets the reference point and resets the best-match state.
    pub fn set_point(&mut self, u: f64, v: f64) {
        self.target_u = u;
        self.target_v = v;
        self.best_index = None;
        self.best_sq_dist = f64::MAX;
    }

    /// Inspects one vertex against the reference point. Deleted vertices are
    /// ignored. Returns `true` when the vertex lies within the tolerance.
    pub fn inspect(&mut self, index: usize, vertex: &DelaunVertex) -> bool {
        if vertex.is_deleted() {
            return false;
        }
        let du = vertex.u - self.target_u;
        let dv = vertex.v - self.target_v;
        let in_tol = if self.use_axis {
            du * du <= self.tolerance_u_sq && dv * dv <= self.tolerance_v_sq
        } else {
            du * du + dv * dv <= self.tolerance_sq
        };
        if in_tol {
            let sq = du * du + dv * dv;
            if sq < self.best_sq_dist {
                self.best_sq_dist = sq;
                self.best_index = Some(index);
            }
        }
        in_tol
    }

    /// Index of the closest matching vertex, if any.
    pub fn best(&self) -> Option<usize> {
        self.best_index
    }
}

/// Spatial index of UV vertices keyed by their 2D coordinates.
///
/// Port of `BRepMesh_VertexTool`: vertices are hashed into a uniform cell grid
/// over a [`Box2d`]; [`VertexTool::add_vertex`] deduplicates nodes coincident
/// within a tolerance so each UV position maps to a single index.
#[derive(Debug, Clone)]
pub struct VertexTool {
    grid: CellGrid,
    vertices: Vec<DelaunVertex>,
    del_nodes: Vec<usize>,
    inspector: VertexInspector,
}

impl VertexTool {
    /// Builds a `cells_u × cells_v` grid covering `box2d`.
    pub fn new(box2d: Box2d, cells_u: usize, cells_v: usize) -> Self {
        Self {
            grid: CellGrid::new(box2d, cells_u, cells_v),
            vertices: Vec::new(),
            del_nodes: Vec::new(),
            inspector: VertexInspector::new(1e-7),
        }
    }

    /// The indexed vertices (deleted slots included, marked [`Movability::Deleted`]).
    pub fn vertices(&self) -> &[DelaunVertex] {
        &self.vertices
    }

    /// Number of indexed vertices (deleted slots included).
    pub fn extent(&self) -> usize {
        self.vertices.len()
    }

    /// True when no vertex has been indexed.
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Vertex with the given index.
    pub fn vertex(&self, index: usize) -> Option<&DelaunVertex> {
        self.vertices.get(index)
    }

    /// Grid cell containing the given UV point.
    pub fn get_cell(&self, u: f64, v: f64) -> (usize, usize) {
        self.grid.cell_of(u, v)
    }

    /// Indices of the vertices hashed into the given cell.
    pub fn select_cell(&self, ci: usize, cj: usize) -> &[usize] {
        self.grid.entries(ci, cj)
    }

    /// The box covered by the grid.
    pub fn box2d(&self) -> Box2d {
        self.grid.box2d
    }

    /// Grid dimensions `(cells_u, cells_v)`.
    pub fn cells(&self) -> (usize, usize) {
        (self.grid.cells_u, self.grid.cells_v)
    }

    /// Index of the vertex closest to `(x, y)` lying within `tol`, if any.
    pub fn find_vertex(&mut self, x: f64, y: f64, tol: f64) -> Option<usize> {
        self.inspector.set_point(x, y);
        self.inspector.set_tolerance(tol);
        let (lo, hi) = self.grid.cell_range(x - tol, y - tol, x + tol, y + tol);
        for cj in lo.1..=hi.1 {
            for ci in lo.0..=hi.0 {
                for idx in self.grid.entries(ci, cj).to_vec() {
                    self.inspector.inspect(idx, &self.vertices[idx]);
                }
            }
        }
        self.inspector.best()
    }

    /// Indices of all vertices lying within `tol` of `(x, y)`.
    pub fn select(&mut self, x: f64, y: f64, tol: f64) -> Vec<usize> {
        self.inspector.set_point(x, y);
        self.inspector.set_tolerance(tol);
        let mut hits = Vec::new();
        let (lo, hi) = self.grid.cell_range(x - tol, y - tol, x + tol, y + tol);
        for cj in lo.1..=hi.1 {
            for ci in lo.0..=hi.0 {
                for idx in self.grid.entries(ci, cj).to_vec() {
                    if self.inspector.inspect(idx, &self.vertices[idx]) {
                        hits.push(idx);
                    }
                }
            }
        }
        hits
    }

    /// Adds a vertex unless one already exists within `tol`; returns the
    /// existing or newly created vertex index.
    pub fn add_vertex(&mut self, vertex: DelaunVertex, tol: f64) -> usize {
        if let Some(existing) = self.find_vertex(vertex.u, vertex.v, tol) {
            return existing;
        }
        let index = if let Some(slot) = self.del_nodes.pop() {
            self.vertices[slot] = vertex;
            slot
        } else {
            self.vertices.push(vertex);
            self.vertices.len() - 1
        };
        self.grid.push(vertex.u, vertex.v, index);
        index
    }

    /// Marks the vertex deleted and removes it from the grid; its slot is
    /// reused by the next [`VertexTool::add_vertex`].
    pub fn delete_vertex(&mut self, index: usize) {
        if let Some(v) = self.vertices.get_mut(index) {
            if v.is_deleted() {
                return;
            }
            let (u, vv) = (v.u, v.v);
            v.movability = Movability::Deleted;
            self.grid.remove(u, vv, index);
            self.del_nodes.push(index);
        }
    }

    /// Consumes the tool and returns the live vertices (deleted slots dropped).
    pub fn into_vertices(self) -> Vec<DelaunVertex> {
        self.vertices.into_iter().filter(|v| !v.is_deleted()).collect()
    }
}

// ---------------------------------------------------------------------------
// CircleInspector / CircleTool
// ---------------------------------------------------------------------------

/// Cell-query predicate for circles shot by a point.
///
/// Port of `BRepMesh_CircleInspector`: a circle is "shot" when the query point
/// lies on or near its circumference, i.e. `dist(point, center)^2 - radius^2`
/// does not exceed the squared tolerance.
#[derive(Debug, Clone)]
pub struct CircleInspector {
    sq_tolerance: f64,
    target_u: f64,
    target_v: f64,
    hits: Vec<usize>,
}

impl CircleInspector {
    /// Builds an inspector with the given linear tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self { sq_tolerance: tolerance * tolerance, target_u: 0.0, target_v: 0.0, hits: Vec::new() }
    }

    /// Sets the linear tolerance (squared internally).
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.sq_tolerance = tolerance * tolerance;
    }

    /// Sets the reference (bullet) point and clears the shot list.
    pub fn set_point(&mut self, u: f64, v: f64) {
        self.target_u = u;
        self.target_v = v;
        self.hits.clear();
    }

    /// Inspects one circle: `true` when the reference point shots it. Deleted
    /// circles (negative radius) are ignored.
    pub fn inspect(&mut self, index: usize, circle: &DelaunCircle) -> bool {
        if circle.is_deleted() {
            return false;
        }
        let du = circle.u - self.target_u;
        let dv = circle.v - self.target_v;
        let hit = du * du + dv * dv - circle.radius * circle.radius <= self.sq_tolerance;
        if hit {
            self.hits.push(index);
        }
        hit
    }

    /// Indices of the circles shot by the reference point.
    pub fn hits(&self) -> &[usize] {
        &self.hits
    }
}

/// Spatial index of circumcircles keyed by their centers.
///
/// Port of `BRepMesh_CircleTool`: each circle's bounding box (center ± radius,
/// clamped to the grid box) is hashed into the cell grid; queries return the
/// circles shot by a point or a circle whose center coincides with a point.
#[derive(Debug, Clone)]
pub struct CircleTool {
    grid: CellGrid,
    circles: Vec<DelaunCircle>,
    inspector: CircleInspector,
}

impl CircleTool {
    /// Builds a `cells_u × cells_v` grid covering `box2d`.
    pub fn new(box2d: Box2d, cells_u: usize, cells_v: usize, tolerance: f64) -> Self {
        Self {
            grid: CellGrid::new(box2d, cells_u, cells_v),
            circles: Vec::new(),
            inspector: CircleInspector::new(tolerance),
        }
    }

    /// The box covered by the grid.
    pub fn box2d(&self) -> Box2d {
        self.grid.box2d
    }

    /// True when no circle has been registered.
    pub fn is_empty(&self) -> bool {
        self.circles.is_empty()
    }

    /// Number of registered circles (deleted slots included).
    pub fn extent(&self) -> usize {
        self.circles.len()
    }

    /// Circle with the given index.
    pub fn circle(&self, index: usize) -> Option<&DelaunCircle> {
        self.circles.get(index)
    }

    /// Registers a circle and returns its index.
    pub fn add_circle(&mut self, circle: DelaunCircle) -> usize {
        let index = self.circles.len();
        self.circles.push(circle);
        self.key_in_grid(index, &circle);
        index
    }

    /// Binds a circle to the given (pre-allocated) index. Mirrors OCCT's
    /// `Bind`; intended for one-shot initialization, not re-keying a moved circle.
    pub fn bind_circle(&mut self, index: usize, circle: DelaunCircle) {
        if index >= self.circles.len() {
            self.circles.resize(index + 1, DelaunCircle::new(0.0, 0.0, -1.0));
        }
        self.circles[index] = circle;
        self.key_in_grid(index, &circle);
    }

    /// Index of the circle whose center lies within `tol` of `(u, v)`, if any.
    pub fn find_circle(&self, u: f64, v: f64, tol: f64) -> Option<usize> {
        let tol_sq = tol * tol;
        let (lo, hi) = self.grid.cell_range(u - tol, v - tol, u + tol, v + tol);
        let mut best: Option<usize> = None;
        let mut best_sq = f64::MAX;
        for cj in lo.1..=hi.1 {
            for ci in lo.0..=hi.0 {
                for &idx in self.grid.entries(ci, cj) {
                    let circle = &self.circles[idx];
                    if circle.is_deleted() {
                        continue;
                    }
                    let du = circle.u - u;
                    let dv = circle.v - v;
                    let sq = du * du + dv * dv;
                    if sq <= tol_sq && sq < best_sq {
                        best_sq = sq;
                        best = Some(idx);
                    }
                }
            }
        }
        best
    }

    /// Indices of the circles whose circumcircle is shot by `(u, v)`.
    ///
    /// A circle's bounding box may cover several cells, so the same index can
    /// be visited more than once; the result is deduplicated.
    pub fn select_circles(&mut self, u: f64, v: f64, tol: f64) -> Vec<usize> {
        self.inspector.set_point(u, v);
        self.inspector.set_tolerance(tol);
        let (lo, hi) = self.grid.cell_range(u - tol, v - tol, u + tol, v + tol);
        for cj in lo.1..=hi.1 {
            for ci in lo.0..=hi.0 {
                for idx in self.grid.entries(ci, cj).to_vec() {
                    self.inspector.inspect(idx, &self.circles[idx]);
                }
            }
        }
        let mut hits = Vec::new();
        for &idx in self.inspector.hits() {
            if !hits.contains(&idx) {
                hits.push(idx);
            }
        }
        hits
    }

    /// Marks the circle deleted (radius set negative, OCCT-style).
    pub fn delete_circle(&mut self, index: usize) {
        if let Some(circle) = self.circles.get_mut(index) {
            if circle.radius > 0.0 {
                circle.radius = -1.0;
            }
        }
    }

    /// Circumcircle of three points, or `None` when they are collinear or
    /// coincident. Returns `(u, v, radius)`. Port of `BRepMesh_CircleTool::MakeCircle`.
    ///
    /// DUPLICATE PORT: the faithful port of `BRepMesh_CircleTool::MakeCircle` is
    /// `delaun/constants.rs:150-183`. This copy is not reachable from the production
    /// mesh pipeline -- its only callers are the tests at `delaun_index.rs:673+`
    /// (`delaun/mod.rs:9-11` still says the sibling file is "under construction").
    /// Documented rather than aligned; do NOT mistake it for a faithful port.
    pub fn make_circle(
        p1: (f64, f64),
        p2: (f64, f64),
        p3: (f64, f64),
    ) -> Option<(f64, f64, f64)> {
        // `BRepMesh_CircleTool.cxx:80-81`: `aPrecision = Precision::PConfusion()`
        // = `Confusion() * 0.01` = 1e-9 (`Precision.hxx:334`),
        // `aSqPrecision = aPrecision * aPrecision`.
        const PRECISION: f64 = 1e-9;
        const SQ_PRECISION: f64 = PRECISION * PRECISION;
        // `BRepMesh_CircleTool.cxx:112`: `std::abs(aD) < gp::Resolution()`,
        // i.e. `RealSmall()` = `DBL_MIN` (`gp.hxx:60`,
        // `Standard_Real.hxx:132-135`). `delaun/constants.rs:167` already uses
        // the same constant for this test; the former port-local 1e-9 guard was
        // not OCCT-derived and rejected near-degenerate circumcircles.
        const DETERMINANT_GUARD: f64 = occt_core::precision::REAL_SMALL;

        let (x1, y1) = p1;
        let (x2, y2) = p2;
        let (x3, y3) = p3;

        let l1 = (x3 - x2, y2 - y3);
        if l1.0 * l1.0 + l1.1 * l1.1 < SQ_PRECISION {
            return None;
        }
        let l2 = (x1 - x3, y3 - y1);
        if l2.0 * l2.0 + l2.1 * l2.1 < SQ_PRECISION {
            return None;
        }
        let l3 = (x2 - x1, y1 - y2);
        if l3.0 * l3.0 + l3.1 * l3.1 < SQ_PRECISION {
            return None;
        }

        let d = 2.0 * (x1 * l1.1 + x2 * l2.1 + x3 * l3.1);
        if d.abs() < DETERMINANT_GUARD {
            return None;
        }
        let inv_d = 1.0 / d;
        let sq1 = x1 * x1 + y1 * y1;
        let sq2 = x2 * x2 + y2 * y2;
        let sq3 = x3 * x3 + y3 * y3;

        let cu = (sq1 * l1.1 + sq2 * l2.1 + sq3 * l3.1) * inv_d;
        let cv = (sq1 * l1.0 + sq2 * l2.0 + sq3 * l3.0) * inv_d;

        let r2 = ((x1 - cu).powi(2) + (y1 - cv).powi(2))
            .max((x2 - cu).powi(2) + (y2 - cv).powi(2))
            .max((x3 - cu).powi(2) + (y3 - cv).powi(2));
        let radius = r2.sqrt() + 2.0 * f64::EPSILON;
        Some((cu, cv, radius))
    }

    /// Hashes the circle's bounding box into the grid cells it covers.
    fn key_in_grid(&mut self, index: usize, circle: &DelaunCircle) {
        let b = self.grid.box2d;
        let min_u = (circle.u - circle.radius).max(b.min_u);
        let max_u = (circle.u + circle.radius).min(b.max_u);
        let min_v = (circle.v - circle.radius).max(b.min_v);
        let max_v = (circle.v + circle.radius).min(b.max_v);
        let (lo, hi) = self.grid.cell_range(min_u, min_v, max_u, max_v);
        for cj in lo.1..=hi.1 {
            for ci in lo.0..=hi.0 {
                let cell = self.grid.flat(ci, cj);
                self.grid.cells[cell].push(index);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_box() -> Box2d {
        Box2d::new(0.0, 0.0, 1.0, 1.0)
    }

    #[test]
    fn vertex_tolerance_find_merges_coincident() {
        let mut vt = VertexTool::new(unit_box(), 4, 4);
        let a = vt.add_vertex(DelaunVertex::new(0.2, 0.3), 1e-6);
        // Exact duplicate returns the existing index.
        let b = vt.add_vertex(DelaunVertex::new(0.2, 0.3), 1e-6);
        assert_eq!(a, b);
        // A point within the tolerance merges into the same index.
        let c = vt.add_vertex(DelaunVertex::new(0.2 + 1e-4, 0.3 + 1e-4), 1e-3);
        assert_eq!(c, a);
        // A distinct point is added.
        let d = vt.add_vertex(DelaunVertex::new(0.8, 0.9), 1e-6);
        assert_ne!(d, a);
        assert_eq!(vt.extent(), 2);

        // find_vertex locates by coordinate within tolerance.
        assert_eq!(vt.find_vertex(0.2, 0.3, 1e-6), Some(a));
        assert_eq!(vt.find_vertex(0.2, 0.3, 0.0), Some(a));
        assert_eq!(vt.find_vertex(0.5, 0.5, 1e-6), None);

        // Deleting frees the slot for reuse by the next addition.
        vt.delete_vertex(d);
        assert!(vt.vertex(d).unwrap().is_deleted());
        let e = vt.add_vertex(DelaunVertex::new(0.7, 0.7), 1e-6);
        assert_eq!(e, d);
    }

    #[test]
    fn vertex_grid_partitions_and_selects() {
        let mut vt = VertexTool::new(unit_box(), 4, 4);
        let a = vt.add_vertex(DelaunVertex::new(0.1, 0.1), 1e-6);
        let b = vt.add_vertex(DelaunVertex::new(0.9, 0.9), 1e-6);
        let c = vt.add_vertex(DelaunVertex::new(0.3, 0.7), 1e-6);

        // Each point lands in its own cell.
        assert_eq!(vt.get_cell(0.1, 0.1), (0, 0));
        assert_eq!(vt.get_cell(0.9, 0.9), (3, 3));
        assert_eq!(vt.get_cell(0.3, 0.7), (1, 2));
        // select_cell returns exactly the hashed entries.
        assert_eq!(vt.select_cell(0, 0), &[a]);
        assert_eq!(vt.select_cell(3, 3), &[b]);
        assert_eq!(vt.select_cell(1, 2), &[c]);
        // Tolerance-based select finds nearby vertices.
        assert_eq!(vt.select(0.1, 0.1, 1e-3), vec![a]);
        assert!(vt.select(0.5, 0.5, 0.05).is_empty());
        // Out-of-box coordinates clamp to the boundary cells.
        assert_eq!(vt.get_cell(-5.0, 2.0), (0, 3));
    }

    #[test]
    fn circle_tool_make_find_and_shoot() {
        let mut ct = CircleTool::new(unit_box(), 4, 4, 1e-6);

        // Right triangle (0,0)-(1,0)-(0,1): circumcenter (0.5,0.5), r = sqrt(0.5).
        let (u, v, r) = CircleTool::make_circle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)).unwrap();
        assert!((u - 0.5).abs() < 1e-9, "center u {u}");
        assert!((v - 0.5).abs() < 1e-9, "center v {v}");
        assert!((r - (0.5f64).sqrt()).abs() < 1e-9, "radius {r}");

        let big = ct.add_circle(DelaunCircle::new(u, v, r));
        let small = ct.add_circle(DelaunCircle::new(0.5, 0.5, 0.05));
        // Center coincidence find.
        assert_eq!(ct.find_circle(0.5, 0.5, 1e-6), Some(big));
        assert_eq!(ct.find_circle(0.1, 0.1, 1e-6), None);

        // Points on / inside the circumcircle are shot by it; the small circle
        // is not (its bbox does not reach these query points).
        assert_eq!(ct.select_circles(0.5, 0.0, 1e-6), vec![big]);
        assert_eq!(ct.select_circles(1.0, 0.5, 1e-6), vec![big]);
        assert_eq!(ct.select_circles(0.25, 0.25, 1e-6), vec![big]);
        // The center of both circles shots both.
        assert_eq!(ct.select_circles(0.5, 0.5, 1e-6), vec![big, small]);

        // Collinear points produce no circle.
        assert!(CircleTool::make_circle((0.0, 0.0), (0.5, 0.0), (1.0, 0.0)).is_none());
    }

    #[test]
    fn circle_grid_partitions_by_center_bbox() {
        let mut ct = CircleTool::new(unit_box(), 4, 4, 1e-6);
        let big = ct.add_circle(DelaunCircle::new(0.5, 0.5, 0.9));
        let small = ct.add_circle(DelaunCircle::new(0.1, 0.1, 0.05));

        // The big circle's bbox clamps to the whole grid box...
        assert_eq!(ct.box2d(), unit_box());
        // ...so it is shot from a corner, while the small circle is not.
        assert_eq!(ct.select_circles(0.95, 0.95, 1e-6), vec![big]);
        assert_eq!(ct.select_circles(0.1, 0.1, 1e-6), vec![big, small]);

        // Deleting a circle removes it from further queries.
        ct.delete_circle(big);
        assert!(ct.circle(big).unwrap().is_deleted());
        assert!(ct.select_circles(0.95, 0.95, 1e-6).is_empty());
    }
}
