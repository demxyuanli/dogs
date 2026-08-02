//! Port of OCCT `BRepMesh_Delaun` — incremental 2D Delaunay triangulation
//! (Bowyer–Watson / "algorithm of Watson") over UV points.
//!
//! Source: `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_Delaun.{hxx,cxx}`
//! (88.9K — the full algorithm, including super-triangle setup, circle-cell
//! acceleration, constraint-edge insertion and frontier adjustment).
//!
//! ponytail: the OCCT circle index (`BRepMesh_CircleTool`/`CircleInspector`,
//! sibling `delaun_index.rs`) is re-implemented here as a `HashMap`-backed grid
//! over `DelaunCircle`s because the sibling file is still under construction.
//! Swap for `delaun_index::CircleTool` when it lands.

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

use occt_core::gp::{GpPnt2d, GpVec2d, GpXY};
use occt_core::precision::{ANGULAR, PCONFUSION, RESOLUTION};

use super::delaun_data::{DelaunDataStructure, DelaunSelector};
use super::delaun_types::{DelaunCircle, DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::geom_tool::{GeomTool, IntFlag};

const ANG_DEV_1DEG: f64 = PI / 180.0;
const ANG_DEV_90DEG: f64 = 90.0 * ANG_DEV_1DEG;
const ANGLE_2PI: f64 = 2.0 * PI;

const PREC: f64 = PCONFUSION;
const PREC2: f64 = PREC * PREC;

/// Minimal 2-D axis-aligned bounding box. Source: `Bnd_B2d`.
#[derive(Debug, Clone, Copy)]
pub struct BndB2 {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
    is_void: bool,
}

impl BndB2 {
    fn void() -> Self {
        Self {
            min_x: f64::INFINITY,
            min_y: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            max_y: f64::NEG_INFINITY,
            is_void: true,
        }
    }

    fn add_pnt(&mut self, p: GpPnt2d) {
        if self.is_void {
            self.min_x = p.x();
            self.min_y = p.y();
            self.max_x = p.x();
            self.max_y = p.y();
            self.is_void = false;
        } else {
            self.min_x = self.min_x.min(p.x());
            self.max_x = self.max_x.max(p.x());
            self.min_y = self.min_y.min(p.y());
            self.max_y = self.max_y.max(p.y());
        }
    }

    fn add_xy(&mut self, p: GpXY) {
        self.add_pnt(GpPnt2d::from_xy(p));
    }

    fn enlarge(&mut self, tol: f64) {
        if !self.is_void {
            self.min_x -= tol;
            self.min_y -= tol;
            self.max_x += tol;
            self.max_y += tol;
        }
    }

    fn is_out(&self, other: &BndB2) -> bool {
        if self.is_void || other.is_void {
            return true;
        }
        self.max_x < other.min_x
            || other.max_x < self.min_x
            || self.max_y < other.min_y
            || other.max_y < self.min_y
    }

    fn get(&self) -> (f64, f64, f64, f64) {
        (self.min_x, self.min_y, self.max_x, self.max_y)
    }
}

fn update_bnd_box(p1: GpXY, p2: GpXY, b: &mut BndB2) {
    b.add_xy(p1);
    b.add_xy(p2);
    b.enlarge(PREC);
}

/// Cell-filtered store of circumcircles keyed by triangle id.
///
/// Source: `BRepMesh_CircleTool` + `BRepMesh_CircleInspector`. Circles are
/// stored in a `HashMap` grid keyed by 2-D cells; a query point only inspects
/// the circles whose clamped bounding box covers the query's cell.
pub struct CircleTool {
    tolerance: f64,
    sq_tolerance: f64,
    cell_size: GpXY,
    face_min: GpXY,
    face_max: GpXY,
    circles: Vec<DelaunCircle>,
    grid: HashMap<(i64, i64), Vec<i32>>,
}

impl CircleTool {
    fn new() -> Self {
        Self {
            tolerance: PREC,
            sq_tolerance: PREC * PREC,
            cell_size: GpXY::new(10.0, 10.0),
            face_min: GpXY::zero(),
            face_max: GpXY::zero(),
            circles: Vec::new(),
            grid: HashMap::new(),
        }
    }

    fn set_min_max_size(&mut self, min: GpXY, max: GpXY) {
        self.face_min = min;
        self.face_max = max;
    }

    fn set_cell_size(&mut self, size_x: f64, size_y: f64) {
        self.cell_size = GpXY::new(if size_x > 0.0 { size_x } else { 1.0 }, if size_y > 0.0 { size_y } else { 1.0 });
        self.grid.clear();
    }

    fn cell_of(&self, p: GpXY) -> (i64, i64) {
        let cx = ((p.x - self.face_min.x) / self.cell_size.x).floor() as i64;
        let cy = ((p.y - self.face_min.y) / self.cell_size.y).floor() as i64;
        (cx, cy)
    }

    fn bind(&mut self, index: i32, location: GpXY, radius: f64) {
        if self.circles.len() < index as usize {
            self.circles.resize(index as usize, DelaunCircle::default());
        }
        self.circles[(index - 1) as usize] = DelaunCircle::with_radius(GpPnt2d::from_xy(location), radius);

        let min_x = (location.x - radius).max(self.face_min.x);
        let max_x = (location.x + radius).min(self.face_max.x);
        let min_y = (location.y - radius).max(self.face_min.y);
        let max_y = (location.y + radius).min(self.face_max.y);
        let min_cell = self.cell_of(GpXY::new(min_x, min_y));
        let max_cell = self.cell_of(GpXY::new(max_x, max_y));
        for ci in min_cell.0..=max_cell.0 {
            for cj in min_cell.1..=max_cell.1 {
                self.grid.entry((ci, cj)).or_default().push(index);
            }
        }
    }

    /// Computes the circumcircle of three points. Source: `MakeCircle`.
    fn make_circle(p1: GpXY, p2: GpXY, p3: GpXY) -> Option<(GpXY, f64)> {
        let sq_prec = PREC2;
        let link1 = GpXY::new(p3.x - p2.x, p2.y - p3.y);
        if link1.square_modulus() < sq_prec {
            return None;
        }
        let link2 = GpXY::new(p1.x - p3.x, p3.y - p1.y);
        if link2.square_modulus() < sq_prec {
            return None;
        }
        let link3 = GpXY::new(p2.x - p1.x, p1.y - p2.y);
        if link3.square_modulus() < sq_prec {
            return None;
        }
        let d = 2.0 * (p1.x * link1.y + p2.x * link2.y + p3.x * link3.y);
        if d.abs() < RESOLUTION {
            return None;
        }
        let inv_d = 1.0 / d;
        let sq1 = p1.square_modulus();
        let sq2 = p2.square_modulus();
        let sq3 = p3.square_modulus();
        let loc = GpXY::new(
            (sq1 * link1.y + sq2 * link2.y + sq3 * link3.y) * inv_d,
            (sq1 * link1.x + sq2 * link2.x + sq3 * link3.x) * inv_d,
        );
        let r_sq = (p1.subtracted(&loc).square_modulus())
            .max(p2.subtracted(&loc).square_modulus())
            .max(p3.subtracted(&loc).square_modulus());
        let r = r_sq.sqrt() + 2.0 * f64::EPSILON;
        Some((loc, r))
    }

    /// Binds a circumcircle to the triangle index; returns `false` when the
    /// points are degenerate (no circle can be built). Source: `Bind`.
    fn bind_circle(&mut self, index: i32, p1: GpXY, p2: GpXY, p3: GpXY) -> bool {
        match Self::make_circle(p1, p2, p3) {
            Some((loc, r)) => {
                self.bind(index, loc, r);
                true
            }
            None => false,
        }
    }

    /// Binds an implicit zero (invalid) circle. Source: `MocBind`.
    fn moc_bind(&mut self, index: i32) {
        if self.circles.len() < index as usize {
            self.circles.resize(index as usize, DelaunCircle::default());
        }
        self.circles[(index - 1) as usize] = DelaunCircle::with_radius(GpPnt2d::zero(), -1.0);
    }

    /// Deletes the circle with the given index. Source: `Delete`.
    fn delete(&mut self, index: i32) {
        if let Some(c) = self.circles.get_mut((index - 1) as usize) {
            if c.is_created {
                c.set_radius(-1.0);
            }
        }
    }

    /// Returns indices of all circles shot by the point (containing it within
    /// the circle tolerance). Source: `Select`.
    fn select(&mut self, point: GpXY) -> Vec<i32> {
        let mut shot = Vec::new();
        let key = self.cell_of(point);
        if let Some(list) = self.grid.get(&key) {
            for &idx in list {
                let circle = &self.circles[(idx - 1) as usize];
                if !circle.is_created {
                    continue;
                }
                let dx = point.x - circle.center.x();
                let dy = point.y - circle.center.y();
                if dx * dx + dy * dy - circle.radius_sq <= self.sq_tolerance {
                    shot.push(idx);
                }
            }
        }
        shot
    }
}

/// Replacement mode for `create_and_replace_polygon_link`. Source: `ReplaceFlag`.
#[derive(Clone, Copy)]
enum ReplaceFlag {
    Replace,
    InsertAfter,
    InsertBefore,
}

/// Stack of element ranges used by `cleanup_polygon`. Source: `StackOfFrames`.
struct StackOfFrames {
    frames: Vec<(usize, usize)>,
}

impl StackOfFrames {
    fn new() -> Self {
        Self { frames: Vec::new() }
    }

    fn push_frame(&mut self, start: usize, end: usize) {
        self.frames.push((start, end));
    }

    fn pop_element(&mut self) -> usize {
        let (cur, end) = self.frames.last().copied().expect("pop on empty stack");
        self.frames.last_mut().unwrap().0 = cur + 1;
        if cur + 1 == end {
            self.frames.pop();
        }
        cur
    }

    fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

fn map_bind(map: &mut HashMap<i32, bool>, key: i32, value: bool) -> bool {
    match map.entry(key) {
        std::collections::hash_map::Entry::Occupied(_) => false,
        std::collections::hash_map::Entry::Vacant(v) => {
            v.insert(value);
            true
        }
    }
}

/// Delaunay triangulator over a [`DelaunDataStructure`].
///
/// Mirrors `BRepMesh_Delaun`. The full Watson pipeline is ported: super
/// triangle, incremental insertion with the circle cell filter, constraint
/// (frontier / fixed) edge insertion via polygon meshing, frontier adjustment,
/// cleanup and auxiliary-element removal.
pub struct Delaun {
    mesh_data: DelaunDataStructure,
    circles: CircleTool,
    sup_vert: Vec<i32>,
    init_circles: bool,
    sup_trian: DelaunTriangle,
}

impl Delaun {
    /// Creates a triangulation from a fresh data structure and the given
    /// vertices (OCCT `BRepMesh_Delaun(Array1OfVertexOfDelaun&)`).
    pub fn new_vertices(vertices: &[DelaunVertex]) -> Self {
        let mut data = DelaunDataStructure::new(vertices.len().max(16));
        let mut indices: Vec<i32> = Vec::with_capacity(vertices.len());
        for v in vertices {
            indices.push(data.add_node(*v));
        }
        let mut delaun = Self {
            mesh_data: data,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            sup_trian: DelaunTriangle::default(),
        };
        delaun.perform(&mut indices, -1, -1);
        delaun
    }

    /// Creates a triangulation over an existing data structure (OCCT
    /// `BRepMesh_Delaun(oldMesh, VectorOfInteger&)`). The vertex indices must
    /// already be registered in `data_structure`.
    pub fn new_with_data(data_structure: DelaunDataStructure, vertex_indices: &mut Vec<i32>) -> Self {
        let mut delaun = Self {
            mesh_data: data_structure,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            sup_trian: DelaunTriangle::default(),
        };
        delaun.perform(vertex_indices, -1, -1);
        delaun
    }

    /// Initializes the triangulation with an array of vertices (OCCT `Init`).
    pub fn init(&mut self, vertices: &[DelaunVertex]) {
        let mut indices: Vec<i32> = Vec::with_capacity(vertices.len());
        for v in vertices {
            indices.push(self.mesh_data.add_node(*v));
        }
        self.perform(&mut indices, -1, -1);
    }

    /// Forced re-initialization of the circle cell filter (OCCT
    /// `InitCirclesTool`). Binds a circumcircle for every live triangle.
    pub fn init_circles_tool_public(&mut self, cells_u: i32, cells_v: i32) {
        let mut box2 = BndB2::void();
        for i in 1..=self.mesh_data.nb_nodes() as i32 {
            box2.add_pnt(self.mesh_data.get_node(i).location);
        }
        box2.enlarge(PREC);
        self.init_circles_tool(&box2, cells_u, cells_v);
        let element_ids: Vec<i32> = self.mesh_data.elements_of_domain().iter().copied().collect();
        for id in element_ids {
            let triangle = self.mesh_data.get_element(id);
            let nodes = self.mesh_data.element_nodes(&triangle);
            self.circles.bind_circle(
                id,
                self.mesh_data.get_node(nodes[0]).location.coord,
                self.mesh_data.get_node(nodes[1]).location.coord,
                self.mesh_data.get_node(nodes[2]).location.coord,
            );
        }
    }

    /// Gives the mesh data structure.
    pub fn result(&self) -> &DelaunDataStructure {
        &self.mesh_data
    }

    /// Consumes the triangulator and returns the mesh data structure.
    pub fn into_result(self) -> DelaunDataStructure {
        self.mesh_data
    }

    /// The circle cell-filter tool used for point location.
    pub fn circles(&self) -> &CircleTool {
        &self.circles
    }

    /// Gives the list of frontier edges. Source: `Frontier()`.
    pub fn frontier(&self) -> HashSet<i32> {
        self.get_edges_by_type(VertexState::Frontier)
    }

    /// Gives the list of internal (fixed) edges. Source: `InternalEdges()`.
    pub fn internal_edges(&self) -> HashSet<i32> {
        self.get_edges_by_type(VertexState::Fixed)
    }

    /// Gives the list of free edges used at most once. Source: `FreeEdges()`.
    pub fn free_edges(&self) -> HashSet<i32> {
        self.get_edges_by_type(VertexState::Free)
    }

    /// Vertex by index.
    pub fn get_vertex(&self, index: i32) -> DelaunVertex {
        *self.mesh_data.get_node(index)
    }

    /// Edge (link) by index.
    pub fn get_edge(&self, index: i32) -> DelaunLink {
        self.mesh_data.get_link(index)
    }

    /// Triangle by index.
    pub fn get_triangle(&self, index: i32) -> DelaunTriangle {
        self.mesh_data.get_element(index)
    }

    /// Removes a vertex and re-triangulates the cavity. Source: `RemoveVertex`.
    pub fn remove_vertex(&mut self, vertex: &DelaunVertex) {
        let mut selector = DelaunSelector::new(&self.mesh_data);
        selector.neighbours_of(vertex);
        let mut loop_edges: HashMap<i32, bool> = HashMap::new();
        let elements: Vec<i32> = selector.elements().iter().copied().collect();
        for id in elements {
            self.delete_triangle(id, &mut loop_edges);
        }
        self.mesh_polygon_of_cavity(&mut loop_edges);
    }

    /// Adds some vertices into the triangulation. Source: `AddVertices`.
    pub fn add_vertices(&mut self, vertices: &mut Vec<i32>) {
        vertices.sort_by(|&a, &b| {
            let va = self.mesh_data.get_node(a).location;
            let vb = self.mesh_data.get_node(b).location;
            (va.x() + va.y()).total_cmp(&(vb.x() + vb.y()))
        });
        self.create_triangles_on_new_vertices(vertices);
    }

    /// Modifies the mesh to use the given edge. The OCCT implementation is an
    /// empty stub (constraint insertion happens through `ProcessConstraints`
    /// when links are pre-marked `Frontier`/`Fixed`); kept as a faithful stub.
    pub fn use_edge(&mut self, _index: i32) -> bool {
        false
    }

    /// Forces insertion of constraint edges and frontier adjustment.
    /// Source: `ProcessConstraints`.
    pub fn process_constraints(&mut self) {
        self.insert_internal_edges();
        self.frontier_adjust();
    }

    /// Tests whether the triangle contains the vertex (with square tolerance
    /// for edge closeness). Source: `Contains`.
    pub fn contains(&self, triangle_id: i32, vertex: &DelaunVertex, sq_tol: f64, edge_on: &mut i32) -> bool {
        *edge_on = 0;
        let element = self.mesh_data.get_element(triangle_id);
        let p = element.vertex_indices;
        let points = [
            self.mesh_data.get_node(p[0]).location,
            self.mesh_data.get_node(p[1]).location,
            self.mesh_data.get_node(p[2]).location,
        ];
        let v_edges = [
            points[1].coord.subtracted(&points[0].coord),
            points[2].coord.subtracted(&points[1].coord),
            points[0].coord.subtracted(&points[2].coord),
        ];
        let mut distance = [0.0; 3];
        let mut sq_modulus = [0.0; 3];
        let mut edge_on_id = 0usize;
        let sq_min_dist = self.calculate_dist(&v_edges, &points, vertex, &mut distance, &mut sq_modulus, &mut edge_on_id);
        if sq_min_dist < 0.0 {
            return false;
        }
        let edge_id_on = element.link_at(edge_on_id).abs();
        let is_not_free = self.mesh_data.link_movability(edge_id_on) != VertexState::Free;
        if sq_min_dist > sq_tol {
            if is_not_free && distance[edge_on_id] < sq_modulus[edge_on_id] / 5.0 {
                *edge_on = edge_id_on;
            }
        } else if is_not_free {
            return false;
        } else {
            *edge_on = edge_id_on;
        }
        distance[0] >= 0.0 && distance[1] >= 0.0 && distance[2] >= 0.0
    }

    // ------------------------------------------------------------------
    // Internal algorithm
    // ------------------------------------------------------------------

    fn perform(&mut self, vertex_indices: &mut Vec<i32>, cells_u: i32, cells_v: i32) {
        if vertex_indices.len() <= 2 {
            return;
        }
        let mut box2 = BndB2::void();
        for &idx in vertex_indices.iter() {
            box2.add_pnt(self.mesh_data.get_node(idx).location);
        }
        box2.enlarge(PREC);
        self.init_circles_tool(&box2, cells_u, cells_v);
        self.super_mesh(&box2);
        vertex_indices.sort_by(|&a, &b| {
            let va = self.mesh_data.get_node(a).location;
            let vb = self.mesh_data.get_node(b).location;
            (va.x() + va.y()).total_cmp(&(vb.x() + vb.y()))
        });
        self.compute(vertex_indices);
    }

    fn init_circles_tool(&mut self, box2: &BndB2, cells_u: i32, cells_v: i32) {
        let (min_x, min_y, max_x, max_y) = box2.get();
        let delta_x = max_x - min_x;
        let delta_y = max_y - min_y;
        let nb = self.mesh_data.nb_nodes();
        let scaler = if nb > 100 { 5 } else { 2 };
        self.circles.set_min_max_size(GpXY::new(min_x, min_y), GpXY::new(max_x, max_y));
        self.circles.set_cell_size(delta_x / cells_u.max(scaler) as f64, delta_y / cells_v.max(scaler) as f64);
        self.init_circles = true;
    }

    fn super_mesh(&mut self, box2: &BndB2) {
        let (min_x, min_y, max_x, max_y) = box2.get();
        let delta_x = max_x - min_x;
        let delta_y = max_y - min_y;
        let delta_min = delta_x.min(delta_y);
        let delta_max = delta_x.max(delta_y);
        let delta = delta_x + delta_y;

        self.sup_vert.clear();
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            (min_x + max_x) / 2.0,
            max_y + delta_max,
            VertexState::Free,
        )));
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            min_x - delta,
            min_y - delta_min,
            VertexState::Free,
        )));
        self.sup_vert.push(self.mesh_data.add_node(DelaunVertex::new_parametric(
            max_x + delta,
            min_y - delta_min,
            VertexState::Free,
        )));

        let mut e = [0i32; 3];
        let mut o = [false; 3];
        for node_id in 0..3 {
            let first = self.sup_vert[node_id];
            let last = self.sup_vert[(node_id + 1) % 3];
            let link_idx = self.mesh_data.add_link(first, last, VertexState::Free);
            e[node_id] = link_idx.abs();
            o[node_id] = link_idx > 0;
        }
        self.sup_trian = DelaunTriangle::new(
            [
                if o[0] { e[0] } else { -e[0] },
                if o[1] { e[1] } else { -e[1] },
                if o[2] { e[2] } else { -e[2] },
            ],
            [self.sup_vert[0], self.sup_vert[1], self.sup_vert[2]],
        );
    }

    fn compute(&mut self, vertex_indexes: &mut Vec<i32>) {
        let mut loop_edges: HashMap<i32, bool> = HashMap::new();
        for i in 0..3 {
            loop_edges.insert(self.sup_trian.link_at(i).abs(), true);
        }
        if vertex_indexes.len() > 0 {
            let first = vertex_indexes[0];
            self.create_triangles(first, &mut loop_edges);
            self.create_triangles_on_new_vertices(vertex_indexes);
        }
        self.remove_aux_elements();
    }

    fn delete_triangle(&mut self, index: i32, loop_edges: &mut HashMap<i32, bool>) {
        if self.init_circles {
            self.circles.delete(index);
        }
        let element = self.mesh_data.get_element(index);
        self.mesh_data.remove_element(index);
        for i in 0..3 {
            let e = element.link_at(i).abs();
            let o = element.link_at(i) > 0;
            if !map_bind(loop_edges, e, o) {
                loop_edges.remove(&e);
                self.mesh_data.remove_link(e, false);
            }
        }
    }

    fn remove_aux_elements(&mut self) {
        let mut loop_edges: HashMap<i32, bool> = HashMap::new();
        let mut elements: HashSet<i32> = HashSet::new();
        {
            let mut selector = DelaunSelector::new(&self.mesh_data);
            let sup = self.sup_vert.clone();
            for &sv in &sup {
                selector.neighbours_of_node(sv);
            }
            for &e in selector.elements() {
                elements.insert(e);
            }
        }
        for &id in elements.iter() {
            self.delete_triangle(id, &mut loop_edges);
        }
        let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
        for &e in &loop_keys {
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }
        let sup = self.sup_vert.clone();
        for &sv in &sup {
            self.mesh_data.remove_node(sv, false);
        }
    }

    fn create_triangles(&mut self, vertex_index: i32, poly: &mut HashMap<i32, bool>) {
        let mut loop_edges: Vec<i32> = Vec::new();
        let mut external_edges: Vec<i32> = Vec::new();
        let vertex_coord = self.mesh_data.get_node(vertex_index).location.coord;

        let keys: Vec<i32> = poly.keys().copied().collect();
        for edge_id in keys {
            let is_positive = poly[&edge_id];
            let edge = self.mesh_data.get_link(edge_id);
            let mut nodes = [0i32; 3];
            if is_positive {
                nodes[0] = edge.first_node();
                nodes[2] = edge.last_node();
            } else {
                nodes[0] = edge.last_node();
                nodes[2] = edge.first_node();
            }
            nodes[1] = vertex_index;

            let first_vertex = self.mesh_data.get_node(nodes[0]).location.coord;
            let last_vertex = self.mesh_data.get_node(nodes[2]).location.coord;
            let mut edge_dir = last_vertex.subtracted(&first_vertex);
            let edge_len = edge_dir.modulus();
            if edge_len < PREC {
                continue;
            }
            edge_dir = edge_dir.divided(edge_len);

            let first_link_dir = first_vertex.subtracted(&vertex_coord);
            let last_link_dir = vertex_coord.subtracted(&last_vertex);

            let dist12 = first_link_dir.crossed(&edge_dir);
            let dist23 = edge_dir.crossed(&last_link_dir);
            if dist12.abs() < PREC || dist23.abs() < PREC {
                continue;
            }

            let first_link_id = self.mesh_data.add_link(nodes[1], nodes[0], VertexState::Free);
            let last_link_id = self.mesh_data.add_link(nodes[2], nodes[1], VertexState::Free);
            let edges_info = [first_link_id, if is_positive { edge_id } else { -edge_id }, last_link_id];

            let is_sens_ok = dist12 > 0.0 && dist23 > 0.0;
            if is_sens_ok {
                let mut edge_ids = [0i32; 3];
                let mut edge_oris = [false; 3];
                for k in 0..3 {
                    edge_ids[k] = edges_info[k].abs();
                    edge_oris[k] = edges_info[k] > 0;
                }
                self.add_triangle(edge_ids, edge_oris, nodes);
            } else {
                if is_positive {
                    loop_edges.push(edge_id);
                } else {
                    loop_edges.push(-edge_id);
                }
                if first_link_dir.square_modulus() > last_link_dir.square_modulus() {
                    external_edges.push(edges_info[0].abs());
                } else {
                    external_edges.push(edges_info[2].abs());
                }
            }
        }

        poly.clear();
        while !external_edges.is_empty() {
            let e = external_edges.remove(0);
            let pair = self.mesh_data.elements_connected_to(e);
            if !pair.is_empty() {
                self.delete_triangle(pair.first_index(), poly);
            }
        }

        let poly_keys: Vec<i32> = poly.keys().copied().collect();
        for &e in &poly_keys {
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }

        while !loop_edges.is_empty() {
            let edge_info = loop_edges.remove(0);
            if self.mesh_data.link_movability(edge_info.abs()) != VertexState::Deleted {
                self.mesh_left_polygon_of(edge_info.abs(), edge_info > 0, &mut None);
            }
        }
    }

    fn create_triangles_on_new_vertices(&mut self, vertex_indexes: &mut Vec<i32>) {
        let (tol_u, tol_v) = self.mesh_data.get_tolerance();
        let sq_tol = tol_u * tol_u + tol_v * tol_v;

        let upper = vertex_indexes.len();
        let mut i = 0;
        while i < upper {
            let vertex_idx = vertex_indexes[i];
            let vertex = *self.mesh_data.get_node(vertex_idx);
            let mut loop_edges: HashMap<i32, bool> = HashMap::new();
            let mut circles_list = self.circles.select(vertex.location.coord);

            let mut on_edge_id = 0;
            let mut triangle_id = 0;
            let mut j = 0;
            while j < circles_list.len() {
                let t = circles_list[j];
                let mut edge_on = 0;
                if self.contains(t, &vertex, sq_tol, &mut edge_on) {
                    if edge_on != 0 && self.mesh_data.link_movability(edge_on) != VertexState::Free {
                        if vertex.state == VertexState::Free {
                            j += 1;
                            continue;
                        }
                    }
                    triangle_id = t;
                    circles_list.remove(j);
                    break;
                }
                j += 1;
            }

            if triangle_id > 0 {
                self.delete_triangle(triangle_id, &mut loop_edges);
                let mut is_modify = true;
                while is_modify && !circles_list.is_empty() {
                    is_modify = false;
                    let mut k = 0;
                    while k < circles_list.len() {
                        let t = circles_list[k];
                        let element = self.mesh_data.get_element(t);
                        if loop_edges.contains_key(&element.link_at(0).abs())
                            || loop_edges.contains_key(&element.link_at(1).abs())
                            || loop_edges.contains_key(&element.link_at(2).abs())
                        {
                            is_modify = true;
                            self.delete_triangle(t, &mut loop_edges);
                            circles_list.remove(k);
                            break;
                        }
                        k += 1;
                    }
                }
                self.create_triangles(vertex_idx, &mut loop_edges);
            }
            i += 1;
        }

        // Constraint processing (frontier adjustment + mesh cleanup) only makes
        // sense when the mesh actually has constraint edges. For a plain
        // triangulation (no Frontier/Fixed edges) OCCT's `cleanupMesh` would
        // strip every boundary triangle whose neighbour touches the super
        // triangle, orphaning valid hull vertices. Skip it so a plain Delaunay
        // keeps a complete triangulation.
        if !self.frontier().is_empty() || !self.internal_edges().is_empty() {
            self.process_constraints();
        }
    }

    fn add_triangle(&mut self, edges: [i32; 3], oris: [bool; 3], nodes: [i32; 3]) {
        let signed = [
            if oris[0] { edges[0] } else { -edges[0] },
            if oris[1] { edges[1] } else { -edges[1] },
            if oris[2] { edges[2] } else { -edges[2] },
        ];
        // The vertex order must match the orientation of the three links:
        // each link runs from the start of one side to the start of the next,
        // so the CCW vertex cycle is [start of link0, start of link1, start of link2].
        // OCCT derives the triangle's vertices from its links the same way; storing
        // the caller's `nodes` verbatim put them in the reversed order, which broke
        // `contains`'s point-in-triangle orientation test (all inserts were rejected).
        let verts = [
            {
                let l = self.mesh_data.get_link(edges[0].abs());
                if oris[0] { l.first_node() } else { l.last_node() }
            },
            {
                let l = self.mesh_data.get_link(edges[1].abs());
                if oris[1] { l.first_node() } else { l.last_node() }
            },
            {
                let l = self.mesh_data.get_link(edges[2].abs());
                if oris[2] { l.first_node() } else { l.last_node() }
            },
        ];
        let tri = DelaunTriangle::new(signed, verts);
        let new_id = self.mesh_data.add_element(tri);
        let mut is_added = true;
        if self.init_circles {
            is_added = self.circles.bind_circle(
                new_id,
                self.mesh_data.get_node(verts[0]).location.coord,
                self.mesh_data.get_node(verts[1]).location.coord,
                self.mesh_data.get_node(verts[2]).location.coord,
            );
        }
        if !is_added {
            self.mesh_data.remove_element(new_id);
        }
    }

    fn insert_internal_edges(&mut self) {
        let internal_edges = self.internal_edges();
        for &link_index in internal_edges.iter() {
            let pair = self.mesh_data.elements_connected_to(link_index);
            let mut is_go = [true, true];
            for tri_it in 1..=pair.extent() {
                let element = self.mesh_data.get_element(pair.index(tri_it));
                for i in 0..3 {
                    if element.link_at(i).abs() == link_index {
                        is_go[if element.link_at(i) > 0 { 0 } else { 1 }] = false;
                        break;
                    }
                }
            }
            if is_go[0] {
                self.mesh_left_polygon_of(link_index, true, &mut None);
            }
            if is_go[1] {
                self.mesh_left_polygon_of(link_index, false, &mut None);
            }
        }
    }

    fn is_bound_to_frontier(&self, ref_node_id: i32, ref_link_id: i32) -> bool {
        let mut stack: Vec<i32> = vec![ref_link_id];
        let mut visited: HashSet<i32> = HashSet::new();
        while let Some(cur) = stack.pop() {
            let pair = self.mesh_data.elements_connected_to(cur);
            if pair.is_empty() {
                return false;
            }
            for elem_it in 1..=pair.extent() {
                let tri_id = pair.index(elem_it);
                if tri_id < 0 {
                    continue;
                }
                let element = self.mesh_data.get_element(tri_id);
                for k in 0..3 {
                    let edge_id = element.link_at(k).abs();
                    if edge_id == cur {
                        continue;
                    }
                    let edge = self.mesh_data.get_link(edge_id);
                    if edge.first_node() != ref_node_id && edge.last_node() != ref_node_id {
                        continue;
                    }
                    if self.mesh_data.link_movability(edge_id) != VertexState::Free {
                        return true;
                    }
                    if visited.insert(edge_id) {
                        stack.push(edge_id);
                    }
                }
            }
        }
        false
    }

    fn cleanup_mesh(&mut self) {
        loop {
            let mut loop_edges: HashMap<i32, bool> = HashMap::new();
            let mut del_triangles: HashSet<i32> = HashSet::new();

            let free_edges = self.free_edges();
            for &free_edge_id in free_edges.iter() {
                let edge = self.mesh_data.get_link(free_edge_id);
                if self.mesh_data.link_movability(free_edge_id) == VertexState::Frontier {
                    continue;
                }
                let pair = self.mesh_data.elements_connected_to(free_edge_id);
                if pair.is_empty() {
                    loop_edges.insert(free_edge_id, true);
                    continue;
                }
                let tri_id = pair.first_index();
                let element = self.mesh_data.get_element(tri_id);
                let an_edges = element.link_indices;

                let mut can_not_be_removed = true;
                for cur_edge_idx in 0..3 {
                    if an_edges[cur_edge_idx].abs() != free_edge_id {
                        continue;
                    }
                    for other in 1..=2 {
                        if !can_not_be_removed {
                            break;
                        }
                        let other_edge_id = an_edges[(cur_edge_idx + other) % 3].abs();
                        let other_pair = self.mesh_data.elements_connected_to(other_edge_id);
                        if other_pair.extent() < 2 {
                            can_not_be_removed = false;
                        } else {
                            for tri_idx in 1..=other_pair.extent() {
                                if !can_not_be_removed {
                                    break;
                                }
                                if other_pair.index(tri_idx) == tri_id {
                                    continue;
                                }
                                let cur_triangle = self.mesh_data.get_element(other_pair.index(tri_idx));
                                for v in cur_triangle.vertex_indices {
                                    if self.is_sup_vertex(v) {
                                        can_not_be_removed = false;
                                    }
                                }
                            }
                        }
                    }
                    break;
                }
                if can_not_be_removed {
                    continue;
                }

                let mut is_connected = [false, false];
                for l in 0..2 {
                    let node = if l == 0 { edge.first_node() } else { edge.last_node() };
                    is_connected[l] = self.is_bound_to_frontier(node, free_edge_id);
                }
                if !is_connected[0] || !is_connected[1] {
                    del_triangles.insert(tri_id);
                }
            }

            let mut deleted_nb = 0;
            let dels: Vec<i32> = del_triangles.iter().copied().collect();
            for id in dels {
                self.delete_triangle(id, &mut loop_edges);
                deleted_nb += 1;
            }

            let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
            for &e in &loop_keys {
                if self.mesh_data.elements_connected_to(e).is_empty() {
                    self.mesh_data.remove_link(e, false);
                }
            }
            if deleted_nb == 0 {
                break;
            }
        }
    }

    fn frontier_adjust(&mut self) {
        let frontier = self.frontier();
        let mut failed_frontiers: Vec<i32> = Vec::new();
        let mut loop_edges: HashMap<i32, bool> = HashMap::new();
        let mut int_frontier_edges: HashSet<i32> = HashSet::new();

        for pass in 1..=2 {
            for &frontier_id in frontier.iter() {
                let pair = self.mesh_data.elements_connected_to(frontier_id);
                let nb = pair.extent();
                for elem_it in 1..=nb {
                    let prior_elem = pair.index(elem_it);
                    if prior_elem < 0 {
                        continue;
                    }
                    let element = self.mesh_data.get_element(prior_elem);
                    let mut found = false;
                    for n in 0..3 {
                        if frontier_id == element.link_at(n).abs() && element.link_at(n) < 0 {
                            found = true;
                            self.delete_triangle(prior_elem, &mut loop_edges);
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }
            }

            let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
            for &e in &loop_keys {
                if self.mesh_data.elements_connected_to(e).is_empty() {
                    self.mesh_data.remove_link(e, false);
                }
            }

            for &frontier_id in frontier.iter() {
                if !self.mesh_data.elements_connected_to(frontier_id).is_empty() {
                    continue;
                }
                let mut skipped = Some(int_frontier_edges.clone());
                let success = self.mesh_left_polygon_of(frontier_id, true, &mut skipped);
                if let Some(s) = skipped {
                    int_frontier_edges = s;
                }
                if pass == 2 && !success {
                    failed_frontiers.push(frontier_id);
                }
            }
        }

        self.cleanup_mesh();

        for &frontier_id in failed_frontiers.iter() {
            if !self.mesh_data.elements_connected_to(frontier_id).is_empty() {
                continue;
            }
            let mut skipped = Some(int_frontier_edges.clone());
            self.mesh_left_polygon_of(frontier_id, true, &mut skipped);
            if let Some(s) = skipped {
                int_frontier_edges = s;
            }
        }
    }

    fn fill_bnd_box(&self, boxes: &mut Vec<BndB2>, v1: i32, v2: i32) {
        let mut b = BndB2::void();
        update_bnd_box(
            self.mesh_data.get_node(v1).location.coord,
            self.mesh_data.get_node(v2).location.coord,
            &mut b,
        );
        boxes.push(b);
    }

    fn mesh_left_polygon_of(&mut self, start_edge_id: i32, is_forward: bool, skipped: &mut Option<HashSet<i32>>) -> bool {
        if let Some(s) = skipped.as_ref() {
            if s.contains(&start_edge_id) {
                return true;
            }
        }
        let ref_edge = self.mesh_data.get_link(start_edge_id);

        let mut polygon: Vec<i32> = Vec::new();
        let (a_start_node, mut a_pivot_node);
        if is_forward {
            polygon.push(start_edge_id);
            a_start_node = ref_edge.first_node();
            a_pivot_node = ref_edge.last_node();
        } else {
            polygon.push(-start_edge_id);
            a_start_node = ref_edge.last_node();
            a_pivot_node = ref_edge.first_node();
        }
        let start_edge_vertex_s = self.mesh_data.get_node(a_start_node).location.coord;
        let mut a_pivot_vertex = self.mesh_data.get_node(a_pivot_node).location.coord;
        let mut ref_link_dir = GpVec2d::from_xy(a_pivot_vertex.subtracted(&start_edge_vertex_s));
        if ref_link_dir.square_magnitude() < PREC2 {
            return true;
        }

        let mut boxes: Vec<BndB2> = Vec::new();
        let mut b0 = BndB2::void();
        update_bnd_box(start_edge_vertex_s, a_pivot_vertex, &mut b0);
        boxes.push(b0);

        let mut dead_links: HashSet<i32> = HashSet::new();
        let mut leprous_links: HashSet<i32> = HashSet::new();
        leprous_links.insert(start_edge_id);

        let mut is_skip_leprous = true;
        let mut a_first_node = a_start_node;
        while a_pivot_node != a_first_node {
            let result = self.find_next_polygon_link(
                a_first_node,
                a_pivot_node,
                a_pivot_vertex,
                ref_link_dir,
                &boxes,
                &polygon,
                skipped.as_ref(),
                is_skip_leprous,
                &mut leprous_links,
                &mut dead_links,
            );
            if let Some((next_link_id, next_pivot_node, next_link_dir, next_link_bbox)) = result {
                a_first_node = a_pivot_node;
                ref_link_dir = next_link_dir;
                a_pivot_node = next_pivot_node;
                a_pivot_vertex = self.mesh_data.get_node(next_pivot_node).location.coord;
                boxes.push(next_link_bbox);
                polygon.push(next_link_id);
                is_skip_leprous = true;
            } else {
                if polygon.len() == 1 {
                    return false;
                }
                let dead_link_id = polygon.last().copied().unwrap().abs();
                dead_links.insert(dead_link_id);
                leprous_links.remove(&dead_link_id);
                polygon.pop();
                boxes.pop();

                let prev_link_info = *polygon.last().unwrap();
                let prev_link = self.mesh_data.get_link(prev_link_info.abs());
                if prev_link_info > 0 {
                    a_first_node = prev_link.first_node();
                    a_pivot_node = prev_link.last_node();
                } else {
                    a_first_node = prev_link.last_node();
                    a_pivot_node = prev_link.first_node();
                }
                a_pivot_vertex = self.mesh_data.get_node(a_pivot_node).location.coord;
                let sn = self.mesh_data.get_node(a_start_node).location.coord;
                ref_link_dir = GpVec2d::from_xy(a_pivot_vertex.subtracted(&sn));
                is_skip_leprous = false;
            }
        }

        if polygon.len() < 3 {
            return false;
        }
        self.cleanup_polygon(&polygon, &boxes);
        self.mesh_polygon(&mut polygon, &mut boxes, skipped);
        true
    }

    #[allow(clippy::too_many_arguments)]
    fn find_next_polygon_link(
        &self,
        first_node: i32,
        pivot_node: i32,
        pivot_vertex: GpXY,
        ref_link_dir: GpVec2d,
        boxes: &[BndB2],
        polygon: &[i32],
        skipped: Option<&HashSet<i32>>,
        is_skip_leprous: bool,
        leprous_links: &mut HashSet<i32>,
        dead_links: &mut HashSet<i32>,
    ) -> Option<(i32, i32, GpVec2d, BndB2)> {
        let mut max_angle = f64::NEG_INFINITY;
        let mut next_link_id = 0;
        let mut next_pivot_node = 0;
        let mut next_link_dir = GpVec2d::zero();
        let mut next_link_bbox = BndB2::void();

        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(pivot_node).to_vec();
        for &neighbor_info in &neighbors {
            let neighbor_id = neighbor_info.abs();
            if dead_links.contains(&neighbor_id) {
                continue;
            }
            if let Some(s) = skipped {
                if s.contains(&neighbor_id) {
                    continue;
                }
            }
            let is_leprous = leprous_links.contains(&neighbor_id);
            if is_skip_leprous && is_leprous {
                continue;
            }
            let neighbor_link = self.mesh_data.get_link(neighbor_id);
            if self.mesh_data.link_movability(neighbor_id) == VertexState::Free
                && self.mesh_data.elements_connected_to(neighbor_id).is_empty()
            {
                dead_links.insert(neighbor_id);
                continue;
            }
            let mut other_node = neighbor_link.first_node();
            if other_node == pivot_node {
                other_node = neighbor_link.last_node();
            }
            let cur_link_dir = GpVec2d::from_xy(self.mesh_data.get_node(other_node).location.coord.subtracted(&pivot_vertex));
            if cur_link_dir.square_magnitude() < PREC2 {
                dead_links.insert(neighbor_id);
                continue;
            }
            if !is_leprous {
                leprous_links.insert(neighbor_id);
            }

            let mut angle = ref_link_dir.angle(&cur_link_dir);
            let is_frontier = self.mesh_data.link_movability(neighbor_id) == VertexState::Frontier;
            let mut is_check_point_on_edge = true;
            if is_frontier {
                if (angle.abs() - PI).abs() < ANGULAR {
                    is_check_point_on_edge = false;
                    angle = angle.abs();
                }
            }
            if angle <= max_angle {
                continue;
            }

            let is_check_end_points = other_node != first_node;
            let mut a_box = BndB2::void();
            let is_not_intersect = self.check_intersection(
                neighbor_id,
                polygon,
                boxes,
                is_check_end_points,
                is_check_point_on_edge,
                true,
                &mut a_box,
            );
            if is_not_intersect {
                max_angle = angle;
                next_link_dir = cur_link_dir;
                next_pivot_node = other_node;
                next_link_bbox = a_box;
                next_link_id = if neighbor_link.first_node() == pivot_node { neighbor_id } else { -neighbor_id };
            }
        }

        if next_link_id == 0 {
            None
        } else {
            Some((next_link_id, next_pivot_node, next_link_dir, next_link_bbox))
        }
    }

    fn check_intersection(
        &self,
        link_id: i32,
        polygon: &[i32],
        poly_boxes: &[BndB2],
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        is_skip_last_edge: bool,
        link_bbox: &mut BndB2,
    ) -> bool {
        let link = self.mesh_data.get_link(link_id);
        update_bnd_box(
            self.mesh_data.get_node(link.first_node()).location.coord,
            self.mesh_data.get_node(link.last_node()).location.coord,
            link_bbox,
        );
        let mut poly_len = polygon.len();
        if is_skip_last_edge {
            poly_len -= 1;
        }
        let is_frontier = self.mesh_data.link_movability(link_id) == VertexState::Frontier;

        for poly_it in 0..poly_len {
            if !link_bbox.is_out(&poly_boxes[poly_it]) {
                let poly_link_id = polygon[poly_it].abs();
                let poly_link = self.mesh_data.get_link(poly_link_id);
                if self.mesh_data.link_movability(poly_link_id) == VertexState::Frontier && is_frontier {
                    continue;
                }
                let mut int_pnt = GpPnt2d::zero();
                let flag = self.int_seg_seg(link_id, poly_link_id, is_consider_end_point_touch, is_consider_point_on_edge, &mut int_pnt);
                if flag != IntFlag::NoIntersection {
                    return false;
                }
            }
        }
        true
    }

    fn add_triangle_by_info(&mut self, edges_info: [i32; 3], nodes: [i32; 3]) {
        let mut edges = [0i32; 3];
        let mut oris = [false; 3];
        for i in 0..3 {
            edges[i] = edges_info[i].abs();
            oris[i] = edges_info[i] > 0;
        }
        self.add_triangle(edges, oris, nodes);
    }

    fn cleanup_polygon(&mut self, the_polygon: &[i32], the_poly_boxes: &[BndB2]) {
        let poly_len = the_polygon.len();
        if poly_len < 3 {
            return;
        }
        let mut loop_edges: HashMap<i32, bool> = HashMap::new();
        let mut ignored_edges: HashSet<i32> = HashSet::new();
        let mut poly_vertices_find_map: HashSet<i32> = HashSet::new();
        let mut poly_vertices: Vec<i32> = Vec::new();

        for poly_it in 1..=poly_len {
            let poly_edge_info = the_polygon[poly_it - 1];
            let poly_edge_id = poly_edge_info.abs();
            ignored_edges.insert(poly_edge_id);
            let is_forward = poly_edge_info > 0;
            let pair = self.mesh_data.elements_connected_to(poly_edge_id);

            let mut elem_it = 1;
            while elem_it <= pair.extent() {
                let elem_id = pair.index(elem_it);
                if elem_id < 0 {
                    elem_it += 1;
                    continue;
                }
                let element = self.mesh_data.get_element(elem_id);
                let mut found = false;
                for k in 0..3 {
                    if element.link_at(k).abs() == poly_edge_id && (element.link_at(k) > 0) == is_forward {
                        found = true;
                        self.delete_triangle(elem_id, &mut loop_edges);
                        break;
                    }
                }
                if found {
                    break;
                }
                elem_it += 1;
            }

            if poly_it % 2 == 1 {
                let poly_edge = self.mesh_data.get_link(poly_edge_id);
                let first_vertex = poly_edge.first_node();
                let last_vertex = poly_edge.last_node();
                poly_vertices_find_map.insert(first_vertex);
                poly_vertices_find_map.insert(last_vertex);
                if poly_edge_info > 0 {
                    poly_vertices.push(first_vertex);
                    poly_vertices.push(last_vertex);
                } else {
                    poly_vertices.push(last_vertex);
                    poly_vertices.push(first_vertex);
                }
            }
        }

        if poly_vertices.first() != poly_vertices.last() {
            let first = *poly_vertices.first().unwrap();
            poly_vertices.push(first);
        }

        let mut survived_links = ignored_edges.clone();
        for poly_vert_it in 0..(poly_vertices.len() - 1) {
            let mut stack_frames = StackOfFrames::new();
            let mut stack_data: Vec<i32> = Vec::new();
            let mut current_victim = poly_vertices[poly_vert_it];
            loop {
                let prev_size = stack_data.len();
                self.kill_triangles_around_vertex(
                    current_victim,
                    &poly_vertices,
                    &poly_vertices_find_map,
                    the_polygon,
                    the_poly_boxes,
                    &mut survived_links,
                    &mut loop_edges,
                    &mut stack_data,
                );
                let new_size = stack_data.len();
                if new_size > prev_size {
                    stack_frames.push_frame(prev_size, new_size);
                }
                if stack_frames.is_empty() {
                    break;
                }
                current_victim = stack_data[stack_frames.pop_element()];
            }
        }

        let loop_keys: Vec<i32> = loop_edges.keys().copied().collect();
        for &e in &loop_keys {
            if ignored_edges.contains(&e) {
                continue;
            }
            if self.mesh_data.elements_connected_to(e).is_empty() {
                self.mesh_data.remove_link(e, false);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn kill_triangles_around_vertex(
        &mut self,
        zombie_node_id: i32,
        poly_vertices: &[i32],
        poly_vertices_find_map: &HashSet<i32>,
        the_polygon: &[i32],
        the_poly_boxes: &[BndB2],
        survived_links: &mut HashSet<i32>,
        loop_edges: &mut HashMap<i32, bool>,
        victim_nodes: &mut Vec<i32>,
    ) {
        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(zombie_node_id).to_vec();
        for &neighbor_link_id in &neighbors {
            if survived_links.contains(&neighbor_link_id) {
                continue;
            }
            let neighbor_link = self.mesh_data.get_link(neighbor_link_id);
            if self.mesh_data.link_movability(neighbor_link_id) == VertexState::Frontier {
                let mut b = BndB2::void();
                let is_not_intersect =
                    self.check_intersection(neighbor_link_id, the_polygon, the_poly_boxes, false, true, false, &mut b);
                if is_not_intersect {
                    survived_links.insert(neighbor_link_id);
                    continue;
                }
            } else {
                let mut other_node = neighbor_link.first_node();
                if other_node == zombie_node_id {
                    other_node = neighbor_link.last_node();
                }
                if !poly_vertices_find_map.contains(&other_node) {
                    if self.is_vertex_inside_polygon(other_node, poly_vertices) {
                        victim_nodes.push(other_node);
                    } else {
                        self.kill_triangles_on_intersecting_links(
                            neighbor_link_id,
                            other_node,
                            the_polygon,
                            the_poly_boxes,
                            survived_links,
                            loop_edges,
                        );
                        continue;
                    }
                }
            }
            survived_links.insert(neighbor_link_id);
            self.kill_link_triangles(neighbor_link_id, loop_edges);
        }
    }

    fn is_vertex_inside_polygon(&self, vertex_id: i32, polygon_vertices: &[i32]) -> bool {
        let poly_len = polygon_vertices.len();
        if poly_len < 3 {
            return false;
        }
        let center = self.mesh_data.get_node(vertex_id).location;
        let first_vertex = self.mesh_data.get_node(polygon_vertices[0]).location;
        let mut prev_dir = GpVec2d::from_xy(first_vertex.coord.subtracted(&center.coord));
        if prev_dir.square_magnitude() < PREC2 {
            return true;
        }
        let mut total_ang = 0.0;
        for poly_it in 1..poly_len {
            let poly_vertex = self.mesh_data.get_node(polygon_vertices[poly_it]).location;
            let cur_dir = GpVec2d::from_xy(poly_vertex.coord.subtracted(&center.coord));
            if cur_dir.square_magnitude() < PREC2 {
                return true;
            }
            total_ang += cur_dir.angle(&prev_dir);
            prev_dir = cur_dir;
        }
        (total_ang.abs() - ANGLE_2PI).abs() <= ANGULAR
    }

    #[allow(clippy::too_many_arguments)]
    fn kill_triangles_on_intersecting_links(
        &mut self,
        link_to_check_id: i32,
        end_point: i32,
        the_polygon: &[i32],
        the_poly_boxes: &[BndB2],
        survived_links: &mut HashSet<i32>,
        loop_edges: &mut HashMap<i32, bool>,
    ) {
        if survived_links.contains(&link_to_check_id) {
            return;
        }
        let mut b = BndB2::void();
        let is_not_intersect =
            self.check_intersection(link_to_check_id, the_polygon, the_poly_boxes, false, false, false, &mut b);
        survived_links.insert(link_to_check_id);
        if is_not_intersect {
            return;
        }
        self.kill_link_triangles(link_to_check_id, loop_edges);
        let neighbors: Vec<i32> = self.mesh_data.links_connected_to(end_point).to_vec();
        for &neighbor_id in &neighbors {
            let neighbor_link = self.mesh_data.get_link(neighbor_id);
            let mut other_node = neighbor_link.first_node();
            if other_node == end_point {
                other_node = neighbor_link.last_node();
            }
            self.kill_triangles_on_intersecting_links(neighbor_id, other_node, the_polygon, the_poly_boxes, survived_links, loop_edges);
        }
    }

    fn kill_link_triangles(&mut self, link_id: i32, loop_edges: &mut HashMap<i32, bool>) {
        let elem_nb = self.mesh_data.elements_connected_to(link_id).extent();
        for _ in 0..elem_nb {
            let elem_id = self.mesh_data.elements_connected_to(link_id).first_index();
            if elem_id < 0 {
                continue;
            }
            self.delete_triangle(elem_id, loop_edges);
        }
    }

    fn get_oriented_nodes(&self, edge_id: i32, is_forward: bool, nodes: &mut [i32]) {
        let edge = self.mesh_data.get_link(edge_id);
        if is_forward {
            nodes[0] = edge.first_node();
            nodes[1] = edge.last_node();
        } else {
            nodes[0] = edge.last_node();
            nodes[1] = edge.first_node();
        }
    }

    fn process_loop(&mut self, link_from: usize, link_to: usize, polygon: &[i32], poly_boxes: &[BndB2]) {
        let nb = link_to - link_from - 1;
        if nb < 3 {
            return;
        }
        let mut sub_polygon: Vec<i32> = Vec::new();
        let mut sub_boxes: Vec<BndB2> = Vec::new();
        for i in 0..nb {
            sub_polygon.push(polygon[link_from + i]);
            sub_boxes.push(poly_boxes[link_from + i]);
        }
        self.mesh_polygon(&mut sub_polygon, &mut sub_boxes, &mut None);
    }

    fn create_and_replace_polygon_link(
        &mut self,
        nodes: &[i32],
        pnts: &[GpPnt2d],
        root_index: usize,
        flag: ReplaceFlag,
        polygon: &mut Vec<i32>,
        poly_boxes: &mut Vec<BndB2>,
    ) -> i32 {
        let new_edge_id = self.mesh_data.add_link(nodes[0], nodes[1], VertexState::Free);
        let mut new_box = BndB2::void();
        update_bnd_box(pnts[0].coord, pnts[1].coord, &mut new_box);
        match flag {
            ReplaceFlag::Replace => {
                polygon[root_index - 1] = new_edge_id;
                poly_boxes[root_index - 1] = new_box;
            }
            ReplaceFlag::InsertAfter => {
                polygon.insert(root_index, new_edge_id);
                poly_boxes.insert(root_index, new_box);
            }
            ReplaceFlag::InsertBefore => {
                polygon.insert(root_index - 1, new_edge_id);
                poly_boxes.insert(root_index - 1, new_box);
            }
        }
        new_edge_id
    }

    fn mesh_polygon(&mut self, the_polygon: &mut Vec<i32>, the_poly_boxes: &mut Vec<BndB2>, skipped: &mut Option<HashSet<i32>>) {
        if self.mesh_elementary_polygon(the_polygon) {
            return;
        }

        let mut poly_len = the_polygon.len();
        let poly_area = self.poly_area(the_polygon, 1, poly_len).abs();
        let small_loop_area = 0.001 * poly_area;

        let mut a_poly_it: usize = 1;
        while a_poly_it < poly_len {
            let cur_edge_info = the_polygon[a_poly_it - 1];
            let mut cur_edge_id = cur_edge_info.abs();
            if self.mesh_data.link_movability(cur_edge_id) != VertexState::Frontier {
                a_poly_it += 1;
                continue;
            }
            let mut cur_nodes = [0i32; 2];
            self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut cur_nodes);
            let mut cur_pnts = [
                self.mesh_data.get_node(cur_nodes[0]).location,
                self.mesh_data.get_node(cur_nodes[1]).location,
            ];

            let mut a_next_poly_it = a_poly_it + 1;
            while a_next_poly_it <= poly_len {
                let next_edge_info = the_polygon[a_next_poly_it - 1];
                let next_edge_id = next_edge_info.abs();
                if self.mesh_data.link_movability(next_edge_id) != VertexState::Frontier {
                    a_next_poly_it += 1;
                    continue;
                }
                let mut next_nodes = [0i32; 2];
                self.get_oriented_nodes(next_edge_id, next_edge_info > 0, &mut next_nodes);
                let next_pnts = [
                    self.mesh_data.get_node(next_nodes[0]).location,
                    self.mesh_data.get_node(next_nodes[1]).location,
                ];

                let mut int_pnt = GpPnt2d::zero();
                let int_flag = self.int_seg_seg(cur_edge_id, next_edge_id, false, true, &mut int_pnt);
                if int_flag == IntFlag::NoIntersection {
                    a_next_poly_it += 1;
                    continue;
                }

                let mut is_remove_from_first = false;
                let mut is_add_replacing_edge = true;
                let mut index_to_remove_to = a_next_poly_it;

                match int_flag {
                    IntFlag::Cross => {
                        let mut loop_area = self.poly_area(the_polygon, a_poly_it + 1, a_next_poly_it);
                        let vec1 = GpVec2d::from_xy(cur_pnts[1].coord.subtracted(&int_pnt.coord));
                        let vec2 = GpVec2d::from_xy(next_pnts[0].coord.subtracted(&int_pnt.coord));
                        loop_area += vec1.crossed(&vec2) / 2.0;
                        if loop_area.abs() > small_loop_area {
                            let next_nodes2 = [next_nodes[0], cur_nodes[0]];
                            let next_pnts2 = [next_pnts[0], cur_pnts[0]];
                            self.create_and_replace_polygon_link(
                                &next_nodes2,
                                &next_pnts2,
                                a_next_poly_it,
                                ReplaceFlag::Replace,
                                the_polygon,
                                the_poly_boxes,
                            );
                            self.process_loop(a_poly_it, a_next_poly_it, the_polygon, the_poly_boxes);
                            return;
                        }
                        let dist1 = int_pnt.square_distance(&next_pnts[0]);
                        let dist2 = int_pnt.square_distance(&next_pnts[1]);
                        let is_close_to_start = dist1 < dist2;
                        let end_point_index = if is_close_to_start { 0 } else { 1 };
                        cur_nodes[1] = next_nodes[end_point_index];
                        cur_pnts[1] = next_pnts[end_point_index];
                        if is_close_to_start {
                            index_to_remove_to -= 1;
                        }
                        if let Some(s) = skipped.as_mut() {
                            for sk in a_poly_it..=index_to_remove_to {
                                s.insert(the_polygon[sk - 1].abs());
                            }
                        }
                    }
                    IntFlag::PointOnSegment => {
                        let mut is_first_chopping = false;
                        let mut check_point_it = 0usize;
                        for cpi in 0..2 {
                            let ref_point = cur_pnts[cpi];
                            let v1 = GpVec2d::from_xy(next_pnts[0].coord.subtracted(&ref_point.coord));
                            let v2 = GpVec2d::from_xy(next_pnts[1].coord.subtracted(&ref_point.coord));
                            if v1.crossed(&v2).abs() < PREC {
                                is_first_chopping = true;
                                check_point_it = cpi;
                                break;
                            }
                        }
                        if is_first_chopping {
                            is_add_replacing_edge = false;
                            is_remove_from_first = check_point_it == 0;
                            let split_link = [next_nodes[0], cur_nodes[check_point_it], next_nodes[1]];
                            let split_pnts = [next_pnts[0], cur_pnts[check_point_it], next_pnts[1]];
                            for sli in 0..2 {
                                let flag = if sli == 0 { ReplaceFlag::Replace } else { ReplaceFlag::InsertAfter };
                                self.create_and_replace_polygon_link(
                                    &split_link[sli..sli + 2],
                                    &split_pnts[sli..sli + 2],
                                    a_next_poly_it,
                                    flag,
                                    the_polygon,
                                    the_poly_boxes,
                                );
                            }
                            self.process_loop(a_poly_it + check_point_it, index_to_remove_to, the_polygon, the_poly_boxes);
                        } else {
                            let split_link_nodes = [next_nodes[1], cur_nodes[1]];
                            let split_link_pnts = [next_pnts[1], cur_pnts[1]];
                            self.create_and_replace_polygon_link(
                                &split_link_nodes,
                                &split_link_pnts,
                                a_poly_it,
                                ReplaceFlag::InsertAfter,
                                the_polygon,
                                the_poly_boxes,
                            );
                            cur_nodes[1] = next_nodes[1];
                            cur_pnts[1] = next_pnts[1];
                            index_to_remove_to += 1;
                            self.process_loop(a_poly_it + 1, index_to_remove_to, the_polygon, the_poly_boxes);
                        }
                    }
                    IntFlag::Glued => {
                        if cur_nodes[1] == next_nodes[0] {
                            cur_nodes[1] = next_nodes[1];
                            cur_pnts[1] = next_pnts[1];
                        }
                    }
                    IntFlag::Same => {
                        self.process_loop(a_poly_it, a_next_poly_it, the_polygon, the_poly_boxes);
                        is_remove_from_first = true;
                        is_add_replacing_edge = false;
                    }
                    _ => {
                        a_next_poly_it += 1;
                        continue;
                    }
                }

                if is_add_replacing_edge {
                    cur_edge_id = self.create_and_replace_polygon_link(
                        &cur_nodes,
                        &cur_pnts,
                        a_poly_it,
                        ReplaceFlag::Replace,
                        the_polygon,
                        the_poly_boxes,
                    );
                }

                let index_to_remove_from = if is_remove_from_first { a_poly_it } else { a_poly_it + 1 };
                the_polygon.drain(index_to_remove_from - 1..=index_to_remove_to - 1);
                the_poly_boxes.drain(index_to_remove_from - 1..=index_to_remove_to - 1);
                poly_len = the_polygon.len();
                if is_remove_from_first {
                    a_poly_it -= 1;
                    break;
                }
                a_next_poly_it = a_poly_it;
            }
            a_poly_it += 1;
        }

        // Decomposition of the (possibly corrected) polygon into triangles.
        let mut cut_poly = std::mem::take(the_polygon);
        let mut cut_boxes = std::mem::take(the_poly_boxes);
        let mut pending: Vec<(Vec<i32>, Vec<BndB2>)> = Vec::new();
        loop {
            let mut poly2: Vec<i32> = Vec::new();
            let mut boxes2: Vec<BndB2> = Vec::new();
            self.decompose_simple_polygon(&mut cut_poly, &mut cut_boxes, &mut poly2, &mut boxes2);
            if !poly2.is_empty() {
                pending.push((poly2, boxes2));
            }
            if cut_poly.is_empty() {
                if pending.is_empty() {
                    break;
                }
                let (p, b) = pending.remove(0);
                cut_poly = p;
                cut_boxes = b;
            }
        }
    }

    fn mesh_elementary_polygon(&mut self, the_polygon: &[i32]) -> bool {
        let poly_len = the_polygon.len();
        if poly_len < 3 {
            return true;
        } else if poly_len > 3 {
            return false;
        }
        let mut edges = [0i32; 3];
        let mut oris = [false; 3];
        for i in 0..3 {
            let info = the_polygon[i];
            edges[i] = info.abs();
            oris[i] = info > 0;
        }
        let edge1 = self.mesh_data.get_link(edges[0]);
        let edge2 = self.mesh_data.get_link(edges[1]);
        let mut nodes = [edge1.first_node(), edge1.last_node(), edge2.first_node()];
        if nodes[2] == nodes[0] || nodes[2] == nodes[1] {
            nodes[2] = edge2.last_node();
        }
        self.add_triangle(edges, oris, nodes);
        true
    }

    fn decompose_simple_polygon(
        &mut self,
        the_polygon: &mut Vec<i32>,
        the_poly_boxes: &mut Vec<BndB2>,
        the_polygon_cut: &mut Vec<i32>,
        the_poly_boxes_cut: &mut Vec<BndB2>,
    ) {
        if self.mesh_elementary_polygon(the_polygon) {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }

        let poly_len = the_polygon.len();
        let first_edge_info = the_polygon[0];
        let first_edge = self.mesh_data.get_link(first_edge_info.abs());
        let mut nodes = [0i32; 3];
        self.get_oriented_nodes(first_edge_info.abs(), first_edge_info > 0, &mut nodes[..2]);

        let mut ref_vertices = [
            self.mesh_data.get_node(nodes[0]).location,
            self.mesh_data.get_node(nodes[1]).location,
            GpPnt2d::zero(),
        ];
        let mut ref_edge_dir = GpVec2d::from_xy(ref_vertices[1].coord.subtracted(&ref_vertices[0].coord));
        let ref_edge_len = ref_edge_dir.magnitude();
        if ref_edge_len < PREC {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }
        ref_edge_dir = ref_edge_dir.divided(ref_edge_len);

        let mut used_link_id = 0usize;
        let mut opt_angle = 0.0;
        let mut min_dist = f64::MAX;
        let mut pivot_node = nodes[1];

        for link_it in 3..=poly_len {
            let link_info = the_polygon[link_it - 1];
            let next_edge = self.mesh_data.get_link(link_info.abs());
            pivot_node = if link_info > 0 { next_edge.first_node() } else { next_edge.last_node() };
            if pivot_node == nodes[1] {
                continue;
            }
            let pivot_vertex = self.mesh_data.get_node(pivot_node).location;
            let distance_dir = GpVec2d::from_xy(pivot_vertex.coord.subtracted(&ref_vertices[1].coord));
            let dist = ref_edge_dir.crossed(&distance_dir);
            let angle = ref_edge_dir.angle(&distance_dir).abs();
            let abs_dist = dist.abs();
            if abs_dist < PREC || dist < 0.0 {
                continue;
            }
            if (abs_dist >= min_dist) && (angle <= opt_angle || angle > ANG_DEV_90DEG) {
                continue;
            }

            let mut is_intersect = false;
            for ref_link_node_it in 0..2 {
                let link_first_node = nodes[ref_link_node_it];
                let link_first_vertex = ref_vertices[ref_link_node_it];
                let mut b = BndB2::void();
                update_bnd_box(link_first_vertex.coord, pivot_vertex.coord, &mut b);
                for check_link_it in 2..=poly_len {
                    if check_link_it == link_it {
                        continue;
                    }
                    if b.is_out(&the_poly_boxes[check_link_it - 1]) {
                        continue;
                    }
                    let poly_link = self.mesh_data.get_link(the_polygon[check_link_it - 1].abs());
                    let check_same = (link_first_node == poly_link.first_node() && pivot_node == poly_link.last_node())
                        || (link_first_node == poly_link.last_node() && pivot_node == poly_link.first_node());
                    if check_same {
                        continue;
                    }
                    let mut int_pnt = GpPnt2d::zero();
                    let flag = self.int_seg_seg_nodes(
                        link_first_node,
                        pivot_node,
                        poly_link.first_node(),
                        poly_link.last_node(),
                        false,
                        false,
                        &mut int_pnt,
                    );
                    if flag != IntFlag::NoIntersection {
                        is_intersect = true;
                        break;
                    }
                }
                if is_intersect {
                    break;
                }
            }
            if is_intersect {
                continue;
            }

            opt_angle = angle;
            min_dist = abs_dist;
            nodes[2] = pivot_node;
            ref_vertices[2] = pivot_vertex;
            used_link_id = link_it;
        }

        if used_link_id == 0 {
            the_polygon.clear();
            the_poly_boxes.clear();
            return;
        }

        let new_edge_0 = self.mesh_data.add_link(nodes[1], nodes[2], VertexState::Free);
        let new_edge_1 = self.mesh_data.add_link(nodes[2], nodes[0], VertexState::Free);
        let new_edges_info = [first_edge_info, new_edge_0, new_edge_1];
        self.add_triangle_by_info(new_edges_info, nodes);

        if used_link_id == 3 {
            the_polygon.remove(0);
            the_poly_boxes.remove(0);
            the_polygon[0] = -new_edges_info[2];
            let mut b = BndB2::void();
            update_bnd_box(ref_vertices[0].coord, ref_vertices[2].coord, &mut b);
            the_poly_boxes[0] = b;
        } else {
            if used_link_id < poly_len {
                the_polygon_cut.clear();
                the_polygon_cut.extend_from_slice(&the_polygon[used_link_id..]);
                the_polygon.truncate(used_link_id);
                the_polygon_cut.insert(0, -new_edges_info[2]);
                the_poly_boxes_cut.clear();
                the_poly_boxes_cut.extend_from_slice(&the_poly_boxes[used_link_id..]);
                the_poly_boxes.truncate(used_link_id);
                let mut b = BndB2::void();
                update_bnd_box(ref_vertices[0].coord, ref_vertices[2].coord, &mut b);
                the_poly_boxes_cut.insert(0, b);
            } else {
                the_polygon.pop();
                the_poly_boxes.pop();
            }
            the_polygon[0] = -new_edges_info[1];
            let mut b = BndB2::void();
            update_bnd_box(ref_vertices[1].coord, ref_vertices[2].coord, &mut b);
            the_poly_boxes[0] = b;
        }
    }

    /// Re-triangulates the cavity left after deleting the triangles of a vertex
    /// (used by `RemoveVertex`). Collects the loop of free edges around the
    /// cavity and meshes it as a polygon.
    fn mesh_polygon_of_cavity(&mut self, loop_edges: &mut HashMap<i32, bool>) {
        let mut boxes: Vec<BndB2> = Vec::new();
        let mut polygon: Vec<i32> = Vec::new();
        let mut loop_edges_count = loop_edges.len();
        let mut keys: Vec<i32> = loop_edges.keys().copied().collect();
        if keys.is_empty() {
            return;
        }
        let an_edge_id = keys[0];
        let edge = self.mesh_data.get_link(an_edge_id);
        let mut first_node = edge.first_node();
        let mut pivot_node = edge.last_node();
        let mut edge_id = an_edge_id;

        let is_positive = loop_edges[&edge_id];
        if !is_positive {
            let tmp = first_node;
            first_node = pivot_node;
            pivot_node = tmp;
            polygon.push(-edge_id);
        } else {
            polygon.push(edge_id);
        }
        let mut b = BndB2::void();
        update_bnd_box(
            self.mesh_data.get_node(first_node).location.coord,
            self.mesh_data.get_node(pivot_node).location.coord,
            &mut b,
        );
        boxes.push(b);
        loop_edges.remove(&edge_id);

        let last_node = first_node;
        while pivot_node != last_node {
            let links: Vec<i32> = self.mesh_data.links_connected_to(pivot_node).to_vec();
            let mut advanced = false;
            for &link_value in &links {
                if link_value != edge_id && loop_edges.contains_key(&link_value) {
                    edge_id = link_value;
                    let e = self.mesh_data.get_link(edge_id);
                    // Traverse the next cavity edge away from the pivot node: if
                    // the link is stored pivot->last use it forward, otherwise
                    // reverse it so the polygon advances to the other endpoint.
                    let current_node = if e.first_node() == pivot_node {
                        polygon.push(edge_id);
                        e.last_node()
                    } else {
                        polygon.push(-edge_id);
                        e.first_node()
                    };
                    let mut b2 = BndB2::void();
                    update_bnd_box(
                        self.mesh_data.get_node(current_node).location.coord,
                        self.mesh_data.get_node(pivot_node).location.coord,
                        &mut b2,
                    );
                    boxes.push(b2);
                    pivot_node = current_node;
                    loop_edges.remove(&edge_id);
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
            if loop_edges_count <= 0 {
                break;
            }
            loop_edges_count -= 1;
        }
        if polygon.len() >= 3 {
            self.mesh_polygon(&mut polygon, &mut boxes, &mut None);
        }
    }

    fn get_edges_by_type(&self, edge_type: VertexState) -> HashSet<i32> {
        let mut result = HashSet::new();
        for &edge in self.mesh_data.links_of_domain() {
            let is_to_add = if edge_type == VertexState::Free {
                self.mesh_data.elements_connected_to(edge).extent() <= 1
            } else {
                self.mesh_data.link_movability(edge) == edge_type
            };
            if is_to_add {
                result.insert(edge);
            }
        }
        result
    }

    fn calculate_dist(
        &self,
        v_edges: &[GpXY; 3],
        points: &[GpPnt2d; 3],
        vertex: &DelaunVertex,
        distance: &mut [f64; 3],
        sq_modulus: &mut [f64; 3],
        edge_on: &mut usize,
    ) -> f64 {
        let mut min_dist = f64::MAX;
        let v = vertex.location.coord;
        for i in 0..3 {
            sq_modulus[i] = v_edges[i].square_modulus();
            if sq_modulus[i] <= PREC2 {
                return -1.0;
            }
            distance[i] = v_edges[i].crossed(&v.subtracted(&points[i].coord));
            let d = distance[i] * distance[i] / sq_modulus[i];
            if d < min_dist {
                *edge_on = i;
                min_dist = d;
            }
        }
        min_dist
    }

    fn int_seg_seg(
        &self,
        edge1_id: i32,
        edge2_id: i32,
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        int_pnt: &mut GpPnt2d,
    ) -> IntFlag {
        let e1 = self.mesh_data.get_link(edge1_id);
        let e2 = self.mesh_data.get_link(edge2_id);
        self.int_seg_seg_nodes(
            e1.first_node(),
            e1.last_node(),
            e2.first_node(),
            e2.last_node(),
            is_consider_end_point_touch,
            is_consider_point_on_edge,
            int_pnt,
        )
    }

    fn int_seg_seg_nodes(
        &self,
        n1: i32,
        n2: i32,
        n3: i32,
        n4: i32,
        is_consider_end_point_touch: bool,
        is_consider_point_on_edge: bool,
        int_pnt: &mut GpPnt2d,
    ) -> IntFlag {
        let p1 = self.mesh_data.get_node(n1).location.coord;
        let p2 = self.mesh_data.get_node(n2).location.coord;
        let p3 = self.mesh_data.get_node(n3).location.coord;
        let p4 = self.mesh_data.get_node(n4).location.coord;
        let (flag, pnt) = GeomTool::int_seg_seg(&p1, &p2, &p3, &p4, is_consider_end_point_touch, is_consider_point_on_edge);
        *int_pnt = GpPnt2d::from_xy(pnt);
        flag
    }

    fn poly_area(&self, the_polygon: &[i32], start_1based: usize, end_1based: usize) -> f64 {
        let mut area = 0.0;
        let poly_len = the_polygon.len();
        if start_1based >= end_1based || start_1based > poly_len {
            return area;
        }
        let mut cur_edge_info = the_polygon[start_1based - 1];
        let mut cur_edge_id = cur_edge_info.abs();
        let mut nodes = [0i32; 2];
        self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut nodes);
        let ref_pnt = self.mesh_data.get_node(nodes[0]).location;
        for poly_it in (start_1based + 1)..=end_1based {
            cur_edge_info = the_polygon[poly_it - 1];
            cur_edge_id = cur_edge_info.abs();
            self.get_oriented_nodes(cur_edge_id, cur_edge_info > 0, &mut nodes);
            let v1 = GpVec2d::from_xy(self.mesh_data.get_node(nodes[0]).location.coord.subtracted(&ref_pnt.coord));
            let v2 = GpVec2d::from_xy(self.mesh_data.get_node(nodes[1]).location.coord.subtracted(&ref_pnt.coord));
            area += v1.crossed(&v2);
        }
        area / 2.0
    }

    fn is_sup_vertex(&self, vertex_idx: i32) -> bool {
        self.sup_vert.contains(&vertex_idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt2d;

    fn v(u: f64, w: f64) -> DelaunVertex {
        DelaunVertex::new_parametric(u, w, VertexState::Free)
    }

    fn tri_verts(ds: &DelaunDataStructure, t: &DelaunTriangle) -> [GpPnt2d; 3] {
        let n = t.vertex_indices;
        [
            ds.get_node(n[0]).location,
            ds.get_node(n[1]).location,
            ds.get_node(n[2]).location,
        ]
    }

    fn signed_area(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
        (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
    }

    /// In-circle test: is `p` strictly inside the circumcircle of (a,b,c)?
    fn in_circumcircle(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
        let orient = signed_area(a, b, c);
        if orient == 0.0 {
            return false;
        }
        let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
        let (bx, by) = (b.x() - p.x(), b.y() - p.y());
        let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
        let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
            - (bx * bx + by * by) * (ax * cy - ay * cx)
            + (cx * cx + cy * cy) * (ax * by - ay * bx);
        det * orient > 0.0
    }

    /// Computes the convex hull of the triangulation's boundary: boundary edges
    /// (exactly one connected triangle) form the hull polygon.
    fn hull_vertices(ds: &DelaunDataStructure) -> Vec<i32> {
        let mut boundary: Vec<(i32, i32)> = Vec::new();
        for &e in ds.links_of_domain() {
            if ds.elements_connected_to(e).extent() == 1 {
                let l = ds.get_link(e);
                boundary.push((l.first_node(), l.last_node()));
            }
        }
        // Walk the boundary polygon, never stepping back across the edge just
        // traversed (the stored link orientation is arbitrary).
        let mut hull = Vec::new();
        if boundary.is_empty() {
            return hull;
        }
        let (start, _) = boundary[0];
        let mut cur = start;
        let mut prev = -1;
        hull.push(cur);
        for _ in 0..boundary.len() {
            let mut next = -1;
            for &(a, b) in &boundary {
                if a == cur && b != prev {
                    next = b;
                    break;
                }
                if b == cur && a != prev {
                    next = a;
                    break;
                }
            }
            if next < 0 || next == start {
                break;
            }
            prev = cur;
            cur = next;
            hull.push(cur);
        }
        hull
    }

    #[test]
    fn grid_3x3_triangle_count() {
        // 9 points, 4 hull corners => 2*9-2-4 = 12 triangles.
        let pts: Vec<DelaunVertex> = (0..3)
            .flat_map(|i| (0..3).map(move |j| v(i as f64, j as f64)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let n = ds.elements_of_domain().len();
        // 9 lattice points: 8 boundary points (4 corners + 4 edge midpoints) on
        // an octagonal hull + 1 interior point => 2*9-2-8 = 8 triangles.
        assert_eq!(n, 8, "3x3 grid must triangulate to 8 triangles, got {n}");
        let h = hull_vertices(ds).len();
        assert_eq!(h, 8, "3x3 grid hull must have 8 boundary vertices, got {h}");
        assert_eq!(n, 2 * 9 - 2 - h);
    }

    #[test]
    fn random_points_triangle_count_and_hull() {
        // Deterministic pseudo-random set covering a spread of sizes.
        let mut rng_state = 0x9E3779B97F4A7C15u64;
        let mut rnd = move || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state as f64 / u64::MAX as f64) * 100.0 - 50.0
        };
        for &n in &[8, 12, 25] {
            let pts: Vec<DelaunVertex> = (0..n).map(|_| v(rnd(), rnd())).collect();
            let delaun = Delaun::new_vertices(&pts);
            let ds = delaun.result();
            let tris = ds.elements_of_domain().len();
            let h = hull_vertices(ds).len();
            assert!(
                tris == 2 * n - 2 - h,
                "N={n}: expected 2N-2-h = {}, got {tris} (h={h})",
                2 * n - 2 - h
            );
        }
    }

    #[test]
    fn all_triangles_non_degenerate() {
        let pts: Vec<DelaunVertex> = (0..4)
            .flat_map(|i| (0..4).map(move |j| v(i as f64 * 0.7, j as f64 * 0.9)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let ids: Vec<i32> = ds.elements_of_domain().iter().copied().collect();
        assert!(ids.len() >= 10);
        for id in ids {
            let t = ds.get_element(id);
            let p = tri_verts(ds, &t);
            let area = signed_area(&p[0], &p[1], &p[2]);
            assert!(area.abs() > 1e-9, "degenerate triangle {id}: area {area}");
            assert!(
                p[0] != p[1] && p[1] != p[2] && p[0] != p[2],
                "repeated vertex in triangle {id}"
            );
        }
    }

    #[test]
    fn empty_circle_property() {
        let mut rng_state = 0x2545F4914F6CDD1Du64;
        let mut rnd = move || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state as f64 / u64::MAX as f64) * 20.0 - 10.0
        };
        let pts: Vec<DelaunVertex> = (0..40).map(|_| v(rnd(), rnd())).collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let tris: Vec<(i32, [i32; 3])> = ds
            .elements_of_domain()
            .iter()
            .map(|&id| (id, ds.get_element(id).vertex_indices))
            .collect();

        // Build undirected edge -> triangle id map.
        let mut edge_tris: HashMap<(i32, i32), Vec<i32>> = HashMap::new();
        for &(id, verts) in &tris {
            for k in 0..3 {
                let (a, b) = (verts[k], verts[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                edge_tris.entry(key).or_default().push(id);
            }
        }

        // For every pair of triangles that do NOT share an edge, assert no
        // vertex of A lies strictly inside B's circumcircle.
        for i in 0..tris.len() {
            for j in (i + 1)..tris.len() {
                let (id_a, verts_a) = tris[i];
                let (id_b, verts_b) = tris[j];
                let mut share_edge = false;
                for ka in 0..3 {
                    let (a, b) = (verts_a[ka], verts_a[(ka + 1) % 3]);
                    let key = if a < b { (a, b) } else { (b, a) };
                    if let Some(list) = edge_tris.get(&key) {
                        if list.len() > 1 && list.contains(&id_a) && list.contains(&id_b) {
                            share_edge = true;
                        }
                    }
                }
                if share_edge {
                    continue;
                }
                let pa = tri_verts(ds, &ds.get_element(id_a));
                let pb = tri_verts(ds, &ds.get_element(id_b));
                for p in &pa {
                    assert!(
                        !in_circumcircle(&pb[0], &pb[1], &pb[2], p),
                        "vertex of triangle {id_a} inside circumcircle of {id_b}"
                    );
                }
                for p in &pb {
                    assert!(
                        !in_circumcircle(&pa[0], &pa[1], &pa[2], p),
                        "vertex of triangle {id_b} inside circumcircle of {id_a}"
                    );
                }
            }
        }
    }

    #[test]
    fn super_triangle_vertices_absent() {
        let pts: Vec<DelaunVertex> = (0..4)
            .flat_map(|i| (0..4).map(move |j| v(i as f64, j as f64)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        // The three super vertices are the last three node slots (1-based
        // indices n-2, n-1, n) and were removed (marked Deleted) after the
        // auxiliary elements were destroyed.
        let n = ds.nb_nodes();
        assert!(n >= pts.len() + 3);
        for i in (n - 2)..=n {
            assert_eq!(
                ds.get_node(i as i32).state,
                VertexState::Deleted,
                "super vertex {i} must be removed"
            );
        }
        // No triangle may reference a super vertex (super vertices are >= n-2).
        for &id in ds.elements_of_domain() {
            let t = ds.get_element(id);
            for w in t.vertex_indices {
                assert!(
                    w < (n - 2) as i32,
                    "triangle {id} references super vertex {w}"
                );
            }
        }
    }

    #[test]
    fn hull_edges_classified_as_free() {
        // A single triangle has three hull edges, each with exactly one
        // connected element, hence classified as Free edges.
        let pts = vec![v(0.0, 0.0), v(3.0, 0.0), v(0.0, 2.0)];
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        assert_eq!(ds.elements_of_domain().len(), 1);
        let free = delaun.free_edges();
        assert_eq!(free.len(), 3, "single triangle hull must expose 3 free edges");
        assert!(delaun.frontier().is_empty());
        assert!(delaun.internal_edges().is_empty());
    }

    #[test]
    fn frontier_edges_are_classified() {
        // Seed the structure with a frontier link and verify classification.
        let mut data = DelaunDataStructure::new(16);
        let a = data.add_node(v(0.0, 0.0));
        let b = data.add_node(v(1.0, 0.0));
        let c = data.add_node(v(0.0, 1.0));
        let d = data.add_node(v(1.0, 1.0));
        let _ab = data.add_link(a, b, VertexState::Frontier);
        let cd = data.add_link(c, d, VertexState::Fixed);
        let _bd = data.add_link(b, d, VertexState::Free);
        let mut delaun = Delaun {
            mesh_data: data,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            sup_trian: DelaunTriangle::default(),
        };
        // Free edges are links with at most one connected element. None of the
        // three links has any element, so all three are classified Free.
        let free = delaun.free_edges();
        assert!(free.contains(&_bd));
        assert_eq!(free.len(), 3);
        let frontier = delaun.frontier();
        assert!(frontier.contains(&_ab));
        assert_eq!(frontier.len(), 1);
        let internal = delaun.internal_edges();
        assert!(internal.contains(&cd));
        assert_eq!(internal.len(), 1);
    }

    #[test]
    fn remove_vertex_retriangulates() {
        // 3x3 grid; removing the center vertex must reduce the triangle count
        // by 2 (the center's link is a triangle fan of 6 -> 4).
        let pts: Vec<DelaunVertex> = (0..3)
            .flat_map(|i| (0..3).map(move |j| v(i as f64, j as f64)))
            .collect();
        let mut delaun = Delaun::new_vertices(&pts);
        let ds0 = delaun.result().clone();
        let n0 = ds0.elements_of_domain().len();
        assert_eq!(n0, 8);
        let center = *ds0.get_node(5); // node (1,1) is 5th 0-based => id 5 (1-based)
        delaun.remove_vertex(&center);
        let ds1 = delaun.result();
        let n1 = ds1.elements_of_domain().len();
        // Removing the interior vertex from the octagon-hull mesh (8 triangles)
        // re-triangulates the octagonal cavity into 6 triangles.
        assert_eq!(n1, 6, "center removal must leave 6 triangles, got {n1}");
    }
}
