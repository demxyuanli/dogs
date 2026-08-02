//! Port of OCCT `BRepMesh_DelaunayDeflectionControlMeshAlgo` — deflection
//! controlled refinement of a 2D Delaunay triangulation (Wave 3).
//!
//! Source:
//! `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx`
//!
//! The OCCT class is a template extending `BRepMesh_DelaunayNodeInsertionMeshAlgo`
//! (the sibling `node_insertion` stub). It walks the generated UV-space mesh and
//! splits every triangle whose linear deflection (chord sagitta) or angular
//! deflection (surface normal change) exceeds the requested tolerances:
//!
//! * `compute_triangle_geometry` evaluates the triangle's 3D normal and 2D area,
//!   rejecting degenerate triangles.
//! * `split_links` checks each link's midpoint: the `LineDeviation` functor
//!   measures how far the surface point sits off the 3D chord, and
//!   `check_link_ends_for_angular_deviation` measures the surface-normal angle
//!   between the two endpoints. Control points are collected and handed to the
//!   caller, which re-triangulates, until convergence or `MinSize` rejection.
//!
//! ponytail: the OCCT class is templated on a `RangeSplitter` that maps real UV
//! to the integer grid of the Delaunay structure; this port operates directly on
//! the raw parametric UV coordinates (the `Delaun` port already does so).
//! `reject_by_min_size` scans all live triangles instead of OCCT's
//! cell-filtered circumcircle tool (`CircleTool`), which is private to `delaun`.

use std::collections::HashSet;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpXY};
use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};
use occt_geom::Surface;

use super::delaun::Delaun;
use super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::geom_tool::GeomTool;
use super::parameters::MeshParameters;

/// Read access to a triangulated mesh required by the deflection control.
///
/// This is the Rust counterpart of OCCT's `getStructure()` accessor on the
/// data structure of Delaun. `Delaun` implements it (all needed queries are
/// already public); a test double may implement it for isolated drives of the
/// refinement loop.
pub trait DeflectionMesh {
    /// Ids of live triangles. Source: `ElementsOfDomain`.
    fn elements_of_domain(&self) -> Vec<i32>;
    /// Triangle by 1-based id. Source: `GetElement`.
    fn triangle(&self, id: i32) -> DelaunTriangle;
    /// Link by 1-based id. Source: `GetLink`.
    fn link(&self, id: i32) -> DelaunLink;
    /// Movability of a link. Source: `BRepMesh_Edge::Movability`.
    fn link_movability(&self, id: i32) -> VertexState;
    /// Vertex by 1-based id. Source: `GetNode`.
    fn vertex(&self, id: i32) -> DelaunVertex;
}

impl DeflectionMesh for Delaun {
    fn elements_of_domain(&self) -> Vec<i32> {
        self.result().elements_of_domain().iter().copied().collect()
    }
    fn triangle(&self, id: i32) -> DelaunTriangle {
        self.get_triangle(id)
    }
    fn link(&self, id: i32) -> DelaunLink {
        self.get_edge(id)
    }
    fn link_movability(&self, id: i32) -> VertexState {
        self.result().link_movability(id)
    }
    fn vertex(&self, id: i32) -> DelaunVertex {
        self.get_vertex(id)
    }
}

/// Functor computing the squared deflection of a point from a reference plane
/// (a reference point plus a unit normal). Source: `NormalDeviation`.
#[derive(Clone, Copy)]
struct NormalDeviation<'a> {
    ref_pnt: &'a GpPnt,
    normal: GpVec,
}

impl<'a> NormalDeviation<'a> {
    fn new(ref_pnt: &'a GpPnt, normal: GpVec) -> Self {
        Self { ref_pnt, normal }
    }
}

/// Functor computing the squared deflection of a point from a chord.
/// Source: `LineDeviation`.
#[derive(Clone, Copy)]
struct LineDeviation<'a> {
    pnt1: &'a GpPnt,
    pnt2: &'a GpPnt,
}

impl<'a> LineDeviation<'a> {
    fn new(pnt1: &'a GpPnt, pnt2: &'a GpPnt) -> Self {
        Self { pnt1, pnt2 }
    }
}

/// Abstraction of the two deflection functors.
trait SquareDeviation {
    fn square_deviation(&self, point: &GpPnt) -> f64;
}

impl SquareDeviation for NormalDeviation<'_> {
    fn square_deviation(&self, point: &GpPnt) -> f64 {
        let v = GpVec::from_pnts(self.ref_pnt, point);
        let d = self.normal.dot(&v).abs();
        d * d
    }
}

impl SquareDeviation for LineDeviation<'_> {
    fn square_deviation(&self, point: &GpPnt) -> f64 {
        GeomTool::square_deflection_of_segment(self.pnt1, self.pnt2, point)
    }
}

/// Geometrical data of a node of a triangle. Source: `TriangleNodeInfo`.
#[derive(Clone, Copy)]
struct TriangleNodeInfo {
    /// Position in parametric (UV) space.
    point_2d: GpXY,
    /// Associated 3D point on the surface.
    point_3d: GpPnt,
    /// Whether the link opposite this node is a frontier (constraint) link.
    is_frontier_link: bool,
}

impl Default for TriangleNodeInfo {
    fn default() -> Self {
        Self {
            point_2d: GpXY::zero(),
            point_3d: GpPnt::zero(),
            is_frontier_link: false,
        }
    }
}

/// Deflection-controlled refinement of a Delaunay triangulation.
///
/// Port of `BRepMesh_DelaunayDeflectionControlMeshAlgo`. The class holds no
/// reference to a `MeshFace`/`MeshSurface`; the surface evaluation and the
/// target face deflection are passed explicitly to [`Self::post_process`]
/// (`ponytail: MeshFace::surface()` + `MeshFace::deflection()` are the two
/// values the OCCT class reads off `getDFace()`).
pub struct DelaunayDeflectionControlMeshAlgo {
    parameters: MeshParameters,
    /// Square of the minimum allowed triangle edge (`MinSize * MinSize`).
    sq_min_size: f64,
    /// Face deflection used as the linear threshold.
    face_deflection: f64,
    /// Largest squared deflection observed in the current optimization run.
    max_sq_deflection: f64,
    /// Becomes `false` as soon as a non-degenerate triangle is processed.
    is_all_degenerated: bool,
    /// Control points (2D) accumulated during the current pass.
    control_nodes: Vec<GpPnt2d>,
    /// Map of oriented links already processed (dedup across passes).
    processed_couples: HashSet<(i32, i32)>,
}

impl DelaunayDeflectionControlMeshAlgo {
    /// Constructs the refinement algorithm from the meshing parameters.
    ///
    /// The effective `MinSize`/`AngleInterior`/`DeflectionInterior` are resolved
    /// exactly like `BRepMesh_IncrementalMesh::initParameters` (negative sentinels
    /// inherit from the base deflection/angle values).
    pub fn new(parameters: &MeshParameters) -> Self {
        let min_size = effective_min_size(parameters);
        Self {
            parameters: parameters.clone(),
            sq_min_size: min_size * min_size,
            face_deflection: 0.0,
            max_sq_deflection: -1.0,
            is_all_degenerated: false,
            control_nodes: Vec::new(),
            processed_couples: HashSet::new(),
        }
    }

    /// Performs deflection-controlled refinement of the mesh.
    ///
    /// Runs up to 11 passes (OCCT's `aIterationsNb`). Each pass scans every live
    /// triangle via `split_triangle_geometry`, collecting 2D control points whose
    /// linear or angular deflection exceeds the tolerances. The collected points
    /// are handed to `insert`, which must insert them into the mesh and
    /// re-triangulate (the OCCT counterpart is `DelaunayNodeInsertionMeshAlgo::
    /// insertNodes`). Iteration stops when no control point is produced, all
    /// triangles are degenerate, or the pass limit is reached.
    ///
    /// `surface` is the face surface (`MeshSurface`), `face_deflection` is the
    /// current face deflection (`MeshFace::deflection`). Returns the largest
    /// deflection achieved (square root of the max squared deflection), 0 when
    /// the mesh is flat or no geometry could be computed.
    pub fn post_process<M: DeflectionMesh>(
        &mut self,
        mesh: &mut M,
        surface: &dyn Surface,
        face_deflection: f64,
        insert: &mut dyn FnMut(&mut M, &[GpPnt2d]) -> bool,
    ) -> f64 {
        const MAX_PASSES: usize = 11;

        self.face_deflection = face_deflection;
        self.processed_couples.clear();
        self.max_sq_deflection = -1.0;
        self.is_all_degenerated = false;

        let mut is_inserted = true;
        for _pass in 0..MAX_PASSES {
            if !is_inserted || self.is_all_degenerated {
                break;
            }

            // Reset per-pass state.
            self.max_sq_deflection = -1.0;
            self.is_all_degenerated = true;
            self.control_nodes.clear();

            let ids: Vec<i32> = mesh.elements_of_domain();
            if ids.is_empty() {
                break;
            }

            for id in ids {
                let triangle = mesh.triangle(id);
                self.split_triangle_geometry(mesh, surface, &triangle);
            }

            if self.control_nodes.is_empty() {
                // No control node produced => nothing to insert; stop refining.
                break;
            }
            is_inserted = insert(mesh, &self.control_nodes);
        }

        if self.max_sq_deflection >= 0.0 {
            self.max_sq_deflection.sqrt()
        } else {
            0.0
        }
    }

    /// Checks geometry of a triangle; if it fails the deflection check, inserts
    /// the control point (triangle center or link midpoints) into the pending list.
    /// Source: `splitTriangleGeometry`.
    fn split_triangle_geometry<M: DeflectionMesh>(
        &mut self,
        mesh: &M,
        surface: &dyn Surface,
        triangle: &DelaunTriangle,
    ) {
        // OCCT skips Deleted triangles here, but `elements_of_domain` only yields
        // live triangles, so the check is redundant for this port.
        let node_indices = triangle.vertex_indices;

        let mut nodes_info = [TriangleNodeInfo::default(); 3];
        self.get_triangle_info(mesh, triangle, node_indices, &mut nodes_info);

        let mut link_vecs = [GpVec::zero(); 3];
        let mut normal = GpVec::zero();
        if self.compute_triangle_geometry(&nodes_info, &mut link_vecs, &mut normal) {
            self.is_all_degenerated = false;

            let center_2d = nodes_info[0]
                .point_2d
                .added(&nodes_info[1].point_2d)
                .added(&nodes_info[2].point_2d)
                .multiplied_scalar(1.0 / 3.0);
            self.use_point(
                mesh,
                surface,
                GpPnt2d::from_xy(center_2d),
                &NormalDeviation::new(&nodes_info[0].point_3d, normal),
            );
            self.split_links(mesh, surface, &nodes_info, &node_indices);
        }
    }

    /// Returns the 2D/3D data of the three nodes of a triangle, together with the
    /// frontier flag of each opposite link. Source: `getTriangleInfo`.
    fn get_triangle_info<M: DeflectionMesh>(
        &self,
        mesh: &M,
        triangle: &DelaunTriangle,
        node_indices: [i32; 3],
        out: &mut [TriangleNodeInfo; 3],
    ) {
        let link_ids = triangle.link_indices;
        for i in 0..3 {
            let vertex = mesh.vertex(node_indices[i]);
            out[i].point_2d = vertex.location.coord;
            out[i].point_3d = vertex.p3d;
            out[i].is_frontier_link =
                mesh.link_movability(link_ids[i].abs()) == VertexState::Frontier;
        }
    }

    /// Computes the 3D link vectors and the triangle normal.
    /// Returns `false` on a degenerate triangle. Source: `computeTriangleGeometry`.
    fn compute_triangle_geometry(
        &self,
        nodes_info: &[TriangleNodeInfo; 3],
        links: &mut [GpVec; 3],
        normal: &mut GpVec,
    ) -> bool {
        self.check_triangle_for_degenerativity_and_get_links(nodes_info, links)
            && self.check_triangle_area_2d(nodes_info)
            && self.compute_normal(&links[0], &links[1], normal)
    }

    /// Computes the three link vectors; rejects degenerate (zero-length) links.
    /// Source: `checkTriangleForDegenerativityAndGetLinks`.
    fn check_triangle_for_degenerativity_and_get_links(
        &self,
        nodes_info: &[TriangleNodeInfo; 3],
        links: &mut [GpVec; 3],
    ) -> bool {
        const MINIMAL_SQ_LENGTH_3D: f64 = 1.0e-12;
        for i in 0..3 {
            links[i] =
                GpVec::from_pnts(&nodes_info[i].point_3d, &nodes_info[(i + 1) % 3].point_3d);
            if links[i].square_magnitude() < MINIMAL_SQ_LENGTH_3D {
                return false;
            }
        }
        true
    }

    /// Checks the triangle area in parametric space for degeneracy.
    /// Source: `checkTriangleArea2d`.
    fn check_triangle_area_2d(&self, nodes_info: &[TriangleNodeInfo; 3]) -> bool {
        const MINIMAL_AREA_2D: f64 = 1.0e-9;
        let link_2d_1 = nodes_info[1].point_2d.subtracted(&nodes_info[0].point_2d);
        let link_2d_2 = nodes_info[2].point_2d.subtracted(&nodes_info[1].point_2d);
        link_2d_1.crossed(&link_2d_2).abs() > MINIMAL_AREA_2D
    }

    /// Computes the unit normal from two link vectors.
    /// Returns `false` if the normal has null magnitude. Source: `computeNormal`.
    fn compute_normal(&self, link1: &GpVec, link2: &GpVec, normal: &mut GpVec) -> bool {
        let cross = link1.crossed(link2);
        if cross.square_magnitude() > RESOLUTION {
            *normal = cross.normalized();
            return true;
        }
        false
    }

    /// Checks the deflection of the three link midpoints; points failing the
    /// linear, MinSize and angular checks are queued for insertion.
    /// Source: `splitLinks`.
    fn split_links<M: DeflectionMesh>(
        &mut self,
        mesh: &M,
        surface: &dyn Surface,
        nodes_info: &[TriangleNodeInfo; 3],
        node_indices: &[i32; 3],
    ) {
        for i in 0..3 {
            if nodes_info[i].is_frontier_link {
                continue;
            }
            let j = (i + 1) % 3;

            // Deduplicate: only process each undirected link once per optimization.
            let (first, last) = if node_indices[i] < node_indices[j] {
                (node_indices[i], node_indices[j])
            } else {
                (node_indices[j], node_indices[i])
            };
            if !self.processed_couples.insert((first, last)) {
                continue;
            }

            let mid_2d = nodes_info[i]
                .point_2d
                .added(&nodes_info[j].point_2d)
                .multiplied_scalar(0.5);
            let mid_2d = GpPnt2d::from_xy(mid_2d);

            let line = LineDeviation::new(&nodes_info[i].point_3d, &nodes_info[j].point_3d);
            if !self.use_point(mesh, surface, mid_2d, &line) {
                let reject_min_size =
                    self.reject_split_links_for_min_size(surface, &nodes_info[i], &nodes_info[j], mid_2d);
                let reject_angular = self.check_link_ends_for_angular_deviation(
                    surface,
                    &nodes_info[i],
                    &nodes_info[j],
                    mid_2d,
                );
                if !reject_min_size && !reject_angular {
                    self.control_nodes.push(mid_2d);
                }
            }
        }
    }

    /// Checks that the two links produced by splitting the given link by its
    /// midpoint satisfy the MinSize requirement. Source:
    /// `rejectSplitLinksForMinSize`.
    fn reject_split_links_for_min_size(
        &self,
        surface: &dyn Surface,
        node1: &TriangleNodeInfo,
        node2: &TriangleNodeInfo,
        mid: GpPnt2d,
    ) -> bool {
        let mid_3d = self.get_point_3d(surface, mid);
        node1.point_3d.square_distance(&mid_3d) < self.sq_min_size
            || node2.point_3d.square_distance(&mid_3d) < self.sq_min_size
    }

    /// Checks the point (located between the given nodes) against the angular
    /// deviation of the surface normals at the two endpoints. Returns `true` when
    /// the angular deflection is acceptable (no split required). Source:
    /// `checkLinkEndsForAngularDeviation`.
    fn check_link_ends_for_angular_deviation(
        &self,
        surface: &dyn Surface,
        node1: &TriangleNodeInfo,
        node2: &TriangleNodeInfo,
        _mid: GpPnt2d,
    ) -> bool {
        let angle_interior = effective_angle_interior(&self.parameters);
        match (
            GeomTool::normal_on_surface(surface, node1.point_2d.x(), node1.point_2d.y()),
            GeomTool::normal_on_surface(surface, node2.point_2d.x(), node2.point_2d.y()),
        ) {
            (Ok((_, normal1)), Ok((_, normal2))) => normal1.angle(&normal2) <= angle_interior,
            // `GeomLib::NormEstim` returning a non-zero status (unable to estimate
            // the normal) makes OCCT accept the link without an angular split.
            _ => true,
        }
    }

    /// Returns the 3D surface point for the given 2D parameter.
    /// Source: `getPoint3d` (via `IMeshData_Face::GetSurface()->D0`).
    fn get_point_3d(&self, surface: &dyn Surface, p2d: GpPnt2d) -> GpPnt {
        surface.d0(p2d.x(), p2d.y())
    }

    /// Computes the deflection of a point; caches it for insertion when it
    /// exceeds the tolerance. Returns `true` if the point has been cached.
    /// Source: `usePoint`.
    fn use_point<M: DeflectionMesh, D: SquareDeviation>(
        &mut self,
        mesh: &M,
        surface: &dyn Surface,
        p2d: GpPnt2d,
        functor: &D,
    ) -> bool {
        let p3d = self.get_point_3d(surface, p2d);
        if !self.check_deflection_of_point_and_update_cache(mesh, p2d, p3d, functor.square_deviation(&p3d)) {
            self.control_nodes.push(p2d);
            return true;
        }
        false
    }

    /// Checks the given point against the linear deflection, updating the total
    /// mesh deflection. Returns `true` when the point fits. Source:
    /// `checkDeflectionOfPointAndUpdateCache`.
    fn check_deflection_of_point_and_update_cache<M: DeflectionMesh>(
        &mut self,
        mesh: &M,
        _p2d: GpPnt2d,
        p3d: GpPnt,
        sq_deflection: f64,
    ) -> bool {
        if sq_deflection > self.max_sq_deflection {
            self.max_sq_deflection = sq_deflection;
        }

        let sq_face = self.face_deflection * self.face_deflection;
        if sq_deflection < sq_face {
            return true;
        }
        self.reject_by_min_size(mesh, p3d)
    }

    /// Checks the distance between the given node and the nodes of the triangles
    /// shot by it against the MinSize criteria. Source: `rejectByMinSize`.
    fn reject_by_min_size<M: DeflectionMesh>(&self, mesh: &M, p3d: GpPnt) -> bool {
        // ponytail: OCCT queries the cell-filtered circumcircle tool (`CircleTool`,
        // private to `delaun`) for triangles shot by the point; this port scans
        // every live triangle's nodes instead — same MinSize semantics, O(N) per
        // query instead of O(1). Swap when `CircleTool::select` becomes public.
        let mut used_nodes: HashSet<i32> = HashSet::new();
        for &id in mesh.elements_of_domain().iter() {
            let triangle = mesh.triangle(id);
            for node_id in triangle.vertex_indices {
                if used_nodes.insert(node_id) {
                    let vertex = mesh.vertex(node_id);
                    if p3d.square_distance(&vertex.p3d) < self.sq_min_size {
                        return true;
                    }
                }
            }
        }
        false
    }
}

/// Resolves `MinSize` exactly like `BRepMesh_IncrementalMesh::initParameters`:
/// a value below `Confusion` inherits `RelMinSize * min(Deflection,
/// DeflectionInterior)` (floored at `Confusion`).
fn effective_min_size(parameters: &MeshParameters) -> f64 {
    if parameters.min_size >= CONFUSION {
        parameters.min_size
    } else {
        let min_deflection = parameters.deflection.min(effective_deflection_interior(parameters));
        (MeshParameters::rel_min_size() * min_deflection).max(CONFUSION)
    }
}

/// Resolves `DeflectionInterior`: a value below `Confusion` inherits the base
/// `Deflection`.
fn effective_deflection_interior(parameters: &MeshParameters) -> f64 {
    if parameters.deflection_interior >= CONFUSION {
        parameters.deflection_interior
    } else {
        parameters.deflection
    }
}

/// Resolves `AngleInterior`: a value below `Angular` inherits the base `Angle`.
fn effective_angle_interior(parameters: &MeshParameters) -> f64 {
    if parameters.angle_interior >= ANGULAR {
        parameters.angle_interior
    } else {
        parameters.angle
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpPln, GpSphere};
    use occt_geom::{GeomPlane, GeomSphere};

    /// A self-contained mesh driving the refinement loop against the real
    /// `Delaun` by batch re-triangulation (OCCT inserts nodes incrementally via
    /// `DelaunayNodeInsertionMeshAlgo`; this wrapper rebuilds the triangulation
    /// with the accumulated vertex list, which exercises the same deflection
    /// control logic).
    struct RefineMesh {
        vertices: Vec<DelaunVertex>,
        delaun: Delaun,
    }

    impl RefineMesh {
        fn from_points(points: &[DelaunVertex]) -> Self {
            Self {
                vertices: points.to_vec(),
                delaun: Delaun::new_vertices(points),
            }
        }

        fn insert_control(&mut self, surface: &dyn Surface, nodes: &[GpPnt2d]) -> bool {
            let before = self.vertices.len();
            for &uv in nodes {
                let p3d = surface.d0(uv.x(), uv.y());
                self.vertices.push(DelaunVertex::new(uv, p3d, self.vertices.len() as i32, VertexState::Free));
            }
            if self.vertices.len() == before {
                return false;
            }
            self.delaun = Delaun::new_vertices(&self.vertices);
            true
        }

        fn triangle_count(&self) -> usize {
            self.delaun.result().elements_of_domain().len()
        }

        fn live_vertex_ids(&self) -> Vec<i32> {
            let mut ids = HashSet::new();
            for &id in self.delaun.result().elements_of_domain() {
                let tri = self.delaun.get_triangle(id);
                for v in tri.vertex_indices {
                    ids.insert(v);
                }
            }
            ids.into_iter().collect()
        }
    }

    impl DeflectionMesh for RefineMesh {
        fn elements_of_domain(&self) -> Vec<i32> {
            self.delaun.result().elements_of_domain().iter().copied().collect()
        }
        fn triangle(&self, id: i32) -> DelaunTriangle {
            self.delaun.get_triangle(id)
        }
        fn link(&self, id: i32) -> DelaunLink {
            self.delaun.get_edge(id)
        }
        fn link_movability(&self, id: i32) -> VertexState {
            self.delaun.result().link_movability(id)
        }
        fn vertex(&self, id: i32) -> DelaunVertex {
            self.delaun.get_vertex(id)
        }
    }

    fn plane() -> Arc<dyn Surface> {
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())))
    }

    fn unit_sphere() -> Arc<dyn Surface> {
        Arc::new(GeomSphere::new(GpSphere::new(GpAx3::standard(), 1.0).unwrap()))
    }

    /// 3x3 grid on the plane z = 0, UV in [0, 2] x [0, 2].
    fn plane_grid(surface: &dyn Surface) -> Vec<DelaunVertex> {
        uv_grid(surface, 3, 3, (0.0, 2.0), (0.0, 2.0))
    }

    /// 3x3 grid over the full parametric domain of a unit sphere.
    fn sphere_grid(surface: &dyn Surface) -> Vec<DelaunVertex> {
        uv_grid(surface, 3, 3, (0.0, 2.0 * PI), (-PI / 2.0, PI / 2.0))
    }

    fn uv_grid(
        surface: &dyn Surface,
        nu: usize,
        nv: usize,
        (u0, u1): (f64, f64),
        (v0, v1): (f64, f64),
    ) -> Vec<DelaunVertex> {
        let du = u1 - u0;
        let dv = v1 - v0;
        let mut points = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let uv = GpPnt2d::new(
                    u0 + du * (i as f64 + 0.5) / nu as f64,
                    v0 + dv * (j as f64 + 0.5) / nv as f64,
                );
                points.push(DelaunVertex::new(uv, surface.d0(uv.x(), uv.y()), 0, VertexState::Free));
            }
        }
        points
    }

    #[test]
    fn normal_deviation_is_distance_to_plane() {
        let ref_pnt = GpPnt::new(0.0, 0.0, 0.0);
        let functor = NormalDeviation::new(&ref_pnt, GpVec::new(0.0, 0.0, 1.0));
        // Points in the plane z = 0 have zero deviation.
        assert!(functor.square_deviation(&GpPnt::new(1.0, 2.0, 0.0)) < 1e-15);
        // Points offset by `d` normal to the plane have deviation d^2.
        assert!((functor.square_deviation(&GpPnt::new(1.0, 2.0, 3.0)) - 9.0).abs() < 1e-12);
        assert!((functor.square_deviation(&GpPnt::new(1.0, 2.0, -0.5)) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn line_deviation_is_distance_to_chord() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(2.0, 0.0, 0.0);
        let functor = LineDeviation::new(&a, &b);
        // On the chord.
        assert!(functor.square_deviation(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-15);
        // Off the chord by 1 unit.
        assert!((functor.square_deviation(&GpPnt::new(1.0, 1.0, 0.0)) - 1.0).abs() < 1e-12);
        // Collinear beyond an endpoint: distance to the supporting line, not the
        // segment, so a point on the line still reads zero.
        assert!(functor.square_deviation(&GpPnt::new(3.0, 0.0, 0.0)) < 1e-15);
    }

    #[test]
    fn plane_triangulation_is_not_refined() {
        // A flat plane has zero linear and angular deflection everywhere: the
        // refinement must leave the 3x3 grid untouched.
        let surface = plane();
        let params = MeshParameters::default();

        let points = plane_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();
        assert_eq!(n0, 8, "3x3 grid triangulates to 8 triangles");

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let max_deflection = algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.001,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert_eq!(mesh.triangle_count(), n0, "plane must not refine");
        assert!(max_deflection < 1e-9, "plane max deflection must be ~0, got {max_deflection}");
    }

    #[test]
    fn sphere_triangulation_is_refined() {
        // A coarse sphere triangulation has huge chord sagitta: the refinement
        // must split links and grow the triangle count.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();
        assert_eq!(n0, 8);

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let max_deflection = algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.01,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert!(
            mesh.triangle_count() > n0,
            "sphere must refine: {n0} -> {}",
            mesh.triangle_count()
        );
        assert!(max_deflection > 0.0, "sphere deflection must be positive");
    }

    #[test]
    fn sphere_refined_vertices_lie_on_surface() {
        // Every node generated by the refinement must sit on the unit sphere.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.01,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        let ids = mesh.live_vertex_ids();
        assert!(ids.len() > points.len(), "refinement must add vertices");
        for id in ids {
            let vertex = mesh.delaun.get_vertex(id);
            let radius = vertex.p3d.distance(&GpPnt::zero());
            assert!(
                (radius - 1.0).abs() < 1e-9,
                "vertex off unit sphere: radius {radius} at ({}, {}, {})",
                vertex.p3d.x(),
                vertex.p3d.y(),
                vertex.p3d.z()
            );
        }
    }

    #[test]
    fn large_min_size_blocks_sphere_refinement() {
        // A MinSize of half the sphere diameter (2.0) makes every candidate
        // control point closer than MinSize to an existing node, so all splits
        // are rejected: the triangle count must stay unchanged.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.min_size = 2.0;

        let points = sphere_grid(surface.as_ref());
        let mut mesh = RefineMesh::from_points(&points);
        let n0 = mesh.triangle_count();

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        algo.post_process(
            &mut mesh,
            surface.as_ref(),
            0.001,
            &mut |m, nodes| m.insert_control(surface.as_ref(), nodes),
        );

        assert_eq!(mesh.triangle_count(), n0, "min_size must reject every split");
    }

    #[test]
    fn post_process_reads_mesh_and_proposes_control_nodes() {
        // Drive the algorithm against the real `Delaun` directly: even with a
        // recording (non-mutating) insert callback the curved sphere mesh must
        // produce control nodes and a positive max deflection.
        let surface = unit_sphere();
        let mut params = MeshParameters::default();
        params.deflection = 0.01;
        params.deflection_interior = 0.01;

        let points = sphere_grid(surface.as_ref());
        let mut delaun = Delaun::new_vertices(&points);
        let n0 = delaun.result().elements_of_domain().len();

        let mut algo = DelaunayDeflectionControlMeshAlgo::new(&params);
        let mut collected: Vec<GpPnt2d> = Vec::new();
        let max_deflection = algo.post_process(
            &mut delaun,
            surface.as_ref(),
            0.01,
            &mut |_mesh, nodes| {
                collected.extend_from_slice(nodes);
                !nodes.is_empty()
            },
        );

        assert!(!collected.is_empty(), "sphere must propose control nodes");
        assert!(max_deflection > 0.0);
        // The recording callback does not mutate the mesh.
        assert_eq!(delaun.result().elements_of_domain().len(), n0);
    }
}
