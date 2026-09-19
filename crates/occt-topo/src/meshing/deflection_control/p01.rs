use super::prelude::*;
use super::*;

pub(super) static DBG_CENTER: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_LINK: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_TESTED: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_SHOT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_SHOT0: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_REJMIN: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_MID: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_BAD3D: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);
pub(super) static DBG_MIN_SQ: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(f64::INFINITY.to_bits());

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
    /// Triangle ids whose circumcircle contains `p` in the same 2D space as
    /// stored vertices. `BRepMesh_CircleTool::Select`. `None` falls back to
    /// scanning every live triangle.
    fn shot_triangles(&self, p: GpXY) -> Option<Vec<i32>>;
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
    fn shot_triangles(&self, p: GpXY) -> Option<Vec<i32>> {
        Some(self.circles().select(p))
    }
}

/// Functor computing the squared deflection of a point from a reference plane
/// (a reference point plus a unit normal). Source: `NormalDeviation`.
#[derive(Clone, Copy)]
pub(super) struct NormalDeviation<'a> {
    pub(super) ref_pnt: &'a GpPnt,
    pub(super) normal: GpVec,
}

impl<'a> NormalDeviation<'a> {
    pub(super) fn new(ref_pnt: &'a GpPnt, normal: GpVec) -> Self {
        Self { ref_pnt, normal }
    }

    pub(super) fn square_deviation(&self, point: &GpPnt) -> f64 {
        SquareDeviation::square_deviation(self, point)
    }
}

/// Functor computing the squared deflection of a point from a chord.
/// Source: `LineDeviation`.
#[derive(Clone, Copy)]
pub(super) struct LineDeviation<'a> {
    pub(super) pnt1: &'a GpPnt,
    pub(super) pnt2: &'a GpPnt,
}

impl<'a> LineDeviation<'a> {
    pub(super) fn new(pnt1: &'a GpPnt, pnt2: &'a GpPnt) -> Self {
        Self { pnt1, pnt2 }
    }

    pub(super) fn square_deviation(&self, point: &GpPnt) -> f64 {
        SquareDeviation::square_deviation(self, point)
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
pub(super) struct TriangleNodeInfo {
    /// Position in parametric (UV) space.
    pub(super) point_2d: GpXY,
    /// Associated 3D point on the surface.
    pub(super) point_3d: GpPnt,
    /// Whether the link opposite this node is a frontier (constraint) link.
    pub(super) is_frontier_link: bool,
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
    pub(super) parameters: MeshParameters,
    /// Square of the minimum allowed triangle edge (`MinSize * MinSize`).
    pub(super) sq_min_size: f64,
    /// Face deflection used as the linear threshold.
    pub(super) face_deflection: f64,
    /// Largest squared deflection observed in the current optimization run.
    pub(super) max_sq_deflection: f64,
    /// Becomes `false` as soon as a non-degenerate triangle is processed.
    pub(super) is_all_degenerated: bool,
    /// Control points (2D) accumulated during the current pass.
    pub(super) control_nodes: Vec<GpPnt2d>,
    /// Map of oriented links already processed (dedup across passes).
    pub(super) processed_couples: HashSet<(i32, i32)>,
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
        pub(super) const MAX_PASSES: usize = 11;

        self.face_deflection = face_deflection;
        self.processed_couples.clear();
        self.max_sq_deflection = -1.0;
        self.is_all_degenerated = false;
        DBG_CENTER.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_LINK.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_TESTED.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_SHOT.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_SHOT0.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_REJMIN.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_MID.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_BAD3D.store(0, std::sync::atomic::Ordering::Relaxed);
        DBG_MIN_SQ.store(f64::INFINITY.to_bits(), std::sync::atomic::Ordering::Relaxed);

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

            for id in &ids {
                let triangle = mesh.triangle(*id);
                self.split_triangle_geometry(mesh, surface, &triangle);
            }

            if self.control_nodes.is_empty() {
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
    pub(super) fn split_triangle_geometry<M: DeflectionMesh>(
        &mut self,
        mesh: &M,
        surface: &dyn Surface,
        triangle: &DelaunTriangle,
    ) {
        // OCCT skips Deleted triangles here, but `elements_of_domain` only yields
        // live triangles, so the check is redundant for this port.
        // `hxx:211`: `ElementNodes`, not cached `vertex_indices`.
        let node_indices = element_nodes_from_links(mesh, triangle);

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
            if self.use_point(
                mesh,
                surface,
                GpPnt2d::from_xy(center_2d),
                &NormalDeviation::new(&nodes_info[0].point_3d, normal),
            ) {
                DBG_CENTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            self.split_links(mesh, surface, &nodes_info, &node_indices);
        }
    }

    /// Returns the 2D/3D data of the three nodes of a triangle, together with the
    /// frontier flag of each opposite link. Source: `getTriangleInfo`.
    pub(super) fn get_triangle_info<M: DeflectionMesh>(
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
    pub(super) fn compute_triangle_geometry(
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
    pub(super) fn check_triangle_for_degenerativity_and_get_links(
        &self,
        nodes_info: &[TriangleNodeInfo; 3],
        links: &mut [GpVec; 3],
    ) -> bool {
        pub(super) const MINIMAL_SQ_LENGTH_3D: f64 = 1.0e-12;
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
    pub(super) fn check_triangle_area_2d(&self, nodes_info: &[TriangleNodeInfo; 3]) -> bool {
        pub(super) const MINIMAL_AREA_2D: f64 = 1.0e-9;
        let link_2d_1 = nodes_info[1].point_2d.subtracted(&nodes_info[0].point_2d);
        let link_2d_2 = nodes_info[2].point_2d.subtracted(&nodes_info[1].point_2d);
        link_2d_1.crossed(&link_2d_2).abs() > MINIMAL_AREA_2D
    }

    /// Computes the unit normal from two link vectors.
    /// Returns `false` if the normal has null magnitude. Source: `computeNormal`.
    pub(super) fn compute_normal(&self, link1: &GpVec, link2: &GpVec, normal: &mut GpVec) -> bool {
        let cross = link1.crossed(link2);
        // `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:284`:
        // `if (aNormal.SquareMagnitude() > gp::Resolution())` with
        // `gp::Resolution()` = `RealSmall()` = `DBL_MIN` (`gp.hxx:60`).
        if cross.square_magnitude() > REAL_SMALL {
            *normal = cross.normalized();
            return true;
        }
        false
    }

    /// Checks the deflection of the three link midpoints; points failing the
    /// linear, MinSize and angular checks are queued for insertion.
    /// Source: `splitLinks`.
    pub(super) fn split_links<M: DeflectionMesh>(
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
            DBG_TESTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

            let mid_2d = nodes_info[i]
                .point_2d
                .added(&nodes_info[j].point_2d)
                .multiplied_scalar(0.5);
            let mid_2d = GpPnt2d::from_xy(mid_2d);

            let line = LineDeviation::new(&nodes_info[i].point_3d, &nodes_info[j].point_3d);
            if self.use_point(mesh, surface, mid_2d, &line) {
                DBG_MID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            } else {
                let reject_min_size =
                    self.reject_split_links_for_min_size(surface, &nodes_info[i], &nodes_info[j], mid_2d);
                let reject_angular = self.check_link_ends_for_angular_deviation(
                    surface,
                    &nodes_info[i],
                    &nodes_info[j],
                    mid_2d,
                );
                if !reject_min_size && !reject_angular {
                    DBG_LINK.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    self.control_nodes.push(mid_2d);
                }
            }
        }
    }

    /// Checks that the two links produced by splitting the given link by its
    /// midpoint satisfy the MinSize requirement. Source:
    /// `rejectSplitLinksForMinSize`.
    pub(super) fn reject_split_links_for_min_size(
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
    pub(super) fn check_link_ends_for_angular_deviation(
        &self,
        surface: &dyn Surface,
        node1: &TriangleNodeInfo,
        node2: &TriangleNodeInfo,
        _mid: GpPnt2d,
    ) -> bool {
        let angle_interior = effective_angle_interior(&self.parameters);
        let (s1, n1) = crate::meshing::geomlib_norm::norm_estim(
            surface,
            node1.point_2d.x(),
            node1.point_2d.y(),
            occt_core::precision::CONFUSION,
        );
        let (s2, n2) = crate::meshing::geomlib_norm::norm_estim(
            surface,
            node2.point_2d.x(),
            node2.point_2d.y(),
            occt_core::precision::CONFUSION,
        );
        match (s1, s2, n1, n2) {
            (0, 0, Some(normal1), Some(normal2)) => normal1.angle(&normal2) <= angle_interior,
            // Non-zero `GeomLib::NormEstim` status: accept the link (cxx:358-368).
            _ => true,
        }
    }

    /// Returns the 3D surface point for the given 2D parameter.
    /// Source: `getPoint3d` (via `IMeshData_Face::GetSurface()->D0`).
    pub(super) fn get_point_3d(&self, surface: &dyn Surface, p2d: GpPnt2d) -> GpPnt {
        surface.d0(p2d.x(), p2d.y())
    }

    /// Computes the deflection of a point; caches it for insertion when it
    /// exceeds the tolerance. Returns `true` if the point has been cached.
    /// Source: `usePoint`.
    pub(super) fn use_point<M: DeflectionMesh, D: SquareDeviation>(
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
    pub(super) fn check_deflection_of_point_and_update_cache<M: DeflectionMesh>(
        &mut self,
        mesh: &M,
        p2d: GpPnt2d,
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
        self.reject_by_min_size(mesh, p2d, p3d)
    }

    /// Checks the distance between the given node and the nodes of the triangles
    /// shot by it against the MinSize criteria. Source: `rejectByMinSize`
    /// (`DelaunayDeflectionControlMeshAlgo.hxx:422-454`).
    pub(super) fn reject_by_min_size<M: DeflectionMesh>(
        &self,
        mesh: &M,
        p2d: GpPnt2d,
        p3d: GpPnt,
    ) -> bool {
        let mut used_nodes: HashSet<i32> = HashSet::new();
        let shot = mesh
            .shot_triangles(p2d.coord)
            .unwrap_or_else(|| mesh.elements_of_domain());
        if shot.is_empty() {
            DBG_SHOT0.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
        DBG_SHOT.fetch_add(shot.len(), std::sync::atomic::Ordering::Relaxed);
        for id in shot {
            let triangle = mesh.triangle(id);
            for node_id in element_nodes_from_links(mesh, &triangle) {
                if used_nodes.insert(node_id) {
                    let vertex = mesh.vertex(node_id);
                    let sq = p3d.square_distance(&vertex.p3d);
                    if !(vertex.p3d.coord.x.is_finite()
                        && vertex.p3d.coord.y.is_finite()
                        && vertex.p3d.coord.z.is_finite())
                    {
                        DBG_BAD3D.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    }
                    DBG_MIN_SQ.fetch_min(sq.to_bits(), std::sync::atomic::Ordering::Relaxed);
                    if sq < self.sq_min_size {
                        DBG_REJMIN.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        return true;
                    }
                }
            }
        }
        false
    }
}

/// `BRepMesh_DataStructureOfDelaun::ElementNodes` (`cxx:259-286`): link 0 both
/// ends and link 2 the remaining corner. Used by `splitTriangleGeometry`
/// (`hxx:211`) and `rejectByMinSize` (`hxx:434-435`).
fn element_nodes_from_links<M: DeflectionMesh>(mesh: &M, triangle: &DelaunTriangle) -> [i32; 3] {
    let mut nodes = [0i32; 3];
    let link1 = mesh.link(triangle.link_at(0).abs());
    if triangle.link_at(0) > 0 {
        nodes[0] = link1.first_node();
        nodes[1] = link1.last_node();
    } else {
        nodes[1] = link1.first_node();
        nodes[0] = link1.last_node();
    }
    let link2 = mesh.link(triangle.link_at(2).abs());
    nodes[2] = if triangle.link_at(2) > 0 {
        link2.first_node()
    } else {
        link2.last_node()
    };
    nodes
}

/// Resolves `MinSize` exactly like `BRepMesh_IncrementalMesh::initParameters`:
/// a value below `Confusion` inherits `RelMinSize * min(Deflection,
/// DeflectionInterior)` (floored at `Confusion`).
pub(super) fn effective_min_size(parameters: &MeshParameters) -> f64 {
    if parameters.min_size >= CONFUSION {
        parameters.min_size
    } else {
        let min_deflection = parameters.deflection.min(effective_deflection_interior(parameters));
        (MeshParameters::rel_min_size() * min_deflection).max(CONFUSION)
    }
}

/// Resolves `DeflectionInterior`: a value below `Confusion` inherits the base
/// `Deflection`.
pub(super) fn effective_deflection_interior(parameters: &MeshParameters) -> f64 {
    if parameters.deflection_interior >= CONFUSION {
        parameters.deflection_interior
    } else {
        parameters.deflection
    }
}

/// Resolves `AngleInterior`: a value below `Angular` inherits the base `Angle`.
pub(super) fn effective_angle_interior(parameters: &MeshParameters) -> f64 {
    if parameters.angle_interior >= ANGULAR {
        parameters.angle_interior
    } else {
        parameters.angle
    }
}
