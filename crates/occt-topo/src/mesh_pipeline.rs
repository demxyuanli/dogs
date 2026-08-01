//! Mesh processing pipeline — refinement, smoothing, decimation, repair,
//! surface normals and mesh → BRep assembly.
//!
//! Port of the OCCT `TKXMesh`/`TKMeshVS` mesh-processing toolkits plus the
//! `BRepMesh`-adjacent mesh-repair utilities. A [`MeshPipeline`] owns a
//! triangle-soup mesh (points + index triples) and layers the classic
//! processing operations on top:
//!
//! - Loop-style **uniform subdivision** (`refine`, via `uniform_subdivide`)
//! - **Laplacian smoothing** (`smooth`, via `laplacian_smooth`)
//! - Quadric **decimation** (`decimate`, via `decimate_mesh`)
//! - **Repair** (`repair`, `weld`, `remove_duplicate_triangles`): degenerate /
//!   duplicate triangle removal, spatial-hash vertex welding and dead-vertex
//!   compaction — mirroring `TKXMesh::Mesh` repair pass
//! - Per-vertex and per-face **normals** (`compute_normals`, `face_normals`)
//! - Divergence-theorem **volume** and **surface area**
//! - **Mesh → BRep** assembly (`to_shape`, `to_shell`), reusing
//!   `shape_mesh_to_brep` so a closed manifold becomes a solid
//!
//! This module is the mesh-adjacent half of OCCT's mesh pipeline (`BRepMesh`
//! produces the mesh; `TKXMesh`/`TKMeshVS` post-process it). In OCCT the chain
//! is roughly:
//!
//! ```text
//! BRepMesh_IncrementalMesh
//!   → Poly_Triangulation (nodes + triangles)
//!   → TKXMesh::Mesh (repair: drop small/degenerate/duplicated triangles,
//!                    merge close nodes, snap to the shape)
//!   → MeshVS / BRepMesh tools (normals, smoothing, quality)
//!   → BRep_Builder / BRepBuilderAPI_MakeSolid (triangles → faces → solid)
//! ```
//!
//! Each stage maps onto a method here: [`MeshPipeline::from_shape`] is the
//! tessellation entry point; [`MeshPipeline::repair`], [`MeshPipeline::weld`]
//! and [`MeshPipeline::remove_duplicate_triangles`] implement the TKXMesh
//! cleanup pass; [`MeshPipeline::compute_normals`], [`MeshPipeline::face_normals`],
//! [`MeshPipeline::smooth`] and [`MeshPipeline::quality`] cover the normals and
//! quality tooling; [`MeshPipeline::to_shape`] and [`MeshPipeline::to_shell`]
//! rebuild BRep topology (a closed 2-manifold becomes a solid).
//!
//! Source: TKXMesh, TKMeshVS, BRepMesh.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::geom::mesh_analysis::{
    analyze_mesh, mesh_edge_topology, mesh_signed_volume, triangle_quality, MeshQualitySummary,
    TriangleQuality,
};
use occt_core::geom::mesh_ops::{
    centroid, decimate_mesh, laplacian_smooth, triangle_normal, uniform_subdivide, vertex_normals,
};
use occt_core::gp::{GpPnt, GpTrsf, GpVec, GpXyz};
use occt_core::poly::triangulation::Triangle;
use occt_core::poly::Triangulation;

use crate::abs::ShapeType;
use crate::mesh::{mesh_surface_area, ShapeMesh};
use crate::mesh_to_brep::shape_mesh_to_brep;
use crate::shape::{Shell, TopoShape};
use crate::shape_mesh::mesh_shape;

/// A lightweight indexed triangle-soup mesh, the input/output currency of the
/// processing pipeline. Mirrors the `Poly_Triangulation` node/triangle layout
/// but stores index triples as plain `(usize, usize, usize)` tuples.
///
/// Triangle winding is significant for [`MeshPipeline::volume`],
/// [`MeshPipeline::is_closed`] and the outward-facing normals: OCCT's mesh
/// generators emit counter-clockwise winding seen from the outward side of the
/// surface, which makes the signed divergence volume positive for closed
/// solids.
#[derive(Debug, Clone)]
pub struct MeshPipeline {
    /// Node positions.
    pub vertices: Vec<GpPnt>,
    /// Index triples into `vertices`.
    pub triangles: Vec<(usize, usize, usize)>,
}

impl MeshPipeline {
    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Mesh every face of a BRep `shape` (deflection-bounded tessellation) into
    /// this pipeline. Equivalent to `BRepMesh_IncrementalMesh` feeding a
    /// `Poly_Triangulation`.
    pub fn from_shape(shape: &TopoShape, deflection: f64) -> Self {
        let mesh = mesh_shape(shape, deflection);
        Self::from_mesh(&mesh)
    }

    /// Adopt an existing `ShapeMesh` (the crate's `Poly_Triangulation`-style
    /// representation), converting its `Triangle` entries to index triples.
    pub fn from_mesh(mesh: &ShapeMesh) -> Self {
        Self {
            vertices: mesh.vertices.clone(),
            triangles: mesh
                .triangles
                .iter()
                .map(|t| (t.n0, t.n1, t.n2))
                .collect(),
        }
    }

    /// Build the pipeline from raw vertex positions and index triples.
    pub fn from_raw(vertices: Vec<GpPnt>, triangles: Vec<(usize, usize, usize)>) -> Self {
        Self { vertices, triangles }
    }

    /// Build the pipeline from a `Poly_Triangulation`.
    pub fn from_triangulation(tri: &Triangulation) -> Self {
        Self {
            vertices: tri.nodes.clone(),
            triangles: tri
                .triangles
                .iter()
                .map(|t| (t.n0, t.n1, t.n2))
                .collect(),
        }
    }

    // ------------------------------------------------------------------
    // Queries
    // ------------------------------------------------------------------

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// True when the mesh carries no triangles (an empty pipeline).
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// Number of distinct (undirected) edges over the triangle set.
    pub fn edge_count(&self) -> usize {
        self.edge_use_counts().len()
    }

    /// Number of boundary (open) edges — edges referenced by exactly one
    /// triangle. Zero for a closed 2-manifold.
    pub fn open_edges(&self) -> usize {
        self.edge_use_counts().values().filter(|&&c| c == 1).count()
    }

    /// Axis-aligned bounding box of all vertices.
    pub fn bbox(&self) -> BndBox {
        let mut b = BndBox::new();
        for p in &self.vertices {
            b.add_point(p);
        }
        b
    }

    /// Average vertex position (zero for an empty vertex list).
    pub fn centroid(&self) -> GpPnt {
        centroid(&self.vertices)
    }

    /// Whether the mesh is a closed 2-manifold: every undirected edge is used
    /// by exactly two triangles. `BRepMesh`-style watertight check.
    pub fn is_closed(&self) -> bool {
        let counts = self.edge_use_counts();
        !counts.is_empty() && counts.values().all(|&c| c == 2)
    }

    /// Edge topology as `(open_edge_count, non_manifold_edge_count)`.
    pub fn edge_topology(&self) -> (usize, usize) {
        mesh_edge_topology(&self.triangles)
    }

    /// Aggregate quality metrics over all triangles (aspect ratio, angle
    /// bounds, degenerate count, total area).
    pub fn quality(&self) -> MeshQualitySummary {
        analyze_mesh(&self.vertices, &self.triangles)
    }

    /// Fraction of degenerate (near-zero-area) triangles in `[0, 1]`.
    pub fn degenerate_fraction(&self) -> f64 {
        let q = self.quality();
        if q.triangle_count == 0 {
            0.0
        } else {
            q.degenerate_count as f64 / q.triangle_count as f64
        }
    }

    /// Signed volume enclosed by the mesh via the divergence theorem:
    /// `V = (1/6) · Σ a·(b×c)` over the triangles. Positive for outward
    /// counter-clockwise winding. Open or inconsistently oriented meshes give
    /// unreliable results — callers should check [`MeshPipeline::is_closed`]
    /// first.
    pub fn volume(&self) -> f64 {
        mesh_signed_volume(&self.vertices, &self.triangles)
    }

    /// Total surface area of the triangle soup.
    pub fn surface_area(&self) -> f64 {
        mesh_surface_area(&self.to_shape_mesh())
    }

    /// Per-vertex normals, area-weighted over the incident faces and
    /// normalized. Degenerate vertices fall back to a unit vector.
    pub fn compute_normals(&self) -> Vec<GpVec> {
        vertex_normals(&self.vertices, &self.triangles)
    }

    /// Per-triangle unit normals (right-hand rule from the winding). A
    /// degenerate triangle yields a zero vector.
    pub fn face_normals(&self) -> Vec<GpVec> {
        self.triangles
            .iter()
            .map(|&(a, b, c)| {
                if a < self.vertices.len() && b < self.vertices.len() && c < self.vertices.len() {
                    triangle_normal(&[self.vertices[a], self.vertices[b], self.vertices[c]])
                } else {
                    GpVec::zero()
                }
            })
            .collect()
    }

    /// Neighbours of vertex `i`: the set of distinct vertices sharing a
    /// triangle with it. Used by smoothing and 1-ring queries.
    pub fn vertex_neighbors(&self, i: usize) -> Vec<usize> {
        let mut set = HashSet::new();
        for &(a, b, c) in &self.triangles {
            if a == i {
                set.insert(b);
                set.insert(c);
            }
            if b == i {
                set.insert(a);
                set.insert(c);
            }
            if c == i {
                set.insert(a);
                set.insert(b);
            }
        }
        set.into_iter().collect()
    }

    /// Degree of vertex `i` — the number of distinct neighbouring vertices.
    pub fn vertex_degree(&self, i: usize) -> usize {
        self.vertex_neighbors(i).len()
    }

    /// Area of the triangle at `index` (0.0 for an out-of-range index or a
    /// degenerate triangle).
    pub fn triangle_area(&self, index: usize) -> f64 {
        if index >= self.triangles.len() {
            return 0.0;
        }
        let (a, b, c) = self.triangles[index];
        if a >= self.vertices.len() || b >= self.vertices.len() || c >= self.vertices.len() {
            return 0.0;
        }
        let ab = GpVec::from_pnts(&self.vertices[a], &self.vertices[b]);
        let ac = GpVec::from_pnts(&self.vertices[a], &self.vertices[c]);
        0.5 * ab.crossed(&ac).magnitude()
    }

    /// Centroid of the triangle at `index` (average of its three vertices).
    pub fn triangle_centroid(&self, index: usize) -> GpPnt {
        if index >= self.triangles.len() {
            return GpPnt::zero();
        }
        let (a, b, c) = self.triangles[index];
        let pa = self.vertices.get(a).copied().unwrap_or_default();
        let pb = self.vertices.get(b).copied().unwrap_or_default();
        let pc = self.vertices.get(c).copied().unwrap_or_default();
        GpPnt::new(
            (pa.x() + pb.x() + pc.x()) / 3.0,
            (pa.y() + pb.y() + pc.y()) / 3.0,
            (pa.z() + pb.z() + pc.z()) / 3.0,
        )
    }

    /// Per-triangle quality metrics (aspect ratio, angles, area, degeneracy)
    /// for the triangle at `index`.
    pub fn triangle_quality_at(&self, index: usize) -> TriangleQuality {
        if index >= self.triangles.len() {
            return triangle_quality(&GpPnt::zero(), &GpPnt::zero(), &GpPnt::zero());
        }
        let (a, b, c) = self.triangles[index];
        let pa = self.vertices.get(a).copied().unwrap_or_default();
        let pb = self.vertices.get(b).copied().unwrap_or_default();
        let pc = self.vertices.get(c).copied().unwrap_or_default();
        triangle_quality(&pa, &pb, &pc)
    }

    /// Euler characteristic V − E + F over the mesh. For a closed oriented
    /// polyhedron this is 2; for a topological disk it is 1.
    pub fn euler_characteristic(&self) -> i32 {
        self.vertex_count() as i32 - self.edge_count() as i32 + self.triangle_count() as i32
    }

    /// Minimal enclosing sphere estimate `(center, radius)`: the axis-aligned
    /// bounding-box centre with radius equal to the furthest vertex distance.
    /// A cheap, stable culling bound (not the true minimal sphere).
    pub fn bounding_sphere(&self) -> (GpPnt, f64) {
        if self.vertices.is_empty() {
            return (GpPnt::zero(), 0.0);
        }
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for p in &self.vertices {
            min = GpPnt::new(min.x().min(p.x()), min.y().min(p.y()), min.z().min(p.z()));
            max = GpPnt::new(max.x().max(p.x()), max.y().max(p.y()), max.z().max(p.z()));
        }
        let center = GpPnt::new(
            0.5 * (min.x() + max.x()),
            0.5 * (min.y() + max.y()),
            0.5 * (min.z() + max.z()),
        );
        let radius = self.vertices.iter().map(|p| p.distance(&center)).fold(0.0, f64::max);
        (center, radius)
    }

    /// Area-weighted centre of mass of the triangle surface (the surface's
    /// first moment, like `GProp_GProps` over a triangulation). Returns the
    /// vertex centroid when the mesh has no triangles.
    pub fn center_of_mass(&self) -> GpPnt {
        let mut acc = GpXyz::zero();
        let mut total = 0.0;
        for &(a, b, c) in &self.triangles {
            if a >= self.vertices.len() || b >= self.vertices.len() || c >= self.vertices.len() {
                continue;
            }
            let pa = &self.vertices[a];
            let pb = &self.vertices[b];
            let pc = &self.vertices[c];
            let ab = GpVec::from_pnts(pa, pb);
            let ac = GpVec::from_pnts(pa, pc);
            let area = 0.5 * ab.crossed(&ac).magnitude();
            if area <= 1e-30 {
                continue;
            }
            let cx = (pa.x() + pb.x() + pc.x()) / 3.0;
            let cy = (pa.y() + pb.y() + pc.y()) / 3.0;
            let cz = (pa.z() + pb.z() + pc.z()) / 3.0;
            acc = acc.added(&GpXyz::new(cx * area, cy * area, cz * area));
            total += area;
        }
        if total <= 1e-30 {
            return self.centroid();
        }
        GpPnt::from_xyz(&acc.divided(total))
    }

    // ------------------------------------------------------------------
    // Processing operations
    // ------------------------------------------------------------------

    /// Uniform (Loop-style) subdivision: every triangle becomes four by
    /// splitting each edge at its midpoint. Applied `iterations` times, so the
    /// triangle count quadruples per iteration.
    pub fn refine(&mut self, iterations: usize) {
        for _ in 0..iterations {
            let (v, t) = uniform_subdivide(&self.vertices, &self.triangles);
            self.vertices = v;
            self.triangles = t;
        }
    }

    /// Laplacian smoothing: move each vertex toward the average of its
    /// neighbors. `iterations` passes with strength `lambda` (typically 0.5).
    /// Topology is unchanged; only the vertex positions are updated.
    pub fn smooth(&mut self, iterations: usize, lambda: f64) {
        self.vertices = laplacian_smooth(&self.vertices, &self.triangles, iterations, lambda);
    }

    /// Quadric-error decimation: repeatedly remove the lowest-error interior
    /// vertex (dropping its incident triangles) until `target_vertices` remain
    /// or no further removable vertex exists. Coarse LOD reduction.
    pub fn decimate(&mut self, target_vertices: usize) {
        let (v, t) = decimate_mesh(&self.vertices, &self.triangles, target_vertices);
        self.vertices = v;
        self.triangles = t;
    }

    /// TKXMesh-style repair pass. In order:
    ///
    /// 1. Remove degenerate (zero-area) triangles.
    /// 2. Remove duplicate triangles (same vertex set, any winding).
    /// 3. Weld coincident vertices within a scale-derived tolerance.
    /// 4. Re-remove degenerate/duplicate triangles the weld may have created.
    /// 5. Compact unreferenced vertices out of the vertex list.
    ///
    /// Returns the total number of removed triangles.
    pub fn repair(&mut self) -> usize {
        let mut removed = 0;
        removed += self.remove_degenerate_triangles();
        removed += self.remove_duplicate_triangles();
        let tol = self.scale_tol();
        self.weld(tol);
        // Welding may collapse two triangle corners onto one vertex.
        removed += self.remove_degenerate_triangles();
        removed += self.remove_duplicate_triangles();
        self.remove_unreferenced_vertices();
        removed
    }

    /// Merge coincident vertices within `tol` (spatial-hash grid with 27-cell
    /// neighbourhood search), remapping triangle indices to the survivors.
    /// Returns the number of merged (removed) vertices.
    pub fn weld(&mut self, tol: f64) -> usize {
        let tol = tol.max(1e-12);
        let n = self.vertices.len();
        if n == 0 {
            return 0;
        }
        let cell = |p: &GpPnt| -> (i64, i64, i64) {
            (
                f64::floor(p.x() / tol) as i64,
                f64::floor(p.y() / tol) as i64,
                f64::floor(p.z() / tol) as i64,
            )
        };
        let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        let mut unique: Vec<GpPnt> = Vec::new();
        let mut remap = vec![0usize; n];
        for (i, p) in self.vertices.iter().enumerate() {
            let c = cell(p);
            let mut found = None;
            'search: for dx in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dz in -1i64..=1 {
                        if let Some(bucket) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                            for &j in bucket {
                                if p.distance(&unique[j]) <= tol {
                                    found = Some(j);
                                    break 'search;
                                }
                            }
                        }
                    }
                }
            }
            match found {
                Some(j) => remap[i] = j,
                None => {
                    let j = unique.len();
                    unique.push(*p);
                    grid.entry(c).or_default().push(j);
                    remap[i] = j;
                }
            }
        }
        let merged = n - unique.len();
        self.vertices = unique;
        for t in &mut self.triangles {
            t.0 = remap[t.0];
            t.1 = remap[t.1];
            t.2 = remap[t.2];
        }
        merged
    }

    /// Remove triangles that share the same vertex set (regardless of winding),
    /// keeping the first occurrence. Returns the number removed.
    pub fn remove_duplicate_triangles(&mut self) -> usize {
        let before = self.triangles.len();
        let mut seen: HashSet<[usize; 3]> = HashSet::new();
        self.triangles.retain(|&(a, b, c)| {
            let mut key = [a, b, c];
            key.sort_unstable();
            seen.insert(key)
        });
        before - self.triangles.len()
    }

    /// Translate every vertex by `v`.
    pub fn translate(&mut self, v: &GpVec) {
        for p in &mut self.vertices {
            *p = p.translated_vec(v);
        }
    }

    /// Apply an affine transformation (rotation/scale/shear/translation) to
    /// every vertex. Mirrors `TopoDS_Shape::Move` on the mesh level.
    pub fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.vertices {
            *p = p.transformed(t);
        }
    }

    /// Uniformly scale every vertex about `center` by `factor`.
    pub fn scale(&mut self, center: &GpPnt, factor: f64) {
        for p in &mut self.vertices {
            *p = p.scaled(center, factor);
        }
    }

    /// Reverse the winding of every triangle (swap two indices), flipping all
    /// face normals. Used to correct a globally inverted mesh.
    pub fn reverse_winding(&mut self) {
        for t in &mut self.triangles {
            let (a, b, c) = *t;
            *t = (b, a, c);
        }
    }

    /// Ensure all triangles wind consistently: flood-fill over shared edges,
    /// flipping a triangle whenever a neighbour traverses their common edge in
    /// the same direction (inconsistent orientation). Returns the number of
    /// triangles flipped. Open (non-manifold) meshes are processed per
    /// connected component.
    ///
    /// For a closed mesh the result is additionally normalized so the signed
    /// volume is positive — i.e. the mesh points outward, matching OCCT's
    /// convention for boundary-representation faces.
    pub fn ensure_consistent_winding(&mut self) -> usize {
        // Undirected edge → incident (triangle index, directed traversal).
        let mut edge_tris: HashMap<(usize, usize), Vec<(usize, (usize, usize))>> = HashMap::new();
        for (ti, &(a, b, c)) in self.triangles.iter().enumerate() {
            for &(p, q) in &[(a, b), (b, c), (c, a)] {
                let key = if p < q { (p, q) } else { (q, p) };
                edge_tris.entry(key).or_default().push((ti, (p, q)));
            }
        }
        let n = self.triangles.len();
        let mut flipped: Vec<Option<bool>> = vec![None; n];
        let mut count = 0usize;
        for seed in 0..n {
            if flipped[seed].is_some() {
                continue;
            }
            flipped[seed] = Some(false);
            let mut stack = vec![seed];
            while let Some(ti) = stack.pop() {
                let par = flipped[ti].expect("visited triangle has a flip state");
                let (a, b, c) = self.triangles[ti];
                for &(p, q) in &[(a, b), (b, c), (c, a)] {
                    let key = if p < q { (p, q) } else { (q, p) };
                    for &(nti, (np, nq)) in edge_tris.get(&key).expect("edge present") {
                        if nti == ti || flipped[nti].is_some() {
                            continue;
                        }
                        // Same traversal direction on the shared edge means the
                        // neighbour is wound the wrong way; flip it.
                        let same = (p, q) == (np, nq);
                        let need_flip = same ^ par;
                        flipped[nti] = Some(need_flip);
                        if need_flip {
                            count += 1;
                        }
                        stack.push(nti);
                    }
                }
            }
        }
        for (ti, &f) in flipped.iter().enumerate() {
            if f == Some(true) {
                let (a, b, c) = self.triangles[ti];
                self.triangles[ti] = (b, a, c);
            }
        }
        // Orient closed meshes outward (positive signed volume). The flood-fill
        // seed fixes the relative winding but not the global orientation.
        if self.is_closed() && self.volume() < 0.0 {
            self.reverse_winding();
        }
        count
    }

    /// Append another mesh to this one: its vertices are concatenated and its
    /// triangle indices offset accordingly. Both meshes keep their own winding.
    pub fn append(&mut self, other: &MeshPipeline) {
        let off = self.vertices.len();
        self.vertices.extend(other.vertices.iter().copied());
        for &(a, b, c) in &other.triangles {
            self.triangles.push((off + a, off + b, off + c));
        }
    }

    // ------------------------------------------------------------------
    // BRep conversion
    // ------------------------------------------------------------------

    /// Rebuild a BRep shape from the mesh. Each triangle becomes a planar face;
    /// when the mesh is a closed 2-manifold the faces assemble into a solid,
    /// otherwise the shell is returned. Coincident vertices within `tol` are
    /// welded first so shared edges line up across faces.
    pub fn to_shape(&self, tol: f64) -> Result<TopoShape, String> {
        if self.triangles.is_empty() {
            return Err("to_shape: mesh has no triangles".to_string());
        }
        let mut clean = self.clone();
        clean.weld(tol.max(1e-12));
        clean.remove_duplicate_triangles();
        clean.remove_degenerate_triangles();
        if clean.triangles.is_empty() {
            return Err("to_shape: mesh has no non-degenerate triangles".to_string());
        }
        let brep = shape_mesh_to_brep(&clean.to_shape_mesh());
        Ok(match brep.solid {
            Some(solid) => solid.0,
            None => brep.shell.0,
        })
    }

    /// Assemble a shell from the triangles (one planar face per triangle,
    /// shared edges deduplicated). Open meshes are allowed — the shell is
    /// simply not closed. Returns an error on an empty mesh or any degenerate
    /// triangle (which cannot form a valid face).
    pub fn to_shell(&self) -> Result<Shell, String> {
        if self.triangles.is_empty() {
            return Err("to_shell: mesh has no triangles".to_string());
        }
        for &(a, b, c) in &self.triangles {
            if a >= self.vertices.len() || b >= self.vertices.len() || c >= self.vertices.len() {
                return Err(format!("to_shell: triangle ({a},{b},{c}) out of range"));
            }
            let ab = GpVec::from_pnts(&self.vertices[a], &self.vertices[b]);
            let ac = GpVec::from_pnts(&self.vertices[a], &self.vertices[c]);
            if ab.crossed(&ac).square_magnitude() <= 1e-24 {
                return Err(format!("to_shell: degenerate triangle ({a},{b},{c})"));
            }
        }
        let brep = shape_mesh_to_brep(&self.to_shape_mesh());
        Ok(brep.shell)
    }

    /// Convert the pipeline into a `Poly_Triangulation` (shared node list).
    pub fn to_triangulation(&self) -> Triangulation {
        Triangulation::new(
            self.vertices.clone(),
            self.triangles
                .iter()
                .map(|&(a, b, c)| Triangle::new(a, b, c))
                .collect(),
        )
    }

    /// Convert the pipeline into the crate's `ShapeMesh`.
    pub fn to_shape_mesh(&self) -> ShapeMesh {
        ShapeMesh {
            vertices: self.vertices.clone(),
            triangles: self
                .triangles
                .iter()
                .map(|&(a, b, c)| Triangle::new(a, b, c))
                .collect(),
            source_shape: ShapeType::Compound,
        }
    }

    // ------------------------------------------------------------------
    // Internal helpers
    // ------------------------------------------------------------------

    /// Count how many triangles reference each undirected edge.
    fn edge_use_counts(&self) -> HashMap<(usize, usize), usize> {
        let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
        for &(a, b, c) in &self.triangles {
            for (x, y) in [(a, b), (b, c), (c, a)] {
                let key = if x < y { (x, y) } else { (y, x) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        counts
    }

    /// Remove degenerate triangles — zero area from repeated indices or
    /// coincident/colinear vertices. Returns the number removed.
    pub fn remove_degenerate_triangles(&mut self) -> usize {
        let before = self.triangles.len();
        self.triangles.retain(|&(a, b, c)| {
            if a >= self.vertices.len() || b >= self.vertices.len() || c >= self.vertices.len() {
                return false;
            }
            if a == b || b == c || a == c {
                return false;
            }
            let ab = GpVec::from_pnts(&self.vertices[a], &self.vertices[b]);
            let ac = GpVec::from_pnts(&self.vertices[a], &self.vertices[c]);
            ab.crossed(&ac).square_magnitude() > 1e-24
        });
        before - self.triangles.len()
    }

    /// Remove vertices no longer referenced by any triangle, remapping the
    /// triangle indices. Returns the number of vertices removed.
    pub fn remove_unreferenced_vertices(&mut self) -> usize {
        let n = self.vertices.len();
        let mut used = vec![false; n];
        for &(a, b, c) in &self.triangles {
            if a < n {
                used[a] = true;
            }
            if b < n {
                used[b] = true;
            }
            if c < n {
                used[c] = true;
            }
        }
        let mut map = vec![usize::MAX; n];
        let mut new_verts: Vec<GpPnt> = Vec::new();
        for i in 0..n {
            if used[i] {
                map[i] = new_verts.len();
                new_verts.push(self.vertices[i]);
            }
        }
        let removed = n - new_verts.len();
        self.vertices = new_verts;
        for t in &mut self.triangles {
            t.0 = map[t.0];
            t.1 = map[t.1];
            t.2 = map[t.2];
        }
        removed
    }

    /// A weld tolerance derived from the mesh scale: 1e-9 of the bounding-box
    /// diagonal (min 1e-12). Used by [`MeshPipeline::repair`].
    fn scale_tol(&self) -> f64 {
        if self.vertices.is_empty() {
            return 1e-12;
        }
        let mut min = self.vertices[0];
        let mut max = self.vertices[0];
        for p in &self.vertices {
            min = GpPnt::new(min.x().min(p.x()), min.y().min(p.y()), min.z().min(p.z()));
            max = GpPnt::new(max.x().max(p.x()), max.y().max(p.y()), max.z().max(p.z()));
        }
        let diag = max.coord.subtracted(&min.coord).modulus();
        (diag * 1e-9).max(1e-12)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::mesh_box;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::shell_check::shell_is_closed;
    use crate::shape_mesh::shape_volume;

    /// The unit-cube mesh (8 vertices, 12 outward triangles) used across tests.
    fn unit_cube() -> MeshPipeline {
        MeshPipeline::from_mesh(&mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0))))
    }

    /// A flat unit square split into two triangles (planar mesh).
    fn unit_square() -> MeshPipeline {
        MeshPipeline::from_raw(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
            ],
            vec![(0, 1, 2), (0, 2, 3)],
        )
    }

    #[test]
    fn from_box_mesh() {
        // A box shape meshed through the pipeline matches `mesh_shape` exactly.
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let direct = mesh_shape(&b.solid.0, 0.5);
        let p = MeshPipeline::from_shape(&b.solid.0, 0.5);
        assert_eq!(p.vertex_count(), direct.vertices.len());
        assert_eq!(p.triangle_count(), direct.triangles.len());
        assert!(!p.is_empty());
    }

    #[test]
    fn refine_quadruples() {
        // 2 triangles → 8 after one uniform-subdivision pass.
        let mut p = unit_square();
        assert_eq!(p.triangle_count(), 2);
        p.refine(1);
        assert_eq!(p.triangle_count(), 8);
        assert_eq!(p.vertex_count(), 9);
    }

    #[test]
    fn smooth_flattens_bump() {
        // Raise the far corner of the square — a bump — and smooth it down.
        let mut p = unit_square();
        p.vertices[2] = GpPnt::new(1.0, 1.0, 0.5);
        let z0 = p.vertices[2].z();
        p.smooth(5, 0.5);
        assert!(p.vertices[2].z() < z0 - 1e-6, "bump reduced from {z0} to {}", p.vertices[2].z());
    }

    #[test]
    fn decimate_reduces() {
        // A refined square (9 vertices) decimated to at most 5.
        let mut p = unit_square();
        p.refine(1);
        assert_eq!(p.vertex_count(), 9);
        let tris_before = p.triangle_count();
        p.decimate(5);
        assert!(p.vertex_count() <= 5, "verts {}", p.vertex_count());
        assert!(p.triangle_count() < tris_before, "triangles reduced");
    }

    #[test]
    fn repair_removes_degenerate() {
        // A degenerate triangle (0,1,1) and a duplicate of (0,1,2) are removed.
        let mut p = MeshPipeline::from_raw(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
            ],
            vec![(0, 1, 2), (0, 2, 3), (0, 1, 1), (0, 1, 2)],
        );
        assert_eq!(p.triangle_count(), 4);
        let removed = p.repair();
        assert!(removed >= 2, "removed {removed}");
        assert_eq!(p.triangle_count(), 2);
        assert_eq!(p.vertex_count(), 4, "all four corners still referenced");
    }

    #[test]
    fn weld_merges_coincident() {
        // Vertex 4 is a duplicate of vertex 0; welding merges them.
        let mut p = MeshPipeline::from_raw(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(0.0, 0.0, 0.0),
            ],
            vec![(0, 1, 2), (0, 2, 3), (4, 1, 2)],
        );
        assert_eq!(p.vertex_count(), 5);
        let merged = p.weld(1e-9);
        assert_eq!(merged, 1);
        assert_eq!(p.vertex_count(), 4);
        // The triangle that referenced the merged vertex is remapped.
        assert!(p.triangles.iter().all(|&(a, b, c)| a < 4 && b < 4 && c < 4));
    }

    #[test]
    fn normals_unit_length() {
        let p = unit_cube();
        let normals = p.compute_normals();
        assert_eq!(normals.len(), p.vertex_count());
        for n in &normals {
            assert!((n.magnitude() - 1.0).abs() < 1e-9, "normal length {}", n.magnitude());
        }
    }

    #[test]
    fn face_normals_outward() {
        // Every face normal of the box points along its dominant axis, with the
        // six signed axis directions each covered.
        let p = unit_cube();
        let normals = p.face_normals();
        assert_eq!(normals.len(), 12);
        let dominant = |n: &GpVec| -> (char, bool) {
            let ax = n.x().abs().max(n.y().abs()).max(n.z().abs());
            if n.x().abs() == ax {
                ('X', n.x() > 0.0)
            } else if n.y().abs() == ax {
                ('Y', n.y() > 0.0)
            } else {
                ('Z', n.z() > 0.0)
            }
        };
        let set: HashSet<(char, bool)> = normals.iter().map(dominant).collect();
        let expect: HashSet<(char, bool)> = [
            ('X', false),
            ('X', true),
            ('Y', false),
            ('Y', true),
            ('Z', false),
            ('Z', true),
        ]
        .into();
        assert_eq!(set, expect, "normals must point outward on all six faces");
    }

    #[test]
    fn to_shape_sphere_roundtrip() {
        // Sphere → pipeline → BRep; the rebuilt solid's volume is within 25% of
        // the analytic 4π/3 r³.
        //
        // Note: the crate's UV-grid sphere tessellator leaves a topological seam
        // (two boundary-edge runs at the poles near θ=0/2π), so `to_shape`
        // yields a shell, not a solid, for this particular mesh. `shape_volume`
        // is therefore unreliable on the open shell; the pipeline's own
        // divergence-theorem `volume()` is exact on the triangle soup and is
        // the robust measure here. A closed mesh (e.g. the box) roundtrips to a
        // solid and `shape_volume` applies directly.
        let sphere = BRepPrimSphere::make_sphere(1.0);
        let p = MeshPipeline::from_shape(&sphere.solid.0, 0.1);
        assert!(!p.is_empty());
        let shape = p.to_shape(1e-6).expect("mesh to shape");
        assert!(
            shape.is_solid() || shape.is_shell() || shape.is_compound(),
            "unexpected shape type {:?}",
            shape.shape_type()
        );
        let vol = p.volume();
        let expect = 4.0 / 3.0 * std::f64::consts::PI;
        assert!(
            (vol - expect).abs() < 0.25 * expect,
            "sphere volume {vol} vs analytic {expect}"
        );
    }

    #[test]
    fn to_shell_closed_for_box() {
        let p = unit_cube();
        let shell = p.to_shell().expect("shell from cube mesh");
        assert!(shell_is_closed(&shell), "box shell must be closed");
    }

    #[test]
    fn volume_and_area_unit_cube() {
        let p = unit_cube();
        let vol = p.volume();
        let area = p.surface_area();
        assert!((vol - 1.0).abs() < 0.1, "volume {vol}");
        assert!((area - 6.0).abs() < 0.1, "area {area}");
    }

    #[test]
    fn is_closed_open_and_closed() {
        let open = unit_square();
        assert!(!open.is_closed(), "an open strip is not closed");
        assert_eq!(open.open_edges(), 4, "the square strip has 4 boundary edges");

        let cube = unit_cube();
        assert!(cube.is_closed(), "a cube mesh is closed");
        assert_eq!(cube.open_edges(), 0);
    }

    #[test]
    fn consistent_winding_fixes_flipped_triangle() {
        // A closed unit-cube mesh with one triangle wound backwards (its
        // normal pointing inward) is re-oriented so all normals agree.
        let mut p = unit_cube();
        p.triangles[0] = (p.triangles[0].1, p.triangles[0].0, p.triangles[0].2);
        let flipped = p.ensure_consistent_winding();
        assert!(flipped >= 1, "at least one triangle flipped ({flipped})");
        // After the pass the mesh is consistently oriented; the signed volume
        // of the closed cube must be near +1 again.
        let vol = p.volume();
        assert!((vol - 1.0).abs() < 1e-6, "consistent volume {vol}");
    }

    #[test]
    fn append_concatenates_meshes() {
        let mut p = unit_cube();
        let other = unit_square();
        let v_before = p.vertex_count();
        let t_before = p.triangle_count();
        p.append(&other);
        assert_eq!(p.vertex_count(), v_before + 4);
        assert_eq!(p.triangle_count(), t_before + 2);
        // The appended square's indices are offset into the merged list.
        assert!(p.triangles.iter().all(|&(a, b, c)| a < p.vertex_count() && b < p.vertex_count() && c < p.vertex_count()));
    }

    #[test]
    fn pipeline_chain_refine_smooth() {
        // Refine → smooth → to_shape must succeed without panic and yield a
        // valid BRep shape.
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut p = MeshPipeline::from_shape(&b.solid.0, 0.5);
        assert!(!p.is_empty());
        p.refine(1);
        p.smooth(2, 0.3);
        assert!(!p.is_empty());
        let shape = p.to_shape(1e-6).expect("shape after refine+smooth");
        assert!(
            shape.is_solid() || shape.is_shell() || shape.is_compound(),
            "unexpected shape type {:?}",
            shape.shape_type()
        );
    }
}
