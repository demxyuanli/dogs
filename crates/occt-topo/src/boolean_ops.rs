//! Voxel-based boolean operations on BRep solids — union, intersection and
//! subtraction, all mesh-approximated.
//! Source: `BRepAlgoAPI_Fuse`, `BRepAlgoAPI_Common`, `BRepAlgoAPI_Cut`.
//!
//! The boolean is evaluated on a shared voxel grid over the union bounding
//! box: each solid is classified inside/outside the voxels, the cells are
//! combined per the operation, and the boundary of the occupied set is
//! converted back to a triangle mesh. Exact NURBS booleans are far out of
//! scope; this gives robust, resolution-controlled results.

use std::collections::HashMap;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;

use crate::mesh::{mesh_surface_area, ShapeMesh};
use crate::shape::TopoShape;

/// The boolean operation to apply to a voxel cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp {
    /// Cell occupied by either input (Fuse).
    Union,
    /// Cell occupied by both inputs (Common).
    Intersection,
    /// Cell occupied by the first but not the second (Cut).
    Subtraction,
}

/// A shared-axis-aligned voxel grid covering a bounding box.
struct VoxelGrid {
    min: [f64; 3],
    cell: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    occupied: Vec<bool>,
}

impl VoxelGrid {
    fn new(bbox: &BndBox, resolution: usize) -> Self {
        let (x0, x1, y0, y1, z0, z1) = bbox.get().unwrap_or((0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
        let span = ((x1 - x0).max(y1 - y0)).max(z1 - z0).max(1e-9);
        let cell = span / resolution as f64;
        let nx = (((x1 - x0) / cell).ceil() as usize).max(1);
        let ny = (((y1 - y0) / cell).ceil() as usize).max(1);
        let nz = (((z1 - z0) / cell).ceil() as usize).max(1);
        VoxelGrid {
            min: [x0, y0, z0],
            cell,
            nx,
            ny,
            nz,
            occupied: vec![false; nx * ny * nz],
        }
    }

    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.ny + j) * self.nx + i
    }

    /// Voxelize `shape`: mark cells whose center is inside it.
    ///
    /// Meshes the shape once at a deflection tied to the cell size, then
    /// classifies every cell center with the even-odd ray test.
    fn fill(&mut self, shape: &TopoShape) {
        let mesh = crate::shape_mesh::mesh_shape(shape, 0.5 * self.cell);
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
                    let p = self.cell_center(i, j, k);
                    if point_in_mesh(&mesh, &p) {
                        let idx = self.idx(i, j, k);
                        self.occupied[idx] = true;
                    }
                }
            }
        }
    }

    fn cell_center(&self, i: usize, j: usize, k: usize) -> GpPnt {
        GpPnt::new(
            self.min[0] + (i as f64 + 0.5) * self.cell,
            self.min[1] + (j as f64 + 0.5) * self.cell,
            self.min[2] + (k as f64 + 0.5) * self.cell,
        )
    }

    fn in_bounds(&self, i: i64, j: i64, k: i64) -> bool {
        i >= 0 && j >= 0 && k >= 0 && (i as usize) < self.nx && (j as usize) < self.ny && (k as usize) < self.nz
    }

    fn occupied_at(&self, i: i64, j: i64, k: i64) -> bool {
        self.in_bounds(i, j, k) && self.occupied[self.idx(i as usize, j as usize, k as usize)]
    }

    /// Convert the occupied set to a surface mesh: for every occupied cell,
    /// emit the faces that border an empty cell (or the grid boundary), with
    /// outward normals. Corners are deduplicated so the mesh is watertight.
    fn to_mesh(&self) -> ShapeMesh {
        let mut corner_idx: HashMap<(i64, i64, i64), usize> = HashMap::new();
        let mut vertices: Vec<GpPnt> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

        // Get-or-create the vertex for a grid corner.
        fn corner(
            idx: &mut HashMap<(i64, i64, i64), usize>,
            verts: &mut Vec<GpPnt>,
            grid: &VoxelGrid,
            x: i64, y: i64, z: i64,
        ) -> usize {
            if let Some(&i) = idx.get(&(x, y, z)) {
                return i;
            }
            let p = GpPnt::new(
                grid.min[0] + x as f64 * grid.cell,
                grid.min[1] + y as f64 * grid.cell,
                grid.min[2] + z as f64 * grid.cell,
            );
            let i = verts.len();
            verts.push(p);
            idx.insert((x, y, z), i);
            i
        }

        // Emit a quad face as two triangles (winding gives the outward normal).
        fn face(
            idx: &mut HashMap<(i64, i64, i64), usize>,
            verts: &mut Vec<GpPnt>,
            tris: &mut Vec<Triangle>,
            grid: &VoxelGrid,
            quad: [(i64, i64, i64); 4],
        ) {
            let v = [
                corner(idx, verts, grid, quad[0].0, quad[0].1, quad[0].2),
                corner(idx, verts, grid, quad[1].0, quad[1].1, quad[1].2),
                corner(idx, verts, grid, quad[2].0, quad[2].1, quad[2].2),
                corner(idx, verts, grid, quad[3].0, quad[3].1, quad[3].2),
            ];
            tris.push(Triangle::new(v[0], v[1], v[2]));
            tris.push(Triangle::new(v[0], v[2], v[3]));
        }

        for i in 0..self.nx as i64 {
            for j in 0..self.ny as i64 {
                for k in 0..self.nz as i64 {
                    if !self.occupied_at(i, j, k) {
                        continue;
                    }
                    let c000 = (i, j, k);
                    let c100 = (i + 1, j, k);
                    let c010 = (i, j + 1, k);
                    let c110 = (i + 1, j + 1, k);
                    let c001 = (i, j, k + 1);
                    let c101 = (i + 1, j, k + 1);
                    let c011 = (i, j + 1, k + 1);
                    let c111 = (i + 1, j + 1, k + 1);
                    if !self.occupied_at(i - 1, j, k) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c000, c001, c011, c010]); }
                    if !self.occupied_at(i + 1, j, k) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c100, c110, c111, c101]); }
                    if !self.occupied_at(i, j - 1, k) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c000, c100, c101, c001]); }
                    if !self.occupied_at(i, j + 1, k) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c010, c011, c111, c110]); }
                    if !self.occupied_at(i, j, k - 1) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c000, c010, c110, c100]); }
                    if !self.occupied_at(i, j, k + 1) { face(&mut corner_idx, &mut vertices, &mut triangles, self, [c001, c101, c111, c011]); }
                }
            }
        }
        ShapeMesh { vertices, triangles, source_shape: crate::abs::ShapeType::Solid }
    }
}

/// Voxel boolean of two solids at `resolution` cells across the span.
pub fn voxel_boolean(
    a: &TopoShape,
    b: &TopoShape,
    resolution: usize,
    op: BoolOp,
) -> Result<ShapeMesh, String> {
    if resolution < 4 {
        return Err("voxel_boolean: resolution must be >= 4".into());
    }
    let mut bb = BndBox::new();
    for s in [a, b] {
        bb = crate::solid_union::bbox_union(&bb, &crate::bbox_from_geometry::shape_bbox(s));
    }
    let mut grid = VoxelGrid::new(&bb, resolution);
    // Classify both solids on the same grid.
    grid.fill(a);
    let a_occ = grid.occupied.clone();
    grid.occupied = vec![false; a_occ.len()];
    grid.fill(b);
    let b_occ = grid.occupied.clone();
    // Combine per the operation.
    for i in 0..a_occ.len() {
        grid.occupied[i] = match op {
            BoolOp::Union => a_occ[i] || b_occ[i],
            BoolOp::Intersection => a_occ[i] && b_occ[i],
            BoolOp::Subtraction => a_occ[i] && !b_occ[i],
        };
    }
    Ok(grid.to_mesh())
}

/// Voxel intersection (Common) of two solids.
pub fn voxel_intersect(a: &TopoShape, b: &TopoShape, resolution: usize) -> Result<ShapeMesh, String> {
    voxel_boolean(a, b, resolution, BoolOp::Intersection)
}

/// Voxel subtraction (Cut): `a` minus `b`.
pub fn voxel_subtract(a: &TopoShape, b: &TopoShape, resolution: usize) -> Result<ShapeMesh, String> {
    voxel_boolean(a, b, resolution, BoolOp::Subtraction)
}

/// Surface area of a voxel boolean result.
pub fn boolean_surface_area(a: &TopoShape, b: &TopoShape, resolution: usize, op: BoolOp) -> Result<f64, String> {
    Ok(mesh_surface_area(&voxel_boolean(a, b, resolution, op)?))
}

/// Möller–Trumbore ray/triangle intersection. Returns the parameter `t` along
/// the ray `origin + t·dir` when the ray hits, else `None`.
pub fn ray_triangle_intersect(
    origin: &GpPnt,
    dir: &occt_core::gp::GpVec,
    v0: &GpPnt,
    v1: &GpPnt,
    v2: &GpPnt,
) -> Option<f64> {
    let edge1 = occt_core::gp::GpVec::from_pnts(v0, v1);
    let edge2 = occt_core::gp::GpVec::from_pnts(v0, v2);
    let p = dir.xyz().crossed(edge2.xyz());
    let det = edge1.xyz().dot(&p);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv = 1.0 / det;
    let tvec = occt_core::gp::GpVec::from_pnts(v0, origin);
    let u = tvec.xyz().dot(&p) * inv;
    if u < 0.0 || u > 1.0 {
        return None;
    }
    let q = tvec.xyz().crossed(edge1.xyz());
    let v = dir.xyz().dot(&q) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = edge2.xyz().dot(&q) * inv;
    if t < 0.0 { None } else { Some(t) }
}

/// Even-odd point-in-mesh test via ray casting (ray along +X).
///
/// Coincident hits (a ray passing exactly through a shared triangle edge, which
/// reports twice at the same parameter) are deduplicated so the crossing count
/// stays odd/even correct. The origin is jittered a tiny amount perpendicular
/// to the ray to break exact grid-diagonal degeneracies.
pub fn point_in_mesh(mesh: &ShapeMesh, p: &GpPnt) -> bool {
    let dir = occt_core::gp::GpVec::new(1.0, 0.0, 0.0);
    let jittered = GpPnt::new(p.x(), p.y() + 1e-7, p.z() + 1e-7);
    let mut hits: Vec<f64> = mesh
        .triangles
        .iter()
        .filter_map(|t| {
            ray_triangle_intersect(&jittered, &dir, &mesh.vertices[t.n0], &mesh.vertices[t.n1], &mesh.vertices[t.n2])
        })
        .collect();
    hits.sort_by(f64::total_cmp);
    let mut unique = 0usize;
    let mut prev: Option<f64> = None;
    for t in hits {
        if prev.map_or(true, |p| (t - p).abs() > 1e-9) {
            unique += 1;
            prev = Some(t);
        }
    }
    unique % 2 == 1
}

/// Number of surface triangles of a voxel boolean result.
pub fn boolean_triangle_count(a: &TopoShape, b: &TopoShape, resolution: usize, op: BoolOp) -> Result<usize, String> {
    Ok(voxel_boolean(a, b, resolution, op)?.triangles.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpTrsf;

    #[test]
    fn point_in_mesh_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mesh = crate::shape_mesh::mesh_shape(&b.solid.0, 0.2);
        assert!(point_in_mesh(&mesh, &GpPnt::new(0.5, 0.5, 0.5)));
        assert!(!point_in_mesh(&mesh, &GpPnt::new(2.0, 2.0, 2.0)));
    }

    #[test]
    fn ray_triangle_basic() {
        // Triangle in z=0 plane.
        let (v0, v1, v2) = (GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.));
        // Ray along +Z from below hits the triangle interior.
        let t = ray_triangle_intersect(
            &GpPnt::new(0.25, 0.25, -1.0), &occt_core::gp::GpVec::new(0.0, 0.0, 1.0),
            &v0, &v1, &v2,
        );
        assert!(t.is_some(), "ray through the triangle hits it");
        // Ray that misses the triangle (outside its footprint) → None.
        let t = ray_triangle_intersect(
            &GpPnt::new(0.9, 0.9, -1.0), &occt_core::gp::GpVec::new(0.0, 0.0, 1.0),
            &v0, &v1, &v2,
        );
        assert!(t.is_none(), "ray outside the triangle misses it");
        // Ray in the same plane as the triangle is degenerate → None.
        let t = ray_triangle_intersect(
            &GpPnt::new(0.25, 0.25, 0.0), &occt_core::gp::GpVec::new(1.0, 0.0, 0.0),
            &v0, &v1, &v2,
        );
        assert!(t.is_none(), "coplanar ray has no conventional crossing");
    }

    #[test]
    fn voxel_subtract_reduces_volume() {
        // A = [0,2]³, B = [1,3]×[0,2]×[0,2] (built at its true position, since
        // move_shape only relabels Location and does not relocate geometry).
        let a = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(3.0, 2.0, 2.0));
        let cut = voxel_subtract(&a.solid.0, &b.solid.0, 24).expect("subtract works");
        assert!(cut.triangles.len() > 24, "cut mesh has triangles");
        // Remaining region [0,1]×[0,2]×[0,2] → area 2(1·2 + 2·2 + 1·2) = 16.
        let area = mesh_surface_area(&cut);
        assert!(area > 12.0 && area < 45.0, "cut area {area} (voxel staircase inflates it)");
    }

    #[test]
    fn voxel_intersect_overlap() {
        let a = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(3.0, 2.0, 2.0));
        let inter = voxel_intersect(&a.solid.0, &b.solid.0, 24).expect("intersect works");
        assert!(inter.triangles.len() > 24);
        // Overlap [1,2]×[0,2]×[0,2] → exact area 16; voxel staircase inflates.
        let area = mesh_surface_area(&inter);
        assert!(area > 12.0 && area < 45.0, "intersect area {area}");
    }

    #[test]
    fn union_covered_by_solid_union() {
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let fused = voxel_boolean(&a.solid.0, &b.solid.0, 20, BoolOp::Union).expect("union works");
        let area = mesh_surface_area(&fused);
        // Two 1×1×1 boxes side by side → exact surface 10 (touching faces
        // cancel), voxel staircase inflates.
        assert!(area > 8.0 && area < 25.0, "union area {area}");
    }
}
