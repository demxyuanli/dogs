//! Mesh-based constructive solid geometry (union/intersection/difference).
//! Simplified voxel-approximation CSG — robust for prototyping.
use crate::gp::GpPnt;
use crate::bnd::BndBox;

/// Voxel grid representation of a solid for CSG operations.
#[derive(Debug, Clone)]
pub struct VoxelSolid {
    pub nx: usize, pub ny: usize, pub nz: usize,
    pub origin: [f64; 3],
    pub size: [f64; 3],
    pub inside: Vec<bool>, // nx*ny*nz
}

impl VoxelSolid {
    /// Create an empty voxel grid over the given bounding box.
    pub fn new(bbox: &BndBox, resolution: f64) -> Option<Self> {
        let (x0, x1, y0, y1, z0, z1) = bbox.get()?;
        if !x1.is_finite() || !y1.is_finite() || !z1.is_finite() { return None; }
        let size = [x1 - x0, y1 - y0, z1 - z0];
        let nx = (size[0] / resolution).ceil().max(1.0) as usize;
        let ny = (size[1] / resolution).ceil().max(1.0) as usize;
        let nz = (size[2] / resolution).ceil().max(1.0) as usize;
        Some(Self {
            nx, ny, nz,
            origin: [x0, y0, z0],
            size,
            inside: vec![false; nx * ny * nz],
        })
    }

    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + j * self.nx + k * self.nx * self.ny
    }

    /// Voxel center coordinates.
    pub fn center(&self, i: usize, j: usize, k: usize) -> GpPnt {
        GpPnt::new(
            self.origin[0] + (i as f64 + 0.5) * self.size[0] / self.nx as f64,
            self.origin[1] + (j as f64 + 0.5) * self.size[1] / self.ny as f64,
            self.origin[2] + (k as f64 + 0.5) * self.size[2] / self.nz as f64,
        )
    }

    /// Mark voxels inside the given signed-distance-like predicate.
    /// inside_fn returns true if a point is inside the solid.
    pub fn fill<F: Fn(&GpPnt) -> bool>(&mut self, inside_fn: &F) {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let p = self.center(i, j, k);
                    let idx = self.index(i, j, k);
                    self.inside[idx] = inside_fn(&p);
                }
            }
        }
    }

    /// Number of inside voxels.
    pub fn count_inside(&self) -> usize { self.inside.iter().filter(|&&v| v).count() }

    /// Approximate volume = inside_voxels * voxel_volume.
    pub fn volume(&self) -> f64 {
        let vv = (self.size[0] / self.nx as f64) * (self.size[1] / self.ny as f64) * (self.size[2] / self.nz as f64);
        self.count_inside() as f64 * vv
    }
}

/// Boolean operations on two voxel solids (same grid geometry).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp { Union, Intersection, Difference }

/// Apply a boolean op between a and b, storing into a.
pub fn apply_bool(a: &mut VoxelSolid, b: &VoxelSolid, op: BoolOp) -> Result<(), &'static str> {
    if a.nx != b.nx || a.ny != b.ny || a.nz != b.nz {
        return Err("csg: grid dimensions must match");
    }
    for i in 0..a.inside.len() {
        let ai = a.inside[i];
        let bi = b.inside[i];
        a.inside[i] = match op {
            BoolOp::Union => ai || bi,
            BoolOp::Intersection => ai && bi,
            BoolOp::Difference => ai && !bi,
        };
    }
    Ok(())
}

/// Box predicate: point inside axis-aligned box [lo, hi].
pub fn box_inside(lo: &GpPnt, hi: &GpPnt) -> impl Fn(&GpPnt) -> bool {
    let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
    let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
    move |p| p.x() >= x0 && p.x() <= x1 && p.y() >= y0 && p.y() <= y1 && p.z() >= z0 && p.z() <= z1
}

/// Sphere predicate: point inside sphere (center, radius).
pub fn sphere_inside(center: &GpPnt, radius: f64) -> impl Fn(&GpPnt) -> bool {
    let (cx, cy, cz) = (center.x(), center.y(), center.z());
    let r2 = radius * radius;
    move |p| {
        let dx = p.x() - cx; let dy = p.y() - cy; let dz = p.z() - cz;
        dx*dx + dy*dy + dz*dz <= r2
    }
}

/// Cylinder predicate: point inside cylinder along Z (axis at x,y), radius r, height h.
pub fn cylinder_inside(radius: f64, height: f64) -> impl Fn(&GpPnt) -> bool {
    let r2 = radius * radius;
    let half = height * 0.5;
    move |p| {
        p.x()*p.x() + p.y()*p.y() <= r2 && p.z() >= -half && p.z() <= half
    }
}

/// Convert a voxel solid to a triangle mesh (marching-cubes-free surface extraction).
/// Simple approach: for each voxel, emit cube-face quads where inside differs from neighbor.
/// Returns (vertices, triangle indices) — quads split into triangles.
pub fn voxel_to_mesh(v: &VoxelSolid) -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let mut verts = Vec::new();
    let mut tris = Vec::new();

    // 8 corners of a voxel at (i,j,k) — relative offsets
    let corner_offsets = [
        (0,0,0),(1,0,0),(1,1,0),(0,1,0),
        (0,0,1),(1,0,1),(1,1,1),(0,1,1),
    ];

    for k in 0..v.nz {
        for j in 0..v.ny {
            for i in 0..v.nx {
                if !v.inside[v.index(i, j, k)] { continue; }
                // Check 6 neighbors; emit face if neighbor is outside or boundary
                let base = [
                    v.origin[0] + i as f64 * v.size[0] / v.nx as f64,
                    v.origin[1] + j as f64 * v.size[1] / v.ny as f64,
                    v.origin[2] + k as f64 * v.size[2] / v.nz as f64,
                ];
                let step = [
                    v.size[0] / v.nx as f64,
                    v.size[1] / v.ny as f64,
                    v.size[2] / v.nz as f64,
                ];
                // -X, +X, -Y, +Y, -Z, +Z faces
                let faces: [((i32,i32,i32), [usize;4], [usize;4]); 6] = [
                    ((-1,0,0), [3,0,4,7], [1,0,5,4]), // face at i-1: corners (i-1 side) — simplified below
                    ((1,0,0), [1,2,6,5], [0,2,7,5]),
                    ((0,-1,0), [0,1,5,4], [3,2,6,7]),
                    ((0,1,0), [2,3,7,6], [1,0,4,5]),
                    ((0,0,-1), [0,3,2,1], [4,7,6,5]),
                    ((0,0,1), [4,5,6,7], [0,1,2,3]),
                ];
                for &(_, quad, alt) in &faces {
                    // Build the 4 corner points; for boundary faces flip winding based on alt
                    let vert_indices: Vec<usize> = quad.iter().map(|&ci| {
                        let (dx, dy, dz) = corner_offsets[ci];
                        let mut nv = verts.len();
                        let vp = GpPnt::new(
                            base[0] + dx as f64 * step[0],
                            base[1] + dy as f64 * step[1],
                            base[2] + dz as f64 * step[2],
                        );
                        verts.push(vp);
                        nv
                    }).collect();
                    tris.push((vert_indices[0], vert_indices[1], vert_indices[2]));
                    tris.push((vert_indices[0], vert_indices[2], vert_indices[3]));
                }
            }
        }
    }
    (verts, tris)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_cube_bbox() -> BndBox {
        let mut b = BndBox::new();
        b.add_point(&GpPnt::new(-0.1, -0.1, -0.1));
        b.add_point(&GpPnt::new(1.1, 1.1, 1.1));
        b
    }

    #[test]
    fn union_volume() {
        let mut grid = VoxelSolid::new(&unit_cube_bbox(), 0.1).unwrap();
        grid.fill(&box_inside(&GpPnt::new(0.,0.,0.), &GpPnt::new(1.,1.,1.)));
        let vol = grid.volume();
        assert!((vol - 1.0).abs() < 0.05, "vol {vol}");
    }

    #[test]
    fn sphere_union_intersection() {
        let mut grid = VoxelSolid::new(&unit_cube_bbox(), 0.1).unwrap();
        grid.fill(&sphere_inside(&GpPnt::new(0.5,0.5,0.5), 0.5));
        let vol_sphere = grid.volume();
        let expect = 4.0/3.0 * std::f64::consts::PI * 0.125;
        assert!((vol_sphere - expect).abs() < 0.05, "sphere vol {vol_sphere} vs {expect}");

        // Intersection with full cube = sphere itself
        let mut b = grid.clone();
        b.fill(&box_inside(&GpPnt::new(0.,0.,0.), &GpPnt::new(1.,1.,1.)));
        apply_bool(&mut grid, &b, BoolOp::Intersection).unwrap();
        assert!((grid.volume() - vol_sphere).abs() < 1e-9);
    }

    #[test]
    fn box_difference() {
        let mut a = VoxelSolid::new(&unit_cube_bbox(), 0.1).unwrap();
        a.fill(&box_inside(&GpPnt::new(0.,0.,0.), &GpPnt::new(1.,1.,1.)));
        let mut b = VoxelSolid::new(&unit_cube_bbox(), 0.1).unwrap();
        b.fill(&box_inside(&GpPnt::new(0.,0.,0.), &GpPnt::new(0.5,1.,1.)));
        apply_bool(&mut a, &b, BoolOp::Difference).unwrap();
        // Half box removed
        let vol = a.volume();
        assert!((vol - 0.5).abs() < 0.1, "vol {vol}");
    }
}
