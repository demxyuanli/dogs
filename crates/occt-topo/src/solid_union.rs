//! Solid union — mesh concatenation, voxel boolean fusion, bbox merge.
//!
//! Source: `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_BooleanOperation` (TKBO).
//!
//! `union_meshes` is a pure triangle-soup merge (no stitching). `voxel_union`
//! rasterizes both solids on a shared grid via even-odd ray casting and
//! produces the fused surface mesh. Both are deliberate approximations —
//! ponytail: `voxel_union` classifies points by ray casting against the faces'
//! (untrimmed) surfaces, so points can be misclassified near a face's phantom
//! boundary; upgrade to a real boundary-aware boolean (BRepAlgoAPI_Fuse) later.

use occt_core::bnd::BndBox;
use occt_core::geom::csg::{apply_bool, voxel_to_mesh, BoolOp, VoxelSolid};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::poly::triangulation::Triangle;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::mesh::ShapeMesh;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;

use crate::face_face::{finite_window, is_plane_like, plane_geometry};

/// Concatenate two meshes, rebasing `b`'s vertex indices. This is a soup
/// merge — coincident vertices are not welded and coplanar triangles are not
/// stitched.
pub fn union_meshes(a: &ShapeMesh, b: &ShapeMesh) -> ShapeMesh {
    let base = a.vertices.len();
    let mut vertices = a.vertices.clone();
    vertices.extend(b.vertices.iter().copied());
    let mut triangles = a.triangles.clone();
    triangles.extend(
        b.triangles
            .iter()
            .map(|t| Triangle::new(t.n0 + base, t.n1 + base, t.n2 + base)),
    );
    ShapeMesh {
        vertices,
        triangles,
        source_shape: ShapeType::Solid,
    }
}

/// Merge two bounding boxes.
pub fn bbox_union(a: &BndBox, b: &BndBox) -> BndBox {
    let mut out = *a;
    out.add_box(b);
    out
}

fn all_faces(shape: &TopoShape) -> Vec<TopoShape> {
    fn walk(s: &TopoShape, out: &mut Vec<TopoShape>) {
        for child in s
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .map(|h| TopoShape::from_handle(h.clone()))
        {
            match child.shape_type() {
                ShapeType::Face => out.push(child),
                ShapeType::Shell | ShapeType::Solid | ShapeType::Compound | ShapeType::CompSolid => walk(&child, out),
                _ => {}
            }
        }
    }
    let mut out = Vec::new();
    walk(shape, &mut out);
    out
}

fn all_edges(shape: &TopoShape) -> Vec<TopoShape> {
    fn walk(s: &TopoShape, out: &mut Vec<TopoShape>) {
        for child in s
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .map(|h| TopoShape::from_handle(h.clone()))
        {
            match child.shape_type() {
                ShapeType::Edge => out.push(child),
                ShapeType::Wire | ShapeType::Face | ShapeType::Shell | ShapeType::Solid
                | ShapeType::Compound | ShapeType::CompSolid => walk(&child, out),
                _ => {}
            }
        }
    }
    let mut out = Vec::new();
    walk(shape, &mut out);
    out
}

/// Bounding box of a shape: exact when it has evaluable edges (their endpoint
/// points), otherwise a fallback sampling of the faces' surfaces.
fn shape_bbox(shape: &TopoShape) -> BndBox {
    let mut bbox = BndBox::new();
    let mut any = false;
    for e in all_edges(shape) {
        let edge = Edge(e);
        if let Some((p1, p2)) = BRepTool::edge_vertices(&edge) {
            bbox.add_point(&p1);
            bbox.add_point(&p2);
            any = true;
        }
    }
    if any {
        return bbox;
    }
    for f in all_faces(shape) {
        if let Some(s) = GeometryRegistry::global().face_surface(&f) {
            let (u0, u1) = finite_window(s.u_range(), 1.0);
            let (v0, v1) = finite_window(s.v_range(), 1.0);
            for i in 0..=8 {
                let u = u0 + (u1 - u0) * i as f64 / 8.0;
                for j in 0..=8 {
                    let v = v0 + (v1 - v0) * j as f64 / 8.0;
                    bbox.add_point(&s.d0(u, v));
                }
            }
        }
    }
    bbox
}

/// Ray–triangle intersection (Möller–Trumbore). Returns the hit distance `t`
/// for a ray `o + t·d`, `t > 0`.
fn ray_triangle(o: &GpPnt, d: &GpVec, a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<f64> {
    let e1 = GpVec::from_pnts(a, b);
    let e2 = GpVec::from_pnts(a, c);
    let pvec = d.crossed(&e2);
    let det = e1.dot(&pvec);
    if det.abs() <= 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    let tvec = GpVec::from_pnts(a, o);
    let u = tvec.dot(&pvec) * inv;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let qvec = tvec.crossed(&e1);
    let v = d.dot(&qvec) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(&qvec) * inv;
    if t > 1e-9 {
        Some(t)
    } else {
        None
    }
}

/// Count crossings of ray `o + t·(+X)` with a surface. Planes are solved
/// exactly; other surfaces are ray-cast against a coarse sampled quad mesh.
fn ray_surface_crossings(o: &GpPnt, dir: &GpVec, s: &dyn occt_geom::Surface) -> usize {
    if is_plane_like(s) {
        let (q, n) = plane_geometry(s);
        let denom = n.dot(dir);
        if denom.abs() <= 1e-12 {
            return 0;
        }
        // n·(q − p) / (n·d) — t > 0 means the crossing is ahead of the ray.
        let t = n.dot(&GpVec::from_pnts(o, &q)) / denom;
        return if t > 1e-9 { 1 } else { 0 };
    }
    // Generic: sample the surface into a grid mesh and ray-cast the triangles.
    let (u0, u1) = finite_window(s.u_range(), 4.0);
    let (v0, v1) = finite_window(s.v_range(), 4.0);
    let (nu, nv) = (24usize, 12usize);
    let mut grid = vec![GpPnt::zero(); (nu + 1) * (nv + 1)];
    for i in 0..=nu {
        let u = u0 + (u1 - u0) * i as f64 / nu as f64;
        for j in 0..=nv {
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            grid[i + j * (nu + 1)] = s.d0(u, v);
        }
    }
    let mut hits = 0usize;
    for i in 0..nu {
        for j in 0..nv {
            let a = i + j * (nu + 1);
            let b = a + 1;
            let c = a + (nu + 1);
            let d = c + 1;
            if ray_triangle(o, dir, &grid[a], &grid[b], &grid[c]).is_some() {
                hits += 1;
            }
            if ray_triangle(o, dir, &grid[b], &grid[d], &grid[c]).is_some() {
                hits += 1;
            }
        }
    }
    hits
}

/// Even-odd ray casting: is `p` inside `shape`? Points outside the shape's
/// bounding box are rejected before casting (avoids phantom crossings from the
/// untrimmed face surfaces).
pub fn point_in_shape(shape: &TopoShape, p: &GpPnt) -> bool {
    let bbox = shape_bbox(shape);
    if !bbox.is_void() && bbox.is_out(p) {
        return false;
    }
    let dir = GpVec::new(1.0, 0.0, 0.0);
    let mut crossings = 0usize;
    for f in all_faces(shape) {
        if let Some(s) = GeometryRegistry::global().face_surface(&f) {
            crossings += ray_surface_crossings(p, &dir, &*s);
        }
    }
    crossings % 2 == 1
}

/// Rasterize both solids on a shared voxel grid (over the union of their
/// bounding boxes) and fuse them with a boolean union. Returns the fused grid.
fn voxel_union_grid(a: &TopoShape, b: &TopoShape, resolution: usize) -> Result<VoxelSolid, String> {
    if resolution == 0 {
        return Err("voxel_union: resolution must be > 0".into());
    }
    let mut bbox = shape_bbox(a);
    bbox.add_box(&shape_bbox(b));
    if bbox.is_void() {
        return Err("voxel_union: empty shape".into());
    }
    let (x0, x1, y0, y1, z0, z1) = bbox
        .get()
        .ok_or_else(|| "voxel_union: empty bounding box".to_string())?;
    if !(x0.is_finite() && x1.is_finite() && y0.is_finite() && y1.is_finite() && z0.is_finite() && z1.is_finite()) {
        return Err("voxel_union: unbounded shape".into());
    }
    let longest = (x1 - x0).max(y1 - y0).max(z1 - z0);
    if longest <= 0.0 {
        return Err("voxel_union: degenerate bounding box".into());
    }
    let cell = longest / resolution as f64;

    let mut va = VoxelSolid::new(&bbox, cell).ok_or("voxel_union: failed to create grid")?;
    let mut vb = VoxelSolid::new(&bbox, cell).ok_or("voxel_union: failed to create grid")?;
    va.fill(&|p: &GpPnt| point_in_shape(a, p));
    vb.fill(&|p: &GpPnt| point_in_shape(b, p));
    apply_bool(&mut va, &vb, BoolOp::Union).map_err(|e| e.to_string())?;
    Ok(va)
}

/// Voxel-approximated boolean union of two solids as a surface mesh.
pub fn voxel_union(a: &TopoShape, b: &TopoShape, resolution: usize) -> Result<ShapeMesh, String> {
    let grid = voxel_union_grid(a, b, resolution)?;
    let (verts, tris) = voxel_to_mesh(&grid);
    Ok(ShapeMesh {
        vertices: verts,
        triangles: tris.into_iter().map(|(i, j, k)| Triangle::new(i, j, k)).collect(),
        source_shape: ShapeType::Solid,
    })
}

/// Volume of the voxel union: inside-voxel count × cell volume.
pub fn union_volume(a: &TopoShape, b: &TopoShape, resolution: usize) -> f64 {
    match voxel_union_grid(a, b, resolution) {
        Ok(grid) => grid.volume(),
        Err(_) => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use crate::builder::TopoBuilder;
    use crate::mesh::{mesh_box, mesh_surface_area};
    use crate::shape::{Face, Solid, Vertex};
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
    use occt_geom::{GeomLine, GeomPlane};

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn corners(lo: &GpPnt, hi: &GpPnt) -> [GpPnt; 8] {
        let (x0, y0, z0) = (lo.x(), lo.y(), lo.z());
        let (x1, y1, z1) = (hi.x(), hi.y(), hi.z());
        [
            GpPnt::new(x0, y0, z0),
            GpPnt::new(x1, y0, z0),
            GpPnt::new(x1, y1, z0),
            GpPnt::new(x0, y1, z0),
            GpPnt::new(x0, y0, z1),
            GpPnt::new(x1, y0, z1),
            GpPnt::new(x1, y1, z1),
            GpPnt::new(x0, y1, z1),
        ]
    }

    const EDGE_PAIRS: [(usize, usize); 12] = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ];

    fn plane_face(b: &TopoBuilder, origin: GpPnt, normal: GpDir) -> GpPln {
        let x_dir = if normal.x().abs() > 0.9 {
            GpDir::new(0.0, 1.0, 0.0).unwrap()
        } else if normal.y().abs() > 0.9 {
            GpDir::new(0.0, 0.0, 1.0).unwrap()
        } else {
            GpDir::new(1.0, 0.0, 0.0).unwrap()
        };
        let ax3 = GpAx3::new(origin, normal, &x_dir).unwrap();
        GpPln::new(ax3)
    }

    /// Axis-aligned box solid: 12 segment edges + 6 planar faces.
    fn make_box_solid(b: &TopoBuilder, lo: &GpPnt, hi: &GpPnt) -> Solid {
        let c = corners(lo, hi);
        let verts: Vec<Vertex> = c.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        let edges: Vec<Edge> = EDGE_PAIRS
            .iter()
            .map(|&(i, j)| {
                let dir = GpDir::from_vec(&GpVec::from_pnts(&c[i], &c[j])).unwrap();
                let lin = GpLin::from_pnt_dir(c[i], dir);
                let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, c[i].distance(&c[j]));
                b.add(&mut e.0, &verts[i].0);
                b.add(&mut e.0, &verts[j].0);
                e
            })
            .collect();
        let normals = [
            GpDir::new(-1.0, 0.0, 0.0).unwrap(),
            GpDir::new(1.0, 0.0, 0.0).unwrap(),
            GpDir::new(0.0, -1.0, 0.0).unwrap(),
            GpDir::new(0.0, 1.0, 0.0).unwrap(),
            GpDir::new(0.0, 0.0, -1.0).unwrap(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
        ];
        let origins = [c[0], c[1], c[0], c[3], c[0], c[4]];
        let face_edges = [
            [8, 7, 11, 3],
            [9, 5, 10, 1],
            [0, 9, 4, 8],
            [2, 10, 6, 11],
            [3, 2, 1, 0],
            [4, 5, 6, 7],
        ];
        let mut faces = Vec::new();
        for i in 0..6 {
            let pln = plane_face(b, origins[i], normals[i]);
            let wire_edges: Vec<Edge> = face_edges[i].iter().map(|&e| edges[e].clone()).collect();
            let wire = b.make_wire(&wire_edges);
            faces.push(b.make_face(Arc::new(GeomPlane::new(pln)), &[wire]));
        }
        let shell = b.make_shell(&faces);
        b.make_solid(&[shell])
    }

    #[test]
    fn union_meshes_concatenates_and_rebases() {
        let a = mesh_box((GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 1.0)));
        let b = mesh_box((GpPnt::new(2.0, 2.0, 2.0), GpPnt::new(3.0, 3.0, 3.0)));
        let u = union_meshes(&a, &b);
        assert_eq!(u.vertices.len(), 16);
        assert_eq!(u.triangles.len(), 24);
        // b's first triangle indices are rebased past a's 8 vertices.
        assert!(u.triangles[12].n0 >= 8);
        assert!(u.triangles[23].n2 < 16);
        // Surface area is additive (disjoint meshes).
        assert!((mesh_surface_area(&u) - 12.0).abs() < 1e-9);
    }

    #[test]
    fn bbox_union_merges() {
        let a = BndBox::from_corners(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b = BndBox::from_corners(&GpPnt::new(2.0, 2.0, 2.0), &GpPnt::new(3.0, 3.0, 3.0));
        let u = bbox_union(&a, &b);
        assert!(u.corner_min().is_equal(&GpPnt::new(0.0, 0.0, 0.0)));
        assert!(u.corner_max().is_equal(&GpPnt::new(3.0, 3.0, 3.0)));

        let c = BndBox::from_corners(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(2.5, 2.5, 2.5));
        let u2 = bbox_union(&a, &c);
        assert!(u2.corner_min().is_equal(&GpPnt::new(0.0, 0.0, 0.0)));
        assert!(u2.corner_max().is_equal(&GpPnt::new(2.5, 2.5, 2.5)));
    }

    #[test]
    fn point_in_shape_inside_and_outside() {
        let b = TopoBuilder::new();
        let solid = make_box_solid(&b, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        assert!(point_in_shape(&solid.0, &GpPnt::new(0.5, 0.5, 0.5)));
        assert!(!point_in_shape(&solid.0, &GpPnt::new(2.0, 2.0, 2.0)));
        assert!(!point_in_shape(&solid.0, &GpPnt::new(-1.0, 0.5, 0.5)));
        clear_tree(&solid.0);
    }

    #[test]
    fn voxel_union_of_overlapping_boxes() {
        let b = TopoBuilder::new();
        let box1 = make_box_solid(&b, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let box2 = make_box_solid(&b, &GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(1.5, 1.5, 1.5));
        let mesh = voxel_union(&box1.0, &box2.0, 24).expect("voxel union ok");
        assert!(!mesh.triangles.is_empty(), "expected a fused surface");
        let vol = union_volume(&box1.0, &box2.0, 24);
        // True union volume is 1.875; voxel estimate over-counts (phantom
        // planes), but must at least contain one unit box and stay within the
        // combined bbox (1.5³ = 3.375).
        assert!(vol > 1.0, "volume {vol}");
        assert!(vol < 3.5, "volume {vol}");
        clear_tree(&box1.0);
        clear_tree(&box2.0);
    }
}
