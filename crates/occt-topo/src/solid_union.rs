//! Mesh soup merge and bounding-box helpers.
//!
//! **UNPORTED (audit A5 / A9 / A14)**: OCCT has no voxel boolean and no triangle-
//! soup "fuse". Real solid boolean is `BRepAlgoAPI_*` -> `BOPAlgo_BOP`
//! (`bop_builder::boolean` / `boolean_via_bopalgo`). The previous file header
//! claiming `BRepAlgoAPI_Fuse` / `BRepAlgoAPI_BooleanOperation` was a false
//! provenance (A9); the voxel ray-cast path is deleted with no production
//! caller left after T-41 / T-79 removed the dispatch fallbacks.

use occt_core::bnd::BndBox;
use occt_core::poly::triangulation::Triangle;

use crate::abs::ShapeType;
use crate::mesh::ShapeMesh;

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

/// Merge two bounding boxes (`Bnd_Box::Add`).
pub fn bbox_union(a: &BndBox, b: &BndBox) -> BndBox {
    let mut out = *a;
    out.add_box(b);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{mesh_box, shape_mesh_surface_area};
    use occt_core::gp::GpPnt;

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
        assert!((shape_mesh_surface_area(&u) - 12.0).abs() < 1e-9);
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
}
