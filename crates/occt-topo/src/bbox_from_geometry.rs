//! Bounding box of a shape computed from its geometry.
//!
//! The box is `BRepBndLib::Add(S, B, /*useTriangulation*/ false)`
//! (`BRepBndLib.cxx:81-215`): the exact port lives in [`crate::brep_bnd_lib`],
//! which walks the faces (`BndLib_AddSurface::Add` restricted to the face's
//! `BRepTools::UVBounds`, and `BndLib_Add3dCurve::Add` on the boundary edges of
//! a plane face), the edges outside faces and the vertices outside edges, using
//! the analytic / control-net extrema of `geom_bnd_lib_curve3d` and
//! `geom_bnd_lib_surface3d`.
//!
//! The three entry points below keep their signatures so the boolean, pave and
//! mesh call sites are unchanged; the previous body sampled 64 points per edge
//! and a 16x16 UV grid per face, a self-invented envelope with no OCCT
//! counterpart (registered as finding B3 in the drift inventory).

use occt_core::bnd::BndBox;

use crate::shape::TopoShape;

/// Test fixture access to the global geometry registry (the tests below set
/// vertex / edge / face geometry directly).
#[cfg(test)]
fn reg() -> &'static crate::tgeometry::GeometryRegistry {
    crate::tgeometry::GeometryRegistry::global()
}

/// Axis-aligned bounding box of `shape` derived from its geometry.
///
/// `BRepBndLib::Add(shape, box, false)`.
pub fn shape_bbox(shape: &TopoShape) -> BndBox {
    crate::brep_bnd_lib::shape_bnd_box(shape)
}

/// Bounding box enlarged by `margin` in every direction (`Bnd_Box::Enlarge`).
pub fn shape_bbox_with_margin(shape: &TopoShape, margin: f64) -> BndBox {
    let mut b = shape_bbox(shape);
    b.enlarge(margin);
    b
}

/// Quick reject: `true` if `shape`'s bounding box does not overlap `other`.
pub fn is_shape_outside(shape: &TopoShape, other: &BndBox) -> bool {
    shape_bbox(shape).is_out_box(other)
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt;
    use crate::shape::Vertex;
    use crate::tgeometry::VertexGeom;
    use crate::wireframe::tests::{segment_edge, square_face};

    #[test]
    fn vertex_bbox_contains_point() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::new(1.0, 2.0, 3.0), tolerance: 0.0 });
        let b = shape_bbox(&v.0);
        assert!(!b.is_out(&GpPnt::new(1.0, 2.0, 3.0)));
        assert!(b.is_out(&GpPnt::new(5.0, 2.0, 3.0)));
    }

    #[test]
    fn edge_bbox_from_curve() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let b = shape_bbox(&e.0);
        let (x0, x1, y0, y1, z0, z1) = b.get().unwrap();
        // `BRepBndLib::Add` -> `BndLib_Add3dCurve::Add(BC, BRep_Tool::Tolerance(E), B)`
        // (`BRepBndLib.cxx:146-149`) and the tolerance is raised to
        // `Precision::Confusion()` (`BRep_Tool.cxx:881-895`), so the box of the
        // segment is `[-t, 1+t] x [-t, t] x [-t, t]`.
        let t = occt_core::precision::CONFUSION;
        assert!((x0 + t).abs() < 1e-12 && (x1 - (1.0 + t)).abs() < 1e-12);
        assert!((y0 + t).abs() < 1e-12 && (y1 - t).abs() < 1e-12);
        assert!((z0 + t).abs() < 1e-12 && (z1 - t).abs() < 1e-12);
    }

    #[test]
    fn face_bbox_from_surface_and_wires() {
        let f = square_face();
        let b = shape_bbox(&f.0);
        let (x0, x1, y0, y1, z0, z1) = b.get().unwrap();
        // A plane face contributes its boundary edge curves plus `Bnd_Box::
        // Enlarge(BRep_Tool::Tolerance(F))` (`BRepBndLib.cxx:127-131`), both of
        // which are `Precision::Confusion()` here.
        let t = occt_core::precision::CONFUSION;
        assert!((x0 + t).abs() < 1e-12 && (x1 - (1.0 + t)).abs() < 1e-12);
        assert!((y0 + t).abs() < 1e-12 && (y1 - (1.0 + t)).abs() < 1e-12);
        assert!((z0 + t).abs() < 1e-12 && (z1 - t).abs() < 1e-12);
    }

    #[test]
    fn bbox_with_margin_enlarges() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::zero(), tolerance: 0.0 });
        let b = shape_bbox_with_margin(&v.0, 0.5);
        let (x0, x1, y0, _y1, z0, _z1) = b.get().unwrap();
        assert!((x0 + 0.5).abs() < 1e-12 && (x1 - 0.5).abs() < 1e-12);
        assert!((y0 + 0.5).abs() < 1e-12 && (z0 + 0.5).abs() < 1e-12);
    }

    #[test]
    fn is_shape_outside_quick_reject() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::zero(), tolerance: 0.0 });
        let far = BndBox::from_corners(&GpPnt::new(10.0, 10.0, 10.0), &GpPnt::new(11.0, 11.0, 11.0));
        assert!(is_shape_outside(&v.0, &far));
        let near = BndBox::from_corners(&GpPnt::new(-1.0, -1.0, -1.0), &GpPnt::new(1.0, 1.0, 1.0));
        assert!(!is_shape_outside(&v.0, &near));
    }
}
