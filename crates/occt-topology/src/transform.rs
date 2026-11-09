//! Topological shape transformation — apply gp_Trsf to shapes.
//! Source: `TopoDS_Shape::Move`, BRep_Tool transforms.
use crate::shape::TopoShape;
use crate::abs::ShapeType;
use occt_core::gp::{GpTrsf, GpVec, GpPnt, GpAx1};
use occt_core::toploc::TopLocLocation;

/// Apply a transform to a shape by composing its location.
/// Mirrors TopoDS_Shape::Move(T) — the TShape is shared, only Location changes.
pub fn move_shape(shape: &TopoShape, t: &GpTrsf) -> TopoShape {
    // Compose the transform on top of the existing location.
    // new_location = T * existing_location (since Move applies T first).
    let mut result = shape.clone();
    result.location = TopLocLocation::composed(t, result.location.clone());
    result
}

/// Copy a shape and apply transform (non-mutating, shares TShape).
pub fn transformed(shape: &TopoShape, t: &GpTrsf) -> TopoShape { move_shape(shape, t) }

/// Translate a shape.
pub fn translated(shape: &TopoShape, v: &GpVec) -> TopoShape {
    let mut t = GpTrsf::identity();
    t.set_translation_vec(v);
    move_shape(shape, &t)
}

/// Rotate a shape around an axis.
pub fn rotated(shape: &TopoShape, axis: &GpAx1, angle: f64) -> TopoShape {
    let mut t = GpTrsf::identity();
    let _ = t.set_rotation_ax1(axis, angle);
    move_shape(shape, &t)
}

/// Scale a shape about a point.
pub fn scaled(shape: &TopoShape, p: &GpPnt, s: f64) -> TopoShape {
    let mut t = GpTrsf::identity();
    let _ = t.set_scale(p, s);
    move_shape(shape, &t)
}

/// Mirror a shape through a point.
pub fn mirrored(shape: &TopoShape, p: &GpPnt) -> TopoShape {
    let mut t = GpTrsf::identity();
    t.set_mirror_pnt(p);
    move_shape(shape, &t)
}

/// Apply transform in place (mutates the TopoShape view).
pub fn move_in_place(shape: &mut TopoShape, t: &GpTrsf) {
    shape.location = TopLocLocation::composed(t, shape.location.clone());
}

/// Compute the accumulated location transform of a shape chain.
pub fn cumulative_transform(shape: &TopoShape) -> GpTrsf {
    shape.location.transformation()
}

/// Transform a vertex's coordinates (if point data were available).
/// Placeholder: returns the location transform for composition.
pub fn transform_data(shape: &TopoShape, point: &GpPnt) -> GpPnt {
    point.transformed(&cumulative_transform(shape))
}

/// Apply the same transform to a list of shapes.
pub fn move_many(shapes: &[TopoShape], t: &GpTrsf) -> Vec<TopoShape> {
    shapes.iter().map(|s| move_shape(s, t)).collect()
}

/// Check if a shape's location is identity.
pub fn is_located(shape: &TopoShape) -> bool { !shape.location.is_identity() }

/// The type after transform is preserved (shape_type is on shared TShape).
pub fn assert_type_preserved(shape: &TopoShape, expected: ShapeType) -> bool {
    shape.shape_type() == expected
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpVec;

    #[test]
    fn translate_vertex() {
        let v = TopoShape::new(ShapeType::Vertex);
        let v2 = translated(&v, &GpVec::new(1., 2., 3.));
        // Location should no longer be identity
        assert!(!v2.location.is_identity());
        // Original untouched (shared TShape, different location)
        assert!(v.location.is_identity());
    }

    #[test]
    fn type_preserved() {
        let f = TopoShape::new(ShapeType::Face);
        let f2 = translated(&f, &GpVec::new(1., 0., 0.));
        assert!(assert_type_preserved(&f2, ShapeType::Face));
    }

    #[test]
    fn cumulative_identity() {
        let v = TopoShape::new(ShapeType::Vertex);
        let t = cumulative_transform(&v);
        assert!(t.form() == occt_core::gp::TrsfForm::Identity);
    }
}
