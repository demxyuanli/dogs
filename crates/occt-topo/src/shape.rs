//! TopoDS shape wrappers — typed views over shared TShape data.
//! Source: `TopoDS_Shape`, `TopoDS_Vertex`, etc.
use std::sync::{Arc, RwLock};
use crate::abs::{ShapeType, Orientation};
use crate::tshape::{HandleTShape, TShape};
use occt_core::toploc::TopLocLocation;

/// Topological shape — shared data + location + orientation.
/// Mirrors OCCT's TopoDS_Shape (handle to TShape + Location + Orientation).
#[derive(Debug, Clone)]
pub struct TopoShape {
    pub tshape: HandleTShape,
    pub location: TopLocLocation,
    pub orientation: Orientation,
}

impl TopoShape {
    /// Wrap an existing TShape handle.
    pub fn from_handle(tshape: HandleTShape) -> Self {
        Self { tshape, location: TopLocLocation::identity(), orientation: Orientation::Forward }
    }

    /// Create a new empty shape of given type.
    pub fn new(shape_type: ShapeType) -> Self {
        Self::from_handle(Arc::new(RwLock::new(TShape::new(shape_type))))
    }

    pub fn shape_type(&self) -> ShapeType {
        self.tshape.read().unwrap().shape_type
    }

    pub fn is_null(&self) -> bool {
        // OCCT considers a shape null if its TShape is null.
        // Here all shapes are non-null; use Option<TopoShape> for null.
        false
    }

    pub fn is_vertex(&self) -> bool { self.shape_type() == ShapeType::Vertex }
    pub fn is_edge(&self) -> bool { self.shape_type() == ShapeType::Edge }
    pub fn is_wire(&self) -> bool { self.shape_type() == ShapeType::Wire }
    pub fn is_face(&self) -> bool { self.shape_type() == ShapeType::Face }
    pub fn is_shell(&self) -> bool { self.shape_type() == ShapeType::Shell }
    pub fn is_solid(&self) -> bool { self.shape_type() == ShapeType::Solid }
    pub fn is_compound(&self) -> bool { self.shape_type() == ShapeType::Compound }

    pub fn free(&self) -> bool { self.tshape.read().unwrap().free() }
    pub fn set_free(&self, v: bool) { self.tshape.write().unwrap().set_free(v); }
    pub fn closed(&self) -> bool { self.tshape.read().unwrap().closed() }
    pub fn set_closed(&self, v: bool) { self.tshape.write().unwrap().set_closed(v); }

    pub fn location(&self) -> &TopLocLocation { &self.location }
    pub fn set_location(&mut self, l: &TopLocLocation) { self.location = l.clone(); }
    pub fn move_location(&mut self, l: &TopLocLocation) {
        self.location = TopLocLocation::composed(&l.transformation(), self.location.clone());
    }

    pub fn orientation(&self) -> Orientation { self.orientation }
    pub fn set_orientation(&mut self, o: Orientation) { self.orientation = o; }

    /// `TopoDS_Shape::Reverse` — swap Forward/Reversed; Internal/External unchanged.
    pub fn reverse(&mut self) {
        self.orientation = self.orientation.reversed();
    }

    /// Copy with reversed orientation.
    pub fn oriented(&self, o: Orientation) -> Self {
        Self { tshape: self.tshape.clone(), location: self.location.clone(), orientation: o }
    }

    /// Copy with same location, no transform.
    pub fn located(&self, l: &TopLocLocation) -> Self {
        Self { tshape: self.tshape.clone(), location: l.clone(), orientation: self.orientation }
    }

    /// Deep-ish copy (shares TShape data, copies view).
    pub fn copied(&self) -> Self { self.clone() }

    /// True if shapes share the same TShape data.
    pub fn same_tshape(&self, other: &Self) -> bool { Arc::ptr_eq(&self.tshape, &other.tshape) }
}

// Typed shape wrappers — thin newtypes delegating to TopoShape.
macro_rules! typed_shape {
    ($name:ident, $ty:ident, $check:ident) => {
        #[derive(Debug, Clone)]
        pub struct $name(pub TopoShape);
        impl $name {
            pub fn new() -> Self { Self(TopoShape::new(ShapeType::$ty)) }
            pub fn wrap(s: TopoShape) -> Option<Self> { if s.$check() { Some(Self(s)) } else { None } }
            pub fn shape(&self) -> &TopoShape { &self.0 }
        }
        impl From<$name> for TopoShape { fn from(v: $name) -> Self { v.0 } }
        impl std::ops::Deref for $name { type Target = TopoShape; fn deref(&self) -> &TopoShape { &self.0 } }
    };
}

typed_shape!(Vertex, Vertex, is_vertex);
typed_shape!(Edge, Edge, is_edge);
typed_shape!(Wire, Wire, is_wire);
typed_shape!(Face, Face, is_face);
typed_shape!(Shell, Shell, is_shell);
typed_shape!(Solid, Solid, is_solid);
typed_shape!(Compound, Compound, is_compound);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_types() {
        let s = TopoShape::new(ShapeType::Vertex);
        assert!(s.is_vertex());
        assert!(!s.is_face());
        assert_eq!(s.shape_type(), ShapeType::Vertex);
    }

    #[test]
    fn orientation_flip() {
        let mut s = TopoShape::new(ShapeType::Face);
        assert!(s.orientation().is_forward());
        s.set_orientation(Orientation::Reversed);
        assert!(s.orientation().is_reversed());
    }

    #[test]
    fn typed_wrappers() {
        let v = Vertex::new();
        assert!(v.is_vertex());
        let s: TopoShape = v.into();
        assert!(s.is_vertex());
    }
}
