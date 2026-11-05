//! `BRepFill_FaceAndOrder` and `BRepFill_EdgeFaceAndOrder` — the two record
//! types `BRepFill_Filling` keeps its constraints in.
//!
//! Source: `TKBool/BRepFill/BRepFill_FaceAndOrder.{hxx,cxx}` and
//! `BRepFill_EdgeFaceAndOrder.{hxx,cxx}` (OCCT 8.0.0). Both are plain
//! "face/edge + continuity" records (`FaceAndOrder.cxx:22-32`,
//! `EdgeFaceAndOrder.cxx:22-35`); the defaulted default constructor leaves
//! `myOrder` uninitialised in the reference, and no `BRepFill_Filling` call site
//! uses it, so the port settles it on `Shape::C0`.

use occt_core::kernel::geomabs::Shape;

use crate::shape::{Edge, Face};

/// `class BRepFill_FaceAndOrder` (`hxx:28-43`): a face plus the continuity the
/// constraint asks for.
#[derive(Debug, Clone)]
pub struct FaceAndOrder {
    /// `myFace` (`hxx:41`).
    pub face: Face,
    /// `myOrder` (`hxx:42`).
    pub order: Shape,
}

impl Default for FaceAndOrder {
    fn default() -> Self {
        Self::new()
    }
}

impl FaceAndOrder {
    /// `BRepFill_FaceAndOrder()` (`cxx:22`).
    pub fn new() -> Self {
        Self {
            face: Face::new(),
            order: Shape::C0,
        }
    }

    /// `BRepFill_FaceAndOrder(aFace, anOrder)` (`cxx:24-28`).
    pub fn with_face(face: &Face, order: Shape) -> Self {
        Self {
            face: face.clone(),
            order,
        }
    }
}

/// `class BRepFill_EdgeFaceAndOrder` (`hxx:30-44`): an edge, the face carrying
/// its tangency/curvature reference (null for a pure G0 bound) and the
/// continuity.
#[derive(Debug, Clone)]
pub struct EdgeFaceAndOrder {
    /// `myEdge` (`hxx:41`).
    pub edge: Edge,
    /// `myFace` (`hxx:42`): `None` is the reference's null `TopoDS_Face`, the
    /// `Add(anEdge, Order, IsBound)` overload (`Filling.cxx:217-233`).
    pub face: Option<Face>,
    /// `myOrder` (`hxx:43`).
    pub order: Shape,
}

impl Default for EdgeFaceAndOrder {
    fn default() -> Self {
        Self::new()
    }
}

impl EdgeFaceAndOrder {
    /// `BRepFill_EdgeFaceAndOrder()` (`cxx:22`).
    pub fn new() -> Self {
        Self {
            edge: Edge::new(),
            face: None,
            order: Shape::C0,
        }
    }

    /// `BRepFill_EdgeFaceAndOrder(anEdge, aFace, anOrder)` (`cxx:24-31`);
    /// `face` stays `None` for the null `TopoDS_Face` of the support-less
    /// overload.
    pub fn with_face(edge: &Edge, face: Option<&Face>, order: Shape) -> Self {
        Self {
            edge: edge.clone(),
            face: face.cloned(),
            order,
        }
    }
}

/// `GeomAbs_Shape` to the `int` the `GeomPlate_CurveConstraint` constructors
/// take (`GeomAbs_Shape.hxx:24-33`: `C0 = 0`, `G1 = 1`, `G2 = 2`, `C1 = 3`,
/// `C2 = 4`, `C3 = 5`, `CN = 6`).
///
/// The Rust enum orders `C1` before `G2`, so a plain `as i32` would turn
/// `GeomAbs_C1` into plate order 2 and silently accept it. OCCT instead pushes
/// 3..6 through and lets the constructor reject anything above `G2`
/// (`GeomPlate_CurveConstraint.cxx:71-75`).
pub fn shape_to_plate_order(order: Shape) -> i32 {
    match order {
        Shape::C0 => 0,
        Shape::G1 => 1,
        Shape::G2 => 2,
        Shape::C1 => 3,
        Shape::C2 => 4,
        Shape::C3 => 5,
        Shape::CN => 6,
    }
}
