//! `TopOpeBRepBuild_ShapeListOfShape` — a shape plus a list of shapes.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ShapeListOfShape.hxx:24-54`,
//! `.cxx:33-77`.
//!
//! `myShape` is a `TopoDS_Shape` that the default constructor leaves null
//! (`cxx:33`), so the port stores [`Option<TopoShape>`].

use crate::shape::TopoShape;

/// `class TopOpeBRepBuild_ShapeListOfShape` (`hxx:27`).
#[derive(Debug, Clone, Default)]
pub struct ShapeListOfShape {
    /// `myList` (`hxx:50`).
    my_list: Vec<TopoShape>,
    /// `myShape` (`hxx:51`), null in the default constructor.
    my_shape: Option<TopoShape>,
}

impl ShapeListOfShape {
    /// `TopOpeBRepBuild_ShapeListOfShape()` (`cxx:33`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepBuild_ShapeListOfShape(S)` (`cxx:37-40`).
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            my_list: Vec::new(),
            my_shape: Some(s.clone()),
        }
    }

    /// `TopOpeBRepBuild_ShapeListOfShape(S, L)` (`cxx:44-49`).
    pub fn with_list(s: &TopoShape, l: Vec<TopoShape>) -> Self {
        Self {
            my_list: l,
            my_shape: Some(s.clone()),
        }
    }

    /// `List()` (`cxx:53-56`).
    pub fn list(&self) -> &[TopoShape] {
        &self.my_list
    }

    /// `ChangeList()` (`cxx:60-63`).
    pub fn change_list(&mut self) -> &mut Vec<TopoShape> {
        &mut self.my_list
    }

    /// `Shape()` (`cxx:67-70`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `ChangeShape()` (`cxx:74-77`).
    pub fn change_shape(&mut self) -> &mut Option<TopoShape> {
        &mut self.my_shape
    }
}
