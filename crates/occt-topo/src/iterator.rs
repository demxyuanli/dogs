//! Shape iterators — walk children of a shape. Source: `TopoDS_Iterator` /
//! `TopExp_Explorer`.
use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;

/// Direct children of `s` with `TopoDS_Iterator` defaults (`cumOri` and
/// `cumLoc` both true): Compose parent orientation, Move parent location.
///
/// Storage on `TShape` keeps the relative view written by `TopoDS_Builder::Add`.
/// Callers that copy children in order to re-`Add` them must read storage
/// directly; every *consumer* walk (MapShapes, TopExp, wire/edge collectors)
/// should go through this helper.
pub fn cumulated_children(s: &TopoShape) -> Vec<TopoShape> {
    let stored = s.tshape.read().expect("poisoned TShape lock").children.clone();
    let parent_ori = s.orientation();
    let parent_loc = s.location().clone();
    let cum_loc = !parent_loc.is_identity();
    stored
        .into_iter()
        .map(|mut child| {
            child.set_orientation(Orientation::compose(parent_ori, child.orientation()));
            if cum_loc {
                child.move_location(&parent_loc);
            }
            child
        })
        .collect()
}

/// A single shape in the flat child list.
#[derive(Debug, Clone)]
pub struct ChildEntry { pub shape: TopoShape }

/// Iterator over the direct children of a shape (`TopoDS_Iterator`).
#[derive(Debug)]
pub struct ShapeIterator {
    children: Vec<TopoShape>,
    pos: usize,
}

impl ShapeIterator {
    pub fn new(children: Vec<TopoShape>) -> Self { Self { children, pos: 0 } }

    /// Direct children of `s`, matching `TopoDS_Iterator` with `cumOri`/`cumLoc`
    /// both true.
    pub fn of_shape(s: &TopoShape) -> Self {
        Self::new(cumulated_children(s))
    }
}

impl Iterator for ShapeIterator {
    type Item = TopoShape;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.children.len() { return None; }
        let item = self.children[self.pos].clone();
        self.pos += 1;
        Some(item)
    }
}

/// Explorer-style traversal: iterate sub-shapes of a given type.
/// Non-matching nodes recurse through [`cumulated_children`].
#[derive(Debug)]
pub struct ShapeExplorer {
    stack: Vec<TopoShape>,
    target: ShapeType,
}

impl ShapeExplorer {
    pub fn new(children: Vec<TopoShape>, target: ShapeType) -> Self {
        Self { stack: children, target }
    }
}

impl Iterator for ShapeExplorer {
    type Item = TopoShape;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(s) = self.stack.pop() {
            if s.shape_type() == self.target {
                return Some(s);
            }
            for child in cumulated_children(&s).into_iter().rev() {
                self.stack.push(child);
            }
        }
        None
    }
}

/// Depth-first count of sub-shapes by type over a provided child tree.
pub fn count_shapes(children: &[TopoShape], target: ShapeType) -> usize {
    children.iter().filter(|s| s.shape_type() == target).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shape::Compound;

    #[test]
    fn iter_children() {
        let c = Compound::new();
        let _ = c.0;
        let children = vec![TopoShape::new(ShapeType::Vertex), TopoShape::new(ShapeType::Edge), TopoShape::new(ShapeType::Vertex)];
        let it = ShapeIterator::new(children);
        assert_eq!(it.count(), 3);
    }

    #[test]
    fn count_by_type() {
        let children = vec![TopoShape::new(ShapeType::Vertex), TopoShape::new(ShapeType::Face), TopoShape::new(ShapeType::Vertex)];
        assert_eq!(count_shapes(&children, ShapeType::Vertex), 2);
        assert_eq!(count_shapes(&children, ShapeType::Face), 1);
    }
}
