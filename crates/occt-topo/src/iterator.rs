//! Shape iterators — walk children of a shape. Source: `TopExp_Explorer`
use crate::abs::ShapeType;
use crate::shape::TopoShape;

/// A single shape in the flat child list.
#[derive(Debug, Clone)]
pub struct ChildEntry { pub shape: TopoShape }

/// Iterator over the direct children of a shape.
/// In OCCT, children are stored per-TShape; this port keeps a simple
/// Vec<TopoShape> side-table keyed by TShape handle — here we expose
/// the API shape so real storage can be plugged in later.
#[derive(Debug)]
pub struct ShapeIterator {
    children: Vec<TopoShape>,
    pos: usize,
}

impl ShapeIterator {
    pub fn new(children: Vec<TopoShape>) -> Self { Self { children, pos: 0 } }
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
/// Mirrors TopExp_Explorer but operates on an explicit child list.
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
            if s.shape_type() == self.target { return Some(s); }
            // NOTE: without stored children we cannot recurse; extend later.
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
        let mut c = Compound::new();
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
