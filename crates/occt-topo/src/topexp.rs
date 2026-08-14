//! Topology exploration — a port of `TopExp_Explorer`.
//!
//! `TShape` does not store children in this port, so the explorer cannot walk
//! the topology tree on its own. The caller seeds children for the root with
//! [`Explorer::set_children`], or supplies a `children_of` callback with
//! [`Explorer::with_children_fn`]; the explorer then walks the tree depth-first,
//! yielding every sub-shape whose type matches the target.

use std::sync::Arc;

use crate::abs::ShapeType;
use crate::shape::TopoShape;

/// A named association between two shapes (TopExp_StackEntry style).
#[derive(Debug, Clone)]
pub struct ShapeAssociation {
    pub parent: TopoShape,
    pub child: TopoShape,
}

/// Depth-first topology explorer.
///
/// `target == ShapeType::Shape` selects *all* sub-shapes; the root itself is
/// not reported, matching TopExp's notion of "all sub-shapes".
pub struct Explorer {
    stack: Vec<TopoShape>,
    target: ShapeType,
    cur: Option<TopoShape>,
    root: Option<TopoShape>,
    children_fn: Option<Box<dyn Fn(&TopoShape) -> Vec<TopoShape>>>,
}

impl Explorer {
    /// Start exploring `s` for sub-shapes of type `target`.
    pub fn new(s: &TopoShape, target: ShapeType) -> Self {
        let mut e = Explorer {
            stack: vec![s.clone()],
            target,
            cur: None,
            root: Some(s.clone()),
            children_fn: None,
        };
        e.advance();
        e
    }

    /// Whether there is another matching sub-shape to yield.
    pub fn more(&self) -> bool {
        self.cur.is_some()
    }

    /// The current sub-shape. Panics once the explorer is exhausted.
    pub fn current(&self) -> &TopoShape {
        self.cur
            .as_ref()
            .expect("Explorer::current called after exhaustion")
    }

    /// Advance to the next matching sub-shape.
    pub fn next(&mut self) {
        self.advance();
    }

    /// Seed the root's children, then restart the walk over them.
    pub fn set_children(&mut self, children: Vec<TopoShape>) {
        let root = self.root.clone();
        self.children_fn = Some(Box::new(move |s: &TopoShape| match &root {
            Some(r) if Arc::ptr_eq(&r.tshape, &s.tshape) => children.clone(),
            _ => Vec::new(),
        }));
        self.restart();
    }

    /// Install a callback that returns the children of a shape during traversal.
    pub fn with_children_fn<F>(&mut self, f: F)
    where
        F: Fn(&TopoShape) -> Vec<TopoShape> + 'static,
    {
        self.children_fn = Some(Box::new(f));
        self.restart();
    }

    /// Rewind to the root and re-run the walk.
    fn restart(&mut self) {
        if let Some(r) = self.root.clone() {
            self.stack = vec![r];
        }
        self.cur = None;
        self.advance();
    }

    fn children_of(&self, s: &TopoShape) -> Vec<TopoShape> {
        match &self.children_fn {
            Some(f) => f(s),
            // Real topology: walk the children stored on the TShape.
            None => s.tshape.read().unwrap().children.clone(),
        }
    }

    fn advance(&mut self) {
        self.cur = None;
        while let Some(s) = self.stack.pop() {
            let is_root = self
                .root
                .as_ref()
                .map_or(false, |r| Arc::ptr_eq(&r.tshape, &s.tshape));
            let matches = if self.target == ShapeType::Shape {
                !is_root
            } else {
                shape_type_of(&s) == self.target
            };
            let kids = self.children_of(&s);
            if matches {
                self.cur = Some(s);
                for k in kids.into_iter().rev() {
                    self.stack.push(k);
                }
                return;
            }
            for k in kids.into_iter().rev() {
                self.stack.push(k);
            }
        }
    }
}

/// Read the topological type of a shape (behind the shared `TShape` lock).
fn shape_type_of(s: &TopoShape) -> ShapeType {
    s.tshape.read().expect("poisoned TShape lock").shape_type()
}

/// Filter `shapes` to those whose type equals `target`.
pub fn map_shapes(shapes: &[TopoShape], target: ShapeType) -> Vec<TopoShape> {
    shapes
        .iter()
        .filter(|s| shape_type_of(s) == target)
        .cloned()
        .collect()
}

/// The first shape of type `target`, if any.
pub fn first_shape(shapes: &[TopoShape], target: ShapeType) -> Option<TopoShape> {
    shapes
        .iter()
        .find(|s| shape_type_of(s) == target)
        .cloned()
}

/// Count the shapes of type `target`.
pub fn nb_shapes(shapes: &[TopoShape], target: ShapeType) -> usize {
    shapes
        .iter()
        .filter(|s| shape_type_of(s) == target)
        .count()
}

/// Collect the distinct vertices referenced by `edges`.
///
/// Reads the real vertex children stored on each edge when present; falls back
/// to uniquifying the input list by `TShape` identity for edges without
/// explicit vertex children.
pub fn vertices_from_edges(edges: &[TopoShape]) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for e in edges {
        let has_children = !e.tshape.read().unwrap().children.is_empty();
        let mut added = false;
        if has_children {
            let kids = e.tshape.read().unwrap().children.clone();
            for v in kids {
                if v.shape_type() != ShapeType::Vertex {
                    continue;
                }
                if !out.iter().any(|o| Arc::ptr_eq(&o.tshape, &v.tshape)) {
                    out.push(v);
                    added = true;
                }
            }
        }
        if !added && !out.iter().any(|o| Arc::ptr_eq(&o.tshape, &e.tshape)) {
            out.push(e.clone());
        }
    }
    out
}

/// Whether `ancestor` is an ancestor of `shape`.
///
/// Parent links are not walkable in this port; the only decidable case is
/// when both arguments name the same `TShape`.
pub fn is_shape_ancestor(ancestor: &TopoShape, shape: &TopoShape) -> bool {
    Arc::ptr_eq(&ancestor.tshape, &shape.tshape)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_shapes_filters_by_type() {
        let shapes = vec![
            TopoShape::new(ShapeType::Vertex),
            TopoShape::new(ShapeType::Edge),
            TopoShape::new(ShapeType::Vertex),
        ];
        let verts = map_shapes(&shapes, ShapeType::Vertex);
        assert_eq!(verts.len(), 2);
        assert_eq!(map_shapes(&shapes, ShapeType::Edge).len(), 1);
        assert!(map_shapes(&shapes, ShapeType::Face).is_empty());
    }

    #[test]
    fn nb_shapes_counts_and_first_shape_finds() {
        let shapes = vec![
            TopoShape::new(ShapeType::Vertex),
            TopoShape::new(ShapeType::Edge),
            TopoShape::new(ShapeType::Vertex),
        ];
        assert_eq!(nb_shapes(&shapes, ShapeType::Vertex), 2);
        assert_eq!(nb_shapes(&shapes, ShapeType::Edge), 1);
        assert_eq!(nb_shapes(&shapes, ShapeType::Face), 0);

        let e = first_shape(&shapes, ShapeType::Edge).expect("edge present");
        assert_eq!(shape_type_of(&e), ShapeType::Edge);
        assert!(first_shape(&shapes, ShapeType::Face).is_none());
    }

    #[test]
    fn explorer_yields_target_children_in_order() {
        let root = TopoShape::new(ShapeType::Compound);
        let v1 = TopoShape::new(ShapeType::Vertex);
        let e1 = TopoShape::new(ShapeType::Edge);
        let v2 = TopoShape::new(ShapeType::Vertex);

        let mut ex = Explorer::new(&root, ShapeType::Vertex);
        ex.set_children(vec![v1.clone(), e1.clone(), v2.clone()]);

        let mut found = Vec::new();
        while ex.more() {
            found.push(ex.current().clone());
            ex.next();
        }

        assert_eq!(found.len(), 2);
        assert!(Arc::ptr_eq(&found[0].tshape, &v1.tshape));
        assert!(Arc::ptr_eq(&found[1].tshape, &v2.tshape));
    }

    #[test]
    fn vertices_from_edges_uniquifies() {
        let a = TopoShape::new(ShapeType::Vertex);
        let b = TopoShape::new(ShapeType::Vertex);
        let edges = vec![a.clone(), b.clone(), a.clone()];
        let vs = vertices_from_edges(&edges);
        assert_eq!(vs.len(), 2);
    }

    #[test]
    fn ancestor_is_identity_for_same_tshape() {
        let s = TopoShape::new(ShapeType::Vertex);
        assert!(is_shape_ancestor(&s, &s));
        let other = TopoShape::new(ShapeType::Vertex);
        assert!(!is_shape_ancestor(&s, &other));
    }
}
