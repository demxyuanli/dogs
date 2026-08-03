//! Connexity block — a connected set of shapes together with the loops
//! (closed shells or wires) extracted from them.
//! Source: `BOPTools_ConnexityBlock.hxx`

use crate::shape::TopoShape;

/// A connexity block is a connected component of shapes of the same kind
/// (faces for the shell splitter, edges for the wire splitter) plus the
/// closed loops built from those shapes.
///
/// The `regular` flag mirrors `BOPTools_ConnexityBlock::IsRegular()`:
/// a block is *regular* when it forms a single closed loop with no
/// multi-connected elements (e.g. every edge of a face block is shared by
/// exactly two faces), so it can be turned into a loop directly without the
/// generic graph traversal.
#[derive(Debug, Clone)]
pub struct ConnexityBlock {
    /// Shapes forming the block (faces / edges).
    shapes: Vec<TopoShape>,
    /// Closed loops (shells / wires) split from `shapes`.
    loops: Vec<TopoShape>,
    /// Whether the block is regular (single simple loop).
    is_regular: bool,
}

impl Default for ConnexityBlock {
    fn default() -> Self {
        Self::new()
    }
}

impl ConnexityBlock {
    /// Empty block, regular by default (matches the OCCT default state).
    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            loops: Vec::new(),
            is_regular: true,
        }
    }

    /// The shapes forming the block.
    pub fn shapes(&self) -> &[TopoShape] {
        &self.shapes
    }

    /// Mutable access to the block shapes (e.g. to re-chain a wire block).
    pub fn change_shapes_mut(&mut self) -> &mut Vec<TopoShape> {
        &mut self.shapes
    }

    /// Set the regularity flag.
    pub fn set_regular(&mut self, flag: bool) {
        self.is_regular = flag;
    }

    /// Whether the block is regular.
    pub fn is_regular(&self) -> bool {
        self.is_regular
    }

    /// The closed loops (shells / wires) built from the block shapes.
    pub fn loops(&self) -> &[TopoShape] {
        &self.loops
    }

    /// Mutable access to the loops.
    pub fn change_loops_mut(&mut self) -> &mut Vec<TopoShape> {
        &mut self.loops
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::ShapeType;

    #[test]
    fn default_state_is_regular_and_empty() {
        let cb = ConnexityBlock::new();
        assert!(cb.is_regular());
        assert!(cb.shapes().is_empty());
        assert!(cb.loops().is_empty());
    }

    #[test]
    fn shapes_loops_roundtrip() {
        let mut cb = ConnexityBlock::new();
        let s = TopoShape::new(ShapeType::Face);
        cb.change_shapes_mut().push(s.clone());
        assert_eq!(cb.shapes().len(), 1);
        assert!(crate::topo_tools_full::is_same(&cb.shapes()[0], &s));

        cb.set_regular(false);
        assert!(!cb.is_regular());

        let w = TopoShape::new(ShapeType::Wire);
        cb.change_loops_mut().push(w.clone());
        assert_eq!(cb.loops().len(), 1);
        assert!(crate::topo_tools_full::is_same(&cb.loops()[0], &w));
    }

    #[test]
    fn default_impl_matches_new() {
        assert_eq!(ConnexityBlock::default().is_regular(), ConnexityBlock::new().is_regular());
    }
}
