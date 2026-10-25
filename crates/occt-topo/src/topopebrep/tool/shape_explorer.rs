//! `TopOpeBRepTool_ShapeExplorer` — `TopExp_Explorer` plus a running index.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_ShapeExplorer.hxx:24-77`.
//!
//! Translation notes:
//! - The class is header-only: it forwards `Init` / `More` / `Current` / `Next`
//!   to `TopExp_Explorer` and only adds `myIndex` (`hxx:85`), which starts at 1
//!   on a non-empty exploration and advances only while `More()` holds
//!   (`hxx:61-66`).
//! - The port materialises the exploration once through
//!   [`crate::topo_tools::explorer_shapes_avoiding`], which is the port's
//!   `TopExp_Explorer(shape, toFind, toAvoid)` (`TopExp_Explorer.cxx:110-172`).

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes_avoiding;

/// `class TopOpeBRepTool_ShapeExplorer` (`hxx:27-87`).
pub struct ShapeExplorer {
    /// The explored sub-shapes, in `TopExp_Explorer` order.
    shapes: Vec<TopoShape>,
    /// The cursor (`TopExp_Explorer::myCurrent`).
    pos: usize,
    /// `myIndex` (`hxx:85`). 1-based, 0 when the exploration is empty.
    index: i32,
}

impl ShapeExplorer {
    /// `TopOpeBRepTool_ShapeExplorer()` (`hxx:30-33`).
    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            pos: 0,
            index: 0,
        }
    }

    /// `TopOpeBRepTool_ShapeExplorer(S, ToFind, ToAvoid)` (`hxx:43-51`).
    pub fn with_shape(s: &TopoShape, to_find: ShapeType, to_avoid: ShapeType) -> Self {
        let mut e = Self::new();
        e.init(s, to_find, to_avoid);
        e
    }

    /// `Init(S, ToFind, ToAvoid)` (`hxx:53-58`).
    pub fn init(&mut self, s: &TopoShape, to_find: ShapeType, to_avoid: ShapeType) {
        self.shapes = explorer_shapes_avoiding(s, to_find, to_avoid);
        self.pos = 0;
        self.index = if self.more() { 1 } else { 0 };
    }

    /// `More()` (`TopExp_Explorer.hxx`, forwarded).
    pub fn more(&self) -> bool {
        self.pos < self.shapes.len()
    }

    /// `Current()` (`TopExp_Explorer.hxx`, forwarded).
    pub fn current(&self) -> Option<&TopoShape> {
        self.shapes.get(self.pos)
    }

    /// `Next()` (`hxx:61-66`): `myIndex++` only while `More()`, then advance.
    pub fn next(&mut self) {
        if self.more() {
            self.index += 1;
        }
        self.pos += 1;
    }

    /// `Index()` (`hxx:69`).
    pub fn index(&self) -> i32 {
        self.index
    }
}

impl Default for ShapeExplorer {
    fn default() -> Self {
        Self::new()
    }
}
