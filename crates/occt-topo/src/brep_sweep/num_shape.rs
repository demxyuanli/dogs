//! `Sweep_NumShape`, `Sweep_NumShapeIterator`, `Sweep_NumShapeTool`.
//!
//! Directing-edge indexation for `BRepSweep_NumLinearRegularSweep`.
//! Sources: `Sweep_NumShape.cxx`, `Sweep_NumShapeIterator.cxx`,
//! `Sweep_NumShapeTool.cxx`.

use crate::abs::{Orientation, ShapeType};

/// Indexed directing edge or one of its vertices (`Sweep_NumShape`).
#[derive(Debug, Clone, Copy)]
pub struct NumShape {
    pub ty: ShapeType,
    pub index: i32,
    pub closed: bool,
    pub beg_inf: bool,
    pub end_inf: bool,
}

impl NumShape {
    pub fn dummy() -> Self {
        Self {
            ty: ShapeType::Shape,
            index: 0,
            closed: false,
            beg_inf: false,
            end_inf: false,
        }
    }

    /// `Sweep_NumShape::Init` (`Sweep_NumShape.cxx:48-59`).
    pub fn init(&mut self, index: i32, ty: ShapeType, closed: bool, beg_inf: bool, end_inf: bool) {
        self.index = index;
        self.ty = ty;
        self.closed = closed;
        self.beg_inf = beg_inf;
        self.end_inf = end_inf;
    }

    /// `Sweep_NumShape::Orientation` (`Sweep_NumShape.cxx:63-79`).
    pub fn orientation(self) -> Orientation {
        if self.ty == ShapeType::Edge {
            Orientation::Forward
        } else if self.index == 2 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }
}

/// `Sweep_NumShapeTool` (`Sweep_NumShapeTool.cxx`).
#[derive(Debug, Clone, Copy)]
pub struct NumShapeTool {
    shape: NumShape,
}

impl NumShapeTool {
    pub fn new(shape: NumShape) -> Self {
        Self { shape }
    }

    pub fn nb_shapes(self) -> i32 {
        if self.shape.ty == ShapeType::Edge {
            if self.shape.closed {
                self.shape.index
            } else {
                self.shape.index + 1
            }
        } else {
            1
        }
    }

    pub fn index(self, shape: NumShape) -> i32 {
        if shape.ty == ShapeType::Edge {
            1
        } else if shape.closed {
            2
        } else {
            shape.index + 1
        }
    }

    pub fn shape(self, an_index: i32) -> NumShape {
        if an_index == 1 {
            self.shape
        } else {
            NumShape {
                ty: ShapeType::Vertex,
                index: an_index - 1,
                closed: self.shape.closed,
                beg_inf: false,
                end_inf: false,
            }
        }
    }

    pub fn has_first_vertex(self) -> bool {
        if self.shape.ty == ShapeType::Edge {
            !self.shape.beg_inf
        } else {
            true
        }
    }

    pub fn has_last_vertex(self) -> bool {
        if self.shape.ty == ShapeType::Edge {
            !self.shape.end_inf
        } else {
            true
        }
    }

    pub fn first_vertex(self) -> NumShape {
        if self.shape.ty == ShapeType::Edge {
            NumShape {
                ty: ShapeType::Vertex,
                index: 1,
                closed: self.shape.closed,
                beg_inf: false,
                end_inf: false,
            }
        } else {
            self.shape
        }
    }

    pub fn last_vertex(self) -> NumShape {
        if self.shape.ty == ShapeType::Edge {
            NumShape {
                ty: ShapeType::Vertex,
                index: self.nb_shapes() - 1,
                closed: self.shape.closed,
                beg_inf: false,
                end_inf: false,
            }
        } else {
            self.shape
        }
    }
}

/// `Sweep_NumShapeIterator` (`Sweep_NumShapeIterator.cxx`).
#[derive(Debug, Clone)]
pub struct NumShapeIterator {
    shape: NumShape,
    current: NumShape,
    range: i32,
    more: bool,
    orientation: Orientation,
}

impl NumShapeIterator {
    pub fn new(shape: NumShape) -> Self {
        let mut it = Self {
            shape,
            current: NumShape::dummy(),
            range: 0,
            more: false,
            orientation: Orientation::Forward,
        };
        it.init(shape);
        it
    }

    fn init(&mut self, shape: NumShape) {
        self.shape = shape;
        if self.shape.ty == ShapeType::Edge {
            let nbvert = self.shape.index;
            self.more = nbvert >= 1;
            if self.more {
                self.range = 1;
                self.current = NumShape {
                    ty: ShapeType::Vertex,
                    index: 1,
                    closed: self.shape.closed,
                    beg_inf: false,
                    end_inf: false,
                };
                if nbvert == 1 {
                    self.orientation = if self.shape.beg_inf {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    };
                } else {
                    self.orientation = Orientation::Forward;
                }
            }
        }
    }

    pub fn more(&self) -> bool {
        self.more
    }

    pub fn value(&self) -> NumShape {
        self.current
    }

    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    pub fn next(&mut self) {
        self.range += 1;
        self.more = self.range <= self.shape.index;
        if self.more && self.shape.ty == ShapeType::Edge {
            self.current = NumShape {
                ty: ShapeType::Vertex,
                index: self.range,
                closed: self.shape.closed,
                beg_inf: false,
                end_inf: false,
            };
            self.orientation = Orientation::Reversed;
        }
    }
}
