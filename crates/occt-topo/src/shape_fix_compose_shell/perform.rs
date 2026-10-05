//! ShapeFix_ComposeShell::Perform / SplitEdges (ShapeFix_ComposeShell.cxx:206-270).

use crate::builder::TopoBuilder;
use crate::shape::Face;

use super::load_wires::load_wires;
use super::shell::{ComposeShell, SHAPEEXTEND_DONE1, SHAPEEXTEND_FAIL6, SHAPEEXTEND_OK};

impl ComposeShell {
    /// `ShapeFix_ComposeShell::Perform()` (`cxx:206-255`).
    pub fn perform(&mut self) -> bool {
        self.status = SHAPEEXTEND_OK;
        self.invert_edge_status = false; // cxx:209

        let face = match self.face() {
            Some(f) => f.clone(),
            None => {
                self.status = SHAPEEXTEND_FAIL6;
                return false;
            }
        };

        // cxx:214-219.
        let mut seqw = load_wires(&face, &self.context);
        if seqw.is_empty() {
            self.status = SHAPEEXTEND_FAIL6;
            return false;
        }

        // cxx:222.
        self.split_by_grid(&mut seqw);
        // cxx:225.
        self.break_wires(&mut seqw);
        // cxx:228-229.
        let mut wires: Vec<super::wire_segment::WireSegment> = Vec::new();
        self.collect_wires(&mut wires, &mut seqw);

        // cxx:232-233.
        let mut faces: Vec<crate::shape::TopoShape> = Vec::new();
        self.dispatch_wires(&mut faces, &mut wires);

        // cxx:236-251.
        let builder = TopoBuilder::new();
        let result = if faces.len() != 1 {
            let fs: Vec<Face> = faces.iter().map(|f| Face(f.clone())).collect();
            builder.make_shell(&fs).0
        } else {
            faces[0].clone()
        };
        let mut result = result;
        result.set_orientation(self.orient);
        self.result = Some(result);

        self.status |= SHAPEEXTEND_DONE1; // cxx:253
        true
    }

    /// `ShapeFix_ComposeShell::SplitEdges()` (`cxx:259-270`).
    pub fn split_edges(&mut self) {
        self.status = SHAPEEXTEND_OK;
        let face = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };
        let mut seqw = load_wires(&face, &self.context);
        self.split_by_grid(&mut seqw);
    }
}
