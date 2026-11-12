//! `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:345-736`).

use std::sync::Arc;

use occt_core::precision::CONFUSION;

use super::prelude::*;
use super::ShapeFixFace;
use crate::shape_fix_compose_shell::ReShape;
use crate::shhealing::{
    least_edge_size, perform_wire_first_round, perform_wire_second_round,
    split_face_after_removed_segment, split_loop_wires,
};

const DONE1: i32 = 0x0001;
const DONE2: i32 = 0x0002;
const DONE3: i32 = 0x0004;
const DONE4: i32 = 0x0008;
const DONE5: i32 = 0x0010;
const DONE6: i32 = 0x0020;
const DONE7: i32 = 0x0040;
const DONE8: i32 = 0x0080;

fn shape_key(shape: &TopoShape) -> usize {
    Arc::as_ptr(&shape.tshape) as usize
}

fn wire_state(wire: &Wire) -> Vec<(usize, Orientation)> {
    wire.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .map(|shape| (shape_key(shape), shape.orientation()))
        .collect()
}

fn edge_counts(wire: &Wire) -> (usize, usize) {
    edges_of_wire(wire).iter().fold((0, 0), |mut counts, edge| {
        match edge.0.orientation() {
            Orientation::Forward | Orientation::Reversed => counts.0 += 1,
            Orientation::Internal | Orientation::External => counts.1 += 1,
        }
        counts
    })
}

fn stored_faces(shape: &TopoShape) -> Vec<Face> {
    if shape.shape_type() == ShapeType::Face {
        vec![Face(shape.clone())]
    } else {
        crate::topo_tools::faces_of(shape)
    }
}

fn stored_wires(face: &Face) -> Vec<Wire> {
    face.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|shape| shape.shape_type() == ShapeType::Wire)
        .cloned()
        .map(Wire)
        .collect()
}

fn apply_face(face: &Face, context: &mut crate::shape_fix_compose_shell::SharedReShape) -> Face {
    let applied = context.apply(&face.0);
    if applied.is_face() {
        Face(applied)
    } else {
        face.clone()
    }
}

fn drop_empty_wires(face: &Face) -> bool {
    let children = face
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let mut removed = false;
    let kept: Vec<TopoShape> = children
        .into_iter()
        .filter(|child| {
            if child.shape_type() != ShapeType::Wire {
                return true;
            }
            let counts = edge_counts(&Wire(child.clone()));
            let keep = counts.0 != 0 || counts.1 != 0;
            removed |= !keep;
            keep
        })
        .collect();
    if removed {
        face.0
            .tshape
            .write()
            .expect("poisoned TShape lock")
            .children = kept;
    }
    removed
}

impl ShapeFixFace {
    /// Runs the complete OCCT 8.0.0 `ShapeFix_Face::Perform` control flow.
    ///
    /// The returned value is true when a `ShapeExtend_DONE*` status was set.
    /// The resulting face, shell, or compound is available in [`Self::result`].
    pub fn perform(&mut self) -> bool {
        self.status = 0;
        self.map_wires.clear();

        let Some(initial_face) = self.face.clone() else {
            self.result = None;
            return false;
        };
        let mut face = apply_face(&initial_face, &mut self.context);
        self.face = Some(face.clone());

        let saved_precision = self.precision;
        let wire_precision = if self.modes.need_fix_wire() && self.modes.auto_correct_precision() {
            let size = least_edge_size(&face.0);
            let corrected = saved_precision.min(size / 2.0) * 1.00001;
            if saved_precision > corrected && corrected > CONFUSION {
                corrected
            } else {
                saved_precision
            }
        } else {
            saved_precision
        };

        if self.modes.need_fix_wire() {
            let first_modes = self.wire_modes.first_wire_round();
            let wires = stored_wires(&face);
            let mut changed = false;
            for mut wire in wires {
                let (manifold, non_manifold) = edge_counts(&wire);
                if manifold == 0 {
                    if non_manifold == 0 {
                        changed = true;
                        self.status |= DONE5;
                    }
                    continue;
                }
                let before = wire_state(&wire);
                perform_wire_first_round(
                    &mut wire,
                    &face,
                    wire_precision,
                    &mut self.context,
                    &first_modes,
                );
                changed |= before != wire_state(&wire);
            }
            if drop_empty_wires(&face) {
                changed = true;
                self.status |= DONE5;
            }
            if changed {
                self.status |= DONE1;
            }
        }

        self.result = Some(face.0.clone());
        if self.modes.need_fix_periodic_degenerated() {
            self.fix_periodic_degenerated();
        }
        if self.modes.need_fix_missing_seam() && self.fix_missing_seam() {
            self.status |= DONE3;
        }

        let root = self.result.clone().unwrap_or_else(|| face.0.clone());
        let second_modes = self.wire_modes.second_wire_round();
        for source_face in stored_faces(&root) {
            face = apply_face(&source_face, &mut self.context);
            self.face = Some(face.clone());
            let mut removed_segment = false;
            let mut wire_changed = false;

            if self.modes.need_fix_wire() {
                for mut wire in stored_wires(&face) {
                    let (manifold, non_manifold) = edge_counts(&wire);
                    if manifold == 0 {
                        if non_manifold == 0 {
                            wire_changed = true;
                            self.status |= DONE5;
                        }
                        continue;
                    }
                    let before = wire_state(&wire);
                    removed_segment |=
                        perform_wire_second_round(&mut wire, &face, wire_precision, &second_modes);
                    wire_changed |= before != wire_state(&wire);
                }
                if self.modes.need_fix_loop_wires() && split_loop_wires(&face) {
                    wire_changed = true;
                    self.status |= DONE7;
                }
                if drop_empty_wires(&face) {
                    wire_changed = true;
                    self.status |= DONE5;
                }
                if wire_changed {
                    self.status |= DONE1;
                }
            }

            if removed_segment && split_face_after_removed_segment(&face, &mut self.context) {
                self.status |= DONE8;
                face = apply_face(&face, &mut self.context);
                self.face = Some(face.clone());
            }

            if self.fix_wires_two_coinc_edges() {
                self.status |= DONE7;
            }
            if self.modes.need_fix_intersecting_wires() && self.fix_intersecting_wires() {
                self.status |= DONE6;
            }
            if self.modes.need_fix_orientation() && self.fix_orientation() {
                self.status |= DONE2;
            }

            let mut need_split = true;
            if self.fix_add_natural_bound() {
                need_split = false;
                self.status |= DONE5;
            }
            if self.modes.need_fix_split_face()
                && need_split
                && self.map_wires.len() > 1
                && self.fix_split_face()
            {
                self.status |= DONE8;
            }
        }

        self.precision = saved_precision;
        let final_root = self.context.apply(&root);
        self.result = Some(final_root.clone());

        if self.modes.need_fix_small_area_wire() {
            for small_face in stored_faces(&final_root) {
                self.face = Some(small_face);
                if self.fix_small_area_wire(self.modes.is_remove_small_face()) {
                    self.status |= DONE4;
                }
            }
            self.result = Some(self.context.apply(&final_root));
        }

        self.status != 0
    }
}
