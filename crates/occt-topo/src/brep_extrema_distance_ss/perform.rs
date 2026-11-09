//! `Perform(S1, S2, Box1, Box2)` and the box-derived accessors.

use super::prelude::*;
use super::DistanceSS;
use super::trim_infinite::*;
use super::statics::*;

impl DistanceSS {
    /// `Perform(S1, S2, Box1, Box2)` (`cxx:656-775`).
    pub fn perform(
        &mut self,
        the_s1: &TopoShape,
        the_s2: &TopoShape,
        the_box1: &BndBox,
        the_box2: &BndBox,
    ) {
        // `cxx:660-662`.
        self.my_seq_sol_shape1.clear();
        self.my_seq_sol_shape2.clear();
        self.my_modif = false;

        // `cxx:664-671`: the box–box gap is the early exit. `Bnd_Box::Distance`
        // and `brep_extrema::shape_bbox_distance` compute the same axis-gap
        // norm for finite boxes.
        let a_bb_dist = crate::brep_extrema::shape_bbox_distance(the_box1, the_box2);
        if a_bb_dist - self.my_dst_ref > self.my_eps {
            return;
        }

        // `cxx:673-772`: the type switch. Every arm borrows the two sequences
        // mutably, so the arms are dispatched explicitly instead of through a
        // closure.
        match (the_s1.shape_type(), the_s2.shape_type()) {
            (ShapeType::Vertex, ShapeType::Vertex) => {
                let (v1, v2) = match (Vertex::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_vertex(&v1, &v2);
            }
            (ShapeType::Vertex, ShapeType::Edge) => {
                let (v1, e2) = match (Vertex::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_edge(&v1, &e2);
            }
            (ShapeType::Vertex, ShapeType::Face) => {
                let (v1, f2) = match (Vertex::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_face(&v1, &f2);
            }
            (ShapeType::Edge, ShapeType::Vertex) => {
                let (e1, v2) = match (Edge::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                vertex_edge_impl(
                    &v2,
                    &e1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Edge, ShapeType::Edge) => {
                let (e1, e2) = match (Edge::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:710-723`: `TRIM_INFINIT_EDGE`.
                let (trim, trim1, trim2) = trim_infinit_edge(&e1, &e2);
                let (t1, t2) = match trim {
                    Some(t) => (
                        if trim1 { t.clone() } else { e1.clone() },
                        if trim2 { t } else { e2.clone() },
                    ),
                    None => (e1.clone(), e2.clone()),
                };
                self.perform_edge_edge(&t1, &t2);
            }
            (ShapeType::Edge, ShapeType::Face) => {
                let (e1, f2) = match (Edge::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:725-734`: `TRIM_INFINIT_FACE(aE1, aF2, ...)`.
                let (trim, is_inf) = trim_infinit_face(&e1.0, &f2.0);
                let f2 = if is_inf {
                    match trim {
                        Some(f) => f,
                        None => f2.clone(),
                    }
                } else {
                    f2.clone()
                };
                self.perform_edge_face(&e1, &f2);
            }
            (ShapeType::Face, ShapeType::Vertex) => {
                let (f1, v2) = match (Face::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                vertex_face_impl(
                    &v2,
                    &f1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    self.my_flag,
                    self.my_algo,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Face, ShapeType::Edge) => {
                let (f1, e2) = match (Face::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:753-759`: `TRIM_INFINIT_FACE(aF1, aE2, ...)`.
                let (trim, is_inf) = trim_infinit_face(&f1.0, &e2.0);
                let f1 = if is_inf {
                    match trim {
                        Some(f) => f,
                        None => f1.clone(),
                    }
                } else {
                    f1.clone()
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                edge_face_impl(
                    &e2,
                    &f1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Face, ShapeType::Face) => {
                let (f1, f2) = match (Face::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_face_face(&f1, &f2);
            }
            // `cxx:771-772`.
            _ => {}
        }
    }

    /// `IsDone()` (`hxx:69`): `myModif`.
    pub fn is_done(&self) -> bool {
        self.my_modif
    }

    /// `DistValue()` (`hxx:72`).
    pub fn dist_value(&self) -> f64 {
        self.my_dst_ref
    }

    /// `Seq1Value()` (`hxx:75`).
    pub fn seq1_value(&self) -> &[SolutionElem] {
        &self.my_seq_sol_shape1
    }

    /// `Seq2Value()` (`hxx:78`).
    pub fn seq2_value(&self) -> &[SolutionElem] {
        &self.my_seq_sol_shape2
    }
}
