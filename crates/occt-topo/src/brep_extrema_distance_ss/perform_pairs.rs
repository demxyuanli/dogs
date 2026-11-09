//! The typed `Perform` pair overloads.

use super::prelude::*;
use super::DistanceSS;
use super::statics::*;

impl DistanceSS {
    /// `Perform(V1, V2, Seq1, Seq2)` (`cxx:780-802`).
    pub fn perform_vertex_vertex(&mut self, s1: &Vertex, s2: &Vertex) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_vertex_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(V, E, Seq1, Seq2)` (`cxx:806-864`).
    pub fn perform_vertex_edge(&mut self, s1: &Vertex, s2: &Edge) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_edge_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(V, F, Seq1, Seq2)` (`cxx:868-925`).
    pub fn perform_vertex_face(&mut self, s1: &Vertex, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_face_impl(
            s1,
            s2,
            &mut seq1,
            &mut seq2,
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

    /// `Perform(E1, E2, Seq1, Seq2)` (`cxx:929-1023`).
    pub fn perform_edge_edge(&mut self, s1: &Edge, s2: &Edge) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        edge_edge_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(E, F, Seq1, Seq2)` (`cxx:1026-1187`).
    pub fn perform_edge_face(&mut self, s1: &Edge, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        edge_face_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(F1, F2, Seq1, Seq2)` (`cxx:1191-1206`).
    pub fn perform_face_face(&mut self, s1: &Face, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        face_face_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }
}
