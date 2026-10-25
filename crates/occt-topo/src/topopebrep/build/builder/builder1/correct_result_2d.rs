//! `TopOpeBRepBuild_Builder1::CorrectResult2d`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1_2.cxx:34-224`
//! (OCCT 8.0.0).
//!
//! `CorrectFace2d` moves the pcurves of an edge when the 2-D contour of a
//! periodic face is corrected; the corrected copies are collected in
//! `myMapOfCorrect2dEdges`. This function rebuilds the result solid so the
//! shell / face / wire / edge tree points at those copies instead of the
//! originals, then rewrites the DS `NewEdges` lists and the split lists of the
//! section edges so the builder's own bookkeeping matches.
//!
//! The map is only filled by `PerformPieceOn2D` / `PerformONParts`
//! (`Builder1.cxx`, UNPORTED), so today the `Extent() == 0` guard at `cxx:38-42`
//! is the exit this function takes. The rest is translated because it defines
//! what those two callers must produce.

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::tools::map_shapes_of;
use crate::topopebrep::ds::curve_explorer::CurveExplorer;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

use super::super::Builder;

impl Builder {
    /// `CorrectResult2d(aResult)` (`cxx:34-224`): `1` when the result was
    /// rebuilt, `0` when there was nothing to correct.
    pub fn correct_result_2d(&mut self, a_result: &mut TopoShape) -> i32 {
        // `cxx:38-42`
        if self.b1.my_map_of_correct_2d_edges.borrow().extent() == 0 {
            return 0;
        }
        // `cxx:44-47`
        if a_result.shape_type() != ShapeType::Solid {
            return 0;
        }

        // `cxx:49-55`: every edge of both arguments, so the copy below can
        // tell an argument edge from a result edge.
        let mut a_source_edges = IndexedShapeSet::new();
        for s in [self.my_shape1.clone(), self.my_shape2.clone()]
            .into_iter()
            .flatten()
        {
            let edges = map_shapes_of(&s, ShapeType::Edge);
            for i in 1..=edges.extent() {
                if let Some(e) = edges.find_key(i) {
                    a_source_edges.add(e);
                }
            }
        }

        // `cxx:57-58`
        let mut ed_map: IndexedShapeData<TopoShape> = IndexedShapeData::new();
        let bb = TopoBuilder::new();

        // `cxx:59-60`
        let mut a_solid = a_result.empty_copied();

        // `cxx:62-161`
        for s in explorer_shapes(a_result, ShapeType::Shell) {
            let mut a_shell = s.empty_copied();
            for f in explorer_shapes(&s, ShapeType::Face) {
                let f_ori = f.orientation();
                let mut a_face = f.empty_copied();
                for w in explorer_shapes(&f, ShapeType::Wire) {
                    let w_ori = w.orientation();
                    let mut a_wire = w.empty_copied();

                    for e in explorer_shapes(&w, ShapeType::Edge) {
                        // `cxx:99-106`: an edge already copied for this result
                        // is reused as is.
                        if let Some(mut an_edge) = ed_map.find_from_key(&e).cloned() {
                            an_edge.set_orientation(e.orientation());
                            bb.add(&mut a_wire, &an_edge);
                            continue;
                        }
                        // `cxx:109-116`: a corrected edge is reused and
                        // remembered.
                        if let Some(corrected) = self
                            .b1
                            .my_map_of_correct_2d_edges
                            .borrow()
                            .find_from_key(&e)
                            .cloned()
                        {
                            let mut an_edge = corrected;
                            an_edge.set_orientation(e.orientation());
                            bb.add(&mut a_wire, &an_edge);
                            ed_map.add(&e, an_edge);
                            continue;
                        }

                        // `cxx:118-142`
                        let mut an_edge;
                        if a_source_edges.contains(&e) {
                            // A copy, so the argument edge is never modified.
                            an_edge = e.empty_copied();
                            ed_map.add(&e, an_edge.clone());

                            let free = an_edge.free();
                            an_edge.set_free(true);
                            for v in explorer_shapes(&e, ShapeType::Vertex) {
                                bb.add(&mut an_edge, &v);
                            }
                            an_edge.set_free(free);
                        } else {
                            an_edge = e.clone();
                        }

                        an_edge.set_orientation(e.orientation());
                        bb.add(&mut a_wire, &an_edge);
                    }

                    // `cxx:147-150`
                    a_wire.set_orientation(w_ori);
                    let closed = BRepTool::is_closed(&a_wire);
                    a_wire.set_closed(closed);
                    bb.add(&mut a_face, &a_wire);
                }

                // `cxx:152-154`
                a_face.set_orientation(f_ori);
                bb.add(&mut a_shell, &a_face);
            }

            // `cxx:156-159`
            a_shell.set_orientation(s.orientation());
            let closed = BRepTool::is_closed(&a_shell);
            a_shell.set_closed(closed);
            bb.add(&mut a_solid, &a_shell);
        }
        *a_result = a_solid;

        // `cxx:163-186`: update the section curves' new edges.
        let curves: Vec<i32> = match self.data_structure() {
            Some(hds) => {
                let mut cex = CurveExplorer::with_ds(hds.ds(), true);
                let mut idx = Vec::new();
                while cex.more() {
                    idx.push(cex.index());
                    cex.next();
                }
                idx
            }
            None => Vec::new(),
        };
        for ic in curves {
            let lse = self.change_new_edges(ic).clone();
            let mut corr = Vec::with_capacity(lse.len());
            for e in &lse {
                match ed_map.find_from_key(e) {
                    Some(new_e) => corr.push(new_e.clone()),
                    None => corr.push(e.clone()),
                }
            }
            let lse = self.change_new_edges(ic);
            lse.clear();
            lse.extend(corr);
        }

        // `cxx:188-222`: update the split lists of the section edges.
        let nes = match self.data_structure() {
            Some(hds) => hds.ds().nb_section_edges(),
            None => 0,
        };
        for i in 1..=nes {
            let es = match self.data_structure() {
                Some(hds) => hds.ds().section_edge(i, true),
                None => None,
            };
            let Some(es) = es else {
                continue;
            };
            if es.is_null() {
                continue;
            }
            // `cxx:200-202`: `TopAbs_State(j)`, `0` is `IN`, `1` `OUT`,
            // `2` `ON`.
            for st in [State::In, State::Out, State::On] {
                let lse = self.change_split(&es, st).clone();
                let mut corr = Vec::with_capacity(lse.len());
                for e in &lse {
                    match ed_map.find_from_key(e) {
                        Some(new_e) => corr.push(new_e.clone()),
                        None => corr.push(e.clone()),
                    }
                }
                let lse = self.change_split(&es, st);
                lse.clear();
                lse.extend(corr);
            }
        }
        1
    }

    /// `myMapOfCorrect2dEdges.Add(E, ECorrected)` — the writers are
    /// `PerformPieceIn2D` / `PerformONParts` (`Builder1.cxx`, UNPORTED).
    #[allow(dead_code)]
    pub fn bind_correct_2d_edge(&mut self, e: &TopoShape, corrected: &TopoShape) {
        self
            .b1
            .my_map_of_correct_2d_edges
            .borrow_mut()
            .add(e, corrected.clone());
    }
}
