//! `BRepOffset_MakeLoops::Build` (`BRepOffset_MakeLoops.cxx:42-179`).

use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_algo_loop::BRepAlgoLoop;
use crate::shape::{Edge, Face, TopoShape};

use super::{substitute_vertices, BRepOffsetMakeLoops};

impl BRepOffsetMakeLoops {
    /// `Build(LF, AsDes, Image, theImageVV, theRange)` (`cxx:42-179`): per face,
    /// feed the already-cut images of its descendant edges (or the edges
    /// themselves when never cut) into a `BRepAlgo_Loop`, rebuild the face and
    /// record both the new faces and the new edges.
    ///
    /// `theImageVV` is only read: `BRepAlgo_Loop::SetImageVV` copies it into the
    /// loop object, so the caller's image is not written back.
    ///
    /// `as_des` is `&mut` here where the reference takes
    /// `const Handle(BRepAlgo_AsDes)&`: the shared handle lets OCCT call
    /// `AsDes->Replace` from a const reference, the port's owned `AsDes` cannot.
    pub fn build(
        &mut self,
        lf: &[TopoShape],
        as_des: &mut AsDes,
        image: &mut BRepAlgoImage,
        the_image_vv: &BRepAlgoImage,
    ) {
        let mut loops = BRepAlgoLoop::new();
        // `Loops.VerticesForSubstitute(myVerVerMap); Loops.SetImageVV(...)`
        // (`cxx:50-51`).
        loops.set_vertices_for_substitute(self.my_ver_ver_map.clone());
        loops.set_image_vv(the_image_vv.clone());

        // `cxx:55-58`: one iteration per face of `LF`.
        for f in lf {
            let face = Face(f.clone());
            // `Loops.Init(F)` (`cxx:64`).
            loops.init(&face);

            // `LE = AsDes->Descendant(F)` (`cxx:72`).
            let le: Vec<TopoShape> = as_des.descendant(&face.0).to_vec();
            // `NCollection_List<TopoDS_Shape> AddedEdges` (`cxx:69`).
            let mut added_edges: Vec<TopoShape> = Vec::new();

            for e_shape in &le {
                let e = Edge(e_shape.clone());
                if image.has_image(&e.0) {
                    // `cxx:74-88`: `E` was already cut on another face; add
                    // every cut edge, oriented like `E`.
                    for ce in image.image(&e.0) {
                        let ce = ce.oriented(e.0.orientation());
                        loops.add_const_edge(&Edge(ce));
                    }
                } else {
                    // `cxx:90-93`: `AddEdge(E, AsDes->Descendant(E))`.
                    let descendants: Vec<TopoShape> = as_des.descendant(&e.0).to_vec();
                    loops.add_edge(&e, &descendants);
                    added_edges.push(e.0.clone());
                }
            }

            // `cxx:98-99`: unwind.
            loops.perform();
            loops.wires_to_faces();

            // `cxx:103-113`: `F => NewFaces`.
            let nf: Vec<TopoShape> = loops.new_faces().to_vec();
            let _ = image.bind_list(&face.0, &nf);

            // `cxx:108-124`: `E => NewEdges`, for the edges handed to
            // `AddEdge` above.
            for e_shape in &added_edges {
                let e = Edge(e_shape.clone());
                let Some(loop_ne) = loops.new_edges(&e).cloned() else {
                    continue;
                };
                if image.has_image(&e.0) {
                    let _ = image.add_list(&e.0, &loop_ne);
                } else {
                    let _ = image.bind_list(&e.0, &loop_ne);
                }
            }
        }

        // `cxx:125-128`: `Loops.GetVerticesForSubstitute(myVerVerMap)`.
        self.my_ver_ver_map = loops.vertices_for_substitute().clone();
        if self.my_ver_ver_map.is_empty() {
            return;
        }

        // `cxx:131-179`: the substitution pass over `LF`.
        let map = self.my_ver_ver_map.clone();
        substitute_vertices(lf, image, as_des, &map);
    }
}
