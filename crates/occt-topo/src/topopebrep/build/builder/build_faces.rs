//! `TopOpeBRepBuild_Builder::BuildFaces` — the surface/face build pass.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuildFaces.cxx:40-107`
//! (OCCT 8.0.0): `BuildFaces(iS, HDS)` `:40-96` and `BuildFaces(HDS)`
//! `:99-107`.
//!
//! The pass builds, for every DS surface, the face of that surface out of the
//! new edges `BuildEdges` left in `myNewEdges`: each edge receives the pcurve
//! the surface/curve interference carries (`BuildTool::PCurve(F, E, CDS, PC)`),
//! is oriented `SCurves.Orientation(TopAbs_IN)` and pushed as a start element of
//! a `WireEdgeSet`, which `FaceBuilder` then turns into face loops and
//! `MakeFaces` into the new faces of the surface.

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topopebrep::build::face_builder::FaceBuilder;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;

use super::Builder;

impl Builder {
    /// `BuildFaces(iS, HDS)` (`cxx:40-96`).
    ///
    /// `myBuildTool.MakeFace(aFace, aTBS)` (`cxx:52`) creates the face the new
    /// edges are attached to; `aTBSTol` (`cxx:49`) is that surface's tolerance
    /// and raises every edge tolerance below it (`cxx:78-81`).
    pub fn build_faces(&mut self, is: i32) {
        // The `cxx:41-89` part only reads the DS, so it runs under the
        // immutable borrow and hands the two owned objects to `cxx:91-93`,
        // which need `&mut self` (`ChangeNewFaces`, `MakeFaces`).
        let Some((a_face, mut fabu)) = (|| {
            let hds = self.data_structure()?;
            let tbs_tol = hds.surface(is).tolerance();
            let mut a_face = TopoShape::new(ShapeType::Face);
            self.my_build_tool.make_face(&mut a_face, hds.surface(is));

            let mut wes = WireEdgeSet::new(&a_face);
            let mut scurves = hds.surface_curves(is);
            while scurves.more() {
                // cxx:64-67
                let ic = scurves.geometry();
                for mut an_edge in self.new_edges(ic).to_vec() {
                    // cxx:77-81: `BRep_Builder::UpdateEdge(aE, aTBSTol)`.
                    if let Some(e) = Edge::wrap(an_edge.clone()) {
                        let tbc_tol = BRepTool::edge_tolerance(&e);
                        if tbc_tol < tbs_tol {
                            GeometryRegistry::global().raise_edge_tolerance(&e.0, tbs_tol);
                        }
                    }
                    // cxx:83-86
                    let ori = scurves.orientation(State::In);
                    self.my_build_tool.set_orientation(&mut an_edge, ori);
                    let pc = scurves.p_curve();
                    self.my_build_tool
                        .p_curve_with_ds(&a_face, &an_edge, hds.curve(ic), pc);
                    wes.add_start_element(&an_edge);
                }
                scurves.next();
            }

            // `TopOpeBRepBuild_FaceBuilder FABU(WES, aFace)` (`cxx:91`): the
            // three-argument constructor defaults `ForceClass` to `false`.
            Some((a_face.clone(), FaceBuilder::with_wire_edge_set(&mut wes, &a_face, false)))
        })() else {
            return;
        };

        // cxx:92-93: `ChangeNewFaces(iS)` is the list `MakeFaces` appends to.
        let mut face_list = self
            .change_new_faces(is)
            .map(std::mem::take)
            .unwrap_or_default();
        self.make_faces(&a_face, &mut fabu, &mut face_list);
        if let Some(slot) = self.change_new_faces(is) {
            *slot = face_list;
        }
    }

    /// `BuildFaces(HDS)` (`cxx:99-107`): one `BuildFaces(iS, HDS)` per DS
    /// surface. `myNewFaces = new HArray1<...>(0, n)` becomes a 1-based array
    /// holding the lists of the surfaces `1..=n`.
    pub fn build_faces_all(&mut self) {
        let n = self.data_structure().map(|h| h.nb_surfaces()).unwrap_or(0);
        self.set_new_faces_len(n as usize);
        for is in 1..=n {
            self.build_faces(is);
        }
    }
}
