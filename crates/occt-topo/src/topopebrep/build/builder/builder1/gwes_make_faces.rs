//! `TopOpeBRepBuild_Builder1::GWESMakeFaces` — the third virtual override.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx`
//! (OCCT 8.0.0, 1373-1412).
//!
//! The base member builds the faces of the wire/edge set; the override then
//! rewrites every built face twice:
//!
//! 1. `myIsKPart == 4` (`KP_ISSOSO`, `k_part/mod.rs:49`) replaces each face by
//!    `TopOpeBRepBuild_Tools::NormalizeFace` (`:1381-1389`), otherwise the list
//!    is copied (`:1391-1394`).
//! 2. Every face of that list is replaced by
//!    `TopOpeBRepBuild_Tools::CorrectFace2d` fed with `mySourceShapes` and
//!    `myMapOfCorrect2dEdges` (`:1397-1411`).
//!
//! Both tools are ported (`build::tools::normalize_face`,
//! `build::tools::correct_face_2d`). OCCT's tools cannot fail — they answer a
//! `TopoDS_Shape&` out-parameter and `NormalizeFace` always assigns it
//! (`Tools.cxx:858-909`) — so a face the port cannot unwrap keeps the shape it
//! was given instead of vanishing from the result.

use crate::shape::{Face, TopoShape};
use crate::topopebrep::build::tools::{correct_face_2d, normalize_face};
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;

use super::super::k_part::KP_ISSOSO;
use super::super::Builder;

impl Builder {
    /// `Builder1::GWESMakeFaces(FF, WES, LOF)` (`cxx:1373-1412`).
    pub(crate) fn b1_gwes_make_faces(
        &mut self,
        ff: &TopoShape,
        wes: &mut WireEdgeSet,
        lof: &mut Vec<TopoShape>,
    ) {
        // `cxx:1378`
        self.base_gwes_make_faces(ff, wes, lof);

        // `cxx:1380-1395`
        let mut corr_lof: Vec<TopoShape> = Vec::with_capacity(lof.len());
        if self.my_is_k_part == KP_ISSOSO {
            for f in lof.iter() {
                match normalize_face(f) {
                    Some(corr) => corr_lof.push(corr),
                    None => corr_lof.push(f.clone()),
                }
            }
        } else {
            corr_lof.extend(lof.iter().cloned());
        }

        // `cxx:1397-1411`: the 2d correction keeps the order of `corr_lof`.
        let mut corr_lof1: Vec<TopoShape> = Vec::with_capacity(corr_lof.len());
        for f in corr_lof.iter() {
            match Face::wrap(f.clone()) {
                Some(face) => corr_lof1.push(correct_face_2d(
                    &face,
                    &self.b1.my_source_shapes,
                    self.b1.my_map_of_correct_2d_edges.clone(),
                )),
                None => corr_lof1.push(f.clone()),
            }
        }

        // `cxx:1410-1411`
        lof.clear();
        lof.extend(corr_lof1);
    }
}
