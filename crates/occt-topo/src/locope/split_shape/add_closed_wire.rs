//! `LocOpe_SplitShape::AddClosedWire` (`LocOpe_SplitShape.cxx:658-732`).
//!
//! A closed splitting wire either becomes a hole of the descendant face it sits
//! in, or it cuts that face in two. The reference decides which by classifying
//! the infinite point of a throwaway face built over the wire alone: a face
//! whose only wire bounds the *outside* means the wire encloses the material,
//! so the wire is reversed before it is used.

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::locope::split_shape::is_inside::{is_inside_wire, is_inside_wire_pair};
use crate::locope::split_shape::{empty_copied_face, SplitShape};
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::explorer_shapes;

impl SplitShape {
    /// `AddClosedWire(W, F)` (`cxx:658-732`).
    pub fn add_closed_wire(&mut self, w: &Wire, f: &Face) -> bool {
        // `cxx:662-677`: the descendant of `F` that contains `W`.
        let lf = self.map_of(&f.0);
        let mut found: Option<usize> = None;
        for (i, s) in lf.iter().enumerate() {
            let Some(fac) = Face::wrap(s.clone()) else {
                continue;
            };
            // `cxx:669-672`: `IsInside(fac, W)`.
            if is_inside_wire(&fac, w) {
                found = Some(i);
                break;
            }
        }
        let Some(found) = found else {
            return false;
        };

        let b = TopoBuilder::new();

        // `cxx:682-688`.
        let mut or_wire = w.0.orientation();
        let mut new_face = empty_copied_face(f);
        new_face.0.set_orientation(Orientation::Forward);
        b.add(&mut new_face.0, &w.0);

        // `cxx:691-702`: a wire that encloses the infinite point bounds a hole,
        // so it is reversed and the face is rebuilt over the reversed wire.
        let classif = BrepTopAdaptorFClass2d::new_with_wires(f, std::slice::from_ref(w), PCONFUSION);
        if classif.perform_infinite_point() == FaceState::In {
            new_face = empty_copied_face(f);
            new_face.0.set_orientation(Orientation::Forward);
            or_wire = or_wire.reversed();
            b.add(&mut new_face.0, &w.0.oriented(or_wire));
        }

        // `cxx:704-706`: `FaceRef` is the descendant, taken out of the list.
        let Some(face_ref) = Face::wrap(lf[found].clone()) else {
            return false;
        };
        let mut face_ref_forward = face_ref.clone();
        face_ref_forward.0.set_orientation(Orientation::Forward);
        let mut new_lf: Vec<TopoShape> = lf[..found].to_vec();
        new_lf.extend(lf[found + 1..].iter().cloned());

        // `cxx:708-710`.
        let mut new_ref = empty_copied_face(&face_ref);
        new_ref.0.set_orientation(Orientation::Forward);

        // `cxx:715-727`: every existing wire of the descendant goes to the side
        // of the split it belongs to; `W` closes `new_ref` with the orientation
        // opposite to the one the new face took.
        for sub in explorer_shapes(&face_ref_forward, ShapeType::Wire) {
            let Some(wir) = Wire::wrap(sub) else {
                continue;
            };
            if is_inside_wire_pair(f, &wir, w) {
                b.add(&mut new_face.0, &wir.0);
            } else {
                b.add(&mut new_ref.0, &wir.0);
            }
        }
        // `cxx:726`.
        b.add(&mut new_ref.0, &w.0.oriented(or_wire.reversed()));

        // `cxx:728-729`.
        new_lf.push(new_ref.0);
        new_lf.push(new_face.0);
        self.map_bind(&f.0, new_lf);
        true
    }
}
