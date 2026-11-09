//! The static `AddPoints` of `LocOpe_CSIntersector`
//! (`LocOpe_CSIntersector.cxx:517-589`).
//!
//! Copies the hits the face intersector holds into the caller's sequence,
//! turning `IntCurveSurface_In` into `TopAbs_FORWARD` (`cxx:530-550`),
//! `IntCurveSurface_Out` into `TopAbs_REVERSED` (`cxx:553-565`) and
//! `IntCurveSurface_Tangent` into `TopAbs_EXTERNAL` (`cxx:566-568`), then
//! inserting each point before the first recorded point with a greater
//! parameter (`cxx:569-587`).

use crate::abs::Orientation;
use crate::int_curves_face::{FaceIntersector, Transition};
use crate::shape::Face;

use super::super::pnt_face::PntFace;

pub(super) fn add_points(the_int: &FaceIntersector, seq: &mut Vec<PntFace>, theface: &Face) {
    let newpnt = the_int.nb_pnt();
    for j in 1..=newpnt {
        let thept = the_int.pnt(j);
        let param = the_int.w_parameter(j);
        let paramu = the_int.u_parameter(j);
        let paramv = the_int.v_parameter(j);

        // `cxx:530-568`.
        let theor = match the_int.transition(j) {
            Transition::In => Orientation::Forward,
            Transition::Out => Orientation::Reversed,
            Transition::Tangent => Orientation::External,
        };
        let newpt = PntFace::new(thept, theface.clone(), theor, param, paramu, paramv);

        // `cxx:570-586`.
        let mut k = 0usize;
        while k < seq.len() {
            if seq[k].parameter() - param > 0.0 {
                break;
            }
            k += 1;
        }
        seq.insert(k, newpt);
    }
}
