//! `ChFi3d_Builder::FindFace`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:4533-4598`
//! (OCCT 8.0.0), the two overloads; the first one (`:4533-4540`) passes a null
//! `Favoid`.
//!
//! The answer is the one face shared by the two arcs of the two common points,
//! skipping `Favoid`. OCCT returns a `bool` and writes the face through an out
//! parameter; the port answers `Option<Face>` (`None` is `Found == false`).
//!
//! Divergence that is not a heuristic: the `OCCT_DEBUG` block that reports
//! "change of face on vertex" and checks `myVFMap(V)` is not translated, so the
//! `V` parameter has no reader - it is kept for the OCCT signature.

use crate::chfi_ds::CommonPoint;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_same;

use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `FindFace(V, P1, P2, Fv, Favoid)` (`cxx:4542-4598`).
    ///
    /// `false` (here `None`) when either common point is not on an arc or when
    /// no face carries both arcs. A missing key in `myEFMap` answers the empty
    /// slice in the port (`Standard_NoSuchObject` in OCCT), so it contributes no
    /// candidate.
    pub fn find_face(
        &self,
        _v: &TopoShape,
        p1: &CommonPoint,
        p2: &CommonPoint,
        f_avoid: Option<&Face>,
    ) -> Option<Face> {
        // `if (!(P1.IsOnArc() && P2.IsOnArc())) return false;`
        if !(p1.is_on_arc() && p2.is_on_arc()) {
            return None;
        }
        let Some(arc1) = p1.arc() else {
            return None;
        };
        let Some(arc2) = p2.arc() else {
            return None;
        };

        for item in self.my_ef_map.get_from_key(arc1) {
            let Some(fv) = Face::wrap(item.clone()) else {
                continue;
            };
            if let Some(avoid) = f_avoid {
                if is_same(&fv.0, &avoid.0) {
                    continue;
                }
            }
            for candidate in self.my_ef_map.get_from_key(arc2) {
                if is_same(candidate, &fv.0) {
                    return Some(fv);
                }
            }
        }
        None
    }
}
