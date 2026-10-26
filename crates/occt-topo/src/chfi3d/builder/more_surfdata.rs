//! `ChFi3d_Builder::MoreSurfdata`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:4601-4705`
//! (OCCT 8.0.0), `const` in `ChFi3d_Builder.hxx:780-801`.
//!
//! Decides whether the intersection at the end concerns more than one
//! `ChFiDS_SurfData`: the last-but-one `SurfData` of the stripe must have a
//! common point that sits on one of the two arcs where the end face meets the
//! fillet, and that point must be within `1e-4` of the corner's own point on
//! the same arc.
//!
//! The OCCT body takes references into the stripe's `SurfData` sequence
//! (`ChFiDS_CommonPoint& CV1 = Fd->ChangeVertex(...)`) and only reads them; the
//! port copies the two `CommonPoint`s, which is the same values without holding
//! the stripe lock across `FindFace`.
//!
//! `num1` is computed by the source (`cxx:4629-4634`) and never read - `num2`
//! is the one the lookups use; the port keeps both to mirror the branch.

use std::sync::Arc;

use crate::chfi_ds::SurfData;
use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::is_same;

use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `MoreSurfdata(Index)` (`cxx:4601-4705`).
    ///
    /// `false` when the vertex has no stripe, when the stripe's `SurfData`
    /// sequence cannot be read, or when neither of the last-but-one common
    /// points matches the corner point (the source's `oksurf`).
    pub fn more_surfdata(&self, index: usize) -> bool {
        // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);`
        let Some(vtx) = self.my_v_data_map.find_key(index).cloned() else {
            return false;
        };
        // `It.Initialize(myVDataMap(Index)); stripe = It.ChangeValue();`
        let Some(stripe) = self
            .my_v_data_map
            .find_from_index(index)
            .and_then(|l| l.first())
            .cloned()
        else {
            return false;
        };
        // `SeqFil = stripe->ChangeSetOfSurfData()->ChangeSequence();`
        let seq: Vec<Arc<SurfData>> = {
            let Ok(guard) = stripe.read() else {
                return false;
            };
            guard.set_of_surf_data().to_vec()
        };
        let nbsurf = seq.len();

        // `int sens = 0; num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
        //  bool isfirst = (sens == 1);`
        let mut sens = 0;
        let Ok(num) = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens) else {
            return false;
        };
        let Ok(num) = usize::try_from(num) else {
            return false;
        };
        let is_first = sens == 1;

        // `ChFiDS_CommonPoint& CV1 = Fd->ChangeVertex(isfirst, 1);` (+ CV2).
        let Some(fd) = num.checked_sub(1).and_then(|i| seq.get(i)) else {
            return false;
        };
        let cv1 = fd.vertex(is_first, 1).clone();
        let cv2 = fd.vertex(is_first, 2).clone();

        // `inters = FindFace(Vtx, CV1, CV2, Fv);` - `Fv` is the face at end.
        let Some(fv) = self.find_face(&vtx, &cv1, &cv2, None) else {
            return false;
        };

        // `if (sens == 1) { num1 = 1; num2 = num1 + 1; }
        //  else { num1 = nbsurf; num2 = num1 - 1; }`
        let _num1 = if sens == 1 { 1usize } else { nbsurf };
        let num2 = if sens == 1 { 2usize } else { nbsurf - 1 };

        if nbsurf == 1 {
            return false;
        }

        // `arc1` / `arc2`: the two edges of `myVEMap(Vtx)` that belong to `Fv`
        // (`cxx:4647-4664`).
        let vertex_edges: Vec<TopoShape> = self.my_ve_map.get_from_key(&vtx).to_vec();

        let mut arc1: Option<Edge> = None;
        for item in &vertex_edges {
            if let Some(e) = Edge::wrap(item.clone()) {
                if crate::chfi3d::builder_c1::contain_e(&fv, &e) {
                    arc1 = Some(e);
                    break;
                }
            }
        }
        let mut arc2: Option<Edge> = None;
        for item in &vertex_edges {
            if let Some(e) = Edge::wrap(item.clone()) {
                if !crate::chfi3d::builder_c1::contain_e(&fv, &e) {
                    continue;
                }
                let same_as_arc1 = match &arc1 {
                    Some(a) => is_same(&a.0, &e.0),
                    None => false,
                };
                if !same_as_arc1 {
                    arc2 = Some(e);
                    break;
                }
            }
        }

        // `Fd1 = SeqFil.ChangeValue(num2); CV3 = Fd1->ChangeVertex(isfirst, 1);
        //  CV4 = Fd1->ChangeVertex(isfirst, 2);`
        let Some(fd1) = num2.checked_sub(1).and_then(|i| seq.get(i)) else {
            return false;
        };
        let cv3 = fd1.vertex(is_first, 1).clone();
        let cv4 = fd1.vertex(is_first, 2).clone();

        // `Arc().IsSame(...)`, with a null arc answering `false` like OCCT's
        // `IsSame` on a null shape.
        let arc_is = |candidate: Option<&TopoShape>, target: Option<&Edge>| match (candidate, target)
        {
            (Some(c), Some(t)) => is_same(c, &t.0),
            _ => false,
        };

        let mut oksurf = false;
        // `if (CV3.IsOnArc())` (`cxx:4670-4688`).
        if cv3.is_on_arc() {
            if arc_is(cv3.arc(), arc1.as_ref()) {
                if cv1.point().distance(cv3.point()) < 1.0e-4 {
                    oksurf = true;
                }
            } else if arc_is(cv3.arc(), arc2.as_ref()) {
                if cv2.point().distance(cv3.point()) < 1.0e-4 {
                    oksurf = true;
                }
            }
        }
        // `if (CV4.IsOnArc())` (`cxx:4689-4702`).
        if cv4.is_on_arc() {
            if cv1.point().distance(cv4.point()) < 1.0e-4 {
                oksurf = true;
            } else if arc_is(cv4.arc(), arc2.as_ref()) {
                if cv2.point().distance(cv4.point()) < 1.0e-4 {
                    oksurf = true;
                }
            }
        }
        oksurf
    }
}
