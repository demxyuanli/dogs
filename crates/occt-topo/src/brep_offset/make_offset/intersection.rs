//! `BRepOffset_MakeOffset::Intersection3D` (`BRepOffset_MakeOffset.cxx:
//! 2855-2939`) and `Intersection2D` (`:2941-2989`).
//!
//! `Intersection3D` runs the 3D face-face phase of the general path:
//!
//! * `:2866-2879` the Complete planar arm returns untouched (that geometry is
//!   built later by `BuildShellsCompleteInter`).
//! * `:2881-2882` `MakeList(OffsetFaces, myInitOffsetFace, myFaces)`.
//! * `:2884-2899` caps with `GeomAbs_Arc`: `Inter.ContextIntByArc(myFaces,
//!   InSide, myAnalyse, myInitOffsetFace, myInitOffsetEdge)` where
//!   `InSide = (myOffset < 0)` ("Temporary", `:2887`).
//! * `:2901-2919` `myInter`: `Inter.CompletInt(OffsetFaces,
//!   myInitOffsetFace)`, then, for `GeomAbs_Intersection` only,
//!   `BRepOffset_Tool::CorrectOrientation(myFaceComp, Inter.NewEdges(),
//!   myAsDes, myInitOffsetFace, myOffset)`.
//! * `:2921-2932` otherwise `Inter.ConnexIntByArc(OffsetFaces, myFaceComp,
//!   myAnalyse, myInitOffsetFace)`.
//!
//! `Intersection2D` (`:2941-2989`) runs `BRepOffset_Inter2d::Compute` on every
//! touched face with `myTol`, then `FuseVertices(aDMVV, myAsDes, myImageVV)`.
//!
//! The port carries `Inter` as the returned value: the reference creates
//! `BRepOffset_Inter3d Inter(myAsDes, Side, myTol)` in `MakeOffsetShape` and
//! reads `Inter.TouchedFaces()` / `Inter.NewEdges()` afterwards, while the
//! port's `BRepOffsetInter3d` owns the `AsDes` it mutates. `Intersection3D`
//! therefore hands back the two maps plus the `AsDes` it ends up holding, and
//! `MakeOffsetShape` keeps them as the `Modif` / `NewEdges` pair.
//!
//! UNPORTED:
//!
//! * `BuildSplitsOfTrimmedFaces` (`_1.cxx:9497-9512`) is not called here; it
//!   only appears in `MakeLoops` / `MakeFaces` (see `loop_rebuild.rs`).
//! * `Message_ProgressScope` / `!More()` -> `BRepOffset_UserBreak`
//!   (`:2874`, `:2907-2911`, `:2927-2931`, `:2966-2971`): no progress
//!   machinery in the port.
//! * the `OCCT_DEBUG` blocks (`:2858-2864`, `:2934-2938`, `:2948-2955`,
//!   `:2979-2986`).

use occt_core::kernel::geomabs::JoinType;

use crate::abs::State;
use crate::brep_offset::BRepOffsetTool;
use crate::brep_offset::{BRepOffsetInter2d, BRepOffsetInter3d, Dmvv};
use crate::shape::Face;

use super::shape_prep::{is_solid, make_list};
use super::BRepOffsetMakeOffset;
use crate::brep_offset::inter3d::ShapeIndexedMap;

/// `Intersection3D(Inter, theRange)` (`cxx:2855-2939`).
///
/// Answers `(Inter.TouchedFaces(), Inter.NewEdges())`; the `AsDes` the
/// reference shares through `Inter` is written back into the object.
pub fn intersection3_d(this: &mut BRepOffsetMakeOffset) -> (ShapeIndexedMap, ShapeIndexedMap) {
    // `cxx:2868-2879`: the Complete planar arm builds its shells elsewhere.
    let planar_complete = this.my_inter
        && this.my_join == JoinType::Intersection
        && this.my_is_planar
        && !this.my_thickening
        && this.my_faces.extent() == 0
        && this.my_shape.as_ref().is_some_and(is_solid);
    if planar_complete {
        return (ShapeIndexedMap::new(), ShapeIndexedMap::new());
    }

    // `cxx:2881-2882`.
    let mut offset_faces: Vec<crate::shape::TopoShape> = Vec::new();
    make_list(&mut offset_faces, &this.my_init_offset_face, &this.my_faces);

    // `cxx:2868`: `TopAbs_State Side` of `MakeOffsetShape` is `TopAbs_IN` for a
    // non-negative offset and `TopAbs_OUT` otherwise (`:872-876`).
    let side = if this.my_offset < 0.0 {
        State::Out
    } else {
        State::In
    };
    let mut inter = BRepOffsetInter3d::new(std::mem::take(&mut this.my_as_des), side, this.my_tol);

    // `cxx:2884-2899`: caps, arc join only.
    if this.my_faces.extent() != 0 {
        // `bool InSide = (myOffset < 0.); // Temporary` (`:2887`).
        let in_side = this.my_offset < 0.0;
        if this.my_join == JoinType::Arc {
            let context_faces: Vec<_> = this.my_faces.iter().cloned().collect();
            inter.context_int_by_arc(
                &context_faces,
                in_side,
                &this.my_analyse,
                &this.my_init_offset_face,
                &mut this.my_init_offset_edge,
            );
        }
    }

    if this.my_inter {
        // `cxx:2904-2906`.
        inter.complet_int(&offset_faces, &this.my_init_offset_face);
        // `cxx:2908-2919`.
        if this.my_join == JoinType::Intersection {
            let new_edges = inter.new_edges().clone();
            BRepOffsetTool::correct_orientation(
                &this.my_face_comp.0,
                &new_edges,
                inter.as_des_mut(),
                &this.my_init_offset_face,
                this.my_offset,
            );
        }
    } else {
        // `cxx:2921-2932`: only between neighbour faces.
        inter.connex_int_by_arc(
            &this.my_face_comp.0,
            &this.my_analyse,
            &this.my_init_offset_face,
        );
    }

    let touched = inter.touched_faces().clone();
    let new_edges = inter.new_edges().clone();
    this.my_as_des = inter.into_as_des();
    (touched, new_edges)
}

impl BRepOffsetMakeOffset {
    /// `Intersection2D(Modif, NewEdges, theRange)` (`cxx:2941-2989`).
    pub fn intersection2_d(
        &mut self,
        modif: &ShapeIndexedMap,
        new_edges: &crate::brep_offset::inter3d::ShapeIndexedMap,
    ) {
        // `cxx:2962-2963`.
        let mut a_dmvv = Dmvv::new();
        let new_edge_keys: std::collections::HashSet<_> = new_edges
            .iter()
            .map(crate::topo_tools::shape_map_key)
            .collect();

        // `cxx:2965-2974`.
        for i in 1..=modif.extent() {
            let Some(shape) = modif.value(i) else {
                continue;
            };
            let Some(f) = Face::wrap(shape.clone()) else {
                continue;
            };
            BRepOffsetInter2d::compute(
                &mut self.my_as_des,
                &f,
                &new_edge_keys,
                self.my_tol,
                &self.my_edge_int_edges,
                &mut a_dmvv,
            );
        }

        // `cxx:2976`.
        let _ = crate::brep_offset::inter2d::fuse_vertices(
            &a_dmvv,
            &mut self.my_as_des,
            &mut self.my_image_vv,
        );
    }
}
