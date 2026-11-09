//! The face bookkeeping of `BRepOffset_MakeOffset`:
//!
//! * `SetFaces` (`BRepOffset_MakeOffset.cxx:795-815`)
//! * `SetFacesWithOffset` (`:818-833`)
//! * `BuildFaceComp` (`:1873-1890`)
//! * `ReplaceRoots` (`:1825-1867`)
//!
//! `SetFaces` fills `myFaces` with the faces to offset and seeds the two
//! history images with a self binding, `SetFacesWithOffset` moves the
//! `SetOffsetOnFace` overrides from the original face to the plane that
//! `IsPlanar` substituted for it, `BuildFaceComp` builds `myFaceComp` (the
//! shape `BRepOffset_Analyse::Perform` analyses), and `ReplaceRoots` puts the
//! original edges / vertices back as image roots once the offset is built.
//!
//! Port guards (no OCCT branch changed):
//!
//! * `TopoDS::Face(myOriginalFaces(ii))` (`:799`) and
//!   `TopoDS::Face(*aPlanface)` (`:804`) force the shape type; the port wraps
//!   through [`Face::wrap`] and skips an item of the wrong type, which OCCT
//!   cannot produce through `AddFace`.
//! * `BRepAlgo_Image::Bind` / `SetRoot` raise `Standard_ConstructionError` on
//!   an already imaged root (`BRepAlgo_Image.cxx`); `SetFaces` runs on a fresh
//!   image after `Clear`, so the port keeps the `Result` and ignores it, as
//!   `BRepOffset_MakeLoops` does.
//! * `myAnalyse.Generated(S)` answers a null shape for an unbound `S`; the
//!   port answers `None` and skips the branch (`:1841`, `:1857`).
//! * `TopExp_Explorer(myFaceComp, TopAbs_EDGE)` at `:1829` really does walk
//!   edges (`aF` is an edge despite its name); the inner `TopExp_Explorer(aF,
//!   TopAbs_EDGE)` therefore re-yields that same edge. Both loops are ported
//!   literally.
//! * `NCollection_Map::Add` answers false for a repeated key (`:1836`,
//!   `:1852`); the port keeps a [`HashSet<ShapeMapKey>`] for the same `View`.
//! * `TopoDS_Iterator(aE)` (`:1848`) is the plain child walk (`cumOri = true`),
//!   the port's [`children_keep_orientation`].

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Face, TopoShape};
use crate::topo_tools::{children_keep_orientation, explorer_shapes, shape_map_key, ShapeMapKey};

use super::BRepOffsetMakeOffset;

/// `SetFaces()` (`cxx:795-815`).
pub fn set_faces(this: &mut BRepOffsetMakeOffset) {
    for ii in 1..=this.my_original_faces.extent() {
        let Some(stored) = this.my_original_faces.value(ii).cloned() else {
            continue;
        };
        let Some(mut a_face) = Face::wrap(stored) else {
            continue;
        };
        if let Some(a_planface) = this
            .my_face_planface_map
            .get(&shape_map_key(&a_face.0))
            .cloned()
        {
            match Face::wrap(a_planface) {
                Some(planface) => a_face = planface,
                None => continue,
            }
        }

        this.my_faces.add(&a_face.0);
        //-------------
        // MAJ SD.
        //-------------
        this.my_init_offset_face.set_root(&a_face.0);
        let _ = this.my_init_offset_face.bind(&a_face.0, &a_face.0);
        this.my_image_offset.set_root(&a_face.0);
    }
}

/// `SetFacesWithOffset()` (`cxx:818-833`): every face that `IsPlanar` replaced
/// by a plane hands its `SetOffsetOnFace` value over to the plane.
///
/// `NCollection_DataMap::Iterator` order is unspecified and the loop body only
/// rebinds independent keys, so the insertion-order walk of
/// [`super::FaceOffsetMap`] is equivalent.
pub fn set_faces_with_offset(this: &mut BRepOffsetMakeOffset) {
    let entries: Vec<(ShapeMapKey, TopoShape)> = this
        .my_face_planface_map
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (key, a_planface) in entries {
        let a_face = match this.face_of_key(&key) {
            Some(f) => f,
            None => continue,
        };
        let Some(an_offset) = this.my_face_offset.get(&a_face.0) else {
            continue;
        };
        this.my_face_offset.un_bind(&a_face.0);
        this.my_face_offset.bind(&a_planface, an_offset);
    }
}

/// `BuildFaceComp()` (`cxx:1873-1890`): a compound of the working faces, with
/// the planar substitutes in place of their originals and the original
/// orientations kept.
pub fn build_face_comp(this: &mut BRepOffsetMakeOffset) {
    let a_builder = TopoBuilder::new();
    let mut a_face_comp = Compound::new();
    let Some(shape) = this.my_shape.clone() else {
        this.my_face_comp = a_face_comp;
        return;
    };
    for a_face in explorer_shapes(&shape, ShapeType::Face) {
        let an_or = a_face.orientation();
        let a_face = match this
            .my_face_planface_map
            .get(&shape_map_key(&a_face))
            .cloned()
        {
            Some(a_planface) => a_planface,
            None => a_face,
        };
        a_builder.add_compound(&mut a_face_comp, &a_face.oriented(an_or));
    }
    this.my_face_comp = a_face_comp;
}

/// `ReplaceRoots()` (`cxx:1825-1867`): replace the artificial faces and edges
/// in the `InitOffset` maps with the original ones.
pub fn replace_roots(this: &mut BRepOffsetMakeOffset) {
    let mut view: HashSet<ShapeMapKey> = HashSet::new();
    for a_f in explorer_shapes(&this.my_face_comp.0, ShapeType::Edge) {
        for a_e in explorer_shapes(&a_f, ShapeType::Edge) {
            if !view.insert(shape_map_key(&a_e)) {
                continue;
            }

            let Some(a_f_gen) = this.my_analyse.generated(&a_e) else {
                continue;
            };

            let _ = this.my_init_offset_face.replace_root(&a_f_gen, &a_e);

            for a_v in children_keep_orientation(&a_e) {
                if !view.insert(shape_map_key(&a_v)) {
                    continue;
                }

                let Some(a_e_gen) = this.my_analyse.generated(&a_v) else {
                    continue;
                };

                let _ = this.my_init_offset_edge.replace_root(&a_e_gen, &a_v);
            }
        }
    }
}
