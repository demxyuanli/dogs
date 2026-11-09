//! `BRepOffset_MakeOffset::MakeOffsetFaces` (`BRepOffset_MakeOffset.cxx:1202-1275`)
//! and `UpdateFaceOffset` (`:2342-2394`).
//!
//! `MakeOffsetFaces` is the bridge between `BRepOffset_Analyse` and
//! `BRepOffset_Offset`: it walks the faces in `BRepLib::SortFaces` order (so
//! the planes and the other analytic surfaces come before the free-form ones),
//! builds one `BRepOffset_Offset` per face - with the per-face
//! `SetOffsetOnFace` value when there is one, `myOffset` otherwise - and fills
//! `ShapeTgt` with the images the next face needs for its shared boundary: the
//! offset of every tangential edge of the face, plus the two endpoint vertex
//! images when those vertices are not shared with a non-tangential neighbour
//! (`JLE.Extent() == JLA.Extent()`, `cxx:1247` / `:1256`). The faces
//! `BRepOffset_Analyse` created on its own (`NewFaces`, `cxx:1265`) get a
//! degenerate offset at `0.0`.
//!
//! `UpdateFaceOffset` propagates a per-face offset over the tangential cloud
//! around that face: `BRepOffset_Analyse::AddFaces` collects the connected
//! faces and every face in that cloud without an explicit offset takes the
//! value of the face the walk started from. `GeomAbs_Arc` also lets the
//! `ChFiDS_Convex` / `ChFiDS_Concave` faces in (`cxx:2369-2374`).
//!
//! Port guards (no OCCT branch changed):
//!
//! * `BRepOffset_Offset::Init` raises `Standard_ConstructionError` /
//!   `Standard_TypeMismatch` / `Standard_NoSuchObject` in several arms
//!   (`BRepOffset_Offset.cxx:466-1094`). OCCT lets that exception leave
//!   `MakeOffsetFaces` and unwind the whole offset; the port records
//!   [`Error::UnknownError`] and returns.
//! * `TopoDS::Edge(aLocalShape)` (`cxx:1235`) and `TopoDS::Face(...)`
//!   (`cxx:1225`, `:1269`) force the type; the port wraps through
//!   [`Edge::wrap`] / [`Face::wrap`] and skips a mismatch, which
//!   `Generated` cannot produce for those inputs.
//! * `OF.Generated(Cur)` on an unbound `Cur` answers a null shape in OCCT and
//!   the following `TopExp::Vertices` then yields two null vertices, whose
//!   `ShapeTgt.Bind(V1, OV1)` binds a null; the port skips the pair instead
//!   (`None`).
//! * `Message_ProgressScope` plus `myError = BRepOffset_UserBreak` on
//!   `!More()` (`cxx:1218-1222`): no progress machinery in the port, same as
//!   the rest of `BRepOffset`.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_offset::{BRepOffsetOffset, Error, MapSf, ShapeMap};
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{Edge, Face};
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key};

use super::BRepOffsetMakeOffset;

/// `MakeOffsetFaces(theMapSF, theRange)` (`cxx:1202-1275`).
///
/// `theMapSF` is `myMapSF` of `BuildOffsetByArc` / `BuildOffsetByInter`; the
/// `theRange` argument and the `Message_ProgressScope` are not ported.
pub fn make_offset_faces(this: &mut BRepOffsetMakeOffset, the_map_sf: &mut MapSf) {
    //
    let offset_outside = this.my_offset > 0.0;
    //
    let a_lf = crate::brep_lib_sort_faces::sort_faces(&this.my_face_comp.0);
    let mut shape_tgt: ShapeMap = HashMap::new();

    for a_f in a_lf {
        let Some(a_f) = Face::wrap(a_f) else {
            continue;
        };
        let a_cur_offset = this
            .my_face_offset
            .get(&a_f.0)
            .unwrap_or(this.my_offset);

        let mut of = BRepOffsetOffset::new();
        if of
            .init_face(
                &a_f,
                a_cur_offset,
                &shape_tgt,
                offset_outside,
                this.my_join,
            )
            .is_err()
        {
            // OCCT lets the exception escape this method; see the module notes.
            this.my_error = Error::UnknownError;
            return;
        }

        for cur in this
            .my_analyse
            .edges_of_face(&a_f, TypeOfConcavity::Tangential)
        {
            if shape_tgt.contains_key(&shape_map_key(&cur.0))
                || this.my_analyse.has_generated(&cur.0)
            {
                continue;
            }
            let Some(a_local_shape) = of.generated(&cur.0) else {
                continue;
            };
            let Some(ote) = Edge::wrap(a_local_shape.clone()) else {
                continue;
            };
            shape_tgt.insert(shape_map_key(&cur.0), a_local_shape);

            let (v1, v2) = edge_vertices(&cur);
            let (ov1, ov2) = edge_vertices(&ote);
            if let Some(v1) = v1.as_ref() {
                if !shape_tgt.contains_key(&shape_map_key(&v1.0)) {
                    let le = this
                        .my_analyse
                        .edges_of_vertex(v1, TypeOfConcavity::Tangential);
                    let la = this.my_analyse.ancestors(&v1.0);
                    if le.len() == la.len() {
                        if let Some(ov1) = ov1.as_ref() {
                            shape_tgt.insert(shape_map_key(&v1.0), ov1.0.clone());
                        }
                    }
                }
            }
            if let Some(v2) = v2.as_ref() {
                if !shape_tgt.contains_key(&shape_map_key(&v2.0)) {
                    let le = this
                        .my_analyse
                        .edges_of_vertex(v2, TypeOfConcavity::Tangential);
                    let la = this.my_analyse.ancestors(&v2.0);
                    if le.len() == la.len() {
                        if let Some(ov2) = ov2.as_ref() {
                            shape_tgt.insert(shape_map_key(&v2.0), ov2.0.clone());
                        }
                    }
                }
            }
        }
        the_map_sf.insert(shape_map_key(&a_f.0), of);
    }
    //
    let a_new_faces = this.my_analyse.new_faces().to_vec();
    for a_f in a_new_faces {
        let Some(a_f) = Face::wrap(a_f) else {
            continue;
        };
        let mut of = BRepOffsetOffset::new();
        if of
            .init_face(&a_f, 0.0, &shape_tgt, offset_outside, this.my_join)
            .is_err()
        {
            this.my_error = Error::UnknownError;
            return;
        }
        the_map_sf.insert(shape_map_key(&a_f.0), of);
    }
}

/// `UpdateFaceOffset()` (`cxx:2342-2394`).
pub fn update_face_offset(this: &mut BRepOffsetMakeOffset) {
    // `NCollection_Map` of the faces already handled in this call.
    let mut m: Vec<crate::topo_tools::ShapeMapKey> = Vec::new();
    // `CopiedMap.Assign(myFaceOffset)`: the walk mutates `myFaceOffset`.
    let copied_map: Vec<(crate::shape::TopoShape, f64)> = this
        .my_face_offset
        .iter()
        .map(|(s, v)| (s.clone(), v))
        .collect();

    let mut rt = TypeOfConcavity::Convex;
    if this.my_offset < 0.0 {
        rt = TypeOfConcavity::Concave;
    }

    let builder = crate::builder::TopoBuilder::new();
    for (f, cur_offset) in copied_map {
        if m.contains(&shape_map_key(&f)) {
            continue;
        }
        m.push(shape_map_key(&f));
        let f = Face(f);

        let mut co = crate::shape::Compound::new();
        let mut dummy: Vec<crate::topo_tools::ShapeMapKey> = Vec::new();
        builder.add_compound(&mut co, &f.0);
        if this.my_join == occt_core::kernel::geomabs::JoinType::Arc {
            this.my_analyse.add_faces_two(
                &f,
                &mut co,
                &mut dummy,
                TypeOfConcavity::Tangential,
                rt,
            );
        } else {
            this.my_analyse
                .add_faces(&f, &mut co, &mut dummy, TypeOfConcavity::Tangential);
        }

        for ff in explorer_shapes(&co.0, ShapeType::Face) {
            let ff_key = shape_map_key(&ff);
            if m.contains(&ff_key) {
                continue;
            }
            m.push(ff_key);
            this.my_face_offset.un_bind(&ff);
            this.my_face_offset.bind(&ff, cur_offset);
        }
    }
}
