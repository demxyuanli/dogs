//! `BOPAlgo_Builder::FillInternalVertices`.
//!
//! Source: `BOPAlgo_Builder_2.cxx:929-1008`. For every source FACE that has
//! images, the DS alone vertices are oriented INTERNAL and paired with each
//! image. `BOPAlgo_VFI::Perform` (`_2.cxx:188-200`) classifies the pair with
//! `IntTools_Context::ComputeVF`; flag 0 means the vertex lies strictly
//! inside the face and `BRep_Builder::Add` attaches it.
//!
//! This is the third step of `FillImagesFaces` (`_2.cxx:215-229`). The
//! previous port in [`crate::bop_build_common::fill_internal_vertices`] used
//! the same `ComputeVF` test but was not wired into `FillImagesFaces`.

use crate::abs::{Orientation, ShapeType};
use crate::bop_build_faces::BopBuilderLike;
use crate::bop_occt_util::{alone_vertices, nb_source};
use crate::bop_pair_sd::{perform_vfi_vector, VertexFaceInternal};
use crate::builder::TopoBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Face, Vertex};

/// `BOPAlgo_Builder::FillInternalVertices` (`_2.cxx:929-1008`).
pub fn fill_internal_vertices_occt<B: BopBuilderLike>(f: &mut B) -> Result<(), String> {
    let n = nb_source(f.ds());
    let fuzzy = f.fuzzy_value();
    let mut vvfi: Vec<VertexFaceInternal> = Vec::new();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Face {
            continue;
        }
        let f_val = si.shape();
        let Some(p_lf_im) = f.history().image(f_val) else { continue };
        if p_lf_im.is_empty() {
            continue;
        }
        let liav = alone_vertices(f.ds(), i);
        for v_idx in liav {
            let Some(v_shape) = f.ds().shape(v_idx).cloned() else { continue };
            let mut v = Vertex(v_shape);
            v.0.set_orientation(Orientation::Internal);
            for f_im in p_lf_im {
                // `BOPAlgo_Builder_2.cxx:876-878`: `SetVertex` / `SetFace` /
                // `SetFuzzyValue` on the freshly constructed pair.
                let mut vfi = VertexFaceInternal::new();
                vfi.set_vertex(v.clone());
                vfi.set_face(Face(f_im.clone()));
                vfi.set_fuzzy_value(fuzzy);
                vvfi.push(vfi);
            }
        }
    }
    let mut ctx = IntToolsContext::new();
    perform_vfi_vector(&mut vvfi, &mut ctx);
    let bld = TopoBuilder::new();
    for vfi in &mut vvfi {
        if !vfi.is_internal() {
            continue;
        }
        // `_2.cxx:904-905` binds `TopoDS_Vertex& aV = aVFI.Vertex()` and
        // `TopoDS_Face& aF = aVFI.Face()`, then `BRep_Builder().Add(aF, aV)`
        // (`_2.cxx:1003-1006`). Rust cannot hold both `&mut` at once, so the
        // vertex view is taken first; the image TShape is shared, so adding to
        // the stored face is visible on every view of it.
        let vertex = vfi.vertex_mut().0.clone();
        let face = vfi.face_mut();
        bld.add(&mut face.0, &vertex);
    }
    Ok(())
}
