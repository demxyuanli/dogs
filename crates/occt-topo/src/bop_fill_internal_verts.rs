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
    let mut a_vvfi: Vec<VertexFaceInternal> = Vec::new();
    for i in 0..n {
        let Some(a_si) = f.ds().shape_info(i) else { continue };
        if a_si.shape_type() != ShapeType::Face {
            continue;
        }
        let a_f = a_si.shape();
        let Some(p_lf_im) = f.history().image(a_f) else { continue };
        if p_lf_im.is_empty() {
            continue;
        }
        let a_liav = alone_vertices(f.ds(), i);
        for v_idx in a_liav {
            let Some(v_shape) = f.ds().shape(v_idx).cloned() else { continue };
            let mut a_v = Vertex(v_shape);
            a_v.0.set_orientation(Orientation::Internal);
            for a_f_im in p_lf_im {
                a_vvfi.push(VertexFaceInternal::from_parts(
                    a_v.clone(),
                    Face(a_f_im.clone()),
                    fuzzy,
                ));
            }
        }
    }
    let mut ctx = IntToolsContext::new();
    perform_vfi_vector(&mut a_vvfi, &mut ctx);
    let bld = TopoBuilder::new();
    for a_vfi in &mut a_vvfi {
        if !a_vfi.is_internal() {
            continue;
        }
        // `BRep_Builder().Add(aF, aV)` — the image TShape is shared, so the
        // add is visible on every view of that face (`_2.cxx:1003-1006`).
        let mut face = a_vfi.face().0.clone();
        let vertex = a_vfi.vertex().0.clone();
        bld.add(&mut face, &vertex);
    }
    Ok(())
}
