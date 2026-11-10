
use super::*;


use occt_topo::builder::TopoBuilder;
use occt_topo::shape::TopoShape;
use occt_topo::tgeometry::GeometryRegistry;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;
use super::super::iges_to_brep_topo_curve::transfer_curve_on_face;

pub(super) fn transfer_trimmed_surface(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `IGESToBRep_TopoSurface::TransferTrimmedSurface` (`cxx:1219-1376`).
    let (surface, outer, inners) = match &session.model.entity(index)?.data {
        EntityData::Trimmed {
            surface,
            outer,
            inners,
        } => (*surface, *outer, inners.clone()),
        _ => return Err("IGES 144 has no trimmed-surface parameters".into()),
    };
    let basis = transfer_topo_surface(session, surface)?;
    if !basis.is_face() {
        return Err(
            "IGES 144 basis surface is not a face: IGESToBRep_TopoSurface.cxx:1254-1291".into(),
        );
    }
    let surf = GeometryRegistry::global()
        .face_surface(&basis)
        .ok_or_else(|| "IGES 144 basis face has no surface".to_string())?;
    let (trans, u_fact) = uv_map(session, surface, &basis)?;
    // `HasOuterContour` is a non-null outer curve (`IGESGeom_TrimmedSurface.cxx:53-55`).
    // `EmptyCopy` then `TransferCurveOnFace` (`cxx:1304-1319`).
    let mut face = make_face(surf, &[]);
    let mut added = false;
    if let Some(outer) = outer {
        let wire = transfer_curve_on_face(session, outer, &face.0, &trans, u_fact, true)?;
        TopoBuilder::new().add_wire(&mut face, &wire);
        added = true;
    }
    for inner in inners.into_iter().flatten() {
        let wire = transfer_curve_on_face(session, inner, &face.0, &trans, u_fact, true)?;
        TopoBuilder::new().add_wire(&mut face, &wire);
        added = true;
    }
    if !added {
        return Ok(basis);
    }
    apply_trimmed_transform(session, index, &mut face.0)?;
    Ok(face.0)
}
