//! Free entity 142. `IGESToBRep_TopoCurve::TransferCurveOnSurface`
//! (`IGESToBRep_TopoCurve.cxx:405-508`).
//!
//! The face is `ParamSurface`, then `EmptyCopy`, then `TransferCurveOnFace`
//! with `isCurveOnSurf` true. The returned shape is the wire (`cxx:577`).
//! A shell with a face count other than 1, or any other shape, returns the
//! 3D curve (`cxx:452-500`). `HasTransf` uses `ConvertLocation` with the
//! default unit 1 (`IGESData_ToolLocation.hxx:151`).

use occt_core::toploc::TopLocLocation;

use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_curve_and_surface::{is_topo_curve, is_topo_surface, Session};
use super::iges_to_brep_topo_curve::{transfer_curve_on_face, transfer_topo_curve};
use super::iges_to_brep_topo_surface::{transfer_topo_surface, uv_map};

pub(in crate::iges_read) fn transfer_curve_on_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    let (surface, curve_3d) = match &session.model.entity(index)?.data {
        EntityData::CurveOnSurface {
            surface, curve_3d, ..
        } => (*surface, *curve_3d),
        _ => return Err("IGES 142 has no curve-on-surface parameters".into()),
    };
    let Some(surface) = surface else {
        return Err("XSTEP_131: IGESToBRep_TopoCurve.cxx:420-425".into());
    };
    if !is_topo_surface(session.model.entity(surface)?.type_number) {
        return Err("XSTEP_131: IGESToBRep_TopoCurve.cxx:420-425".into());
    }
    let myshape = transfer_topo_surface(session, surface)?;
    let face = if myshape.is_face() {
        myshape
    } else if myshape.is_shell() {
        let kids = myshape.tshape.read().expect("lock").children.clone();
        let faces: Vec<TopoShape> = kids.iter().filter(|child| child.is_face()).cloned().collect();
        if kids.len() == 1 && faces.len() == 1 {
            faces.into_iter().next().expect("one face")
        } else {
            return curve_3d_fallback(session, index, curve_3d);
        }
    } else {
        return curve_3d_fallback(session, index, curve_3d);
    };
    let surf = GeometryRegistry::global()
        .face_surface(&face)
        .ok_or_else(|| "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1851".to_string())?;
    let (trans, u_fact) = uv_map(session, surface, &face)?;
    // `face.EmptyCopy` (`cxx:505`) keeps the surface and the location.
    let mut blank = TopoBuilder::new().make_face(surf, &[]);
    if !face.location().is_identity() {
        blank.0.move_location(face.location());
    }
    // `isCurveOnSurf` is true, so `Check` does not test closure (`cxx:506`, `:545`).
    let mut wire = transfer_curve_on_face(session, index, &blank.0, &trans, u_fact, false)?;
    move_wire_unit1(session, index, &mut wire.0)?;
    Ok(wire.0)
}

fn curve_3d_fallback(
    session: &mut Session<'_>,
    index: usize,
    curve_3d: Option<usize>,
) -> Result<TopoShape, String> {
    // `cxx:452-474` and the default arm `cxx:478-500`. `IGES_1062` is a warning
    // and is not stored. A missing 3D curve is `IGES_1061`.
    let Some(curve) = curve_3d else {
        return Err(format!(
            "IGES_1061: IGESToBRep_TopoCurve.cxx:469-473 entity {index}"
        ));
    };
    if !is_topo_curve(session.model.entity(curve)?.type_number) {
        return Err("IGES 142 3D curve is not a topo curve: IGESToBRep_TopoCurve.cxx:455".into());
    }
    transfer_topo_curve(session, curve)
}

fn move_wire_unit1(session: &Session<'_>, index: usize, shape: &mut TopoShape) -> Result<(), String> {
    // `cxx:558-571`. The 3-argument `ConvertLocation` uses unit 1.
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, 1.0) {
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(())
}
