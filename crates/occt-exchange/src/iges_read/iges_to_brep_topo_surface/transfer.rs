
use super::*;


use occt_topo::shape::TopoShape;

use super::super::iges_to_brep_curve_and_surface::{is_basic_surface, Session};

pub(in crate::iges_read) fn transfer_topo_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `IGESToBRep_TopoSurface::TransferTopoSurface` (`cxx:184-232`).
    // `TheULength = 1` is reset even when the shape is already cached (`cxx:191-202`).
    session.u_length = 1.0;
    if let Some(hit) = session.cache.get(&index) {
        return Ok(hit.clone());
    }
    let type_number = session.model.entity(index)?.type_number;
    let shape = if is_basic_surface(type_number) {
        super::super::iges_to_brep_basic_surface::transfer_topo_basic_surface(session, index)?
    } else if type_number == 144 {
        transfer_trimmed_surface(session, index)?
    } else if type_number == 120 {
        transfer_revolution(session, index)?
    } else if type_number == 122 {
        super::super::tabulated_cylinder::transfer_tabulated_cylinder(session, index)?
    } else if type_number == 118 {
        super::super::ruled_surface::transfer_ruled_surface(session, index)?
    } else if type_number == 108 {
        transfer_plane(session, index)?
    } else if type_number == 143 {
        super::super::bounded_surface::transfer_bounded_surface(session, index)?
    } else if type_number == 140 {
        super::super::offset_surface::transfer_offset_surface(session, index)?
    } else if type_number == 402 {
        // `cxx:245-249`: `res = TransferPerforate(st402_9)` (limit: planes only).
        transfer_perforate(session, index)?
    } else {
        return Err(format!(
            "UNPORTED topo surface {type_number}: IGESToBRep_TopoSurface.cxx:184"
        ));
    };
    session.cache.insert(index, shape.clone());
    Ok(shape)
}
