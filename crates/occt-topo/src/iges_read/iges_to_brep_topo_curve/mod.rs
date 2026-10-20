//! `IGESToBRep_TopoCurve` for entities 110, 100, 126, and 102, and
//! `TransferCurveOnFace` (`IGESToBRep_TopoCurve.cxx:512-578`).
//!
//! The 142 UV arm is `IGESControl_IGESBoundary::Transfer`
//! (`IGESControl_IGESBoundary.cxx:146-536`). `read.stdsameparameter.mode`
//! is 0 (`XSAlgo.cxx:46`), so the copied pcurve goes through `FixReversed2d`
//! then `FixSameParameter(edge)`. A failed `ConnectNextWire` rebuilds both
//! wires with `ShapeFix_Wire::Perform` (`cxx:441-508`). `Check`
//! (`cxx:60-109`) then keeps that wire, or the 3D or 2D wire when closure fails.

mod basic_curve;
mod boundary;
mod composite_curve;
mod curve_on_face;
mod pcurve_edge;
mod pcurve_transfer;
mod wire_data_boundary;

pub(in crate::iges_read) use self::boundary::{
    boundary_after_check, boundary_both, boundary_max_tol, BoundaryPiece,
};
pub(in crate::iges_read) use self::curve_on_face::transfer_curve_on_face;
pub(in crate::iges_read) use self::pcurve_transfer::{
    curve_sequence, edges_from_2d, edges_from_2d_basic, edges_of, transfer_pcurve,
};
pub(in crate::iges_read) use self::wire_data_boundary::WireDataBoundary;

use crate::shape::TopoShape;

use super::iges_to_brep_basic_curve::basic_curve_3d;
use super::iges_to_brep_curve_and_surface::{is_basic_curve, Session};

use basic_curve::{transfer_point, transfer_topo_basic_curve};
use composite_curve::transfer_composite;

pub(in crate::iges_read) fn transfer_topo_curve(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `IGESToBRep_TopoCurve::TransferTopoCurve` (`cxx:859-905`).
    if let Some(hit) = session.cache.get(&index) {
        return Ok(hit.clone());
    }
    let ent = session.model.entity(index)?;
    let type_number = ent.type_number;
    let shape = match type_number {
        110 | 100 | 104 | 106 | 112 | 126 => {
            // `TransferTopoBasicCurve` sets `ModeTransfer(false)` (`cxx:967-968`),
            // so a directory 124 is applied inside `TransferBasicCurve`.
            let span = basic_curve_3d(session.model, index)?;
            transfer_topo_basic_curve(&span)?
        }
        102 => transfer_composite(session, index)?,
        142 => super::curve_on_surface::transfer_curve_on_surface(session, index)?,
        141 => super::bounded_surface::transfer_boundary(session, index)?,
        116 => transfer_point(session, index)?,
        130 => super::offset_curve::transfer_offset_curve(session, index)?,
        _ if is_basic_curve(type_number) => {
            return Err(format!(
                "UNPORTED basic curve {type_number}: IGESToBRep_TopoCurve.cxx:871-873"
            ));
        }
        _ => {
            return Err(format!(
                "UNPORTED topo curve {type_number}: IGESToBRep_TopoCurve.cxx:859"
            ));
        }
    };
    session.cache.insert(index, shape.clone());
    Ok(shape)
}

/// `read.iges.bspline.continuity` (`IGESData.cxx:66`) -- the
/// `IGESToBRep_CurveAndSurface::myContinuity` that `GetContinuity()` returns
/// (`IGESToBRep_CurveAndSurface.cxx:58-64`).
const CONTINUITY: i32 = 1;
