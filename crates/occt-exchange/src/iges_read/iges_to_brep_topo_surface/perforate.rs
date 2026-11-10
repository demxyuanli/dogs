
use super::*;

use occt_core::toploc::TopLocLocation;

use occt_topo::builder::TopoBuilder;
use occt_topo::shape::{Face, TopoShape, Wire};

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;

/// `TransferPerforate` (`IGESToBRep_TopoSurface.cxx:1526-1604`), the 402 form 9
/// arm of `TransferTopoSurface` (`cxx:245-249`, "limit: Planes only"): the
/// parent plane gives the face and its external contour, every child plane a
/// hole.
pub(super) fn transfer_perforate(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    let (parent, children) = match &session.model.entity(index)?.data {
        EntityData::SingleParent { parent, children } => (*parent, children.clone()),
        _ => return Err(format!("IGES {index} is not a 402 form 9 SingleParent")),
    };
    // `DeclareAndCast(IGESGeom_Plane, p0, st->SingleParent()); if (p0.IsNull())`
    // (`cxx:1543-1551`): the SingleParent does not describe a holed face.
    let Some(parent) = parent else {
        return Err("XSTEP_206: IGESToBRep_TopoSurface.cxx:1543-1551".into());
    };
    // `res = TransferPlaneParts(p0, pln, trsf, true)` (`cxx:1553`), which starts
    // with the face and its external contour.
    let (shape, _pln, trsf) = transfer_plane_parts(session, parent, true)?;
    let shape = shape.ok_or_else(|| "IGES 1005: IGESToBRep_TopoSurface.cxx:1526-1538".to_string())?;
    let mut face = Face(shape);
    // `int nb = st->NbChildren(); for (int i = 1; i <= nb; i++)` (`cxx:1555-1598`).
    for child in children {
        let Some(child) = child else {
            continue;
        };
        // `DeclareAndCast(IGESGeom_Plane, pi, st->Child(i)); if (pi.IsNull())`
        // (`cxx:1557-1566`): `IGES_1285` "A child is not a Plane, skipped".
        if !matches!(session.model.entity(child), Ok(e) if e.type_number == 108) {
            continue;
        }
        let (wire, _pli, trsi) = transfer_plane_parts(session, child, false)?;
        // `cxx:1582-1588` compares the two `gp_Ax3` with `IsCoplanar` and only
        // sends `IGES_1295`; the wire is added either way and the port has no
        // messenger.
        // `if (wire.ShapeType() != TopAbs_WIRE)` (`cxx:1569-1581`): `IGES_1156`
        // "A hole could not be transferred, skipped".
        let Some(wire) = wire.filter(|w| w.is_wire()) else {
            continue;
        };
        // `wire.Location(locw, false)` (`cxx:1589-1593`): the child's own local
        // transformation is composed onto the hole contour.
        let mut wire = Wire(wire);
        wire.0
            .move_location(&TopLocLocation::composed(&trsi, TopLocLocation::identity()));
        TopoBuilder::new().add_wire(&mut face, &wire);
    }
    // `res.Location(loc, false)` (`cxx:1597-1603`).
    Ok(apply_plane_trsf(face.0, &trsf))
}
