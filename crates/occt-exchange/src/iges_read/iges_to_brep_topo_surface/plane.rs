use std::sync::Arc;

use super::*;

use occt_core::gp::{
    GpPln, GpPnt, GpTrsf, GpTrsf2d,
};
use occt_core::precision::CONFUSION;
use occt_geom::GeomPlane;

use occt_topo::brep_gprop_full::surface_properties;
use occt_topo::builder::TopoBuilder;
use occt_topo::shape::{Edge, TopoShape, Wire};
use occt_topo::shhealing::fix_connected_all;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_basic_curve::eps_geom;
use super::super::iges_to_brep_curve_and_surface::{is_topo_curve, Session};
use super::super::iges_to_brep_topo_curve::{transfer_curve_on_face, transfer_topo_curve};

pub(super) fn transfer_plane(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferPlane` (`cxx:1493-1516`) calls `TransferPlaneParts` with
    // `first = true` and then applies the matrix that function computed
    // (`cxx:1511-1515`).
    let (shape, _pln, trsf) = transfer_plane_parts(session, index, true)?;
    let shape = shape.ok_or_else(|| "IGES 108 plane transfer returned no shape".to_string())?;
    Ok(apply_plane_trsf(shape, &trsf))
}

/// `TransferPlaneParts` (`IGESToBRep_TopoSurface.cxx:1608-1785`). `first = true`
/// returns the face together with its external contour, `first = false` the
/// contour alone, which is what `TransferPerforate` adds as a hole
/// (`cxx:1567`). `None` is OCCT's null `TopoDS_Shape`.
pub(super) fn transfer_plane_parts(
    session: &mut Session<'_>,
    index: usize,
    first: bool,
) -> Result<(Option<TopoShape>, GpPln, GpTrsf), String> {
    // Geom equation is `ax+by+cz+d = 0` with `d = -D` (`cxx:1624-1628`).
    let ent = session.model.entity(index)?;
    let form = ent.form_number;
    let (a, b, c, d, bounding) = match &ent.data {
        EntityData::Plane {
            a,
            b,
            c,
            d,
            bounding,
        } => (*a, *b, *c, *d, *bounding),
        _ => return Err("IGES 108 has no plane parameters".into()),
    };
    let pln_raw = GpPln::from_coefficients(a, b, c, -d)
        .map_err(|_| "IGES 108 plane normal is null: gp_Pln.cxx:86".to_string())?;
    let mut pln = pln_raw.clone();
    let unit = session.model.unit_factor;
    if unit != 1.0 {
        // `geomPln->Scale(origin, GetUnitFactor)` (`cxx:1637`).
        let loc = pln.location();
        pln.set_location(&GpPnt::new(loc.x() * unit, loc.y() * unit, loc.z() * unit));
    }
    // `cxx:1647-1658`: the matrix is computed here and applied by the caller.
    let trsf = plane_trsf(session, index)?;
    // `cxx:1631-1644`: only `first` creates the face. Otherwise `plane` stays
    // null and only `wire` can be returned.
    let mut face = if first {
        Some(make_face(Arc::new(GeomPlane::new(pln.clone())), &[]))
    } else {
        None
    };
    // `HasBoundingCurve` is a non-null curve pointer (`IGESGeom_Plane.cxx:70-72`).
    // A null pointer leaves the plane untrimmed.
    let wire = match bounding {
        Some(curve) => plane_bounding_wire(
            session,
            face.as_mut(),
            curve,
            &pln_raw,
            // `IGESGeom_Plane::HasBoundingCurveHole` (`IGESGeom_Plane.cxx:45,75`):
            // form number `-1` is a hole, `+1` an external bound.
            form == -1,
            unit,
        )?,
        None => None,
    };
    // `cxx:1756-1762`: only `first` attaches the contour to the face.
    if let (Some(f), Some(w)) = (face.as_mut(), wire.as_ref()) {
        TopoBuilder::new().add_wire(f, w);
    }
    let res = if first {
        face.map(|f| f.0)
    } else {
        wire.map(|w| w.0)
    };
    Ok((res, pln, trsf))
}

/// The contour half of `TransferPlaneParts` (`IGESToBRep_TopoSurface.cxx:1660-1775`).
/// Returns the transferred contour, or `None` where OCCT leaves `wire` null.
pub(super) fn plane_bounding_wire(
    session: &mut Session<'_>,
    face: Option<&mut occt_topo::shape::Face>,
    curve: usize,
    pln_raw: &GpPln,
    hole: bool,
    unit: f64,
) -> Result<Option<Wire>, String> {
    if session.model.entity(curve).is_err() {
        // `cxx:1665-1671`: null bounding curve, keep the infinite plane.
        return Ok(None);
    }
    let curve_type = session.model.entity(curve)?.type_number;
    if !is_topo_curve(curve_type) {
        // `cxx:1765-1773`: not a topo curve. The warning leaves the plane untrimmed.
        return Ok(None);
    }
    if curve_type == 142 {
        // `cxx:1677-1693`. The local `gp_Trsf` is identity, so the `gp_Trsf2d`
        // conversion (`gp_Trsf2d.hxx:232-241`) is identity. `TheULength` is the
        // value reset at `TransferTopoSurface` (`cxx:191`). `cxx:1682` discards
        // the transferred shape, so `wire` stays null and only the face is
        // trimmed.
        if let Some(f) = face {
            let trans = GpTrsf2d::identity();
            if let Ok(w) =
                transfer_curve_on_face(session, curve, &f.0, &trans, session.u_length, true)
            {
                TopoBuilder::new().add_wire(f, &w);
            }
        }
        return Ok(None);
    }
    let shape = match transfer_topo_curve(session, curve) {
        Ok(shape) => shape,
        // `cxx:1694-1726`: a failed transfer has no shape type and takes the
        // `default` arm (the plane for `first`, the null wire otherwise).
        Err(_) => return Ok(None),
    };
    let mut wire = if shape.is_edge() {
        TopoBuilder::new().make_wire(&[Edge(shape)])
    } else if shape.is_wire() {
        Wire(shape)
    } else {
        // `cxx:1710-1725`, anything else warns and returns the plane (`first`) or
        // the null wire.
        return Ok(None);
    };
    // Naked `ShapeFix_Wire` (`cxx:1729-1732`). Default precision is
    // `Precision::Confusion` (`ShapeFix_Root.cxx:26`). `FixConnected(GetMaxTol)`
    // passes that value as the check precision (`ShapeFix_Wire.cxx:557-569`).
    let eps = eps_geom(session.model.resolution) * unit;
    let _ = fix_connected_all(&mut wire, CONFUSION, eps.max(1.0));
    // `BRepLib_MakeFace(pln, wire, Inside=false)` (`cxx:1733`,
    // `BRepLib_MakeFace.cxx:262-270`) uses the unscaled plane and only feeds
    // `BRepGProp::SurfaceProperties`.
    //
    // `cxx:1734-1738` returns the null shape when `BRepLib_MakeFace::IsDone()`
    // is false. That arm is dead for this call, so the port runs straight into
    // the properties stage, which is OCCT's own behaviour: the constructor is
    // `Init(Pl, false, Confusion()); Add(W);` (`BRepLib_MakeFace.cxx:262-270`),
    // `Init(S, false, tol)` only sets `myError = BRepLib_FaceDone` and builds
    // the natural face without ever calling `Done()`
    // (`BRepLib_MakeFace.cxx:366-388`), and `Add(const TopoDS_Wire&)` ends with
    // `Done()` unconditionally (`BRepLib_MakeFace.cxx:872-877`). `IsDone()` is
    // `BRepBuilderAPI_Command::myDone` (`BRepBuilderAPI_Command.hxx:42` :58,
    // set by `Done()` at `:52`), so it is true by the time the test is reached.
    let probe = TopoBuilder::new().make_face(Arc::new(GeomPlane::new(pln_raw.clone())), &[wire.clone()]);
    // `BRepGProp::SurfaceProperties` (`cxx:1741-1752`). A missing mass leaves
    // the wire as transferred; the sign test is the only reverse.
    if let Ok((props, _)) = surface_properties(&probe.0) {
        if props.mass() < 0.0 {
            if !hole {
                wire.0.reverse();
            }
        } else if hole {
            wire.0.reverse();
        }
    }
    Ok(Some(wire))
}
