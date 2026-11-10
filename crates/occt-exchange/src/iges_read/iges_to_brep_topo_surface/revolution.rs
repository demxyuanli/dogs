use std::f64::consts::PI;
use std::sync::Arc;

use super::*;

use occt_core::gp::{
    GpAx1, GpDir, GpPnt, GpTrsf,
    GpVec,
};
use occt_core::precision::PCONFUSION;
use occt_core::toploc::TopLocLocation;
use occt_geom::{GeomSurfaceOfRevolution, Surface};

use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::{Edge, TopoShape};

use super::super::iges_data_iges_model::EntityData;
use super::super::basis_parameter_curve::rotate_span;
use super::super::iges_to_brep_basic_curve::{
    convert_location, entity_location, line_transformed_ends, CurveSpan,
};
use super::super::iges_to_brep_curve_and_surface::Session;
use super::super::iges_to_brep_topo_curve::transfer_topo_curve;

pub(super) fn transfer_revolution(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferSurfaceOfRevolution` (`IGESToBRep_TopoSurface.cxx:699-879`).
    let (axis, generatrix, start_iges, end_iges) = match &session.model.entity(index)?.data {
        EntityData::Revolution {
            axis,
            generatrix,
            start_angle,
            end_angle,
        } => (*axis, *generatrix, *start_angle, *end_angle),
        _ => return Err("IGES 120 has no revolution parameters".into()),
    };
    if !super::super::iges_to_brep_curve_and_surface::is_topo_curve(
        session.model.entity(generatrix)?.type_number,
    ) {
        return Err(
            "IGES 120 generatrix is not a topo curve: IGESToBRep_TopoSurface.cxx:718-724".into(),
        );
    }
    if session.model.entity(generatrix)?.type_number == 110 {
        if let EntityData::Line { start, end } = &session.model.entity(generatrix)?.data {
            let p1 = GpPnt::new(start[0], start[1], start[2]);
            let p2 = GpPnt::new(end[0], end[1], end[2]);
            session.u_length = p1.distance(&p2) * session.model.unit_factor;
        }
    }
    if session.model.entity(axis)?.type_number != 110 {
        return Err("IGES 120 axis is not a line: IGESToBRep_TopoSurface.cxx:735-739".into());
    }
    // `TransferTopoCurve` then `extractCurve3d` (`cxx:744`, `cxx:155-179`).
    // One edge builds `Geom_SurfaceOfRevolution` (`cxx:793-821`). Any other
    // count calls `BRepPrimAPI_MakeRevol` (`cxx:837-848`).
    let generated = transfer_topo_curve(session, generatrix)?;
    let span = match extract_curve_3d(&generated) {
        Some(span) => span,
        None => return revol_profile(session, index, &generated, axis),
    };
    // `igesAxis->TransformedStartPoint/EndPoint` then unit scale, not the
    // transferred line's parameter ends (`IGESToBRep_TopoSurface.cxx:758-768`).
    // `gp_Vec(pt2, pt1)` is `pt1 - pt2`.
    let (pt1, pt2) = line_transformed_ends(session.model, axis)?;
    let axis_dir = GpDir::from_vec(&GpVec::from_pnts(&pt2, &pt1))
        .map_err(|_| "IGES 120 axis is null: IGESToBRep_TopoSurface.cxx:768".to_string())?;
    let revol_axis = GpAx1::new(pt1, axis_dir);
    let start_angle = 2.0 * PI - end_iges;
    let end_angle = 2.0 * PI - start_iges;
    let mut delta = end_angle - start_angle;
    let full = delta > 2.0 * PI - PCONFUSION;
    if full {
        delta = 2.0 * PI;
    }
    let mut generatrix_curve = rotate_span(&span, &revol_axis, start_angle)?;
    let mut axis = revol_axis;
    // `Shape.Move` (`cxx:860-870`) stores the location. The mesh reads the
    // surface with no `TopLoc`, so the same `gp_Trsf` is applied to the
    // generatrix and axis. `gp_Dir::Transform` uses the linear part
    // (`Geom_SurfaceOfRevolution.cxx:419-424`). A negative determinant
    // reverses that direction so `S'(u) = T(S(u))` at the original UV.
    // `UReverse` (`cxx:425-428`) also remaps `U` to `2*pi-U`; `Shape.Move`
    // does not rewrite pcurves, so that remap is not applied.
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            let mut gen = generatrix_curve.clone_dyn();
            gen.transform(&tr);
            generatrix_curve = Arc::from(gen);
            axis = moved_revolution_axis(&axis, &tr)?;
        }
    }
    // `cxx:796-820`: the native path stores the plain `Geom_SurfaceOfRevolution`
    // on the face for the full *and* the partial angle. The `!IsFullAngle`
    // restriction goes through `BRepBuilderAPI_MakeFace::Init(aResultSurf, 0,
    // endAngle - startAngle, VF, VL, Confusion)`, i.e. `BRepLib_MakeFace::Init`
    // (`BRepLib_MakeFace.cxx:463-546`), which keeps the whole basis surface and
    // only wraps in `Geom_RectangularTrimmedSurface` for an *offset* basis
    // (`BRepLib_MakeFace.cxx:505-517`). `TransferFace` then calls `EmptyCopy()`
    // (`IGESToBRep_BRepEntity.cxx:506`), dropping the isoparametric wire, so the
    // shape keeps the untrimmed surface: the probe reads it back as
    // `surf=Geom_SurfaceOfRevolution u=[0,6.28319] uclosed=1`, not as a trim
    // with `u=[0,endAngle-startAngle]`.
    let _ = (delta, full);
    let surf: Arc<dyn Surface> = Arc::new(GeomSurfaceOfRevolution::new(generatrix_curve, axis));
    Ok(make_face(surf, &[]).0)
}

/// `BRepPrimAPI_MakeRevol` after `generatrix.Move(startLoc)` (`cxx:779-848`).
/// The entity location is `Shape.Move` (`cxx:860-870`), not baked into the
/// surfaces.
pub(super) fn revol_profile(
    session: &mut Session<'_>,
    index: usize,
    generated: &TopoShape,
    axis_entity: usize,
) -> Result<TopoShape, String> {
    let (pt1, pt2) = line_transformed_ends(session.model, axis_entity)?;
    let axis_dir = GpDir::from_vec(&GpVec::from_pnts(&pt2, &pt1))
        .map_err(|_| "IGES 120 axis is null: IGESToBRep_TopoSurface.cxx:768".to_string())?;
    let revol_axis = GpAx1::new(pt1, axis_dir);
    let (start_iges, end_iges) = match &session.model.entity(index)?.data {
        EntityData::Revolution {
            start_angle,
            end_angle,
            ..
        } => (*start_angle, *end_angle),
        _ => return Err("IGES 120 has no revolution parameters".into()),
    };
    let start_angle = 2.0 * PI - end_iges;
    let end_angle = 2.0 * PI - start_iges;
    let mut delta = end_angle - start_angle;
    if delta > 2.0 * PI - PCONFUSION {
        delta = 2.0 * PI;
    }
    let mut profile = generated.clone();
    let mut start = GpTrsf::identity();
    start
        .set_rotation_ax1(&revol_axis, start_angle)
        .map_err(|_| "IGES 120 start rotation: IGESToBRep_TopoSurface.cxx:779".to_string())?;
    profile.move_location(&TopLocLocation::composed(
        &start,
        TopLocLocation::identity(),
    ));
    let mut shape = occt_topo::brep_sweep::make_revol(&profile, &revol_axis, delta)?;
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(shape)
}

/// `extractCurve3d` (`IGESToBRep_TopoSurface.cxx:155-179`).
pub(in crate::iges_read) fn extract_curve_3d(shape: &TopoShape) -> Option<CurveSpan> {
    let mut edges = Vec::new();
    collect_edges(shape, &mut edges);
    if edges.len() != 1 {
        return None;
    }
    let edge = Edge(edges.remove(0));
    let curve = BRepTool::edge_curve(&edge)?;
    let (first, last) = BRepTool::edge_parameters(&edge);
    if !first.is_finite() || !last.is_finite() {
        return None;
    }
    Some(CurveSpan { curve, first, last })
}
