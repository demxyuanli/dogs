//! `IGESToBRep_TopoSurface::TransferTopoSurface` for entities 144, 108, and 120
//! (`IGESToBRep_TopoSurface.cxx:184-232`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::{clib2d, slib};
use occt_core::gp::{
    GpAx1, GpAx22d, GpAx2d, GpCirc2d, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d, GpTrsf, GpTrsf2d,
    GpVec, GpVec2d, TrsfForm,
};
use occt_core::precision::{Precision, ANGULAR, CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::{GeomPlane, GeomRectangularTrimmedSurface, GeomSurfaceOfRevolution, Surface};
use occt_geom2d::Curve2d;

use crate::abs::ShapeType;
use crate::brep_gprop_full::surface_properties;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Wire};
use crate::shhealing::fix_connected_all;
use crate::tgeometry::GeometryRegistry;

use super::iges_data_iges_model::EntityData;
use super::basis_parameter_curve::rotate_span;
use super::iges_to_brep_basic_curve::{
    convert_location, entity_location, eps_geom, line_transformed_ends, CurveSpan,
};
use super::iges_to_brep_curve_and_surface::{is_basic_surface, is_topo_curve, Session};
use super::iges_to_brep_topo_curve::{transfer_curve_on_face, transfer_topo_curve};

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
        super::iges_to_brep_basic_surface::transfer_topo_basic_surface(session, index)?
    } else if type_number == 144 {
        transfer_trimmed_surface(session, index)?
    } else if type_number == 120 {
        transfer_revolution(session, index)?
    } else if type_number == 122 {
        super::tabulated_cylinder::transfer_tabulated_cylinder(session, index)?
    } else if type_number == 118 {
        super::ruled_surface::transfer_ruled_surface(session, index)?
    } else if type_number == 108 {
        transfer_plane(session, index)?
    } else if type_number == 143 {
        super::bounded_surface::transfer_bounded_surface(session, index)?
    } else if type_number == 140 {
        super::offset_surface::transfer_offset_surface(session, index)?
    } else {
        return Err(format!(
            "UNPORTED topo surface {type_number}: IGESToBRep_TopoSurface.cxx:184"
        ));
    };
    session.cache.insert(index, shape.clone());
    Ok(shape)
}

fn transfer_trimmed_surface(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
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

fn transfer_revolution(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
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
    if !super::iges_to_brep_curve_and_surface::is_topo_curve(
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
    let v0 = generatrix_curve.first_parameter();
    let v1 = generatrix_curve.last_parameter();
    let revol = Arc::new(GeomSurfaceOfRevolution::new(generatrix_curve, axis));
    let surf: Arc<dyn Surface> = if full {
        revol
    } else {
        Arc::new(GeomRectangularTrimmedSurface::uv(
            revol,
            0.0,
            delta,
            v0,
            v1,
        ))
    };
    Ok(make_face(surf, &[]).0)
}

/// `BRepPrimAPI_MakeRevol` after `generatrix.Move(startLoc)` (`cxx:779-848`).
/// The entity location is `Shape.Move` (`cxx:860-870`), not baked into the
/// surfaces.
fn revol_profile(
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
    let mut shape = crate::brep_sweep::make_revol(&profile, &revol_axis, delta)?;
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

fn collect_edges(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    if shape.shape_type() == ShapeType::Edge {
        out.push(shape.clone());
        return;
    }
    let children = shape.tshape.read().expect("lock").children.clone();
    for child in children {
        collect_edges(&child, out);
    }
}

fn moved_revolution_axis(axis: &GpAx1, tr: &GpTrsf) -> Result<GpAx1, String> {
    let loc = axis.location().transformed(tr);
    let mut d = *axis.direction().xyz();
    tr.transforms_xyz_dir(&mut d);
    let mut dir = GpDir::from_xyz(&d)
        .map_err(|_| "IGES 120 axis is null after location: gp_Dir::Transform".to_string())?;
    let det = tr.scale_factor() * tr.h_vectorial_part().determinant();
    if det < 0.0 {
        dir.reverse();
    }
    Ok(GpAx1::new(loc, dir))
}

/// Parametric map applied to a 142 or 141 UV curve.
/// `ParamSurface` (`IGESToBRep_TopoSurface.cxx:1846-2070`).
///
/// Entity 120 stored as `Geom_SurfaceOfRevolution` does not enter the
/// cylinder/cone/torus/sphere shift (`cxx:1859-1862`). `TheULength` is reset
/// to 1 at `TransferTopoSurface` (`cxx:191`) and replaced by the line length
/// only for a line generatrix (`cxx:732`).
pub(in crate::iges_read) fn uv_map(
    session: &Session<'_>,
    basis: usize,
    face: &TopoShape,
) -> Result<(GpTrsf2d, f64), String> {
    let mut trans = GpTrsf2d::identity();
    let mut cscale = session.u_length;
    let mut uscale = 1.0;
    let surf = BRepTool::face_surface_world(&crate::shape::Face(face.clone())).ok_or_else(|| {
        "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1851".to_string()
    })?;
    // One rectangular unwrap (`cxx:1853-1857`). A nested trim stays a trim.
    let surf = surf.rectangular_trimmed_basis().unwrap_or(surf);
    let (mut paramu, mut paramv) = elementary_param_shift(session, basis, face, surf.as_ref());
    if paramu.abs() <= CONFUSION {
        paramu = 0.0;
    }
    if paramv.abs() <= CONFUSION {
        paramv = 0.0;
    }
    // Offset unwrap is one level (`cxx:1988-1991`).
    let mut isrf = basis;
    if session.model.entity(isrf)?.type_number == 140 {
        if let EntityData::OffsetSurface { surface, .. } = &session.model.entity(isrf)?.data {
            isrf = *surface;
        }
    }
    let isrf_type = session.model.entity(isrf)?.type_number;
    if isrf_type == 120 {
        if let EntityData::Revolution { end_angle, .. } = &session.model.entity(isrf)?.data {
            let mut tmp = GpTrsf2d::identity();
            tmp.set_translation_vec(&GpVec2d::new(0.0, -2.0 * PI));
            trans.pre_multiply(&tmp);
            tmp = GpTrsf2d::identity();
            tmp.set_mirror_ax2d(&GpAx2d::new(
                GpPnt2d::new(0.0, 0.0),
                GpDir2d::new(1.0, 0.0).expect("ox"),
            ));
            trans.pre_multiply(&tmp);
            tmp = GpTrsf2d::identity();
            tmp.set_mirror_ax2d(&GpAx2d::new(
                GpPnt2d::new(0.0, 0.0),
                GpDir2d::new(1.0, 1.0).expect("diag"),
            ));
            trans.pre_multiply(&tmp);
            uscale = 1.0 / cscale;
            paramu = -(2.0 * PI - *end_angle);
        }
    } else {
        paramu = 0.0; // cxx:2010-2012
    }
    if isrf_type == 118 {
        uscale = 1.0 / cscale; // cxx:2015-2018
    }
    if isrf_type == 122 {
        // `Surf->Bounds` after the rectangular unwrap (`cxx:2021-2032`).
        let (umin, umax) = surf.u_range();
        let uln = (umax - umin).abs();
        uscale = uln / cscale;
        paramu = umin / uln;
    }
    if isrf_type == 192 || isrf_type == 194 {
        uscale = PI / 180.0; // cxx:2035-2038
    }
    if isrf_type == 196 {
        cscale = PI / 180.0; // cxx:2041-2044
        uscale = 1.0;
    }
    if isrf_type == 198 {
        let mut tmp = GpTrsf2d::identity();
        tmp.set_translation_vec(&GpVec2d::new(0.0, -360.0));
        trans.pre_multiply(&tmp);
        tmp = GpTrsf2d::identity();
        tmp.set_mirror_ax2d(&GpAx2d::new(
            GpPnt2d::new(0.0, 0.0),
            GpDir2d::new(1.0, 0.0).expect("ox"),
        ));
        trans.pre_multiply(&tmp);
        tmp = GpTrsf2d::identity();
        tmp.set_mirror_ax2d(&GpAx2d::new(
            GpPnt2d::new(0.0, 0.0),
            GpDir2d::new(1.0, 1.0).expect("diag"),
        ));
        trans.pre_multiply(&tmp);
        if paramv > 0.0 {
            paramv = paramv * 180.0 / PI; // cxx:2056-2058
        }
        cscale = PI / 180.0;
        uscale = 1.0;
    }
    let mut tmp = GpTrsf2d::identity();
    tmp.set_translation_pnts(&GpPnt2d::new(0.0, 0.0), &GpPnt2d::new(paramu, paramv));
    trans.pre_multiply(&tmp);
    tmp = GpTrsf2d::identity();
    tmp.set_scale(&GpPnt2d::new(0.0, 0.0), cscale)
        .map_err(|e| e.to_string())?;
    trans.pre_multiply(&tmp);
    Ok((trans, uscale))
}

/// Cylinder, cone, torus, or sphere parameter origin (`cxx:1859-1974`).
/// A missing V-direction `Geom2d_Line` leaves the shift at 0: `Curve` on a
/// null edge at `cxx:1882` is not given a stand-in point.
fn elementary_param_shift(
    session: &Session<'_>,
    entity: usize,
    face: &TopoShape,
    surf: &dyn Surface,
) -> (f64, f64) {
    if surf.rectangular_trimmed_basis().is_some() {
        return (0.0, 0.0);
    }
    let kind = if surf.gp_cylinder().is_some() {
        1
    } else if surf.gp_cone().is_some() {
        2
    } else if surf.gp_torus().is_some() {
        3
    } else if surf.gp_sphere().is_some() {
        4
    } else {
        0
    };
    if kind == 0 {
        return (0.0, 0.0);
    }
    let dy = GpDir2d::new(0.0, 1.0).expect("dy");
    let face_key = GeometryRegistry::shape_key(face);
    let reg = GeometryRegistry::global();
    let mut edges = Vec::new();
    collect_edges(face, &mut edges);
    let mut chosen: Option<crate::shape::Edge> = None;
    for edge in &edges {
        let Some(curve2d) = reg.edge_pcurve(edge, face_key) else {
            // `Curve2d->IsKind` on a null pcurve (`cxx:1870`).
            continue;
        };
        if curve2d.is_geom2d_line() {
            if let Some(line) = curve2d.gp_lin2d() {
                if line.direction().is_parallel(&dy, ANGULAR) {
                    chosen = Some(crate::shape::Edge(edge.clone()));
                    break;
                }
            }
        }
    }
    let Some(edge) = chosen else {
        return (0.0, 0.0);
    };
    let Some(curve3d) = BRepTool::edge_curve_world(&edge) else {
        return (0.0, 0.0);
    };
    let (mut first, _) = BRepTool::edge_parameters(&edge);
    if Precision::is_negative_infinite(first) {
        first = 0.0; // cxx:1883-1885
    }
    let point = curve3d.d0(first);
    let mut paramu = 0.0;
    let mut paramv = 0.0;
    if let Some(cyl) = surf.gp_cylinder() {
        let (u, v) = slib::cylinder_parameters(&cyl.position(), &point);
        paramu = u;
        paramv = v;
    } else if let Some(cone) = surf.gp_cone() {
        let (u, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &point);
        paramu = u;
        paramv = v;
    } else if let Some(torus) = surf.gp_torus() {
        let (u, v) = slib::torus_parameters(
            torus.position(),
            torus.major_radius(),
            torus.minor_radius(),
            &point,
        );
        paramu = u;
        paramv = v;
    } else if surf.gp_sphere().is_some()
        && session
            .model
            .entity(entity)
            .map(|e| e.type_number != 196)
            .unwrap_or(true)
    {
        if let Some(sphere) = surf.gp_sphere() {
            let (u, v) = slib::sphere_parameters(sphere.position(), &point);
            paramu = u;
            paramv = v;
        }
    }
    if let Ok(ent) = session.model.entity(entity) {
        if ent.type_number == 120 {
            if let EntityData::Revolution { generatrix, .. } = &ent.data {
                if let Ok(gen) = session.model.entity(*generatrix) {
                    if let EntityData::CircularArc { center, start, .. } = &gen.data {
                        paramv -= circ_arc_parameter(center, start); // cxx:1943-1949
                        if surf.gp_sphere().is_some() {
                            paramv += crate::shhealing::adjust_by_period(paramv, 0.0, PI * 2.0);
                        } else if surf.gp_torus().is_some() {
                            paramv += crate::shhealing::adjust_by_period(paramv, PI, PI * 2.0);
                        }
                    }
                }
            }
        } else if ent.type_number == 122 {
            if let EntityData::Tabulated { directrix, .. } = &ent.data {
                if let Ok(dir) = session.model.entity(*directrix) {
                    if let EntityData::CircularArc { center, start, .. } = &dir.data {
                        paramu -= circ_arc_parameter(center, start); // cxx:1964-1970
                        paramu += crate::shhealing::adjust_by_period(paramu, PI, PI * 2.0);
                    }
                }
            }
        }
    }
    (paramu, paramv)
}

fn circ_arc_parameter(center: &[f64; 2], start: &[f64; 2]) -> f64 {
    let dx = center[0] - start[0];
    let dy = center[1] - start[1];
    let radius = (dx * dx + dy * dy).sqrt();
    let ax = GpAx22d::from_xdir(
        GpPnt2d::new(center[0], center[1]),
        GpDir2d::new(1.0, 0.0).expect("ox"),
    );
    clib2d::parameter_circ2d(&GpCirc2d::new(ax, radius), &GpPnt2d::new(start[0], start[1]))
}

fn apply_trimmed_transform(
    session: &Session<'_>,
    index: usize,
    face: &mut TopoShape,
) -> Result<(), String> {
    // `IGESToBRep_TopoSurface.cxx:1321-1371`.
    let Some(g) = entity_location(session.model, index)? else {
        return Ok(());
    };
    // `s1`/`s2`/`s3` are squared column lengths (`cxx:1328-1334`), not moduli.
    let s1 = col_sq(&g, 1);
    let s2 = col_sq(&g, 2);
    let s3 = col_sq(&g, 3);
    if (s1 - s2).abs() > CONFUSION || (s1 - s3).abs() > CONFUSION {
        // `BRepBuilderAPI_GTransform::Perform(face, true)` (`cxx:1336-1343`).
        return super::gtransform::gtransform_face(face, &g);
    }
    let tx = g.value(1, 4);
    let ty = g.value(2, 4);
    let tz = g.value(3, 4);
    let tmp = (g.value(1, 1) - 1.0).abs()
        + g.value(1, 2).abs()
        + g.value(1, 3).abs()
        + g.value(2, 1).abs()
        + (g.value(2, 2) - 1.0).abs()
        + g.value(2, 3).abs()
        + g.value(3, 1).abs()
        + g.value(3, 2).abs()
        + (g.value(3, 3) - 1.0).abs()
        + (tx * tx + ty * ty + tz * tz).sqrt();
    if tmp <= CONFUSION {
        return Ok(());
    }
    let mut tr = GpTrsf::identity();
    tr.set_values(
        g.value(1, 1),
        g.value(1, 2),
        g.value(1, 3),
        tx,
        g.value(2, 1),
        g.value(2, 2),
        g.value(2, 3),
        ty,
        g.value(3, 1),
        g.value(3, 2),
        g.value(3, 3),
        tz,
    )
    .map_err(|e| e.to_string())?;
    move_face(face, &tr);
    Ok(())
}

fn col_sq(g: &occt_core::gp::GpGTrsf, col: usize) -> f64 {
    let x = g.value(1, col);
    let y = g.value(2, col);
    let z = g.value(3, col);
    x * x + y * y + z * z
}

fn move_face(face: &TopoShape, tr: &GpTrsf) {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(face);
    let old_key = reg.debug_repr_key(face_key);
    if let Some(surf) = reg.face_surface(face) {
        let mut moved = surf.clone_dyn();
        moved.transform(tr);
        reg.set_face(face, crate::tgeometry::FaceGeom::new(Arc::from(moved)));
    }
    let new_key = reg.debug_repr_key(face_key);
    let wires = face.tshape.read().expect("lock").children.clone();
    for wire in wires {
        let edges = wire.tshape.read().expect("lock").children.clone();
        for edge in edges {
            if !edge.is_edge() {
                continue;
            }
            if old_key != new_key {
                let mut ts = edge.tshape.write().expect("lock");
                let slot = ts.edge_pcurves_mut();
                if let Some(curves) = slot.curves.remove(&old_key) {
                    slot.curves.insert(new_key, curves);
                }
                if let Some(range) = slot.ranges.remove(&old_key) {
                    slot.ranges.insert(new_key, range);
                }
            }
            if let Some(curve) = reg.edge_curve(&edge) {
                let moved = Arc::from(curve.transformed(tr));
                edge.tshape.write().expect("lock").edge_core_mut().curve = Some(moved);
            }
            let verts = edge.tshape.read().expect("lock").children.clone();
            for vert in verts {
                if vert.is_vertex() {
                    let p = reg.vertex_point(&vert).transformed(tr);
                    crate::shape::Vertex(vert).set_point(p);
                }
            }
        }
    }
}

fn transfer_plane(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferPlane` / `TransferPlaneParts` (`cxx:1493-1516`, `:1608-1785`).
    // Geom equation is `ax+by+cz+d = 0` with `d = -D` (`cxx:1624-1628`).
    // `TransferPlane` calls `TransferPlaneParts` with `first = true`.
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
    let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln));
    let mut face = make_face(surface, &[]);
    // `HasBoundingCurve` is a non-null curve pointer (`IGESGeom_Plane.cxx:70-72`).
    // A null pointer leaves the face without a wire.
    let Some(curve) = bounding else {
        return finish_plane(session, index, face);
    };
    if session.model.entity(curve).is_err() {
        // `cxx:1665-1671`: null bounding curve, keep the infinite plane.
        return finish_plane(session, index, face);
    }
    let curve_type = session.model.entity(curve)?.type_number;
    if !is_topo_curve(curve_type) {
        // `cxx:1765-1773`: not a topo curve. The warning leaves the plane untrimmed.
        return finish_plane(session, index, face);
    }
    if curve_type == 142 {
        // `cxx:1677-1691`. The local `gp_Trsf` is identity, so the `gp_Trsf2d`
        // conversion (`gp_Trsf2d.hxx:232-241`) is identity. `TheULength` is the
        // value reset at `TransferTopoSurface` (`cxx:191`).
        let trans = GpTrsf2d::identity();
        if let Ok(wire) =
            transfer_curve_on_face(session, curve, &face.0, &trans, session.u_length, true)
        {
            TopoBuilder::new().add_wire(&mut face, &wire);
        }
        return finish_plane(session, index, face);
    }
    let shape = match transfer_topo_curve(session, curve) {
        Ok(shape) => shape,
        Err(_) => return finish_plane(session, index, face),
    };
    let mut wire = if shape.is_edge() {
        TopoBuilder::new().make_wire(&[Edge(shape)])
    } else if shape.is_wire() {
        Wire(shape)
    } else {
        // `cxx:1710-1725`, `first` is true: return the plane, not the curve.
        return finish_plane(session, index, face);
    };
    // Naked `ShapeFix_Wire` (`cxx:1729-1732`). Default precision is
    // `Precision::Confusion` (`ShapeFix_Root.cxx:26`). `FixConnected(GetMaxTol)`
    // passes that value as the check precision (`ShapeFix_Wire.cxx:557-569`).
    let eps = eps_geom(session.model.resolution) * unit;
    let _ = fix_connected_all(&mut wire, CONFUSION, eps.max(1.0));
    // `BRepLib_MakeFace(pln, wire, Inside=false)` (`cxx:1733`, `BRepLib_MakeFace.cxx:262-270`)
    // uses the unscaled plane and only feeds `BRepGProp::SurfaceProperties`.
    let probe = TopoBuilder::new().make_face(Arc::new(GeomPlane::new(pln_raw)), &[wire.clone()]);
    // `BRepGProp::SurfaceProperties` (`cxx:1741-1752`). A missing mass leaves
    // the wire as transferred; the sign test is the only reverse.
    if let Ok((props, _)) = surface_properties(&probe.0) {
        let hole = form == -1;
        if props.mass() < 0.0 {
            if !hole {
                wire.0.reverse();
            }
        } else if hole {
            wire.0.reverse();
        }
    }
    TopoBuilder::new().add_wire(&mut face, &wire);
    finish_plane(session, index, face)
}

/// `TransferPlane` applies the matrix computed in `TransferPlaneParts`
/// (`IGESToBRep_TopoSurface.cxx:1511-1515`, `:1647-1658`). A failed
/// `ConvertLocation` leaves the transform identity and does not move the face.
fn finish_plane(
    session: &Session<'_>,
    index: usize,
    mut face: crate::shape::Face,
) -> Result<TopoShape, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if tr.form() != TrsfForm::Identity {
                face.0
                    .set_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
            }
        }
    }
    Ok(face.0)
}

/// `BRep_Builder::MakeFace` + `B.UpdateFace(plane, surf, TopLoc_Location(),
/// Precision::Confusion())` (`IGESToBRep_TopoSurface.cxx:302-303`, `:1641`),
/// `BRepLib_MakeFace(surf, Precision::Confusion())` (`cxx:308`, `:1186`) and
/// `BRepBuilderAPI_MakeFace::Init(surf, ..., Precision::Confusion())`
/// (`cxx:811-815`, `:961-966`): every face this module creates carries
/// `Precision::Confusion()`, not the `BRep_TFace` default `RealEpsilon()`.
/// `BRepBndLib::Add` feeds that tolerance to `BndLib_AddSurface::Add`
/// (`BRepBndLib.cxx:109`), so it sets `Bnd_Box::Gap` and therefore
/// `Prs3d::GetDeflection`.
fn make_face(surface: Arc<dyn Surface>, wires: &[Wire]) -> crate::shape::Face {
    TopoBuilder::new().make_face_with_tolerance(surface, wires, CONFUSION)
}
