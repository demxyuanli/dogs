use std::f64::consts::PI;

use super::*;

use occt_core::elib::{clib2d, slib};
use occt_core::gp::{
    GpAx22d, GpCirc2d, GpDir2d, GpPnt2d,
};
use occt_core::precision::{Precision, ANGULAR};
use occt_geom::Surface;

use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::TopoShape;
use occt_topo::tgeometry::GeometryRegistry;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;

/// Cylinder, cone, torus, or sphere parameter origin (`cxx:1859-1974`).
/// A missing V-direction `Geom2d_Line` leaves the shift at 0: `Curve` on a
/// null edge at `cxx:1882` is not given a stand-in point.
pub(super) fn elementary_param_shift(
    session: &Session<'_>,
    entity: usize,
    face: &TopoShape,
    surf: &dyn Surface,
) -> (f64, f64) {
    if surf.rectangular_trimmed_basis().is_some() {
        return (0.0, 0.0);
    }
    // `cxx:1848-1974`: `paramu`/`paramv` are filled inside a single
    // `if (Surf->IsKind(Cylindrical|Conical|Toroidal|Spherical))` block
    // (`cxx:1859-1974`, the revolution and tabulated arms at `cxx:1939-1974` are
    // nested in it, four spaces deep), so a `Geom_SurfaceOfRevolution`
    // (`cxx:1993`) gets no shift from this pass.
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
    let mut chosen: Option<occt_topo::shape::Edge> = None;
    for edge in &edges {
        let Some(curve2d) = reg.edge_pcurve(edge, face_key) else {
            // `Curve2d->IsKind` on a null pcurve (`cxx:1870`).
            continue;
        };
        if curve2d.is_geom2d_line() {
            if let Some(line) = curve2d.gp_lin2d() {
                if line.direction().is_parallel(&dy, ANGULAR) {
                    chosen = Some(occt_topo::shape::Edge(edge.clone()));
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
                            paramv += occt_topo::shhealing::adjust_by_period(paramv, 0.0, PI * 2.0);
                        } else if surf.gp_torus().is_some() {
                            paramv += occt_topo::shhealing::adjust_by_period(paramv, PI, PI * 2.0);
                        }
                    }
                }
            }
        } else if ent.type_number == 122 {
            if let EntityData::Tabulated { directrix, .. } = &ent.data {
                if let Ok(dir) = session.model.entity(*directrix) {
                    if let EntityData::CircularArc { center, start, .. } = &dir.data {
                        paramu -= circ_arc_parameter(center, start); // cxx:1964-1970
                        paramu += occt_topo::shhealing::adjust_by_period(paramu, PI, PI * 2.0);
                    }
                }
            }
        }
    }
    (paramu, paramv)
}

pub(super) fn circ_arc_parameter(center: &[f64; 2], start: &[f64; 2]) -> f64 {
    let dx = center[0] - start[0];
    let dy = center[1] - start[1];
    let radius = (dx * dx + dy * dy).sqrt();
    let ax = GpAx22d::from_xdir(
        GpPnt2d::new(center[0], center[1]),
        GpDir2d::new(1.0, 0.0).expect("ox"),
    );
    clib2d::parameter_circ2d(&GpCirc2d::new(ax, radius), &GpPnt2d::new(start[0], start[1]))
}
