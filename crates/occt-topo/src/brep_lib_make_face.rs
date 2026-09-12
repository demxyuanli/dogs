//! `BRepLib_MakeFace::Init` natural bounds for a surface UV rectangle.
//!
//! Source: `BRepLib_MakeFace.cxx` (`Init(S, Bound)` at 366, `Init` UV at 463,
//! `IsDegenerated` at 392). `BRepBuilderAPI_MakeFace(S, TolDegen)` uses this
//! path. Offset SameParameter (`cxx:860-866`) is `brep_lib_same_parameter`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpDir2d, GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::{surface_ops, Curve, GeomCircle, GeomRectangularTrimmedSurface, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dLine;

use crate::abs::Orientation;
use crate::builder::TopoBuilder;
use crate::geom_bnd_lib_elclib2d::adjust_periodic;
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};

/// `BRepLib_MakeFace::IsDegenerated` (`cxx:392-458`).
fn is_degenerated(curve: &dyn Curve, max_tol: f64) -> (bool, f64) {
    let confusion = Precision::CONFUSION;
    if let Some(radius) = curve.circle_radius() {
        if radius > max_tol {
            return (false, confusion);
        }
        return (true, radius.max(confusion));
    }
    if let Some(poles) = curve.bspline_poles() {
        return poles_collapsed(poles, max_tol, confusion);
    }
    if let Some(poles) = curve.bezier_poles() {
        return poles_collapsed(poles, max_tol, confusion);
    }
    (false, confusion)
}

/// Shared BSpline/Bezier pole-collapse (`cxx:411-456`).
fn poles_collapsed(poles: &[GpPnt], max_tol: f64, confusion: f64) -> (bool, f64) {
    if poles.is_empty() {
        return (false, confusion);
    }
    let p1 = poles[0];
    let max_tol2 = max_tol * max_tol;
    let mut max_pole_dist2 = 0.0;
    for p2 in poles.iter().skip(1) {
        let d2 = p1.square_distance(p2);
        if d2 > max_tol2 {
            return (false, confusion);
        }
        if d2 > max_pole_dist2 {
            max_pole_dist2 = d2;
        }
    }
    (true, (1.000001 * max_pole_dist2.sqrt()).max(confusion))
}

/// OCCT `UIso(u)`: fixed U, vary V. Sphere uses `ElSLib::SphereUIso`.
fn surface_u_iso(s: &dyn Surface, u: f64) -> Arc<dyn Curve> {
    if let Some(sp) = s.gp_sphere() {
        Arc::new(GeomCircle::new(slib::sphere_u_iso(
            sp.position(),
            sp.radius(),
            u,
        )))
    } else if let Some(c) = s.u_iso_curve(u) {
        c
    } else {
        surface_ops::iso_v(s, u)
    }
}

/// OCCT `VIso(v)`: fixed V, vary U. Sphere uses `ElSLib::SphereVIso`.
fn surface_v_iso(s: &dyn Surface, v: f64) -> Arc<dyn Curve> {
    if let Some(sp) = s.gp_sphere() {
        Arc::new(GeomCircle::new(slib::sphere_v_iso(
            sp.position(),
            sp.radius(),
            v,
        )))
    } else if let Some(c) = s.v_iso_curve(v) {
        c
    } else {
        surface_ops::iso_u(s, v)
    }
}

fn dir2d_x() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("unit x")
}

fn dir2d_y() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("unit y")
}

fn make_bound_edge(curve: Arc<dyn Curve>, first: f64, last: f64, degen: bool, tol: f64) -> Edge {
    let e = Edge::new();
    let mut g = EdgeGeom::new(curve, first, last);
    g.tolerance = tol;
    g.degenerated = degen;
    GeometryRegistry::global().set_edge(&e.0, g);
    e
}

fn add_vertices(b: &TopoBuilder, e: &mut Edge, first: Option<&Vertex>, last: Option<&Vertex>) {
    if let Some(v) = first {
        let mut vf = v.0.clone();
        vf.set_orientation(Orientation::Forward);
        b.add(&mut e.0, &vf);
    }
    if let Some(v) = last {
        let mut vl = v.0.clone();
        vl.set_orientation(Orientation::Reversed);
        b.add(&mut e.0, &vl);
    }
}

/// `BRepLib_MakeFace::Init(S, Bound=true, TolDegen)` (`cxx:366-376`, `463-858`).
pub fn make_face_from_surface(surface: Arc<dyn Surface>, tol_degen: f64) -> Face {
    let (u_min, u_max) = surface.u_range();
    let (v_min, v_max) = surface.v_range();
    make_face_uv(surface, u_min, u_max, v_min, v_max, tol_degen)
}

/// `BRepLib_MakeFace::Init(S, Um, UM, Vm, VM, TolDegen)` (`cxx:463`).
pub fn make_face_uv(
    surface: Arc<dyn Surface>,
    um: f64,
    u_m: f64,
    vm: f64,
    v_m: f64,
    tol_degen: f64,
) -> Face {
    let mut u_min = um;
    let mut u_max = u_m;
    let mut v_min = vm;
    let mut v_max = v_m;

    // `BRepLib_MakeFace.cxx:479-518` — unwrap RTS, recut Offset of extrusion/revolution.
    let mut surface = surface;
    let bs = surface
        .rectangular_trimmed_basis()
        .unwrap_or_else(|| surface.clone());
    let offset_surface = bs.is_offset_surface();
    let epsilon = Precision::PCONFUSION;
    let (nat_umin, nat_umax) = bs.u_range();
    let (nat_vmin, nat_vmax) = bs.v_range();
    if offset_surface {
        if let Some(base) = bs.offset_basis_surface() {
            if base.is_surface_of_linear_extrusion() {
                if Precision::is_infinite(nat_umin) || Precision::is_infinite(nat_umax) {
                    surface = Arc::new(GeomRectangularTrimmedSurface::uv(
                        bs.clone(),
                        u_min,
                        u_max,
                        v_min,
                        v_max,
                    ));
                } else {
                    surface = Arc::new(GeomRectangularTrimmedSurface::one_param(
                        bs.clone(),
                        v_min,
                        v_max,
                        false,
                    ));
                }
            } else if base.is_surface_of_revolution() {
                if Precision::is_infinite(nat_vmin) || Precision::is_infinite(nat_vmax) {
                    surface = Arc::new(GeomRectangularTrimmedSurface::one_param(
                        bs.clone(),
                        v_min,
                        v_max,
                        false,
                    ));
                }
            }
        }
    }

    if surface.is_u_periodic() {
        adjust_periodic(nat_umin, nat_umax, epsilon, &mut u_min, &mut u_max);
    } else if u_min > u_max {
        std::mem::swap(&mut u_min, &mut u_max);
        if (nat_umin - u_min > epsilon) || (u_max - nat_umax > epsilon) {
            return TopoBuilder::new().make_face(surface, &[]);
        }
    }
    if surface.is_v_periodic() {
        adjust_periodic(nat_vmin, nat_vmax, epsilon, &mut v_min, &mut v_max);
    } else if v_min > v_max {
        std::mem::swap(&mut v_min, &mut v_max);
        if (nat_vmin - v_min > epsilon) || (v_max - nat_vmax > epsilon) {
            return TopoBuilder::new().make_face(surface, &[]);
        }
    }

    let umininf = Precision::is_negative_infinite(u_min);
    let umaxinf = Precision::is_positive_infinite(u_max);
    let vmininf = Precision::is_negative_infinite(v_min);
    let vmaxinf = Precision::is_positive_infinite(v_max);

    let uclosed = surface.is_u_closed()
        && (u_min - nat_umin).abs() < epsilon
        && (u_max - nat_umax).abs() < epsilon;
    let vclosed = surface.is_v_closed()
        && (v_min - nat_vmin).abs() < epsilon
        && (v_max - nat_vmax).abs() < epsilon;

    let max_tol = tol_degen;
    let mut umin_tol = Precision::CONFUSION;
    let mut umax_tol = Precision::CONFUSION;
    let mut vmin_tol = Precision::CONFUSION;
    let mut vmax_tol = Precision::CONFUSION;
    let mut dumin = false;
    let mut dumax = false;
    let mut dvmin = false;
    let mut dvmax = false;
    let mut c_umin: Option<Arc<dyn Curve>> = None;
    let mut c_umax: Option<Arc<dyn Curve>> = None;
    let mut c_vmin: Option<Arc<dyn Curve>> = None;
    let mut c_vmax: Option<Arc<dyn Curve>> = None;

    if !umininf {
        let c = surface_u_iso(surface.as_ref(), u_min);
        let (d, t) = is_degenerated(c.as_ref(), max_tol);
        dumin = d;
        umin_tol = t;
        c_umin = Some(c);
    }
    if !umaxinf {
        let c = surface_u_iso(surface.as_ref(), u_max);
        let (d, t) = is_degenerated(c.as_ref(), max_tol);
        dumax = d;
        umax_tol = t;
        c_umax = Some(c);
    }
    if !vmininf {
        let c = surface_v_iso(surface.as_ref(), v_min);
        let (d, t) = is_degenerated(c.as_ref(), max_tol);
        dvmin = d;
        vmin_tol = t;
        c_vmin = Some(c);
    }
    if !vmaxinf {
        let c = surface_v_iso(surface.as_ref(), v_max);
        let (d, t) = is_degenerated(c.as_ref(), max_tol);
        dvmax = d;
        vmax_tol = t;
        c_vmax = Some(c);
    }

    let b = TopoBuilder::new();
    let mut v00: Option<Vertex> = None;
    let mut v10: Option<Vertex> = None;
    let mut v11: Option<Vertex> = None;
    let mut v01: Option<Vertex> = None;
    if !umininf {
        if !vmininf {
            v00 = Some(b.make_vertex(
                surface.value(u_min, v_min),
                umin_tol.max(vmin_tol),
            ));
        }
        if !vmaxinf {
            v01 = Some(b.make_vertex(
                surface.value(u_min, v_max),
                umin_tol.max(vmax_tol),
            ));
        }
    }
    if !umaxinf {
        if !vmininf {
            v10 = Some(b.make_vertex(
                surface.value(u_max, v_min),
                umax_tol.max(vmin_tol),
            ));
        }
        if !vmaxinf {
            v11 = Some(b.make_vertex(
                surface.value(u_max, v_max),
                umax_tol.max(vmax_tol),
            ));
        }
    }
    if uclosed {
        v10 = v00.clone();
        v11 = v01.clone();
    }
    if vclosed {
        v01 = v00.clone();
        v11 = v10.clone();
    }
    if dumin {
        v00 = v01.clone();
    }
    if dumax {
        v10 = v11.clone();
    }
    if dvmin {
        v00 = v10.clone();
    }
    if dvmax {
        v01 = v11.clone();
    }

    let l_umin: Option<Arc<dyn Curve2d>> = if !umininf {
        Some(Arc::new(Geom2dLine::from_pnt_dir(
            GpPnt2d::new(u_min, 0.0),
            dir2d_y(),
        )))
    } else {
        None
    };
    let l_umax: Option<Arc<dyn Curve2d>> = if !umaxinf {
        Some(Arc::new(Geom2dLine::from_pnt_dir(
            GpPnt2d::new(u_max, 0.0),
            dir2d_y(),
        )))
    } else {
        None
    };
    let l_vmin: Option<Arc<dyn Curve2d>> = if !vmininf {
        Some(Arc::new(Geom2dLine::from_pnt_dir(
            GpPnt2d::new(0.0, v_min),
            dir2d_x(),
        )))
    } else {
        None
    };
    let l_vmax: Option<Arc<dyn Curve2d>> = if !vmaxinf {
        Some(Arc::new(Geom2dLine::from_pnt_dir(
            GpPnt2d::new(0.0, v_max),
            dir2d_x(),
        )))
    } else {
        None
    };

    let mut face = b.make_face(surface, &[]);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();

    let mut e_umin: Option<Edge> = None;
    let mut e_umax: Option<Edge> = None;
    let mut e_vmin: Option<Edge> = None;
    let mut e_vmax: Option<Edge> = None;

    if !umininf {
        if let Some(c) = c_umin {
            let mut e = make_bound_edge(c, v_min, v_max, dumin, umin_tol);
            if uclosed {
                if let (Some(lumax), Some(lumin)) = (l_umax.clone(), l_umin.clone()) {
                    reg.set_edge_pcurves(&e.0, face_key, vec![lumax, lumin]);
                    if let Some(mut g) = reg.edge_geom(&e.0) {
                        g.tolerance = umin_tol.max(umax_tol);
                        reg.set_edge(&e.0, g);
                    }
                }
            } else if let Some(lumin) = l_umin.clone() {
                reg.set_edge_pcurve(&e.0, face_key, lumin);
            }
            add_vertices(&b, &mut e, v00.as_ref(), v01.as_ref());
            e_umin = Some(e);
        }
    }
    if !umaxinf {
        if uclosed {
            e_umax = e_umin.clone();
        } else if let Some(c) = c_umax {
            let mut e = make_bound_edge(c, v_min, v_max, dumax, umax_tol);
            if let Some(lumax) = l_umax.clone() {
                reg.set_edge_pcurve(&e.0, face_key, lumax);
            }
            add_vertices(&b, &mut e, v10.as_ref(), v11.as_ref());
            e_umax = Some(e);
        }
    }
    if !vmininf {
        if let Some(c) = c_vmin {
            let mut e = make_bound_edge(c, u_min, u_max, dvmin, vmin_tol);
            if vclosed {
                if let (Some(lvmin), Some(lvmax)) = (l_vmin.clone(), l_vmax.clone()) {
                    reg.set_edge_pcurves(&e.0, face_key, vec![lvmin, lvmax]);
                    if let Some(mut g) = reg.edge_geom(&e.0) {
                        g.tolerance = vmin_tol.max(vmax_tol);
                        reg.set_edge(&e.0, g);
                    }
                }
            } else if let Some(lvmin) = l_vmin.clone() {
                reg.set_edge_pcurve(&e.0, face_key, lvmin);
            }
            add_vertices(&b, &mut e, v00.as_ref(), v10.as_ref());
            e_vmin = Some(e);
        }
    }
    if !vmaxinf {
        if vclosed {
            e_vmax = e_vmin.clone();
        } else if let Some(c) = c_vmax {
            let mut e = make_bound_edge(c, u_min, u_max, dvmax, vmax_tol);
            if let Some(lvmax) = l_vmax.clone() {
                reg.set_edge_pcurve(&e.0, face_key, lvmax);
            }
            add_vertices(&b, &mut e, v01.as_ref(), v11.as_ref());
            e_vmax = Some(e);
        }
    }

    if let Some(e) = e_umin.as_mut() {
        e.0.set_orientation(Orientation::Reversed);
    }
    if let Some(e) = e_vmax.as_mut() {
        e.0.set_orientation(Orientation::Reversed);
    }

    if !umininf && !umaxinf && vmininf && vmaxinf {
        if let Some(e) = e_umin {
            let w = b.make_wire(&[e]);
            b.add_wire(&mut face, &w);
        }
        if let Some(e) = e_umax {
            let w = b.make_wire(&[e]);
            b.add_wire(&mut face, &w);
        }
        face.0.set_closed(uclosed);
    } else if umininf && umaxinf && !vmininf && !vmaxinf {
        if let Some(e) = e_vmin {
            let w = b.make_wire(&[e]);
            b.add_wire(&mut face, &w);
        }
        if let Some(e) = e_vmax {
            let w = b.make_wire(&[e]);
            b.add_wire(&mut face, &w);
        }
        face.0.set_closed(vclosed);
    } else if !umininf || !umaxinf || !vmininf || !vmaxinf {
        let mut edges: Vec<Edge> = Vec::new();
        if let Some(e) = e_umin {
            edges.push(e);
        }
        if let Some(e) = e_vmin {
            edges.push(e);
        }
        if let Some(e) = e_umax {
            edges.push(e);
        }
        if let Some(e) = e_vmax {
            edges.push(e);
        }
        let w = b.make_wire(&edges);
        w.0.set_closed(!umininf && !umaxinf && !vmininf && !vmaxinf);
        b.add_wire(&mut face, &w);
        face.0.set_closed(uclosed && vclosed);
    }

    if offset_surface {
        crate::brep_lib_same_parameter::same_parameter_face(
            &face,
            Precision::CONFUSION,
            true,
        );
    }

    face
}
