//! `ComputeCurve3d` (`BRepOffset_Offset.cxx:190-385`) and the three local
//! `UpdateEdge` helpers of `BRepOffset_Offset.cxx` (`:115-183`).
//!
//! `ComputeCurve3d` only builds a 3D curve for the iso cases of an analytic
//! surface and for a p-curve on a plane; every other case is deliberately left
//! without a 3D curve (`cxx:376-384`, "les courbes 3d ... ne sont calculees
//! que si necessaire"). That is what makes the offset face's edges carry
//! exact circles and lines.
//!
//! The `TopLoc_Location` argument of OCCT's helpers is not reproduced: OCCT
//! stores it on the `BRep_GCurve` representation (`BRep_Builder::UpdateEdge(E,
//! C, L, Tol)`, `BRep_Builder.cxx:1180-1191`), and this crate stores the
//! representation location on the edge view instead. `init_face` puts the
//! source `L` there when it creates the edge, so the local curve written by
//! these helpers is evaluated in the same frame as OCCT's.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx1, GpDir, GpDir2d, GpLin};
use occt_core::precision::{ANGULAR, PCONFUSION};
use occt_geom::{Curve, GeomCircle, GeomLine, Surface};
use occt_geom2d::Curve2d;

use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `gp::DX2d()` (`gp.hxx`).
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };
/// `gp::DY2d()` (`gp.hxx`).
const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };

/// `GeomAdaptor_Surface::GetType()` for the five surface types
/// `ComputeCurve3d` switches on. The queries reach through a
/// `Geom_RectangularTrimmedSurface`, which is exactly what
/// `GeomAdaptor_Surface::load` does (`GeomAdaptor_Surface.cxx:423-425`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceType {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

fn surface_type(s: &Arc<dyn Surface>) -> SurfaceType {
    if s.gp_pln().is_some() {
        SurfaceType::Plane
    } else if s.gp_sphere().is_some() {
        SurfaceType::Sphere
    } else if s.gp_cylinder().is_some() {
        SurfaceType::Cylinder
    } else if s.gp_cone().is_some() {
        SurfaceType::Cone
    } else if s.gp_torus().is_some() {
        SurfaceType::Torus
    } else {
        SurfaceType::Other
    }
}

/// The axis direction `gp_Dir` of `Dir1.Crossed(Dir2)` (`cxx:227`, `:245`):
/// the cross product of two `gp_Ax3` directions. A null product cannot come
/// out of an `gp_Ax3`, so the fallback keeps the first direction.
fn crossed_dir(a: &GpDir, b: &GpDir) -> GpDir {
    a.crossed(b).unwrap_or(*a)
}

/// `BRep_Builder::UpdateEdge(E, C, L, Tol)` (`cxx:115-131`): a
/// `Geom_TrimmedCurve` is replaced by its basis curve, then the 3D curve and
/// the raised tolerance are stored.
pub fn update_edge_curve3d(edge: &Edge, curve: &Arc<dyn Curve>, tol: f64) {
    let stored = match curve.trim_basis() {
        Some(basis) => basis,
        None => curve.clone(),
    };
    GeometryRegistry::global().set_edge_curve3d(&edge.0, stored, tol);
}

fn unwrap_2d(curve: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    match curve.trimmed_basis() {
        Some(basis) => Arc::from(basis.clone_dyn()),
        None => curve.clone(),
    }
}

/// `BRep_Builder::UpdateEdge(E, C2d, F, Tol)` (`cxx:135-151`).
pub fn update_edge_pcurve(edge: &Edge, curve: &Arc<dyn Curve2d>, face: &Face, tol: f64) {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    reg.set_edge_pcurve(&edge.0, face_key, unwrap_2d(curve));
    reg.raise_edge_tolerance(&edge.0, tol);
}

/// `BRep_Builder::UpdateEdge(E, C1, C2, F, Tol)` (`cxx:155-183`), the seam
/// overload: `C1` is the forward representation and `C2` the reversed one.
pub fn update_edge_pcurves(
    edge: &Edge,
    c1: &Arc<dyn Curve2d>,
    c2: &Arc<dyn Curve2d>,
    face: &Face,
    tol: f64,
) {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    reg.set_edge_pcurves(&edge.0, face_key, vec![unwrap_2d(c1), unwrap_2d(c2)]);
    reg.raise_edge_tolerance(&edge.0, tol);
}

/// `ComputeCurve3d(Edge, Curve, Surf, Loc, Tol)`
/// (`BRepOffset_Offset.cxx:190-385`).
pub fn compute_curve3d(edge: &Edge, curve: &Arc<dyn Curve2d>, surf: &Arc<dyn Surface>, tol: f64) {
    let sty = surface_type(surf);
    let line = curve.gp_lin2d();

    // `cxx:208-369`: only a line p-curve on a non-plane analytic surface
    // produces a 3D curve here.
    if sty != SurfaceType::Plane {
        let Some(line) = line else {
            return;
        };
        let d = *line.direction();
        let p = line.location();

        if d.is_parallel(&DX2D, ANGULAR) {
            // `cxx:213-290`: Iso V.
            match sty {
                SurfaceType::Sphere => {
                    // `cxx:217-221`.
                    if (p.y().abs() - std::f64::consts::FRAC_PI_2).abs() < PCONFUSION {
                        GeometryRegistry::global().set_degenerated(&edge.0, true);
                        return;
                    }
                    let Some(sp) = surf.gp_sphere() else {
                        return;
                    };
                    let axis = *sp.position();
                    let mut ci = slib::sphere_v_iso(&axis, sp.radius(), p.y());
                    let drev = crossed_dir(axis.x_direction(), axis.y_direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.x());
                    update_circle(edge, ci, &d, &DX2D, tol);
                }
                SurfaceType::Cylinder => {
                    let Some(cyl) = surf.gp_cylinder() else {
                        return;
                    };
                    let axis = cyl.position();
                    let mut ci = slib::cylinder_v_iso(&axis, cyl.radius(), p.y());
                    let drev = crossed_dir(axis.x_direction(), axis.y_direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.x());
                    update_circle(edge, ci, &d, &DX2D, tol);
                }
                SurfaceType::Cone => {
                    let Some(cone) = surf.gp_cone() else {
                        return;
                    };
                    let axis = cone.position();
                    let mut ci = slib::cone_v_iso(
                        &axis,
                        cone.radius(),
                        cone.semi_angle(),
                        p.y(),
                    );
                    let drev = crossed_dir(axis.x_direction(), axis.y_direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.x());
                    update_circle(edge, ci, &d, &DX2D, tol);
                }
                SurfaceType::Torus => {
                    let Some(tore) = surf.gp_torus() else {
                        return;
                    };
                    let axis = tore.position();
                    let mut ci = slib::torus_v_iso(
                        &axis,
                        tore.major_radius(),
                        tore.minor_radius(),
                        p.y(),
                    );
                    let drev = crossed_dir(axis.x_direction(), axis.y_direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.x());
                    update_circle(edge, ci, &d, &DX2D, tol);
                }
                _ => {}
            }
        } else if d.is_parallel(&DY2D, ANGULAR) {
            // `cxx:291-367`: Iso U.
            match sty {
                SurfaceType::Sphere => {
                    let Some(sp) = surf.gp_sphere() else {
                        return;
                    };
                    let axis = *sp.position();
                    // `cxx:299`: the iso-0 circle.
                    let mut ci = slib::sphere_u_iso(&axis, sp.radius(), 0.0);
                    // `cxx:302-304`: rotation about `XDirection ^ Direction`.
                    let drev = crossed_dir(axis.x_direction(), &axis.direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.y());
                    // `cxx:307-309`: then about `XDirection ^ YDirection`.
                    let drev = crossed_dir(axis.x_direction(), axis.y_direction());
                    ci.rotate(&GpAx1::new(axis.location(), drev), p.x());
                    update_circle(edge, ci, &d, &DY2D, tol);
                }
                SurfaceType::Cylinder => {
                    let Some(cyl) = surf.gp_cylinder() else {
                        return;
                    };
                    let mut l = slib::cylinder_u_iso(&cyl.position(), cyl.radius(), p.x());
                    translate_lin(&mut l, p.y());
                    update_line(edge, l, &d, &DY2D, tol);
                }
                SurfaceType::Cone => {
                    let Some(cone) = surf.gp_cone() else {
                        return;
                    };
                    let mut l = slib::cone_u_iso(
                        &cone.position(),
                        cone.radius(),
                        cone.semi_angle(),
                        p.x(),
                    );
                    translate_lin(&mut l, p.y());
                    update_line(edge, l, &d, &DY2D, tol);
                }
                SurfaceType::Torus => {
                    let Some(tore) = surf.gp_torus() else {
                        return;
                    };
                    let axis = tore.position();
                    let mut ci = slib::torus_u_iso(
                        &axis,
                        tore.major_radius(),
                        tore.minor_radius(),
                        p.x(),
                    );
                    let cax = ci.axis();
                    ci.rotate(&cax, p.y());
                    update_circle(edge, ci, &d, &DY2D, tol);
                }
                _ => {}
            }
        }
        return;
    }

    // `cxx:370-375`: the plane arm, `GeomAPI::To3d(Curve, S.Plane())`.
    let Some(pln) = surf.gp_pln() else {
        return;
    };
    let axes = pln.position().ax2();
    let Some(c3d) = occt_geom::geom_lib::to_3d(&axes, curve.as_ref()) else {
        return;
    };
    update_edge_curve3d(edge, &c3d, tol);
}

/// `cxx:230-235`: wrap the circle, reverse it when the p-curve line is
/// opposite to the iso direction, and store it.
fn update_circle(
    edge: &Edge,
    ci: occt_core::gp::GpCirc,
    line_dir: &GpDir2d,
    iso_dir: &GpDir2d,
    tol: f64,
) {
    let mut circle = GeomCircle::new(ci);
    if line_dir.is_opposite(iso_dir, ANGULAR) {
        circle.reverse();
    }
    let curve: Arc<dyn Curve> = Arc::new(circle);
    update_edge_curve3d(edge, &curve, tol);
}

/// `cxx:327-332`: wrap the line, reverse it when the p-curve line is opposite
/// to the iso direction, and store it.
fn update_line(edge: &Edge, l: GpLin, line_dir: &GpDir2d, iso_dir: &GpDir2d, tol: f64) {
    let mut line = GeomLine::new(l);
    if line_dir.is_opposite(iso_dir, ANGULAR) {
        line.reverse();
    }
    let curve: Arc<dyn Curve> = Arc::new(line);
    update_edge_curve3d(edge, &curve, tol);
}

/// `L.Translate(Tr)` with `Tr = L.Direction() * amount` (`cxx:324-326`).
fn translate_lin(l: &mut GpLin, amount: f64) {
    let step = l.direction().xyz().multiplied(amount);
    let p = occt_core::gp::GpPnt::from_xyz(&l.location().xyz().added(&step));
    l.set_location(p);
}
