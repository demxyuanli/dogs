//! `BRepLib::SameParameter` for a shape (forced, mutable).
//! Source: `BRepLib.cxx` `InternalSameParameter` 913-1008, edge
//! `SameParameter` 1251-1739, `ComputeTol` 1070-1188.

use std::collections::HashSet;

use occt_core::precision::{CONFUSION, INFINITE, Precision};
use occt_geom::{ApproxSameParameter, Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools_full::{edge_vertices, edges_of, faces_of};

const NCONTROL: i32 = 22;

fn u_resolution(s: &dyn Surface, tol: f64) -> f64 {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let um = 0.5 * (u0 + u1);
    let vm = 0.5 * (v0 + v1);
    let (_, du, _) = s.d1(um, vm);
    let m = du.magnitude();
    if m > 1.0e-12 {
        tol / m
    } else {
        tol
    }
}

fn v_resolution(s: &dyn Surface, tol: f64) -> f64 {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let um = 0.5 * (u0 + u1);
    let vm = 0.5 * (v0 + v1);
    let (_, _, dv) = s.d1(um, vm);
    let m = dv.magnitude();
    if m > 1.0e-12 {
        tol / m
    } else {
        tol
    }
}

/// `BRepLib.cxx:1070-1188`.
fn compute_tol(
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    nbp: i32,
) -> f64 {
    let nbp = nbp.max(1) as usize;
    let mut dist = vec![-1.0; nbp + 10];
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let du = 0.01 * (ul - uf);
    let dv = 0.01 * (vl - vf);
    let is_u_per = surf.is_u_periodic();
    let is_v_per = surf.is_v_periodic();
    let dsdu = 1.0 / u_resolution(surf, 1.0).max(1.0e-30);
    let dsdv = 1.0 / v_resolution(surf, 1.0).max(1.0e-30);
    let mut d2: f64 = 0.0;
    let mut dapp: f64 = -1.0;
    for i in 0..=nbp {
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.d0(u);
        let puv = c2d.d0(u);
        if !is_u_per {
            if puv.x() < uf - du {
                dapp = dapp.max(dsdu * (uf - puv.x()));
                continue;
            } else if puv.x() > ul + du {
                dapp = dapp.max(dsdu * (puv.x() - ul));
                continue;
            }
        }
        if !is_v_per {
            if puv.y() < vf - dv {
                dapp = dapp.max(dsdv * (vf - puv.y()));
                continue;
            } else if puv.y() > vl + dv {
                dapp = dapp.max(dsdv * (puv.y() - vl));
                continue;
            }
        }
        let pcons = surf.d0(puv.x(), puv.y());
        if Precision::is_infinite(pcons.x())
            || Precision::is_infinite(pcons.y())
            || Precision::is_infinite(pcons.z())
        {
            return INFINITE;
        }
        let temp = pc3d.square_distance(&pcons);
        if i + 1 < dist.len() {
            dist[i + 1] = temp;
        }
        d2 = d2.max(temp);
    }
    if Precision::is_infinite(d2) {
        return d2;
    }
    d2 = d2.sqrt();
    if dapp > d2 {
        return dapp;
    }
    let mut n1 = 0i32;
    let mut n2 = 0i32;
    for &di in dist.iter().skip(1) {
        if di > 0.0 {
            if di < 1.0 {
                n1 += 1;
            } else {
                n2 += 1;
            }
        }
    }
    let mut n3 = 0i32;
    if n1 > n2 && n2 != 0 {
        n3 = 100 * n2 / (n1 + n2);
    }
    let mut ana = false;
    let mut d2_ana: f64 = 0.0;
    if n3 < 10 && n3 != 0 {
        ana = true;
        for &di in dist.iter().skip(1) {
            if di > 0.0 && di < 1.0 {
                d2_ana = d2_ana.max(di);
            }
        }
    }
    d2 = if !ana { 1.5 * d2 } else { 1.5 * d2_ana.sqrt() };
    d2.max(1.0e-7)
}

fn same_parameter_edge(edge: &Edge, face: &Face, the_tol: f64) -> f64 {
    let reg = GeometryRegistry::global();
    let Some(mut g) = reg.edge_geom(&edge.0) else {
        return -1.0;
    };
    if g.same_parameter {
        return -1.0;
    }
    if g.degenerated {
        return -1.0;
    }
    let c3d = g.curve.clone();
    let mut f3d = g.first;
    let mut l3d = g.last;
    if !c3d.is_periodic() {
        let udeb = c3d.first_parameter();
        let ufin = c3d.last_parameter();
        if udeb > f3d {
            f3d = udeb;
        }
        if l3d > ufin {
            l3d = ufin;
        }
    }
    let Some(surf) = BRepTool::face_surface(face) else {
        return -1.0;
    };
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcs = reg.edge_pcurves(&edge.0, face_key);
    if pcs.is_empty() {
        g.first = f3d;
        g.last = l3d;
        g.same_range = true;
        reg.set_edge(&edge.0, g);
        return -1.0;
    }
    let mut is_same_p = true;
    let mut maxdist = 0.0;
    let an_edge_tol = g.tolerance;
    let same_range = g.same_range;
    let mut ya_pcu = false;
    const BIG_ERROR: f64 = 1.0e10;
    for pc in &pcs {
        ya_pcu = true;
        let _ = same_range;
        // Unported: `GeomLib::SameRange` when `!SameRange` (`cxx:1387-1392`).
        let error = compute_tol(c3d.as_ref(), pc.as_ref(), surf.as_ref(), f3d, l3d, NCONTROL);
        if error > BIG_ERROR {
            maxdist = error;
            break;
        }
        let same_p = ApproxSameParameter::new(
            c3d.as_ref(),
            pc.as_ref(),
            surf.as_ref(),
            f3d,
            l3d,
            the_tol,
        );
        if same_p.same_parameter {
            maxdist = maxdist.max(same_p.tol_reached);
        } else if same_p.done {
            maxdist = maxdist.max(same_p.tol_reached.min(error));
        } else {
            is_same_p = false;
        }
        if !is_same_p {
            // `BRepCheck::PrecCurve` / `PrecSurface` (`BRepCheck.cxx:70-129`):
            // ellipse/cone epsilon specials are not ported; default is `RealEpsilon`.
            let prec = f64::EPSILON;
            let cur_tol = an_edge_tol + prec;
            if cur_tol >= error {
                maxdist = maxdist.max(an_edge_tol);
                is_same_p = true;
            }
        }
    }
    g.first = f3d;
    g.last = l3d;
    g.same_range = true;
    let mut new_tol = -1.0;
    if is_same_p {
        if ya_pcu {
            maxdist = maxdist.max(CONFUSION);
            new_tol = maxdist;
            g.tolerance = maxdist;
        }
        g.same_parameter = true;
    }
    reg.set_edge(&edge.0, g);
    new_tol
}

fn update_v_tol(edge: &Edge, new_tol: f64) {
    if new_tol <= 0.0 {
        return;
    }
    let reg = GeometryRegistry::global();
    let (v1, v2) = edge_vertices(edge);
    for v in [v1, v2].into_iter().flatten() {
        let mut vg = reg.vertex_geom(&v.0).unwrap_or(VertexGeom {
            point: Default::default(),
            tolerance: 0.0,
        });
        if new_tol > vg.tolerance {
            vg.tolerance = new_tol;
            reg.set_vertex(&v.0, vg);
        }
    }
}

fn update_tolerances(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    for face in faces_of(shape) {
        let ft = BRepTool::face_surface(&face)
            .map(|_| reg.face_tolerance(&face.0))
            .unwrap_or(0.0);
        for e in edges_of(&face.0) {
            if let Some(mut g) = reg.edge_geom(&e.0) {
                if g.tolerance < ft {
                    g.tolerance = ft;
                    reg.set_edge(&e.0, g);
                }
            }
            let et = reg.edge_geom(&e.0).map(|g| g.tolerance).unwrap_or(0.0);
            let (v1, v2) = edge_vertices(&e);
            for v in [v1, v2].into_iter().flatten() {
                let mut vg = reg.vertex_geom(&v.0).unwrap_or(VertexGeom {
                    point: Default::default(),
                    tolerance: 0.0,
                });
                if vg.tolerance < et {
                    vg.tolerance = et;
                    reg.set_vertex(&v.0, vg);
                }
            }
        }
    }
}

/// `BRepLib::SameParameter(S, Tolerance, forced)` (`cxx:1016`).
pub fn same_parameter_shape(shape: &TopoShape, tolerance: f64, forced: bool) {
    let reg = GeometryRegistry::global();
    let mut done: HashSet<usize> = HashSet::new();
    for e in edges_of(shape) {
        let key = GeometryRegistry::shape_key(&e.0);
        if !done.insert(key) {
            continue;
        }
        if let Some(mut g) = reg.edge_geom(&e.0) {
            if forced && (g.same_range || g.same_parameter) {
                g.same_range = false;
                g.same_parameter = false;
                reg.set_edge(&e.0, g);
            }
        }
        let faces = faces_of(shape);
        let face = faces.into_iter().find(|f| {
            edges_of(&f.0)
                .iter()
                .any(|fe| GeometryRegistry::shape_key(&fe.0) == key)
        });
        let Some(face) = face else {
            continue;
        };
        let new_tol = same_parameter_edge(&e, &face, tolerance);
        update_v_tol(&e, new_tol);
    }
    // Unported: planar `GetEdgeTol` (`cxx:973-1002`). Offset faces skip it.
    update_tolerances(shape);
}

/// `BRepLib::SameParameter(F, Confusion, true)` from `BRepLib_MakeFace.cxx:865`.
pub fn same_parameter_face(face: &Face, tolerance: f64, forced: bool) {
    same_parameter_shape(&face.0, tolerance, forced);
}
