//! `BRepTools::UVBounds` / `BRepTools::AddUVBounds` on the parametric plane of a
//! face.
//!
//! Source: `TKBRep/BRepTools/BRepTools.cxx:64-367` (OCCT 8.0.0): `UVBounds(F)`
//! `:64-80`, `UVBounds(F, E)` `:105-122`, `AddUVBounds(FF, B)` `:126-157`,
//! `AddUVBounds(aF, aE, aB)` `:172-367`.
//!
//! These are the *p-curve* bounds: the 2d box of the face's boundary p-curves,
//! clamped to the surface's parametric span. They are not the surface's natural
//! range ([`crate::brep_tool::BRepTool::uv_bounds_world`]), which keeps the
//! surface bounds of a trimmed face.
//!
//! `BndLib_Add2dCurve::Add(C, T1, T2, Tol, Box)` is
//! [`crate::bnd_lib_add2d::add_geom2d_range`]; the periodic verification of
//! `BRepTools::AddUVBounds` (`cxx:201-361`) is transcribed with the same point
//! probes and the same `100 * Precision::Confusion()^2` threshold.

use occt_core::bnd::BndBox2d;
use occt_geom2d::curve::Curve2d;

use crate::bnd_lib_add2d;
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topexp::Explorer;
use crate::abs::ShapeType;

/// `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)` (`cxx:64-81`).
///
/// When the box is void the four outputs are set to `0.0`, exactly as
/// `cxx:76-79`.
pub fn uv_bounds_face(f: &Face) -> (f64, f64, f64, f64) {
    let mut b = BndBox2d::new();
    add_uv_bounds_face(f, &mut b);
    match b.get() {
        Some((umin, vmin, umax, vmax)) => (umin, umax, vmin, vmax),
        None => (0.0, 0.0, 0.0, 0.0),
    }
}

/// `BRepTools::UVBounds(F, E, UMin, UMax, VMin, VMax)` (`cxx:105-122`).
pub fn uv_bounds_face_edge(f: &Face, e: &Edge) -> (f64, f64, f64, f64) {
    let mut b = BndBox2d::new();
    add_uv_bounds_face_edge(f, e, &mut b);
    match b.get() {
        Some((umin, vmin, umax, vmax)) => (umin, umax, vmin, vmax),
        None => (0.0, 0.0, 0.0, 0.0),
    }
}

/// `BRepTools::AddUVBounds(FF, B)` (`cxx:126-157`): the union over the face's
/// edges, falling back to the surface's natural bounds when the face carries no
/// p-curve at all (`cxx:139-154`).
pub fn add_uv_bounds_face(f: &Face, b: &mut BndBox2d) {
    let mut a_box = BndBox2d::new();
    let mut ex = Explorer::new(&f.0, ShapeType::Edge);
    while ex.more() {
        let e = ex.current().clone();
        ex.next();
        if let Some(edge) = Edge::wrap(e) {
            add_uv_bounds_face_edge(f, &edge, &mut a_box);
        }
    }
    if a_box.is_void() {
        let Some(s) = BRepTool::face_surface_world(f) else {
            return;
        };
        let (umin, umax) = s.u_range();
        let (vmin, vmax) = s.v_range();
        a_box.update(umin, vmin, umax, vmax);
    }
    b.add_box(&a_box);
}

/// `BRepTools::AddUVBounds(aF, aE, aB)` (`cxx:172-367`).
pub fn add_uv_bounds_face_edge(f: &Face, e: &Edge, b: &mut BndBox2d) {
    let Some((pc, t1, t2)) = boptools_2d::curve_on_surface_range(e, f) else {
        // `cxx:179-182`: a null `CurveOnSurface` returns without touching `aB`.
        return;
    };
    add_uv_bounds_pcurve(f, pc.as_ref(), t1, t2, b);
}

/// The part of `BRepTools::AddUVBounds` after `CurveOnSurface`
/// (`cxx:184-367`), taking the p-curve the caller already resolved.
///
/// `funks.cxx` needs this arm for an edge that is not on the face: it projects
/// the 3d curve first and then boxes the resulting p-curve.
pub fn add_uv_bounds_pcurve(
    f: &Face,
    pc: &dyn Curve2d,
    t1: f64,
    t2: f64,
    b: &mut BndBox2d,
) {
    let mut box_c = BndBox2d::new();
    bnd_lib_add2d::add_geom2d_range(pc, t1, t2, 0.0, &mut box_c);
    let (a_xmin, a_ymin, a_xmax, a_ymax) = match box_c.get() {
        Some(v) => v,
        // `cxx:187-190`: a void p-curve box leaves the four doubles at `0.0`.
        None => (0.0, 0.0, 0.0, 0.0),
    };

    let Some(surf) = BRepTool::face_surface_world(f) else {
        return;
    };
    let (a_umin, a_umax) = surf.u_range();
    let (a_vmin, a_vmax) = surf.v_range();

    // `cxx:194-199`: a `Geom_RectangularTrimmedSurface` is unwrapped for the
    // periodicity probes only; the bounds above stay those of the trimmed one.
    let basis = match surf.rectangular_trimmed_basis() {
        Some(base) => base,
        None => surf,
    };

    let (mut a_xmin, mut a_ymin, mut a_xmax, mut a_ymax) =
        (a_xmin, a_ymin, a_xmax, a_ymax);

    if !basis.is_u_periodic() {
        let mut is_u_periodic = false;
        // `cxx:208-211`
        if basis.is_bspline_surface() && (a_xmin < a_umin || a_xmax > a_umax) {
            let tol2 = 100.0 * occt_core::precision::CONFUSION * occt_core::precision::CONFUSION;
            is_u_periodic = true;
            // 1. verify that the surface is U-closed (`cxx:213-229`)
            if !basis.is_u_closed() {
                let a_v_step = a_vmax - a_vmin;
                let mut a_v = a_vmin;
                while a_v <= a_vmax {
                    let p1 = basis.value(a_umin, a_v);
                    let p2 = basis.value(a_umax, a_v);
                    if p1.square_distance(&p2) > tol2 {
                        is_u_periodic = false;
                        break;
                    }
                    if a_v_step <= 0.0 {
                        break;
                    }
                    a_v += a_v_step;
                }
            }
            // 2. verify periodicity of the surface inside the UV bounds of the
            // edge (`cxx:231-267`)
            if is_u_periodic {
                let a_v = (a_vmin + a_vmax) * 0.5;
                let mut a_u = [0.0f64; 6];
                let mut a_upp = [0.0f64; 6];
                let mut a_nb_pnt = 0;
                if a_xmin < a_umin {
                    a_u[0] = a_xmin;
                    a_u[1] = (a_xmin + a_umin) * 0.5;
                    a_u[2] = a_umin;
                    a_upp[0] = a_u[0] + a_umax - a_umin;
                    a_upp[1] = a_u[1] + a_umax - a_umin;
                    a_upp[2] = a_u[2] + a_umax - a_umin;
                    a_nb_pnt += 3;
                }
                if a_xmax > a_umax {
                    a_u[a_nb_pnt] = a_umax;
                    a_u[a_nb_pnt + 1] = (a_xmax + a_umax) * 0.5;
                    a_u[a_nb_pnt + 2] = a_xmax;
                    a_upp[a_nb_pnt] = a_u[a_nb_pnt] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 1] = a_u[a_nb_pnt + 1] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 2] = a_u[a_nb_pnt + 2] - a_umax + a_umin;
                    a_nb_pnt += 3;
                }
                for i in 0..a_nb_pnt {
                    let p1 = basis.value(a_u[i], a_v);
                    let p2 = basis.value(a_upp[i], a_v);
                    if p1.square_distance(&p2) > tol2 {
                        is_u_periodic = false;
                        break;
                    }
                }
            }
        }
        // `cxx:270-280`
        if !is_u_periodic {
            if a_xmin < a_umin && a_umin < a_xmax {
                a_xmin = a_umin;
            }
            if a_xmin < a_umax && a_umax < a_xmax {
                a_xmax = a_umax;
            }
        }
    }

    if !basis.is_v_periodic() {
        let mut is_v_periodic = false;
        // `cxx:290-292`
        if basis.is_bspline_surface() && (a_ymin < a_vmin || a_ymax > a_vmax) {
            let tol2 = 100.0 * occt_core::precision::CONFUSION * occt_core::precision::CONFUSION;
            is_v_periodic = true;
            // 1. verify that the surface is V-closed (`cxx:295-309`)
            if !basis.is_v_closed() {
                let a_u_step = a_umax - a_umin;
                let mut a_u = a_umin;
                while a_u <= a_umax {
                    let p1 = basis.value(a_u, a_vmin);
                    let p2 = basis.value(a_u, a_vmax);
                    if p1.square_distance(&p2) > tol2 {
                        is_v_periodic = false;
                        break;
                    }
                    if a_u_step <= 0.0 {
                        break;
                    }
                    a_u += a_u_step;
                }
            }
            // 2. verify periodicity of the surface inside the UV bounds of the
            // edge (`cxx:311-348`)
            if is_v_periodic {
                let a_u = (a_umin + a_umax) * 0.5;
                let mut a_v = [0.0f64; 6];
                let mut a_vpp = [0.0f64; 6];
                let mut a_nb_pnt = 0;
                if a_ymin < a_vmin {
                    a_v[0] = a_ymin;
                    a_v[1] = (a_ymin + a_vmin) * 0.5;
                    a_v[2] = a_vmin;
                    a_vpp[0] = a_v[0] + a_vmax - a_vmin;
                    a_vpp[1] = a_v[1] + a_vmax - a_vmin;
                    a_vpp[2] = a_v[2] + a_vmax - a_vmin;
                    a_nb_pnt += 3;
                }
                if a_ymax > a_vmax {
                    a_v[a_nb_pnt] = a_vmax;
                    a_v[a_nb_pnt + 1] = (a_ymax + a_vmax) * 0.5;
                    a_v[a_nb_pnt + 2] = a_ymax;
                    a_vpp[a_nb_pnt] = a_v[a_nb_pnt] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 1] = a_v[a_nb_pnt + 1] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 2] = a_v[a_nb_pnt + 2] - a_vmax + a_vmin;
                    a_nb_pnt += 3;
                }
                for i in 0..a_nb_pnt {
                    let p1 = basis.value(a_u, a_v[i]);
                    let p2 = basis.value(a_u, a_vpp[i]);
                    if p1.square_distance(&p2) > tol2 {
                        is_v_periodic = false;
                        break;
                    }
                }
            }
        }
        // `cxx:350-362`
        if !is_v_periodic {
            if a_ymin < a_vmin && a_vmin < a_ymax {
                a_ymin = a_vmin;
            }
            if a_ymin < a_vmax && a_vmax < a_ymax {
                a_ymax = a_vmax;
            }
        }
    }

    // `cxx:364-366`
    let mut a_box_s = BndBox2d::new();
    a_box_s.update(a_xmin, a_ymin, a_xmax, a_ymax);
    b.add_box(&a_box_s);
}
