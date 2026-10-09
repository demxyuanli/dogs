//! UV bounding boxes of faces, wires and edges.
//!
//! Source: `BRepTools.cxx` `AddUVBounds` (face at 126, wire at 161, edge at
//! 172). PerformAreas builds a `BOPTools_Box2dTree` of hole-face UV boxes
//! with `BRepTools::AddUVBounds(aHFace, aBox)` and queries it with the
//! growth-face UV box.
//!
//! Edge overload:
//! 1. `BRep_Tool::CurveOnSurface` — missing pcurve → return;
//! 2. `BndLib_Add2dCurve::Add(C2D, T1, T2, 0., BoxC)` then `Get` the UV range;
//! 3. `Bounds()` of the STORED surface first (`BRepTools.cxx:191-192`), THEN
//!    peel `Geom_RectangularTrimmedSurface` down to its basis surface, which the
//!    reference uses only for the periodicity/closedness tests
//!    (`BRepTools.cxx:194-199`). The port keeps the stored surface for both
//!    steps: the ranges still agree (`Geom_RectangularTrimmedSurface::Bounds`
//!    returns the trim, `Geom_RectangularTrimmedSurface.cxx:433-440`, exactly as
//!    our `u_range`/`v_range` do), but a trimmed wrapper answers
//!    `IsUPeriodic() == false` for a periodic basis whose trim is not a whole
//!    number of periods (`Geom_RectangularTrimmedSurface.cxx:509-523`). OCCT
//!    tests the BASIS there (periodic -> no clamping), we test the wrapper and
//!    take the clamp below, so the two differ for such faces;
//! 4. if the surface is not U-periodic, clamp the U range to surface bounds
//!    unless a B-spline extra-periodicity check succeeds;
//! 5. same for V;
//! 6. `BoxS.Update(Xmin, Ymin, Xmax, Ymax)` then `B.Add(BoxS)`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::bnd_lib_add2d::add_geom2d_range;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Wire};
use crate::topo_tools_full::{edges_of, edges_of_wire};

/// `BRepTools::AddUVBounds(const TopoDS_Face& FF, Bnd_Box2d& B)`.
pub fn add_uv_bounds_face(ff: &Face, b: &mut BndBox2d) {
    let mut f = ff.clone();
    f.0.set_orientation(Orientation::Forward);
    let mut a_box = BndBox2d::new();
    for e in edges_of(&f.0) {
        add_uv_bounds_edge(&f, &e, &mut a_box);
    }
    if a_box.is_void() {
        let Some(a_surf) = BRepTool::face_surface(&f) else {
            return;
        };
        let (u_min, u_max) = a_surf.u_range();
        let (v_min, v_max) = a_surf.v_range();
        a_box.update(u_min, v_min, u_max, v_max);
    }
    b.add_box(&a_box);
}

/// `BRepTools::AddUVBounds(Face, Wire, Box)`.
pub fn add_uv_bounds_wire(f: &Face, w: &Wire, b: &mut BndBox2d) {
    for e in edges_of_wire(w) {
        add_uv_bounds_edge(f, &e, b);
    }
}

fn surface_is_bspline(a_s: &dyn Surface) -> bool {
    a_s.is_bspline_surface()
}

/// Extra U-periodicity verification for B-spline surfaces
/// (`BRepTools.cxx:202-268`).
///
/// 1. Verify the surface is U-closed (2 points) when `IsUClosed` is false.
/// 2. Verify periodicity inside the edge UV-bounds (3 or 6 points).
fn verify_u_periodic_bspline(
    a_s: &dyn Surface,
    a_xmin: f64,
    a_xmax: f64,
    a_umin: f64,
    a_umax: f64,
    a_vmin: f64,
    a_vmax: f64,
) -> bool {
    if !surface_is_bspline(a_s) {
        return false;
    }
    if !(a_xmin < a_umin || a_xmax > a_umax) {
        return false;
    }
    let a_tol2 = 100.0 * CONFUSION * CONFUSION;
    let mut is_u_periodic = true;
    let a_v_step = a_vmax - a_vmin;
    let mut a_v = a_vmin;
    while a_v <= a_vmax {
        let p1 = a_s.d0(a_umin, a_v);
        let p2 = a_s.d0(a_umax, a_v);
        if p1.square_distance(&p2) > a_tol2 {
            is_u_periodic = false;
            break;
        }
        if a_v_step.abs() <= f64::EPSILON {
            break;
        }
        a_v += a_v_step;
    }
    if is_u_periodic {
        let a_v = 0.5 * (a_vmin + a_vmax);
        let mut a_u = [0.0; 6];
        let mut a_upp = [0.0; 6];
        let mut a_nb_pnt = 0usize;
        if a_xmin < a_umin {
            a_u[0] = a_xmin;
            a_u[1] = 0.5 * (a_xmin + a_umin);
            a_u[2] = a_umin;
            a_upp[0] = a_u[0] + a_umax - a_umin;
            a_upp[1] = a_u[1] + a_umax - a_umin;
            a_upp[2] = a_u[2] + a_umax - a_umin;
            a_nb_pnt += 3;
        }
        if a_xmax > a_umax {
            a_u[a_nb_pnt] = a_umax;
            a_u[a_nb_pnt + 1] = 0.5 * (a_xmax + a_umax);
            a_u[a_nb_pnt + 2] = a_xmax;
            a_upp[a_nb_pnt] = a_u[a_nb_pnt] - a_umax + a_umin;
            a_upp[a_nb_pnt + 1] = a_u[a_nb_pnt + 1] - a_umax + a_umin;
            a_upp[a_nb_pnt + 2] = a_u[a_nb_pnt + 2] - a_umax + a_umin;
            a_nb_pnt += 3;
        }
        for an_ind in 0..a_nb_pnt {
            let p1 = a_s.d0(a_u[an_ind], a_v);
            let p2 = a_s.d0(a_upp[an_ind], a_v);
            if p1.square_distance(&p2) > a_tol2 {
                is_u_periodic = false;
                break;
            }
        }
    }
    is_u_periodic
}

/// Extra V-periodicity verification (`BRepTools.cxx:283-349`).
fn verify_v_periodic_bspline(
    a_s: &dyn Surface,
    a_ymin: f64,
    a_ymax: f64,
    a_umin: f64,
    a_umax: f64,
    a_vmin: f64,
    a_vmax: f64,
) -> bool {
    if !surface_is_bspline(a_s) {
        return false;
    }
    if !(a_ymin < a_vmin || a_ymax > a_vmax) {
        return false;
    }
    let a_tol2 = 100.0 * CONFUSION * CONFUSION;
    let mut is_v_periodic = true;
    let a_u_step = a_umax - a_umin;
    let mut a_u = a_umin;
    while a_u <= a_umax {
        let p1 = a_s.d0(a_u, a_vmin);
        let p2 = a_s.d0(a_u, a_vmax);
        if p1.square_distance(&p2) > a_tol2 {
            is_v_periodic = false;
            break;
        }
        if a_u_step.abs() <= f64::EPSILON {
            break;
        }
        a_u += a_u_step;
    }
    if is_v_periodic {
        let a_u = 0.5 * (a_umin + a_umax);
        let mut a_v = [0.0; 6];
        let mut a_vpp = [0.0; 6];
        let mut a_nb_pnt = 0usize;
        if a_ymin < a_vmin {
            a_v[0] = a_ymin;
            a_v[1] = 0.5 * (a_ymin + a_vmin);
            a_v[2] = a_vmin;
            a_vpp[0] = a_v[0] + a_vmax - a_vmin;
            a_vpp[1] = a_v[1] + a_vmax - a_vmin;
            a_vpp[2] = a_v[2] + a_vmax - a_vmin;
            a_nb_pnt += 3;
        }
        if a_ymax > a_vmax {
            a_v[a_nb_pnt] = a_vmax;
            a_v[a_nb_pnt + 1] = 0.5 * (a_ymax + a_vmax);
            a_v[a_nb_pnt + 2] = a_ymax;
            a_vpp[a_nb_pnt] = a_v[a_nb_pnt] - a_vmax + a_vmin;
            a_vpp[a_nb_pnt + 1] = a_v[a_nb_pnt + 1] - a_vmax + a_vmin;
            a_vpp[a_nb_pnt + 2] = a_v[a_nb_pnt + 2] - a_vmax + a_vmin;
            a_nb_pnt += 3;
        }
        for an_ind in 0..a_nb_pnt {
            let p1 = a_s.d0(a_u, a_v[an_ind]);
            let p2 = a_s.d0(a_u, a_vpp[an_ind]);
            if p1.square_distance(&p2) > a_tol2 {
                is_v_periodic = false;
                break;
            }
        }
    }
    is_v_periodic
}

/// `BRepTools::AddUVBounds(Face, Edge, Box)`.
pub fn add_uv_bounds_edge(a_f: &Face, a_e: &Edge, a_b: &mut BndBox2d) {
    let Some((a_c2d, a_t1, a_t2)) = curve_on_surface_range(a_e, a_f) else {
        return;
    };
    let mut a_box_c = BndBox2d::new();
    add_geom2d_range(a_c2d.as_ref(), a_t1, a_t2, 0.0, &mut a_box_c);
    let mut a_xmin = 0.0;
    let mut a_ymin = 0.0;
    let mut a_xmax = 0.0;
    let mut a_ymax = 0.0;
    if !a_box_c.is_void() {
        if let Some((x, y, xx, yy)) = a_box_c.get() {
            a_xmin = x;
            a_ymin = y;
            a_xmax = xx;
            a_ymax = yy;
        }
    }
    // `BRepTools.cxx:191-199`: the ranges come from the **stored** surface
    // (`Geom_RectangularTrimmedSurface::Bounds` returns the trim), but the
    // periodicity / closedness tests and every `Value` call are made on the
    // BASIS surface the trim wrapper is peeled down to. Testing the wrapper
    // itself answers `IsUPeriodic() == false` for a periodic basis trimmed to
    // less than a whole period (`Geom_RectangularTrimmedSurface.cxx:509-523`)
    // and fails the `BSplineSurface` type check, so the clamp below would fire
    // where OCCT skips it.
    let stored = BRepTool::face_surface(a_f);
    let Some(stored) = stored else {
        return;
    };
    let (a_umin, a_umax) = stored.u_range();
    let (a_vmin, a_vmax) = stored.v_range();
    let basis = stored.rectangular_trimmed_basis();
    let a_s: &dyn Surface = match basis.as_deref() {
        Some(b) => b,
        None => stored.as_ref(),
    };

    if !a_s.is_u_periodic() {
        let is_u_periodic = verify_u_periodic_bspline(
            a_s,
            a_xmin,
            a_xmax,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );
        if !is_u_periodic {
            if a_xmin < a_umin && a_umin < a_xmax {
                a_xmin = a_umin;
            }
            if a_xmin < a_umax && a_umax < a_xmax {
                a_xmax = a_umax;
            }
        }
    }

    if !a_s.is_v_periodic() {
        let is_v_periodic = verify_v_periodic_bspline(
            a_s,
            a_ymin,
            a_ymax,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );
        if !is_v_periodic {
            if a_ymin < a_vmin && a_vmin < a_ymax {
                a_ymin = a_vmin;
            }
            if a_ymin < a_vmax && a_vmax < a_ymax {
                a_ymax = a_vmax;
            }
        }
    }

    let mut a_box_s = BndBox2d::new();
    a_box_s.update(a_xmin, a_ymin, a_xmax, a_ymax);
    a_b.add_box(&a_box_s);
    let _ = CONFUSION;
}

/// Convenience: UV box of a face (empty → natural bounds).
pub fn uv_box_of_face(face: &Face) -> BndBox2d {
    let mut b = BndBox2d::new();
    add_uv_bounds_face(face, &mut b);
    b
}

/// Convenience: UV box of an edge on a face.
pub fn uv_box_of_edge_on_face(face: &Face, edge: &Edge) -> BndBox2d {
    let mut b = BndBox2d::new();
    add_uv_bounds_edge(face, edge, &mut b);
    b
}
