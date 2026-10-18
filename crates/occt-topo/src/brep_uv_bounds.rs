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
    let mut bnd_box = BndBox2d::new();
    for e in edges_of(&f.0) {
        add_uv_bounds_edge(&f, &e, &mut bnd_box);
    }
    if bnd_box.is_void() {
        let Some(surf) = BRepTool::face_surface(&f) else {
            return;
        };
        let (u_min, u_max) = surf.u_range();
        let (v_min, v_max) = surf.v_range();
        bnd_box.update(u_min, v_min, u_max, v_max);
    }
    b.add_box(&bnd_box);
}

/// `BRepTools::AddUVBounds(Face, Wire, Box)`.
pub fn add_uv_bounds_wire(f: &Face, w: &Wire, b: &mut BndBox2d) {
    for e in edges_of_wire(w) {
        add_uv_bounds_edge(f, &e, b);
    }
}

fn surface_is_bspline(s: &dyn Surface) -> bool {
    s.is_bspline_surface()
}

/// Extra U-periodicity verification for B-spline surfaces
/// (`BRepTools.cxx:202-268`).
///
/// 1. Verify the surface is U-closed (2 points) when `IsUClosed` is false.
/// 2. Verify periodicity inside the edge UV-bounds (3 or 6 points).
fn verify_u_periodic_bspline(
    s: &dyn Surface,
    xmin: f64,
    xmax: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> bool {
    if !surface_is_bspline(s) {
        return false;
    }
    if !(xmin < umin || xmax > umax) {
        return false;
    }
    let tol2 = 100.0 * CONFUSION * CONFUSION;
    let mut is_u_periodic = true;
    let v_step = vmax - vmin;
    let mut v = vmin;
    while v <= vmax {
        let p1 = s.d0(umin, v);
        let p2 = s.d0(umax, v);
        if p1.square_distance(&p2) > tol2 {
            is_u_periodic = false;
            break;
        }
        if v_step.abs() <= f64::EPSILON {
            break;
        }
        v += v_step;
    }
    if is_u_periodic {
        let v = 0.5 * (vmin + vmax);
        let mut u = [0.0; 6];
        let mut upp = [0.0; 6];
        let mut nb_pnt = 0usize;
        if xmin < umin {
            u[0] = xmin;
            u[1] = 0.5 * (xmin + umin);
            u[2] = umin;
            upp[0] = u[0] + umax - umin;
            upp[1] = u[1] + umax - umin;
            upp[2] = u[2] + umax - umin;
            nb_pnt += 3;
        }
        if xmax > umax {
            u[nb_pnt] = umax;
            u[nb_pnt + 1] = 0.5 * (xmax + umax);
            u[nb_pnt + 2] = xmax;
            upp[nb_pnt] = u[nb_pnt] - umax + umin;
            upp[nb_pnt + 1] = u[nb_pnt + 1] - umax + umin;
            upp[nb_pnt + 2] = u[nb_pnt + 2] - umax + umin;
            nb_pnt += 3;
        }
        for an_ind in 0..nb_pnt {
            let p1 = s.d0(u[an_ind], v);
            let p2 = s.d0(upp[an_ind], v);
            if p1.square_distance(&p2) > tol2 {
                is_u_periodic = false;
                break;
            }
        }
    }
    is_u_periodic
}

/// Extra V-periodicity verification (`BRepTools.cxx:283-349`).
fn verify_v_periodic_bspline(
    s: &dyn Surface,
    ymin: f64,
    ymax: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> bool {
    if !surface_is_bspline(s) {
        return false;
    }
    if !(ymin < vmin || ymax > vmax) {
        return false;
    }
    let tol2 = 100.0 * CONFUSION * CONFUSION;
    let mut is_v_periodic = true;
    let u_step = umax - umin;
    let mut u = umin;
    while u <= umax {
        let p1 = s.d0(u, vmin);
        let p2 = s.d0(u, vmax);
        if p1.square_distance(&p2) > tol2 {
            is_v_periodic = false;
            break;
        }
        if u_step.abs() <= f64::EPSILON {
            break;
        }
        u += u_step;
    }
    if is_v_periodic {
        let u = 0.5 * (umin + umax);
        let mut v = [0.0; 6];
        let mut vpp = [0.0; 6];
        let mut nb_pnt = 0usize;
        if ymin < vmin {
            v[0] = ymin;
            v[1] = 0.5 * (ymin + vmin);
            v[2] = vmin;
            vpp[0] = v[0] + vmax - vmin;
            vpp[1] = v[1] + vmax - vmin;
            vpp[2] = v[2] + vmax - vmin;
            nb_pnt += 3;
        }
        if ymax > vmax {
            v[nb_pnt] = vmax;
            v[nb_pnt + 1] = 0.5 * (ymax + vmax);
            v[nb_pnt + 2] = ymax;
            vpp[nb_pnt] = v[nb_pnt] - vmax + vmin;
            vpp[nb_pnt + 1] = v[nb_pnt + 1] - vmax + vmin;
            vpp[nb_pnt + 2] = v[nb_pnt + 2] - vmax + vmin;
            nb_pnt += 3;
        }
        for an_ind in 0..nb_pnt {
            let p1 = s.d0(u, v[an_ind]);
            let p2 = s.d0(u, vpp[an_ind]);
            if p1.square_distance(&p2) > tol2 {
                is_v_periodic = false;
                break;
            }
        }
    }
    is_v_periodic
}

/// `BRepTools::AddUVBounds(Face, Edge, Box)`.
pub fn add_uv_bounds_edge(f: &Face, e: &Edge, a_b: &mut BndBox2d) {
    let Some((c2d, t1, t2)) = curve_on_surface_range(e, f) else {
        return;
    };
    let mut box_c = BndBox2d::new();
    add_geom2d_range(c2d.as_ref(), t1, t2, 0.0, &mut box_c);
    let mut xmin = 0.0;
    let mut ymin = 0.0;
    let mut xmax = 0.0;
    let mut ymax = 0.0;
    if !box_c.is_void() {
        if let Some((x, y, xx, yy)) = box_c.get() {
            xmin = x;
            ymin = y;
            xmax = xx;
            ymax = yy;
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
    let stored = BRepTool::face_surface(f);
    let Some(stored) = stored else {
        return;
    };
    let (umin, umax) = stored.u_range();
    let (vmin, vmax) = stored.v_range();
    let basis = stored.rectangular_trimmed_basis();
    let s: &dyn Surface = match basis.as_deref() {
        Some(b) => b,
        None => stored.as_ref(),
    };

    if !s.is_u_periodic() {
        let is_u_periodic = verify_u_periodic_bspline(
            s,
            xmin,
            xmax,
            umin,
            umax,
            vmin,
            vmax,
        );
        if !is_u_periodic {
            if xmin < umin && umin < xmax {
                xmin = umin;
            }
            if xmin < umax && umax < xmax {
                xmax = umax;
            }
        }
    }

    if !s.is_v_periodic() {
        let is_v_periodic = verify_v_periodic_bspline(
            s,
            ymin,
            ymax,
            umin,
            umax,
            vmin,
            vmax,
        );
        if !is_v_periodic {
            if ymin < vmin && vmin < ymax {
                ymin = vmin;
            }
            if ymin < vmax && vmax < ymax {
                ymax = vmax;
            }
        }
    }

    let mut box_s = BndBox2d::new();
    box_s.update(xmin, ymin, xmax, ymax);
    a_b.add_box(&box_s);
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
