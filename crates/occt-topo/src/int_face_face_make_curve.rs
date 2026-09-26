//! Walking-line `MakeCurve` path. Source: `IntTools_FaceFace.cxx:695` (Walking
//! branch), `GeomInt_IntSS::MakeBSpline` / `MakeBSpline2d`,
//! `ComputeTolReached3d`, `PrepareLines3D`, `ClassifyLin2d`, `ApproxWithPCurves`.
//!
//! Walking approximation goes through [`crate::geom_int::WlApprox`]. When
//! `IsDone` is false the same OCCT fallback applies: degree-1 B-splines
//! through the walking points.

use std::sync::Arc;

use occt_core::gp::{GpCylinder, GpDir2d, GpLin, GpLin2d, GpPnt, GpPnt2d, GpSphere};
use occt_core::precision::{CONFUSION, Precision};
use occt_geom::{Curve, GeomBSplineCurve, Surface};
use occt_geom::geom_api::project_point_on_surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::geom_int::{GeomIntLine, LineConstructor, ParametrizationType, TopolTool, WlApprox};
use super::bounds::parameter_out_of_boundary;
use crate::int_tools_lines::{is_closed_curve, reject_lines, split_curve};
use crate::int_tools_segpln::compute_tolerance;
use crate::int_tools_wline::{
    decomposition_of_wline, not_use_surfaces_for_approx_wline, WLine,
};
use crate::inttools_data::{CurveKind};
use crate::shape::Face;

use super::{curve_range, cylinder_from_surface, sphere_from_surface, FaceFace, FaceFaceCurve};

/// Flattened knot vector for OCCT `Geom_BSplineCurve(poles, knots, mults, 1)`
/// with unique knots `0..n-1` and end multiplicities 2.
fn degree1_polyline_knots(n: usize) -> Vec<f64> {
    let mut k = Vec::with_capacity(n + 2);
    k.push(0.0);
    for i in 0..n {
        k.push(i as f64);
    }
    k.push((n - 1) as f64);
    k
}

/// `GeomInt_IntSS::MakeBSpline` (`GeomInt_IntSS_1.cxx:1452`).
pub fn make_bspline(wl: &WLine, ideb: i32, ifin: i32) -> Option<Arc<dyn Curve>> {
    if ifin < ideb {
        return None;
    }
    let nbpnt = (ifin - ideb + 1) as usize;
    if nbpnt < 2 {
        return None;
    }
    let mut poles = Vec::with_capacity(nbpnt);
    for i in 0..nbpnt {
        poles.push(wl.point(ideb + i as i32).value());
    }
    let knots = degree1_polyline_knots(nbpnt);
    GeomBSplineCurve::new(poles, knots, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve>)
}

/// `GeomInt_IntSS::MakeBSpline2d` (`GeomInt_IntSS_1.cxx:1473`).
pub fn make_bspline2d(wl: &WLine, ideb: i32, ifin: i32, on_first: bool) -> Option<Arc<dyn Curve2d>> {
    if ifin < ideb {
        return None;
    }
    let nbpnt = (ifin - ideb + 1) as usize;
    if nbpnt < 2 {
        return None;
    }
    let mut xs = Vec::with_capacity(nbpnt);
    let mut ys = Vec::with_capacity(nbpnt);
    for i in 0..nbpnt {
        let p = wl.point(ideb + i as i32);
        let (u, v) = if on_first {
            p.parameters_on_s1()
        } else {
            p.parameters_on_s2()
        };
        xs.push(u);
        ys.push(v);
    }
    let knots = degree1_polyline_knots(nbpnt);
    Geom2dBSplineCurve::new(xs, ys, knots, 1)
        .ok()
        .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
}

/// `ApproxWithPCurves` (`IntTools_FaceFace.cxx:2357`).
pub fn approx_with_pcurves(the_cyl: &GpCylinder, the_sph: &GpSphere) -> bool {
    let b_res = true;
    let r1 = the_cyl.radius();
    let r2 = the_sph.radius();
    {
        let a_eps = 1.0e-7;
        let a_rc2 = r1 * r1;
        let a_ax3_sph = the_sph.position();
        let a_loc_sph = a_ax3_sph.location();
        let a_dir_sph = a_ax3_sph.direction();
        let a_ax1_cyl = the_cyl.axis();
        let a_lin_cyl = GpLin::new(a_ax1_cyl);
        let mut a_apex = GpPnt::new(
            a_loc_sph.x() + r2 * a_dir_sph.x(),
            a_loc_sph.y() + r2 * a_dir_sph.y(),
            a_loc_sph.z() + r2 * a_dir_sph.z(),
        );
        let a_d2 = a_lin_cyl.square_distance(&a_apex);
        if (a_d2 - a_rc2).abs() < a_eps {
            return !b_res;
        }
        a_apex = GpPnt::new(
            a_loc_sph.x() - r2 * a_dir_sph.x(),
            a_loc_sph.y() - r2 * a_dir_sph.y(),
            a_loc_sph.z() - r2 * a_dir_sph.z(),
        );
        let a_d2 = a_lin_cyl.square_distance(&a_apex);
        if (a_d2 - a_rc2).abs() < a_eps {
            return !b_res;
        }
    }
    if r1 < 2.0 * r2 {
        return b_res;
    }
    let an_cyl_ax = GpLin::new(the_cyl.axis());
    let a_dist = an_cyl_ax.distance(&the_sph.location());
    let a_d_rel = (a_dist - r1).abs() / r2;
    if a_d_rel > 0.2 {
        return b_res;
    }
    let loc = an_cyl_ax.location();
    let dir = an_cyl_ax.direction();
    let par = (the_sph.location().x() - loc.x()) * dir.x()
        + (the_sph.location().y() - loc.y()) * dir.y()
        + (the_sph.location().z() - loc.z()) * dir.z();
    let a_p = GpPnt::new(
        loc.x() + par * dir.x(),
        loc.y() + par * dir.y(),
        loc.z() + par * dir.z(),
    );
    let a_v = occt_core::gp::GpVec::from_pnts(&a_p, &the_sph.location());
    let dd = a_v.dot(&occt_core::gp::GpVec::from_xyz(the_sph.position().x_direction().xyz()));
    if a_dist < r1 && dd > 0.0 {
        return false;
    }
    if a_dist > r1 && dd < 0.0 {
        return false;
    }
    b_res
}

fn lin2d_coefficients(lin: &GpLin2d) -> (f64, f64, f64) {
    let d = lin.direction();
    let loc = lin.location();
    let a = d.y();
    let b = -d.x();
    let c = -(a * loc.x() + b * loc.y());
    (a, b, c)
}

fn lin2d_parameter(lin: &GpLin2d, p: &GpPnt2d) -> f64 {
    let loc = lin.location();
    let d = lin.direction();
    (p.x() - loc.x()) * d.x() + (p.y() - loc.y()) * d.y()
}

fn inter(d1: f64, d2: f64, tol: f64) -> bool {
    (d1 > tol && d2 < -tol)
        || (d1 < -tol && d2 > tol)
        || ((d1 <= tol && d1 >= -tol) && (d2 > tol || d2 < -tol))
        || ((d2 <= tol && d2 >= -tol) && (d1 > tol || d1 < -tol))
}

fn coinc(d1: f64, d2: f64, tol: f64) -> bool {
    (d1 <= tol && d1 >= -tol) && (d2 <= tol && d2 >= -tol)
}

/// `ClassifyLin2d` (`IntTools_FaceFace.cxx:2574`).
pub fn classify_lin2d(
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    the_lin2d: &GpLin2d,
    the_tol: f64,
) -> Option<(f64, f64)> {
    let (a, b, c) = lin2d_coefficients(the_lin2d);
    let mut par = [0.0f64; 2];
    let mut nbi = 0usize;
    let xmin = umin;
    let xmax = umax;
    let ymin = vmin;
    let ymax = vmax;

    let mut d1 = a * xmin + b * ymin + c;
    let mut d2 = a * xmin + b * ymax + c;
    if inter(d1, d2, the_tol) {
        let y = -(c + a * xmin) / b;
        par[nbi] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmin, y));
        nbi += 1;
    } else if coinc(d1, d2, the_tol) {
        par[0] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmin, ymin));
        par[1] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmin, ymax));
        nbi = 2;
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > the_tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmax + b * ymax + c;
    if d1 > the_tol || d1 < -the_tol {
        if inter(d1, d2, the_tol) {
            let x = -(c + b * ymax) / a;
            par[nbi] = lin2d_parameter(the_lin2d, &GpPnt2d::new(x, ymax));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmin, ymax));
            par[1] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmax, ymax));
            nbi = 2;
        }
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > the_tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmax + b * ymin + c;
    if d1 > the_tol || d1 < -the_tol {
        if inter(d1, d2, the_tol) {
            let y = -(c + a * xmax) / b;
            par[nbi] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmax, y));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmax, ymax));
            par[1] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmax, ymin));
            nbi = 2;
        }
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > the_tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmin + b * ymin + c;
    if d1 > the_tol || d1 < -the_tol {
        if inter(d1, d2, the_tol) {
            let x = -(c + b * ymin) / a;
            par[nbi] = lin2d_parameter(the_lin2d, &GpPnt2d::new(x, ymin));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmax, ymin));
            par[1] = lin2d_parameter(the_lin2d, &GpPnt2d::new(xmin, ymin));
            nbi = 2;
        }
    }
    if nbi == 2 {
        if (par[0] - par[1]).abs() > the_tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        }
    } else {
        None
    }
}

/// `ApproxParameters` (`IntTools_FaceFace.cxx:2736`).
pub fn approx_parameters(kind_a: SurfaceKind, kind_b: SurfaceKind) -> (i32, i32, i32) {
    let mut i_nb_iter = 0;
    let i_deg_min = 4;
    let mut i_deg_max = 8;
    if (kind_a == SurfaceKind::Cylinder && kind_b == SurfaceKind::Torus)
        || (kind_b == SurfaceKind::Cylinder && kind_a == SurfaceKind::Torus)
    {
        i_deg_max = 6;
    }
    if kind_a == SurfaceKind::Cylinder && kind_b == SurfaceKind::Cylinder {
        i_nb_iter = 1;
    }
    (i_deg_min, i_deg_max, i_nb_iter)
}

/// `Tolerances` (`IntTools_FaceFace.cxx:2787`).
pub fn tolerances(kind_a: SurfaceKind, kind_b: SurfaceKind, a_tol_tang: &mut f64) {
    if (kind_a == SurfaceKind::Cylinder && kind_b == SurfaceKind::Torus)
        || (kind_b == SurfaceKind::Cylinder && kind_a == SurfaceKind::Torus)
    {
        *a_tol_tang *= 0.1;
    }
}

/// Golden-section max of curve-to-surface distance on `[first, last]`.
fn find_max_distance_interval(
    the_c: &dyn Curve,
    the_first: f64,
    the_last: f64,
    surf: &dyn Surface,
    the_eps: f64,
) -> f64 {
    let a_cf = 0.6180339887498948;
    let mut a_a = the_first;
    let mut a_b = the_last;
    let mut a_x1 = a_b - a_cf * (a_b - a_a);
    let mut a_f1 = max_distance(the_c, a_x1, surf);
    let mut a_x2 = a_a + a_cf * (a_b - a_a);
    let mut a_f2 = max_distance(the_c, a_x2, surf);
    while (a_x1 - a_x2).abs() > the_eps {
        if a_f1 > a_f2 {
            a_b = a_x2;
            a_x2 = a_x1;
            a_f2 = a_f1;
            a_x1 = a_b - a_cf * (a_b - a_a);
            a_f1 = max_distance(the_c, a_x1, surf);
        } else {
            a_a = a_x1;
            a_x1 = a_x2;
            a_f1 = a_f2;
            a_x2 = a_a + a_cf * (a_b - a_a);
            a_f2 = max_distance(the_c, a_x2, surf);
        }
    }
    let mut a_f = max_distance(the_c, 0.5 * (a_x1 + a_x2), surf);
    if a_f1 > a_f {
        a_f = a_f1;
    }
    if a_f2 > a_f {
        a_f = a_f2;
    }
    a_f
}

fn max_distance(the_c: &dyn Curve, a_t: f64, surf: &dyn Surface) -> f64 {
    let a_p = the_c.d0(a_t);
    project_point_on_surface(surf, &a_p, 0.0)
        .map(|p| p.distance)
        .unwrap_or(0.0)
}

/// `FindMaxDistance` over a face (`IntTools_FaceFace.cxx:2899`).
pub fn find_max_distance_face(
    the_curve: &dyn Curve,
    the_first: f64,
    the_last: f64,
    the_face: &Face,
) -> f64 {
    let Some(surf) = BRepTool::face_surface(the_face) else {
        return 0.0;
    };
    let a_nb_s = 11.0;
    let a_dt = (the_last - the_first) / a_nb_s;
    let mut a_d_max = 0.0;
    let an_eps = 1.0e-4 * a_dt;
    let mut a_t2 = the_first;
    loop {
        let a_t1 = a_t2;
        a_t2 += a_dt;
        if a_t2 > the_last {
            break;
        }
        let a_d = find_max_distance_interval(the_curve, a_t1, a_t2, surf.as_ref(), an_eps);
        if a_d > a_d_max {
            a_d_max = a_d;
        }
    }
    a_d_max
}

/// `CheckPCurve` (`IntTools_FaceFace.cxx:3010`) with a single CN interval.
pub fn check_pcurve(a_pc: &dyn Curve2d, a_face: &Face) -> bool {
    const N_POINTS: i32 = 23;
    let (mut umin, mut umax, mut vmin, mut vmax) = crate::brep_surface::face_uv_bounds(a_face);
    let tol_u = ((umax - umin) * 0.01).max(CONFUSION);
    let tol_v = ((vmax - vmin) * 0.01).max(CONFUSION);
    let fp = a_pc.first_parameter();
    let lp = a_pc.last_parameter();
    if let Some(a_surf) = BRepTool::face_surface(a_face) {
        let pnt = a_pc.d0(0.5 * (fp + lp));
        let (u, v) = (pnt.x(), pnt.y());
        if a_surf.is_u_periodic() {
            let a_per = crate::int_tools_wline::u_period(a_surf.as_ref()).unwrap_or(0.0);
            if a_per > 0.0 {
                let mut nshift = ((u - umin) / a_per) as i32;
                if u < umin + a_per * nshift as f64 {
                    nshift -= 1;
                }
                umin += a_per * nshift as f64;
                umax += a_per * nshift as f64;
            }
        }
        if a_surf.is_v_periodic() {
            let a_per = crate::int_tools_wline::v_period(a_surf.as_ref()).unwrap_or(0.0);
            if a_per > 0.0 {
                let mut nshift = ((v - vmin) / a_per) as i32;
                if v < vmin + a_per * nshift as f64 {
                    nshift -= 1;
                }
                vmin += a_per * nshift as f64;
                vmax += a_per * nshift as f64;
            }
        }
    }
    let d_t = (lp - fp) / N_POINTS as f64;
    let mut a_t = fp;
    for _i in 1..N_POINTS {
        a_t += d_t;
        let a_p2d = a_pc.d0(a_t);
        let (u, v) = (a_p2d.x(), a_p2d.y());
        if umin - u > tol_u || u - umax > tol_u || vmin - v > tol_v || v - vmax > tol_v {
            return false;
        }
    }
    true
}

fn wline_to_curve(
    ff: &FaceFace,
    wl: &WLine,
    ifprm: i32,
    ilprm: i32,
    _fa: &Face,
    _fb: &Face,
) -> Option<FaceFaceCurve> {
    let curve = make_bspline(wl, ifprm, ilprm)?;
    let pcurve1 = if ff.approx1 {
        make_bspline2d(wl, ifprm, ilprm, true)
    } else {
        None
    };
    let pcurve2 = if ff.approx2 {
        make_bspline2d(wl, ifprm, ilprm, false)
    } else {
        None
    };
    let range = curve_range(curve.as_ref());
    Some(FaceFaceCurve {
        kind: CurveKind::BSpline,
        curve,
        range,
        face1_idx: 0,
        face2_idx: 1,
        pcurve1,
        pcurve2,
        tolerance: 0.0,
        tangential_tolerance: 0.0,
    })
}

fn line_constructor_parts(lc: &LineConstructor) -> Vec<(i32, i32)> {
    let mut out = Vec::new();
    for i in 1..=lc.nb_parts() {
        if let Some((fprm, lprm)) = lc.part(i) {
            out.push((fprm as i32, lprm as i32));
        }
    }
    out
}

fn curve_from_approx(mbspc: &crate::geom_int::MultiBSpCurve) -> Option<FaceFaceCurve> {
    let curve = mbspc.curve3d.clone()?;
    let range = curve_range(curve.as_ref());
    Some(FaceFaceCurve {
        kind: CurveKind::BSpline,
        curve,
        range,
        face1_idx: 0,
        face2_idx: 1,
        pcurve1: mbspc.curve2d_s1.clone(),
        pcurve2: mbspc.curve2d_s2.clone(),
        tolerance: 0.0,
        tangential_tolerance: 0.0,
    })
}

impl FaceFace {
    /// Walking-line `MakeCurve` (`IntTools_FaceFace.cxx` Walking case).
    ///
    /// LineConstructor splits the WLine into in-domain parts, then
    /// `DecompositionOfWLine` cuts periodic-boundary crossings. Each piece is
    /// approximated with `GeomInt_WLApprox`; if that is not done, OCCT's
    /// `MakeBSpline` fallback is used.
    pub(crate) fn make_curve_walking(
        &self,
        wl: &WLine,
        sa: &dyn Surface,
        sb: &dyn Surface,
        fa: &Face,
        fb: &Face,
    ) -> Vec<FaceFaceCurve> {
        if wl.nb_pnts() < 2 {
            return Vec::new();
        }

        let sa_arc: Arc<dyn Surface> = Arc::from(sa.clone_dyn());
        let sb_arc: Arc<dyn Surface> = Arc::from(sb.clone_dyn());
        let mut lc = LineConstructor::new();
        lc.load(
            TopolTool::from_face(fa, sa_arc.as_ref()),
            TopolTool::from_face(fb, sb_arc.as_ref()),
            sa_arc,
            sb_arc,
        );
        lc.perform(&GeomIntLine::Walking(wl.clone()));
        if !lc.is_done() || lc.nb_parts() <= 0 {
            return Vec::new();
        }
        let line_parts = line_constructor_parts(&lc);
        let b_avoid = false;

        if !self.approx {
            let mut out = Vec::new();
            for &(ifprm, ilprm) in &line_parts {
                if let Some(c) = wline_to_curve(self, wl, ifprm, ilprm, fa, fb) {
                    out.push(c);
                }
            }
            return out;
        }

        let typs1 = classify_surface(sa);
        let typs2 = classify_surface(sb);
        let mut an_approx1 = self.approx1;
        let mut an_approx2 = self.approx2;
        let mut my_tol_approx = self.tol_approx;
        let mut an_with_pc = true;
        if typs1 == SurfaceKind::Cylinder && typs2 == SurfaceKind::Sphere {
            if let (Some(cyl), Some(sph)) = (cylinder_from_surface(sa), sphere_from_surface(sb)) {
                an_with_pc = approx_with_pcurves(&cyl, &sph);
            }
        } else if typs1 == SurfaceKind::Sphere && typs2 == SurfaceKind::Cylinder {
            if let (Some(sph), Some(cyl)) = (sphere_from_surface(sa), cylinder_from_surface(sb)) {
                an_with_pc = approx_with_pcurves(&cyl, &sph);
            }
        }
        if !an_with_pc {
            my_tol_approx = 1.0e-5;
            an_approx1 = false;
            an_approx2 = false;
        }

        let mut a_seq_of_l: Vec<WLine> = Vec::new();
        let b_is_decomposited = decomposition_of_wline(
            wl,
            sa,
            sb,
            fa,
            fb,
            &line_parts,
            b_avoid,
            self.tol,
            &mut a_seq_of_l,
        );

        let nbiter = if b_is_decomposited {
            a_seq_of_l.len() as i32
        } else {
            line_parts.len() as i32
        };

        let (i_deg_min, i_deg_max, i_nb_iter) = approx_parameters(typs1, typs2);
        let mut _tol_tang = self.tol;
        tolerances(typs1, typs2, &mut _tol_tang);

        let mut out = Vec::new();
        for i in 1..=nbiter {
            let (piece, ifprm, ilprm): (&WLine, i32, i32) = if b_is_decomposited {
                let w = &a_seq_of_l[(i as usize).saturating_sub(1)];
                (w, 1, w.nb_pnts())
            } else {
                let (f, l) = line_parts[(i as usize).saturating_sub(1)];
                (wl, f, l)
            };
            let mut an_approx = self.approx;
            let mut a1 = an_approx1;
            let mut a2 = an_approx2;
            if typs1 == SurfaceKind::Plane {
                an_approx = false;
                a1 = true;
            } else if typs2 == SurfaceKind::Plane {
                an_approx = false;
                a2 = true;
            }
            let b_reject = not_use_surfaces_for_approx_wline(fa, fb, piece, ifprm, ilprm);
            let mut theapp3d = WlApprox::new();
            theapp3d.set_parameters(
                my_tol_approx,
                my_tol_approx,
                i_deg_min,
                i_deg_max,
                i_nb_iter,
                30,
                !b_reject,
                ParametrizationType::ChordLength,
            );
            if typs1 == SurfaceKind::Plane {
                theapp3d.perform(sa, sb, piece, false, true, a2, ifprm, ilprm);
            } else if typs2 == SurfaceKind::Plane {
                theapp3d.perform(sa, sb, piece, false, a1, true, ifprm, ilprm);
            } else {
                theapp3d.perform(sa, sb, piece, an_approx, a1, a2, ifprm, ilprm);
            }
            let mut pushed = false;
            if theapp3d.is_done() {
                for j in 1..=theapp3d.nb_multi_curves() {
                    if let Some(mbspc) = theapp3d.value(j) {
                        if let Some(c) = curve_from_approx(mbspc) {
                            out.push(c);
                            pushed = true;
                        }
                    }
                }
            }
            if !pushed {
                if let Some(c) = wline_to_curve(self, piece, ifprm, ilprm, fa, fb) {
                    out.push(c);
                }
            }
        }
        out
    }

    /// `ComputeTolReached3d` (`IntTools_FaceFace.cxx:613`).
    pub(crate) fn compute_tol_reached_3d(&mut self, fa: &Face, fb: &Face) {
        let a_tol_f_max = BRepTool::face_tolerance(fa).max(BRepTool::face_tolerance(fb));
        let Some(a_s1) = BRepTool::face_surface(fa) else {
            return;
        };
        let Some(a_s2) = BRepTool::face_surface(fb) else {
            return;
        };
        for c in &mut self.result.curves {
            let a_first = c.curve.first_parameter();
            let a_last = c.curve.last_parameter();
            let mut a_tol_c = c.tolerance;
            for j in 0..2 {
                let a_c2d = if j == 0 { &c.pcurve1 } else { &c.pcurve2 };
                if let Some(pc) = a_c2d {
                    let a_s = if j == 0 { a_s1.as_ref() } else { a_s2.as_ref() };
                    if let Some((a_d, _a_t)) =
                        compute_tolerance(c.curve.as_ref(), pc.as_ref(), a_s, a_first, a_last)
                    {
                        if a_d > a_tol_c {
                            a_tol_c = a_d;
                        }
                    }
                } else {
                    let a_f = if j == 0 { fa } else { fb };
                    let a_d = find_max_distance_face(c.curve.as_ref(), a_first, a_last, a_f);
                    if a_d > a_tol_c {
                        a_tol_c = a_d;
                    }
                }
            }
            c.tolerance = a_tol_c;
            if c.tangential_tolerance < a_tol_f_max {
                c.tangential_tolerance = a_tol_f_max;
            }
        }
    }

    /// `PrepareLines3D` (`IntTools_FaceFace.cxx:1932`).
    pub fn prepare_lines_3d(&mut self, b_to_split: bool) {
        let mut a_new_cvs: Vec<FaceFaceCurve> = Vec::new();
        for c in &self.result.curves {
            if b_to_split {
                let split = split_curve(c);
                if split.is_empty() {
                    a_new_cvs.push(c.clone());
                } else {
                    a_new_cvs.extend(split);
                }
            } else {
                a_new_cvs.push(c.clone());
            }
        }
        let a_type1 = self.face1.as_ref().and_then(BRepTool::face_surface);
        let a_type2 = self.face2.as_ref().and_then(BRepTool::face_surface);
        let k1 = a_type1.as_ref().map(|s| classify_surface(s.as_ref()));
        let k2 = a_type2.as_ref().map(|s| classify_surface(s.as_ref()));
        let plane_cone = matches!(
            (k1, k2),
            (Some(SurfaceKind::Plane), Some(SurfaceKind::Cone))
                | (Some(SurfaceKind::Cone), Some(SurfaceKind::Plane))
        );
        if plane_cone && a_new_cvs.len() == 4 && a_new_cvs[0].kind == CurveKind::Line {
            a_new_cvs = reject_lines(&a_new_cvs);
        }
        self.result.curves = a_new_cvs;
    }
}

/// Used by `ParameterOutOfBoundary` callers in the circle/ellipse MakeCurve path.
#[allow(dead_code)]
pub fn adjust_closed_part(
    newc: &dyn Curve,
    fprm: f64,
    lprm: f64,
    fa: &Face,
    fb: &Face,
    my_tol: f64,
) -> (f64, f64) {
    let mut f = fprm;
    let mut l = lprm;
    if let Some(anew) = parameter_out_of_boundary(fprm, newc, fa, fb, lprm, false, my_tol) {
        f = anew;
    }
    if let Some(anew) = parameter_out_of_boundary(lprm, newc, fa, fb, fprm, true, my_tol) {
        l = anew;
    }
    (f, l)
}

#[allow(dead_code)]
fn _is_closed_hint(c: &dyn Curve) -> bool {
    is_closed_curve(c)
}

#[allow(dead_code)]
fn _dir2d_x() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).unwrap_or_default()
}

#[allow(dead_code)]
fn _precision_infinite(x: f64) -> bool {
    Precision::is_infinite(x)
}
