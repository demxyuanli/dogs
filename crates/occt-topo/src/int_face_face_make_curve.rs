//! Walking-line `MakeCurve` path. Source: `IntTools_FaceFace.cxx:695` (Walking
//! branch), `GeomInt_IntSS::MakeBSpline` / `MakeBSpline2d`,
//! `ComputeTolReached3d`, `PrepareLines3D`, `ClassifyLin2d`, `ApproxWithPCurves`.
//!
//! Walking approximation goes through [`crate::geom_int::WlApprox`]. When
//! `IsDone` is false the same OCCT fallback applies: degree-1 B-splines
//! through the walking points.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

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

use super::{curve_range, gp_cylinder_from_surface, gp_sphere_from_surface, FaceFace, FaceFaceCurve};

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
pub fn approx_with_pcurves(cyl: &GpCylinder, sph: &GpSphere) -> bool {
    let b_res = true;
    let r1 = cyl.radius();
    let r2 = sph.radius();
    {
        let eps = 1.0e-7;
        let rc2 = r1 * r1;
        let ax3_sph = sph.position();
        let loc_sph = ax3_sph.location();
        let dir_sph = ax3_sph.direction();
        let ax1_cyl = cyl.axis();
        let lin_cyl = GpLin::new(ax1_cyl);
        let mut apex = GpPnt::new(
            loc_sph.x() + r2 * dir_sph.x(),
            loc_sph.y() + r2 * dir_sph.y(),
            loc_sph.z() + r2 * dir_sph.z(),
        );
        let d2 = lin_cyl.square_distance(&apex);
        if (d2 - rc2).abs() < eps {
            return !b_res;
        }
        apex = GpPnt::new(
            loc_sph.x() - r2 * dir_sph.x(),
            loc_sph.y() - r2 * dir_sph.y(),
            loc_sph.z() - r2 * dir_sph.z(),
        );
        let d2 = lin_cyl.square_distance(&apex);
        if (d2 - rc2).abs() < eps {
            return !b_res;
        }
    }
    if r1 < 2.0 * r2 {
        return b_res;
    }
    let an_cyl_ax = GpLin::new(cyl.axis());
    let dist = an_cyl_ax.distance(&sph.location());
    let d_rel = (dist - r1).abs() / r2;
    if d_rel > 0.2 {
        return b_res;
    }
    let loc = an_cyl_ax.location();
    let dir = an_cyl_ax.direction();
    let par = (sph.location().x() - loc.x()) * dir.x()
        + (sph.location().y() - loc.y()) * dir.y()
        + (sph.location().z() - loc.z()) * dir.z();
    let p = GpPnt::new(
        loc.x() + par * dir.x(),
        loc.y() + par * dir.y(),
        loc.z() + par * dir.z(),
    );
    let v = occt_core::gp::GpVec::from_pnts(&p, &sph.location());
    let dd = v.dot(&occt_core::gp::GpVec::from_xyz(sph.position().x_direction().xyz()));
    if dist < r1 && dd > 0.0 {
        return false;
    }
    if dist > r1 && dd < 0.0 {
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
    lin2d: &GpLin2d,
    tol: f64,
) -> Option<(f64, f64)> {
    let (a, b, c) = lin2d_coefficients(lin2d);
    let mut par = [0.0f64; 2];
    let mut nbi = 0usize;
    let xmin = umin;
    let xmax = umax;
    let ymin = vmin;
    let ymax = vmax;

    let mut d1 = a * xmin + b * ymin + c;
    let mut d2 = a * xmin + b * ymax + c;
    if inter(d1, d2, tol) {
        let y = -(c + a * xmin) / b;
        par[nbi] = lin2d_parameter(lin2d, &GpPnt2d::new(xmin, y));
        nbi += 1;
    } else if coinc(d1, d2, tol) {
        par[0] = lin2d_parameter(lin2d, &GpPnt2d::new(xmin, ymin));
        par[1] = lin2d_parameter(lin2d, &GpPnt2d::new(xmin, ymax));
        nbi = 2;
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmax + b * ymax + c;
    if d1 > tol || d1 < -tol {
        if inter(d1, d2, tol) {
            let x = -(c + b * ymax) / a;
            par[nbi] = lin2d_parameter(lin2d, &GpPnt2d::new(x, ymax));
            nbi += 1;
        } else if coinc(d1, d2, tol) {
            par[0] = lin2d_parameter(lin2d, &GpPnt2d::new(xmin, ymax));
            par[1] = lin2d_parameter(lin2d, &GpPnt2d::new(xmax, ymax));
            nbi = 2;
        }
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmax + b * ymin + c;
    if d1 > tol || d1 < -tol {
        if inter(d1, d2, tol) {
            let y = -(c + a * xmax) / b;
            par[nbi] = lin2d_parameter(lin2d, &GpPnt2d::new(xmax, y));
            nbi += 1;
        } else if coinc(d1, d2, tol) {
            par[0] = lin2d_parameter(lin2d, &GpPnt2d::new(xmax, ymax));
            par[1] = lin2d_parameter(lin2d, &GpPnt2d::new(xmax, ymin));
            nbi = 2;
        }
    }
    if nbi == 2 {
        return if (par[0] - par[1]).abs() > tol {
            Some((par[0].min(par[1]), par[0].max(par[1])))
        } else {
            None
        };
    }

    d1 = d2;
    d2 = a * xmin + b * ymin + c;
    if d1 > tol || d1 < -tol {
        if inter(d1, d2, tol) {
            let x = -(c + b * ymin) / a;
            par[nbi] = lin2d_parameter(lin2d, &GpPnt2d::new(x, ymin));
            nbi += 1;
        } else if coinc(d1, d2, tol) {
            par[0] = lin2d_parameter(lin2d, &GpPnt2d::new(xmax, ymin));
            par[1] = lin2d_parameter(lin2d, &GpPnt2d::new(xmin, ymin));
            nbi = 2;
        }
    }
    if nbi == 2 {
        if (par[0] - par[1]).abs() > tol {
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
pub fn tolerances(kind_a: SurfaceKind, kind_b: SurfaceKind, tol_tang: &mut f64) {
    if (kind_a == SurfaceKind::Cylinder && kind_b == SurfaceKind::Torus)
        || (kind_b == SurfaceKind::Cylinder && kind_a == SurfaceKind::Torus)
    {
        *tol_tang *= 0.1;
    }
}

/// Golden-section max of curve-to-surface distance on `[first, last]`.
fn find_max_distance_interval(
    c: &dyn Curve,
    first: f64,
    last: f64,
    surf: &dyn Surface,
    eps: f64,
) -> f64 {
    let cf = 0.6180339887498948;
    let mut a = first;
    let mut b = last;
    let mut x1 = b - cf * (b - a);
    let mut f1 = max_distance(c, x1, surf);
    let mut x2 = a + cf * (b - a);
    let mut f2 = max_distance(c, x2, surf);
    while (x1 - x2).abs() > eps {
        if f1 > f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - cf * (b - a);
            f1 = max_distance(c, x1, surf);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + cf * (b - a);
            f2 = max_distance(c, x2, surf);
        }
    }
    let mut f = max_distance(c, 0.5 * (x1 + x2), surf);
    if f1 > f {
        f = f1;
    }
    if f2 > f {
        f = f2;
    }
    f
}

fn max_distance(c: &dyn Curve, t: f64, surf: &dyn Surface) -> f64 {
    let a_p = c.d0(t);
    project_point_on_surface(surf, &a_p, 0.0)
        .map(|p| p.distance)
        .unwrap_or(0.0)
}

/// `FindMaxDistance` over a face (`IntTools_FaceFace.cxx:2899`).
pub fn find_max_distance_face(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    face: &Face,
) -> f64 {
    let Some(surf) = BRepTool::face_surface(face) else {
        return 0.0;
    };
    let nb_s = 11.0;
    let dt = (last - first) / nb_s;
    let mut d_max = 0.0;
    let an_eps = 1.0e-4 * dt;
    let mut t2 = first;
    loop {
        let t1 = t2;
        t2 += dt;
        if t2 > last {
            break;
        }
        let d = find_max_distance_interval(curve, t1, t2, surf.as_ref(), an_eps);
        if d > d_max {
            d_max = d;
        }
    }
    d_max
}

/// `CheckPCurve` (`IntTools_FaceFace.cxx:3010`) with a single CN interval.
pub fn face_pcurve_check(pc: &dyn Curve2d, face: &Face) -> bool {
    const N_POINTS: i32 = 23;
    let (mut umin, mut umax, mut vmin, mut vmax) = crate::brep_surface::face_uv_bounds(face);
    let tol_u = ((umax - umin) * 0.01).max(CONFUSION);
    let tol_v = ((vmax - vmin) * 0.01).max(CONFUSION);
    let fp = pc.first_parameter();
    let lp = pc.last_parameter();
    if let Some(surf) = BRepTool::face_surface(face) {
        let pnt = pc.d0(0.5 * (fp + lp));
        let (u, v) = (pnt.x(), pnt.y());
        if surf.is_u_periodic() {
            let per = crate::int_tools_wline::u_period(surf.as_ref()).unwrap_or(0.0);
            if per > 0.0 {
                let mut nshift = ((u - umin) / per) as i32;
                if u < umin + per * nshift as f64 {
                    nshift -= 1;
                }
                umin += per * nshift as f64;
                umax += per * nshift as f64;
            }
        }
        if surf.is_v_periodic() {
            let per = crate::int_tools_wline::v_period(surf.as_ref()).unwrap_or(0.0);
            if per > 0.0 {
                let mut nshift = ((v - vmin) / per) as i32;
                if v < vmin + per * nshift as f64 {
                    nshift -= 1;
                }
                vmin += per * nshift as f64;
                vmax += per * nshift as f64;
            }
        }
    }
    let d_t = (lp - fp) / N_POINTS as f64;
    let mut t = fp;
    for _i in 1..N_POINTS {
        t += d_t;
        let p2d = pc.d0(t);
        let (u, v) = (p2d.x(), p2d.y());
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
            if let (Some(cyl), Some(sph)) = (gp_cylinder_from_surface(sa), gp_sphere_from_surface(sb)) {
                an_with_pc = approx_with_pcurves(&cyl, &sph);
            }
        } else if typs1 == SurfaceKind::Sphere && typs2 == SurfaceKind::Cylinder {
            if let (Some(sph), Some(cyl)) = (gp_sphere_from_surface(sa), gp_cylinder_from_surface(sb)) {
                an_with_pc = approx_with_pcurves(&cyl, &sph);
            }
        }
        if !an_with_pc {
            my_tol_approx = 1.0e-5;
            an_approx1 = false;
            an_approx2 = false;
        }

        let mut seq_of_l: Vec<WLine> = Vec::new();
        let b_is_decomposited = decomposition_of_wline(
            wl,
            sa,
            sb,
            fa,
            fb,
            &line_parts,
            b_avoid,
            self.tol,
            &mut seq_of_l,
        );

        let nbiter = if b_is_decomposited {
            seq_of_l.len() as i32
        } else {
            line_parts.len() as i32
        };

        let (i_deg_min, i_deg_max, i_nb_iter) = approx_parameters(typs1, typs2);
        let mut _tol_tang = self.tol;
        tolerances(typs1, typs2, &mut _tol_tang);

        let mut out = Vec::new();
        for i in 1..=nbiter {
            let (piece, ifprm, ilprm): (&WLine, i32, i32) = if b_is_decomposited {
                let w = &seq_of_l[(i as usize).saturating_sub(1)];
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
        let tol_f_max = BRepTool::face_tolerance(fa).max(BRepTool::face_tolerance(fb));
        let Some(s1) = BRepTool::face_surface(fa) else {
            return;
        };
        let Some(s2) = BRepTool::face_surface(fb) else {
            return;
        };
        for c in &mut self.result.curves {
            let first = c.curve.first_parameter();
            let last = c.curve.last_parameter();
            let mut tol_c = c.tolerance;
            for j in 0..2 {
                let c2d = if j == 0 { &c.pcurve1 } else { &c.pcurve2 };
                if let Some(pc) = c2d {
                    let s = if j == 0 { &s1 } else { &s2 };
                    if let Some((d, _a_t)) =
                        compute_tolerance(&c.curve, pc, s, first, last)
                    {
                        if d > tol_c {
                            tol_c = d;
                        }
                    }
                } else {
                    let f = if j == 0 { fa } else { fb };
                    let d = find_max_distance_face(c.curve.as_ref(), first, last, f);
                    if d > tol_c {
                        tol_c = d;
                    }
                }
            }
            c.tolerance = tol_c;
            if c.tangential_tolerance < tol_f_max {
                c.tangential_tolerance = tol_f_max;
            }
        }
    }

    /// `PrepareLines3D` (`IntTools_FaceFace.cxx:1932`).
    pub fn prepare_lines_3d(&mut self, b_to_split: bool) {
        let mut new_cvs: Vec<FaceFaceCurve> = Vec::new();
        for c in &self.result.curves {
            if b_to_split {
                let split = split_curve(c);
                if split.is_empty() {
                    new_cvs.push(c.clone());
                } else {
                    new_cvs.extend(split);
                }
            } else {
                new_cvs.push(c.clone());
            }
        }
        let type1 = self.face1.as_ref().and_then(BRepTool::face_surface);
        let type2 = self.face2.as_ref().and_then(BRepTool::face_surface);
        let k1 = type1.as_ref().map(|s| classify_surface(s.as_ref()));
        let k2 = type2.as_ref().map(|s| classify_surface(s.as_ref()));
        let plane_cone = matches!(
            (k1, k2),
            (Some(SurfaceKind::Plane), Some(SurfaceKind::Cone))
                | (Some(SurfaceKind::Cone), Some(SurfaceKind::Plane))
        );
        if plane_cone && new_cvs.len() == 4 && new_cvs[0].kind == CurveKind::Line {
            new_cvs = reject_lines(&new_cvs);
        }
        self.result.curves = new_cvs;
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
    tol: f64,
) -> (f64, f64) {
    let mut f = fprm;
    let mut l = lprm;
    if let Some(anew) = parameter_out_of_boundary(fprm, newc, fa, fb, lprm, false, tol) {
        f = anew;
    }
    if let Some(anew) = parameter_out_of_boundary(lprm, newc, fa, fb, fprm, true, tol) {
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
