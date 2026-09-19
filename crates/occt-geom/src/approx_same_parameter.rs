//! `Approx_SameParameter`.
//! Source: `Approx_SameParameter.cxx` (`Build` 318-539, `BuildInitialDistribution`
//! 547-579, `IncreaseInitialNbSamples` 588-648, `CheckSameParameter` 653-766,
//! `ComputeTangents` 770-809, `Interpolate` 813-853, `Check` 192-266,
//! Evaluator 64-102, `ComputeTolReached` 153-188, `IncreaseNbPoles` 857-991).
//!
//! `CheckSameParameter` uses `Extrema_LocateExtPC` then `Extrema_ExtPC`
//! whole-space search (`cxx:700-758`).

use std::sync::Arc;

use occt_core::bspl::banded_interp::{interpolate_contact, knot_sequence};
use occt_core::bspl::eval::{eval_curve, eval_curve_d1};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};

use crate::adv_approx::ApproxAFunction1dPair;
use crate::curve::Curve;
use crate::surface::Surface;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

const NB_SAMPLES: usize = 22;
const MAX_ARRAY: usize = 1000;

/// Result of `Approx_SameParameter`.
pub struct ApproxSameParameter {
    pub done: bool,
    pub same_parameter: bool,
    pub tol_reached: f64,
    /// Rebuilt 2d curve when `IsDone && !IsSameParameter` (`cxx:1648-1654`).
    pub curve2d: Option<Arc<dyn Curve2d>>,
}

struct DistData {
    pc3d: [f64; MAX_ARRAY],
    pc2d: [f64; MAX_ARRAY],
    nb: usize,
    c3d_f: f64,
    c3d_l: f64,
    c2d_f: f64,
    c2d_l: f64,
    tol: f64,
}

fn cons_value(c2d: &dyn Curve2d, s: &dyn Surface, t: f64) -> GpPnt {
    let uv = c2d.d0(t);
    s.d0(uv.x(), uv.y())
}

fn cons_d1(c2d: &dyn Curve2d, s: &dyn Surface, t: f64) -> (GpPnt, GpVec) {
    let (uv, duv) = c2d.d1(t);
    let (p, du, dv) = s.d1(uv.x(), uv.y());
    let d = du.multiplied_scalar(duv.x()).added(&dv.multiplied_scalar(duv.y()));
    (p, d)
}

fn compute_tol_reached(c3d: &dyn Curve, c2d: &dyn Curve2d, s: &dyn Surface, first: f64, last: f64) -> f64 {
    let nbp = 2 * NB_SAMPLES;
    let mut d2: f64 = 0.0;
    for i in 0..=nbp {
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.d0(u);
        let pcons = cons_value(c2d, s, u);
        if !pcons.x().is_finite() || !pcons.y().is_finite() || !pcons.z().is_finite() {
            return INFINITE;
        }
        d2 = d2.max(pc3d.square_distance(&pcons));
    }
    (1.05 * d2.sqrt()).max(CONFUSION)
}

fn project_point_on_curve(
    init: f64,
    point: &GpPnt,
    tolerance: f64,
    curve: &dyn Curve,
    first: f64,
    last: f64,
) -> Option<f64> {
    let mut param = init;
    for _ in 0..30 {
        let (a_point, d1, d2) = curve.d2(param);
        let vector = GpVec::from_pnts(&a_point, point);
        let func = vector.dot(&d1);
        if func.abs() < tolerance * d1.magnitude() {
            return Some(param);
        }
        let func_derivative = vector.dot(&d2) - d1.dot(&d1);
        if func_derivative.abs() > 1.0e-12 {
            param -= func / func_derivative;
        }
        param = param.max(first).min(last);
    }
    None
}

fn increase_initial_nb_samples(data: &mut DistData, c2d: &dyn Curve2d) -> bool {
    let intervals = c2d.parameter_intervals(2);
    let length = intervals.len();
    if length == 0 {
        return true;
    }
    let at = |i1: usize| intervals[i1 - 1];
    let mut inter = 1usize;
    let mut nb_int = length;
    while inter <= nb_int && at(inter) <= data.c3d_f + PCONFUSION {
        inter += 1;
    }
    while nb_int > 0 && at(nb_int) >= data.c3d_l - PCONFUSION {
        nb_int -= 1;
    }
    let mut new_par: Vec<f64> = Vec::new();
    new_par.push(data.c3d_f);
    let mut ii = 1usize;
    while inter <= nb_int || (ii < NB_SAMPLES && inter <= length) {
        if at(inter) < data.pc2d[ii] {
            new_par.push(at(inter));
            if data.pc2d[ii] - at(inter) <= PCONFUSION {
                ii += 1;
                if ii > NB_SAMPLES {
                    ii = NB_SAMPLES;
                }
            }
            inter += 1;
        } else {
            if (at(inter) - data.pc2d[ii]) > PCONFUSION {
                new_par.push(data.pc2d[ii]);
            }
            ii += 1;
        }
    }
    data.nb = new_par.len();
    if data.nb > MAX_ARRAY - 1 {
        return false;
    }
    for i in 1..data.nb {
        let v = new_par[i];
        data.pc2d[i] = v;
        data.pc3d[i] = v;
    }
    data.pc3d[data.nb] = data.c3d_l;
    data.pc2d[data.nb] = data.c2d_l;
    true
}

fn build_initial(data: &mut DistData, c2d: &dyn Curve2d) -> bool {
    let deltacons = (data.c2d_l - data.c2d_f) / NB_SAMPLES as f64;
    let deltac3d = (data.c3d_l - data.c3d_f) / NB_SAMPLES as f64;
    let mut wcons = data.c2d_f;
    let mut wc3d = data.c3d_f;
    for ii in 0..NB_SAMPLES {
        data.pc2d[ii] = wcons;
        data.pc3d[ii] = wc3d;
        wcons += deltacons;
        wc3d += deltac3d;
    }
    data.nb = NB_SAMPLES;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    if c2d.continuity() < 2 {
        return increase_initial_nb_samples(data, c2d);
    }
    true
}

fn check_same_parameter(
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    sq_dist: &mut f64,
) -> bool {
    let tol2 = data.tol * data.tol;
    let mut is_same = true;
    let pcons = cons_value(c2d, s, data.c2d_f);
    let pc3d = c3d.d0(data.c3d_f);
    let mut dmax2 = pcons.square_distance(&pc3d);
    let pcons = cons_value(c2d, s, data.c2d_l);
    let pc3d = c3d.d0(data.c3d_l);
    dmax2 = dmax2.max(pcons.square_distance(&pc3d));

    let mut count = 1usize;
    let mut previousp = data.c3d_f;
    let mut initp = 0.0;
    let bornesup = data.c3d_l - PCONFUSION;
    let mut is_proj_ok = false;
    for ii in 1..data.nb {
        let pcons = cons_value(c2d, s, data.pc2d[ii]);
        let pc3d_pt = c3d.d0(data.pc3d[ii]);
        let dist2 = pcons.square_distance(&pc3d_pt);
        let is_use = dist2 <= tol2 && data.pc3d[ii] > data.pc3d[count - 1] + PCONFUSION;
        if is_use {
            if dmax2 < dist2 {
                dmax2 = dist2;
            }
            initp = data.pc3d[ii];
            previousp = initp;
            data.pc3d[count] = data.pc3d[ii];
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
            continue;
        }
        if !is_proj_ok {
            initp = data.pc3d[ii];
        }
        is_proj_ok = false;
        is_same = false;
        let mut curp = initp;
        if let Some(p) = project_point_on_curve(initp, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
        {
            curp = p;
            is_proj_ok = true;
        }
        is_proj_ok = is_proj_ok && curp > previousp + PCONFUSION && curp < bornesup;
        if is_proj_ok {
            initp = curp;
            previousp = curp;
            data.pc3d[count] = curp;
            data.pc2d[count] = data.pc2d[ii];
            count += 1;
            continue;
        }
        // Whole parameter space search (`Extrema_ExtPC`, cxx:728-758).
        let extrema = crate::extrema_pc::point_curve_extrema_all(c3d, &pcons);
        if extrema.is_empty() {
            continue;
        }
        let mut best_u: Option<f64> = None;
        let mut best_d2 = f64::MAX;
        for e in &extrema {
            if e.u1 < data.c3d_f - PCONFUSION || e.u1 > data.c3d_l + PCONFUSION {
                continue;
            }
            let d2 = e.p2.square_distance(&pcons);
            if d2 < best_d2 {
                best_d2 = d2;
                best_u = Some(e.u1);
            }
        }
        if let Some(curp) = best_u {
            if curp > previousp + PCONFUSION && curp < bornesup {
                initp = curp;
                previousp = curp;
                data.pc3d[count] = curp;
                data.pc2d[count] = data.pc2d[ii];
                count += 1;
                is_proj_ok = true;
            }
        }
    }
    data.nb = count;
    data.pc2d[data.nb] = data.c2d_l;
    data.pc3d[data.nb] = data.c3d_l;
    *sq_dist = dmax2;
    is_same
}

fn compute_tangents(c3d: &dyn Curve, c2d: &dyn Curve2d, s: &dyn Surface) -> Option<(f64, f64)> {
    const SMALL: f64 = 1.0e-12;
    let a_param_first = c3d.first_parameter();
    let (_, a_vec_cons) = cons_d1(c2d, s, a_param_first);
    let (_, a_vec) = c3d.d1(a_param_first);
    let mag = a_vec_cons.magnitude();
    if mag <= SMALL {
        return None;
    }
    let first_t = a_vec.magnitude() / mag;
    let a_param_last = c3d.last_parameter();
    let (_, a_vec_cons) = cons_d1(c2d, s, a_param_last);
    let (_, a_vec) = c3d.d1(a_param_last);
    let mag = a_vec_cons.magnitude();
    if mag <= SMALL {
        return None;
    }
    let last_t = a_vec.magnitude() / mag;
    Some((first_t, last_t))
}

fn interpolate_reparam(data: &DistData, tang_first: f64, tang_last: f64) -> Option<(Vec<f64>, Vec<f64>)> {
    let num_poles = data.nb + 3;
    let num_knots = data.nb + 7;
    let mut poles = vec![0.0; num_poles];
    let mut flat = vec![0.0; num_knots];
    let mut contact = vec![0i32; num_poles];
    let mut parameters = vec![0.0; num_poles];
    contact[1] = 1;
    contact[num_poles - 2] = 1;
    for i in 0..4 {
        flat[i] = data.c3d_f;
        flat[num_poles + i] = data.c3d_l;
    }
    poles[0] = data.c2d_f;
    poles[num_poles - 1] = data.c2d_l;
    poles[1] = tang_first;
    poles[num_poles - 2] = tang_last;
    parameters[0] = data.c3d_f;
    parameters[1] = data.c3d_f;
    parameters[num_poles - 2] = data.c3d_l;
    parameters[num_poles - 1] = data.c3d_l;
    for ii in 3..=(num_poles - 2) {
        poles[ii - 1] = data.pc2d[ii - 2];
        parameters[ii - 1] = data.pc3d[ii - 2];
        flat[ii + 1] = data.pc3d[ii - 2];
    }
    interpolate_contact(3, &flat, &parameters, &contact, &mut poles, 1).ok()?;
    Some((poles, flat))
}

fn eval_1d_d0(flat: &[f64], poles: &[f64], u: f64) -> f64 {
    let pts: Vec<GpPnt> = poles.iter().map(|&x| GpPnt::new(x, 0.0, 0.0)).collect();
    eval_curve(&pts, flat, 3, u).x()
}

fn eval_1d_d1(flat: &[f64], poles: &[f64], u: f64) -> (f64, f64) {
    let pts: Vec<GpPnt> = poles.iter().map(|&x| GpPnt::new(x, 0.0, 0.0)).collect();
    let (p, d) = eval_curve_d1(&pts, flat, 3, u);
    (p.x(), d.x())
}

fn check_reparam(
    flat: &[f64],
    poles: &[f64],
    nbp: usize,
    pc3d: &[f64],
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    tol: &mut f64,
    oldtol: f64,
) -> bool {
    let mut a_param_first = 3.0 * pc3d[0] - 2.0 * pc3d[nbp - 1];
    let mut a_param_last = 3.0 * pc3d[nbp - 1] - 2.0 * pc3d[0];
    let first_par = c2d.first_parameter();
    let last_par = c2d.last_parameter();
    if a_param_first < first_par {
        a_param_first = first_par;
    }
    if a_param_last > last_par {
        a_param_last = last_par;
    }
    let d = *tol;
    let nn = 2 * nbp;
    let unsurnn = 1.0 / nn as f64;
    let mut tprev = a_param_first;
    let mut d2: f64 = 0.0;
    for i in 0..=nn {
        let t = unsurnn * i as f64;
        let tc3d = pc3d[0] * (1.0 - t) + pc3d[nbp - 1] * t;
        let pc3d_pt = c3d.d0(tc3d);
        let tcons = eval_1d_d0(flat, poles, tc3d);
        if tcons < tprev || tcons > a_param_last {
            *tol = INFINITE;
            return false;
        }
        tprev = tcons;
        let pcons = cons_value(c2d, s, tcons);
        d2 = d2.max(pc3d_pt.square_distance(&pcons));
    }
    *tol = d2.sqrt();
    for i in 1..poles.len() {
        if poles[i - 1] > poles[i] {
            return false;
        }
    }
    *tol <= d || *tol > 0.8 * oldtol
}

/// `Precision::Parametric(R3d)` (`Precision.hxx:328`, `= R3d * 0.01`), the
/// `default` arm of `GeomAdaptor_Surface::UResolution` / `VResolution`.
#[inline]
fn parametric_resolution(r3d: f64) -> f64 {
    r3d * PCONFUSION / CONFUSION
}

/// `GeomAdaptor_Surface::UResolution` (`GeomAdaptor_Surface.cxx:1818-1896`).
/// `r3d` is OCCT's `R3d` argument.
///
/// Shared by every OCCT call site that asks a surface for a parametric
/// resolution: `Approx_SameParameter.cxx:429` / `:504`
/// (`mySurf->UResolution(besttol)`), `BRepLib.cxx:1085`
/// (`DSdu = 1./surf->UResolution(1.)`) and `IntTools_WLineTool.cxx:153`
/// (`aGAS1.UResolution(aDelta)`).
pub fn u_resolution(s: &dyn Surface, r3d: f64) -> f64 {
    // `GeomAdaptor_Surface::load` peels a `Geom_RectangularTrimmedSurface`
    // down to its basis surface (`cxx:423-430`), so the switch below sees
    // the basis type. NOTE: OCCT keeps the trimmed range as `myUFirst..myVLast`;
    // this recursion uses the basis range, which only matters for the
    // unbounded-V Cone guard (`cxx:1855-1858`).
    if let Some(b) = s.rectangular_trimmed_basis() {
        return u_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_SurfaceOfExtrusion: BasisCurve->Resolution(R3d)`
    // (`cxx:1824-1827`), live through `GeomSurfaceOfLinearExtrusion` and the
    // trimmed/offset wrappers that forward `extrusion_basis_curve`.
    if let Some(c) = s.extrusion_basis_curve() {
        return c.resolution(r3d);
    }
    // `case GeomAbs_OffsetSurface: BasisAdaptor->UResolution(R3d)`
    // (`cxx:1882-1885`).
    if let Some(b) = s.offset_basis_surface() {
        return u_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_Cone` (`cxx:1853-1866`). The guard is the infinite-domain
    // one: `myVLast - myVFirst > 1.e10` means "not truly bounded".
    if let Some((radius, semi_angle)) = s.cone_ref() {
        let (v_first, v_last) = s.v_range();
        if v_last - v_first > 1.0e10 {
            return parametric_resolution(r3d);
        }
        // `S->VIso(V)` radius is `Radius + V*sin(Angle)`, taken absolute
        // (`ElSLib.cxx:1793-1811`).
        let rayon1 = (radius + v_last * semi_angle.sin()).abs();
        let rayon2 = (radius + v_first * semi_angle.sin()).abs();
        let r = if rayon1 > rayon2 { rayon1 } else { rayon2 };
        return if r > CONFUSION { r3d / r } else { 0.0 };
    }
    // `case GeomAbs_Plane: return R3d` (`cxx:1867-1869`).
    if s.gp_pln().is_some() {
        return r3d;
    }
    // `GeomAbs_Torus` / `GeomAbs_Sphere` / `GeomAbs_Cylinder` compute
    // `Res = R3d / (2. * R)` and then fall into the `2.*asin(Res)` tail
    // (`cxx:1828-1852`, `:1890-1895`).
    let res = if let Some(t) = s.gp_torus() {
        let r = t.major_radius + t.minor_radius;
        if r > CONFUSION { r3d / (2.0 * r) } else { 0.0 }
    } else if let Some(sp) = s.gp_sphere() {
        if sp.radius > CONFUSION { r3d / (2.0 * sp.radius) } else { 0.0 }
    } else if let Some(cy) = s.gp_cylinder() {
        if cy.radius > CONFUSION { r3d / (2.0 * cy.radius) } else { 0.0 }
    } else if let Some((ures, _)) = s.uv_resolution(r3d) {
        // `GeomAbs_BSplineSurface` (`cxx:1875-1881`) and the offset/revolution
        // wrappers. UNPORTED: `Geom_BezierSurface::Resolution`
        // (`GeomAdaptor_Surface.cxx:1870-1874`) has no Rust implementation, so
        // a Bezier surface falls through to the `default` arm below.
        return ures;
    } else {
        // `default: return Precision::Parametric(R3d)` (`cxx:1886-1887`).
        return parametric_resolution(r3d);
    };
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * std::f64::consts::PI
    }
}

/// `GeomAdaptor_Surface::VResolution` (`GeomAdaptor_Surface.cxx:1900-1958`).
/// `r3d` is OCCT's `R3d` argument.
pub fn v_resolution(s: &dyn Surface, r3d: f64) -> f64 {
    // Trimmed-surface unwrap, as in `u_resolution` above.
    if let Some(b) = s.rectangular_trimmed_basis() {
        return v_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_SurfaceOfRevolution: BasisCurve->Resolution(R3d)`
    // (`cxx:1905-1908`).
    if let Some(c) = s.revolution_basis_curve() {
        return c.resolution(r3d);
    }
    // `case GeomAbs_OffsetSurface: BasisAdaptor->VResolution(R3d)`
    // (`cxx:1944-1947`).
    if let Some(b) = s.offset_basis_surface() {
        return v_resolution(b.as_ref(), r3d);
    }
    // `GeomAbs_Torus` (MinorRadius, `cxx:1910-1918`) and `GeomAbs_Sphere`
    // (`cxx:1919-1927`) share the `2.*asin(Res)` tail (`cxx:1952-1957`).
    let res = if let Some(t) = s.gp_torus() {
        let r = t.minor_radius;
        if r > CONFUSION { r3d / (2.0 * r) } else { 0.0 }
    } else if let Some(sp) = s.gp_sphere() {
        if sp.radius > CONFUSION { r3d / (2.0 * sp.radius) } else { 0.0 }
    } else if s.gp_pln().is_some()
        || s.gp_cylinder().is_some()
        || s.cone_ref().is_some()
        || s.extrusion_basis_curve().is_some()
    {
        // `case GeomAbs_SurfaceOfExtrusion / GeomAbs_Cylinder / GeomAbs_Cone /
        // GeomAbs_Plane: return R3d` (`cxx:1928-1933`).
        return r3d;
    } else if let Some((_, vres)) = s.uv_resolution(r3d) {
        return vres;
    } else {
        // `default: return Precision::Parametric(R3d)` (`cxx:1948-1949`).
        return parametric_resolution(r3d);
    };
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * std::f64::consts::PI
    }
}

fn curve2d_from_approx(approx: &ApproxAFunction1dPair) -> Option<Arc<dyn Curve2d>> {
    if !approx.has_result || approx.poles_u.is_empty() {
        return None;
    }
    let deg = approx.degree.max(1) as usize;
    let flat = knot_sequence(&approx.knots, &approx.mults, approx.degree);
    let c = Geom2dBSplineCurve::new(approx.poles_u.clone(), approx.poles_v.clone(), flat, deg).ok()?;
    Some(Arc::new(c))
}

/// `Approx_SameParameter::IncreaseNbPoles` (`cxx:857-991`).
fn increase_nb_poles(
    poles: &[f64],
    flat: &[f64],
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    best_sq_tol: &mut f64,
) -> bool {
    let mut new_pc2d = [0.0; MAX_ARRAY];
    let mut new_pc3d = [0.0; MAX_ARRAY];
    let mut newcount = 0usize;
    for ii in 0..data.nb {
        new_pc2d[newcount] = data.pc2d[ii];
        new_pc3d[newcount] = data.pc3d[ii];
        newcount += 1;
        if data.nb - ii + newcount == MAX_ARRAY {
            continue;
        }
        let mid_c3d = 0.5 * (data.pc3d[ii] + data.pc3d[ii + 1]);
        let eval_result = eval_1d_d0(flat, poles, mid_c3d);
        if eval_result < data.pc2d[ii] || eval_result > data.pc2d[ii + 1] {
            let ucons = 0.5 * (data.pc2d[ii] + data.pc2d[ii + 1]);
            let uc3d = 0.5 * (data.pc3d[ii] + data.pc3d[ii + 1]);
            let pcons = cons_value(c2d, s, ucons);
            if let Some(curp) =
                project_point_on_curve(uc3d, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
            {
                let dist_2 = c3d.d0(curp).square_distance(&pcons);
                if dist_2 > *best_sq_tol {
                    *best_sq_tol = dist_2;
                }
                if curp > data.pc3d[ii] + PCONFUSION && curp < data.pc3d[ii + 1] - PCONFUSION {
                    new_pc3d[newcount] = curp;
                    new_pc2d[newcount] = ucons;
                    newcount += 1;
                }
            }
        }
    }
    new_pc3d[newcount] = data.pc3d[data.nb];
    new_pc2d[newcount] = data.pc2d[data.nb];
    if data.nb != newcount && newcount < MAX_ARRAY - 1 {
        data.pc2d = new_pc2d;
        data.pc3d = new_pc3d;
        data.nb = newcount;
        return true;
    }

    newcount = 0;
    for n in 0..data.nb {
        new_pc3d[newcount] = data.pc3d[n];
        new_pc2d[newcount] = data.pc2d[n];
        newcount += 1;
        if data.nb - n + newcount == MAX_ARRAY {
            continue;
        }
        let ucons = 0.5 * (data.pc2d[n] + data.pc2d[n + 1]);
        let uc3d = 0.5 * (data.pc3d[n] + data.pc3d[n + 1]);
        let pcons = cons_value(c2d, s, ucons);
        if let Some(curp) =
            project_point_on_curve(uc3d, &pcons, data.tol, c3d, data.c3d_f, data.c3d_l)
        {
            let dist_2 = c3d.d0(curp).square_distance(&pcons);
            if dist_2 > *best_sq_tol {
                *best_sq_tol = dist_2;
            }
            if curp > data.pc3d[n] + PCONFUSION && curp < data.pc3d[n + 1] - PCONFUSION {
                new_pc3d[newcount] = curp;
                new_pc2d[newcount] = ucons;
                newcount += 1;
            }
        }
    }
    new_pc3d[newcount] = data.pc3d[data.nb];
    new_pc2d[newcount] = data.pc2d[data.nb];
    if data.nb != newcount {
        data.pc2d = new_pc2d;
        data.pc3d = new_pc3d;
        data.nb = newcount;
        return true;
    }
    false
}

fn approx_evaluator_eval<'a>(
    flat: &'a [f64],
    poles: &'a [f64],
    c2d: &'a dyn Curve2d,
) -> impl Fn(f64, i32, &mut [f64]) -> i32 + 'a {
    move |param: f64, der: i32, out: &mut [f64]| -> i32 {
        if der == 0 {
            let tcons = eval_1d_d0(flat, poles, param);
            let uv = c2d.d0(tcons);
            if out.len() >= 2 {
                out[0] = uv.x();
                out[1] = uv.y();
            }
            0
        } else if der == 1 {
            let (tcons, dt) = eval_1d_d1(flat, poles, param);
            let (_, duv) = c2d.d1(tcons);
            if out.len() >= 2 {
                out[0] = duv.x() * dt;
                out[1] = duv.y() * dt;
            }
            0
        } else {
            1
        }
    }
}

fn rebuild_pcurve(
    data: &mut DistData,
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    tang_first: f64,
    tang_last: f64,
) -> (bool, f64, Option<Arc<dyn Curve2d>>) {
    let cont = if c2d.continuity() > 2 {
        Shape::C1
    } else if c2d.continuity() == 0 {
        Shape::C0
    } else {
        Shape::C1
    };
    let mut besttol2 = data.tol * data.tol;
    let mut tolsov = INFINITE;
    let mut done = false;
    let mut curve2d: Option<Arc<dyn Curve2d>> = None;
    let mut tol_reached = compute_tol_reached(c3d, c2d, s, data.c3d_f, data.c3d_l);
    let mut interpolok = false;
    let mut has_count_changed = false;
    // `cxx:400-480`: densify loop over IncreaseNbPoles until AdvApprox accepts.
    loop {
        let Some((poles, flat)) = interpolate_reparam(data, tang_first, tang_last) else {
            return (false, tol_reached, None);
        };
        let mut algtol = besttol2.sqrt();
        interpolok = check_reparam(
            &flat,
            &poles,
            data.nb + 1,
            &data.pc3d,
            c3d,
            c2d,
            s,
            &mut algtol,
            tolsov,
        );
        tolsov = algtol;
        if interpolok {
            let besttol = besttol2.sqrt();
            let tol_u = u_resolution(s, besttol);
            let tol_v = v_resolution(s, besttol);
            let eval = approx_evaluator_eval(&flat, &poles, c2d);
            if let Ok(approx) = ApproxAFunction1dPair::approx(
                data.c3d_f,
                data.c3d_l,
                cont,
                11,
                1000,
                tol_u,
                tol_v,
                &eval,
            ) {
                if approx.done || approx.has_result {
                    if let Some(new_c2d) = curve2d_from_approx(&approx) {
                        let new_tol =
                            compute_tol_reached(c3d, new_c2d.as_ref(), s, data.c3d_f, data.c3d_l);
                        const MULT: f64 = 250.0;
                        if new_tol < MULT * besttol {
                            done = true;
                            tol_reached = new_tol;
                            curve2d = Some(new_c2d);
                            break;
                        } else if data.nb < MAX_ARRAY - 1 {
                            interpolok = false;
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        if !interpolok {
            has_count_changed =
                increase_nb_poles(&poles, &flat, data, c3d, c2d, s, &mut besttol2);
        }
        if interpolok || !has_count_changed {
            break;
        }
    }

    // `cxx:482-539`: post-loop fallback — accept AdvApprox if better than original.
    if !done {
        tol_reached = compute_tol_reached(c3d, c2d, s, data.c3d_f, data.c3d_l);
        let Some((poles, flat)) = interpolate_reparam(data, tang_first, tang_last) else {
            return (false, tol_reached, None);
        };
        let besttol = besttol2.sqrt();
        let tol_u = u_resolution(s, besttol);
        let tol_v = v_resolution(s, besttol);
        let eval = approx_evaluator_eval(&flat, &poles, c2d);
        let Ok(approx) = ApproxAFunction1dPair::approx(
            data.c3d_f,
            data.c3d_l,
            cont,
            11,
            40,
            tol_u,
            tol_v,
            &eval,
        ) else {
            return (false, tol_reached, None);
        };
        if !approx.done && !approx.has_result {
            return (false, tol_reached, None);
        }
        let Some(new_c2d) = curve2d_from_approx(&approx) else {
            return (false, tol_reached, None);
        };
        let approx_tol = compute_tol_reached(c3d, new_c2d.as_ref(), s, data.c3d_f, data.c3d_l);
        if approx_tol < tol_reached {
            tol_reached = approx_tol;
            curve2d = Some(new_c2d);
        }
        done = true;
    }
    (done, tol_reached, curve2d)
}

impl ApproxSameParameter {
    /// Adaptor constructor used by `BRepLib::SameParameter` (`cxx:1631`).
    pub fn new(
        c3d: &dyn Curve,
        c2d: &dyn Curve2d,
        surf: &dyn Surface,
        first: f64,
        last: f64,
        tol: f64,
    ) -> Self {
        let mut data = DistData {
            pc3d: [0.0; MAX_ARRAY],
            pc2d: [0.0; MAX_ARRAY],
            nb: 0,
            c3d_f: first,
            c3d_l: last,
            c2d_f: first,
            c2d_l: last,
            tol,
        };
        if !build_initial(&mut data, c2d) {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        let keep_min = data.nb - ((0.3 * data.nb as f64) as usize);
        let mut sq = 0.0;
        let same = check_same_parameter(&mut data, c3d, c2d, surf, &mut sq);
        if same {
            return Self {
                done: true,
                same_parameter: true,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        if data.nb < keep_min {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        }
        let Some((tang_first, tang_last)) = compute_tangents(c3d, c2d, surf) else {
            return Self {
                done: false,
                same_parameter: false,
                tol_reached: compute_tol_reached(c3d, c2d, surf, first, last),
                curve2d: None,
            };
        };
        let (done, tol_reached, curve2d) =
            rebuild_pcurve(&mut data, c3d, c2d, surf, tang_first, tang_last);
        Self {
            done,
            same_parameter: false,
            tol_reached,
            curve2d,
        }
    }
}
