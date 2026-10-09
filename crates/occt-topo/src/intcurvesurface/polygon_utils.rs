//! `IntCurveSurface_TheHCurveTool::SamplePars` / `NbSamples` and
//! `Adaptor3d_HSurfaceTool::NbSamplesU/V`.
//!
//! Source: `IntCurveSurface_TheHCurveTool.cxx:31-70, 72-299` and
//! `Adaptor3d_HSurfaceTool.cxx:27-116` (TKGeomAlgo / TKG3d). These decide the
//! sample parameters from which `IntCurveSurface_HInter` builds its curve
//! polygon and surface polyhedron, so the interference engine and the
//! polyhedron parameter recovery both inherit the OCCT sampling counts.

use occt_geom::{Curve, Surface};

/// `myMinPnts` (`IntCurveSurface_TheHCurveTool.cxx:32`).
pub const MY_MIN_PNTS: usize = 5;

/// `IntCurveSurface_TheHCurveTool::NbSamples` (`...cxx:31-70`).
pub fn nb_samples(curve: &dyn Curve, u0: f64, u1: f64) -> usize {
    const NBS_OTHER: f64 = 10.0;
    let mut nbs = NBS_OTHER;

    if curve.is_line() {
        nbs = 2.0;
    } else if curve.bezier_poles().is_some() {
        nbs = 3.0 + curve.bezier_poles().map_or(0, |p| p.len()) as f64;
    } else if let Some(degree) = curve.nurbs_degree() {
        let nb_knots = curve.bspline_knots().map_or(0, unique_knot_count) as f64;
        nbs = nb_knots * degree as f64 * (curve.last_parameter() - curve.first_parameter());
        let a_range = u1 - u0;
        if a_range.abs() > occt_core::precision::PCONFUSION {
            nbs /= a_range;
        }
        if nbs < 2.0 {
            nbs = 2.0;
        }
    }

    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as usize
}

/// Number of distinct knots in a flat knot vector (OCCT `NbKnots`).
pub fn unique_knot_count(knots: &[f64]) -> usize {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    for &k in knots {
        if !(prev == k) {
            n += 1;
            prev = k;
        }
    }
    n
}

/// `IntCurveSurface_TheHCurveTool::SamplePars` (`...cxx:72-299`). Returns the
/// sampled parameters in order (OCCT's 1-based `Pars(1..NbSamples)`).
pub fn sample_pars(curve: &dyn Curve, u0: f64, u1: f64, defl: f64, nb_min: usize) -> Vec<f64> {
    // `...cxx:80-99`: uniform counts for lines / Bezier / other curves.
    if curve.nurbs_degree().is_none() {
        let mut nbs = if curve.is_line() {
            2.0
        } else if let Some(poles) = curve.bezier_poles() {
            3.0 + poles.len() as f64
        } else {
            10.0
        };
        if nbs > 50.0 {
            nbs = 50.0;
        }
        let nnbs = nbs as usize;
        if nnbs < 2 {
            return vec![u0, u1];
        }
        let du = (u1 - u0) / (nnbs - 1) as f64;
        let mut pars = vec![0.0; nnbs];
        pars[0] = u0;
        pars[nnbs - 1] = u1;
        let mut u = u0 + du;
        for p in pars.iter_mut().take(nnbs - 1).skip(1) {
            *p = u;
            u += du;
        }
        return pars;
    }

    // `...cxx:115-196`: B-spline branch, sample per knot span at degree
    // subdivisions.
    let Some(knots) = curve.bspline_knots() else {
        return vec![u0, u1];
    };
    let degree = curve.nurbs_degree().unwrap_or(1);
    if knots.len() < degree + 1 || degree == 0 {
        return vec![u0, u1];
    }
    let mut ui1 = degree;
    let mut ui2 = knots.len() - 1 - degree;

    for i in ui1..ui2 {
        if u0 >= knots[i] && u0 < knots[i + 1] {
            ui1 = i;
            break;
        }
    }
    let mut i = ui2;
    while i > ui1 {
        if u1 <= knots[i] && u1 > knots[i - 1] {
            ui2 = i;
            break;
        }
        i -= 1;
    }

    let mut nbsu = ui2 - ui1 + 1;
    nbsu += (nbsu - 1) * (degree - 1);
    let mut b_uniform = false;
    if nbsu < nb_min {
        nbsu = nb_min;
        b_uniform = true;
    }

    // 1-based arrays to mirror the OCCT index arithmetic.
    let mut a_pars = vec![0.0; nbsu + 1];
    let mut a_flg = vec![false; nbsu + 1];
    if nbsu >= 2 {
        if b_uniform {
            let t1 = u0;
            let t2 = u1;
            let dt = (t2 - t1) / (nbsu - 1) as f64;
            a_pars[1] = t1;
            a_pars[nbsu] = t2;
            let mut t = t1 + dt;
            for p in a_pars.iter_mut().take(nbsu).skip(2) {
                *p = t;
                t += dt;
            }
        } else {
            let nbi = degree;
            let mut k = 0usize;
            let mut t1 = u0;
            for ii in (ui1 + 1)..=ui2 {
                let t2 = if ii == ui2 { u1 } else { knots[ii] };
                let dt = (t2 - t1) / nbi as f64;
                let mut j = 1usize;
                loop {
                    k += 1;
                    if k <= nbsu {
                        a_pars[k] = t1;
                        a_flg[k] = false;
                    }
                    t1 += dt;
                    j += 1;
                    if j > nbi {
                        break;
                    }
                }
                t1 = t2;
            }
            k += 1;
            if k <= nbsu {
                a_pars[k] = t1;
            }
        }
    } else {
        return vec![u0, u1];
    }

    // `...cxx:197-268`: deflection analysis over the raw samples.
    let a_defl2 = (defl * defl).max(1e-9);
    let tol = (0.01 * a_defl2).max(1e-9);

    let mut nb_samples = 2usize;
    a_flg[1] = true;
    a_flg[nbsu] = true;
    let mut j = 1usize;
    let mut b_cont = true;
    while j < nbsu - 1 && b_cont {
        if a_flg[j + 1] {
            j += 1;
            continue;
        }
        let p1 = curve.d0(a_pars[j]);
        let mut k = j + 2;
        let mut last_k = k;
        while k <= nbsu {
            last_k = k;
            let p2 = curve.d0(a_pars[k]);
            if p1.square_distance(&p2) <= tol {
                k += 1;
                continue;
            }
            let mut ok = true;
            if let Ok(dir) = occt_core::gp::GpDir::from_vec(&occt_core::gp::GpVec::from_pnts(&p1, &p2))
            {
                let lin = occt_core::gp::GpLin::from_pnt_dir(p1, dir);
                for l in (j + 1)..k {
                    if a_flg[l] {
                        ok = false;
                        break;
                    }
                    let pp = curve.d0(a_pars[l]);
                    let d = lin.square_distance(&pp);
                    if d <= a_defl2 {
                        continue;
                    }
                    ok = false;
                    break;
                }
            } else {
                ok = false;
            }
            if !ok {
                j = k - 1;
                a_flg[j] = true;
                nb_samples += 1;
                break;
            }
            if a_flg[k] {
                j = k;
                break;
            }
            k += 1;
        }
        if last_k >= nbsu {
            b_cont = false;
        }
    }

    // `...cxx:270-287`: too few kept points -> uniform fallback.
    if nb_samples < nb_min {
        let n = nb_min.max(2);
        let dt = (u1 - u0) / (n - 1) as f64;
        let mut pars = vec![0.0; n];
        pars[0] = u0;
        pars[n - 1] = u1;
        let mut t = u0 + dt;
        for p in pars.iter_mut().take(n - 1).skip(1) {
            *p = t;
            t += dt;
        }
        return pars;
    }

    // `...cxx:289-297`: collect the flagged samples.
    let mut pars = Vec::with_capacity(nb_samples);
    for i in 1..=nbsu {
        if a_flg[i] {
            pars.push(a_pars[i]);
        }
    }
    pars
}

/// `Adaptor3d_HSurfaceTool::NbSamplesU` (`...cxx:27-45`).
pub fn surface_nb_samples_u(surface: &dyn Surface) -> usize {
    if surface.gp_pln().is_some() {
        return 2;
    }
    if surface.is_bezier_surface() {
        return 3 + surface.nb_u_poles().max(0) as usize;
    }
    if surface.is_bspline_surface() {
        let nb_knots = surface
            .bspline_surface_uknots()
            .map_or(0, unique_knot_count);
        let nbs = nb_knots * surface.u_degree().max(0) as usize;
        return if nbs < 2 { 2 } else { nbs };
    }
    if surface.gp_torus().is_some() {
        return 20;
    }
    10
}

/// `Adaptor3d_HSurfaceTool::NbSamplesV` (`...cxx:47-70`).
pub fn surface_nb_samples_v(surface: &dyn Surface) -> usize {
    if surface.gp_pln().is_some() {
        return 2;
    }
    if surface.is_bezier_surface() {
        return 3 + surface.nb_v_poles().max(0) as usize;
    }
    if surface.is_bspline_surface() {
        let nb_knots = surface
            .bspline_surface_vknots()
            .map_or(0, unique_knot_count);
        let nbs = nb_knots * surface.v_degree().max(0) as usize;
        return if nbs < 2 { 2 } else { nbs };
    }
    if surface.gp_cylinder().is_some()
        || surface.gp_cone().is_some()
        || surface.gp_sphere().is_some()
        || surface.gp_torus().is_some()
        || surface.is_surface_of_revolution()
        || surface.is_surface_of_linear_extrusion()
    {
        return 15;
    }
    10
}

/// `Adaptor3d_HSurfaceTool::NbSamplesU(S, u1, u2)` (`...cxx:72-93`).
pub fn surface_nb_samples_u_range(surface: &dyn Surface, u1: f64, u2: f64) -> usize {
    let nbs = surface_nb_samples_u(surface);
    let mut n = nbs;
    if nbs > 10 {
        let (uf, ul) = surface.u_range();
        n = (n as f64 * ((u2 - u1) / (ul - uf))) as usize;
        if n > nbs || n > 50 {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}

/// `Adaptor3d_HSurfaceTool::NbSamplesV(S, v1, v2)` (`...cxx:95-116`).
pub fn surface_nb_samples_v_range(surface: &dyn Surface, v1: f64, v2: f64) -> usize {
    let nbs = surface_nb_samples_v(surface);
    let mut n = nbs;
    if nbs > 10 {
        let (vf, vl) = surface.v_range();
        n = (n as f64 * ((v2 - v1) / (vl - vf))) as usize;
        if n > nbs || n > 50 {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}
