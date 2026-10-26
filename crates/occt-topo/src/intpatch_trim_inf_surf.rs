//! `FUN_GetMinMaxXYZPnt` / `FUN_TrimInfSurf`
//! (`IntPatch_Intersection.cxx:276-448`).
//!
//! Both are `static` helpers of `IntPatch_Intersection.cxx` used by
//! `ParamParamPerfom` when exactly one of the two domains is infinite
//! (`:1687-1708`): the infinite surface is trimmed to a finite window derived
//! from the extrema of the point cloud of the finite one, because
//! `IntPatch_PrmPrmIntersection` cannot sample an infinite domain.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::Precision;
use occt_geom::extrema_surf::ExtPs;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::Surface;

/// `FUN_GetMinMaxXYZPnt(S, pMin, pMax)` (`cxx:278-304`).
///
/// Samples the surface every quarter of its parametric span and keeps the
/// points with the smallest / largest `XYZ().Modulus()`.
pub(crate) fn get_min_max_xyz_pnt(s: &dyn Surface) -> (GpPnt, GpPnt) {
    let (u_first, u_last) = s.u_range();
    let (v_first, v_last) = s.v_range();
    let du = 0.25 * (u_last - u_first).abs();
    let dv = 0.25 * (v_last - v_first).abs();
    let mut t_min = f64::MAX;
    let mut t_max = -t_min;
    let mut p_min = GpPnt::zero();
    let mut p_max = GpPnt::zero();

    // `for (double U = FirstU; U <= LastU; U += DU)` (`cxx:286-287`). A zero
    // step would spin forever, so it is only reachable for a degenerate range;
    // OCCT has the same exposure and relies on `DU > 0` for a real surface.
    if du <= 0.0 || dv <= 0.0 {
        return (p_min, p_max);
    }
    let mut u = u_first;
    while u <= u_last {
        let mut v = v_first;
        while v <= v_last {
            let p = s.d0(u, v);
            let c = p.xyz().modulus();
            if c > t_max {
                t_max = c;
                p_max = p;
            }
            if c < t_min {
                t_min = c;
                p_min = p;
            }
            v += dv;
        }
        u += du;
    }
    (p_min, p_max)
}

/// `FUN_TrimInfSurf(Pmin, Pmax, InfSurf, AlternativeTrimPrm, TrimS)`
/// (`cxx:310-448`).
///
/// `Extrema_ExtPS` over the (infinite) surface locates the parameters of the
/// two cloud corners; `TP` becomes the largest of their absolute parameters,
/// and the infinite directions are cut at `1000 * TP` (`:401`).
pub(crate) fn trim_inf_surf(
    p_min: &GpPnt,
    p_max: &GpPnt,
    inf_surf: &dyn Surface,
    alternative_trim_prm: f64,
) -> Arc<dyn Surface> {
    let (u0, u1) = inf_surf.u_range();
    let (v0, v1) = inf_surf.v_range();

    let mut ext = ExtPs::new();
    ext.initialize(inf_surf, u0, u1, v0, v1, 1.0e-7, 1.0e-7);

    let mut u_max = -1.0e100;
    let mut u_min = 1.0e100;
    let mut v_max = -1.0e100;
    let mut v_min = 1.0e100;

    // `// Process Pmin` (`cxx:327-354`).
    ext.perform(p_min);
    let is_done_1 = ext.is_done();
    if is_done_1 {
        for i in 1..=ext.nb_ext() {
            let (cu, cv, _) = ext.point(i);
            if cu > u_max {
                u_max = cu;
            }
            if cu < u_min {
                u_min = cu;
            }
            if cv > v_max {
                v_max = cv;
            }
            if cv < v_min {
                v_min = cv;
            }
        }
    }

    // `// Process Pmax` (`cxx:356-382`).
    ext.perform(p_max);
    let is_done_2 = ext.is_done();
    if is_done_2 {
        for i in 1..=ext.nb_ext() {
            let (cu, cv, _) = ext.point(i);
            if cu > u_max {
                u_max = cu;
            }
            if cu < u_min {
                u_min = cu;
            }
            if cv > v_max {
                v_max = cv;
            }
            if cv < v_min {
                v_min = cv;
            }
        }
    }

    let mut tp = alternative_trim_prm;
    if is_done_1 || is_done_2 {
        tp = u_min
            .abs()
            .max(u_max.abs().max(v_min.abs().max(v_max.abs())));
    }
    if tp == 0.0 {
        // `TrimS = InfSurf;` (`cxx:390-392`): the original handle is kept.
        return Arc::from(inf_surf.clone_dyn());
    }

    let u_inf = Precision::is_negative_infinite(u0);
    let u_sup = Precision::is_positive_infinite(u1);
    let v_inf = Precision::is_negative_infinite(v0);
    let v_sup = Precision::is_positive_infinite(v1);
    let tp1000 = 1000.0 * tp;

    // The V pass (`cxx:403-418`). The flags are independent `if`s in the
    // `.cxx`, so at most one of them fires for a well-formed range.
    let mut is_trimed = false;
    let mut trimmed: Option<Arc<dyn Surface>> = None;
    if v_inf && v_sup {
        trimmed = v_trim(inf_surf, -tp1000, tp1000);
        is_trimed = true;
    }
    if v_inf && !v_sup {
        trimmed = v_trim(inf_surf, -tp1000, v1);
        is_trimed = true;
    }
    if !v_inf && v_sup {
        trimmed = v_trim(inf_surf, v0, tp1000);
        is_trimed = true;
    }

    // The U pass (`cxx:419-448`): over the V-trimmed patch when one was made,
    // over the original surface otherwise.
    if is_trimed {
        let tmp = match trimmed {
            Some(t) => t,
            None => return Arc::from(inf_surf.clone_dyn()),
        };
        if u_inf && u_sup {
            if let Some(t) = u_trim(&*tmp, -tp1000, tp1000) {
                return t;
            }
        } else if u_inf && !u_sup {
            if let Some(t) = u_trim(&*tmp, -tp1000, u1) {
                return t;
            }
        } else if !u_inf && u_sup {
            if let Some(t) = u_trim(&*tmp, u0, tp1000) {
                return t;
            }
        }
        tmp
    } else if u_inf && u_sup {
        u_trim(inf_surf, -tp1000, tp1000).unwrap_or_else(|| Arc::from(inf_surf.clone_dyn()))
    } else if u_inf && !u_sup {
        u_trim(inf_surf, -tp1000, u1).unwrap_or_else(|| Arc::from(inf_surf.clone_dyn()))
    } else if !u_inf && u_sup {
        u_trim(inf_surf, u0, tp1000).unwrap_or_else(|| Arc::from(inf_surf.clone_dyn()))
    } else {
        Arc::from(inf_surf.clone_dyn())
    }
}

/// `GeomAdaptor_Surface::UTrim(First, Last, Tol)` (`GeomAdaptor_Surface.cxx`):
/// `new Geom_RectangularTrimmedSurface(mySurf, First, Last, true)` wrapped in
/// an adaptor. The trimmed surface carries the range, so the port returns it
/// directly.
fn u_trim(s: &dyn Surface, first: f64, last: f64) -> Option<Arc<dyn Surface>> {
    GeomRectangularTrimmedSurface::one_param(Arc::from(s.clone_dyn()), first, last, true)
        .ok()
        .map(|t| Arc::new(t) as Arc<dyn Surface>)
}

/// `GeomAdaptor_Surface::VTrim(First, Last, Tol)`: the V mirror of
/// [`u_trim`].
fn v_trim(s: &dyn Surface, first: f64, last: f64) -> Option<Arc<dyn Surface>> {
    GeomRectangularTrimmedSurface::one_param(Arc::from(s.clone_dyn()), first, last, false)
        .ok()
        .map(|t| Arc::new(t) as Arc<dyn Surface>)
}
