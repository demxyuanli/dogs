//! Port of ShapeCustom_Surface (ShapeCustom_Surface.cxx) — only the
//! ConvertToPeriodic entry the STEP face path uses
//! (StepToTopoDS_TranslateFace.cxx:557-567).
//!
//! OCCT forces every B-spline STEP face surface that is geometrically closed
//! but not parametrically periodic into the periodic representation. Without it
//! a pcurve that runs once around the closed direction leaves the surface's
//! parameter window (the ported projection removes the seam jump by one
//! period), so evaluating it extrapolates and ShapeFix_Edge::FixSameParameter
//! records a huge deviation as a vertex tolerance. See
//! specs/_a3n00_gap_analysis.md 3.5.

use std::sync::Arc;

use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::Surface;

use crate::pcurve_full::{sa_is_u_closed, sa_is_v_closed};

/// Split a flat knot vector into (unique knots, multiplicities).
fn unique_km(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    GeomBSplineSurface::unique_knots_mults(flat)
}

/// ShapeCustom_Surface::ConvertToPeriodic(substitute=false, preci)
/// (ShapeCustom_Surface.cxx:475-616). Returns None when nothing changed
/// (cxx:606-609) or the surface is not an exact Geom_BSplineSurface.
pub fn convert_to_periodic(surf: &Arc<dyn Surface>, preci: f64) -> Option<Arc<dyn Surface>> {
    // down_cast<Geom_BSplineSurface>(mySurf) (cxx:479-483).
    if !surf.is_bspline_surface() {
        return None;
    }
    let mut bs = surf.osculating_bspline()?;

    // ShapeAnalysis_Surface sas(mySurf); sas.IsU/VClosed(preci) (cxx:485-492).
    let u_closed = sa_is_u_closed(surf.as_ref(), preci);
    let v_closed = sa_is_v_closed(surf.as_ref(), preci);
    if !u_closed && !v_closed {
        return None;
    }

    let mut converted = false;

    if u_closed && !bs.is_u_periodic() && bs.nb_poles_u() > 3 {
        let mut set = true;
        let (uk, um) = unique_km(&bs.knots_u);
        // cxx:500-541: when both end multiplicities are degree + 1, first
        // rearrange the knots so they become 1, then SetUPeriodic.
        if um.len() >= 2 && um[0] == bs.deg_u as i32 + 1 && um[um.len() - 1] == bs.deg_u as i32 + 1 {
            let n = uk.len();
            let a = 0.5 * ((uk[1] - uk[0]) + (uk[n - 1] - uk[n - 2]));
            let mut new_knots = vec![0.0f64; n + 2];
            let mut new_mults = vec![0i32; n + 2];
            new_knots[0] = uk[0] - a;
            new_knots[n + 1] = uk[n - 1] + a;
            new_mults[0] = 1;
            new_mults[n + 1] = 1;
            for i in 1..=n {
                new_knots[i] = uk[i - 1];
                new_mults[i] = um[i - 1];
            }
            new_mults[1] = bs.deg_u as i32;
            new_mults[n] = bs.deg_u as i32;
            let (vk, vm) = unique_km(&bs.knots_v);
            bs = rebuild(&bs, new_knots, new_mults, vk, vm)?;
        } else if um.len() >= 2
            && (um[0] > bs.deg_u as i32 || um[um.len() - 1] > bs.deg_u as i32 + 1)
        {
            set = false;
        }
        if set {
            bs.set_u_periodic();
            converted = true;
        }
    }

    if v_closed && !bs.is_v_periodic() && bs.poles.first().map_or(0, |r| r.len()) > 3 {
        let mut set = true;
        let (vk, vm) = unique_km(&bs.knots_v);
        // cxx:553-594.
        if vm.len() >= 2 && vm[0] == bs.deg_v as i32 + 1 && vm[vm.len() - 1] == bs.deg_v as i32 + 1 {
            let n = vk.len();
            let a = 0.5 * ((vk[1] - vk[0]) + (vk[n - 1] - vk[n - 2]));
            let mut new_knots = vec![0.0f64; n + 2];
            let mut new_mults = vec![0i32; n + 2];
            new_knots[0] = vk[0] - a;
            new_knots[n + 1] = vk[n - 1] + a;
            new_mults[0] = 1;
            new_mults[n + 1] = 1;
            for i in 1..=n {
                new_knots[i] = vk[i - 1];
                new_mults[i] = vm[i - 1];
            }
            new_mults[1] = bs.deg_v as i32;
            new_mults[n] = bs.deg_v as i32;
            let (uk2, um2) = unique_km(&bs.knots_u);
            bs = rebuild(&bs, uk2, um2, new_knots, new_mults)?;
        } else if vm.len() >= 2
            && (vm[0] > bs.deg_v as i32 || vm[vm.len() - 1] > bs.deg_v as i32 + 1)
        {
            set = false;
        }
        if set {
            bs.set_v_periodic();
            converted = true;
        }
    }

    if !converted {
        return None;
    }
    Some(Arc::new(bs))
}

/// Rebuild with new unique knot/multiplicity arrays in both directions; poles,
/// weights and degrees are kept (cxx:525-535 / 578-587).
fn rebuild(
    bs: &GeomBSplineSurface,
    u_knots: Vec<f64>,
    u_mults: Vec<i32>,
    v_knots: Vec<f64>,
    v_mults: Vec<i32>,
) -> Option<GeomBSplineSurface> {
    let ku = occt_core::bspl::banded_interp::knot_sequence(&u_knots, &u_mults, bs.deg_u as i32);
    let kv = occt_core::bspl::banded_interp::knot_sequence(&v_knots, &v_mults, bs.deg_v as i32);
    match &bs.weights {
        Some(w) => {
            GeomBSplineSurface::rational(bs.poles.clone(), w.clone(), ku, kv, bs.deg_u, bs.deg_v)
                .ok()
        }
        None => GeomBSplineSurface::new(bs.poles.clone(), ku, kv, bs.deg_u, bs.deg_v).ok(),
    }
}
