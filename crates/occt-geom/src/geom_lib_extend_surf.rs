//! `GeomLib::ExtendSurfByLength` and its private helper `ExtendKPart`.
//!
//! Source: `ModelingData/TKGeomBase/GeomLib/GeomLib.cxx:1416-1975`
//! (OCCT 8.0.0).
//!
//! Companions ported elsewhere:
//! * `BSplCLib::TangExtendToConstraint` (`BSplCLib.cxx:3794-4238`) is
//!   [`crate::bspl_tang_extend::tang_extend_to_constraint`],
//! * `BSplCLib::Eval` (`cxx:3640-3790`) is
//!   [`occt_core::bspl::eval_extrap::eval_dim`],
//! * `BSplCLib::Knots` (`cxx:572-598`) is
//!   [`occt_core::bspl::knots::unique_knots_mults`],
//! * `BSplCLib::KnotSequence` (`cxx:488-516`) is
//!   [`occt_core::bspl::banded_interp::knot_sequence`],
//! * `Geom_BSplineSurface::IncreaseDegree` / `Segment` /
//!   `KnotSequence` are the `GeomBSplineSurface` methods of the same names.
//!
//! UNPORTED (`GeomLib.cxx:1509-1527`): when the incoming surface is neither a
//! `Geom_BSplineSurface` nor a `Geom_RectangularTrimmedSurface` over one, OCCT
//! first runs `GeomConvert_ApproxSurface` (`Tol = Precision::Confusion()`,
//! `C1/C1`, degrees 14/14, segments 16, prec 1) and falls back to
//! `GeomConvert::SurfaceToBSplineSurface`. Neither converter is available in
//! `occt-geom`, so [`extend_surf_by_length`] reports that arm as `None`.

use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::eval_extrap::eval_dim;
use occt_core::bspl::knots::unique_knots_mults;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::PCONFUSION;

use crate::bspl_tang_extend::tang_extend_to_constraint;
use crate::bspline_surface::GeomBSplineSurface;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;

/// `ExtendKPart(Surface, Length, InU, After)` (`GeomLib.cxx:1416-1483`).
///
/// The OCCT parameter is the `Geom_RectangularTrimmedSurface` down-cast of the
/// incoming surface; a null handle is the "not a trimmed patch" answer, which
/// the port expresses with `rectangular_trimmed_basis() == None`. On success
/// OCCT replaces the handle with a *new* trimmed patch over the same basis
/// surface and the shifted bounds (`cxx:1449`, `cxx:1472`).
fn extend_k_part(
    surface: &Arc<dyn Surface>,
    length: f64,
    in_u: bool,
    after: bool,
) -> Option<Arc<dyn Surface>> {
    let support = surface.rectangular_trimmed_basis()?;
    let (uf, ul) = surface.u_range();
    let (vf, vl) = surface.v_range();

    // `GeomAdaptor_Surface AS(Surface); Type = AS.GetType()` (`cxx:1433-1434`).
    // The adaptor of a trimmed patch reports the *basis* type, which is what
    // the `gp_*` hooks of `GeomRectangularTrimmedSurface` forward.
    let supported = if in_u {
        // `cxx:1436-1456`: planes only.
        surface.gp_pln().is_some()
    } else {
        // `cxx:1457-1479`: planes, cylinders and surfaces of linear extrusion.
        surface.gp_pln().is_some()
            || surface.gp_cylinder().is_some()
            || surface.is_surface_of_linear_extrusion()
    };
    if !supported {
        return None;
    }

    let (uf, ul, vf, vl) = if in_u {
        if after {
            (uf, ul + length, vf, vl)
        } else {
            (uf - length, ul, vf, vl)
        }
    } else if after {
        (uf, ul, vf, vl + length)
    } else {
        (uf, ul, vf - length, vl)
    };
    GeomRectangularTrimmedSurface::uv(support, uf, ul, vf, vl)
        .ok()
        .map(|t| Arc::new(t) as Arc<dyn Surface>)
}

/// The `Geom_BSplineSurface` OCCT extends in place: the surface itself when it
/// is one, else the basis of the trimmed patch (`GeomLib.cxx:1505-1508`).
fn bspline_target(surface: &Arc<dyn Surface>) -> Option<GeomBSplineSurface> {
    if surface.is_bspline_surface() {
        return surface.osculating_bspline();
    }
    let basis = surface.rectangular_trimmed_basis()?;
    if basis.is_bspline_surface() {
        basis.osculating_bspline()
    } else {
        None
    }
}

/// `GeomLib::ExtendSurfByLength(Surface, Length, Continuity, InU, After)`
/// (`GeomLib.cxx:1485-1975`).
///
/// `None` is used for the two arms the port cannot reproduce: the
/// `Continuity` guard is a no-op there too (`cxx:1492-1495` returns the
/// unchanged surface, so the caller gets a clone), while the
/// `GeomConvert` conversion (`cxx:1509-1527`) and any construction failure are
/// reported as `None`.
pub fn extend_surf_by_length(
    surface: &Arc<dyn Surface>,
    length: f64,
    continuity: i32,
    in_u: bool,
    after: bool,
) -> Option<Arc<dyn Surface>> {
    // `if (Continuity < 0 || Continuity > 3) return;` (`cxx:1492-1495`).
    if !(0..=3).contains(&continuity) {
        return Some(surface.clone());
    }

    // Kpart ? (`cxx:1497-1503`).
    if let Some(trimmed) = extend_k_part(surface, length, in_u, after) {
        return Some(trimmed);
    }

    // The b-spline the extension works on (`cxx:1505-1527`).
    let bs = bspline_target(surface)?;
    let extended = extend_bspline_surface(bs, length, continuity, in_u, after)?;
    Some(Arc::new(extended) as Arc<dyn Surface>)
}

/// The b-spline part of `GeomLib::ExtendSurfByLength`
/// (`GeomLib.cxx:1528-1975`), from the degree elevation down to the
/// `new Geom_BSplineSurface(...)` at the end.
///
/// `bs` is the surface OCCT mutates in place (`IncreaseDegree`, `Segment`);
/// the returned surface is OCCT's `Res`, which the caller installs. `None`
/// covers the arms whose OCCT exception leaves no result.
pub fn extend_bspline_surface(
    mut bs: GeomBSplineSurface,
    length: f64,
    continuity: i32,
    in_u: bool,
    after: bool,
) -> Option<GeomBSplineSurface> {
    // Degree sufficient for the requested continuity (`cxx:1528-1536`).
    if in_u && bs.deg_u() < continuity as usize + 1 {
        bs.increase_degree(continuity as usize + 1, bs.deg_v()).ok()?;
    }
    if !in_u && bs.deg_v() < continuity as usize + 1 {
        bs.increase_degree(bs.deg_u(), continuity as usize + 1).ok()?;
    }

    // If BS was periodic in the extension direction, it will no longer be
    // (`cxx:1538-1543`).
    if (in_u && bs.is_u_periodic()) || (!in_u && bs.is_v_periodic()) {
        let (u0, u1) = bs.u_range();
        let (v0, v1) = bs.v_range();
        bs.segment_uv(u0, u1, v0, v1, 0.0, 0.0).ok()?;
    }

    // `rational` covers both directions (OCC bug 0022694, `cxx:1545-1552`).
    let rational = bs.is_u_rational() || bs.is_v_rational();
    let eps_w = 10.0 * PCONFUSION;
    let gap = if rational { 4usize } else { 3usize };

    let mut cdeg = 0usize;
    let mut cdim = 0usize;
    let mut nbp = 0usize;
    let mut psize = 1usize;
    let mut fknots: Vec<f64> = Vec::new();
    let mut poles: Vec<f64> = Vec::new();
    let mut point: Vec<f64> = Vec::new();
    let mut tgte: Vec<f64> = Vec::new();
    let mut lambda: Vec<f64> = Vec::new();
    let mut lambmin = length;
    let mut ok = false;

    // `for (Kount = 0, Ok = false; Kount <= 2 && !Ok; Kount++)`
    // (`cxx:1562-1753`).
    let mut kount = 0i32;
    while kount <= 2 && !ok {
        // The surface written as a one-variable b-spline (`cxx:1570-1581`).
        if in_u {
            cdeg = bs.deg_u();
            nbp = bs.nb_poles_u();
            cdim = bs.nb_poles_v() * gap;
        } else {
            cdeg = bs.deg_v();
            nbp = bs.nb_poles_v();
            cdim = bs.nb_poles_u() * gap;
        }

        // The flat knots (`cxx:1583-1592`).
        fknots = if in_u {
            bs.knots_u().to_vec()
        } else {
            bs.knots_v().to_vec()
        };

        // The parameter of the connection knot (`cxx:1594-1603`).
        let tbord = if after {
            *fknots.get(fknots.len() - 1 - cdeg)?
        } else {
            *fknots.get(cdeg)?
        };

        // The poles (`cxx:1605-1642`).
        psize = cdim * nbp;
        poles = vec![0.0; psize];
        let mut ipole = 0usize;
        if in_u {
            for ii in 0..nbp {
                for jj in 0..bs.nb_poles_v() {
                    let p = bs.pole(ii, jj);
                    poles[ipole] = p.x();
                    poles[ipole + 1] = p.y();
                    poles[ipole + 2] = p.z();
                    if rational {
                        poles[ipole + 3] = bs.weight(ii, jj);
                    }
                    ipole += gap;
                }
            }
        } else {
            for jj in 0..nbp {
                for ii in 0..bs.nb_poles_u() {
                    let p = bs.pole(ii, jj);
                    poles[ipole] = p.x();
                    poles[ipole + 1] = p.y();
                    poles[ipole + 2] = p.z();
                    if rational {
                        poles[ipole + 3] = bs.weight(ii, jj);
                    }
                    ipole += gap;
                }
            }
        }

        // `BSplCLib::Eval(Tbord, periodic_flag = false, derivative_request,
        // extrap_mode[0] = Cdeg, Cdeg, FKnots, Cdim, *Padr, *Radr)`
        // (`cxx:1660-1677`).
        let derivative_request = continuity.max(1);
        let res = eval_dim(
            tbord,
            false,
            derivative_request,
            [cdeg as i32; 2],
            cdeg as i32,
            &fknots,
            cdim,
            &poles,
        )?;
        ok = true;
        point = res[0..cdim].to_vec();
        tgte = res[cdim..2 * cdim].to_vec();
        lambda = vec![0.0f64; cdim];

        // Constraint to reach (`cxx:1679-1741`).
        let tgtol = 1.0e-12;
        let mut old_t = GpVec::new(0.0, 0.0, 0.0);
        let mut old_n = 0.0;
        if rational {
            for k in 0..cdim / gap {
                tgte[k * gap + 3] = 0.0;
            }
            for k in 0..cdim / gap {
                let base = k * gap;
                let cur_t = GpVec::new(tgte[base], tgte[base + 1], tgte[base + 2]);
                let n_tgte = cur_t.magnitude();
                if n_tgte > tgtol {
                    let val = length / n_tgte;
                    // The segment given by the poles can run opposite to the
                    // derivative (parts of a torus, `cxx:1694-1700`).
                    if old_n > tgtol && cur_t.angle(&old_t) > 2.0 {
                        ok = false;
                    }
                    lambda[base] = val;
                    lambda[base + 1] = val;
                    lambda[base + 2] = val;
                    lambda[base + 3] = 0.0;
                    lambmin = lambmin.min(val);
                } else {
                    lambda[base] = 0.0;
                    lambda[base + 1] = 0.0;
                    lambda[base + 2] = 0.0;
                    lambda[base + 3] = 0.0;
                }
                old_t = cur_t;
                old_n = n_tgte;
            }
        } else {
            for k in 0..cdim / gap {
                let base = k * gap;
                let cur_t = GpVec::new(tgte[base], tgte[base + 1], tgte[base + 2]);
                let n_tgte = cur_t.magnitude();
                if n_tgte > tgtol {
                    let val = length / n_tgte;
                    if old_n > tgtol && cur_t.angle(&old_t) > 2.0 {
                        ok = false;
                    }
                    lambda[base] = val;
                    lambda[base + 1] = val;
                    lambda[base + 2] = val;
                    lambmin = lambmin.min(val);
                } else {
                    lambda[base] = 0.0;
                    lambda[base + 1] = 0.0;
                    lambda[base + 2] = 0.0;
                }
                old_t = cur_t;
                old_n = n_tgte;
            }
        }

        if !ok && kount < 2 {
            // Bring the poles of the border iso closer and retry
            // (`cxx:1742-1752`).
            if in_u {
                bs.increase_degree(bs.deg_u(), bs.deg_v() + 2).ok()?;
            } else {
                bs.increase_degree(bs.deg_u() + 2, bs.deg_v()).ok()?;
            }
        }
        kount += 1;
    }

    // The constraint point (`cxx:1755-1771`).
    let mut constraint = vec![0.0f64; cdim];
    for ii in 0..cdim {
        constraint[ii] = if after {
            point[ii] + lambda[ii] * tgte[ii]
        } else {
            point[ii] - lambda[ii] * tgte[ii]
        };
    }

    // Special case of rational (`cxx:1773-1788`).
    if rational {
        let mut ipole = 0usize;
        while ipole < psize {
            poles[ipole] *= poles[ipole + 3];
            poles[ipole + 1] *= poles[ipole + 3];
            poles[ipole + 2] *= poles[ipole + 3];
            ipole += gap;
        }
        for ii in (0..cdim).step_by(gap) {
            constraint[ii] *= constraint[ii + 3];
            constraint[ii + 1] *= constraint[ii + 3];
            constraint[ii + 2] *= constraint[ii + 3];
        }
    }

    // The extension (`cxx:1790-1918`).
    let mut ext_ok = false;
    let mut fk_res: Vec<f64> = Vec::new();
    let mut new_poles: Vec<Vec<GpPnt>> = Vec::new();
    let mut new_weights: Vec<Vec<f64>> = Vec::new();

    let mut kount = 1i32;
    while kount <= 5 && !ext_ok {
        let ext = tang_extend_to_constraint(
            &fknots,
            lambmin,
            nbp,
            &poles,
            cdim,
            cdeg,
            &constraint,
            continuity as usize,
            after,
        )?;
        // `FKRes(1..NbKnots) = FK(1..NbKnots)` (`cxx:1934-1938`): `flat_knots`
        // is exactly the `NbKnots` values the extension reported.
        fk_res = ext.flat_knots;
        let pres = ext.poles;

        // Copy the result poles as 3D points and weights (`cxx:1821-1906`).
        let (nu, nv) = if in_u {
            (ext.nb_poles, bs.nb_poles_v())
        } else {
            (bs.nb_poles_u(), ext.nb_poles)
        };
        let mut np = vec![vec![GpPnt::zero(); nv]; nu];
        // `NewW.Init(1.)` for the non-rational surface (`cxx:1839-1842`).
        let mut nw = vec![vec![1.0f64; nv]; nu];
        let mut null_weight = false;
        if in_u {
            'outer_u: for ii in 0..nu {
                for jj in 0..nv {
                    let indice = ii * cdim + jj * gap;
                    np[ii][jj] = GpPnt::new(pres[indice], pres[indice + 1], pres[indice + 2]);
                    if rational {
                        let mut ww = pres[indice + 3];
                        if (ww - 1.0).abs() < eps_w {
                            ww = 1.0;
                        }
                        if ww < eps_w {
                            null_weight = true;
                            break 'outer_u;
                        }
                        nw[ii][jj] = ww;
                        let q = np[ii][jj];
                        np[ii][jj] = GpPnt::new(q.x() / ww, q.y() / ww, q.z() / ww);
                    }
                }
            }
        } else {
            'outer_v: for jj in 0..nv {
                for ii in 0..nu {
                    let indice = ii * gap + jj * cdim;
                    np[ii][jj] = GpPnt::new(pres[indice], pres[indice + 1], pres[indice + 2]);
                    if rational {
                        let mut ww = pres[indice + 3];
                        if (ww - 1.0).abs() < eps_w {
                            ww = 1.0;
                        }
                        if ww < eps_w {
                            null_weight = true;
                            break 'outer_v;
                        }
                        nw[ii][jj] = ww;
                        let q = np[ii][jj];
                        np[ii][jj] = GpPnt::new(q.x() / ww, q.y() / ww, q.z() / ww);
                    }
                }
            }
        }

        new_poles = np;
        new_weights = nw;
        if null_weight {
            lambmin /= 3.0;
        } else {
            ext_ok = true;
        }
        kount += 1;
    }

    // Flat knots -> knots with multiplicities, and the degrees of the result
    // (`cxx:1920-1958`).
    let mut usize_ = bs.uknots().len();
    let mut vsize = bs.vknots().len();
    if in_u {
        usize_ += 1;
    } else {
        vsize += 1;
    }
    let (udeg, vdeg, uknots, umults, vknots, vmults) = if in_u {
        let (uk, mut um) = unique_knots_mults(&fk_res);
        // `UMults(Usize) = UDeg + 1` (`cxx:1945`). OCCT writes the slot past
        // the derived multiplicities when the extension did not add a distinct
        // knot; that slot is uninitialised storage there, so the port keeps
        // the derived set in that case.
        if um.len() == usize_ {
            if let Some(last) = um.last_mut() {
                *last = cdeg as i32 + 1;
            }
        }
        (
            cdeg,
            bs.deg_v(),
            uk,
            um,
            bs.vknots().to_vec(),
            bs.vmults().to_vec(),
        )
    } else {
        let (vk, mut vm) = unique_knots_mults(&fk_res);
        // `VMults(Vsize) = VDeg + 1` (`cxx:1955`), see the U arm above.
        if vm.len() == vsize {
            if let Some(last) = vm.last_mut() {
                *last = cdeg as i32 + 1;
            }
        }
        (
            bs.deg_u(),
            cdeg,
            bs.uknots().to_vec(),
            bs.umults().to_vec(),
            vk,
            vm,
        )
    };

    // The knot sequence of the extended direction. The direction is
    // non-periodic (see the `Segment` above), the other one keeps the flat
    // sequence the source b-spline stores.
    let flat_u = if in_u {
        knot_sequence(&uknots, &umults, udeg as i32)
    } else {
        bs.knots_u().to_vec()
    };
    let flat_v = if in_u {
        bs.knots_v().to_vec()
    } else {
        knot_sequence(&vknots, &vmults, vdeg as i32)
    };

    // `new Geom_BSplineSurface(NewPoles, NewWeights, UKnots, VKnots, UMults,
    // VMults, UDeg, VDeg, IsUPeriodic, IsVPeriodic)` (`cxx:1961-1972`).
    //
    // OCCT always hands the weight grid over; all weights are 1 for the
    // non-rational case, which is exactly the port's `weights == None` state.
    let weights = if rational { Some(new_weights) } else { None };
    Some(GeomBSplineSurface::from_parts(
        new_poles,
        weights,
        flat_u,
        flat_v,
        udeg,
        vdeg,
        bs.is_u_periodic(),
        bs.is_v_periodic(),
    ))
}
