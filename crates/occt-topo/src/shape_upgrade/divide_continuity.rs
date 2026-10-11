//! `ShapeAlgo_AlgoContainer::C0ShapeToC1Shape`
//! (`ShapeAlgo_AlgoContainer.cxx:390-399`) for one face.
//!
//! `ShapeUpgrade_ShapeDivideContinuity::GetSplitFaceTool` installs a C1
//! `ShapeUpgrade_SplitSurfaceContinuity`. `FaceDivide::SplitSurface`
//! (`ShapeUpgrade_FaceDivide.cxx:98-202`) runs that tool, then
//! `ShapeFix_ComposeShell` when the surface changed.
//!
//! Knot removal that leaves a single patch (`ShapeUpgrade_SplitSurfaceContinuity.cxx:269-401`
//! and the swept-basis arm `:113-178`) is applied here. Extra split values
//! go through `ShapeUpgrade_SplitSurface::Build` (`ShapeUpgrade_SplitSurface.cxx:240`)
//! and `ShapeFix_ComposeShell` (`ShapeUpgrade_FaceDivide.cxx:181-198`).
//!
//! `FaceDivide::SplitCurves` (`ShapeUpgrade_FaceDivide.cxx:206-253`) is not
//! called. Entity 140 rebuilds the offset with `BRepLib_MakeFace` on the
//! corrected surface (`IGESToBRep_TopoSurface.cxx:1183-1194`), so the wires
//! on the intermediate face are not the result.

use std::sync::Arc;

use occt_core::precision::{Precision, PCONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::{Curve, GeomSurfaceOfLinearExtrusion, GeomSurfaceOfRevolution, Surface};

use super::split_surface::{split_bspline, split_extrusion, split_revolution, C0Basis};
use crate::brep_tools;
use crate::shape::Face;
use crate::topo_tools_full::edges_of;

const CONT_C1: i32 = 1;

/// Correct a C0 basis the way `C0ShapeToC1Shape` does. A single patch is the
/// corrected surface. A patch grid is composed; a shell is returned as a shape
/// because `TransferOffsetSurface` does not offset a non-face
/// (`IGESToBRep_TopoSurface.cxx:1136-1140`).
pub fn correct_c0_basis(
    face: &Face,
    surf: Arc<dyn Surface>,
    tolerance: f64,
) -> Result<C0Basis, String> {
    let (mut uf, mut ul, mut vf, mut vl) = face_uv_bounds(face, surf.as_ref());
    if Precision::is_infinite(uf)
        || Precision::is_infinite(ul)
        || Precision::is_infinite(vf)
        || Precision::is_infinite(vl)
    {
        // `ShapeUpgrade_FaceDivide::SplitSurface` (`cxx:121-125`) returns
        // before `Compute`. The face is unchanged.
        return Ok(C0Basis::Surface(surf));
    }
    extend_to_surface(surf.as_ref(), &mut uf, &mut ul, &mut vf, &mut vl);
    let (uf, ul, vf, vl) = clamp_init(surf.as_ref(), uf, ul, vf, vl);

    if let Some(basis) = surf.rectangular_trimmed_basis() {
        // `ShapeUpgrade_SplitSurface::Build` unwraps a trim and splits the
        // basis (`ShapeUpgrade_SplitSurface.cxx:413-426`).
        return correct_c0_basis(face, basis, tolerance);
    }
    if surf.is_bspline_surface() {
        let Some(mut bs) = surf.osculating_bspline() else {
            return Ok(C0Basis::Surface(surf));
        };
        let mut u_splits = vec![uf, ul];
        let mut v_splits = vec![vf, vl];
        raise_surface_knots(&mut bs, &mut u_splits, &mut v_splits, CONT_C1, tolerance);
        if u_splits.len() > 2 || v_splits.len() > 2 {
            return split_bspline(face, &bs, &mut u_splits, &mut v_splits);
        }
        return Ok(C0Basis::Surface(Arc::new(bs)));
    }

    if surf.is_surface_of_revolution() {
        return correct_revolution(face, surf, vf, vl, tolerance);
    }
    if surf.is_surface_of_linear_extrusion() {
        return correct_extrusion(face, surf, uf, ul, tolerance);
    }
    Ok(C0Basis::Surface(surf))
}

fn correct_revolution(
    face: &Face,
    surf: Arc<dyn Surface>,
    v_first: f64,
    v_last: f64,
    tolerance: f64,
) -> Result<C0Basis, String> {
    let Some(basis) = surf.revolution_basis_curve() else {
        return Ok(C0Basis::Surface(surf));
    };
    let Some(axis) = surf.revolution_axis() else {
        return Ok(C0Basis::Surface(surf));
    };
    let Some(mut curve) = clone_bspline(basis.as_ref()) else {
        return Ok(C0Basis::Surface(surf));
    };
    let mut splits = vec![v_first, v_last];
    // `ShapeUpgrade_SplitSurfaceContinuity` revolution arm (`cxx:113-143`).
    raise_curve_knots(&mut curve, &mut splits, CONT_C1, tolerance);
    if splits.len() > 2 {
        let (u0, u1) = surf.u_range();
        return split_revolution(face, &curve, axis, u0, u1, &mut splits);
    }
    Ok(C0Basis::Surface(Arc::new(GeomSurfaceOfRevolution::new(
        Arc::new(curve),
        axis,
    ))))
}

fn correct_extrusion(
    face: &Face,
    surf: Arc<dyn Surface>,
    u_first: f64,
    u_last: f64,
    tolerance: f64,
) -> Result<C0Basis, String> {
    let Some(basis) = surf.extrusion_basis_curve() else {
        return Ok(C0Basis::Surface(surf));
    };
    let Some(dir) = surf.extrusion_direction() else {
        return Ok(C0Basis::Surface(surf));
    };
    let Some(mut curve) = clone_bspline(basis.as_ref()) else {
        return Ok(C0Basis::Surface(surf));
    };
    let mut splits = vec![u_first, u_last];
    // `ShapeUpgrade_SplitSurfaceContinuity` extrusion arm (`cxx:145-177`).
    raise_curve_knots(&mut curve, &mut splits, CONT_C1, tolerance);
    if splits.len() > 2 {
        let (v0, v1) = surf.v_range();
        return split_extrusion(face, &curve, dir, v0, v1, &mut splits);
    }
    Ok(C0Basis::Surface(Arc::new(GeomSurfaceOfLinearExtrusion::new(
        Arc::new(curve),
        dir,
    ))))
}

/// `ShapeAnalysis::GetFaceUVBounds` (`ShapeAnalysis.cxx:268-299`).
///
/// A face with no edges uses the surface bounds (`cxx:277-281`). A face with
/// edges uses the pcurve box. `FillBndBox` (`cxx:296`) is the 20-sample
/// curve box; the port uses `BRepTools::UVBounds` (`brep_tools::uv_bounds`).
fn face_uv_bounds(face: &Face, surf: &dyn Surface) -> (f64, f64, f64, f64) {
    if edges_of(&face.0).is_empty() {
        let (u0, u1) = surf.u_range();
        let (v0, v1) = surf.v_range();
        (u0, u1, v0, v1)
    } else {
        brep_tools::uv_bounds(face)
    }
}

/// The 1% extension in `ShapeUpgrade_FaceDivide::SplitSurface` (`cxx:127-153`).
fn extend_to_surface(surf: &dyn Surface, uf: &mut f64, ul: &mut f64, vf: &mut f64, vl: &mut f64) {
    let (su0, su1) = surf.u_range();
    let (sv0, sv1) = surf.v_range();
    if !surf.is_u_periodic() {
        let du = (*ul - *uf) * 0.01;
        if *uf > su0 {
            *uf -= du.min(*uf - su0);
        }
        if *ul < su1 {
            *ul += du.min(su1 - *ul);
        }
    }
    if !surf.is_v_periodic() {
        let dv = (*vl - *vf) * 0.01;
        if *vf > sv0 {
            *vf -= dv.min(*vf - sv0);
        }
        if *vl < sv1 {
            *vl += dv.min(sv1 - *vl);
        }
    }
}

/// `ShapeUpgrade_SplitSurface::Init` clamp (`cxx:101-134`). `theArea` is 0,
/// so the area branch (`cxx:136`) is not entered.
fn clamp_init(surf: &dyn Surface, u_first: f64, u_last: f64, v_first: f64, v_last: f64) -> (f64, f64, f64, f64) {
    let (mut u1, mut u2) = surf.u_range();
    let (mut v1, mut v2) = surf.v_range();
    let precision = PCONFUSION;
    if surf.is_u_periodic() && u_last - u_first <= u2 - u1 + precision {
        u1 = u_first;
        u2 = u1 + surf.u_period();
    }
    if surf.is_v_periodic() && v_last - v_first <= v2 - v1 + precision {
        v1 = v_first;
        v2 = v1 + surf.v_period();
    }
    let (uf, ul) = if u_first > u2 - precision || u_last < u1 - precision {
        (u1, u2)
    } else {
        (u1.max(u_first), u2.min(u_last))
    };
    let (vf, vl) = if v_first > v2 - precision || v_last < v1 - precision {
        (v1, v2)
    } else {
        (v1.max(v_first), v2.min(v_last))
    };
    (uf, ul, vf, vl)
}

/// U then V interior knots of `ShapeUpgrade_SplitSurfaceContinuity::Compute`
/// (`cxx:269-401`). `myCont` is 1 (`SetCriterion(GeomAbs_C1)`).
fn raise_surface_knots(
    bs: &mut GeomBSplineSurface,
    u_splits: &mut Vec<f64>,
    v_splits: &mut Vec<f64>,
    cont: i32,
    tolerance: f64,
) {
    raise_surface_dir(bs, u_splits, true, cont, tolerance);
    raise_surface_dir(bs, v_splits, false, cont, tolerance);
}

fn raise_surface_dir(
    bs: &mut GeomBSplineSurface,
    splits: &mut Vec<f64>,
    u_direction: bool,
    cont: i32,
    tolerance: f64,
) {
    let precision = PCONFUSION;
    let (knots, mults) = if u_direction {
        bs.distinct_knots_and_mults_u()
    } else {
        bs.distinct_knots_and_mults_v()
    };
    if knots.len() <= 2 {
        return;
    }
    let degree = if u_direction { bs.deg_u } else { bs.deg_v };
    let periodic = if u_direction { bs.u_periodic } else { bs.v_periodic };
    let (first_k, last_k) = if periodic {
        (1, mults.len() as i32)
    } else {
        (
            occt_core::bspl::locate::first_u_knot_index(degree as i32, &mults),
            occt_core::bspl::locate::last_u_knot_index(degree as i32, &mults),
        )
    };
    let first_ind = first_k + 1;
    let mut last_ind = last_k - 1;
    let mut window_first = splits.first().copied().unwrap_or(0.0);
    let mut j = 2i32;
    while j <= splits.len() as i32 {
        let window_last = splits[(j as usize) - 1];
        let mut iknot = first_ind;
        while iknot <= last_ind {
            let (knots, mults) = if u_direction {
                bs.distinct_knots_and_mults_u()
            } else {
                bs.distinct_knots_and_mults_v()
            };
            if iknot < 1 || iknot > knots.len() as i32 {
                break;
            }
            let valknot = knots[(iknot as usize) - 1];
            if valknot <= window_first + precision {
                iknot += 1;
                continue;
            }
            if valknot >= window_last - precision {
                break;
            }
            let mult = mults[(iknot as usize) - 1];
            let mut continuity = degree as i32 - mult;
            if continuity < cont {
                let new_mult = degree as i32 - cont;
                let mut corrected = false;
                if new_mult >= 0 {
                    let removed = if u_direction {
                        bs.remove_u_knot(iknot, new_mult, tolerance)
                    } else {
                        bs.remove_v_knot(iknot, new_mult, tolerance)
                    };
                    corrected = removed.unwrap_or(false);
                }
                if corrected && new_mult > 0 {
                    let (_, mults) = if u_direction {
                        bs.distinct_knots_and_mults_u()
                    } else {
                        bs.distinct_knots_and_mults_v()
                    };
                    if (iknot as usize) <= mults.len() {
                        continuity = degree as i32 - mults[(iknot as usize) - 1];
                        corrected = continuity >= cont;
                    }
                }
                if corrected {
                    if new_mult == 0 {
                        iknot -= 1;
                        last_ind -= 1;
                    }
                } else {
                    let (knots, _) = if u_direction {
                        bs.distinct_knots_and_mults_u()
                    } else {
                        bs.distinct_knots_and_mults_v()
                    };
                    if (iknot as usize) <= knots.len() {
                        let u = knots[(iknot as usize) - 1];
                        splits.insert((j as usize) - 1, u);
                        j += 1;
                    }
                }
            }
            iknot += 1;
        }
        window_first = window_last;
        j += 1;
    }
}

/// `ShapeUpgrade_SplitCurve3dContinuity::Compute` B-spline arm
/// (`ShapeUpgrade_SplitCurve3dContinuity.cxx:146-235`).
fn raise_curve_knots(curve: &mut GeomBSplineCurve, splits: &mut Vec<f64>, cont: i32, tolerance: f64) {
    let precision = PCONFUSION;
    let (knots, mults) = curve.distinct_knots_and_mults();
    if knots.len() <= 2 {
        return;
    }
    let degree = curve.degree() as i32;
    let (first_k, last_k) = if curve.is_periodic() {
        (1, mults.len() as i32)
    } else {
        (
            occt_core::bspl::locate::first_u_knot_index(degree, &mults),
            occt_core::bspl::locate::last_u_knot_index(degree, &mults),
        )
    };
    let mut last_ind = last_k - 1;
    let first_ind = first_k + 1;
    let mut window_first = splits.first().copied().unwrap_or(0.0);
    let mut j = 2i32;
    while j <= splits.len() as i32 {
        let window_last = splits[(j as usize) - 1];
        let mut iknot = first_ind;
        while iknot <= last_ind {
            let (knots, mults) = curve.distinct_knots_and_mults();
            if iknot < 1 || iknot > knots.len() as i32 {
                break;
            }
            let valknot = knots[(iknot as usize) - 1];
            if valknot <= window_first + precision {
                iknot += 1;
                continue;
            }
            if valknot > window_last - precision {
                break;
            }
            let mult = mults[(iknot as usize) - 1];
            let mut continuity = degree - mult;
            if continuity < cont {
                let mut new_mult = degree - cont;
                if new_mult < 0 {
                    new_mult = 0;
                }
                let corrected = match curve.remove_knot(iknot, new_mult, tolerance) {
                    Ok(flag) => flag,
                    Err(_) => false,
                };
                let mut corrected = corrected;
                if corrected && new_mult > 0 {
                    let (_, mults) = curve.distinct_knots_and_mults();
                    if (iknot as usize) <= mults.len() {
                        continuity = degree - mults[(iknot as usize) - 1];
                        corrected = continuity >= cont;
                    }
                }
                if corrected {
                    if new_mult == 0 {
                        let (_, mults) = curve.distinct_knots_and_mults();
                        let last_k = if curve.is_periodic() {
                            mults.len() as i32
                        } else {
                            occt_core::bspl::locate::last_u_knot_index(degree, &mults)
                        };
                        last_ind = last_k - 1;
                        iknot -= 1;
                    }
                } else {
                    let (knots, _) = curve.distinct_knots_and_mults();
                    if (iknot as usize) <= knots.len() {
                        let u = knots[(iknot as usize) - 1];
                        splits.insert((j as usize) - 1, u);
                        j += 1;
                    }
                }
            }
            iknot += 1;
        }
        window_first = window_last;
        j += 1;
    }
}

fn clone_bspline(curve: &dyn Curve) -> Option<GeomBSplineCurve> {
    let poles = curve.bspline_poles()?.to_vec();
    let knots = curve.bspline_knots()?.to_vec();
    let degree = curve.nurbs_degree()?;
    let mut built = if let Some(weights) = curve.bspline_weights() {
        GeomBSplineCurve::rational(poles, weights.to_vec(), knots, degree).ok()?
    } else {
        GeomBSplineCurve::new(poles, knots, degree).ok()?
    };
    built.periodic = curve.is_periodic();
    Some(built)
}
