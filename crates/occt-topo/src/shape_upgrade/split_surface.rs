//! `ShapeUpgrade_SplitSurface::Build` (`ShapeUpgrade_SplitSurface.cxx:240-635`)
//! for a B-spline grid and for the revolution / extrusion generatrix split
//! (`ShapeUpgrade_SplitCurve3d::Build`, `ShapeUpgrade_SplitCurve3d.cxx:246-315`).
//!
//! More than one patch is handed to `ShapeFix_ComposeShell`
//! (`ShapeUpgrade_FaceDivide.cxx:181-198`). A shell is not a face, so the
//! IGES offset caller returns it (`IGESToBRep_TopoSurface.cxx:1136-1140`).

use std::sync::Arc;

use occt_core::precision::{PCONFUSION, CONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::{Curve, GeomRectangularTrimmedSurface, Surface};

use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::shape_fix_compose_shell::{ComposeShell, CompositeSurface};

/// One corrected surface, or the shell `ComposeShell` built from a patch grid.
pub enum C0Basis {
    Surface(Arc<dyn Surface>),
    Shape(TopoShape),
}

/// `ShapeUpgrade_SplitSurface::Build` B-spline arm (`cxx:505-618`).
pub fn split_bspline(
    face: &Face,
    bs: &GeomBSplineSurface,
    u_splits: &mut [f64],
    v_splits: &mut [f64],
) -> Result<C0Basis, String> {
    let (uk, _) = bs.distinct_knots_and_mults_u();
    let (vk, _) = bs.distinct_knots_and_mults_v();
    snap_splits(&uk, u_splits);
    snap_splits(&vk, v_splits);
    let mut rows = Vec::new();
    let mut u1 = u_splits[0];
    for irow in 1..u_splits.len() {
        let u2 = u_splits[irow];
        let mut v1 = v_splits[0];
        let mut row = Vec::new();
        for icol in 1..v_splits.len() {
            let v2 = v_splits[icol];
            row.push(segment_bspline(bs, u1, u2, v1, v2));
            v1 = v2;
        }
        rows.push(row);
        u1 = u2;
    }
    finish_grid(face, rows, u_splits, v_splits)
}

/// Revolution arm (`ShapeUpgrade_SplitSurface.cxx:255-331`) when only V is split.
pub fn split_revolution(
    face: &Face,
    curve: &GeomBSplineCurve,
    axis: occt_core::gp::GpAx1,
    u_first: f64,
    u_last: f64,
    v_splits: &mut [f64],
) -> Result<C0Basis, String> {
    let (knots, _) = curve.distinct_knots_and_mults();
    snap_splits(&knots, v_splits);
    let pieces = split_bspline_curve(curve, v_splits);
    let mut row = Vec::with_capacity(pieces.len());
    for piece in pieces {
        row.push(Arc::new(occt_geom::GeomSurfaceOfRevolution::new(piece, axis)) as Arc<dyn Surface>);
    }
    finish_grid(face, vec![row], &[u_first, u_last], v_splits)
}

/// Extrusion arm (`ShapeUpgrade_SplitSurface.cxx:334-410`) when only U is split.
pub fn split_extrusion(
    face: &Face,
    curve: &GeomBSplineCurve,
    direction: occt_core::gp::GpDir,
    v_first: f64,
    v_last: f64,
    u_splits: &mut [f64],
) -> Result<C0Basis, String> {
    let (knots, _) = curve.distinct_knots_and_mults();
    snap_splits(&knots, u_splits);
    let pieces = split_bspline_curve(curve, u_splits);
    let mut rows = Vec::with_capacity(pieces.len());
    for piece in pieces {
        rows.push(vec![
            Arc::new(occt_geom::GeomSurfaceOfLinearExtrusion::new(piece, direction))
                as Arc<dyn Surface>,
        ]);
    }
    finish_grid(face, rows, u_splits, &[v_first, v_last])
}

fn segment_bspline(bs: &GeomBSplineSurface, u1: f64, u2: f64, v1: f64, v2: f64) -> Arc<dyn Surface> {
    let mut copy = bs.clone();
    if copy.segment_uv(u1, u2, v1, v2, PCONFUSION, PCONFUSION).is_ok() {
        return Arc::new(copy);
    }
    // `cxx:593-602`: a failed `Segment` becomes a rectangular trim.
    Arc::new(GeomRectangularTrimmedSurface::uv(
        Arc::new(bs.clone()),
        u1,
        u2,
        v1,
        v2,
    ))
}

/// `ShapeUpgrade_SplitCurve3d::Build` B-spline span loop (`cxx:273-315`).
fn split_bspline_curve(curve: &GeomBSplineCurve, splits: &[f64]) -> Vec<Arc<dyn Curve>> {
    let mut out = Vec::new();
    for w in splits.windows(2) {
        let mut copy = curve.clone();
        if copy.segment(w[0], w[1], PCONFUSION).is_ok() {
            out.push(Arc::new(copy) as Arc<dyn Curve>);
        } else {
            out.push(Arc::new(occt_geom::trimmed::GeomTrimmedCurve::new(
                Arc::new(curve.clone()),
                w[0],
                w[1],
            )) as Arc<dyn Curve>);
        }
    }
    out
}

fn finish_grid(
    face: &Face,
    patches: Vec<Vec<Arc<dyn Surface>>>,
    u_joints: &[f64],
    v_joints: &[f64],
) -> Result<C0Basis, String> {
    let one = patches.len() == 1 && patches.first().map_or(false, |r| r.len() == 1);
    if one {
        return Ok(C0Basis::Surface(patches.into_iter().next().unwrap().into_iter().next().unwrap()));
    }
    // `ShapeUpgrade_FaceDivide::SplitSurface` (`cxx:181-198`).
    let mut grid = CompositeSurface::new();
    grid.init_joints(patches, u_joints, v_joints);
    let mut shell = ComposeShell::new();
    shell.init(grid, face, CONFUSION);
    let _ = shell.perform();
    let Some(shape) = shell.result().cloned() else {
        return Err("ShapeUpgrade_FaceDivide::SplitSurface: ShapeUpgrade_FaceDivide.cxx:193".into());
    };
    if shape.is_face() {
        let face = Face(shape);
        let surf = BRepTool::face_surface_world(&face).ok_or_else(|| {
            "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1160".to_string()
        })?;
        return Ok(C0Basis::Surface(surf));
    }
    Ok(C0Basis::Shape(shape))
}

/// Snap a split value onto a knot within `PConfusion`
/// (`ShapeUpgrade_SplitSurface.cxx:511-529`).
fn snap_splits(knots: &[f64], splits: &mut [f64]) {
    if knots.is_empty() {
        return;
    }
    let last = knots.len() - 1;
    let mut j = 0usize;
    for sp in splits.iter_mut() {
        while j <= last {
            if *sp > knots[j] + PCONFUSION {
                j += 1;
                continue;
            }
            if *sp < knots[j] - PCONFUSION {
                break;
            }
            *sp = knots[j];
            j += 1;
            break;
        }
        if j == last {
            break;
        }
    }
}
