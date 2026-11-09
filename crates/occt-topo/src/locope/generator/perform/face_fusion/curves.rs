//! p-curve placement on the rebuilt faces (`LocOpe_Generator.cxx:638-684`,
//! `cxx:864-921`, `cxx:1021-1065`).
//!
//! When the merged face does not sit on a plane, every edge put back on it needs
//! a p-curve on the new face. The reference projects the 3d curve onto the
//! stored surface, then walks the projected curve by whole `2*pi` steps until it
//! lands inside the face's U range (`cxx:657-676`) - the "recalage dans la
//! facette" - because a projected curve may come back shifted by a period.

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_core::precision::{epsilon, ANGULAR};
use occt_geom::geom_proj_lib;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_offset::update_edge_pcurve;
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::shape::{Edge, Face, TopoShape};

/// `Geom_TrimmedCurve` peel of the p-curve call sites (`cxx:650-653`).
pub(super) fn untrimmed(c: Arc<dyn Curve>) -> Arc<dyn Curve> {
    if c.is_geom_trimmed() {
        c.trim_basis().unwrap_or(c)
    } else {
        c
    }
}

/// `gp_Geom2d_Curve::Translated(gp_Vec2d)` (`cxx:662`, `cxx:670`).
fn translated(curve: &mut Box<dyn Curve2d>, dx: f64) {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, 0.0));
    curve.transform(&t);
}

/// `cxx:638-684`: project `edg` onto `s` and store the result as the p-curve of
/// `new_face`.
///
/// `closed_face` is the `facbis` the reference passes to
/// `BRepTools::IsReallyClosed`; a really closed edge gets no p-curve here.
pub(super) fn put_curve_2d(
    new_face: &TopoShape,
    closed_face: &TopoShape,
    edg: &TopoShape,
    s: &Arc<dyn Surface>,
    u_minc: f64,
    u_maxc: f64,
) {
    let (Some(face), Some(target), Some(nf)) = (
        Face::wrap(closed_face.clone()),
        Edge::wrap(edg.clone()),
        Face::wrap(new_face.clone()),
    ) else {
        return;
    };
    // `cxx:642-643`.
    let tol = BRepTool::edge_tolerance(&target);
    let Some(curve) = BRepTool::edge_curve_world(&target) else {
        return;
    };
    let curve = untrimmed(curve);
    let (f, l) = BRepTool::edge_parameters(&target);
    // `cxx:656`: `GeomProjLib::Curve2d(C, f, l, S, tol)`. The reference keeps
    // `tol` as the edge tolerance when it stores the p-curve afterwards, so the
    // projection's own `Precision::Confusion()` clamp stays local.
    let mut proj_tol = tol;
    let Some(c2d) = geom_proj_lib::curve_2d_range(&curve, f, l, s, &mut proj_tol) else {
        return;
    };
    // `cxx:657-676`: shift by whole periods until the curve sits inside the
    // face's U range.
    let mut c2d = c2d.clone_dyn();
    let tttol = ANGULAR;
    loop {
        let (pf, pl) = (c2d.value(f), c2d.value(l));
        if pf.x().min(pl.x()) >= u_maxc - tttol {
            translated(&mut c2d, -2.0 * std::f64::consts::PI);
        } else {
            break;
        }
    }
    loop {
        let (pf, pl) = (c2d.value(f), c2d.value(l));
        if pf.x().max(pl.x()) <= u_minc + tttol {
            translated(&mut c2d, 2.0 * std::f64::consts::PI);
        } else {
            break;
        }
    }
    // `cxx:677-683`.
    if !is_really_closed(&target, &face) {
        update_edge_pcurve(&target, &Arc::from(c2d), &nf, tol);
    }
}

/// `cxx:864-921`: the p-curve of a twin edge appended after the wire was
/// rebuilt, realigned on the U period of a periodic surface.
#[allow(clippy::too_many_arguments)]
pub(super) fn put_twin_curve_2d(
    new_face: &TopoShape,
    fac: &TopoShape,
    edg: &TopoShape,
    twin: &TopoShape,
    s: &Arc<dyn Surface>,
    tol: f64,
    dont_fuse: Option<&TopoShape>,
    edg_shortened: bool,
) {
    let Some(twin_edge) = Edge::wrap(twin.clone()) else {
        return;
    };
    let Some(twin_shape) = Face::wrap(new_face.clone()) else {
        return;
    };
    let Some(curve) = BRepTool::edge_curve_world(&twin_edge) else {
        return;
    };
    let curve = untrimmed(curve);
    let (f, l) = BRepTool::edge_parameters(&twin_edge);
    let mut proj_tol = tol;
    let Some(ptc) = geom_proj_lib::curve_2d_range(&curve, f, l, s, &mut proj_tol) else {
        return;
    };
    let mut ptc = ptc.clone_dyn();
    // `cxx:869-918`.
    if s.is_u_periodic() {
        let u_ref = match (dont_fuse, Face::wrap(fac.clone()), Edge::wrap(edg.clone())) {
            (Some(dont), Some(fac_face), Some(edge)) => {
                let oredge = dont.orientation();
                // `cxx:874-902`.
                let probe = if edg_shortened {
                    edge.0.oriented(oredge)
                } else {
                    edge.0.oriented(oredge.reversed())
                };
                let Some(probe) = Edge::wrap(probe) else {
                    return;
                };
                BRepTool::uv_points(&probe, &fac_face).0.x()
            }
            (None, Some(fac_face), Some(edge)) => {
                // `cxx:903-906`.
                BRepTool::uv_points(&edge, &fac_face).0.x()
            }
            _ => return,
        };
        let new_u = ptc.value(f).x();
        if (new_u - u_ref).abs() > epsilon(s.u_period()) {
            translated(&mut ptc, u_ref - new_u);
        }
    }
    update_edge_pcurve(&twin_edge, &Arc::from(ptc), &twin_shape, tol);
}
