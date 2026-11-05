//! The spine-parameter statics of `ChFi3d_Builder_SpKP.cxx`:
//! `ComputeAbscissa` (`:202-213`), `ParamOnSpine` (`:216-333`),
//! `YaUnVoisin` (`:336-368`) and `ResetProl` (`:510-545`).
//!
//! Sources: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0).

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln};
use occt_geom::plane::GeomPlane;
use occt_geom::surface::Surface;

use crate::chfi3d::builder_0_helpers::inter_plane_edge::inter_plane_edge;
use crate::chfi_ds::{SpineHandle, SurfData};
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::Edge;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

/// `Spine->CurrentElementarySpine(iedge)` (`ChFiDS_Spine.hxx`).
pub(super) fn elementary_spine(
    spine: &SpineHandle,
    iedge: i32,
) -> Result<CurveAdaptor, String> {
    let mut guard = spine
        .write()
        .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
    guard
        .current_elementary_spine(iedge)
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: no elementary spine".to_string())
}

/// `Spine->Edges(iedge)`.
pub(super) fn spine_edge(spine: &SpineHandle, iedge: i32) -> Result<Edge, String> {
    let guard = spine
        .read()
        .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
    let e = guard
        .edges(iedge as usize)
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: no edge on the spine".to_string())?;
    Edge::wrap(e.clone())
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: not an edge".to_string())
}

/// `ComputeAbscissa(C, U)` (`cxx:202-213`).
pub(super) fn compute_abscissa(c: &CurveAdaptor, u: f64) -> f64 {
    match c.get_type() {
        GeomAbsCurveType::Line => u,
        GeomAbsCurveType::Circle => c
            .curve()
            .and_then(|base| base.gp_circ())
            .map(|circle| circle.radius() * u)
            .unwrap_or(0.),
        _ => 0.,
    }
}

/// `ParamOnSpine(DStr, ptg, CD, Spine, iedge, intf, intl, tol, pok)`
/// (`cxx:216-333`).
///
/// Answers `(Nl, pok)`.
#[allow(clippy::too_many_arguments)]
pub(super) fn param_on_spine(
    ds: &DataStructure,
    ptg: f64,
    cd: &SurfData,
    spine: &SpineHandle,
    iedge: i32,
    intf: bool,
    intl: bool,
    tol: f64,
) -> Result<(f64, bool), String> {
    let (f, l) = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
        (
            guard
                .first_parameter_of(iedge)
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: FirstParameter".to_string())?,
            guard
                .last_parameter_of(iedge)
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: LastParameter".to_string())?,
        )
    };

    let elementary = elementary_spine(spine, iedge)?;
    let mut nl = compute_abscissa(&elementary, ptg) + f;
    if (nl >= f - tol || intf) && (nl <= l + tol || intl) {
        return Ok((nl, true));
    }

    // Construction of the plane containing the section of CD with parameter
    // ptg.
    let mut c3d = None;
    if cd.interference_on_s1().line_index() != 0 {
        c3d = ds
            .curve(cd.interference_on_s1().line_index())
            .curve()
            .cloned();
    }
    if c3d.is_none() {
        c3d = ds
            .curve(cd.interference_on_s2().line_index())
            .curve()
            .cloned();
    }
    let Some(c3d) = c3d else {
        return Ok((nl, false));
    };
    let (pp, vv) = c3d.d1(ptg);
    let dir = GpDir::new(vv.x(), vv.y(), vv.z())
        .map_err(|_| "ChFi3d_Builder::SplitKPart: null section direction".to_string())?;
    let ax3 = GpAx3::from_point_dir(pp, dir)
        .map_err(|_| "ChFi3d_Builder::SplitKPart: null section plane".to_string())?;
    let plan: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(ax3)));

    let (periodic, nb_edges) = {
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
        (guard.is_periodic(), guard.nb_edges() as i32)
    };

    // Intersection plane spine.
    let mut sens = 1;
    if nl <= f {
        sens = -1;
    }
    let mut ii = iedge + sens;
    if periodic {
        if ii <= 0 {
            ii += nb_edges;
        }
        if ii > nb_edges {
            ii -= nb_edges;
        }
    } else if ii < 1 || ii > nb_edges {
        return Ok((nl, true));
    }

    let mut found = false;
    let mut fini = false;
    while !found && !fini {
        let o = spine_edge(spine, ii)?.orientation();
        let first = (o.is_forward() && sens == 1) || (o.is_reversed() && sens == -1);
        let ce = elementary_spine(spine, ii)?;
        let tolc = ce.resolution(tol);
        let curve = ce
            .curve()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::SplitKPart: empty elementary spine".to_string())?;
        found = inter_plane_edge(plan.as_ref(), curve.as_ref(), &mut nl, first, tolc);
        if found {
            nl = spine
                .write()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
                .absc_on(nl, ii)
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: Absc".to_string())?;
        }

        ii += sens;
        if periodic {
            if ii <= 0 {
                ii += nb_edges;
            }
            if ii > nb_edges {
                ii -= nb_edges;
            }
            fini = ii == iedge;
        } else {
            fini = ii < 1 || ii > nb_edges;
        }
    }
    Ok((nl, found))
}

/// `YaUnVoisin(Spine, iedge, ivois, isfirst)` (`cxx:336-368`): answers the
/// neighbour index when there is one.
pub(super) fn ya_un_voisin(
    spine: &SpineHandle,
    iedge: i32,
    is_first: bool,
) -> Result<Option<i32>, String> {
    let guard = spine
        .read()
        .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
    let nbed = guard.nb_edges() as i32;
    if nbed == 1 {
        return Ok(None);
    }
    let periodic = guard.is_periodic();
    let mut ivois = if is_first { iedge - 1 } else { iedge + 1 };
    if periodic {
        if ivois == 0 {
            ivois = nbed;
        }
        if ivois == nbed + 1 {
            ivois = 1;
        }
    }
    Ok(if ivois > 0 && ivois <= nbed {
        Some(ivois)
    } else {
        None
    })
}

/// `ResetProl(DStr, CD, Spine, iedge, isfirst)` (`cxx:510-545`).
pub(super) fn reset_prol(
    ds: &DataStructure,
    cd: &mut SurfData,
    spine: &SpineHandle,
    iedge: i32,
    is_first: bool,
) -> Result<f64, String> {
    let bc = elementary_spine(spine, iedge)?;
    let edglen = bc.last_parameter() - bc.first_parameter();
    let surf: Arc<dyn Surface> = ds
        .surface(cd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null surface".to_string())?;
    let par = if is_first { 0. } else { edglen };
    let mut sppar = 0.;
    for i in 1..=2 {
        cd.interference_mut(i).set_parameter(par, is_first);
        let pc = cd
            .interference(i)
            .p_curve_on_surf()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null pcurve".to_string())?;
        let p = pc.value(par);
        cd.vertex_mut(is_first, i).reset();
        cd.vertex_mut(is_first, i)
            .set_point(&surf.value(p.x(), p.y()));
        if is_first {
            sppar = spine
                .read()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
                .first_parameter_of(iedge)
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: FirstParameter".to_string())?;
            cd.set_first_spine_param(sppar);
        } else {
            sppar = spine
                .read()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
                .last_parameter_of(iedge)
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: LastParameter".to_string())?;
            cd.set_last_spine_param(sppar);
        }
    }
    Ok(sppar)
}

/// `ChFi3d_KParticular(Spine, ivois, S1, S2)` reached from `Trunc`, kept here
/// so `trunc.rs` only names one entry point per helper.
pub(super) fn k_particular(
    spine: &SpineHandle,
    ivois: i32,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
) -> bool {
    crate::chfi3d::builder_0_helpers::k_particular::ch_fi3d_k_particular(
        spine, ivois, s1, s2,
    )
}
