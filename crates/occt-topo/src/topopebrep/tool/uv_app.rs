//! `TopOpeBRepTool_TOOL` uv-application family.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `FUN_tool_isobounds(F, u1, u2, v1, v2)` (`TopOpeBRepTool_TOPOLOGY.cxx:155-220`)
//! - `FUN_nearestISO(F, xpar, isoU, xinf, xsup)` (`TOPOLOGY.cxx:585-643`)
//! - `FUN_tool_getdxx(F, E, dxx)` (`TopOpeBRepTool_matter.cxx:150-230`)
//! - `TopOpeBRepTool_TOOL::uvApp(f, e, pare, eps, uvapp)`
//!   (`TopOpeBRepTool_TOOL.cxx:1491-1512`)
//!
//! `uvApp` moves the uv point of `pare` on `f` by `eps` along the direction that
//! stays inside `f`2d, which is the coarse sample `MatterKPtg`
//! (`TOOL.cxx:1710-1771`) needs.

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::abs::Orientation;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::tool::c2d;
use crate::topopebrep::tool::geometry;
use crate::topopebrep::tool::iso;
use crate::topopebrep::tool::topology;

/// `FUN_tool_isobounds(F, u1, u2, v1, v2)` (`cxx:155-220`): the uv box spanned
/// by the iso lines that bound `F`.
///
/// Every edge must have a p-curve on `F` (`cxx:189-193` -> `false`) and that
/// p-curve must be isoparametric (`cxx:213-216` -> `false`). `uclosed`/`vclosed`
/// come from `FUN_tool_closedS` (`cxx:176-178`) but no branch reads them: the
/// `if (!uvclosed)` arm is commented out (`cxx:180-184`).
pub fn fun_tool_isobounds(f: &TopoShape) -> Option<(f64, f64, f64, f64)> {
    let face = Face(f.clone());
    let mut u1 = 1.0e7f64;
    let mut v1 = 1.0e7f64;
    let mut u2 = -1.0e7f64;
    let mut v2 = -1.0e7f64;

    if crate::brep_tool::BRepTool::face_surface_world(&face).is_none() {
        return None;
    }
    // `cxx:176-178`: the call is unconditional; its result is unused.
    let _ = geometry::fun_tool_closed_s(&face);

    for e in explorer_shapes(f, crate::abs::ShapeType::Edge) {
        let c2df = c2d::fc2d_curve_on_surface(&e, f, false)?;
        let (curve, first, last, _tol) = c2df.pc();
        let curve = curve?;
        // `cxx:198-216`: not an iso line -> `false`.
        iso::uv_iso_curve(curve.as_ref())?;
        let p_first = curve.d0(first);
        let p_last = curve.d0(last);
        u1 = u1.min(p_first.x());
        u2 = u2.max(p_first.x());
        v1 = v1.min(p_first.y());
        v2 = v2.max(p_first.y());
        u1 = u1.min(p_last.x());
        u2 = u2.max(p_last.x());
        v1 = v1.min(p_last.y());
        v2 = v2.max(p_last.y());
    }
    Some((u1, u2, v1, v2))
}

/// `FUN_nearestISO(F, xpar, isoU, xinf, xsup)` (`cxx:585-643`): narrow
/// `[xinf, xsup]` onto the iso lines of `F` that bracket `xpar`.
///
/// `xpar` on either bound answers `true` immediately (`cxx:597-603`). Every
/// edge of `F` must be an iso line of the requested axis (`cxx:614-622`), else
/// `false`.
pub fn fun_nearest_iso(f: &TopoShape, xpar: f64, iso_u: bool, xinf: &mut f64, xsup: &mut f64) -> bool {
    let tol2d = 1.0e-6;
    let df = xpar - *xinf;
    let onf = df.abs() < tol2d;
    let dl = xpar - *xsup;
    let onl = dl.abs() < tol2d;
    if onf || onl {
        return true;
    }

    let iso_v = !iso_u;
    for e in explorer_shapes(f, crate::abs::ShapeType::Edge) {
        let Some(iso) = iso::uv_iso_edge_face(&e, f) else {
            return false;
        };
        let compare = (iso_u && iso.iso_u) || (iso_v && iso.iso_v);
        if !compare {
            return false;
        }
        let xvalue = if iso.iso_u { iso.origin.x() } else { iso.origin.y() };
        if xvalue < xpar && *xinf < xvalue {
            *xinf = xvalue;
        }
        if xpar < xvalue && xvalue < *xsup {
            *xsup = xvalue;
        }
    }
    true
}

/// `FUN_tool_getdxx(F, E, dxx)` (`matter.cxx:150-231`): the uv vector along
/// which `F`2d has matter next to `E`.
///
/// OCCT names the third parameter `parE` but comments it out in the signature
/// (`matter.cxx:152-153`); the answer only depends on `E`'s iso position.
/// `None` covers the four `return false` arms: orientation unknown
/// (`cxx:161-166`), `INTERNAL`/`EXTERNAL` (`cxx:167-170`), not an iso line
/// (`cxx:172-177`), bounds/nearest-iso failure (`cxx:179-193`).
pub fn fun_tool_getdxx(f: &TopoShape, e: &TopoShape) -> Option<GpVec2d> {
    let ori = topology::fun_tool_orient_ein_f_forward(e, f)?;
    if ori == Orientation::Internal || ori == Orientation::External {
        return None;
    }

    let iso = iso::uv_iso_edge_face(e, f)?;
    let (u1, u2, v1, v2) = fun_tool_isobounds(f)?;

    let xpar = if iso.iso_u { iso.origin.x() } else { iso.origin.y() };
    let mut xinf = if iso.iso_u { u1 } else { v1 };
    let mut xsup = if iso.iso_u { u2 } else { v2 };
    if !fun_nearest_iso(f, xpar, iso.iso_u, &mut xinf, &mut xsup) {
        return None;
    }

    let ypar = if iso.iso_u { iso.dir.y } else { iso.dir.x };
    // `cxx:197-204`: the two `if`s (not `else if`) of OCCT.
    let mut matter_after_xpar = false;
    if iso.iso_u {
        matter_after_xpar = ypar < 0.0;
    }
    if iso.iso_v {
        matter_after_xpar = ypar > 0.0;
    }
    if ori == Orientation::Reversed {
        matter_after_xpar = !matter_after_xpar;
    }

    let dx = if matter_after_xpar {
        xsup - xpar
    } else {
        xinf - xpar
    };
    // `cxx:221-228`: both `if`s write `dxx`; exactly one axis is an iso.
    let mut dxx = GpVec2d::new(0.0, 0.0);
    if iso.iso_u {
        dxx = GpVec2d::new(dx, 0.0);
    }
    if iso.iso_v {
        dxx = GpVec2d::new(0.0, dx);
    }
    Some(dxx)
}

/// `TopOpeBRepTool_TOOL::uvApp(f, e, pare, eps, uvapp)` (`cxx:1491-1512`):
/// `UV(pare on e / f)` translated by `eps * getdxx`.
pub fn uv_app(f: &TopoShape, e: &TopoShape, pare: f64, eps: f64) -> Option<GpPnt2d> {
    let uv = topology::fun_tool_par_on_e_f(e, pare, f)?;
    let dxx = fun_tool_getdxx(f, e)?;
    Some(uv.translated_vec(&dxx.multiplied_scalar(eps)))
}
