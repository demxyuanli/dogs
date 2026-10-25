//! `TopOpeBRepTool_PROJECT` — the surface-projection helper the orientation
//! comparisons read.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_PROJECT.hxx:98`,
//! `TopOpeBRepTool_PROJECT.cxx:225-256`.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::surface::Surface;

use crate::brep_tool::BRepTool;
use crate::shape::Face;

/// `FUN_tool_projPonS(P, S, UV, dist, anExtFlag, anExtAlgo)`
/// (`PROJECT.cxx:225-256`): `GeomAPI_ProjectPointOnSurf::Init(S, S->Bounds(),
/// anExtAlgo)` then `Perform(P)`, reporting `UV` and `LowerDistance()`.
///
/// `Extrema_ExtPS::SetFlag` (`cxx:238`) selects which extremum the search
/// keeps; the port's `project_point_on_surface` (`Extrema_ExtPS`) always
/// returns the lowest-distance solution, which is `Extrema_ExtFlag_MIN` — the
/// `Extrema_ExtFlag` default (`Extrema_ExtFlag.hxx:24`) and the flag the only
/// caller in `ShapeTool.cxx:488` uses.
///
/// `false` covers both OCCT failure branches: `!IsDone()` (`cxx:242-245`) and
/// `NbPoints() == 0` (`cxx:246-249`).
pub fn fun_tool_proj_pon_s(p: &GpPnt, s: &Arc<dyn Surface>) -> Option<(f64, f64)> {
    let proj = occt_geom::geom_api::project_point_on_surface(
        s.as_ref(),
        p,
        occt_core::precision::CONFUSION,
    )?;
    Some((proj.u, proj.v))
}

/// `FUN_tool_projPonF(P, F, UV, dist, anExtFlag, anExtAlgo)`
/// (`TopOpeBRepTool_PROJECT.cxx:350-361`): `BRep_Tool::Surface(F)` then
/// [`fun_tool_proj_pon_s`]. `dist` starts at `1.` in the C++ (`cxx:358`) and is
/// only overwritten by a successful projection.
pub fn fun_tool_proj_pon_f(p: &GpPnt, f: &Face) -> Option<(f64, f64, f64)> {
    let s = BRepTool::face_surface_world(f)?;
    let proj = occt_geom::geom_api::project_point_on_surface(
        s.as_ref(),
        p,
        occt_core::precision::CONFUSION,
    )?;
    Some((proj.u, proj.v, proj.distance))
}

/// `FUN_tool_projPonboundedF(P, F, UV, dist)` (`PROJECT.cxx:304-347`):
/// `BRepExtrema_ExtPF` of the point's vertex against the face, keeping the
/// smallest `SquareDistance` solution (`cxx:325-340`) and answering
/// `dist = sqrt(dmin2)` (`cxx:344`).
///
/// Unlike [`fun_tool_proj_pon_f`], the face restriction is taken into account.
pub fn fun_tool_proj_pon_bounded_f(p: &GpPnt, f: &Face) -> Option<(f64, f64, f64)> {
    let ext_pf = crate::brep_extrema_pf::ExtPF::perform(p, f);
    // `cxx:315-323`
    if !ext_pf.is_done() {
        return None;
    }
    let npt = ext_pf.nb_ext();
    if npt == 0 {
        return None;
    }
    // `cxx:325-340`: index of the smallest square distance.
    let mut dmin2 = 1.0e14;
    let mut imin = 0usize;
    for i in 1..=npt {
        let Some(d2) = ext_pf.square_distance(i) else {
            continue;
        };
        if d2 < dmin2 {
            imin = i;
            dmin2 = d2;
        }
    }
    if imin == 0 {
        return None;
    }
    // `cxx:342-346`
    let (u, v) = ext_pf.parameter(imin)?;
    Some((u, v, dmin2.sqrt()))
}
