//! `TopOpeBRepTool_TOOL` approximate-normal family.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `FUN_tool_nC2dINSIDES(tgC2d)` (`TopOpeBRepTool_matter.cxx:45-62`)
//! - `TopOpeBRepTool_TOOL::outUVbounds(uv, F)` (`TopOpeBRepTool_TOOL.cxx:1915-1921`)
//! - `TopOpeBRepTool_TOOL::tryOriEinF(par, e, f)` (`TOOL.cxx:779-814`)
//! - `TopOpeBRepTool_TOOL::NgApp(par, e, f, tola, ngApp)` (`TOOL.cxx:818-918`)
//! - `TopOpeBRepTool_TOOL::tryNgApp(par, e, f, tola, Ng)` (`TOOL.cxx:922-944`)
//!
//! `Matter` needs a second, *approximate* sample of the normal when the two
//! faces are tangent on the edge (the matter angle collapses to `0` or `2*pi`);
//! `tryNgApp` is the entry point the matter helpers call for that.

use occt_core::gp::{GpDir, GpDir2d, GpPnt2d, GpVec2d};

use crate::shape::{Face, TopoShape};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::closed::ori_in_sor;
use crate::topopebrep::tool::curvature::min_duv;
use crate::topopebrep::tool::geometry::fun_tool_ng_s;
use crate::topopebrep::tool::topology::fun_tool_ng_geom_f;
use crate::topopebrep::tool::topology::fun_tool_par_on_e_f;

/// `FUN_tool_nC2dINSIDES(tgC2d)` (`matter.cxx:45-62`): `Y = Z ^ X` for
/// `X = (tgC2d, 0)`, `Z = (0, 0, 1)`, i.e. the in-plane left normal.
pub fn n_c2d_inside_s(tg: &GpDir2d) -> Option<GpDir2d> {
    // `gp_Dir(tgC2d.X(), tgC2d.Y(), 0.)` normalises; a null tangent throws in
    // OCCT, the port answers `None`.
    let l = (tg.x() * tg.x() + tg.y() * tg.y()).sqrt();
    if l <= occt_core::precision::REAL_SMALL {
        return None;
    }
    // `Z ^ X = (-X.y, X.x, 0)`
    GpDir2d::new(-tg.y() / l, tg.x() / l).ok()
}

/// `TopOpeBRepTool_TOOL::outUVbounds(uv, F)` (`TOOL.cxx:1915-1921`): `uv` lies
/// outside the surface bounds of `F`.
pub fn out_uv_bounds(uv: &GpPnt2d, f: &TopoShape) -> bool {
    let bs = SurfaceAdaptor::of_face(&Face(f.clone()));
    let out_u = uv.x() > bs.last_u_parameter() || uv.x() < bs.first_u_parameter();
    let out_v = uv.y() > bs.last_v_parameter() || uv.y() < bs.first_v_parameter();
    out_u || out_v
}

/// `TopOpeBRepTool_TOOL::tryOriEinF(par, e, f)` (`TOOL.cxx:779-814`): `1` when
/// `e` is `FORWARD` in `f`, `2` when `REVERSED`, `0` when undecided.
///
/// `checkclo = true` (`cxx:786`); the stored p-curve must already exist
/// (`cxx:795-799`), and the normal side is sampled at `minDUV(f) * 1e-1`
/// (`cxx:808-812`).
pub fn try_ori_e_in_f(par: f64, e: &TopoShape, f: &TopoShape) -> i32 {
    let oeinf = ori_in_sor(e, f, true);
    if oeinf != 0 {
        return oeinf;
    }

    let Some(c2df) = crate::topopebrep::tool::c2d::fc2d_has_old_curve_on_surface(e, f) else {
        return 0;
    };
    let (curve, _first, _last, _tol) = c2df.pc();
    let Some(curve) = curve else {
        return 0;
    };

    // `pc->D1(par, uv, tg2d)`; `gp_Dir2d(-tg2d.Y(), tg2d.X())` is the
    // left normal, oriented inside `f`2d for a `FORWARD` edge.
    let (uv, tg2d) = curve.d1(par);
    let Some(n2d) = GpDir2d::new(-tg2d.y(), tg2d.x()).ok() else {
        return 0;
    };

    let mut delta = min_duv(&Face(f.clone()));
    delta *= 1.0e-1;
    let ouv = uv.translated_vec(&GpVec2d::from_dir2d(&n2d).multiplied_scalar(delta));
    if out_uv_bounds(&ouv, f) {
        2
    } else {
        1
    }
}

/// `TopOpeBRepTool_TOOL::NgApp(par, e, f, tola, ngApp)` (`TOOL.cxx:818-918`):
/// the normal of `f` shifted off the edge, so that it differs from the normal
/// *on* the edge (plane faces and non-quadrics answer `None`, `cxx:840-852`).
pub fn ng_app(par: f64, e: &TopoShape, f: &TopoShape, tola: f64) -> Option<GpDir> {
    let surface = crate::brep_tool::BRepTool::face_surface_world(&Face(f.clone()))?;
    if crate::topopebrep::tool::geometry::fun_tool_plane(&Face(f.clone())) {
        return None;
    }
    if !crate::topopebrep::tool::geometry::fun_tool_quad_surface(surface.as_ref()) {
        return None;
    }
    // `cxx:854-861`: the p-curve must already be stored.
    crate::topopebrep::tool::c2d::fc2d_has_old_curve_on_surface(e, f)?;

    // `cxx:863-870`: the orientation is taken on `f.Oriented(FORWARD)`.
    let fwd = f.clone().oriented(crate::abs::Orientation::Forward);
    let orieinf = try_ori_e_in_f(par, e, &fwd);
    if orieinf == 0 {
        return None;
    }

    let uv0 = fun_tool_par_on_e_f(e, par, f)?;
    let mut ng = fun_tool_ng_s(&uv0, surface.as_ref());

    let c2df = crate::topopebrep::tool::c2d::fc2d_has_old_curve_on_surface(e, f)?;
    let (curve, _first, _last, _tol) = c2df.pc();
    let curve = curve?;
    let (uv, tg2d) = curve.d1(par);

    let mut n2d_inside_s = n_c2d_inside_s(&GpDir2d::new(tg2d.x(), tg2d.y()).ok()?)?;
    if orieinf == 2 {
        n2d_inside_s.reverse();
    }

    // `cxx:894-916`: five samples of growing length, the first one whose
    // normal both differs from `ng` and is not degenerated wins.
    let nmax = 5;
    let mut same = false;
    let mut delta = 0.45678;
    let mut duv = GpVec2d::from_dir2d(&n2d_inside_s).multiplied_scalar(0.45678);
    for _i in 1..=nmax {
        let newuv = uv.translated_vec(&duv);
        let newng = fun_tool_ng_s(&newuv, surface.as_ref());
        // `ng.IsEqual(newng, tola)` is `|Angle| <= tola` (`gp_Dir::IsEqual`);
        // `gp_Vec newng = FUN_tool_ngS(...)` normalises, so `newng.Magnitude()`
        // is `1` and `cxx:907` reduces to `1 > tola`.
        same = ng.angle(&newng).abs() <= tola;
        let okk = 1.0 > tola;
        if !same && okk {
            ng = newng;
            break;
        }
        delta *= 1.25; // `cxx:913`: NYI
        duv = GpVec2d::from_dir2d(&n2d_inside_s).multiplied_scalar(delta);
    }
    if same {
        None
    } else {
        Some(ng)
    }
}

/// `TopOpeBRepTool_TOOL::tryNgApp(par, e, f, tola, Ng)` (`TOOL.cxx:922-944`):
/// [`ng_app`], falling back to the plain geometric normal.
pub fn try_ng_app(par: f64, e: &TopoShape, f: &TopoShape, tola: f64) -> Option<GpDir> {
    let uv = fun_tool_par_on_e_f(e, par, f)?;
    let ng = fun_tool_ng_geom_f(&uv, &Face(f.clone()))?;
    Some(ng_app(par, e, f, tola).unwrap_or(ng))
}
