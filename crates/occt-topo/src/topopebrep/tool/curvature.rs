//! Curvature, frame (`XX` / `Nt`) and parameter-space helpers of
//! `TopOpeBRepTool_TOOL` — the slice `TopOpeBRepTool_makeTransition` reads.
//!
//! Sources (OCCT 8.0.0, `TKBool/TopOpeBRepTool`):
//! - `TopOpeBRepTool_TOOL::IsQuad(E)` / `IsQuad(F)` (`TOOL.cxx:948-959`)
//! - `TopOpeBRepTool_TOOL::CurvE` (`TOOL.cxx:964-1016`)
//! - `FUN_analyticcS` (`TOOL.cxx:1037-1125`) and `CurvF`
//!   (`TOOL.cxx:1131-1187`)
//! - `TopOpeBRepTool_TOOL::Getduv` (`TOOL.cxx:1430-1486`)
//! - `TopOpeBRepTool_TOOL::XX` (`TOOL.cxx:1515-1544`)
//! - `TopOpeBRepTool_TOOL::Nt` (`TOOL.cxx:1548-1561`)
//! - `TopOpeBRepTool_TOOL::Getstp3dF` (`TOOL.cxx:1775-1814`)
//! - `TopOpeBRepTool_TOOL::minDUV` (`TOOL.cxx:1852-1859`)
//! - `TopOpeBRepTool_TOOL::stuvF` (`TOOL.cxx:1867-1911`)
//! - `TopOpeBRepTool_TOOL::EdgeONFace` (`TOOL.cxx:1978-2112`)
//! - `FUN_tool_getgeomxx` (`TopOpeBRepTool_matter.cxx:66-84`)
//! - `FUN_nullprodv` (`TOOL.cxx:85-90`)
//!
//! `FUN_tool_nggeomF` — which `Nt` and `XX` bottom out on — is
//! [`super::topology::fun_tool_ng_geom_f`].

use occt_core::gp::{GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec};

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape};

use super::adaptor::{CurveAdaptor, SurfaceAdaptor};
use super::closed::basis_surface;
use super::geometry::fun_quad_ct;
use super::lprop::CLProps;
use super::project::fun_tool_proj_pon_f;
use super::project_curve::fun_tool_bounds;
use super::tangent::tg_geom_e;
use super::topology::{
    fun_tool_ng_geom_f, fun_tool_orient_ein_f_forward, fun_tool_par_f_tol, fun_tool_value_uv,
};

/// `INFFIRST (-1)` / `SUPLAST (-2)` / `ONFIRST (1)` / `ONLAST (2)`
/// (`TOOL.cxx:1864-1867`), the `stuvF` answer.
pub const INFFIRST: i32 = -1;
pub const SUPLAST: i32 = -2;
pub const ONFIRST: i32 = 1;
pub const ONLAST: i32 = 2;

/// `FUN_nullprodv(prodv)` (`TOOL.cxx:85-90`): `|prodv| < 1.e-6`
/// (the `NYITOLXPU` constant is 1.e-6, not `Precision::Angular()`).
fn fun_nullprodv(prodv: f64) -> bool {
    prodv.abs() < 1.0e-6
}

/// `TopOpeBRepTool_TOOL::IsQuad(E)` (`TOOL.cxx:948-951`).
pub fn is_quad_e(e: &TopoShape) -> bool {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    fun_quad_ct(bc.get_type())
}

/// `TopOpeBRepTool_TOOL::IsQuad(F)` (`TOOL.cxx:954-959`):
/// `FUN_tool_quad(BASISSURFACE(F))`.
pub fn is_quad_f(f: &Face) -> bool {
    let Some(s) = BRepTool::face_surface_world(f) else {
        return false;
    };
    let s = basis_surface(&s);
    fun_tool_quad_surface_checked(s.as_ref())
}

/// `FUN_tool_quad(const Handle(Geom_Surface)&)` (`GEOMETRY.cxx:294-325`)
/// reached from `IsQuad(F)`; [`super::geometry::fun_tool_quad_surface`] is the
/// same predicate.
fn fun_tool_quad_surface_checked(s: &dyn occt_geom::Surface) -> bool {
    super::geometry::fun_tool_quad_surface(s)
}

/// `TopOpeBRepTool_TOOL::CurvE(E, par, tg0, curv)` (`TOOL.cxx:964-1016`).
///
/// `None` is the C++ `false` (parameter not usable / the `NYI` general case);
/// `0.0` is attached to the "null curvature" arms, which the C++ reports as
/// `curv = 0.; return true`.
pub fn curv_e(e: &TopoShape, par: f64, tg0: &GpDir) -> Option<f64> {
    let bac = CurveAdaptor::of_edge(&Edge(e.clone()));
    let ct = bac.get_type();
    let line = ct == GeomAbsCurveType::Line;
    let tola = occt_core::precision::ANGULAR * 1.0e3; // NYITOLXPU
    if line {
        let dir = bac.curve()?.gp_line()?.direction();
        let dot = dir.dot(tg0);
        return if (1.0 - dot).abs() >= tola {
            Some(0.0)
        } else {
            None
        };
    }

    let conf = occt_core::precision::CONFUSION;
    let mut clp = CLProps::new(&bac, par, 2, conf).ok()?;
    let mut der_order = 2;
    if !clp.is_tangent_defined(&bac, par, &mut der_order) {
        return None;
    }
    let mut curv = clp.curvature_of(&bac, par, &mut der_order).ok()?.abs();

    let tol = conf * 1.0e2; // NYITOLXPU
    if curv < tol {
        return Some(0.0);
    }
    let n = clp.normal(&bac, par, &mut der_order).ok()?;
    let t = clp.tangent(&bac, par, &mut der_order).ok()?;
    let axis = n.crossed(&t).ok()?;
    let dot = axis.dot(tg0).abs();
    if dot < tola {
        curv = 0.0;
        return Some(curv);
    }
    if (1.0 - dot).abs() < tola {
        return Some(curv);
    }
    None // nyi general case
}

/// `FUN_analyticcS(uv0, S, ngS, tg0, curv, direct)` (`TOOL.cxx:1037-1125`):
/// `Some((curv, direct))` when the apparent contour is analytic.
///
/// The `plane` arm gives `curv = 0`; cylinder / cone / torus read the axis and
/// the `Direct()` flag, sphere reads `1 / R`. The `isMaxAcurv` arm asks
/// `GeomLProp_SLProps` for the curvature whose direction is the axis.
fn fun_analyticc_s(
    uv0: &GpPnt2d,
    s: &dyn occt_geom::Surface,
    ng_s: &GpDir,
    tg0: &GpDir,
    direct: &mut bool,
) -> Option<f64> {
    *direct = true;
    let st = surface_kind(s);
    let plane = st == GeomAbsSurfaceType::Plane;
    let cyl = st == GeomAbsSurfaceType::Cylinder;
    let cone = st == GeomAbsSurfaceType::Cone;
    let sphe = st == GeomAbsSurfaceType::Sphere;
    let torus = st == GeomAbsSurfaceType::Torus;

    let mut curv = 0.0;
    let mut curvdone = false;
    if plane {
        curv = 0.0;
        curvdone = true;
    }
    if cyl || cone || torus {
        let mut axis: Option<GpDir> = None;
        if cyl {
            if let Some(c) = s.gp_cylinder() {
                axis = Some(c.axis().direction().clone());
                *direct = c.is_direct();
            }
        }
        if cone {
            if let Some(c) = s.gp_cone() {
                axis = Some(c.axis().direction().clone());
                *direct = c.is_direct();
            }
        }
        if torus {
            if let Some(t) = s.gp_torus() {
                axis = Some(t.axis().direction().clone());
                *direct = t.is_direct();
            }
        }
        let prod = match &axis {
            Some(a) => a.dot(tg0),
            None => 0.0,
        };
        let is_max_acurv = fun_nullprodv(1.0 - prod.abs());
        let mut nullcurv = fun_nullprodv(prod);

        let prod2 = ng_s.dot(tg0);
        if cyl || cone {
            nullcurv = nullcurv || fun_nullprodv(1.0 - prod2.abs());
        }

        if nullcurv {
            curv = 0.0;
            curvdone = true;
        }
        if is_max_acurv {
            let conf = occt_core::precision::CONFUSION;
            if let Ok(mut slp) =
                occt_geom::geom_lprop::SlProps::new(s, uv0.x(), uv0.y(), 2, conf)
            {
                if slp.is_curvature_defined(s) {
                    let min_acurv = slp.min_curvature(s).unwrap_or(0.0).abs();
                    let max_acurv = slp.max_curvature(s).unwrap_or(0.0).abs();
                    let is_amax = max_acurv > min_acurv;
                    curv = if is_amax { max_acurv } else { min_acurv };
                }
            }
            curvdone = true;
        }
    }
    if sphe {
        if let Some(sp) = s.gp_sphere() {
            curv = 1.0 / sp.radius();
            curvdone = true;
            *direct = sp.is_direct();
        }
    }

    if curvdone {
        Some(curv)
    } else {
        None
    }
}

/// `TopOpeBRepTool_TOOL::CurvF(F, uv, tg0, curv, direct)` (`TOOL.cxx:1131-1187`).
///
/// `None` is the C++ `false` — the basis surface is null, the analytic contour
/// is not recognised and `GeomLProp_SLProps` has no defined curvature.
pub fn curv_f(f: &TopoShape, uv: &GpPnt2d, tg0: &GpDir) -> Option<(f64, bool)> {
    let face = Face(f.clone());
    let ng_s = fun_tool_ng_geom_f(uv, &face)?;
    let s = BRepTool::face_surface_world(&face)?;
    let s = basis_surface(&s);

    let mut direct = true;
    if let Some(curv) = fun_analyticc_s(uv, s.as_ref(), &ng_s, tg0, &mut direct) {
        return Some((curv, direct));
    }

    let tola = 1.0e-6; // NYITOLXPU
    let conf = occt_core::precision::CONFUSION;
    let mut slp = occt_geom::geom_lprop::SlProps::new(s.as_ref(), uv.x(), uv.y(), 2, conf).ok()?;
    if !slp.is_curvature_defined(s.as_ref()) {
        return None;
    }
    let curdef;
    let (max_d, min_d) = slp.curvature_directions(s.as_ref()).ok()?;
    let min_curv = slp.min_curvature(s.as_ref()).ok()?;
    let max_curv = slp.max_curvature(s.as_ref()).ok()?;

    let ng_vec = GpVec::new(ng_s.x(), ng_s.y(), ng_s.z());
    let dmax = GpVec::new(max_d.x(), max_d.y(), max_d.z()).crossed(&ng_vec);
    let dmin = GpVec::new(min_d.x(), min_d.y(), min_d.z()).crossed(&ng_vec);
    let npl = GpVec::new(tg0.x(), tg0.y(), tg0.z());
    let dotmax = dmax.dot(&npl);
    let iscurmax = (1.0 - dotmax).abs() < tola;
    let mut curv = 0.0;
    if iscurmax {
        direct = max_curv < 0.0;
        curv = max_curv.abs();
    }
    let dotmin = dmin.dot(&npl);
    let iscurmin = (1.0 - dotmin).abs() < tola;
    if iscurmin {
        direct = min_curv < 0.0;
        curv = min_curv.abs();
    }
    curdef = iscurmax || iscurmin;
    if curdef {
        Some((curv, direct))
    } else {
        None
    }
}

/// `TopOpeBRepTool_TOOL::minDUV(F)` (`TOOL.cxx:1852-1859`).
pub fn min_duv(f: &Face) -> f64 {
    let bs = SurfaceAdaptor::of_face(f);
    let mut delta = bs.last_u_parameter() - bs.first_u_parameter();
    let tmp = bs.last_v_parameter() - bs.first_v_parameter();
    if tmp < delta {
        delta = tmp;
    }
    delta
}

/// `TopOpeBRepTool_TOOL::stuvF(uv, f, onU, onV)` (`TOOL.cxx:1867-1911`):
/// `INFFIRST` / `ONFIRST` / `ONLAST` / `SUPLAST` / `0` per direction.
pub fn stuv_f(uv: &GpPnt2d, f: &Face) -> (i32, i32) {
    let bs = SurfaceAdaptor::of_face(f);
    let tolf = BRepTool::face_tolerance(f);
    let (tolu, tolv) = match bs.surface() {
        Some(s) => (
            occt_geom::approx_same_parameter::u_resolution(s.as_ref(), tolf),
            occt_geom::approx_same_parameter::v_resolution(s.as_ref(), tolf),
        ),
        None => (tolf, tolf),
    };
    let (u, v) = (uv.x(), uv.y());
    let uf = bs.first_u_parameter();
    let ul = bs.last_u_parameter();
    let vf = bs.first_v_parameter();
    let vl = bs.last_v_parameter();
    let onuf = (uf - u).abs() < tolu;
    let onul = (ul - u).abs() < tolu;
    let onvf = (vf - v).abs() < tolv;
    let onvl = (vl - v).abs() < tolv;
    let mut on_u = 0;
    let mut on_v = 0;
    if onuf {
        on_u = ONFIRST;
    }
    if onul {
        on_u = ONLAST;
    }
    if onvf {
        on_v = ONFIRST;
    }
    if onvl {
        on_v = ONLAST;
    }
    if u < (uf - tolu) {
        on_u = INFFIRST;
    }
    if u > (ul + tolu) {
        on_u = SUPLAST;
    }
    if v < (vf - tolv) {
        on_v = INFFIRST;
    }
    if v > (vl + tolv) {
        on_v = SUPLAST;
    }
    (on_u, on_v)
}

/// `TopOpeBRepTool_TOOL::Getduv(f, uv, dir, factor, duv)` (`TOOL.cxx:1430-1486`).
///
/// `dir` is a 3-D vector translated from `P(uv)` and reprojected on `f`; the
/// uv difference is unwrapped across the periodic direction. `None` is the
/// C++ `false` (non-quadric face, or the reprojection misses `tol3d`).
pub fn getduv(f: &Face, uv: &GpPnt2d, dir: &GpVec, factor: f64) -> Option<GpDir2d> {
    if !is_quad_f(f) {
        return None;
    }
    let p0 = fun_tool_value_uv(uv, f)?;
    let p = p0.translated_vec(&dir.multiplied_scalar(factor));
    let (u_proj, v_proj, d) = fun_tool_proj_pon_f(&p, f)?;
    let tolf = BRepTool::face_tolerance(f) * 1.0e2; // NYIXPUTOL
    if d > tolf {
        return None;
    }
    let uvtr = GpPnt2d::new(u_proj, v_proj);
    let mut dx = uvtr.x() - uv.x();
    let mut dy = uvtr.y() - uv.y();
    let s = BRepTool::face_surface_world(f)?;
    let s = basis_surface(&s);
    let pconf = occt_core::precision::PCONFUSION;
    if s.is_u_periodic() && dx.abs() > s.u_period() / 2.0 {
        let period = s.u_period();
        let mut u1 = uv.x();
        let mut u2 = uvtr.x();
        occt_core::elib::clib::adjust_periodic(0.0, period, pconf, &mut u1, &mut u2);
        dx = u2 - u1;
        if dx > period / 2.0 {
            dx -= period;
        }
    }
    if s.is_v_periodic() && dy.abs() > s.v_period() / 2.0 {
        let period = s.v_period();
        let mut v1 = uv.y();
        let mut v2 = uvtr.y();
        occt_core::elib::clib::adjust_periodic(0.0, period, pconf, &mut v1, &mut v2);
        dy = v2 - v1;
        if dy > period / 2.0 {
            dy -= period;
        }
    }
    GpDir2d::new(dx, dy).ok()
}

/// `TopOpeBRepTool_TOOL::XX(uv, f, par, e, XX)` (`TOOL.cxx:1515-1544`): the
/// in-face normal to `e` at `par`, oriented inside `f`.
///
/// `FUN_tool_getgeomxx(f, e, par, ng)` is `N(uv) ^ tg(e, par)`; a null vector
/// is `false` (`cxx:1528`), and a `REVERSED` occurrence of `e` in `f` flips it.
pub fn xx(uv: &GpPnt2d, f: &TopoShape, par: f64, e: &TopoShape) -> Option<GpDir> {
    let ng = fun_tool_ng_geom_f(uv, &Face(f.clone()))?;
    let geomxx = get_geom_xx(e, par, &ng)?;

    let tol = occt_core::precision::CONFUSION * 1.0e2; // NYITOL
    if geomxx.magnitude() < tol {
        return None;
    }
    let oef = fun_tool_orient_ein_f_forward(e, f)?;
    let mut result = GpDir::from_vec(&geomxx).ok()?;
    if oef == Orientation::Reversed {
        result.reverse();
    }
    Some(result)
}

/// `FUN_tool_getxx(Fi, Ei, parEi, ngFi, XX)` (`TopOpeBRepTool_matter.cxx:106-127`):
/// [`get_geom_xx`] with the `gp::Resolution()` guard, turned into a direction
/// and flipped when `Ei` is `REVERSED` in `Fi` (`FUN_tool_orientEinFFORWARD`).
pub fn get_xx(e: &TopoShape, f: &TopoShape, par_e: f64, ng_f: &GpDir) -> Option<GpDir> {
    let v = get_geom_xx(e, par_e, ng_f)?;
    if v.magnitude() < occt_core::precision::REAL_SMALL {
        return None;
    }
    let mut result = GpDir::from_vec(&v).ok()?;
    if fun_tool_orient_ein_f_forward(e, f)? == Orientation::Reversed {
        result.reverse();
    }
    Some(result)
}

/// `FUN_tool_getgeomxx(Fi, Ei, parEi, ngFi)` (`matter.cxx:66-84`):
/// `FUN_tool_nCinsideS(tgE, ngS) = ngS ^ tgE`.
///
/// The `FUN_tool_nCinsideS` `NYIRAISE` arm (`matter.cxx:80`) is `None`.
pub fn get_geom_xx(e: &TopoShape, par_e: f64, ng_f: &GpDir) -> Option<GpVec> {
    let tg = tg_geom_e(par_e, e)?;
    Some(GpVec::new(ng_f.x(), ng_f.y(), ng_f.z()).crossed(&tg))
}

/// `TopOpeBRepTool_TOOL::Nt(uv, f, normt)` (`TOOL.cxx:1548-1561`): the
/// *topological* normal — the geometric one, flipped for a `REVERSED` face.
pub fn nt(uv: &GpPnt2d, f: &Face) -> Option<GpDir> {
    let n = fun_tool_ng_geom_f(uv, f)?;
    if f.orientation() == Orientation::Reversed {
        Some(n.reversed())
    } else {
        Some(n)
    }
}

/// `TopOpeBRepTool_TOOL::Getstp3dF(p, f, uv, st)` (`TOOL.cxx:1775-1814`).
///
/// `ON` when the projection distance is below the face tolerance, otherwise
/// `IN`/`OUT` from the sign of `(p - proj) . Nt(uv, f)`.
pub fn get_stp_3d_f(p: &GpPnt, f: &Face) -> Option<(GpPnt2d, State)> {
    let tol3d = BRepTool::face_tolerance(f);
    let (u, v, d) = fun_tool_proj_pon_f(p, f)?;
    let uv = GpPnt2d::new(u, v);
    if d < tol3d {
        return Some((uv, State::On));
    }
    let ppr = fun_tool_value_uv(&uv, f)?;
    let ntf = nt(&uv, f)?;
    let dppr = GpVec::from_pnts(p, &ppr);
    let dot = dppr.dot(&GpVec::new(ntf.x(), ntf.y(), ntf.z()));
    let is_out = dot < 0.0;
    Some((uv, if is_out { State::Out } else { State::In }))
}

/// `TopOpeBRepTool_TOOL::EdgeONFace(par, ed, uv, fa, isonfa)`
/// (`TOOL.cxx:1978-2112`): whether the point of `ed` at `par` is on `fa`.
///
/// The prequesitory of the C++ is `pt(par, ed) == pt(uv, fa)`. A degenerated
/// edge is `ON` by definition (`cxx:1986-1991`); tangency is then resolved by
/// plane / cylinder quadric cases and, failing those, by projecting a second
/// point of `ed` on `fa`.
pub fn edge_on_face(par: f64, e: &TopoShape, uv: &GpPnt2d, f: &Face) -> Option<bool> {
    let face = f.clone();
    if BRepTool::is_degenerated(&Edge(e.clone())) {
        return Some(true);
    }
    let tola = occt_core::precision::ANGULAR * 1.0e2; // NYITOLXPU
    let tge = tg_geom_e(par, e)?;
    let ngf = fun_tool_ng_geom_f(uv, &face)?;

    let a_prod_dot = tge.dot(&GpVec::new(ngf.x(), ngf.y(), ngf.z()));
    let etgf = a_prod_dot.abs() < tola;
    if !etgf {
        return Some(false); // "not tangent, thus not ON" is the C++ `return true` with isonfa false
    }

    let bs = SurfaceAdaptor::of_face(&face);
    let st = bs.get_type();
    let plane = st == GeomAbsSurfaceType::Plane;
    let cylinder = st == GeomAbsSurfaceType::Cylinder;

    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let ct = bc.get_type();
    let line = ct == GeomAbsCurveType::Line;
    let circle = ct == GeomAbsCurveType::Circle;

    let tole = bc.tolerance();
    let tol1de = bc.resolution(tole);
    let tolf = BRepTool::face_tolerance(&face);
    let tol3d = tole.max(tolf) * 1.0e2; // NYITOLXPU

    // NYIxpu100299 : for other analytic geometries
    if plane && line {
        return Some(true);
    }
    if plane {
        let mut det = true;
        let ne = if circle {
            bc.curve()?.gp_circ()?.axis().direction().clone()
        } else if ct == GeomAbsCurveType::Ellipse {
            bc.curve()?.gp_ellipse()?.axis().direction().clone()
        } else if ct == GeomAbsCurveType::Hyperbola {
            bc.curve()?.gp_hyperbola()?.axis().direction().clone()
        } else if ct == GeomAbsCurveType::Parabola {
            bc.curve()?.gp_parabola()?.axis().direction().clone()
        } else {
            det = false;
            GpDir::default_dir()
        };
        if det {
            let prod = ne.dot(&ngf);
            return Some((1.0 - prod.abs()).abs() < tola);
        }
    } else if cylinder {
        let mut det = true;
        let ne = if line {
            GpDir::from_vec(&tge).ok()?
        } else if circle {
            bc.curve()?.gp_circ()?.axis().direction().clone()
        } else {
            det = false;
            GpDir::default_dir()
        };
        let axicy = bs.surface()?.gp_cylinder()?.axis().direction().clone();
        if det {
            let prod = ne.dot(&axicy);
            let mut isonfa = (1.0 - prod.abs()).abs() < tola;
            if isonfa && circle {
                let radci = bc.curve()?.gp_circ()?.radius();
                let radcy = bs.surface()?.gp_cylinder()?.radius();
                isonfa = (radci - radcy).abs() < tol3d;
            }
            return Some(isonfa);
        }
    }

    // !!!!!!!!!!!!!!!! NOT STILL OK !!!!!!!!!!!!!!
    // projecting point of <ed> on <fa>
    let x = 0.12345;
    let (fb, lb) = fun_tool_bounds(e);
    let onf = (par - fb).abs() < tol1de;
    let opar = if onf {
        (1.0 - x) * fb + x * lb
    } else {
        (1.0 - x) * fb + x * par
    };
    let opc = bc.d0(opar)?;

    let ouv = fun_tool_par_f_tol(e, opar, &face, tolf)?;
    let ops = fun_tool_value_uv(&ouv, &face)?;
    let dd = opc.distance(&ops);
    Some(dd < tol3d)
}
