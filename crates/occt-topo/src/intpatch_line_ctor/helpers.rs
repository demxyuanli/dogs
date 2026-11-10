//! File-local statics of `IntPatch_LineConstructor`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_LineConstructor.cxx`
//! - `Recadre` (`cxx:44-153`)
//! - `Parameters` (`cxx:156-210`)
//! - `LocalFirstParameter` (`cxx:214-292`)
//! - `LocalLastParameter` (`cxx:296-377`)
//! - `ComputeParametricTolerance` (`cxx:381-406`)
//! - `AppendSameVertexA` (`cxx:416-447`)
//! - `AppendSameVertexG` (`cxx:450-503`)
//! - `AppendSameVertexW` (`cxx:506-546`)
//! - `AppendSameVertexR` (`cxx:549-589`)
//! - `GetVertexTolerance` (`cxx:896-910`)
//! - `IsSegmentSmall` (`cxx:914-946`)
//!
//! Translation notes:
//! - `Recadre` classifies both surfaces with `classify_surface`
//!   (`brep_surface.rs`) instead of `Adaptor3d_Surface::GetType()`, which the
//!   port's `Surface` trait does not expose; `classify_surface` follows the
//!   same `GeomAdaptor_Surface::Load` order, so a rectangular-trimmed surface
//!   still answers with its basis type.
//! - `Recadre` hard-codes `lmf = M_PI + M_PI` (`cxx:97`, `:113`, `:128`,
//!   `:140`); the torus V period is therefore 2π in OCCT too (`cxx:113`), not
//!   `VPeriod()`.
//! - `Parameters` throws `Standard_ConstructionError` on a non-quadric
//!   (`cxx:184`, `:204`); the port answers `Err`.
//! - `LocalFirstParameter` / `LocalLastParameter` read
//!   `MyCurve.IsFirstOpen()` through `IntPatch_ALine::FirstParameter`
//!   (`IntPatch_ALine.lxx:38-50`); the port's `Curve` trait exposes
//!   `is_periodic()` instead, which is how `GeomAdaptor_Curve` answers
//!   `IsFirstOpen` (`!IsClosed()`).
//! - `GetVertexTolerance` is the bare `vtx.Tolerance()` of `cxx:900`: the
//!   `Adaptor3d_TopolTool` lookups are commented out in the source.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{epsilon, INFINITE, PCONFUSION};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::geom_int::{GLineKind, GeomIntLine};
use crate::int_surf_quadric::IntSurfQuadric;
use crate::int_tools_wline::{PatchPoint, WLine};

/// `Recadre(myHS1, myHS2, u1, v1, u2, v2)` (`cxx:44-153`): shift each periodic
/// parameter into the surface's `[First, Last]` range by whole 2π turns.
pub(crate) fn recadre(s1: &dyn Surface, s2: &dyn Surface, uv: &mut (f64, f64, f64, f64)) {
    let (u1_periodic, v1_periodic) = periodic_flags(s1);
    let (u2_periodic, v2_periodic) = periodic_flags(s2);
    let lmf = std::f64::consts::PI * 2.0;

    if u1_periodic {
        let (f, l) = s1.u_range();
        while uv.0 < f {
            uv.0 += lmf;
        }
        while uv.0 > l {
            uv.0 -= lmf;
        }
    }
    if v1_periodic {
        let (f, l) = s1.v_range();
        while uv.1 < f {
            uv.1 += lmf;
        }
        while uv.1 > l {
            uv.1 -= lmf;
        }
    }
    if u2_periodic {
        let (f, l) = s2.u_range();
        while uv.2 < f {
            uv.2 += lmf;
        }
        while uv.2 > l {
            uv.2 -= lmf;
        }
    }
    if v2_periodic {
        let (f, l) = s2.v_range();
        while uv.3 < f {
            uv.3 += lmf;
        }
        while uv.3 > l {
            uv.3 -= lmf;
        }
    }
}

/// The `(U periodic, V periodic)` switch of `Recadre` (`cxx:58-94`): only
/// cylinder / cone / sphere (U) and torus (U and V) count as periodic here.
fn periodic_flags(s: &dyn Surface) -> (bool, bool) {
    match classify_surface(s) {
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere => (true, false),
        SurfaceKind::Torus => (true, true),
        // `//-- Le cas de biparametrees periodiques est gere en amont`
        // (`cxx:70`, `:90`).
        SurfaceKind::Plane | SurfaceKind::Other => (false, false),
    }
}

/// `Parameters(myHS1, myHS2, Ptref, U1, V1, U2, V2)` (`cxx:156-210`): UV of a
/// 3D point on the two quadrics.
pub(crate) fn parameters(
    s1: &dyn Surface,
    s2: &dyn Surface,
    p: &GpPnt,
) -> Result<(f64, f64, f64, f64), String> {
    let quad1 = quadric_of(s1)?;
    let quad2 = quadric_of(s2)?;
    let (u1, v1) = quad1.parameters(p);
    let (u2, v2) = quad2.parameters(p);
    Ok((u1, v1, u2, v2))
}

/// The `IntSurf_Quadric::SetValue` switch of `Parameters` (`cxx:164-208`).
fn quadric_of(s: &dyn Surface) -> Result<IntSurfQuadric, String> {
    let mut q = IntSurfQuadric::new();
    match classify_surface(s) {
        SurfaceKind::Plane => match s.gp_pln() {
            Some(p) => q.set_value_plane(&p),
            None => return Err("IntPatch_LineConstructor::Parameters".to_string()),
        },
        SurfaceKind::Cylinder => match s.gp_cylinder() {
            Some(c) => q.set_value_cylinder(&c),
            None => return Err("IntPatch_LineConstructor::Parameters".to_string()),
        },
        SurfaceKind::Cone => match s.gp_cone() {
            Some(c) => q.set_value_cone(&c),
            None => return Err("IntPatch_LineConstructor::Parameters".to_string()),
        },
        SurfaceKind::Sphere => match s.gp_sphere() {
            Some(sp) => q.set_value_sphere(&sp),
            None => return Err("IntPatch_LineConstructor::Parameters".to_string()),
        },
        SurfaceKind::Torus => match s.gp_torus() {
            Some(t) => q.set_value_torus(&t),
            None => return Err("IntPatch_LineConstructor::Parameters".to_string()),
        },
        SurfaceKind::Other => return Err("IntPatch_IntSS::MakeCurve".to_string()),
    }
    Ok(q)
}

/// `LocalFirstParameter(L)` (`cxx:214-292`).
pub(crate) fn local_first_parameter(l: &GeomIntLine) -> f64 {
    match l {
        GeomIntLine::Analytic(a) => {
            if a.has_first_point {
                return a.first_point().parameter_on_line();
            }
            let f = a.curve.first_parameter();
            let included = a.curve.is_periodic();
            if !included {
                f + epsilon(f)
            } else {
                f
            }
        }
        GeomIntLine::Restriction(r) => {
            if r.has_first_point {
                r.first_point().parameter_on_line()
            } else {
                -INFINITE
            }
        }
        GeomIntLine::Walking(w) => {
            if w.has_first_point {
                w.first_point().parameter_on_line()
            } else {
                1.0
            }
        }
        GeomIntLine::Geometric(g) => {
            if g.has_first_point {
                return g.first_point().parameter_on_line();
            }
            match g.kind {
                GLineKind::Lin(_) | GLineKind::Parab(_) | GLineKind::Hypr(_) => -INFINITE,
                GLineKind::Circ(_) | GLineKind::Elips(_) => 0.0,
            }
        }
    }
}

/// `LocalLastParameter(L)` (`cxx:296-377`).
pub(crate) fn local_last_parameter(l: &GeomIntLine) -> f64 {
    match l {
        GeomIntLine::Analytic(a) => {
            if a.has_last_point {
                return a.last_point().parameter_on_line();
            }
            let last = a.curve.last_parameter();
            let included = a.curve.is_periodic();
            if !included {
                last - epsilon(last)
            } else {
                last
            }
        }
        GeomIntLine::Restriction(r) => {
            if r.has_last_point {
                r.last_point().parameter_on_line()
            } else {
                INFINITE
            }
        }
        GeomIntLine::Walking(w) => {
            if w.has_last_point {
                w.last_point().parameter_on_line()
            } else {
                w.nb_pnts() as f64
            }
        }
        GeomIntLine::Geometric(g) => {
            if g.has_last_point {
                return g.last_point().parameter_on_line();
            }
            match g.kind {
                GLineKind::Lin(_) | GLineKind::Parab(_) | GLineKind::Hypr(_) => INFINITE,
                GLineKind::Circ(_) | GLineKind::Elips(_) => std::f64::consts::PI * 2.0,
            }
        }
    }
}

/// `ComputeParametricTolerance(theTol3d, theD1u, theD1v)` (`cxx:381-406`).
pub(crate) fn compute_parametric_tolerance(tol3d: f64, d1u: &GpVec, d1v: &GpVec) -> f64 {
    let nad1u = d1u.magnitude();
    let nad1v = d1v.magnitude();
    let tolu = if nad1u > 1e-12 { tol3d / nad1u } else { 0.1 };
    let tolv = if nad1v > 1e-12 { tol3d / nad1v } else { 0.1 };
    if tolu > tolv {
        tolu
    } else {
        tolv
    }
}

/// `AppendSameVertexA(alig, L, index, TabIndex)` (`cxx:416-447`). Returns the
/// number of vertices appended; `tab_index` is OCCT's 1-based `TabIndex`.
pub(crate) fn append_same_vertex_a(
    alig: &mut crate::geom_int::ALine,
    l: &crate::geom_int::ALine,
    index: i32,
    tab_index: &mut [i32],
) -> i32 {
    let mut a = 0;
    let n = l.nb_vertex();
    let vtx_index = l.vertex(index).clone();
    let tol1 = vtx_index.tolerance;
    for i in 1..=n {
        if i != index {
            let vtxi = l.vertex(i).clone();
            let mut tol2 = vtxi.tolerance;
            if tol2 < tol1 {
                tol2 = tol1;
            }
            let d = vtx_index.p.distance(&vtxi.p);
            if d <= tol2 {
                alig.add_vertex(vtxi);
                a += 1;
                tab_index[i as usize] = tab_index[index as usize];
            }
        }
    }
    a
}

/// `AppendSameVertexG(glig, L, index, decal, TabIndex)` (`cxx:450-503`).
pub(crate) fn append_same_vertex_g(
    glig: &mut crate::geom_int::GLine,
    l: &crate::geom_int::GLine,
    index: i32,
    decal: f64,
    tab_index: &mut [i32],
) -> i32 {
    let mut a = 0;
    let n = l.nb_vertex();
    let vtx_index = l.vertex(index).clone();
    let tol1 = vtx_index.tolerance;
    for i in 1..=n {
        if i != index {
            let vtxi = l.vertex(i).clone();
            let mut add = false;
            let mut tol2 = vtxi.tolerance;
            if tol2 < tol1 {
                tol2 = tol1;
            }
            let d = vtx_index.p.distance(&vtxi.p);
            if d <= tol2 {
                add = true;
            } else {
                // `//-- Le test suivant a ete ajoute le 20 aout 98`
                // (`cxx:480-489`).
                let p1 = vtx_index.param_on_line;
                let p2 = vtxi.param_on_line;
                if (p1 - p2).abs() < PCONFUSION {
                    add = true;
                }
            }
            if add {
                let p1 = vtx_index.param_on_line;
                let mut a_vtx = vtxi;
                a_vtx.param_on_line = p1 + decal;
                glig.add_vertex(a_vtx);
                a += 1;
                tab_index[i as usize] = tab_index[index as usize];
            }
        }
    }
    a
}

/// `AppendSameVertexW(wlig, L, index, par, TabIndex)` (`cxx:506-546`).
pub(crate) fn append_same_vertex_w(
    wlig: &mut WLine,
    l: &WLine,
    index: i32,
    par: f64,
    tab_index: &mut [i32],
) -> i32 {
    let mut a = 0;
    let n = l.nb_vertex();
    let vtx_index = l.vertex(index).clone();
    let pnt_index = vtx_index.p;
    let tol1 = vtx_index.tolerance;
    let (u1, v1, u2, v2) = vtx_index.parameters();
    for i in 1..=n {
        if i != index {
            let mut vtxi = l.vertex(i).clone();
            let d = pnt_index.distance(&vtxi.p);
            let mut tol2 = vtxi.tolerance;
            if tol2 < tol1 {
                tol2 = tol1;
            }
            if d <= tol2 {
                vtxi.param_on_line = par;
                vtxi.u1 = u1;
                vtxi.v1 = v1;
                vtxi.u2 = u2;
                vtxi.v2 = v2;
                vtxi.p = pnt_index;
                wlig.vertices.push(vtxi);
                a += 1;
                tab_index[i as usize] = tab_index[index as usize];
            }
        }
    }
    a
}

/// `AppendSameVertexR(rlig, L, index, TabIndex)` (`cxx:549-589`).
pub(crate) fn append_same_vertex_r(
    rlig: &mut crate::geom_int::RLine,
    l: &crate::geom_int::RLine,
    index: i32,
    tab_index: &mut [i32],
) -> i32 {
    let mut a = 0;
    let n = l.nb_vertex();
    let vtx_index = l.vertex(index).clone();
    let tol1 = vtx_index.tolerance;
    for i in 1..=n {
        if i != index {
            let vtxi = l.vertex(i).clone();
            let d = vtx_index.p.distance(&vtxi.p);
            let mut tol2 = vtxi.tolerance;
            if tol2 < tol1 {
                tol2 = tol1;
            }
            if d < tol2 {
                if vtxi.param_on_line != vtx_index.param_on_line {
                    let mut vtxi_cop = l.vertex(i).clone();
                    vtxi_cop.param_on_line = vtx_index.param_on_line;
                    rlig.add_vertex(vtxi_cop);
                } else {
                    rlig.add_vertex(vtxi);
                }
                a += 1;
                tab_index[i as usize] = tab_index[index as usize];
            }
        }
    }
    a
}

/// `GetVertexTolerance(vtx)` (`cxx:896-910`).
pub(crate) fn get_vertex_tolerance(vtx: &PatchPoint) -> f64 {
    vtx.tolerance
}

/// `IsSegmentSmall(WLine, ivFirst, ivLast)` (`cxx:914-946`).
pub(crate) fn is_segment_small(wline: &WLine, iv_first: i32, iv_last: i32) -> bool {
    let vtx_f = wline.vertex(iv_first);
    let vtx_l = wline.vertex(iv_last);
    let ip_f = vtx_f.param_on_line as i32;
    let ip_l = vtx_l.param_on_line as i32;
    if ip_f >= ip_l {
        return true;
    }

    let tol_f = get_vertex_tolerance(vtx_f);
    let tol_l = get_vertex_tolerance(vtx_l);
    let tol = if tol_f > tol_l { tol_f } else { tol_l };

    let mut len = 0.0;
    let mut p1 = wline.point(ip_f).value();
    for i in (ip_f + 1)..=ip_l {
        let p2 = wline.point(i).value();
        len += p1.distance(&p2);
        if len > tol {
            break;
        }
        p1 = p2;
    }
    len <= tol
}

