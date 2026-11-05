//! The same-parameter family of `BRepFill_Sweep`:
//! `CheckSameParameter` (`cxx:251-283`), `CheckSameParameterExact`
//! (`cxx:288-311`), `SameParameter` (`cxx:319-391`), `CorrectSameParameter`
//! (`cxx:394-448`) and `UpdateEdgeOnPlane` (`cxx:494-517`).
//!
//! Translation notes:
//! - `Adaptor3d_Curve` / `Adaptor3d_Surface` are the port's `Curve` / `Surface`
//!   traits. `GeomAdaptor_Curve(C3d, f, l)` only rebinds the reported range, so
//!   the sampling loops here take the range explicitly and evaluate the curve
//!   itself (`GeomAdaptor_Curve::Value` forwards to `Geom_Curve::Value`).
//! - `Adaptor3d_CurveOnSurface` is [`occt_geom::curve_on_surface::CurveOnSurface`].
//! - `BRep_Builder::UpdateEdge(E, C, Tol)` is
//!   [`GeometryRegistry::set_edge_curve3d`]; `UpdateEdge(E, Tol)` is
//!   [`GeometryRegistry::raise_edge_tolerance`]; `UpdateEdge(E, C2, S, L, Tol)`
//!   is [`GeometryRegistry::update_edge_pcurve_on_surface`].
//! - `BRep_Builder::UpdateVertex(V, Tol)` is
//!   [`GeometryRegistry::raise_vertex_tolerance`].
//! - `BRepAdaptor_Curve::Initialize(E, F)` (`BRepAdaptor_Curve.cxx:117-137`)
//!   builds an `Adaptor3d_CurveOnSurface` from the pcurve and the face's
//!   surface, with `myTrsf` = the *face* location; `Value` applies it, while the
//!   3D curve of `BRep_Tool::Curve` is evaluated without any location.

use std::sync::Arc;

use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::approx::curve_on_surface::ApproxCurveOnSurface;
use occt_geom::approx_same_parameter::ApproxSameParameter;
use occt_geom::check_curve_on_surface::CheckCurveOnSurface;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_check::edge_tolerance;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

use super::helpers;

/// `CheckSameParameter` (`cxx:251-283`): 45 samples of the distance between the
/// 3D curve and the pcurve lifted to the surface.
///
/// `first` / `last` are the range the reference reads off the `C3d` adaptor,
/// i.e. `BRep_Tool::Curve(E, f, l)`'s `f` / `l`. The returned flag is the
/// reference's return value and the second member is `tolreached` (doubled on
/// both exits, `cxx:276-280`).
pub fn check_same_parameter(
    c3d: &dyn Curve,
    first: f64,
    last: f64,
    pcurv: &dyn Curve2d,
    s: &dyn Surface,
    tol3d: f64,
) -> (bool, f64) {
    let nbp = 45;
    let step = 1.0 / (nbp - 1) as f64;
    let mut tolreached = 0.0_f64;
    for i in 0..nbp {
        let t = step * i as f64;
        let t = (1.0 - t) * first + t * last;
        let uv = pcurv.d0(t);
        let p_s = s.d0(uv.x(), uv.y());
        let p_c = c3d.d0(t);
        let d2 = p_s.square_distance(&p_c);
        tolreached = tolreached.max(d2);
    }
    tolreached = tolreached.sqrt();
    if tolreached > tol3d {
        // `cxx:274-278`.
        tolreached *= 2.0;
        return (false, tolreached);
    }
    // `cxx:279-280`.
    tolreached *= 2.0;
    (true, tolreached.max(CONFUSION))
}

/// `CheckSameParameterExact` (`cxx:288-311`): `GeomLib_CheckCurveOnSurface`
/// over the `Adaptor3d_CurveOnSurface`.
///
/// `curve_on_surface` is the `Adaptor3d_CurveOnSurface`; `curve2d` is its
/// `GetCurve()` handle, which [`CheckCurveOnSurface::perform`] takes beside it.
pub fn check_same_parameter_exact(
    c3d: &dyn Curve,
    curve_on_surface: &dyn Curve,
    curve2d: &dyn Curve2d,
    tol3d: f64,
) -> (bool, f64) {
    let mut check = CheckCurveOnSurface::with_curve(c3d, PCONFUSION);
    check.set_parallel(false);
    check.perform(curve_on_surface, curve2d);

    let mut tolreached = check.max_distance();
    if tolreached > tol3d {
        return (false, tolreached);
    }
    // `cxx:305-308`.
    tolreached = tolreached.max(CONFUSION);
    tolreached *= 1.05;
    (true, tolreached)
}

/// `SameParameter` (`cxx:319-391`): makes the pcurve same-parameter with the 3D
/// curve, answering `false` on the two failure exits.
///
/// `pcurv` is the in/out `Handle(Geom2d_Curve)&`: the approximation's rebuilt
/// 2D curve replaces it exactly where `cxx:358-361` / `cxx:384-388` assign
/// `sp.Curve2d()`. The returned flag is the reference's return value and the
/// second member is `tolreached`.
pub fn same_parameter(
    e: &Edge,
    pcurv: &mut Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    tol3d: f64,
) -> (bool, f64) {
    let reg = GeometryRegistry::global();
    let (f, l) = BRepTool::edge_parameters(e);
    let Some(c3d) = reg.edge_curve(&e.0) else {
        // `BRep_Tool::Curve` returned a null handle: the reference dereferences
        // it at `cxx:327`, which the port cannot reproduce.
        return (false, 0.0);
    };

    // `cxx:334-337`.
    let (ok, mut tolreached) =
        check_same_parameter(c3d.as_ref(), f, l, pcurv.as_ref(), surf.as_ref(), tol3d);
    if ok {
        return (true, tolreached);
    }

    // `cxx:339-356`: the pcurve is the only representation of `E`.
    if !helpers::has_pcurves(e) {
        let mut app_curve = ApproxCurveOnSurface::new(
            pcurv.clone(),
            surf.clone(),
            pcurv.first_parameter(),
            pcurv.last_parameter(),
            CONFUSION,
        );
        // `AppCurve.Perform(10, 10, GeomAbs_C1, true)` (`cxx:346`).
        let _ = app_curve.perform(10, 10, 2, true, false);
        if app_curve.is_done() && app_curve.has_result() {
            if let Some(new3d) = app_curve.curve3d() {
                tolreached = app_curve.max_error3d();
                // `B.UpdateEdge(E, C3d, tolreached)` (`cxx:352-353`).
                reg.set_edge_curve3d(&e.0, Arc::new(new3d.clone()), tolreached);
                return (true, tolreached);
            }
        }
    }

    // `cxx:358-370`.
    let sp = ApproxSameParameter::new(c3d.as_ref(), pcurv.as_ref(), surf.as_ref(), f, l, tol3d);
    if sp.done && !sp.same_parameter {
        if let Some(c2d) = &sp.curve2d {
            *pcurv = c2d.clone();
        }
    } else if !sp.done && !sp.same_parameter {
        return (false, tolreached);
    }

    // `cxx:372-383`: `sp.Curve3d()` is the adaptor over the input 3D curve and
    // `sp.CurveOnSurface()` the one over the (possibly rebuilt) pcurve.
    let curve_on_surface = CurveOnSurface::new(pcurv.clone(), surf.clone());
    let (exact_ok, res_tol) =
        check_same_parameter_exact(c3d.as_ref(), &curve_on_surface, pcurv.as_ref(), tol3d);
    if !exact_ok && res_tol > tolreached {
        return (false, tolreached);
    }
    tolreached = 1.1 * res_tol;
    // `cxx:384-388`.
    if sp.done && !sp.same_parameter {
        if let Some(c2d) = &sp.curve2d {
            *pcurv = c2d.clone();
        }
    }
    (true, tolreached)
}

/// `CorrectSameParameter` (`cxx:394-448`): raise the edge tolerance to the
/// largest distance between the 3D curve and either face's pcurve.
///
/// `the_edge` is the reference's `TopoDS_Edge&`: `cxx:415-418` reverses it in
/// place when both faces are the same, so the second `BRepAdaptor_Curve` samples
/// the reversed pcurve. That reversal is a property of the caller's `TopoDS_Edge`
/// view, hence `&mut Edge` here; a null face is `None`.
pub fn correct_same_parameter(
    the_edge: &mut Edge,
    the_face1: Option<&Face>,
    the_face2: Option<&Face>,
) {
    let reg = GeometryRegistry::global();
    if reg.is_degenerated_edge(&the_edge.0) {
        // `cxx:398-401`.
        return;
    }
    let (fpar, lpar) = BRepTool::edge_parameters(the_edge);
    let Some(a_curve) = reg.edge_curve(&the_edge.0) else {
        // `BRep_Tool::Curve` returned a null handle; `cxx:434` would dereference
        // it.
        return;
    };

    // `cxx:407-422`: `BRepAdaptor_Curve::Initialize(theEdge, theFace)`.
    let mut pcurve_exists = [false, false];
    let mut adaptors: [Option<CurveOnSurface>; 2] = [None, None];
    let mut face_trsf = [
        occt_core::gp::GpTrsf::identity(),
        occt_core::gp::GpTrsf::identity(),
    ];
    if let Some(f1) = the_face1 {
        pcurve_exists[0] = true;
        if let Some(a) = face_adaptor(&the_edge.0, f1) {
            adaptors[0] = Some(a.0);
            face_trsf[0] = a.1;
        }
    }
    if let (Some(f1), Some(f2)) = (the_face1, the_face2) {
        if crate::topo_tools::is_same(&f1.0, &f2.0) {
            // `cxx:415-418`: `theEdge.Reverse()`.
            the_edge.0.reverse();
        }
    }
    if let Some(f2) = the_face2 {
        pcurve_exists[1] = true;
        if let Some(a) = face_adaptor(&the_edge.0, f2) {
            adaptors[1] = Some(a.0);
            face_trsf[1] = a.1;
        }
    }

    // `cxx:424-443`.
    const NCONTROL: usize = 23;
    let delta = (lpar - fpar) / NCONTROL as f64;
    let mut max_sq_dist = 0.0_f64;
    for i in 0..=NCONTROL {
        let a_param = fpar + i as f64 * delta;
        let a_pnt = a_curve.d0(a_param);
        for j in 0..2 {
            if !pcurve_exists[j] {
                continue;
            }
            let Some(ad) = adaptors[j].as_ref() else {
                continue;
            };
            let from_face = ad.d0(a_param).transformed(&face_trsf[j]);
            let sq_dist = a_pnt.square_distance(&from_face);
            if sq_dist > max_sq_dist {
                max_sq_dist = sq_dist;
            }
        }
    }

    // `cxx:446-448`: `BB.UpdateEdge(theEdge, sqrt(MaxSqDist))`.
    let a_tol = max_sq_dist.sqrt();
    reg.raise_edge_tolerance(&the_edge.0, a_tol);
}

/// `BRepAdaptor_Curve::Initialize(E, F)` (`BRepAdaptor_Curve.cxx:117-137`): the
/// `Adaptor3d_CurveOnSurface` over `(pcurve on F, surface of F)` plus the face
/// location `myTrsf`.
fn face_adaptor(e: &crate::shape::TopoShape, f: &Face) -> Option<(CurveOnSurface, occt_core::gp::GpTrsf)> {
    let reg = GeometryRegistry::global();
    let surf = BRepTool::face_surface(f)?;
    let pc = reg.edge_pcurve(e, GeometryRegistry::shape_key(&f.0))?;
    Some((
        CurveOnSurface::new(pc, surf),
        f.0.location().transformation(),
    ))
}

/// `UpdateEdgeOnPlane` (`cxx:494-517`, inside the `OCC500(apo)` block).
///
/// `BRep_Tool::CurveOnSurface(E, F, f, l)` is the pcurve of `E` on the surface
/// of `F`; re-attaching it with `BB.UpdateEdge(E, C2d, S, Loc, Tol)` is
/// [`GeometryRegistry::update_edge_pcurve_on_surface`], which removes the
/// representation when the pcurve is null exactly like `cxx:149`.
pub fn update_edge_on_plane(f: &Face, e: &Edge) {
    let reg = GeometryRegistry::global();
    let c2d = reg.edge_pcurve(&e.0, GeometryRegistry::shape_key(&f.0));
    let s = BRepTool::face_surface(f);
    let mut tol = BRepTool::edge_tolerance(e);
    // `BB.UpdateEdge(E, C2d, S, Loc, Tol)` (`cxx:501`).
    if let Some(surf) = &s {
        reg.update_edge_pcurve_on_surface(&e.0, surf, c2d, tol);
    }
    // `Tol = max(Tol, BRepCheck_Edge(E).Tolerance())` (`cxx:502-504`).
    tol = tol.max(edge_tolerance(e));
    reg.raise_edge_tolerance(&e.0, tol);
    tol *= 1.01;
    // `TopExp::FirstVertex` / `LastVertex` (`cxx:506-516`).
    let (first_v, last_v) = crate::topo_tools::edge_vertices(e);
    for v in [first_v, last_v].into_iter().flatten() {
        if BRepTool::vertex_tolerance(&v) < tol {
            reg.raise_vertex_tolerance(&v.0, tol);
        }
    }
}

