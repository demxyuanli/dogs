//! `BRepBuilderAPI_GTransform::Perform` (`BRepBuilderAPI_GTransform.cxx:47-57`)
//! on one face: `BRepBuilderAPI_NurbsConvert` then
//! `BRepTools_GTrsfModification` pole updates.
//!
//! Surface conversion follows `GeomConvert::SurfaceToBSplineSurface` for a
//! trimmed plane (`GeomConvert_1.cxx:410-432`), a cylinder, cone, sphere, or
//! torus (`:433-543`), a surface of revolution (`:544-644`), a trimmed
//! linear extrusion (`:646-691`), a trimmed Bezier (`:693-730`), and a
//! trimmed B-spline (`:732-768`). A shifted full-period seam uses
//! `CheckAndSegment` (`Geom_BSplineSurface.cxx:880`). A non-periodic edge
//! whose knots miss the edge range is rescaled by `BSplCLib::Reparametrize`
//! (`BRepTools_NurbsConvertModification.cxx:469-476`). An infinite basis that
//! is not elementary and not an extrusion reaches `Iso::MakeApprox`.
//! `Discretise` cuts with `UpdateInU` and `UpdateInV`. A finished square is
//! stored. `Compute3DErrors` runs, then `ConvertBS` returns the B-spline.
//! The trimmed fallthrough (`GeomConvert_1.cxx:771`) and the direct
//! fallthrough (`:935`) call `GeomConvert_ApproxSurface` with tolerance
//! `1e-4`, degree 14, and `PrecisCode` 1. A thrown approximation decrements
//! `GeomAbs_Shape` and tries again (`ShapeConstruct.cxx:301`). A failed iso
//! cuts with `PrefAndRec::Value`. `CutSense` 0 keeps a stored square
//! (`AdvApp2Var_ApproxAFunc2Var.cxx:442`). `D3` stores `D3UUV` when
//! `UOrder == 2` and `D3UVV` when `UOrder == 1`
//! (`GeomConvert_ApproxSurface.cxx:193`). `DN` stores `EvalDN` (`:212`).
//! Offset uses `EvaluateD3` and `EvaluateDN`. `mmapptt_` orders -1, 0, 1,
//! and 2 read `MMAPGSS`, `MMAPGS0`, `MMAPGS1`, and `MMAPGS2`. Jacobi order
//! -1 truncates with `mmtrpj0_` (`AdvApp2Var_MathBase.cxx:9006`).
//! `mmfmtb1_` transposes the IsoU tables
//! (`AdvApp2Var_ApproxF2var.cxx:5634`). `MakeApprox` stays on one 3D
//! subspace (`AdvApp2Var_Iso.cxx:122`). The rational Bezier `Segment` arm
//! runs when weights are stored (`Geom_BezierSurface.cxx:961`).
//! `EvalD0`, `EvalD1`, `EvalD2`, and `EvalDN` use the rational B-spline
//! evaluator
//! (`Geom_BezierSurface.cxx:1425`). `UIso` and `VIso` call `BSplSLib::Iso`
//! and keep curve weights when that direction is rational
//! (`Geom_BezierSurface.cxx:1776`). A trimmed rational Bezier becomes a
//! rational B-spline (`GeomConvert_1.cxx:714`). Curve `Segment` trims
//! weights with the poles (`Geom_BezierCurve.cxx:395`). `SetWeight` turns a
//! curve or surface rational (`Geom_BezierCurve.cxx:455`,
//! `Geom_BezierSurface.cxx:1179`). `SetWeightCol` and `SetWeightRow`
//! rewrite one column or a row span (`Geom_BezierSurface.cxx:1219`,
//! `:1262`). `UReverse` and `VReverse` swap poles and weights
//! (`:1305`, `:1352`). `ExchangeUV` transposes poles and weights
//! (`Geom_BezierSurface.cxx:490`). `Increase` elevates degree through
//! `BSplSLib::IncreaseDegree` (`:532`). Insert and remove pole columns
//! and rows (`Geom_BezierSurface.cxx:652-937`). `SetPole` and column/row
//! setters (`:1047-1175`). Weighted construction (`:416`). Curve
//! `Increase`, insert/remove poles, and `Reverse`
//! (`Geom_BezierCurve.cxx:141`, `:210`, `:298`, `:352`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpGTrsf, GpPnt, GpTrsf, GpVec, GpXyz};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{PCONFUSION, CONFUSION};
use occt_geom::bspline_curve::{GeomBSplineCurve, MAX_DEGREE};
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::convert_approx_curve::GeomConvertApproxCurve;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::surface_of_revolution::GeomSurfaceOfRevolution;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom::{Curve, Surface};

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::{FaceGeom, GeometryRegistry};
use crate::topo_tools::is_same;

pub(in crate::iges_read) fn gtransform_face(face: &mut TopoShape, g: &GpGTrsf) -> Result<(), String> {
    let scale = g_scale(g);
    let reg = GeometryRegistry::global();
    let surf = reg
        .face_surface(face)
        .ok_or_else(|| "GTransform face has no surface: BRepTools_GTrsfModification.cxx:82-85".to_string())?;
    let mut bspline = surface_as_bspline(face, surf.as_ref())?;
    transform_surface_poles(&mut bspline, g);
    let tol = reg.face_tolerance(face) * scale;
    reg.set_face(face, FaceGeom::new(Arc::new(bspline)));
    reg.set_face_tolerance(face, tol);
    let wires = face.tshape.read().expect("lock").children.clone();
    let mut seen_edges = Vec::new();
    let mut seen_verts = Vec::new();
    for wire in wires {
        let edges = wire.tshape.read().expect("lock").children.clone();
        for edge in edges {
            if !edge.is_edge() || seen_edges.iter().any(|old: &TopoShape| is_same(old, &edge)) {
                continue;
            }
            seen_edges.push(edge.clone());
            transform_edge(&Edge(edge), g, scale, &mut seen_verts)?;
        }
    }
    // `RevFace = myGTrsf.IsNegative()` (`BRepTools_GTrsfModification.cxx:92`).
    if g.is_negative() {
        face.reverse();
    }
    Ok(())
}

fn g_scale(g: &GpGTrsf) -> f64 {
    // `BRepTools_GTrsfModification.cxx:47-59`.
    let mut loc1 = g.value(1, 1).abs().max(g.value(1, 2).abs());
    let loc2 = g.value(2, 1).abs().max(g.value(2, 2).abs());
    let loc3 = g.value(3, 1).abs().max(g.value(3, 2).abs());
    let loc4 = g.value(1, 3).abs().max(g.value(2, 3).abs());
    loc1 = loc1.max(loc2);
    loc1 = loc1.max(loc3.max(loc4));
    loc1.max(g.value(3, 3).abs())
}

fn transform_xyz(g: &GpGTrsf, p: &GpPnt) -> GpPnt {
    let mut xyz = *p.xyz();
    g.transforms(&mut xyz);
    GpPnt::from_xyz(&xyz)
}

fn transform_surface_poles(surf: &mut GeomBSplineSurface, g: &GpGTrsf) {
    for row in &mut surf.poles {
        for pole in row {
            *pole = transform_xyz(g, pole);
        }
    }
}

fn surface_as_bspline(face: &TopoShape, surf: &dyn Surface) -> Result<GeomBSplineSurface, String> {
    if surf.is_bspline_surface() || surf.is_bezier_surface() {
        return surf.osculating_bspline().ok_or_else(|| {
            "GTransform BSpline copy failed: BRepTools_GTrsfModification.cxx:96-127".into()
        });
    }
    let (basis, u1, u2, v1, v2) = trim_window(face, surf);
    if !u1.is_finite() || !u2.is_finite() || !v1.is_finite() || !v2.is_finite() {
        return Err(
            "GeomConvert::SurfaceToBSplineSurface infinite: GeomConvert_1.cxx:354-357".into(),
        );
    }
    if basis.gp_pln().is_some() && !basis.is_surface_of_revolution() {
        return plane_patch(basis.as_ref(), u1, u2, v1, v2);
    }
    if let Some(built) = super::elementary_bspline::elementary_patch(basis.as_ref(), u1, u2, v1, v2)?
    {
        return Ok(built);
    }
    if basis.is_surface_of_revolution() {
        return revolution_bspline(basis.as_ref(), u1, u2, v1, v2);
    }
    if surf.rectangular_trimmed_basis().is_some() && basis.is_surface_of_linear_extrusion() {
        return extrusion_bspline(surf, basis.as_ref(), u1, u2, v1, v2);
    }
    if surf.rectangular_trimmed_basis().is_some() && basis.is_bezier_surface() {
        return trimmed_bezier_bspline(basis.as_ref(), u1, u2, v1, v2);
    }
    if surf.rectangular_trimmed_basis().is_some() && basis.is_bspline_surface() {
        return trimmed_bspline(basis.as_ref(), u1, u2, v1, v2);
    }
    if surf.rectangular_trimmed_basis().is_some() {
        // `GeomConvert_1.cxx:771-787`: `Tol3d` 1e-4, degree 14, `PrecisCode` 1.
        // C2 intervals decide one continuity for both directions.
        const C2: u8 = 4;
        let order = if surf.nb_u_intervals(C2) > 1 || surf.nb_v_intervals(C2) > 1 {
            1
        } else {
            2
        };
        return geomconvert_approx(surf, u1, u2, v1, v2, order, order, max_seg_cn(surf)?);
    }
    // `GeomConvert_1.cxx:935-961`: same tolerance and segment formula.
    // `IsCNu` and `IsCNv` choose the two orders independently.
    let order_u = if surf.is_cnu(2)? {
        2
    } else if surf.is_cnu(1)? {
        1
    } else {
        0
    };
    let order_v = if surf.is_cnv(2)? {
        2
    } else if surf.is_cnv(1)? {
        1
    } else {
        0
    };
    geomconvert_approx(surf, u1, u2, v1, v2, order_u, order_v, max_seg_cn(surf)?)
}

/// `GeomConvert_1.cxx:785` and `:959`: `4 * (NbUIntervals(CN) + 1) * (NbVIntervals(CN) + 1)`.
fn max_seg_cn(surf: &dyn Surface) -> Result<i32, String> {
    const CN: u8 = 6;
    let nu = surf.nb_u_intervals(CN).checked_add(1).ok_or_else(|| {
        "GeomConvert_1.cxx:785".to_string()
    })?;
    let nv = surf.nb_v_intervals(CN).checked_add(1).ok_or_else(|| {
        "GeomConvert_1.cxx:785".to_string()
    })?;
    nu.checked_mul(nv)
        .and_then(|p| p.checked_mul(4))
        .ok_or_else(|| "GeomConvert_1.cxx:785".to_string())
}

/// `GeomConvert_ApproxSurface` at `GeomConvert_1.cxx:786` and `:960`.
fn geomconvert_approx(
    surf: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    order_u: i32,
    order_v: i32,
    max_seg: i32,
) -> Result<GeomBSplineSurface, String> {
    approximate_prefix(surf, u1, u2, v1, v2, order_u, order_v, 14, max_seg, 1, 1.0e-4)
}

/// `ShapeConstruct::ConvertSurfaceToBSpline` (`ShapeConstruct.cxx:159-309`).
///
/// An elementary basis is re-trimmed and sent through
/// `GeomConvert::SurfaceToBSplineSurface` (`:181-185`).
pub(in crate::iges_read) fn convert_surface_to_bspline(
    surf: &dyn Surface,
    uf: f64,
    ul: f64,
    vf: f64,
    vl: f64,
) -> Result<GeomBSplineSurface, String> {
    let basis = match surf.rectangular_trimmed_basis() {
        Some(basis) => basis,
        None => Arc::from(surf.clone_dyn()),
    };
    let s = basis.as_ref();
    let elementary = (s.gp_pln().is_some() && !s.is_surface_of_revolution())
        || s.gp_cylinder().is_some()
        || s.gp_cone().is_some()
        || s.gp_sphere().is_some()
        || s.gp_torus().is_some();
    if elementary {
        if s.gp_pln().is_some() && !s.is_surface_of_revolution() {
            return plane_patch(s, uf, ul, vf, vl);
        }
        return super::elementary_bspline::elementary_patch(s, uf, ul, vf, vl)?.ok_or_else(|| {
            "ConvertSurfaceToBSpline elementary: ShapeConstruct.cxx:185".to_string()
        });
    }
    if s.is_surface_of_linear_extrusion() {
        // `ShapeAlgo_AlgoContainer.cxx:410-418`: Confusion, C1, 100 segments,
        // `Geom_BSplineSurface::MaxDegree`.
        return infinite_extrusion_bspline(s, uf, ul, vf, vl);
    }
    // `ShapeConstruct.cxx:238-255`: a revolution whose generatrix is an offset
    // is rebuilt from `ConvertCurveToBSpline` before the approximation.
    let mut current = basis;
    if current.is_surface_of_revolution() {
        if let Some(generatrix) = current.revolution_basis_curve() {
            if generatrix.offset_curve().is_some() {
                let mut cnt = generatrix.continuity();
                if cnt > Shape::C2 as u8 {
                    cnt = Shape::C2 as u8;
                }
                let axis = current.revolution_axis().ok_or_else(|| {
                    "revolution has no axis: ShapeConstruct.cxx:248".to_string()
                })?;
                let bspl = convert_curve_to_bspline(
                    generatrix.as_ref(),
                    vf,
                    vl,
                    CONFUSION,
                    shape_from_continuity(cnt),
                    100,
                    MAX_DEGREE as i32,
                )?;
                current = Arc::new(GeomSurfaceOfRevolution::new(Arc::new(bspl), axis));
            }
        }
    }
    // `ShapeConstruct.cxx:257-309`: start at C1, and on a thrown approximation
    // decrement the `GeomAbs_Shape` and try again (`:301-305`).
    shapeconstruct_approx(current.as_ref(), uf, ul, vf, vl)
}

/// `ShapeConstruct.cxx:259-309`. `Continuity` is `GeomAbs_C1` (value 2).
/// `G1` and `G2` throw `UContinuity Error` (`AdvApp2Var_ApproxAFunc2Var.cxx:191`)
/// and the catch decrements again. `cnt == 0` does not decrement
/// (`ShapeConstruct.cxx:301-305`); the error is returned so that spin ends.
fn shapeconstruct_approx(
    surf: &dyn Surface,
    uf: f64,
    ul: f64,
    vf: f64,
    vl: f64,
) -> Result<GeomBSplineSurface, String> {
    let mut cnt: i32 = Shape::C1 as i32;
    loop {
        let order = match cnt {
            x if x == Shape::C0 as i32 => 0,
            x if x == Shape::C1 as i32 => 1,
            x if x == Shape::C2 as i32 => 2,
            _ => {
                if cnt > 0 {
                    cnt -= 1;
                    continue;
                }
                return Err("AdvApp2Var_ApproxAFunc2Var : UContinuity Error".into());
            }
        };
        match approximate_prefix(
            surf,
            uf,
            ul,
            vf,
            vl,
            order,
            order,
            MAX_DEGREE as i32,
            100,
            0,
            CONFUSION * 0.5,
        ) {
            Ok(built) => return Ok(built),
            Err(_) if cnt > 0 => {
                cnt -= 1;
                continue;
            }
            Err(err) => return Err(err),
        }
    }
}

fn shape_from_continuity(value: u8) -> Shape {
    match value {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        _ => Shape::CN,
    }
}

/// `AdvApp2Var_ApproxAFunc2Var::Init` (`AdvApp2Var_ApproxAFunc2Var.cxx:164-234`).
/// `order_u` / `order_v` are 0, 1, or 2 (`C0` / `C1` / `C2`). IsoV (`ifav = 2`).
fn approximate_prefix(
    surf: &dyn Surface,
    uf: f64,
    ul: f64,
    vf: f64,
    vl: f64,
    order_u: i32,
    order_v: i32,
    max_degree: i32,
    max_segments: i32,
    precis: i32,
    tol: f64,
) -> Result<GeomBSplineSurface, String> {
    let iu = order_u;
    let iv = order_v;
    let ndu = (max_degree + 1).max(2 * iu + 2);
    let ndv = (max_degree + 1).max(2 * iv + 2);
    if ndu < 2 * iu + 2 {
        return Err("AdvApp2Var_ApproxAFunc2Var : UMaxDegree Error".into());
    }
    if ndv < 2 * iv + 2 {
        return Err("AdvApp2Var_ApproxAFunc2Var : VMaxDegree Error".into());
    }
    let precis = precis.clamp(0, 3);
    // `GeomConvert_ApproxSurface::Approximate` (`GeomConvert_ApproxSurface.cxx:360-376`):
    // one 3D subspace. The caller passes `Tol3d / 2` or the `GeomConvert_1` `1e-4`.
    let frontier = [tol, tol, tol, tol];
    let ctx = occt_geom::adv_app2var::context_until_roots(
        2, iu, iv, ndu, ndv, precis, 0, 0, 1, &[], &[], &[tol], &[], &[], &[frontier],
    )?;
    occt_geom::adv_app2var::approximate_until_cdi(surf, uf, ul, vf, vl, &ctx, max_segments)
}

fn infinite_extrusion_bspline(
    extru: &dyn Surface,
    uf: f64,
    ul: f64,
    vf: f64,
    vl: f64,
) -> Result<GeomBSplineSurface, String> {
    // `ShapeConstruct.cxx:188-235`. The generatrix is converted on `[UF, UL]`,
    // then each pole row is translated from `Value(UF, 0)` to `Value(UF, V)`.
    let meridian = extru.extrusion_basis_curve().ok_or_else(|| {
        "extrusion has no generatrix: ShapeConstruct.cxx:192".to_string()
    })?;
    let curve = convert_curve_to_bspline(meridian.as_ref(), uf, ul, CONFUSION, Shape::C1, 100, MAX_DEGREE as i32)?;
    let origin = extru.d0(uf, 0.0);
    let at_first = extru.d0(uf, vf);
    let at_last = extru.d0(uf, vl);
    let mut shift_first = GpTrsf::identity();
    shift_first.set_translation_pnts(&origin, &at_first);
    let mut shift_last = GpTrsf::identity();
    shift_last.set_translation_pnts(&origin, &at_last);
    let (u_knots, u_mults) = curve.distinct_knots_and_mults();
    let nu = curve.poles.len();
    let mut poles = Vec::with_capacity(nu);
    let mut weights = Vec::with_capacity(nu);
    for i in 0..nu {
        let w = curve.weight((i as i32) + 1);
        poles.push(vec![
            curve.poles[i].transformed(&shift_first),
            curve.poles[i].transformed(&shift_last),
        ]);
        weights.push(vec![w, w]);
    }
    let period = u_knots.last().copied().unwrap_or(0.0) - u_knots.first().copied().unwrap_or(0.0);
    let flat_u = flat_knots(&u_knots, &u_mults, curve.degree as i32, curve.periodic, period);
    let flat_v = knot_sequence(&[vf, vl], &[2, 2], 1);
    let mut built = GeomBSplineSurface::rational(poles, weights, flat_u, flat_v, curve.degree, 1)
        .map_err(|e| format!("infinite extrusion BSpline: {e}"))?;
    built.u_periodic = curve.periodic;
    built.v_periodic = false;
    Ok(built)
}

/// `ShapeConstruct::ConvertCurveToBSpline` (`ShapeConstruct.cxx:61-111`).
fn convert_curve_to_bspline(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    tol: f64,
    continuity: Shape,
    max_segments: i32,
    mut max_degree: i32,
) -> Result<GeomBSplineCurve, String> {
    if curve.bspline_poles().is_some() && !curve.is_geom_trimmed() {
        return Ok(GeomBSplineCurve {
            poles: curve.bspline_poles().unwrap().to_vec(),
            weights: curve.bspline_weights().map(|w| w.to_vec()),
            knots: curve.bspline_knots().unwrap_or(&[]).to_vec(),
            degree: curve.nurbs_degree().unwrap_or(1),
            periodic: curve.is_periodic(),
        });
    }
    if curve.gp_circ().is_some()
        || curve.gp_ellipse().is_some()
        || curve.gp_hyperbola().is_some()
        || curve.gp_parabola().is_some()
    {
        max_degree = max_degree.min(6);
    }
    let trimmed = GeomTrimmedCurve::new(Arc::from(curve.clone_dyn()), first, last);
    let approx = GeomConvertApproxCurve::new(&trimmed, tol, continuity, max_segments, max_degree);
    if approx.has_result() {
        return approx.curve().cloned().ok_or_else(|| {
            "GeomConvert_ApproxCurve::Curve: ShapeConstruct.cxx:92".to_string()
        });
    }
    curve_to_bspline_curve(curve, ParameterisationType::QuasiAngular)
        .map_err(|_| "CurveToBSplineCurve: ShapeConstruct.cxx:96".to_string())
}

fn trim_window(face: &TopoShape, surf: &dyn Surface) -> (Arc<dyn Surface>, f64, f64, f64, f64) {
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        let (u1, u2) = surf.u_range();
        let (v1, v2) = surf.v_range();
        return (basis, u1, u2, v1, v2);
    }
    let (mut u1, mut u2, mut v1, mut v2) = crate::brep_tools::uv_bounds(&Face(face.clone()));
    let (su1, su2) = surf.u_range();
    let (sv1, sv2) = surf.v_range();
    let tol_par = 0.1 * GeometryRegistry::global().face_tolerance(face);
    if (u1 - su1).abs() <= tol_par {
        u1 = su1;
    }
    if (u2 - su2).abs() <= tol_par {
        u2 = su2;
    }
    if (v1 - sv1).abs() <= tol_par {
        v1 = sv1;
    }
    if (v2 - sv2).abs() <= tol_par {
        v2 = sv2;
    }
    if !surf.is_u_periodic() {
        if su1.is_finite() {
            u1 = u1.max(su1);
        }
        if su2.is_finite() {
            u2 = u2.min(su2);
        }
    }
    if !surf.is_v_periodic() {
        if sv1.is_finite() {
            v1 = v1.max(sv1);
        }
        if sv2.is_finite() {
            v2 = v2.min(sv2);
        }
    }
    (Arc::from(surf.clone_dyn()), u1, u2, v1, v2)
}

fn plane_patch(surf: &dyn Surface, u1: f64, u2: f64, v1: f64, v2: f64) -> Result<GeomBSplineSurface, String> {
    // `GeomConvert_1.cxx:410-432`.
    let poles = vec![
        vec![surf.d0(u1, v1), surf.d0(u1, v2)],
        vec![surf.d0(u2, v1), surf.d0(u2, v2)],
    ];
    let knots_u = knot_sequence(&[u1, u2], &[2, 2], 1);
    let knots_v = knot_sequence(&[v1, v2], &[2, 2], 1);
    GeomBSplineSurface::new(poles, knots_u, knots_v, 1, 1)
        .map_err(|e| format!("plane BSpline: {e}"))
}

fn trimmed_bezier_bspline(
    basis: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    // `GeomConvert_1.cxx:693-730`. `Segment` reparameterizes the patch onto
    // `[0, 1] x [0, 1]`, then the B-spline knots are 0 and 1.
    let (poles, n_u, n_v) = basis.bezier_control_grid().ok_or_else(|| {
        "trimmed Bezier has no poles: GeomConvert_1.cxx:696".to_string()
    })?;
    let mut bez = occt_geom::bezier_surface::GeomBezierSurface::new(poles, n_u, n_v)
        .map_err(|e| format!("trimmed Bezier: {e}"))?;
    if let Some((weights, u_rat, v_rat)) = basis.bezier_weights() {
        bez.adopt_weights(weights, u_rat, v_rat);
    }
    bez.segment(u1, u2, v1, v2);
    bez.converted_bspline()
        .ok_or_else(|| "trimmed Bezier BSpline: GeomConvert_1.cxx:714-729".to_string())
}

fn trimmed_bspline(
    basis: &dyn Surface,
    mut u1: f64,
    mut u2: f64,
    mut v1: f64,
    mut v2: f64,
) -> Result<GeomBSplineSurface, String> {
    // `GeomConvert_1.cxx:732-768`.
    let mut bs = basis.osculating_bspline().ok_or_else(|| {
        "trimmed BSpline copy failed: GeomConvert_1.cxx:734".to_string()
    })?;
    let (umin, umax) = basis.u_range();
    let (vmin, vmax) = basis.v_range();
    if !basis.is_u_periodic() {
        if u1 < umin {
            u1 = umin;
        }
        if u2 > umax {
            u2 = umax;
        }
    }
    if !basis.is_v_periodic() {
        if v1 < vmin {
            v1 = vmin;
        }
        if v2 > vmax {
            v2 = vmax;
        }
    }
    if basis.is_u_periodic() || basis.is_v_periodic() {
        bs.check_and_segment(u1, u2, v1, v2)?;
    } else {
        bs.segment_uv(u1, u2, v1, v2, PCONFUSION, PCONFUSION)?;
    }
    Ok(bs)
}

fn extrusion_bspline(
    trimmed: &dyn Surface,
    extru: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    // `GeomConvert_1.cxx:646-691`. U is the basis curve, V is the extrusion.
    let meridian = extru.extrusion_basis_curve().ok_or_else(|| {
        "extrusion has no generatrix: GeomConvert_1.cxx:651".to_string()
    })?;
    let direction = extru
        .extrusion_direction()
        .ok_or_else(|| "extrusion has no direction: GeomConvert_1.cxx:672".to_string())?;
    let curve = if trimmed.is_u_closed() {
        curve_to_bspline_curve(meridian.as_ref(), ParameterisationType::TgtThetaOver2)
    } else {
        let span = GeomTrimmedCurve::new(meridian, u1, u2);
        curve_to_bspline_curve(&span, ParameterisationType::TgtThetaOver2)
    }
    .map_err(|_| "CurveToBSplineCurve meridian: GeomConvert_1.cxx:655-660".to_string())?;
    let (u_knots, u_mults) = curve.distinct_knots_and_mults();
    let along = GpVec::from_xyz(direction.xyz());
    let dv1 = along.multiplied_scalar(v1);
    let dv2 = along.multiplied_scalar(v2);
    let nu = curve.poles.len();
    let mut poles = Vec::with_capacity(nu);
    let mut weights = Vec::with_capacity(nu);
    for i in 0..nu {
        let w = curve.weight((i as i32) + 1);
        poles.push(vec![
            curve.poles[i].translated_vec(&dv1),
            curve.poles[i].translated_vec(&dv2),
        ]);
        weights.push(vec![w, w]);
    }
    let period = u_knots.last().copied().unwrap_or(0.0) - u_knots.first().copied().unwrap_or(0.0);
    let flat_u = flat_knots(&u_knots, &u_mults, curve.degree as i32, curve.periodic, period);
    let flat_v = knot_sequence(&[v1, v2], &[2, 2], 1);
    let mut built = GeomBSplineSurface::rational(poles, weights, flat_u, flat_v, curve.degree, 1)
        .map_err(|e| format!("extrusion BSpline: {e}"))?;
    built.u_periodic = curve.periodic;
    built.v_periodic = false;
    Ok(built)
}

fn revolution_bspline(
    surf: &dyn Surface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<GeomBSplineSurface, String> {
    // `GeomConvert_1.cxx:544-644`.
    let meridian = surf
        .revolution_basis_curve()
        .ok_or_else(|| "revolution has no generatrix: GeomConvert_1.cxx:548".to_string())?;
    let axis = surf
        .revolution_axis()
        .ok_or_else(|| "revolution has no axis: GeomConvert_1.cxx:607".to_string())?;
    let v_closed = v_span_closed(surf, v1, v2);
    let curve = if v_closed {
        curve_to_bspline_curve(meridian.as_ref(), ParameterisationType::TgtThetaOver2)
    } else {
        let trimmed = GeomTrimmedCurve::new(meridian, v1, v2);
        curve_to_bspline_curve(&trimmed, ParameterisationType::TgtThetaOver2)
    }
    .map_err(|_| "CurveToBSplineCurve meridian: GeomConvert_1.cxx:552-557".to_string())?;
    let u_periodic = surf.is_u_periodic() && (u2 - u1 - surf.u_period()).abs() <= PCONFUSION;
    let (nb_u, nb_knots, alfa, periodic) = if u_periodic {
        (6usize, 4usize, PI / 3.0, true)
    } else {
        let spans = (1.2 * (u2 - u1) / PI).trunc() as usize + 1;
        let alfa = (u2 - u1) / (spans as f64 * 2.0);
        (2 * spans + 1, spans + 1, alfa, false)
    };
    let mut u_knots = Vec::with_capacity(nb_knots);
    let mut u_mults = vec![2i32; nb_knots];
    for i in 0..nb_knots {
        u_knots.push(u1 + i as f64 * 2.0 * alfa);
    }
    if !periodic {
        u_mults[0] += 1;
        u_mults[nb_knots - 1] += 1;
    }
    let nv = curve.poles.len();
    let mut meridian_poles = curve.poles.clone();
    let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); nv]; nb_u];
    let mut weights = vec![vec![1.0; nv]; nb_u];
    let mut iu = 0usize;
    while iu < nb_u {
        let mut rot = GpTrsf::identity();
        rot.set_rotation_ax1(&axis, u1 + iu as f64 * alfa)
            .map_err(|e| e.to_string())?;
        for j in 0..nv {
            poles[iu][j] = meridian_poles[j].transformed(&rot);
            weights[iu][j] = curve_weight(&curve, j);
        }
        iu += 2;
    }
    let mut aff = GpGTrsf::identity();
    aff.set_affinity_ax1(&axis, 1.0 / alfa.cos());
    for pole in &mut meridian_poles {
        *pole = transform_xyz(&aff, pole);
    }
    let mut iu = 1usize;
    while iu < nb_u {
        let mut rot = GpTrsf::identity();
        rot.set_rotation_ax1(&axis, u1 + iu as f64 * alfa)
            .map_err(|e| e.to_string())?;
        for j in 0..nv {
            poles[iu][j] = meridian_poles[j].transformed(&rot);
            weights[iu][j] = curve_weight(&curve, j) * alfa.cos();
        }
        iu += 2;
    }
    let knots_u = flat_knots(&u_knots, &u_mults, 2, periodic, u2 - u1);
    let knots_v = curve.knots.clone();
    let mut built = GeomBSplineSurface::rational(poles, weights, knots_u, knots_v, 2, curve.degree)
        .map_err(|e| format!("revolution BSpline: {e}"))?;
    built.u_periodic = periodic;
    built.v_periodic = curve.periodic;
    Ok(built)
}

fn curve_weight(curve: &GeomBSplineCurve, index: usize) -> f64 {
    curve.weights.as_ref().map(|w| w[index]).unwrap_or(1.0)
}

fn v_span_closed(surf: &dyn Surface, v1: f64, v2: f64) -> bool {
    if !surf.is_v_periodic() {
        return false;
    }
    let length = v2 - v1;
    let period = surf.v_period();
    length > PCONFUSION && (length - (length / period).round() * period).abs() <= PCONFUSION
}

pub(super) fn flat_knots(knots: &[f64], mults: &[i32], degree: i32, periodic: bool, period: f64) -> Vec<f64> {
    // `BSplCLib::KnotSequence` (`BSplCLib.cxx:488-539`).
    let body = knot_sequence(knots, mults, degree);
    if !periodic || knots.len() < 2 {
        return body;
    }
    let m1 = (degree + 1 - mults[0]).max(0);
    let mut prefix = Vec::new();
    let mut j = knots.len() - 2;
    let mut m = 1i32;
    for _ in 0..m1 {
        prefix.push(knots[j] - period);
        m += 1;
        if m > mults[j] && j > 0 {
            j -= 1;
            m = 1;
        }
    }
    prefix.reverse();
    prefix.extend(body);
    prefix
}

fn transform_edge(
    edge: &Edge,
    g: &GpGTrsf,
    scale: f64,
    seen_verts: &mut Vec<TopoShape>,
) -> Result<(), String> {
    let reg = GeometryRegistry::global();
    // `NewCurve` reads the edge tolerance before the later scale
    // (`BRepTools_NurbsConvertModification.cxx:409`).
    let edge_tol = reg.edge_tolerance(&edge.0);
    reg.set_edge_tolerance(&edge.0, edge_tol * scale);
    if !BRepTool::is_degenerated(edge) {
        if let Some(curve) = reg.edge_curve(&edge.0) {
            let (first, last) = BRepTool::edge_parameters(edge);
            // `NewCurve` trims a non-NURBS curve to the edge range before
            // `CurveToBSplineCurve` (`NurbsConvertModification.cxx:430-463`).
            let source: Arc<dyn Curve> =
                if curve.bspline_poles().is_some()
                    || curve.bezier_poles().is_some()
                    || curve.is_geom_trimmed()
                    || !first.is_finite()
                    || !last.is_finite()
                {
                    curve
                } else {
                    Arc::new(GeomTrimmedCurve::new(curve, first, last))
                };
            let mut moved = curve_to_bspline_curve(source.as_ref(), ParameterisationType::TgtThetaOver2)
                .map_err(|_| "CurveToBSplineCurve: BRepTools_NurbsConvertModification.cxx:463".to_string())?;
            for pole in &mut moved.poles {
                *pole = transform_xyz(g, pole);
            }
            let knot_first = moved.first_parameter();
            let knot_last = moved.last_parameter();
            // `Resolution` then `BSplCLib::Reparametrize` + `SetKnots`
            // (`BRepTools_NurbsConvertModification.cxx:469-476`).
            if !moved.periodic {
                let u_tol = moved.resolution(edge_tol);
                if (first - knot_first).abs() > u_tol || (last - knot_last).abs() > u_tol {
                    let (mut knots, _) = moved.distinct_knots_and_mults();
                    occt_core::bspl::knots::reparametrize(first, last, &mut knots);
                    moved.set_knots(&knots).map_err(|_| {
                        "Geom_BSplineCurve::SetKnots: BRepTools_NurbsConvertModification.cxx:476"
                    })?;
                }
            }
            edge.tshape.write().expect("lock").edge_core_mut().curve = Some(Arc::new(moved));
        }
    }
    let verts = edge.0.tshape.read().expect("lock").children.clone();
    for vert in verts {
        if !vert.is_vertex() || seen_verts.iter().any(|old| is_same(old, &vert)) {
            continue;
        }
        seen_verts.push(vert.clone());
        let p = reg.vertex_point(&vert);
        crate::shape::Vertex(vert.clone()).set_point(transform_xyz(g, &p));
        let tol = reg.vertex_tolerance(&vert) * scale;
        crate::shape::Vertex(vert).set_tolerance(tol);
    }
    Ok(())
}
