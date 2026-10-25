//! `TopOpeBRepTool_ShapeTool` — the orientation-comparison family and
//! `EdgeData`.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_ShapeTool.hxx:34-96`,
//! `TopOpeBRepTool_ShapeTool.cxx:394-609`, `:613-653`.
//!
//! Only the members this batch needs are translated: `ShapesSameOriented`,
//! `FacesSameOriented`, `SurfacesSameOriented`, `EdgesSameOriented`,
//! `CurvesSameOriented`, both `EdgeData` overloads, `Pnt`, `Tolerance`,
//! `AdjustOnPeriodic` and the `Resolution3d` family. The rest of the class
//! (`Baseline`, `BASISCURVE` / `BASISSURFACE`, `Bounds`, `UVRep`, ...) is
//! still UNPORTED — see the note at the end of this file.

use occt_core::gp::GpDir;
use occt_core::precision;

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

use super::adaptor::{CurveAdaptor, SurfaceAdaptor};
use super::lprop::CLProps;
use super::project::fun_tool_proj_pon_s;

/// `TopOpeBRepTool_ShapeTool::ShapesSameOriented(S1, S2)`
/// (`ShapeTool.cxx:394-440`).
pub fn shapes_same_oriented(s1: &TopoShape, s2: &TopoShape) -> bool {
    let mut so = true;

    if is_same(s1, s2) {
        let o1 = s1.orientation();
        let o2 = s2.orientation();
        if is_forward_or_reversed(o1) && is_forward_or_reversed(o2) {
            so = o1 == o2;
            return so;
        }
    }

    let t1 = s1.shape_type();
    let t2 = s2.shape_type();
    if t1 == ShapeType::Solid && t2 == ShapeType::Solid {
        so = true;
    } else if t1 == ShapeType::Face && t2 == ShapeType::Face {
        so = faces_same_oriented(s1, s2);
    } else if t1 == ShapeType::Edge && t2 == ShapeType::Edge {
        so = edges_same_oriented(s1, s2);
    } else if t1 == ShapeType::Vertex && t2 == ShapeType::Vertex {
        let o1 = s1.orientation();
        let o2 = s2.orientation();
        if is_external_or_internal(o1) || is_external_or_internal(o2) {
            so = true;
        } else {
            so = o1 == o2;
        }
    }

    so
}

/// `TopAbs_FORWARD || TopAbs_REVERSED`.
fn is_forward_or_reversed(o: Orientation) -> bool {
    matches!(o, Orientation::Forward | Orientation::Reversed)
}

/// `TopAbs_EXTERNAL || TopAbs_INTERNAL`.
fn is_external_or_internal(o: Orientation) -> bool {
    matches!(o, Orientation::External | Orientation::Internal)
}

/// `TopOpeBRepTool_ShapeTool::FacesSameOriented(S1, Sref)`
/// (`ShapeTool.cxx:522-546`).
pub fn faces_same_oriented(s1: &TopoShape, s_ref: &TopoShape) -> bool {
    let s2 = s_ref;
    let (Some(f1), Some(f2)) = (Face::wrap(s1.clone()), Face::wrap(s2.clone())) else {
        return true;
    };
    let o1 = f1.0.orientation();
    let o2 = f2.0.orientation();
    if is_external_or_internal(o1) || is_external_or_internal(o2) {
        return true;
    }

    // `const bool computerestriction = false;` (`cxx:535`)
    let bas1 = SurfaceAdaptor::of_face(&f1);
    let bas2 = SurfaceAdaptor::of_face(&f2);
    let so = is_same(&f1.0, &f2.0) || surfaces_same_oriented(&bas1, &bas2);
    let mut b = so;
    if o1 != o2 {
        b = !so;
    }
    b
}

/// `TopOpeBRepTool_ShapeTool::SurfacesSameOriented(S1, Sref)`
/// (`ShapeTool.cxx:444-518`).
///
/// The `else` arm (`cxx:504-513`) only prints a NYI trace under
/// `OCCT_DEBUG`; `so` keeps its initial `true`, as here.
pub fn surfaces_same_oriented(s1: &SurfaceAdaptor, s_ref: &SurfaceAdaptor) -> bool {
    use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;

    let s2 = s_ref;
    let st1 = s1.get_type();
    let st2 = s2.get_type();

    let mut so = true;

    if st1 == GeomAbsSurfaceType::Plane && st2 == GeomAbsSurfaceType::Plane {
        let n1 = normal_at_first_parameters(s1);
        let n2 = normal_at_first_parameters(s2);
        if let (Some(n1), Some(n2)) = (n1, n2) {
            let d = n1.dot(&n2);
            so = d > 0.0;
        }
    } else if st1 == GeomAbsSurfaceType::Cylinder && st2 == GeomAbsSurfaceType::Cylinder {
        // "On peut projeter n'importe quel point. prenons donc l'origine"
        let (Some((p1, d1u, d1v)), Some(hs2)) = (s1.d1(0.0, 0.0), s2.surface().cloned()) else {
            return so;
        };
        let n1 = d1u.crossed(&d1v);
        let Some((u2, v2)) = fun_tool_proj_pon_s(&p1, &hs2) else {
            // `if (!ok) { return so; // NYI : raise }` (`cxx:489-492`)
            return so;
        };
        let Some((_p2, d2u, d2v)) = s2.d1(u2, v2) else {
            return so;
        };
        let n2 = d2u.crossed(&d2v);
        let d = n1.dot(&n2);
        so = d > 0.0;
    }
    // `else` (`cxx:504-513`): NYI trace only, `so` unchanged.

    so
}

/// `SurfacesSameOriented`'s plane arm (`cxx:453-472`):
/// `D1(FirstUParameter, FirstVParameter)` then `d1u.Crossed(d1v)`.
fn normal_at_first_parameters(s: &SurfaceAdaptor) -> Option<occt_core::gp::GpVec> {
    let (_, d1u, d1v) = s.d1(s.first_u_parameter(), s.first_v_parameter())?;
    Some(d1u.crossed(&d1v))
}

/// `TopOpeBRepTool_ShapeTool::EdgesSameOriented(S1, Sref)`
/// (`ShapeTool.cxx:589-609`).
pub fn edges_same_oriented(s1: &TopoShape, s_ref: &TopoShape) -> bool {
    let s2 = s_ref;
    let (Some(e1), Some(e2)) = (Edge::wrap(s1.clone()), Edge::wrap(s2.clone())) else {
        return true;
    };
    let o1 = e1.0.orientation();
    let o2 = e2.0.orientation();
    if is_external_or_internal(o1) || is_external_or_internal(o2) {
        return true;
    }
    let bac1 = CurveAdaptor::of_edge(&e1);
    let bac2 = CurveAdaptor::of_edge(&e2);
    let so = curves_same_oriented(&bac1, &bac2);
    let mut b = so;
    if o1 != o2 {
        b = !so;
    }
    b
}

/// `TopOpeBRepTool_ShapeTool::CurvesSameOriented(C1, Cref)`
/// (`ShapeTool.cxx:549-587`).
///
/// The `else` arm (`cxx:569-582`) only prints a NYI trace under
/// `OCCT_DEBUG`; `so` keeps its initial `true`, as here.
pub fn curves_same_oriented(c1: &CurveAdaptor, c_ref: &CurveAdaptor) -> bool {
    use crate::geom_abs_curve_type::GeomAbsCurveType;

    let c2 = c_ref;
    let ct1 = c1.get_type();
    let ct2 = c2.get_type();
    let mut so = true;

    if ct1 == GeomAbsCurveType::Line && ct2 == GeomAbsCurveType::Line {
        let p1 = c1.first_parameter();
        let (t1, _c1) = match edge_data_tangent_and_curvature(c1, p1) {
            Some(v) => v,
            // `EdgeData` raises `LProp_NotDefined` through `BL.Tangent(T)` when
            // the tangent is undefined; OCCT does not catch it here either.
            None => return so,
        };
        let p2 = c2.first_parameter();
        let (t2, _c2) = match edge_data_tangent_and_curvature(c2, p2) {
            Some(v) => v,
            None => return so,
        };
        let d = t1.dot(&t2);
        so = d > 0.0;
    }
    // `else` (`cxx:569-582`): NYI trace only, `so` unchanged.

    so
}

/// `TopOpeBRepTool_ShapeTool::EdgeData(BAC, P, T, N, C)`
/// (`ShapeTool.cxx:613-635`), with the `N` output collapsed into the return
/// value.
///
/// The `tolm` clamp (`cxx:626-630`) takes `max(tol, max(Epsilon(0.),
/// RealLast()))`, i.e. `RealLast()`, so `|C| > tolm` is never true and `N` is
/// never written — the xpu150399 workaround as it stands in OCCT 8.0.0. The
/// tangent is still computed and can fail with `LProp_NotDefined`.
pub fn edge_data(c: &CurveAdaptor, p: f64) -> Result<EdgeDataResult, String> {
    let tol = precision::ANGULAR;
    let mut der_order = 0;
    let mut props = CLProps::new(c, p, 2, tol)?;
    let tangent = props.tangent(c, p, &mut der_order)?;
    let curvature = props.curvature_of(c, p, &mut der_order)?;

    let tol1 = f64::EPSILON; // Epsilon(0.)
    let tol2 = f64::MAX; // RealLast()
    let tolm = tol.max(tol1.max(tol2));

    let normal = if curvature.abs() > tolm {
        props.normal(c, p, &mut der_order).ok()
    } else {
        None
    };
    Ok(EdgeDataResult {
        tangent,
        normal,
        curvature,
        tolerance: tol,
    })
}

/// The out-parameters of `EdgeData` (`ShapeTool.cxx:613-625`):
/// `T` (tangent), `N` (normal), `C` (curvature) and the returned `tol`.
#[derive(Debug, Clone)]
pub struct EdgeDataResult {
    /// `gp_Dir& T`
    pub tangent: GpDir,
    /// `gp_Dir& N` — `None` when `|C| <= tolm` (`cxx:630-633`), which is
    /// always the case with OCCT 8.0.0's `tolm`.
    pub normal: Option<GpDir>,
    /// `double& C`
    pub curvature: f64,
    /// The return value `Precision::Angular()`.
    pub tolerance: f64,
}

/// The `(T, C)` pair `CurvesSameOriented` reads (`cxx:561-568`); the `N` and
/// `tol` outputs of `EdgeData` are not read there.
fn edge_data_tangent_and_curvature(c: &CurveAdaptor, p: f64) -> Option<(GpDir, f64)> {
    let r = edge_data(c, p).ok()?;
    Some((r.tangent, r.curvature))
}

/// `TopOpeBRepTool_ShapeTool::Tolerance(S)` (`ShapeTool.cxx:46-69`):
/// face / edge / vertex tolerances; every other shape type raises
/// `Standard_ProgramError`, which is `Err` here.
pub fn tolerance_of(s: &TopoShape) -> Result<f64, String> {
    match s.shape_type() {
        ShapeType::Face => Ok(crate::brep_tool::BRepTool::face_tolerance(&Face(s.clone()))),
        ShapeType::Edge => Ok(crate::brep_tool::BRepTool::edge_tolerance(&Edge(s.clone()))),
        ShapeType::Vertex => Ok(crate::brep_tool::BRepTool::vertex_tolerance(&crate::shape::Vertex(
            s.clone(),
        ))),
        _ => Err("TopOpeBRepTool_ShapeTool : Shape has no tolerance".into()),
    }
}

/// `TopOpeBRepTool_ShapeTool::Pnt(S)` (`ShapeTool.cxx:73-80`): the vertex
/// point; any other shape type raises `Standard_ProgramError`.
pub fn pnt_of(s: &TopoShape) -> Result<occt_core::gp::GpPnt, String> {
    if s.shape_type() != ShapeType::Vertex {
        return Err("TopOpeBRepTool_ShapeTool::Pnt".into());
    }
    Ok(crate::brep_tool::BRepTool::vertex_point(&crate::shape::Vertex(
        s.clone(),
    )))
}

/// `TopOpeBRepTool_ShapeTool::AdjustOnPeriodic(F, u, v)`
/// (`ShapeTool.cxx:207-253`). Maps `u` / `v` into the period of the face's
/// UV bounds when the supporting surface is periodic in that direction.
pub fn adjust_on_periodic(f: &Face, u: &mut f64, v: &mut f64) {
    let Some(surf) = crate::brep_tool::BRepTool::face_surface(f) else {
        // `BRep_Tool::Surface` returning a null handle is a null dereference in
        // OCCT; the port leaves the parameters untouched instead.
        return;
    };
    let is_uperio = surf.is_u_periodic();
    let is_vperio = surf.is_v_periodic();

    // exit if surface supporting F is not periodic on U or V
    if !is_uperio && !is_vperio {
        return;
    }

    let (uf_first, _uf_last, vf_first, _vf_last) = crate::brep_surface::face_uv_bounds(f);
    let tol = precision::PCONFUSION;

    if is_uperio {
        let u_period = surf.u_period();
        if (*u - uf_first - u_period).abs() > tol {
            *u = occt_core::elib::clib::in_period(*u, uf_first, uf_first + u_period);
        }
    }
    if is_vperio {
        let v_period = surf.v_period();
        if (*v - vf_first - v_period).abs() > tol {
            *v = occt_core::elib::clib::in_period(*v, vf_first, vf_first + v_period);
        }
    }
}

/// `TopOpeBRepTool_ShapeTool::Resolution3dU(SU, Tol2d)`
/// (`ShapeTool.cxx:655-661`): `r3dunit * (Tol2d / GAS.UResolution(r3dunit))`,
/// with the fixed `r3dunit = 0.00001` OCCT uses to avoid the `RangeError` of a
/// torus adaptor at `R3d = 1`.
pub fn resolution3d_u(su: &dyn occt_geom::surface::Surface, tol2d: f64) -> f64 {
    let r3dunit = 0.00001;
    let ru = occt_geom::approx_same_parameter::u_resolution(su, r3dunit);
    r3dunit * (tol2d / ru)
}

/// `TopOpeBRepTool_ShapeTool::Resolution3dV(SU, Tol2d)`
/// (`ShapeTool.cxx:665-673`).
pub fn resolution3d_v(su: &dyn occt_geom::surface::Surface, tol2d: f64) -> f64 {
    let r3dunit = 0.00001;
    let rv = occt_geom::approx_same_parameter::v_resolution(su, r3dunit);
    r3dunit * (tol2d / rv)
}

/// `TopOpeBRepTool_ShapeTool::Resolution3d(SU, Tol2d)`
/// (`ShapeTool.cxx:677-684`): the larger of the two directions.
pub fn resolution3d(su: &dyn occt_geom::surface::Surface, tol2d: f64) -> f64 {
    let ru = resolution3d_u(su, tol2d);
    let rv = resolution3d_v(su, tol2d);
    ru.max(rv)
}

/// `TopOpeBRepTool_ShapeTool::Resolution3d(F, Tol2d)` (`ShapeTool.cxx:688-694`):
/// `Resolution3d(BRep_Tool::Surface(F), Tol2d)`.
///
/// `None` stands for a face without geometry, where OCCT dereferences the null
/// `Geom_Surface` handle inside `GeomAdaptor_Surface`.
pub fn resolution3d_face(f: &TopoShape, tol2d: f64) -> Option<f64> {
    let su = crate::brep_tool::BRepTool::face_surface(&Face(f.clone()))?;
    Some(resolution3d(su.as_ref(), tol2d))
}

// UNPORTED: the rest of `TopOpeBRepTool_ShapeTool` (`hxx:34-96`,
// `ShapeTool.cxx:82-390`, `:696-`):
// - `BASISCURVE` / `BASISSURFACE` / `Bounds` / `UVRep` / `AnyFace` /
//   `Faces` / `ShapeType`-driven helpers (`cxx:82-390`), which
//   `TopOpeBRep` and `TopOpeBRepBuild` drive.
// They are ported with their first consumer, not speculatively.

/// `TopOpeBRepTool_ShapeTool::Closed(S1, S2)` (`ShapeTool.cxx:82-104`): `E`
/// is a closed edge of `F` when it is closed on `F` and occurs at least twice
/// among `F`'s edges.
pub fn closed(s1: &TopoShape, s2: &TopoShape) -> bool {
    let e = Edge(s1.clone());
    let f = Face(s2.clone());
    let brepclosed = crate::brep_tool::BRepTool::is_closed_edge_face(&e, &f);
    if !brepclosed {
        return false;
    }
    let mut n = 0;
    let mut ex = crate::topexp::Explorer::new(s2, ShapeType::Edge);
    while ex.more() {
        if is_same(ex.current(), s1) {
            n += 1;
        }
        ex.next();
    }
    n >= 2
}

