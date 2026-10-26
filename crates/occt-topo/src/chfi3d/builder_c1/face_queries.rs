//! The face / edge / vertex probes of `ChFi3d_Builder_C1.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:1654-1825`
//! (OCCT 8.0.0): `cherche_face` (`:1660-1694`), `cherche_edge1` (`:1701-1730`),
//! `containV` (`:1737-1758`), `containE` (`:1764-1785`), `IsShrink`
//! (`:1792-1825`).
//!
//! The three `TopExp::MapShapes`-based probes (`cherche_face`, `cherche_edge1`,
//! `containE` / `containV`) walk the sub-shapes of a face and compare with
//! `TopoDS_Shape::IsSame`. The port uses
//! [`crate::topo_tools::map_shapes`] (the crate's `TopExp::MapShapes`) and
//! [`crate::topo_tools::is_same`] (`TShape` + `Location`, the
//! `TopTools_ShapeMapHasher` key).
//!
//! `IsShrink` samples a B-spline / Bezier pcurve with `math_FunctionSample(Pf,
//! Pl, 10)`; `GetParameter(i)` is `((n - i) * Pf + (i - 1) * Pl) / (n - 1)`
//! (`math_FunctionSample.cxx:44-48`), i.e. eleven points from `Pf` to `Pl`
//! inclusive - reproduced literally rather than re-spaced.

use occt_core::gp::GpPnt2d;
use occt_core::kernel::geomabs::CurveType;
use occt_geom2d::adaptor_curve::type2d;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_same, map_shapes};

/// The sub-shapes of `shape` of one type, in the crate's `MapShapes` order.
fn sub_shapes(shape: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    map_shapes(shape, &[t])
        .remove(&t)
        .unwrap_or_default()
}

/// `cherche_face(map, E, F1, F2, F3, F)` (`cxx:1660-1694`): the first face of
/// `map` that is none of `f1` / `f2` / `f3` and contains `e`.
///
/// OCCT raises `Standard_ConstructionError("Failed to find face.")` when `F` is
/// still null; that is the `Err` here.
pub fn cherche_face(
    map: &[TopoShape],
    e: &Edge,
    f1: &Face,
    f2: &Face,
    f3: &Face,
) -> Result<Face, String> {
    for item in map {
        let Some(fcur) = Face::wrap(item.clone()) else {
            continue;
        };
        if is_same(&fcur.0, &f1.0) || is_same(&fcur.0, &f2.0) || is_same(&fcur.0, &f3.0) {
            continue;
        }
        for e_cur in sub_shapes(&fcur.0, ShapeType::Edge) {
            if is_same(&e_cur, &e.0) {
                return Ok(fcur);
            }
        }
    }
    Err("Failed to find face.".to_string())
}

/// `cherche_edge1(F1, F2, Edge)` (`cxx:1701-1730`): the first edge that belongs
/// to both faces.
///
/// OCCT raises `Standard_ConstructionError("Failed to find edge.")` when `Edge`
/// is null; that is the `Err` here.
pub fn cherche_edge1(f1: &Face, f2: &Face) -> Result<Edge, String> {
    let edges1 = sub_shapes(&f1.0, ShapeType::Edge);
    let edges2 = sub_shapes(&f2.0, ShapeType::Edge);
    for e1 in &edges1 {
        for e2 in &edges2 {
            if is_same(e2, e1) {
                if let Some(edge) = Edge::wrap(e1.clone()) {
                    return Ok(edge);
                }
            }
        }
    }
    Err("Failed to find edge.".to_string())
}

/// `containV(F1, V)` (`cxx:1737-1758`): whether `v` is one of the vertices of
/// `f1`.
pub fn contain_v(f1: &Face, v: &TopoShape) -> bool {
    sub_shapes(&f1.0, ShapeType::Vertex)
        .iter()
        .any(|v_cur| is_same(v_cur, v))
}

/// `containE(F1, E)` (`cxx:1764-1785`): whether `e` is one of the edges of
/// `f1`.
pub fn contain_e(f1: &Face, e: &Edge) -> bool {
    sub_shapes(&f1.0, ShapeType::Edge)
        .iter()
        .any(|e_cur| is_same(e_cur, &e.0))
}

/// `gp_Pnt2d::Coord(isU ? 1 : 2)`.
fn coord2d(p: &GpPnt2d, is_u: bool) -> f64 {
    if is_u {
        p.x()
    } else {
        p.y()
    }
}

/// `IsShrink(PC, Pf, Pl, Param, isU, tol)` (`cxx:1792-1825`): whether the
/// `isU`-th coordinate of every sampled point of the pcurve between `Pf` and
/// `Pl` is within `tol` of `Param`.
///
/// The `GeomAbs_Line` arm tests the two ends only; the Bezier / B-spline arm
/// tests the eleven `math_FunctionSample(Pf, Pl, 10)` parameters; every other
/// curve type answers `false` (the source's `default:;`).
pub fn is_shrink(pc: &dyn Curve2d, pf: f64, pl: f64, param: f64, is_u: bool, tol: f64) -> bool {
    match type2d(pc) {
        // `case GeomAbs_Line:` (`cxx:1801-1807`).
        CurveType::Line => {
            let p1 = pc.d0(pf);
            let p2 = pc.d0(pl);
            (coord2d(&p1, is_u) - param).abs() <= tol
                && (coord2d(&p2, is_u) - param).abs() <= tol
        }
        // `case GeomAbs_BezierCurve: case GeomAbs_BSplineCurve:` (`cxx:1808-1821`).
        CurveType::BezierCurve | CurveType::BSplineCurve => {
            // `math_FunctionSample aSample(Pf, Pl, 10);`
            const N: usize = 10;
            for i in 1..=N {
                let u = ((N - i) as f64 * pf + (i - 1) as f64 * pl) / (N - 1) as f64;
                let p = pc.d0(u);
                if (coord2d(&p, is_u) - param).abs() > tol {
                    return false;
                }
            }
            true
        }
        // `default:;` (`cxx:1821-1822`).
        _ => false,
    }
}

