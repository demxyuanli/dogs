//! `HLRBRep_EdgeFaceTool` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_EdgeFaceTool.hxx:25-49`,
//! `.cxx:24-104`): the two edge/face queries `HLRBRep_Data::Update` needs.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `CurvatureValue(F, U, V, Tg)` | `cxx:26-57` | [`curvature_value`] |
//! | `UVPoint(Par, E, F, U, V)` | `cxx:61-104` | [`uv_point`] |
//!
//! Translation notes:
//! * `BRep_Tool::CurveOnSurface(E, F, f, l).IsNull()` (`cxx:73-78`) is
//!   [`Curve2dAdaptor::curve`] answering `None`: the port's `BRepAdaptor_Curve2d`
//!   resolves the same pcurve (`brep_adaptor_curve2d.rs:62-96`).
//! * `BRepLib_MakeVertex(P)` (`cxx:81`) is [`TopoBuilder::make_vertex`] with
//!   `Precision::Confusion()`, the tolerance `BRepLib_MakeVertex` gives the
//!   vertex (`BRepLib_MakeVertex.cxx:30-34`).
//! * `BRepExtrema_ExtPF(vertex, face)` is [`ExtPF::new`] with the
//!   `Extrema_ExtFlag_MIN` / `Extrema_ExtAlgo_Grad` defaults of the
//!   two-argument constructor (`BRepExtrema_ExtPF.hxx:31-34`,
//!   `.cxx:32-37`).
//! * `CurvatureValue` answers the double directly; the two `return 0.` arms
//!   (`cxx:50`, `cxx:56`) are the single `0.0` return.

use occt_core::gp::{GpDir, GpPnt};
use occt_core::precision::{Precision, REAL_SMALL};
use occt_geom::extrema_surf::ExtremaExtAlgo;
use occt_geom::extrema_surf::ExtremaExtFlag;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_extrema_ext_pf::ExtPF;
use crate::builder::TopoBuilder;
use crate::hlrbrep_curve::HlrBRepCurve;
use crate::hlrbrep_surface::HlrBRepSurface;
use crate::shape::{Edge, Face};

/// `HLRBRep_EdgeFaceTool::CurvatureValue(F, U, V, Tg)` (`cxx:26-57`): the
/// normal curvature of `F` at `(U, V)` along the direction `Tg`.
pub fn curvature_value(f: &HlrBRepSurface, u: f64, v: f64, tg: &GpDir) -> f64 {
    // `D2(U, V, P, D1U, D1V, D2U, D2V, D2UV)` (`cxx:30`). OCCT dereferences the
    // unported projection without a test; the port's `d2` answers `None` for a
    // surface that has no analytic D2, which is the `det <= Resolution` arm's
    // zero (`cxx:56`).
    let Some((_p, d1u, d1v, d2u, d2v, d2uv)) = f.d2(u, v) else {
        return 0.0;
    };

    let t = tg.xyz();
    let d1ut = d1u.coord.dot(t);
    let d1vt = d1v.coord.dot(t);
    let d1ud1v = d1u.coord.dot(&d1v.coord);
    let nmu2 = d1u.coord.dot(&d1u.coord);
    let nmv2 = d1v.coord.dot(&d1v.coord);
    let det = nmu2 * nmv2 - d1ud1v * d1ud1v;

    // `if (det > gp::Resolution())` (`cxx:38`).
    if det > REAL_SMALL {
        let alfa = (d1ut * nmv2 - d1vt * d1ud1v) / det;
        let beta = (d1vt * nmu2 - d1ut * d1ud1v) / det;
        let alfa2 = alfa * alfa;
        let beta2 = beta * beta;
        let alfabeta = alfa * beta;

        // `gp_Vec Nm = D1U ^ D1V; Nm.Normalize();` (`cxx:44-45`).
        let nm = d1u.coord.crossed(&d1v.coord).normalized();

        let n = nm.dot(&d2u.coord) * alfa2
            + 2.0 * nm.dot(&d2uv.coord) * alfabeta
            + nm.dot(&d2v.coord) * beta2;
        let d = nmu2 * alfa2 + 2.0 * d1ud1v * alfabeta + nmv2 * beta2;
        return n / d;
    }
    0.0 // `cxx:56`
}

/// `HLRBRep_EdgeFaceTool::UVPoint(Par, E, F, U, V)` (`cxx:61-104`): the surface
/// parameters of the point of edge `e` at parameter `par`, or `None` where OCCT
/// answers `false` (`cxx:97-100`: no extrema solution).
pub fn uv_point(par: f64, e: &HlrBRepCurve, f: &HlrBRepSurface) -> Option<(f64, f64)> {
    // `Curve().Edge()` and `Surface().Face()` (`cxx:73-74`).
    let edge: Edge = e.edge()?.clone();
    let face: Face = f.surface_adaptor().face()?.clone();

    // `BRep_Tool::CurveOnSurface(E, F, f, l).IsNull()` (`cxx:73-78`).
    let pc = Curve2dAdaptor::of_edge_face(&edge, &face);
    if pc.curve().is_none() {
        // `BRepExtrema_ExtPF proj(BRepLib_MakeVertex(Value3D(Par)), Face)`
        // (`cxx:81-90`).
        let p: GpPnt = e.value3d(par);
        let vertex = TopoBuilder.make_vertex(p, Precision::CONFUSION);
        let proj = ExtPF::new(
            &vertex,
            &face,
            ExtremaExtFlag::Min,
            ExtremaExtAlgo::Grad,
        );
        let n = proj.nb_ext();
        let mut index = 0usize;
        let mut dist2 = f64::MAX; // `RealLast()` (`cxx:87`)
        for i in 1..=n {
            let newdist2 = proj.square_distance(i);
            if newdist2 < dist2 {
                dist2 = newdist2;
                index = i;
            }
        }
        if index == 0 {
            return None; // `cxx:97-100`
        }
        // `proj.Parameter(index, U, V)` (`cxx:102`).
        return Some(proj.parameter(index));
    }

    // `BRepAdaptor_Curve2d PC(E, F); gp_Pnt2d P2d; PC.D0(Par, P2d);`
    // (`cxx:99-102`).
    let p2d = pc.d0(par)?;
    Some((p2d.x(), p2d.y()))
}
