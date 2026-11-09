//! `BRepFeat_RibSlot::EdgeExtention`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:582-638`).
//!
//! Lengthens one end of `e` by `bnd / 10` along its own curve. An elementary
//! curve (line, circle, ellipse, hyperbola, parabola) simply re-runs
//! `BRepLib_MakeEdge` on the widened parameter range; every other curve is
//! extended geometrically with `GeomLib::ExtendCurveToPoint` towards a point
//! `bnd / 1000` down the tangent line.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::elib::clib;
use occt_core::gp::{GpDir, GpLin};
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::curve::Curve;
use occt_geom::geom_lib_extend_curve::extend_curve_to_point;
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurveBasis;

use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d, init_curve_points, EdgeBuild};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Edge;
use crate::shhealing::{first_vertex, last_vertex};
use crate::topo_tools::vertex_position;

/// `BRepLib_MakeEdge::Shape()`: the edge plus the endpoint vertices `Init`
/// decided on.
fn make_edge(curve: Arc<dyn Curve>, build: &EdgeBuild) -> Option<Edge> {
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve, build.first, build.last);
    add_endpoint_vertices(&b, &mut e, build);
    Some(e)
}

/// `EdgeExtention(TopoDS_Edge& e, bnd, FirstLast)` (`cxx:584-638`).
///
/// `None` is the null edge the reference leaves behind when `BRepLib_MakeEdge`
/// is not done, and the `Standard_ConstructionError` of a null tangent
/// direction (`cxx:621` / `cxx:630`: `gp_Dir(-vct)` on a null vector).
pub fn edge_extension(e: &Edge, bnd: f64, first_last: bool) -> Option<Edge> {
    let (f, l) = BRepTool::edge_parameters(e);
    let cu = BRepTool::edge_curve_world(e)?;

    // `if (cu->DynamicType() == STANDARD_TYPE(Geom_Line) || ... Geom_Parabola)`
    // (`cxx:596-599`).
    if cu.gp_line().is_some()
        || cu.gp_circ().is_some()
        || cu.gp_ellipse().is_some()
        || cu.gp_hyperbola().is_some()
        || cu.gp_parabola().is_some()
    {
        // `BRepLib_MakeEdge Edg(cu, f - bnd/10., l)` / `(cu, f, l + bnd/10.)`
        // (`cxx:602-609`).
        let build = if first_last {
            init_3d(&cu, f - bnd / 10.0, l)?
        } else {
            init_3d(&cu, f, l + bnd / 10.0)?
        };
        return make_edge(cu, &build);
    }

    // `C = new Geom_TrimmedCurve(cu, f, l)` (`cxx:592`), in basis parameters
    // like the reference.
    let trimmed = GeomTrimmedCurveBasis::new(cu, f, l);

    // `C->D1(f, pnt, vct); ln = new Geom_Line(pnt, -vct)` and the mirrored
    // `C->D1(l, pnt, vct); ln = new Geom_Line(pnt, vct)` (`cxx:618-630`).
    let (pnt, vct) = if first_last {
        let (pnt, vct) = trimmed.d1(f);
        (pnt, vct.reversed())
    } else {
        trimmed.d1(l)
    };
    let dir = GpDir::from_vec(&vct).ok()?;
    let line = GeomLine::new(GpLin::from_pnt_dir(pnt, dir));
    // `ln->D0(bnd / 1000., Pt)`.
    let pt = clib::line_value(line.lin(), bnd / 1000.0);

    // `GeomLib::ExtendCurveToPoint(C, Pt, GeomAbs_G1, false)` before the first
    // parameter, `(C, Pt, GeomAbs_G1, true)` after the last (`cxx:624`, `:632`).
    let mut extended = curve_to_bspline_curve(&trimmed, ParameterisationType::QuasiAngular).ok()?;
    extend_curve_to_point(&mut extended, &pt, 1, !first_last).ok()?;
    let extended: Arc<dyn Curve> = Arc::new(extended);

    // `BRepLib_MakeEdge Edg(C, Pt, BRep_Tool::Pnt(TopExp::LastVertex(e, true)))`
    // and its `FirstVertex` mirror (`cxx:625`, `:633`).
    let build = if first_last {
        let last = last_vertex(e)?;
        init_curve_points(&extended, &pt, &vertex_position(&last))?
    } else {
        let first = first_vertex(e)?;
        init_curve_points(&extended, &vertex_position(&first), &pt)?
    };
    make_edge(extended, &build)
}
