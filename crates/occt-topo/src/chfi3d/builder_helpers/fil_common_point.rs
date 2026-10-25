//! `ChFi3d_TrsfTrans` / `ChFi3d_FilCommonPoint` — the transition mapping and
//! the loading of a `ChFiDS_CommonPoint` out of a blending extremity.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2174-2296`
//! (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - `occ::down_cast<BRepTopAdaptor_HVertex>(SP.Vertex())->Vertex()` becomes
//!   [`Extremity::vertex`], which answers `Result<Option<&HVertex>, String>`;
//!   the port skips the vertex when it is absent instead of dereferencing the
//!   null handle,
//! - `occ::down_cast<BRepAdaptor_Curve2d>(PR.Arc())` is a plain read of
//!   [`PointOnRst::arc`]: that member is already the `Curve2dAdaptor` the
//!   downcast would yield (see `point_on_rst.rs`),
//! - `TopExp::Vertices(E, V1, V2)` is [`crate::topo_tools::edge_vertices`],
//! - `BRep_Tool::Tolerance` / `Pnt` / `Parameter` are the
//!   [`BRepTool`] queries of the same names.

use crate::abs::Orientation;
use crate::brep_blend::extremity::Extremity;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::common_point::CommonPoint;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools;

/// `ChFi3d_TrsfTrans(T1)` (`cxx:2174-2186`): `IntSurf_In` is `FORWARD`,
/// `IntSurf_Out` is `REVERSED`, everything else `INTERNAL`.
pub fn trsf_trans(t1: IntSurfTypeTrans) -> Orientation {
    match t1 {
        IntSurfTypeTrans::In => Orientation::Forward,
        IntSurfTypeTrans::Out => Orientation::Reversed,
        _ => Orientation::Internal,
    }
}

/// `ChFi3d_FilCommonPoint(SP, TransLine, Start, CP, Tol)` (`cxx:2193-2296`):
/// fills `cp` from the blending extremity, snapping onto a pre-existing vertex
/// when the point sits on one.
pub fn fil_common_point(
    sp: &Extremity,
    trans_line: IntSurfTypeTrans,
    start: bool,
    cp: &mut CommonPoint,
    tol: f64,
) {
    let mut maxtol = tol.max(cp.tolerance());

    cp.set_point(sp.value());
    if sp.has_tangent() {
        if let Ok(tangent) = sp.tangent() {
            if start {
                cp.set_vector(&tangent.reversed());
            } else {
                cp.set_vector(tangent);
            }
        }
    }

    // and the parameter of the spine.
    cp.set_parameter(sp.parameter_on_guide());

    if sp.is_vertex() {
        // The Vertex is loaded if required (inside of a face).
        if let Ok(Some(hvertex)) = sp.vertex() {
            let v = hvertex.vertex();
            cp.set_vertex(&v.0);
            let dist = sp.value().distance(&BRepTool::vertex_point(v));
            // `cxx:2223-2228` (OCC571): `maxtol += Dist` was replaced by
            // `maxtol = max(Dist, maxtol)`.
            maxtol = dist.max(maxtol);
            cp.set_point(&BRepTool::vertex_point(v));
        }
    }

    if sp.nb_point_on_rst() != 0 {
        // An arc, and/or a vertex is loaded.
        if let Some(pr) = sp.point_on_rst(1) {
            if let Some(harc) = pr.arc() {
                let e: Edge = harc.edge().clone();
                let (v0, v1) = topo_tools::edge_vertices(&e);
                let pnt = |v: &Option<crate::shape::Vertex>| match v {
                    Some(v) => BRepTool::vertex_point(v),
                    None => *sp.value(),
                };
                let distf = sp.value().distance(&pnt(&v0));
                let distl = sp.value().distance(&pnt(&v1));
                let (index_min, dist) = if distf < distl {
                    (0usize, distf)
                } else {
                    (1usize, distl)
                };
                let vmin = if index_min == 0 { &v0 } else { &v1 };

                let le_param_amoi;
                let snap = match vmin {
                    Some(v) => dist <= maxtol + BRepTool::vertex_tolerance(v),
                    None => false,
                };
                if let Some(v) = vmin {
                    if snap {
                        // A pre-existing vertex has been met.
                        cp.set_vertex(&v.0);
                        cp.set_point(&BRepTool::vertex_point(v));
                        maxtol = BRepTool::vertex_tolerance(v).max(maxtol);
                        // `cxx:2267-2271` (OCC571).
                        maxtol = dist.max(maxtol);
                        le_param_amoi = BRepTool::parameter_of_vertex(&e.0, &v.0);
                    } else {
                        // Creation of an arc only.
                        maxtol = BRepTool::edge_tolerance(&e).max(maxtol);
                        maxtol = sp.tolerance().max(maxtol);
                        le_param_amoi = pr.parameter_on_arc();
                    }
                } else {
                    maxtol = BRepTool::edge_tolerance(&e).max(maxtol);
                    maxtol = sp.tolerance().max(maxtol);
                    le_param_amoi = pr.parameter_on_arc();
                }

                // Definition of the arc.
                let or = e.orientation();
                let tr = if start {
                    Orientation::compose(trsf_trans(trans_line), or).reversed()
                } else {
                    Orientation::compose(trsf_trans(trans_line), or)
                };
                cp.set_arc(maxtol, &e.0, le_param_amoi, tr);
            }
        }
    }
    // Finally, the tolerance.
    cp.set_tolerance(maxtol);
}

/// The `TopoDS_Vertex` view the caller needs when only the shape is at hand.
pub fn vertex_shape(v: &crate::shape::Vertex) -> &TopoShape {
    &v.0
}
