//! `ChFi3d_Builder::ComputeCurve2d` - the p-curve of a curve on a face.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:456-483`
//! (OCCT 8.0.0).
//!
//! ```text
//! BRepLib_MakeEdge Bedge(Ct); TopoDS_Edge edg = Bedge.Edge();
//! BRepAlgo_NormalProjection OrtProj;
//! OrtProj.Init(Face); OrtProj.Add(edg);
//! OrtProj.SetParams(1.e-6, 1.e-6, GeomAbs_C1, 14, 16);
//! OrtProj.SetLimit(false); OrtProj.Compute3d(false); OrtProj.Build();
//! if (OrtProj.IsDone()) {
//!   TopExp::MapShapes(OrtProj.Projection(), TopAbs_EDGE, MapE1);
//!   if (MapE1.Extent() != 0) {
//!     E1 = TopoDS::Edge(MapE1(1));
//!     C2d = BRep_Tool::CurveOnSurface(E1, Face, up1, up2);
//!   }
//! }
//! ```
//!
//! `BRepAlgo_NormalProjection` is [`crate::normal_projection`]. The call
//! passes `SetLimit(false)` (`cxx:469`), i.e. it takes exactly the branch of
//! `BRepAlgo_NormalProjection::Build` the port covers; the face-bounds branch
//! that the port does not cover is never reached from here.
//!
//! `BRep_Tool::CurveOnSurface(E1, Face, up1, up2)` answers the p-curve and its
//! range; [`crate::boptools_2d::curve_on_surface_range`] is that call.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::brep_lib_make_edge;
use crate::builder::TopoBuilder;
use crate::normal_projection::NormalProjection;
use crate::shape::Face;
use crate::topo_tools::edges_of;

/// `ComputeCurve2d(Ct, Face, C2d)` (`cxx:456-483`). Answers the 2D curve of
/// `ct` on `face`, or `None` where OCCT leaves `C2d` untouched
/// (`!OrtProj.IsDone()` or `MapE1.Extent() == 0`).
pub fn compute_curve2d(ct: &Arc<dyn Curve>, face: &Face) -> Option<Arc<dyn Curve2d>> {
    // `cxx:462`: `BRepLib_MakeEdge Bedge(Ct)` - the whole curve.
    let builder = TopoBuilder::new();
    let first = ct.first_parameter();
    let last = ct.last_parameter();
    let mut edg = builder.make_edge(ct.clone(), first, last);
    if let Some(build) = brep_lib_make_edge::init_3d(ct, first, last) {
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edg, &build);
    }

    // `cxx:463-470`.
    let mut proj = NormalProjection::new();
    proj.init(face.0.clone());
    proj.add(&edg.0);
    proj.set_params(1e-6, 1e-6, GeomAbsShape::C1 as u8, 14, 16);
    proj.set_limit(false);
    proj.compute3d(false);
    if proj.build().is_err() {
        return None;
    }

    // `cxx:471-478`: `TopExp::MapShapes(Projection(), TopAbs_EDGE, MapE1)`
    // and the first entry of that map.
    if !proj.is_done() {
        return None;
    }
    let projected = edges_of(proj.projection());
    let e1 = projected.first()?;

    // `cxx:480`: `BRep_Tool::CurveOnSurface(E1, Face, up1, up2)`.
    crate::boptools_2d::curve_on_surface_range(e1, face).map(|(c, _, _)| c)
}
