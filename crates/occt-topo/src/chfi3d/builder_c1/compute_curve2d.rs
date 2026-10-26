//! `ComputeCurve2d` - UNPORTED.
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
//! The body is a single dependency: `BRepAlgo_NormalProjection`
//! (`ModelingAlgorithms/TKTopAlgo/BRepAlgo/BRepAlgo_NormalProjection.cxx`),
//! which is not ported - it builds a projected wire through `BRepAlgo_Approx`
//! and `GeomProjLib`'s general-surface projector. Neither is ported, so this
//! entry answers `None` (OCCT's "no edge was projected": `MapE1.Extent() == 0`
//! or `!OrtProj.IsDone()`), which is the same shape of answer the OCCT body
//! gives on failure. The two callers that need a real projection
//! (`PerformOneCorner`'s `cxx:1338` / `:1342`, `PerformIntersectionAtEnd`'s
//! `cxx:3434` / `:3504`) are UNPORTED together with the rest of their bodies.

use std::sync::Arc;

use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::shape::Face;

/// `ComputeCurve2d(Ct, Face, C2d)` (`cxx:456-483`) - UNPORTED, see the module
/// note. Answers the 2D curve of `ct` on `face`, or `None`.
pub fn compute_curve2d(_ct: &Arc<dyn Curve>, _face: &Face) -> Option<Arc<dyn Curve2d>> {
    None
}
