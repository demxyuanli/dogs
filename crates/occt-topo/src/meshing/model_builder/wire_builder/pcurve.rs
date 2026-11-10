use super::super::prelude::*;

/// `BRep_Tool::CurveOnSurface(edge, face, f, l)` (`BRep_Tool.cxx:327-361`): the
/// edge's pcurve on the face together with the **`BRep_CurveOnSurface`
/// representation's own** parameter range (`GC->Range(First, Last)`, cxx:353),
/// not the edge's 3D range. `ShapeAnalysis_Edge::PCurve` (cxx:192-208) hands
/// exactly these `f`/`l` to its callers (the 2D branch of
/// `ShapeAnalysis_Wire::CheckOrder`, cxx:627-649, reads `c2d->Value(f)` /
/// `c2d->Value(l)`). Falls back to the 3D range when the edge carries no stored
/// COS range for this face (what the port's earlier code always used).
///
/// The distinction only matters when the two differ: after a split of a
/// periodic edge the child keeps the parent's pcurve but its COS range is the
/// sliced window (`ShapeBuild_Edge::CopyRanges`, `ShapeBuild_Edge.cxx:206-334`),
/// so evaluating at the 3D range would move every 2D endpoint by one period.
/// `ShapeAnalysis_Wire::CheckOrder` (cxx:593-651) and `ShapeAnalysis::TotCross2D`
/// both read the pcurve *through* this accessor, so the port has to as well.
pub(super) fn pcurve_and_range(e: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    if let Some(v) = crate::boptools_2d::curve_on_surface_oriented(e, face, false) {
        return Some(v);
    }
    let pc = make_pcurve_full(e, face).ok()?;
    let fwd = Edge(e.0.oriented(Orientation::Forward));
    let (a, b) = BRepTool::edge_parameters(&fwd);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    Some((pc, a, b))
}

/// Uniform sample of a pcurve over `[a, b]` (endpoints + interior) with the
/// point count OCCT uses for both `ShapeAnalysis::TotCross2D`
/// (`ShapeAnalysis_Curve::GetSamplePoints`, `ShapeAnalysis_Curve.cxx:1317-1337`)
/// and the `BRepTopAdaptor_FClass2d` boundary polygon
/// (`BRepTopAdaptor_FClass2d.cxx:179-185`).
pub(crate) fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    let n = sample_count_2d(pc, a, b);
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = a + (b - a) * i as f64 / (n - 1) as f64;
        pts.push(pc.d0(t));
    }
    pts
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:73-91`) plus the `nbs *= 4` applied by the
/// two callers that share this count:
///
/// * `ShapeAnalysis_Curve.cxx:1325-1331` (`GetSamplePoints` for `TotCross2D`)
/// * `BRepTopAdaptor_FClass2d.cxx:182-184`
///
/// Both apply the `*4` at their own call site, so it stays here rather than in
/// [`crate::curve_sampling_2d::geom2d_nb_samples`].
pub(crate) fn sample_count_2d(pc: &dyn Curve2d, a: f64, b: f64) -> usize {
    let mut nbs = crate::curve_sampling_2d::geom2d_nb_samples(pc, a, b);
    if nbs > 2 {
        nbs *= 4;
    }
    // `nbPoints` never returns less than 2 (the `Geom2d_Line` arm).
    nbs.max(2)
}
