use super::prelude::*;
use super::statics::*;


/// `BRep_Builder::UpdateEdge(anE, aCurve, max_deviation)` +
/// `UpdateEdge(anE, aLine2d, theSurf, Loc, max_deviation)` +
/// `Range(anE, FirstPar, LastPar)` (`cxx:3564-3574`, `:3589-3599`).
///
/// The 3D curve comes from `GeomLib::BuildCurve3d(Confusion, ConS, FirstPar,
/// LastPar, aCurve, maxdev, avgdev)` with the header defaults `GeomAbs_C1`,
/// `MaxDegree = 14`, `MaxSegment = 30` (`GeomLib.hxx:1181-1189`).
pub(super) fn rebuild_connecting_edge(edge: &Edge, surf: &Arc<dyn Surface>, p1: GpPnt2d, p2: GpPnt2d) {
    let reg = GeometryRegistry::global();
    let Some(dir) = GpDir2d::new(p2.x() - p1.x(), p2.y() - p1.y()).ok() else {
        return;
    };
    let a_line2d: Arc<dyn Curve2d> = Arc::new(line2d_at(p1, dir));
    let first_par = 0.0;
    let last_par = p1.distance(&p2);
    let cos: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(
        a_line2d.clone(),
        surf.clone(),
        first_par,
        last_par,
    ));
    let built = geom_lib_build_curve3d(
        CONFUSION,
        cos.as_ref(),
        surf,
        &a_line2d,
        first_par,
        last_par,
        Continuity::C1,
        14,
        30,
    );
    match built.curve {
        Some(curve) => reg.set_edge_curve3d(&edge.0, curve, built.max_deviation),
        None => reg.remove_edge_curve3d(&edge.0),
    }
    reg.update_edge_pcurve_on_surface(
        &edge.0,
        surf,
        Some(a_line2d),
        built.max_deviation,
    );
    reg.set_ranges_all(&edge.0, first_par, last_par);
}
