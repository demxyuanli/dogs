use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, TrsfForm};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::builder::TopoBuilder;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, TopoShape};
use crate::shhealing::ShapeBuildEdge;
use crate::tgeometry::GeometryRegistry;

use super::super::iges_to_brep_basic_curve::{basic_curve_2d, Curve2dSpan};
use super::super::iges_to_brep_curve_and_surface::Session;

/// An edge carrying only a pcurve on `face`, the shape
/// `ShapeBuild_Edge::MakeEdge(edge, C2d, face, first, last)`
/// (`ShapeBuild_Edge.cxx:666-702`) creates. The port's `EdgeGeom::curve` is
/// non-nullable, so the pcurve is wrapped in a `CurveOnSurface` 3D
/// representation, and `EdgeGeomCore::curve3d_null` records that this wrapper
/// is only `BRepAdaptor_Curve(E, F)`'s evaluator, not a `BRep_Curve3D` curve.
///
/// The representation the wrapper stands in for is OCCT's edge that holds a
/// `BRep_Curve3D` whose curve is null plus a `BRep_CurveOnSurface`. The null
/// `BRep_Curve3D` is what `BRep_Builder::Degenerated(E, true)` installs
/// (`BRep_Builder.cxx:1073-1085` `UpdateCurves(curves, nullCurve, loc)`), and
/// its range is the `BRep_GCurve` default `[-RealLast(), RealLast()]`
/// (`BRep_GCurve` ctor). The two are equivalent for every range query the port
/// makes: `BRep_Tool::Range` skips that null-curve representation and falls
/// through to the curve-on-surface one (`BRep_Tool.cxx:931-960`), and
/// `BRep_Tool::Curve` returns null with `[-RealLast(), RealLast()]`
/// (`BRep_Tool.cxx:172-196`), which is `curve3d_null` here. Verified
/// against the oracle: on the IGES sphere face 0 the port's two pcurve-only
/// edges carry the same `pcrange = [0, 2*PI]` and the same `BRep_Tool::Range`
/// value OCCT reports (`occt_probe --igespc`, `erange`), and on the cone face 0
/// edge 0 the oracle reports `null3d=0 c3=Geom_BSplineCurve` because
/// `ShapeFix_Wire::FixEdgeCurves`'s `FixAddCurve3d` arm
/// (`ShapeFix_Wire.cxx:881-918`) rebuilt it from that pcurve.
pub(super) fn pcurve_edge(
    face: &TopoShape,
    pc: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
) -> Edge {
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .expect("IGES face has no surface");
    // `ShapeBuild_Edge().MakeEdge(newedge, curve->Reversed(), face, ...)`
    // (`IGESToBRep_IGESBoundary.cxx:408-421`) -> `BRepBuilderAPI_MakeEdge(C2d,
    // S, p1, p2)` -> `BRepLib_MakeEdge::Init(C2d, S, null, null, p1, p2)`
    // (`BRepLib_MakeEdge.cxx:905-1094`). `ShapeBuild_Edge::MakeEdge` only
    // prints in `OCCT_DEBUG` and leaves `newedge` null when `!ME.IsDone()`
    // (`ShapeBuild_Edge.cxx:859-865`), and `ReverseCurves2d` then stores that
    // null edge; the caller here keeps the un-reversed pcurve range instead.
    let Some(build) = crate::brep_lib_make_edge::init_2d(&pc, &surf, first, last) else {
        let mut plain = TopoBuilder::new().make_edge(
            Arc::new(CurveOnSurface::new(pc.clone(), surf, first, last)),
            first,
            last,
        );
        GeometryRegistry::global().set_curve3d_null(&plain.0, true);
        attach_pcurve(&plain, face, pc, first, last);
        return plain;
    };
    let curve3d: Arc<dyn occt_geom::Curve> =
        Arc::new(CurveOnSurface::new(pc.clone(), surf, build.first, build.last));
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve3d, build.first, build.last);
    // `BRep_Tool::Curve` answers null on this edge: `ShapeBuild_Edge::MakeEdge`
    // leaves a `BRep_Curve3D` holding a null curve (`BRep_Builder::UpdateCurves`,
    // `BRep_Builder.cxx:86-92`). `curve` above is only what
    // `BRepAdaptor_Curve(E, F)` builds from the pcurve.
    GeometryRegistry::global().set_curve3d_null(&edge.0, true);
    crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    if build.reverse {
        edge.0.set_orientation(crate::abs::Orientation::Reversed);
    }
    attach_pcurve(&edge, face, pc, build.first, build.last);
    edge
}

/// `IGESToBRep_IGESBoundary::ReverseCurves2d` (`cxx:394-426`).
///
/// `ShapeExtend_WireData::Reverse(face)` (`ShapeExtend_WireData.cxx:545-572`)
/// reverses the edge order and every edge orientation before the pcurves are
/// reversed (`Reverse`, `ShapeExtend_WireData.cxx:483-506`), then runs
/// `ComputeSeams(true)` (`:163-222`) and `SwapSeam` on the seam occurrences
/// (`:511-543`). That seam half is not carried here. `SwapSeam` exchanges the
/// two PCurve representations an edge holds for the face (`:530-542`), while
/// every edge of this list is rebuilt by `pcurve_edge`, which gives it its own
/// TShape and a single pcurve registration for the face, so `ComputeSeams` (it
/// pairs a FORWARD and a REVERSED occurrence of the SAME TShape,
/// `ShapeExtend_WireData.cxx:190-215`) has nothing to bind. The port's seam
/// bookkeeping lives in `shape_fix_compose_shell::wire_data::swap_seam`.
///
/// Reachability: the `toreverse2d` arm is only fed by the 508 Loop path
/// (`IGESToBRep_BRepEntity.cxx:394-427`), and none of the fourteen
/// `_igesdump/oracle` models carries an entity 502, 504, 508, 510, 514 or 186.
pub(super) fn reverse_curves_2d(edges: &[Edge], face: &TopoShape) -> Vec<Edge> {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let mut built = Vec::with_capacity(edges.len());
    for edge in edges.iter().rev() {
        let orientation = edge.0.orientation().reversed();
        let Some(pc) = reg.edge_pcurves(&edge.0, key).into_iter().next() else {
            let mut plain = edge.clone();
            plain.0.set_orientation(orientation);
            built.push(plain);
            continue;
        };
        let (p1, p2) = reg
            .pcurve_range(&edge.0, key)
            .unwrap_or_else(|| reg.edge_parameters(&edge.0));
        let (first, last) = if pc.is_periodic() {
            (pc.reversed_parameter(p2), pc.reversed_parameter(p1))
        } else {
            (
                pc.first_parameter().max(pc.reversed_parameter(p2)),
                pc.last_parameter().min(pc.reversed_parameter(p1)),
            )
        };
        let rev: Arc<dyn Curve2d> = Arc::from(pc.reversed());
        let mut neu = pcurve_edge(face, rev, first, last);
        neu.0.set_orientation(orientation);
        built.push(neu);
    }
    built
}

pub(super) fn mapped_pcurve(
    session: &Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<(Arc<dyn Curve2d>, f64, f64), String> {
    // `Transfer2dTopoBasicCurve` (`cxx:1217-1230`).
    let span = basic_curve_2d(session.model, index)?;
    let (curve, mut first, mut last) = apply_uv_trsf(span, trans);
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .ok_or_else(|| "IGES face has no surface".to_string())?;
    let mut ntrsf = GpTrsf2d::identity();
    if surf.gp_pln().is_some() {
        ntrsf
            .set_scale(&GpPnt2d::new(0.0, 0.0), session.model.unit_factor)
            .map_err(|e| e.to_string())?;
    }
    let pc = ShapeBuildEdge.transform_pcurve(&curve, &ntrsf, u_fact, &mut first, &mut last);
    Ok((pc, first, last))
}

pub(super) fn apply_uv_trsf(span: Curve2dSpan, trans: &GpTrsf2d) -> (Arc<dyn Curve2d>, f64, f64) {
    if trans.form() == TrsfForm::Identity {
        return (span.curve, span.first, span.last);
    }
    let nf = span.curve.transformed_parameter(span.first, trans);
    let nl = span.curve.transformed_parameter(span.last, trans);
    let mut cloned = span.curve.clone_dyn();
    cloned.transform(trans);
    let basis = match cloned.trimmed_basis() {
        Some(b) => Arc::from(b.clone_dyn()),
        None => Arc::from(cloned),
    };
    let (first, last) = if nf <= nl { (nf, nl) } else { (nl, nf) };
    (
        Arc::new(Geom2dTrimmedCurve::new(basis, first, last)),
        first,
        last,
    )
}

pub(super) fn attach_pcurve(
    edge: &Edge,
    face: &TopoShape,
    pc: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
) {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    reg.set_edge_pcurve(&edge.0, key, pc);
    reg.set_pcurve_range(&edge.0, key, first, last);
    let (f3, l3) = reg.edge_parameters(&edge.0);
    if (f3 - first).abs() > occt_core::precision::PCONFUSION
        || (l3 - last).abs() > occt_core::precision::PCONFUSION
    {
        reg.set_same_parameter(&edge.0, false);
        reg.set_same_range(&edge.0, false);
    }
}
