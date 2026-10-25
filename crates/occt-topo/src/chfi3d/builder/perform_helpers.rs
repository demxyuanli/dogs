//! File statics of `ChFi3d_Builder_1.cxx` that only `PerformExtremity` /
//! `PerformElement` use.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_1.cxx`
//! (OCCT 8.0.0):
//!
//! - `ReorderFaces` (`:65-123`),
//! - `ConcatCurves` (`:125-147`),
//! - `MakeOffsetEdge` (`:148-304`),
//! - `TangentExtremity` (`:625-711`),
//! - `TangentOnVertex` (`:713-728`).
//!
//! Divergences, each also stated at the call site:
//! - `ReorderFaces` loops in OCCT (`for (;;)`) until it meets one of the two
//!   faces; the port returns when a pass finds no candidate, which is the point
//!   where OCCT would spin with a stale `CurEdge` / `CurFace`,
//! - `ConcatCurves` appends a null `Geom_Curve` and loops forever when every
//!   `Add` fails; the port stops instead,
//! - `Geom_OffsetSurface::Surface()` canonicalisation is UNPORTED, so
//!   `MakeOffsetEdge` always takes OCCT's null-handle fallback (see
//!   [`make_offset_edge`]).

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::geom_lprop::SlProps;
use occt_geom::GeomOffsetSurface;
use occt_geom::GeomRectangularTrimmedSurface;
use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::make_2d;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi3d::builder_helpers::conex_faces;
use crate::chfi_ds::map::ChFiDSMap;
use crate::geom_int::IntSS;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_same, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

/// `ReorderFaces(F1, F2, FirstFace, PrevEdge, CommonVertex, EFMap)`
/// (`cxx:65-123`): makes `F1` the face reached from `FirstFace` through
/// `CommonVertex`, walking around the vertex while both candidates are neither
/// of the two.
pub(super) fn reorder_faces(
    f1: &mut Face,
    f2: &mut Face,
    first_face: &Face,
    prev_edge: &Edge,
    common_vertex: &TopoShape,
    ef_map: &ChFiDSMap,
) {
    if is_same(&f1.0, &first_face.0) {
        return;
    } else if is_same(&f2.0, &first_face.0) {
        std::mem::swap(f1, f2);
        return;
    }

    // Loop until find <theF1> or <theF2>
    let mut prev_edge = prev_edge.clone();
    let mut prev_face = first_face.clone();
    loop {
        let map = map_shapes_and_ancestors(&prev_face.0, ShapeType::Vertex, ShapeType::Edge);
        let mut cur_edge: Option<Edge> = None;
        let mut cur_face: Option<Face> = None;
        if let Some(edges) = map.get(&shape_map_key(common_vertex)) {
            for value in edges {
                let Some(edge) = Edge::wrap(value.clone()) else {
                    continue;
                };
                if is_same(&edge.0, &prev_edge.0) {
                    continue;
                }
                let flist = ef_map.find_from_key(&edge.0);
                // `Flist.First()` / `Flist.Last()` on an empty list raise
                // `Standard_NoSuchObject`; the EF map holds the edge's faces,
                // so an empty list cannot come from `ChFiDS_Map::Fill`.
                let Some(first) = flist.first().and_then(|f| Face::wrap(f.clone())) else {
                    continue;
                };
                let Some(last) = flist.last().and_then(|f| Face::wrap(f.clone())) else {
                    continue;
                };
                let face = if is_same(&prev_face.0, &first.0) { last } else { first };
                if is_same(&face.0, &f1.0) {
                    return;
                } else if is_same(&face.0, &f2.0) {
                    std::mem::swap(f1, f2);
                    return;
                }
                cur_edge = Some(edge);
                cur_face = Some(face);
            }
        }
        // OCCT's `for (;;)` reuses the last `CurEdge` / `CurFace` and would walk
        // the same ring again (or null handles) when a pass finds no candidate;
        // the port stops there.
        let (Some(edge), Some(face)) = (cur_edge, cur_face) else {
            return;
        };
        prev_edge = edge;
        prev_face = face;
    }
}

/// `ConcatCurves(theCurves, theNewCurves)` (`cxx:125-147`): glues as many
/// curves as `GeomConvert_CompCurveToBSplineCurve` can absorb into one B-spline
/// per pass.
///
/// `Err` carries the `GeomConvert` exception OCCT would let escape.
pub(super) fn concat_curves(
    curves: &mut Vec<Arc<dyn Curve>>,
    new_curves: &mut Vec<Arc<dyn Curve>>,
) -> Result<(), String> {
    while !curves.is_empty() {
        let mut concat = CompCurveToBSplineCurve::new(ParameterisationType::TgtThetaOver2);
        let mut i = 0;
        while i < curves.len() {
            let curve = curves[i].clone();
            // `Concat.Add(aBoundedCurve, 1.e-5, true)`, then the `false`
            // (prepend) retry; the `WithRatio` / `MinM` defaults of `Add` are
            // `true` / `0` (`GeomConvert_CompCurveToBSplineCurve.hxx:56-60`).
            let mut success = concat
                .add(curve.as_ref(), 1.0e-5, true, true, 0)
                .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve::Add({e:?})"))?;
            if !success {
                success = concat
                    .add(curve.as_ref(), 1.0e-5, false, true, 0)
                    .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve::Add({e:?})"))?;
            }
            if success {
                curves.remove(i);
            } else {
                i += 1;
            }
        }
        // OCCT appends `Concat.BSplineCurve()` even when it is null, and its
        // outer `while` then never terminates; the port stops instead.
        let Some(bspline) = concat.bspline_curve() else {
            return Ok(());
        };
        new_curves.push(Arc::new(bspline.clone()));
    }
    Ok(())
}

/// `MakeOffsetEdge(theEdge, Distance, S1, S2)` (`cxx:148-304`): the edge
/// obtained by intersecting the surface offset from `S1` with `S2`, trimmed to
/// the ends of `theEdge`.
///
/// `None` is OCCT's null `TopoDS_Edge` (no intersection line, a failed trim, or
/// an out-of-range concatenation).
///
/// UNPORTED — `Geom_OffsetSurface::Surface()` (`Geom_OffsetSurface.cxx:867+`)
/// returns an *equivalent canonical* surface (translated plane, wider cylinder,
/// ...) when one exists, so OCCT usually feeds an exact analytic surface to
/// `GeomInt_IntSS`. That canonicalisation is not ported, so this function takes
/// OCCT's own `if (OffsetTrGS1.IsNull()) OffsetTrGS1 = MakeOffsetSurf;` arm
/// (`cxx:169-174`) unconditionally and intersects the general offset surface.
pub(super) fn make_offset_edge(
    edge: &Edge,
    distance: f64,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
) -> Option<Edge> {
    let f1 = s1.face()?.clone();
    let f2 = s2.face()?.clone();
    let gs1 = BRepTool::face_surface(&f1)?;
    let tr_gs1 = GeomRectangularTrimmedSurface::uv(
        gs1,
        s1.first_u_parameter(),
        s1.last_u_parameter(),
        s1.first_v_parameter(),
        s1.last_v_parameter(),
    )
    .ok()?;
    let mut offset = -distance;
    if f1.0.orientation() == Orientation::Reversed {
        offset = distance;
    }
    let make_offset_surf = GeomOffsetSurface::new(Arc::new(tr_gs1), offset);
    let offset_tr_gs1: Arc<dyn Surface> = Arc::new(make_offset_surf);
    let gs2 = BRepTool::face_surface(&f2)?;
    let tr_gs2 = GeomRectangularTrimmedSurface::uv(
        gs2,
        s2.first_u_parameter(),
        s2.last_u_parameter(),
        s2.first_v_parameter(),
        s2.last_v_parameter(),
    )
    .ok()?;
    let mut intersector = IntSS::new();
    // `GeomInt_IntSS(OffsetTrGS1, TrGS2, Precision::Confusion())` with the
    // header defaults `Approx = true`, `ApproxS1 = ApproxS2 = false`
    // (`GeomInt_IntSS.hxx`).
    intersector.perform(
        offset_tr_gs1,
        Arc::new(tr_gs2),
        CONFUSION,
        true,
        false,
        false,
    );
    if !intersector.is_done() || intersector.nb_lines() == 0 {
        return None;
    }

    let mut int_curve: Arc<dyn Curve> = intersector.line(1)?.curve.clone();
    let a_ba_curve = CurveAdaptor::of_edge(edge);
    let ends: [GpPnt; 2] = [
        a_ba_curve.d0(a_ba_curve.first_parameter())?,
        a_ba_curve.d0(a_ba_curve.last_parameter())?,
    ];

    if intersector.nb_lines() > 1 {
        let mut curves: Vec<Arc<dyn Curve>> = (1..=intersector.nb_lines())
            .filter_map(|i| intersector.line(i).map(|l| l.curve.clone()))
            .collect();
        let mut new_curves: Vec<Arc<dyn Curve>> = Vec::new();
        concat_curves(&mut curves, &mut new_curves).ok()?;

        let mut min_dist = f64::MAX;
        let mut imin = 0usize;
        for (i, curve) in new_curves.iter().enumerate() {
            // `Extrema_ExtPC Projector(Ends[0], GAcurve)`.
            let projector = point_curve_extrema_all(curve.as_ref(), &ends[0]);
            if projector.is_empty() {
                continue;
            }
            for solution in &projector {
                let a_dist = solution.distance * solution.distance;
                if a_dist < min_dist {
                    min_dist = a_dist;
                    imin = i;
                }
            }
        }
        int_curve = new_curves.get(imin)?.clone();
    }

    // Projection of extremities onto <IntCurve>
    let mut params = [0.0f64; 2];
    for ind_end in 0..2 {
        if ind_end == 1 && a_ba_curve.is_closed() {
            break;
        }
        let projector = point_curve_extrema_all(int_curve.as_ref(), &ends[ind_end]);
        let mut param = [0.0f64; 4];
        let mut dist = [0.0f64; 4];
        // `TrimmedSquareDistances(dist[1], dist[2], Pnt[1], Pnt[2])`.
        param[1] = int_curve.first_parameter();
        param[2] = int_curve.last_parameter();
        dist[1] = ends[ind_end].square_distance(&int_curve.d0(param[1]));
        dist[2] = ends[ind_end].square_distance(&int_curve.d0(param[2]));
        dist[3] = f64::MAX;
        if !projector.is_empty() {
            let mut best = 0usize;
            for i in 1..projector.len() {
                let di = projector[i].distance * projector[i].distance;
                let db = projector[best].distance * projector[best].distance;
                if di < db {
                    best = i;
                }
            }
            // `Projector.Point(imin).Parameter()` / `.SquareDistance(imin)`.
            param[3] = projector[best].u2;
            dist[3] = projector[best].distance * projector[best].distance;
        }

        let mut imin = 1usize;
        for i in 2..=3 {
            if dist[i] < dist[imin] {
                imin = i;
            }
        }
        params[ind_end] = param[imin];
    }
    if a_ba_curve.is_closed() {
        // temporary
        params[1] = int_curve.last_parameter();
    }
    if params[0] > params[1] {
        let fpnt = int_curve.d0(int_curve.first_parameter());
        let lpnt = int_curve.d0(int_curve.last_parameter());
        let is_closed = fpnt.square_distance(&lpnt) <= SQUARE_CONFUSION;
        if is_closed {
            params[1] = int_curve.last_parameter();
        } else {
            let new_first_par = int_curve.reversed_parameter(params[0]);
            let new_last_par = int_curve.reversed_parameter(params[1]);
            int_curve = Arc::from(int_curve.reversed());
            params[0] = new_first_par;
            params[1] = new_last_par;
        }
    }
    if a_ba_curve.is_closed() {
        // check the direction of closed curve
        let (_, tangent) = a_ba_curve.d1(a_ba_curve.first_parameter())?;
        let (_, offset_tangent) = int_curve.d1(params[0]);
        if tangent.dot(&offset_tangent) < 0.0 {
            int_curve = Arc::from(int_curve.reversed());
        }
    }

    // `BRepLib_MakeEdge(IntCurve, Params[0], Params[1])` (`cxx:301`).
    let build = init_3d(&int_curve, params[0], params[1])?;
    let builder = TopoBuilder::new();
    let mut offset_edge = builder.make_edge(int_curve, build.first, build.last);
    add_endpoint_vertices(&builder, &mut offset_edge, &build);
    if build.closed {
        offset_edge.0.set_closed(true);
    }
    Some(offset_edge)
}

/// `TangentExtremity(V, E, hs1, hs2, tang)` (`cxx:625-711`): the two faces of
/// `E` are tangent at `V`, both normals taken on the p-curve of `E`.
///
/// `Ok(false)` covers OCCT's two `return false` arms for an undefined normal;
/// `Err` carries the p-curve / `GeomLProp` failure OCCT would raise.
pub(super) fn tangent_extremity(
    v: &TopoShape,
    e: &Edge,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    tang: f64,
) -> Result<bool, String> {
    let face1 = hs1
        .face()
        .ok_or_else(|| "TangentExtremity: null face on hs1".to_string())?;
    let o1 = face1.0.orientation();
    // `f1.Orientation(TopAbs_FORWARD)` works on the local copy only.
    let f1 = Face(face1.0.oriented(Orientation::Forward));
    let face2 = hs2
        .face()
        .ok_or_else(|| "TangentExtremity: null face on hs2".to_string())?;
    let o2 = face2.0.orientation();
    let f2 = Face(face2.0.oriented(Orientation::Forward));

    let e1 = Edge(e.0.oriented(Orientation::Forward));
    let mut e2 = Edge(e.0.oriented(Orientation::Forward));
    if is_same(&f1.0, &f2.0) && BRepTool::is_closed_edge_face(&e1, &f1) {
        e2 = Edge(e.0.oriented(Orientation::Reversed));
    }
    let p1 = BRepTool::parameter_of_vertex_on_face(v, &e1, &f1);
    let p2 = BRepTool::parameter_of_vertex_on_face(v, &e2, &f2);
    let eps = 1.0e-9;

    let surf1 = hs1
        .surface()
        .ok_or_else(|| "TangentExtremity: null surface on hs1".to_string())?;
    let pc1 = make_2d(&e1, &f1)?;
    let uv1 = pc1.d0(p1);
    let mut prop1 = SlProps::new(surf1.as_ref(), uv1.x(), uv1.y(), 1, eps)?;
    if !prop1.is_normal_defined(surf1.as_ref()) {
        // It is not known...
        return Ok(false);
    }
    let n1 = prop1.normal(surf1.as_ref())?;
    let mut n1 = GpVec::new(n1.x(), n1.y(), n1.z());
    if o1 == Orientation::Reversed {
        n1.reverse();
    }

    let surf2 = hs2
        .surface()
        .ok_or_else(|| "TangentExtremity: null surface on hs2".to_string())?;
    let pc2 = make_2d(&e2, &f2)?;
    let uv2 = pc2.d0(p2);
    let mut prop2 = SlProps::new(surf2.as_ref(), uv2.x(), uv2.y(), 1, eps)?;
    if !prop2.is_normal_defined(surf2.as_ref()) {
        // It is not known...
        return Ok(false);
    }
    let n2 = prop2.normal(surf2.as_ref())?;
    let mut n2 = GpVec::new(n2.x(), n2.y(), n2.z());
    if o2 == Orientation::Reversed {
        n2.reverse();
    }

    Ok(n1.angle(&n2) < tang)
}

/// `TangentOnVertex(V, E, EFMap, tang)` (`cxx:713-728`): the support faces of
/// `E` are tangent at `V`.
pub(super) fn tangent_on_vertex(
    v: &TopoShape,
    e: &Edge,
    ef_map: &ChFiDSMap,
    tang: f64,
) -> Result<bool, String> {
    let (ff1, ff2) = conex_faces(e, ef_map);
    let (Some(ff1), Some(ff2)) = (ff1, ff2) else {
        return Ok(false);
    };
    let s1 = SurfaceAdaptor::of_face(&ff1);
    let s2 = SurfaceAdaptor::of_face(&ff2);
    tangent_extremity(v, e, &s1, &s2, tang)
}
