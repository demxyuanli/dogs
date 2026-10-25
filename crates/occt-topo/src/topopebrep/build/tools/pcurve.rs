//! `TopOpeBRepBuild_Tools` — p-curve recomputation and face rebuilding.
//!
//! Sources: `TopOpeBRepBuild_Tools.cxx:664-677` (`UpdatePCurves`), `:679-714`
//! (`UpdateEdgeOnFace`), `:716-836` (`UpdateEdgeOnPeriodicalFace`), `:855-906`
//! (`NormalizeFace`), `:908-924` (`CorrectFace2d`) — OCCT 8.0.0.
//!
//! `BRep_Builder::UpdateEdge(E, C2d, F, Tol)` and the seam overload
//! `UpdateEdge(E, C1, C2, F, Tol)` map onto
//! [`GeometryRegistry::set_edge_pcurve`] / [`GeometryRegistry::set_edge_pcurves`];
//! the tolerance write is `BRep_TEdge::UpdateTolerance` (`BRep_TEdge.hxx:50-52`),
//! i.e. a maximum, so it is written as `max(current, tol)`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::correct_face_2d::CorrectFace2d;
use crate::topopebrep::build::tools::geom::get_normal_to_face_on_edge;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

/// `FC2D_CurveOnSurface(E, F, f2, l2, tolpc, true)`
/// (`TopOpeBRepTool_2d.cxx:356-376`): the p-curve and its `BRep_GCurve` range.
///
/// UNPORTED: `TolReached2d` from the `FC2D_make2d` projection arm
/// (`TopOpeBRepTool_2d.cxx:301`, `TopOpeBRepTool_CurveTool::MakePCurveOnFace`)
/// — the port's projection has no tolerance output (see the same note in
/// `topopebrep/interf/hctxee2d.rs:23-31`). Every caller folds it into
/// `max(tolE, tolpc)`, and the stored-p-curve arm sets `tol = Tolerance(E)`
/// (`2d.cxx:165`), so the tolerance written below is `Tolerance(E)`.
fn fc2d_curve_on_surface(e: &Edge, f: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    crate::boptools_2d::build_pcurve_for_edge_on_face(e, f).ok()?;
    crate::boptools_2d::curve_on_surface_range(e, f)
}

/// `UpdatePCurves(aWire, fromFace, toFace)` (`cxx:664-677`): recomputes the
/// p-curve of every edge of `aWire` on `toFace`.
pub fn update_pcurves(a_wire: &TopoShape, from_face: &Face, to_face: &Face) {
    for e in explorer_shapes(a_wire, crate::abs::ShapeType::Edge) {
        if let Some(edge) = Edge::wrap(e) {
            update_edge_on_face(&edge, from_face, to_face);
        }
    }
}

/// `BRep_Builder::UpdateEdge(E, C, F, Tol)` (`BRep_Builder.cxx:1180-1210`):
/// writes the p-curve on `F` and raises the edge tolerance to `tol`.
pub fn write_edge_pcurve(e: &Edge, c: Arc<dyn Curve2d>, f: &Face, tol: f64) {
    let reg = GeometryRegistry::global();
    reg.set_edge_pcurve(&e.0, GeometryRegistry::shape_key(&f.0), c);
    update_tolerance(e, tol);
}

/// The seam overload `BRep_Builder::UpdateEdge(E, C1, C2, F, Tol)`.
pub fn write_edge_pcurves(e: &Edge, c1: Arc<dyn Curve2d>, c2: Arc<dyn Curve2d>, f: &Face, tol: f64) {
    let reg = GeometryRegistry::global();
    reg.set_edge_pcurves(&e.0, GeometryRegistry::shape_key(&f.0), vec![c1, c2]);
    update_tolerance(e, tol);
}

/// `BRep_TEdge::UpdateTolerance(T)` (`BRep_TEdge.hxx:50-52`): the max of `T`
/// and the current tolerance.
fn update_tolerance(e: &Edge, tol: f64) {
    let reg = GeometryRegistry::global();
    let current = BRepTool::edge_tolerance(e);
    if tol > current {
        reg.set_edge_tolerance(&e.0, tol);
    }
}

/// `UpdateEdgeOnFace(aEdgeToUpdate, fromFace, toFace)` (`cxx:679-714`).
pub fn update_edge_on_face(a_edge_to_update: &Edge, from_face: &Face, to_face: &Face) {
    let tol_e = BRepTool::edge_tolerance(a_edge_to_update);

    if BRepTool::is_degenerated(a_edge_to_update) {
        // A degenerated edge has no computable p-curve on `toFace`; keep the one
        // of `fromFace` as a trimmed copy (`cxx:689-699`).
        let Some((c2d, f2, l2)) = fc2d_curve_on_surface(a_edge_to_update, from_face) else {
            return;
        };
        // `tolpc` — see the UNPORTED note on [`fc2d_curve_on_surface`].
        let tol = tol_e;
        let c2dn: Arc<dyn Curve2d> = Arc::from(c2d.transformed(&GpTrsf2d::identity()));
        let new_c2d: Arc<dyn Curve2d> =
            Arc::new(occt_geom2d::trimmed::Geom2dTrimmedCurve::new(c2dn, f2, l2));
        write_edge_pcurve(a_edge_to_update, new_c2d, to_face, tol);
    } else if BRepTool::is_closed_edge_face(a_edge_to_update, from_face) {
        update_edge_on_periodical_face(a_edge_to_update, from_face, to_face);
    } else {
        let Some((c2d, _f2, _l2)) = fc2d_curve_on_surface(a_edge_to_update, to_face) else {
            return;
        };
        write_edge_pcurve(a_edge_to_update, c2d, to_face, tol_e);
    }
}

/// `UpdateEdgeOnPeriodicalFace(aEdgeToUpdate, fromFace, toFace)`
/// (`cxx:716-836`): recomputes both seam p-curves, choosing their order from
/// the 2D direction of the new curve and the 2D offset between the two sides.
pub fn update_edge_on_periodical_face(a_edge_to_update: &Edge, from_face: &Face, to_face: &Face) {
    // `BRep_Tool::CurveOnSurface(newE, tFace, fc, lc)` non-null -> nothing to do
    // (`cxx:729-735`).
    if let Some((_pc, _fc, _lc)) =
        crate::boptools_2d::curve_on_surface_oriented(a_edge_to_update, to_face, false)
    {
        return;
    }

    // `cxx:738-749`: opposite normals mean the two faces see the edge reversed.
    let a_n1 = get_normal_to_face_on_edge(from_face, a_edge_to_update);
    let a_n2 = get_normal_to_face_on_edge(to_face, a_edge_to_update);
    let diff_oriented = match (a_n1, a_n2) {
        (Some(n1), Some(n2)) => n1.dot(&n2) < 0.0,
        _ => false,
    };

    let tol_e = BRepTool::edge_tolerance(a_edge_to_update);
    // `tol = max(tolpc, tolE)` with `tolpc` unported — see [`fc2d_curve_on_surface`].
    let tol = tol_e;

    // First p-curve (`cxx:751-757`).
    let Some((c2d, f2, l2)) = fc2d_curve_on_surface(a_edge_to_update, to_face) else {
        return;
    };

    // `cxx:766-778`: the p-curves of both edge orientations on `fFace`.
    let Some((old_c2d_for, ff, _lf)) =
        crate::boptools_2d::curve_on_surface_oriented(a_edge_to_update, from_face, false)
    else {
        return;
    };
    let mut new_e_rev = a_edge_to_update.0.clone();
    new_e_rev.set_orientation(Orientation::Reversed);
    let Some((old_c2d_rev, fr, _lr)) = crate::boptools_2d::curve_on_surface_oriented(
        &Edge(new_e_rev),
        from_face,
        false,
    ) else {
        return;
    };

    let a_uvf = old_c2d_for.d0(ff);
    let a_uvr = old_c2d_rev.d0(fr);

    let a_tr_v = if !diff_oriented {
        GpVec2d::new(a_uvr.x() - a_uvf.x(), a_uvr.y() - a_uvf.y())
    } else {
        GpVec2d::new(a_uvf.x() - a_uvr.x(), a_uvf.y() - a_uvr.y())
    };

    // `cxx:788-789`: `aux = gp_Vec2d((0,0), (1,1))`.
    let aux = GpVec2d::new(1.0, 1.0);
    let dir = aux.dot(&a_tr_v) >= 0.0;

    // `cxx:791-809`: the order of the two seam p-curves.
    let a_uvfv = c2d.d0(f2);
    let a_uvlv = c2d.d0(l2);
    let c2d_vec = GpVec2d::new(a_uvlv.x() - a_uvfv.x(), a_uvlv.y() - a_uvfv.y());

    let scalar = GpVec2d::new(0.0, 1.0).dot(&c2d_vec);
    let first_order = if scalar.abs() <= 1e-10 {
        // Compute along the X axis instead.
        GpVec2d::new(1.0, 0.0).dot(&c2d_vec) >= 0.0
    } else {
        scalar <= 0.0
    };

    // `aTrC = C2D->Copy(); aTrC->Translate(aTrV)` (`cxx:810-812`).
    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(&a_tr_v);
    let a_tr_c: Arc<dyn Curve2d> = Arc::from(c2d.transformed(&tr));

    if dir {
        if first_order {
            write_edge_pcurves(a_edge_to_update, c2d, a_tr_c, to_face, tol);
        } else {
            write_edge_pcurves(a_edge_to_update, a_tr_c, c2d, to_face, tol);
        }
    } else if !first_order {
        write_edge_pcurves(a_edge_to_update, c2d, a_tr_c, to_face, tol);
    } else {
        write_edge_pcurves(a_edge_to_update, a_tr_c, c2d, to_face, tol);
    }
}

/// `NormalizeFace(oldFace, corrFace)` (`cxx:855-906`): rebuilds the face with
/// one forward wire per source wire, dropping the `INTERNAL` / `EXTERNAL`
/// edges. A wire that ends up empty is dropped (`cxx:897-902`).
pub fn normalize_face(old_face: &TopoShape) -> Option<TopoShape> {
    let Some(mut a_f1) = Face::wrap(old_face.clone()) else {
        return None;
    };
    a_f1.0.set_orientation(Orientation::Forward);

    let surf = BRepTool::face_surface_world(&a_f1)?;
    let tol_f1 = BRepTool::face_tolerance(&a_f1);
    let builder = TopoBuilder::new();
    let mut a_new_face = builder.make_face_with_tolerance(surf, &[], tol_f1);

    for wire in explorer_shapes(&a_f1.0, crate::abs::ShapeType::Wire) {
        let mut a_wire = wire;
        let wire_orientation = a_wire.orientation();
        a_wire.set_orientation(Orientation::Forward);

        let mut good_edges: Vec<Edge> = Vec::new();
        for e in explorer_shapes(&a_wire, crate::abs::ShapeType::Edge) {
            let o = e.orientation();
            if o == Orientation::External || o == Orientation::Internal {
                continue;
            }
            if let Some(edge) = Edge::wrap(e) {
                good_edges.push(edge);
            }
        }
        if good_edges.is_empty() {
            continue;
        }
        let mut a_new_wire = builder.make_wire(&good_edges);
        a_new_wire.0.set_orientation(wire_orientation);
        builder.add_wire(&mut a_new_face, &a_new_wire);
    }

    a_new_face.0.set_orientation(old_face.orientation());
    Some(a_new_face.0)
}

/// `CorrectFace2d(oldFace, corrFace, aSourceShapes, aMapOfCorrect2dEdges)`
/// (`cxx:908-924`).
///
/// OCCT runs the corrector on `oldFace` and then assigns `corrFace = oldFace`:
/// the moves are written in place on the original edges recorded in
/// `a_map_of_correct2d_edges`, so the caller's `old_face` is the corrected one.
pub fn correct_face_2d(
    old_face: &Face,
    a_source_shapes: &IndexedShapeSet,
    a_map_of_correct2d_edges: Rc<RefCell<IndexedShapeData<TopoShape>>>,
) -> TopoShape {
    let mut a_correct_face2d =
        CorrectFace2d::with_face(old_face, a_source_shapes, a_map_of_correct2d_edges);
    a_correct_face2d.perform();
    old_face.0.clone()
}
