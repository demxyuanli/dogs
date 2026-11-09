//! The two file-local `Store` overloads of `BRepOffset_Inter2d.cxx`
//! (`:276-390` and `:396-412`): record the intersection vertices of one edge
//! (or of an edge pair) in the 2D `AsDes`, and the coincident pairs in
//! `theDMVV`.
//!
//! Vertex identity throughout is `TopoDS_Shape::IsSame`, i.e.
//! [`is_same`]/[`ShapeMapKey`] (TShape plus location, orientation ignored).
//! `Store` deliberately keeps *every* vertex it is given instead of merging:
//! `FuseVertices` (`cxx:2335-2393`, `fuse_vertices.rs`) fuses the chains
//! afterwards, and `theDMVV` is how it finds them.

use std::collections::HashSet;

use occt_core::gp::GpPnt;
use occt_geom::extrema_pc::point_curve_extrema_range;

use crate::brep_algo::AsDes;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

/// `Store(theEdge, theLV, theTol, IsToUpdate, theAsDes2d, theDMVV)`
/// (`cxx:276-390`).
///
/// Port notes:
/// * `BRep_Builder::UpdateVertex(V, theTol)` (`cxx:290`) is
///   [`GeometryRegistry::raise_vertex_tolerance`] (raise only), the same call
///   `brep_builder_update_vertex` documents.
/// * `GeomAPI_ProjectPointOnCurve(aC, aT1, aT2)` (`cxx:308-310`) followed by
///   `Perform(aP)` (`cxx:351`) is [`point_curve_extrema_range`]; its `None` is
///   `NbPoints() == 0`, the `continue` of `cxx:352-355`. A null
///   `BRep_Tool::Curve` (`cxx:308`) cannot be initialized in OCCT either, so
///   the port skips the projection (same `continue`).
/// * `theAsDes2d->Descendant(theEdge)` is captured once at `cxx:294` and read
///   while the loop appends to the same list; the port clones it, which is what
///   the `NCollection_List` iterator sees.
pub(crate) fn store_on_edge(
    edge: &TopoShape,
    lv: &[TopoShape],
    tol: f64,
    is_to_update: bool,
    as_des: &mut AsDes,
    dmvv: &mut crate::brep_offset::Dmvv,
) {
    let reg = GeometryRegistry::global();

    // `cxx:287-292`: raise every vertex tolerance to `theTol`.
    for v in lv {
        reg.raise_vertex_tolerance(v, tol);
    }

    // `cxx:294`: `aLVEx` is a *reference* to the AsDes list. A bound key keeps
    // growing as the loop below adds to it (`cxx:388`), so every iteration
    // re-scans the list as it stands then; an unbound key answers the static
    // empty list of `Descendant`, which never sees those later `Add`s.
    let live_ex = as_des.has_descendant(edge);

    // `cxx:295-303`: nothing to compare against and nothing to add.
    if !is_to_update && !live_ex {
        if !lv.is_empty() {
            as_des.add_list(edge, lv);
        }
        return;
    }

    // `cxx:305-312`: the projection is only primed for the update arm.
    let mut proj: Option<(std::sync::Arc<dyn occt_geom::Curve>, f64, f64)> = None;
    let mut tol_e = 0.0;
    if is_to_update {
        if let Some(c) = reg.edge_curve(edge) {
            let (t1, t2) = BRepTool::edge_parameters(&crate::shape::Edge(edge.clone()));
            proj = Some((c, t1, t2));
        }
        tol_e = reg.edge_tolerance(edge);
    }

    // `cxx:315-318`: a `TopTools_ShapeMapHasher` set, so `IsSame` repeats are
    // skipped (the second occurrence of the same vertex).
    let mut mv: HashSet<ShapeMapKey> = HashSet::new();
    for v in lv {
        if !mv.insert(shape_map_key(v)) {
            continue;
        }
        // `cxx:322-323`.
        let p: GpPnt = reg.vertex_point(v);
        let a_tol = reg.vertex_tolerance(v);

        // `cxx:327-341`: the already-stored vertices coincident with this one,
        // scanned over the live list (see `live_ex` above).
        let lv_ex: Vec<TopoShape> = if live_ex {
            as_des.descendant(edge).to_vec()
        } else {
            Vec::new()
        };
        let mut lvc: Vec<TopoShape> = Vec::new();
        let mut same_found = false;
        for v_ex in &lv_ex {
            if is_same(v, v_ex) {
                // `cxx:331-334`: `break` leaves the iterator `More()`, the
                // `continue` of `cxx:345-348`.
                same_found = true;
                break;
            }
            let p_ex = reg.vertex_point(v_ex);
            let tol_v_ex = reg.vertex_tolerance(v_ex);
            if p.distance(&p_ex) <= a_tol + tol_v_ex {
                lvc.push(v_ex.clone());
            }
        }
        if same_found {
            continue;
        }

        // `cxx:349-362`.
        if is_to_update {
            let Some((c, t1, t2)) = proj.as_ref() else {
                continue;
            };
            let Some(pair) = point_curve_extrema_range(c.as_ref(), &p, *t1, *t2) else {
                continue;
            };
            if pair.distance > a_tol + tol_e {
                continue;
            }
            // `aV.Oriented(TopAbs_INTERNAL)` then `UpdateVertex(.., aT, theEdge,
            // aTol)`.
            let v_int = v.oriented(crate::abs::Orientation::Internal);
            let _ = update_vertex_on_edge(&v_int, pair.u2, edge, a_tol);
        }

        // `cxx:366-379`: bind the coincident pairs both ways. The first loop
        // appends the *new* vertex under every stored one, the second appends
        // the stored list under the new one.
        if !lvc.is_empty() {
            for vc in &lvc {
                dmvv.add(vc, v);
            }
            dmvv.add_list(v, &lvc);
        }
        // `cxx:388`.
        as_des.add(edge, v);
    }
}

/// `Store(theE1, theE2, theLV1, theLV2, theTol, theAsDes2d, theDMVV)`
/// (`cxx:396-412`): the per-edge overload, `theE1`/`theLV1` first.
pub(crate) fn store_on_two_edges(
    e1: &TopoShape,
    e2: &TopoShape,
    lv1: &[TopoShape],
    lv2: &[TopoShape],
    tol: f64,
    as_des: &mut AsDes,
    dmvv: &mut crate::brep_offset::Dmvv,
) {
    for (e, lv) in [(e1, lv1), (e2, lv2)] {
        store_on_edge(e, lv, tol, false, as_des, dmvv);
    }
}

/// The `TopoDS::Vertex` wrap of a vertex `TopoShape` (`Store` casts its list
/// entries this way).
pub(crate) fn as_vertex(s: &TopoShape) -> Vertex {
    Vertex(s.clone())
}

