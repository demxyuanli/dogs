//! `LocOpe_Generator.cxx:526-1090`: the rebuild of the merged faces.
//!
//! For every base face that has faces fused into it, the wires that still carry
//! a removed or rebuilt edge are rebuilt from
//!
//! 1. the untouched wires of the faces being fused (`cxx:591-688`),
//! 2. the surviving edges of the wire itself, with the p-curves of a really
//!    closed edge paired as the seam overload of `UpdateEdge` wants
//!    (`cxx:697-925`),
//! 3. the edges of the fused faces that no `EdgAdded` entry covers yet, paired
//!    with their twin from `theEEMap` and oriented by the tangent directions at
//!    their shared vertex (`cxx:810-921`, `cxx:929-1068`).
//!
//! `BRepAlgo_Loop` then turns the edge list into wires and faces, and every face
//! of the fusion group answers that same new face list (`cxx:1072-1088`).

mod curves;

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_algo_loop::BRepAlgoLoop;
use crate::brep_offset::{update_edge_pcurve, update_edge_pcurves};
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::brep_tools_uv::uv_bounds_face;
use crate::iterator::cumulated_children;
use crate::locope::shape_maps::{ShapeSet, ShapeValueMap};
use crate::locope::split_shape::empty_copied_face;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

use super::Stage;

/// `cxx:526-1090`.
pub(super) fn fuse_faces(stage: &mut Stage, mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>) {
    let entries: Vec<(TopoShape, Vec<TopoShape>)> = stage
        .ff_map
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    for (fac, list) in entries {
        rebuild_face(stage, mod_shapes, &fac, &list);
    }
}

/// One iteration of the `cxx:531-1089` face loop.
fn rebuild_face(
    stage: &mut Stage,
    mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>,
    fac: &TopoShape,
    list: &[TopoShape],
) {
    let mut mod_face = false;
    let mut list_of_edg: Vec<TopoShape> = Vec::new();
    let mut edg_added = ShapeSet::new();

    // `cxx:537-546`: `orface` is the orientation `fac` carries inside
    // `myShape`. The reference reads a finished explorer when the face is not
    // found, which raises; every key here comes from `myShape`.
    let Some(orface) = explorer_shapes(&stage.shape, ShapeType::Face)
        .iter()
        .find(|f| is_same(f, fac))
        .map(|f| f.orientation())
    else {
        return;
    };
    // `cxx:547-553`.
    let Some(face) = Face::wrap(fac.clone()) else {
        return;
    };
    let mut new_face = empty_copied_face(&face);
    new_face.0.set_orientation(Orientation::Forward);
    let Some(mut s) = BRepTool::face_surface(&face) else {
        return;
    };
    if let Some(basis) = s.rectangular_trimmed_basis() {
        s = basis;
    }
    // `P` (`cxx:452`): the surface only needs p-curves when it is not a plane.
    let p_is_null = s.gp_pln().is_none();
    let mut u_minc = 0.0;
    let mut u_maxc = 0.0;

    for wir in explorer_shapes(&fac.oriented(Orientation::Forward), ShapeType::Wire) {
        let wire_edges = explorer_shapes(&wir, ShapeType::Edge);
        // `cxx:559-568`: does the wire carry a removed or rebuilt edge?
        let wire_touched = wire_edges
            .iter()
            .any(|e| stage.to_remove.contains(e) || mod_shapes.is_bound(e));
        if !wire_touched {
            // `cxx:569-576`: an untouched wire keeps all of its edges.
            list_of_edg.extend(wire_edges.iter().cloned());
            continue;
        }
        if !mod_face {
            // `cxx:580-587`: the UV box of the base face, only for the
            // period realignment of a non-planar surface.
            if p_is_null {
                let (umin, umax, _, _) = uv_bounds_face(&face);
                u_minc = umin;
                u_maxc = umax;
            }
            // `cxx:591-688`: the untouched wires of every face being fused.
            for facbis in list {
                let (facbis, _) = oriented_face(stage, facbis);
                for w in explorer_shapes(&facbis, ShapeType::Wire) {
                    let edges = explorer_shapes(&w, ShapeType::Edge);
                    if edges.iter().any(|e| stage.to_remove.contains(e)) {
                        continue;
                    }
                    for e in &edges {
                        // `cxx:634-635`.
                        list_of_edg.push(e.oriented(Orientation::compose(orface, e.orientation())));
                        edg_added.add(e);
                        // `cxx:639-684`.
                        if p_is_null {
                            curves::put_curve_2d(&new_face.0, &facbis, e, &s, u_minc, u_maxc);
                        }
                    }
                }
            }
            mod_face = true;
        }

        // `cxx:697-925`: rebuild the wire itself.
        let mut c2d1: Option<Arc<dyn Curve2d>> = None;
        for edg in &wire_edges {
            // `cxx:699-706`.
            let orient = edg.orientation();
            if !stage.to_remove.contains(edg) && !stage.ee_map.is_bound(edg) {
                list_of_edg.push(edg.clone());
                continue;
            }
            if !(mod_shapes.is_bound(edg) || stage.ee_map.is_bound(edg)) {
                continue;
            }
            // `cxx:709-715`: the rebuilt edge, else the edge itself.
            let new_edg = if mod_shapes.is_bound(edg) {
                mod_shapes
                    .find(edg)
                    .and_then(|l| l.first().cloned())
                    .unwrap_or_else(|| edg.clone())
            } else {
                edg.clone()
            };
            list_of_edg.push(new_edg.oriented(orient));
            let new_edge_wrapped = Edge::wrap(new_edg.clone());
            let cur_new = new_edge_wrapped.as_ref().and_then(|e| BRepTool::edge_curve_world(e)).map(curves::untrimmed);
            let tol_new = new_edge_wrapped
                .as_ref()
                .map(|e| BRepTool::edge_tolerance(e))
                .unwrap_or(0.0);
            if p_is_null {
                // `cxx:730-765`: the p-curve of the edge on the base face.
                if let Some(c2d) = curve_on_surface_oriented(edg, fac, orient) {
                    let closed = really_closed(edg, fac);
                    if !closed {
                        update_edge_pcurve(&edge_of(&new_edg), &c2d, &new_face, tol_new);
                    } else if c2d1.is_none() {
                        c2d1 = Some(c2d);
                    } else if orient == Orientation::Forward {
                        update_edge_pcurves(
                            &edge_of(&new_edg),
                            &c2d,
                            c2d1.as_ref().expect("bound above"),
                            &new_face,
                            tol_new,
                        );
                    } else {
                        update_edge_pcurves(
                            &edge_of(&new_edg),
                            c2d1.as_ref().expect("bound above"),
                            &c2d,
                            &new_face,
                            tol_new,
                        );
                    }
                }
            }
            // `cxx:767-801`: `AddPart` - the twin edges of `edg` belong to the
            // new wire as well.
            let add_part = add_part(stage, mod_shapes, edg, fac, orient);
            if !add_part {
                continue;
            }
            let Some(twins) = stage.ee_map.find(edg).cloned() else {
                continue;
            };
            for twin in twins {
                // `cxx:809-811`.
                if edg_added.contains(&twin) {
                    continue;
                }
                // `cxx:812-848`: the tangent directions at the shared vertex.
                let shared = cumulated_children(&new_edg)
                    .into_iter()
                    .find(|c| cumulated_children(&twin).iter().any(|d| is_same(c, d)));
                let (dir1, dir2) = match &shared {
                    Some(v) => {
                        let d1 = cur_new.as_ref().map(|c| {
                            let prm = BRepTool::parameter_of_vertex(&new_edg, v);
                            c.d1(prm).1
                        });
                        let d2 = Edge::wrap(twin.clone())
                            .and_then(|e| BRepTool::edge_curve_world(&e))
                            .map(curves::untrimmed)
                            .map(|c| {
                                let prm = BRepTool::parameter_of_vertex(&twin, v);
                                c.d1(prm).1
                            });
                        (d1, d2)
                    }
                    // `cxx:849-851`.
                    None => (None, None),
                };
                edg_added.add(&twin);
                // `cxx:856-863`: reverse the twin when the tangents disagree.
                // A missing curve where the reference dereferences null keeps
                // the twin forward.
                let dot = match (&dir1, &dir2) {
                    (Some(a), Some(b)) => a.dot(b),
                    _ => 1.0,
                };
                let twin_orient = if dot < 0.0 { orient.reversed() } else { orient };
                list_of_edg.push(twin.oriented(twin_orient));
                // `cxx:864-921`.
                if p_is_null {
                    curves::put_twin_curve_2d(
                        &new_face.0,
                        fac,
                        edg,
                        &twin,
                        &s,
                        tol_new,
                        stage.dont_fuse.find(edg),
                        mod_shapes.is_bound(edg),
                    );
                }
            }
        }

        // `cxx:929-1068`: the edges of the fused faces that are still missing.
        let mut include_in_w = false;
        for facbis in list {
            let (facbis, genface) = oriented_face(stage, facbis);
            for w in explorer_shapes(&facbis, ShapeType::Wire) {
                let edges = explorer_shapes(&w, ShapeType::Edge);
                let touched = edges.iter().any(|e| stage.to_remove.contains(e));
                if genface {
                    include_in_w = true;
                } else {
                    // The reference computes `isIncludedInW` from the removed
                    // edges' vertices here (`cxx:976-1004`) and never reads it;
                    // `includeinw` keeps whatever the `genface` arm above left.
                    let _ = touched;
                }
                if !include_in_w {
                    continue;
                }
                for e in &edges {
                    if stage.to_remove.contains(e) || edg_added.contains(e) {
                        continue;
                    }
                    // `cxx:1016-1020`.
                    list_of_edg.push(e.oriented(Orientation::compose(orface, e.orientation())));
                    edg_added.add(e);
                    // `cxx:1021-1065`.
                    if p_is_null {
                        curves::put_curve_2d(&new_face.0, &facbis, e, &s, u_minc, u_maxc);
                    }
                }
            }
        }
    }

    // `cxx:1072-1088`.
    if list_of_edg.is_empty() {
        return;
    }
    let mut loop_builder = BRepAlgoLoop::new();
    loop_builder.init(&new_face);
    loop_builder.add_const_edges(&list_of_edg);
    loop_builder.perform();
    loop_builder.wires_to_faces();
    let list_of_faces: Vec<TopoShape> = loop_builder.new_faces().to_vec();
    stage.to_remove.add(fac);
    mod_shapes.bind(fac, list_of_faces.clone());
    for facbis in list {
        mod_shapes.bind(facbis, list_of_faces.clone());
    }
}

/// `cxx:933-958`: the face with the orientation `G->OrientedFaces()` gives it,
/// else the orientation it carries in `myShape`. The flag tells whether it came
/// from the generated list (`genface`).
fn oriented_face(stage: &Stage, facbis: &TopoShape) -> (TopoShape, bool) {
    for f in stage.g.oriented_faces() {
        if is_same(f, facbis) {
            return (facbis.oriented(f.orientation()), true);
        }
    }
    for f in explorer_shapes(&stage.shape, ShapeType::Face) {
        if is_same(&f, facbis) {
            return (facbis.oriented(f.orientation()), false);
        }
    }
    (facbis.clone(), false)
}

/// `cxx:767-801`.
fn add_part(
    stage: &Stage,
    mod_shapes: &ShapeValueMap<Vec<TopoShape>>,
    edg: &TopoShape,
    fac: &TopoShape,
    orient: Orientation,
) -> bool {
    let Some(dont) = stage.dont_fuse.find(edg) else {
        return false;
    };
    if !really_closed(edg, fac) {
        if is_same(dont, fac) {
            mod_shapes.is_bound(edg)
        } else {
            !stage.to_remove.contains(edg)
        }
    } else if mod_shapes.is_bound(edg) {
        orient == dont.orientation()
    } else {
        orient == dont.orientation().reversed()
    }
}

/// `BRep_Tool::CurveOnSurface(Edge(E.Oriented(Ori)), Face(F.Oriented(FORWARD)))`
/// (`cxx:736-739`), without the range the reference overwrites anyway.
fn curve_on_surface_oriented(
    edg: &TopoShape,
    fac: &TopoShape,
    orient: Orientation,
) -> Option<Arc<dyn Curve2d>> {
    let e = Edge::wrap(edg.oriented(orient))?;
    let f = Face::wrap(fac.oriented(Orientation::Forward))?;
    curve_on_surface_range(&e, &f).map(|(c, _, _)| c)
}

/// `BRepTools::IsReallyClosed(const TopoDS_Edge&, const TopoDS_Face&)`.
fn really_closed(edg: &TopoShape, fac: &TopoShape) -> bool {
    match (Edge::wrap(edg.clone()), Face::wrap(fac.clone())) {
        (Some(e), Some(f)) => is_really_closed(&e, &f),
        _ => false,
    }
}

/// `TopoDS::Edge` on a raw shape.
fn edge_of(e: &TopoShape) -> Edge {
    Edge::wrap(e.clone()).expect("LocOpe_Generator: expected an edge")
}
