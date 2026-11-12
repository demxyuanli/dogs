//! `Build` phase 2: the set of intermediate sections.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MiddlePath.cxx`
//! (OCCT 8.0.0): `cxx:522-843`.
//!
//! The loop walks the two neighbouring paths `j-1` and `j` (wrapping to path 1
//! on the last step of a closed section) one segment at a time, pads a path
//! that has already reached its end with the last vertex of its last edge
//! (`cxx:547-558`), and then looks for the edge of `myInitialShape` that joins
//! `PrevVertex` and `CurVertex` on the face shared by `E1`/`E2` and the previous
//! section's edge (`cxx:596-642`). A missing connection is bridged by a line in
//! the face's 2D parametric space (`cxx:660-812`); three candidate lines are
//! built and `IsValidEdge` / `TypeOfEdge` pick one (`cxx:749-788`).
//!
//! Faithfully kept, including the two "restart" devices of the reference: the
//! `ChooseEdge == 1` arm nullifies `EdgeSeq(1..j-2)`, inserts the first vertex
//! of `myPaths(k)(i)` into every path `k <= j-1` and rewinds `j` to 1
//! (`cxx:790-806`).
//!
//! Null-dereference arms: the reference indexes `myPaths`, `SectionsEdges` and
//! `EdgeSeq` without bounds checks and dereferences a null `theFace` / null
//! p-curve / null `Geom2d_Line` inside `BRepLib_MakeEdge`, which is
//! `Standard_OutOfRange` / `Standard_NullObject` / `Standard_ConstructionError`
//! in OCCT. The port answers `Err` for exactly those points and does not add
//! any rule the reference does not have.

use std::collections::HashSet;

use crate::abs::{Orientation, ShapeType};
use crate::brep_offset_api_middle_path::statics::{
    build_curve3d_default, is_valid_edge, make_edge_2d, make_line_2d, pcurve_on_face, type_of_edge,
};
use crate::brep_offset_api_middle_path::{BRepOffsetApiMiddlePath, Ctx};
use crate::brep_tools_wire_explorer::{first_vertex_cum_ori, last_vertex_cum_ori};
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key};

/// `TopoDS::Edge(aShape)`.
fn as_edge(s: &TopoShape) -> Result<Edge, String> {
    Edge::wrap(s.clone()).ok_or_else(|| "Standard_TypeMismatch: TopoDS::Edge".to_string())
}

/// `TopoDS::Vertex(aShape)`.
fn as_vertex(s: &TopoShape) -> Result<Vertex, String> {
    Vertex::wrap(s.clone()).ok_or_else(|| "Standard_TypeMismatch: TopoDS::Vertex".to_string())
}

/// `IsSame` where either side may be a null shape (`IsSame(nullptr)` is false).
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        _ => false,
    }
}

fn null_object(what: &str) -> String {
    format!("Standard_NullObject: {what}")
}

/// `for (;;) { ... }` over `i` (`cxx:530-843`).
pub(super) fn build_sections(
    mp: &mut BRepOffsetApiMiddlePath,
    ctx: &mut Ctx,
) -> Result<(), String> {
    // `cxx:522-528`.
    let nb_e = ctx.edge_seq.len();
    let nb_paths = mp.my_paths.len();
    let nb_ver = if mp.my_closed_section {
        mp.my_paths.len() + 1
    } else {
        mp.my_paths.len()
    };

    // `cxx:529`.
    let mut i: usize = 1;
    loop {
        // `cxx:532-535`: `EdgeSeq(j).Nullify()` for every j.
        for item in ctx.edge_seq.iter_mut() {
            *item = None;
        }
        let mut to_insert_vertex = false;

        // `cxx:539-827`: `for (j = 2; j <= NbVer; j++)`. The reference rewinds
        // `j` to 1 in the `ChooseEdge == 1` arm (`cxx:805`), so this is a
        // `while` loop with the increment written out at every exit.
        let mut j: usize = 2;
        while j <= nb_ver {
            // `cxx:541-544`.
            if ctx.edge_seq[j - 2].is_some() {
                j += 1;
                continue;
            }
            // `myPaths(j-1)` (`cxx:547`, `:559`, ...) and
            // `myPaths((j <= NbPaths) ? j : 1)` (`cxx:555`, ...).
            let p1 = j - 2;
            let p2 = if j <= nb_paths { j - 1 } else { 0 };

            // `cxx:547-558`: pad a path that already reached its end.
            if mp.my_paths[p1].len() < i {
                let a_e1 = as_edge(
                    mp.my_paths[p1]
                        .get(i - 2)
                        .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())?,
                )?;
                let last_ver =
                    last_vertex_cum_ori(&a_e1).ok_or_else(|| null_object("TopExp::LastVertex"))?;
                mp.my_paths[p1].push(last_ver.0);
            }
            if mp.my_paths[p2].len() < i {
                let a_e2 = as_edge(
                    mp.my_paths[p2]
                        .get(i - 2)
                        .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())?,
                )?;
                let last_ver =
                    last_vertex_cum_ori(&a_e2).ok_or_else(|| null_object("TopExp::LastVertex"))?;
                mp.my_paths[p2].push(last_ver.0);
            }

            // `cxx:561-577`: `ToInsertVertex` inserts the first vertex of the
            // segment at `i` in front of it.
            if to_insert_vertex {
                if mp.my_paths[p1][i - 1].shape_type() == ShapeType::Edge {
                    let a_e1 = as_edge(&mp.my_paths[p1][i - 1])?;
                    let fver = first_vertex_cum_ori(&a_e1)
                        .ok_or_else(|| null_object("TopExp::FirstVertex"))?;
                    mp.my_paths[p1].insert(i - 1, fver.0);
                }
                if mp.my_paths[p2][i - 1].shape_type() == ShapeType::Edge {
                    let a_e2 = as_edge(&mp.my_paths[p2][i - 1])?;
                    let fver = first_vertex_cum_ori(&a_e2)
                        .ok_or_else(|| null_object("TopExp::FirstVertex"))?;
                    mp.my_paths[p2].insert(i - 1, fver.0);
                }
                to_insert_vertex = false;
            }

            // `cxx:579-587`.
            let mut e1: Option<Edge> = None;
            let mut e2: Option<Edge> = None;
            if mp.my_paths[p1][i - 1].shape_type() == ShapeType::Edge {
                e1 = Some(as_edge(&mp.my_paths[p1][i - 1])?);
            }
            if mp.my_paths[p2][i - 1].shape_type() == ShapeType::Edge {
                e2 = Some(as_edge(&mp.my_paths[p2][i - 1])?);
            }
            // `cxx:588`: `E12 = TopoDS::Edge(SectionsEdges(i)(j-1))`.
            let e12 = as_edge(
                ctx.sections_edges
                    .get(i - 1)
                    .and_then(|seq| seq.get(j - 2))
                    .and_then(|s| s.as_ref())
                    .ok_or_else(|| "Standard_OutOfRange: SectionsEdges".to_string())?,
            )?;

            // `cxx:590-595`: both null => the proper edge is the previous
            // section's edge between the same two vertices.
            let Some(e1_or_e2) = e1.clone().or_else(|| e2.clone()) else {
                ctx.edge_seq[j - 2] = Some(e12.0);
                j += 1;
                continue;
            };

            // `cxx:596-614`: the face carrying both `E1orE2` and `E12`.
            let lf: Vec<TopoShape> = ctx
                .efmap
                .get(&shape_map_key(&e1_or_e2.0))
                .cloned()
                .unwrap_or_default();
            let lf2: Vec<TopoShape> = ctx
                .efmap
                .get(&shape_map_key(&e12.0))
                .cloned()
                .unwrap_or_default();
            let mut the_face: Option<Face> = None;
            'faces: for a_face in &lf {
                for a_face2 in &lf2 {
                    if is_same(a_face, a_face2) {
                        the_face = Face::wrap(a_face.clone());
                        break 'faces;
                    }
                }
            }
            // `EFmap` is never empty for a section edge the reference built, so
            // a missing face is the null `theFace` dereference of `cxx:629`.
            let the_face = the_face.ok_or_else(|| null_object("theFace"))?;

            // `cxx:616-621`.
            let prev_vertex = match &e1 {
                None => Some(as_vertex(&mp.my_paths[p1][i - 1])?),
                Some(e) => last_vertex_cum_ori(e),
            };
            let cur_vertex = match &e2 {
                None => Some(as_vertex(&mp.my_paths[p2][i - 1])?),
                Some(e) => last_vertex_cum_ori(e),
            };

            // `cxx:626-642`: the edge of the face that joins the two vertices.
            let edges_of_the_face: HashSet<_> = explorer_shapes(&the_face.0, ShapeType::Edge)
                .iter()
                .map(shape_map_key)
                .collect();
            let mut proper_edge: Option<Edge> = None;
            if let Some(pv) = &prev_vertex {
                let le: Vec<TopoShape> = ctx
                    .vemap
                    .find_from_key(&pv.0)
                    .map(|s| s.to_vec())
                    .unwrap_or_default();
                for it in le {
                    let an_edge = as_edge(&it)?;
                    let (v1, v2) = edge_vertices(&an_edge);
                    let joined = (same_vertex(&v1, &prev_vertex) && same_vertex(&v2, &cur_vertex))
                        || (same_vertex(&v1, &cur_vertex) && same_vertex(&v2, &prev_vertex));
                    if !joined {
                        continue;
                    }
                    // `EdgesOfTheFace.Contains(anEdge)`.
                    if !edges_of_the_face.contains(&shape_map_key(&an_edge.0)) {
                        continue;
                    }
                    // `!anEdge.IsSame(E1)`.
                    if e1.as_ref().map_or(false, |e| is_same(&an_edge.0, &e.0)) {
                        continue;
                    }
                    proper_edge = Some(an_edge);
                    break;
                }
            }

            // `cxx:644-649`: both segments are punctual.
            if mp.my_paths[p1][i - 1].shape_type() == ShapeType::Vertex
                && mp.my_paths[p2][i - 1].shape_type() == ShapeType::Vertex
            {
                ctx.edge_seq[j - 2] = proper_edge.map(|e| e.0);
                j += 1;
                continue;
            }

            // `cxx:650-651`.
            let prev_prev_ver = match &e1 {
                None => prev_vertex.clone(),
                Some(e) => first_vertex_cum_ori(e),
            };
            let prev_cur_ver = match &e2 {
                None => cur_vertex.clone(),
                Some(e) => first_vertex_cum_ori(e),
            };

            if let Some(pe) = proper_edge {
                // `cxx:815-826`.
                ctx.edge_seq[j - 2] = Some(pe.0);
                j += 1;
                continue;
            }

            // `cxx:655-658`: no connection between the two vertices.
            let list_one_face: Vec<TopoShape> = vec![the_face.0.clone()];

            if e1.is_none() || e2.is_none() {
                // `cxx:660-700`: one of the two segments is punctual, so the
                // previous element of that path is the edge to use.
                let ee1 = match &e1 {
                    Some(e) => e.clone(),
                    None => as_edge(
                        mp.my_paths[p1]
                            .get(i - 2)
                            .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())?,
                    )?,
                };
                let ee2 = match &e2 {
                    Some(e) => e.clone(),
                    None => as_edge(
                        mp.my_paths[p2]
                            .get(i - 2)
                            .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())?,
                    )?,
                };
                // `cxx:670-690`.
                let (pc1, the_surf, fpar1, lpar1) = pcurve_on_face(&ee1, &the_face)
                    .ok_or_else(|| null_object("BRep_Tool::CurveOnSurface"))?;
                let (pc2, _, fpar2, lpar2) = pcurve_on_face(&ee2, &the_face)
                    .ok_or_else(|| null_object("BRep_Tool::CurveOnSurface"))?;
                let last_par1 = if ee1.0.orientation() == Orientation::Forward {
                    lpar1
                } else {
                    fpar1
                };
                let last_par2 = if ee2.0.orientation() == Orientation::Forward {
                    lpar2
                } else {
                    fpar2
                };
                // `cxx:690-699`.
                let first_pnt2d = pc1.value(last_par1);
                let last_pnt2d = pc2.value(last_par2);
                let the_line = make_line_2d(&first_pnt2d, &last_pnt2d)
                    .ok_or_else(|| "gce_ConfusedPoints: GC_MakeLine2d".to_string())?;
                let len_ne = first_pnt2d.distance(&last_pnt2d);
                let new_edge = make_edge_2d(
                    &the_line,
                    &the_surf,
                    &the_face,
                    prev_vertex.as_ref(),
                    cur_vertex.as_ref(),
                    0.0,
                    len_ne,
                )
                .ok_or_else(|| null_object("BRepLib_MakeEdge"))?;
                // `cxx:698`.
                build_curve3d_default(&new_edge);
                ctx.edge_seq[j - 2] = Some(new_edge.0.clone());
                ctx.efmap
                    .insert(shape_map_key(&new_edge.0), list_one_face.clone());
            }

            if let (Some(ee1), Some(ee2)) = (e1.clone(), e2.clone()) {
                // `cxx:702-812`: both segments are edges.
                let (pc1, the_surf, fpar1, lpar1) = pcurve_on_face(&ee1, &the_face)
                    .ok_or_else(|| null_object("BRep_Tool::CurveOnSurface"))?;
                let (pc2, _, fpar2, lpar2) = pcurve_on_face(&ee2, &the_face)
                    .ok_or_else(|| null_object("BRep_Tool::CurveOnSurface"))?;
                // `cxx:707-726`.
                let (first_par1, last_par1) = if ee1.0.orientation() == Orientation::Forward {
                    (fpar1, lpar1)
                } else {
                    (lpar1, fpar1)
                };
                let (first_par2, last_par2) = if ee2.0.orientation() == Orientation::Forward {
                    (fpar2, lpar2)
                } else {
                    (lpar2, fpar2)
                };
                // `cxx:728-735`.
                let first_pnt2d = pc1.value(last_par1);
                let last_pnt2d = pc2.value(last_par2);
                let the_line = make_line_2d(&first_pnt2d, &last_pnt2d)
                    .ok_or_else(|| "gce_ConfusedPoints: GC_MakeLine2d".to_string())?;
                let len_ne = first_pnt2d.distance(&last_pnt2d);
                let new_edge = make_edge_2d(
                    &the_line,
                    &the_surf,
                    &the_face,
                    prev_vertex.as_ref(),
                    cur_vertex.as_ref(),
                    0.0,
                    len_ne,
                )
                .ok_or_else(|| null_object("BRepLib_MakeEdge"))?;
                // `cxx:736`.
                build_curve3d_default(&new_edge);

                // `cxx:737-747`: the two alternative connections.
                let prev_first_pnt2d = pc1.value(first_par1);
                let prev_last_pnt2d = pc2.value(first_par2);
                let line1 = make_line_2d(&prev_first_pnt2d, &last_pnt2d)
                    .ok_or_else(|| "gce_ConfusedPoints: GC_MakeLine2d".to_string())?;
                let line2 = make_line_2d(&first_pnt2d, &prev_last_pnt2d)
                    .ok_or_else(|| "gce_ConfusedPoints: GC_MakeLine2d".to_string())?;
                let len_ne1 = prev_first_pnt2d.distance(&last_pnt2d);
                let new_edge1 = make_edge_2d(
                    &line1,
                    &the_surf,
                    &the_face,
                    prev_prev_ver.as_ref(),
                    cur_vertex.as_ref(),
                    0.0,
                    len_ne1,
                )
                .ok_or_else(|| null_object("BRepLib_MakeEdge"))?;
                // `cxx:744`.
                build_curve3d_default(&new_edge1);
                let len_ne2 = first_pnt2d.distance(&prev_last_pnt2d);
                let new_edge2 = make_edge_2d(
                    &line2,
                    &the_surf,
                    &the_face,
                    prev_vertex.as_ref(),
                    prev_cur_ver.as_ref(),
                    0.0,
                    len_ne2,
                )
                .ok_or_else(|| null_object("BRepLib_MakeEdge"))?;
                // `cxx:747`.
                build_curve3d_default(&new_edge2);

                // `cxx:749-788`.
                let good_ne = is_valid_edge(&new_edge, &the_face);
                let good_ne1 = is_valid_edge(&new_edge1, &the_face);
                let type_e1 = type_of_edge(&ee1);
                let type_e2 = type_of_edge(&ee2);

                let mut choose_edge = 0u8;
                if !good_ne || type_e1 != type_e2 {
                    if type_e1 == type_e2 {
                        // `cxx:758-768`: `!good_ne`.
                        choose_edge = if good_ne1 { 1 } else { 2 };
                    } else if type_e1 == GeomAbsCurveType::Line {
                        // `cxx:770-778`.
                        choose_edge = 1;
                    } else if type_e2 == GeomAbsCurveType::Line {
                        choose_edge = 2;
                    }
                    // `cxx:779-781`: "to be developed later" leaves
                    // `ChooseEdge == 0`.
                }

                if choose_edge == 0 {
                    // `cxx:784-789`.
                    ctx.edge_seq[j - 2] = Some(new_edge.0.clone());
                    ctx.efmap
                        .insert(shape_map_key(&new_edge.0), list_one_face.clone());
                } else if choose_edge == 1 {
                    // `cxx:790-806`.
                    ctx.edge_seq[j - 2] = Some(new_edge1.0.clone());
                    ctx.efmap
                        .insert(shape_map_key(&new_edge1.0), list_one_face.clone());
                    for k in 0..(j - 2) {
                        ctx.edge_seq[k] = None;
                    }
                    for k in 0..=(j - 2) {
                        let a_last_edge = as_edge(
                            mp.my_paths[k]
                                .get(i - 1)
                                .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())?,
                        )?;
                        let vertex_as_edge = first_vertex_cum_ori(&a_last_edge)
                            .ok_or_else(|| null_object("TopExp::FirstVertex"))?;
                        mp.my_paths[k].insert(i - 1, vertex_as_edge.0);
                    }
                    // `cxx:805`: `j = 1`, i.e. 2 on the next turn.
                    j = 1;
                } else {
                    // `cxx:807-812`.
                    ctx.edge_seq[j - 2] = Some(new_edge2.0.clone());
                    ctx.efmap
                        .insert(shape_map_key(&new_edge2.0), list_one_face.clone());
                    to_insert_vertex = true;
                }
            }

            j += 1;
        }

        // `cxx:828`.
        ctx.sections_edges.push(ctx.edge_seq.clone());

        // `cxx:830-840`: every section edge has to be an end edge.
        let nb_end_edges = ctx
            .edge_seq
            .iter()
            .flatten()
            .filter(|s| ctx.end_edges.contains(&shape_map_key(s)))
            .count();
        if nb_end_edges == nb_e {
            break;
        }

        // `cxx:842`.
        i += 1;
    }

    Ok(())
}
