//! Remainder of `BOPAlgo_BuilderFace.cxx`: `PerformShapesToAvoid`,
//! `PerformLoops`, `PerformAreas`, `PerformInternalShapes`, `IsGrowthWire`,
//! `IsInside`, `MakeInternalWires`.
//!
//! Source: `BOPAlgo_BuilderFace.cxx` (Perform at 117, ShapesToAvoid at 152,
//! Loops at 239, Areas at 387, InternalShapes at 618, MakeInternalWires at
//! 782, IsInside at 842).

use std::collections::{HashMap, HashSet};

use occt_core::gp::GpPnt2d;

use crate::abs::{Orientation, ShapeType};
use crate::bbox_from_geometry::shape_bbox;
use crate::bop_box2d_tree::{Box2dTree, Box2dTreeSelector};
use crate::boptools_2d::curve_on_surface;
use crate::brep_tool::BRepTool;
use crate::brep_uv_bounds::{add_uv_bounds_edge, add_uv_bounds_face};
use crate::builder::TopoBuilder;
use crate::builder_face::FaceBuilder;
use crate::fclass2d::{FaceState, FClass2d};
use crate::int_tools_full::IntToolsContext;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::{FaceGeom, GeometryRegistry};
use crate::topo_tools_full::{edge_vertices, edges_of, vertices_of};
use crate::wire_splitter::{WireEdgeSet, WireSplitter};

fn sk(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

fn map_vertex_edges(edges: &[TopoShape]) -> HashMap<usize, Vec<TopoShape>> {
    let mut mve: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for e in edges {
        for v in vertices_of(e) {
            mve.entry(sk(&v.0)).or_default().push(e.clone());
        }
    }
    mve
}

/// `BOPAlgo_BuilderFace::PerformShapesToAvoid`.
pub fn perform_shapes_to_avoid(fb: &mut FaceBuilder) -> Result<(), String> {
    fb.base.avoid.clear();
    loop {
        let mut b_found = false;
        let live: Vec<TopoShape> = fb
            .base
            .shapes
            .iter()
            .filter(|s| s.shape_type() == ShapeType::Edge && !fb.base.is_avoided(s))
            .cloned()
            .collect();
        let mve = map_vertex_edges(&live);
        for (_vk, a_le) in &mve {
            let a_nb_e = a_le.len();
            if a_nb_e == 0 {
                continue;
            }
            let a_e1 = &a_le[0];
            if a_nb_e == 1 {
                if BRepTool::is_degenerated(&Edge(a_e1.clone())) {
                    continue;
                }
                let internal_v = vertices_of(a_e1).iter().any(|v| v.0.orientation() == Orientation::Internal);
                if internal_v {
                    continue;
                }
                b_found = true;
                if !fb.base.avoid.iter().any(|a| a.same_tshape(a_e1)) {
                    fb.base.avoid.push(a_e1.clone());
                }
            } else if a_nb_e == 2 {
                let a_e2 = &a_le[1];
                if a_e2.same_tshape(a_e1) {
                    let (v1x, v2x) = edge_vertices(&Edge(a_e1.clone()));
                    let same_ends = match (v1x, v2x) {
                        (Some(a), Some(b)) => a.0.same_tshape(&b.0),
                        _ => false,
                    };
                    if same_ends {
                        continue;
                    }
                    b_found = true;
                    if !fb.base.avoid.iter().any(|a| a.same_tshape(a_e1)) {
                        fb.base.avoid.push(a_e1.clone());
                    }
                    if !fb.base.avoid.iter().any(|a| a.same_tshape(a_e2)) {
                        fb.base.avoid.push(a_e2.clone());
                    }
                }
            }
        }
        if !b_found {
            break;
        }
    }
    Ok(())
}

/// `BOPAlgo_BuilderFace::PerformLoops`.
pub fn perform_loops(fb: &mut FaceBuilder) -> Result<(), String> {
    fb.base.loops.clear();
    let mut a_wes = WireEdgeSet::new();
    if let Some(face) = fb.face.clone() {
        a_wes.set_face(face);
    }
    for s in &fb.base.shapes {
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        if fb.base.is_avoided(s) {
            continue;
        }
        a_wes.add_edge(Edge(s.clone()));
    }
    if a_wes.is_empty() {
        return Ok(());
    }
    let mut a_wsp = WireSplitter::new();
    a_wsp.set_wes(a_wes);
    if let Err(e) = a_wsp.perform() {
        fb.base.warnings.push(e);
        return Ok(());
    }
    for w in a_wsp.wires() {
        fb.base.loops.push(w.clone());
    }
    let mut a_mep: HashSet<usize> = HashSet::new();
    for w in &fb.base.loops {
        for e in edges_of(w) {
            a_mep.insert(sk(&e.0));
        }
    }
    for a in &fb.base.avoid {
        a_mep.insert(sk(a));
    }
    for s in &fb.base.shapes {
        if s.shape_type() == ShapeType::Edge && !a_mep.contains(&sk(s)) {
            if !fb.base.avoid.iter().any(|a| a.same_tshape(s)) {
                fb.base.avoid.push(s.clone());
            }
        }
    }
    fb.base.internal.clear();
    let a_nb_ea = fb.base.avoid.len();
    if a_nb_ea == 0 {
        return Ok(());
    }
    let mve = map_vertex_edges(&fb.base.avoid);
    let mut a_m_added: HashSet<usize> = HashSet::new();
    let mut b_flag = true;
    let avoid = fb.base.avoid.clone();
    for a_ee in &avoid {
        if !b_flag {
            break;
        }
        if !a_m_added.insert(sk(a_ee)) {
            continue;
        }
        let mut edges: Vec<Edge> = vec![Edge(a_ee.clone())];
        let mut growing = true;
        while growing && b_flag {
            growing = false;
            let current = edges.clone();
            for e in &current {
                for v in vertices_of(&e.0) {
                    if let Some(a_le) = mve.get(&sk(&v.0)) {
                        for a_ex in a_le {
                            if a_m_added.insert(sk(a_ex)) {
                                edges.push(Edge(a_ex.clone()));
                                growing = true;
                                if a_m_added.len() == a_nb_ea {
                                    b_flag = false;
                                }
                            }
                        }
                    }
                }
            }
        }
        let a_w = TopoBuilder::new().make_wire(&edges);
        fb.base.internal.push(a_w.0);
    }
    Ok(())
}

/// `IsGrowthWire` (`BOPAlgo_BuilderFace.cxx`).
pub fn is_growth_wire(the_wire: &TopoShape, a_mhe: &HashSet<usize>) -> bool {
    for e in edges_of(the_wire) {
        if a_mhe.contains(&sk(&e.0)) {
            return true;
        }
    }
    false
}

/// `IsInside` (`BOPAlgo_BuilderFace.cxx:842`).
///
/// Shared-edge rule: if the face already contains an edge of the wire, OCCT
/// **returns false immediately** (the wire cannot lie inside that face).
pub fn is_inside(
    the_wire_or_edge: &TopoShape,
    the_f: &Face,
    the_ctx: &mut IntToolsContext,
) -> bool {
    let mut a_face_edges: HashSet<usize> = HashSet::new();
    for e in edges_of(&the_f.0) {
        a_face_edges.insert(sk(&e.0));
    }
    let edges: Vec<Edge> = if the_wire_or_edge.shape_type() == ShapeType::Edge {
        vec![Edge(the_wire_or_edge.clone())]
    } else {
        edges_of(the_wire_or_edge)
    };
    let Ok(a_clsf) = the_ctx.fclass2d(the_f) else {
        return false;
    };
    let mut is_inside = false;
    for e in &edges {
        if BRepTool::is_degenerated(e) {
            continue;
        }
        if a_face_edges.contains(&sk(&e.0)) {
            return is_inside;
        }
        let Some((c2d, a_t1, a_t2)) = crate::boptools_2d::curve_on_surface_range(e, the_f) else {
            continue;
        };
        let p2 = c2d.d0(0.5 * (a_t1 + a_t2));
        let st = a_clsf.perform(GpPnt2d::new(p2.x(), p2.y()));
        is_inside = st == FaceState::In;
        break;
    }
    is_inside
}

fn first_wire_of(shape: &TopoShape) -> Option<Wire> {
    for w in ShapeIterator::of_shape(shape) {
        if w.shape_type() == ShapeType::Wire {
            return Some(Wire(w));
        }
    }
    None
}

/// `BRep_Builder::MakeFace(face, S, Loc, Tol)` — empty face on the generatrix
/// surface. Orientation stays FORWARD (`myFace` is already FORWARD).
fn make_empty_face(fb: &FaceBuilder, natural: bool) -> Result<Face, String> {
    let Some(gf) = fb.face.as_ref() else {
        return Err("builder_face: no generatrix face".into());
    };
    let Some(surf) = BRepTool::face_surface(gf) else {
        return Err("builder_face: no surface".into());
    };
    let a_tol = BRepTool::face_tolerance(gf);
    let mut geom = FaceGeom::new(surf);
    geom.tolerance = a_tol;
    geom.natural_restriction = natural;
    let face = Face::new();
    GeometryRegistry::global().set_face(&face.0, geom);
    Ok(face)
}

fn make_face_on(fb: &FaceBuilder, wire: &Wire) -> Result<Face, String> {
    let Some(gf) = fb.face.as_ref() else {
        return Err("builder_face: no generatrix face".into());
    };
    let Some(surf) = BRepTool::face_surface(gf) else {
        return Err("builder_face: no surface".into());
    };
    let a_tol = BRepTool::face_tolerance(gf);
    let mut face = TopoBuilder::new().make_face(surf, &[wire.clone()]);
    GeometryRegistry::global().set_face_tolerance(&face.0, a_tol);
    GeometryRegistry::global().set_natural_restriction(&face.0, false);
    // OCCT CurveOnSurface is keyed by (edge, surface). Copy the generatrix
    // pcurves onto the draft face so FClass2d Init sees the same UV as
    // WireSplitter Coord2d.
    for e in edges_of(&wire.0) {
        if curve_on_surface(&e, &face).is_some() {
            continue;
        }
        if let Some(pc) = curve_on_surface(&e, gf) {
            let _ = crate::boptools_2d::attach_existing_pcurve(&e, &face, pc);
        }
    }
    Ok(face)
}

/// `BOPAlgo_BuilderFace::PerformAreas` (`BOPAlgo_BuilderFace.cxx:387-614`).
pub fn perform_areas(fb: &mut FaceBuilder) -> Result<(), String> {
    fb.base.areas.clear();
    let Some(gf) = fb.face.clone() else {
        return Ok(());
    };
    let a_tol = BRepTool::face_tolerance(&gf);
    let mut ctx = IntToolsContext::new();

    if fb.base.loops.is_empty() {
        if ctx.is_infinite_face(&gf) {
            let natural = BRepTool::natural_restriction(&gf);
            let a_face = make_empty_face(fb, natural)?;
            fb.base.areas.push(a_face.0);
        }
        return Ok(());
    }

    let mut a_new_faces: Vec<TopoShape> = Vec::new();
    let mut a_hole_faces: Vec<TopoShape> = Vec::new();
    let mut a_mhe: HashSet<usize> = HashSet::new();

    for a_wire in fb.base.loops.clone() {
        let wire = Wire(a_wire.clone());
        let a_face = make_face_on(fb, &wire)?;
        let mut b_is_growth = is_growth_wire(&a_wire, &a_mhe);
        if !b_is_growth {
            if let Ok(a_clsf) = FClass2d::new(&a_face, a_tol) {
                b_is_growth = !a_clsf.is_hole();
            }
        }
        if b_is_growth {
            a_new_faces.push(a_face.0);
        } else {
            a_hole_faces.push(a_face.0);
            for e in edges_of(&a_wire) {
                a_mhe.insert(sk(&e.0));
            }
        }
    }

    if a_hole_faces.is_empty() {
        fb.base.areas.extend(a_new_faces);
        return Ok(());
    }

    let a_nb_h = a_hole_faces.len();
    let mut a_box_tree = Box2dTree::new();
    a_box_tree.set_size(a_nb_h);
    for (i, a_h_face) in a_hole_faces.iter().enumerate() {
        let mut a_box = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&Face(a_h_face.clone()), &mut a_box);
        a_box_tree.add(i as i32, a_box);
    }
    a_box_tree.build();

    let mut a_hole_face_map: HashMap<usize, TopoShape> = HashMap::new();
    let mut a_selector = Box2dTreeSelector::new();
    a_selector.set_bvh_set(&a_box_tree);

    for a_face in &a_new_faces {
        let af = Face(a_face.clone());
        let mut a_box = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&af, &mut a_box);
        a_selector.clear();
        a_selector.set_box(a_box);
        let a_li = a_selector.select();
        for k in a_li {
            let idx = k as usize;
            if idx >= a_hole_faces.len() {
                continue;
            }
            let a_hole = &a_hole_faces[idx];
            if !is_inside(a_hole, &af, &mut ctx) {
                continue;
            }
            let hk = sk(a_hole);
            if let Some(was) = a_hole_face_map.get(&hk).cloned() {
                if is_inside(a_face, &Face(was), &mut ctx) {
                    a_hole_face_map.insert(hk, a_face.clone());
                }
            } else {
                a_hole_face_map.insert(hk, a_face.clone());
            }
        }
    }

    let mut a_face_holes: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for (hk, a_face) in &a_hole_face_map {
        let Some(hole) = a_hole_faces.iter().find(|h| sk(h) == *hk) else {
            continue;
        };
        a_face_holes
            .entry(sk(a_face))
            .or_default()
            .push(hole.clone());
    }

    if a_hole_faces.len() != a_hole_face_map.len() {
        let a_box_f = shape_bbox(&gf.0);
        if a_box_f.is_open_xmin()
            || a_box_f.is_open_xmax()
            || a_box_f.is_open_ymin()
            || a_box_f.is_open_ymax()
            || a_box_f.is_open_zmin()
            || a_box_f.is_open_zmax()
        {
            let a_face = make_empty_face(fb, false)?;
            let mut unused: Vec<TopoShape> = Vec::new();
            for a_hole in &a_hole_faces {
                if !a_hole_face_map.contains_key(&sk(a_hole)) {
                    unused.push(a_hole.clone());
                }
            }
            a_face_holes.insert(sk(&a_face.0), unused);
            a_new_faces.push(a_face.0);
        }
    }

    let b = TopoBuilder::new();
    for a_face in a_new_faces {
        let holes = a_face_holes.get(&sk(&a_face)).cloned().unwrap_or_default();
        let mut face = Face(a_face);
        if !holes.is_empty() {
            for a_f_hole in &holes {
                if let Some(w) = first_wire_of(a_f_hole) {
                    b.add_wire(&mut face, &w);
                }
            }
            ctx.fclass2d_invalidate(&face);
            let _ = ctx.fclass2d(&face);
            let _ = a_tol;
        }
        fb.base.areas.push(face.0);
    }
    Ok(())
}

/// `MakeInternalWires`.
pub fn make_internal_wires(the_me: &[TopoShape]) -> Vec<TopoShape> {
    let mut the_wires: Vec<TopoShape> = Vec::new();
    let mve = map_vertex_edges(the_me);
    let mut a_added: HashSet<usize> = HashSet::new();
    for a_ee in the_me {
        if !a_added.insert(sk(a_ee)) {
            continue;
        }
        let mut edges: Vec<Edge> = Vec::new();
        let mut e0 = a_ee.clone();
        e0.set_orientation(Orientation::Internal);
        edges.push(Edge(e0));
        let mut growing = true;
        while growing {
            growing = false;
            let current = edges.clone();
            for e in &current {
                for v in vertices_of(&e.0) {
                    if let Some(a_le) = mve.get(&sk(&v.0)) {
                        for a_el in a_le {
                            if a_added.insert(sk(a_el)) {
                                let mut el = a_el.clone();
                                el.set_orientation(Orientation::Internal);
                                edges.push(Edge(el));
                                growing = true;
                            }
                        }
                    }
                }
            }
        }
        the_wires.push(TopoBuilder::new().make_wire(&edges).0);
    }
    the_wires
}

/// `BOPAlgo_BuilderFace::PerformInternalShapes` (`BOPAlgo_BuilderFace.cxx:618`).
pub fn perform_internal_shapes(fb: &mut FaceBuilder) -> Result<(), String> {
    if fb.base.avoid_internal_shapes {
        return Ok(());
    }
    if fb.base.internal.is_empty() {
        return Ok(());
    }
    let Some(gf) = fb.face.clone() else {
        return Ok(());
    };

    let mut a_box_tree = Box2dTree::new();
    let mut an_edges_map: Vec<TopoShape> = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();
    for w in &fb.base.internal {
        for e in edges_of(w) {
            if !seen.insert(sk(&e.0)) {
                continue;
            }
            let mut a_box_e = occt_core::bnd::BndBox2d::new();
            add_uv_bounds_edge(&gf, &e, &mut a_box_e);
            a_box_tree.add(an_edges_map.len() as i32, a_box_e);
            an_edges_map.push(e.0);
        }
    }
    a_box_tree.build();

    let mut a_me_done: HashSet<i32> = HashSet::new();
    let mut ctx = IntToolsContext::new();
    let areas = fb.base.areas.clone();
    for a_f in &areas {
        let mut face = Face(a_f.clone());
        let mut a_box_f = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&face, &mut a_box_f);
        let mut a_selector = Box2dTreeSelector::new();
        a_selector.set_bvh_set(&a_box_tree);
        a_selector.set_box(a_box_f);
        let a_li = a_selector.select();
        if a_li.is_empty() {
            continue;
        }
        let mut an_edges_inside: Vec<TopoShape> = Vec::new();
        for n_e in a_li {
            if a_me_done.contains(&n_e) {
                continue;
            }
            let idx = n_e as usize;
            if idx >= an_edges_map.len() {
                continue;
            }
            let a_e = &an_edges_map[idx];
            if is_inside(a_e, &face, &mut ctx) {
                an_edges_inside.push(a_e.clone());
                a_me_done.insert(n_e);
            }
        }
        if an_edges_inside.is_empty() {
            continue;
        }
        let a_lsi = make_internal_wires(&an_edges_inside);
        let b = TopoBuilder::new();
        for wi in a_lsi {
            b.add_wire(&mut face, &Wire(wi));
        }
        if let Some(pos) = fb.base.areas.iter().position(|a| a.same_tshape(a_f)) {
            fb.base.areas[pos] = face.0;
        }
        if a_me_done.len() == an_edges_map.len() {
            return Ok(());
        }
    }

    let unused: Vec<TopoShape> = an_edges_map
        .iter()
        .enumerate()
        .filter(|(i, _)| !a_me_done.contains(&(*i as i32)))
        .map(|(_, e)| e.clone())
        .collect();
    if !unused.is_empty() {
        let a_lsi = make_internal_wires(&unused);
        let b = TopoBuilder::new();
        let mut a_w_shape = b.make_compound_of(&[]);
        b.add_compound(&mut a_w_shape, &gf.0);
        if a_lsi.len() == 1 {
            b.add_compound(&mut a_w_shape, &a_lsi[0]);
        } else {
            let mut a_ce = b.make_compound_of(&[]);
            for w in &a_lsi {
                b.add_compound(&mut a_ce, w);
            }
            b.add_compound(&mut a_w_shape, &a_ce.0);
        }
        let _ = a_w_shape;
        fb.base.warnings.push(format!(
            "BOPAlgo_AlertFaceBuilderUnusedEdges: {} unused internal edge group(s)",
            a_lsi.len()
        ));
    }
    Ok(())
}

/// `BOPAlgo_BuilderFace::CheckData`.
pub fn check_data(fb: &FaceBuilder) -> Result<(), String> {
    if fb.face.is_none() {
        return Err("BOPAlgo_AlertNullInputShapes".into());
    }
    Ok(())
}

/// `BOPAlgo_BuilderFace::Perform`.
pub fn perform(fb: &mut FaceBuilder) -> Result<(), String> {
    fb.base.errors.clear();
    check_data(fb)?;
    perform_shapes_to_avoid(fb)?;
    perform_loops(fb)?;
    perform_areas(fb)?;
    perform_internal_shapes(fb)?;
    Ok(())
}

/// Keep `Vertex` referenced for MapShapesAndAncestors helpers.
pub fn vertex_key(v: &Vertex) -> usize {
    sk(&v.0)
}
