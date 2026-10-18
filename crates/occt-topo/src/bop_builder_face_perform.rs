//! Remainder of `BOPAlgo_BuilderFace.cxx`: `PerformShapesToAvoid`,
//! `PerformLoops`, `PerformAreas`, `PerformInternalShapes`, `IsGrowthWire`,
//! `IsInside`, `MakeInternalWires`.
//!
//! Source: `BOPAlgo_BuilderFace.cxx` (Perform at 117, ShapesToAvoid at 152,
//! Loops at 239, Areas at 387, InternalShapes at 618, MakeInternalWires at
//! 782, IsInside at 842).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

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
use crate::topo_tools::{edge_vertices, edges_of, is_equal, is_same, vertices_of};
use crate::wire_splitter::{WireEdgeSet, WireSplitter};

fn sk(s: &TopoShape) -> usize {
    GeometryRegistry::shape_key(s)
}

/// `TopExp::MapShapesAndAncestors(theShape, TopAbs_VERTEX, TopAbs_EDGE, theMap)`
/// (`TopExp.cxx:80-120`): the ancestor edge is appended once per *vertex
/// occurrence*, not once per distinct vertex. `TopExp_Explorer` walks with
/// `TopoDS_Iterator` and never uniquifies (`TopExp_Explorer.cxx:110-170`), so an
/// edge whose two children are the same vertex -- a closed edge -- is appended
/// twice. `BOPAlgo_BuilderFace::PerformShapesToAvoid` counts on that
/// (`BOPAlgo_BuilderFace.cxx:211-225`: `aNbE == 2` with `aE1.IsSame(aE2)` and
/// `aV1x.IsSame(aV2x)` -> `continue`, the closed edge is *not* avoided).
fn map_vertex_edges(edges: &[TopoShape]) -> HashMap<usize, Vec<TopoShape>> {
    let mut mve: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for e in edges {
        for v in stored_vertex_children(e) {
            mve.entry(sk(&v)).or_default().push(e.clone());
        }
    }
    mve
}

/// The stored child vertices of a shape, in stored order, without uniquifying
/// (`TopoDS_Iterator`).
fn stored_vertex_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
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
        for (_vk, le) in &mve {
            let nb_e = le.len();
            if nb_e == 0 {
                continue;
            }
            let e1 = &le[0];
            if nb_e == 1 {
                if BRepTool::is_degenerated(&Edge(e1.clone())) {
                    continue;
                }
                let internal_v = vertices_of(e1).iter().any(|v| v.0.orientation() == Orientation::Internal);
                if internal_v {
                    continue;
                }
                b_found = true;
                if !fb.base.avoid.iter().any(|a| is_equal(a, e1)) {
                    fb.base.avoid.push(e1.clone());
                }
            } else if nb_e == 2 {
                let e2 = &le[1];
                if is_same(e2, e1) {
                    let (v1x, v2x) = edge_vertices(&Edge(e1.clone()));
                    let same_ends = match (v1x, v2x) {
                        (Some(a), Some(b)) => is_same(&a.0, &b.0),
                        _ => false,
                    };
                    if same_ends {
                        continue;
                    }
                    b_found = true;
                    if !fb.base.avoid.iter().any(|a| is_equal(a, e1)) {
                        fb.base.avoid.push(e1.clone());
                    }
                    if !fb.base.avoid.iter().any(|a| is_equal(a, e2)) {
                        fb.base.avoid.push(e2.clone());
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
    let mut wes = WireEdgeSet::new();
    if let Some(face) = fb.face.clone() {
        wes.set_face(face);
    }
    for s in &fb.base.shapes {
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        if fb.base.is_avoided(s) {
            continue;
        }
        wes.add_edge(Edge(s.clone()));
    }
    if wes.is_empty() {
        return Ok(());
    }
    let mut wsp = WireSplitter::new();
    wsp.set_wes(wes);
    if let Err(e) = wsp.perform() {
        fb.base.warnings.push(e);
        return Ok(());
    }
    for w in wsp.wires() {
        fb.base.loops.push(w.clone());
    }
    let mut mep: HashSet<usize> = HashSet::new();
    for w in &fb.base.loops {
        for e in edges_of(w) {
            mep.insert(sk(&e.0));
        }
    }
    for a in &fb.base.avoid {
        mep.insert(sk(a));
    }
    for s in &fb.base.shapes {
        if s.shape_type() == ShapeType::Edge && !mep.contains(&sk(s)) {
            if !fb.base.avoid.iter().any(|a| is_equal(a, s)) {
                fb.base.avoid.push(s.clone());
            }
        }
    }
    fb.base.internal.clear();
    let nb_ea = fb.base.avoid.len();
    if nb_ea == 0 {
        return Ok(());
    }
    let mve = map_vertex_edges(&fb.base.avoid);
    let mut m_added: HashSet<usize> = HashSet::new();
    let mut b_flag = true;
    let avoid = fb.base.avoid.clone();
    for ee in &avoid {
        if !b_flag {
            break;
        }
        if !m_added.insert(sk(ee)) {
            continue;
        }
        let mut edges: Vec<Edge> = vec![Edge(ee.clone())];
        let mut growing = true;
        while growing && b_flag {
            growing = false;
            let current = edges.clone();
            for e in &current {
                for v in vertices_of(&e.0) {
                    if let Some(le) = mve.get(&sk(&v.0)) {
                        for ex in le {
                            if m_added.insert(sk(ex)) {
                                edges.push(Edge(ex.clone()));
                                growing = true;
                                if m_added.len() == nb_ea {
                                    b_flag = false;
                                }
                            }
                        }
                    }
                }
            }
        }
        let w_val = TopoBuilder::new().make_wire(&edges);
        fb.base.internal.push(w_val.0);
    }
    Ok(())
}

/// `IsGrowthWire` (`BOPAlgo_BuilderFace.cxx`).
pub fn is_growth_wire(wire: &TopoShape, mhe: &HashSet<usize>) -> bool {
    for e in edges_of(wire) {
        if mhe.contains(&sk(&e.0)) {
            return true;
        }
    }
    false
}

/// `IsInside` (`BOPAlgo_BuilderFace.cxx:842`).
///
/// Shared-edge rule: if the face already contains an edge of the wire, OCCT
/// **returns false immediately** (the wire cannot lie inside that face).
pub fn face_perform_is_inside(
    wire_or_edge: &TopoShape,
    f: &Face,
    ctx: &mut IntToolsContext,
) -> bool {
    let mut face_edges: HashSet<usize> = HashSet::new();
    for e in edges_of(&f.0) {
        face_edges.insert(sk(&e.0));
    }
    let edges: Vec<Edge> = if wire_or_edge.shape_type() == ShapeType::Edge {
        vec![Edge(wire_or_edge.clone())]
    } else {
        edges_of(wire_or_edge)
    };
    let Ok(clsf) = ctx.fclass2d(f) else {
        return false;
    };
    let mut is_inside = false;
    for e in &edges {
        if BRepTool::is_degenerated(e) {
            continue;
        }
        if face_edges.contains(&sk(&e.0)) {
            return is_inside;
        }
        let Some((c2d, t1, t2)) = crate::boptools_2d::curve_on_surface_range(e, f) else {
            continue;
        };
        let p2 = c2d.d0(0.5 * (t1 + t2));
        let st = clsf.perform(GpPnt2d::new(p2.x(), p2.y()));
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
    let tol = BRepTool::face_tolerance(gf);
    let mut geom = FaceGeom::new(surf);
    geom.tolerance = tol;
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
    let tol = BRepTool::face_tolerance(gf);
    let face = TopoBuilder::new().make_face(surf, &[wire.clone()]);
    GeometryRegistry::global().set_face_tolerance(&face.0, tol);
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
    let tol = BRepTool::face_tolerance(&gf);
    let mut ctx = IntToolsContext::new();

    if fb.base.loops.is_empty() {
        if ctx.is_infinite_face(&gf) {
            let natural = BRepTool::natural_restriction(&gf);
            let face_val = make_empty_face(fb, natural)?;
            fb.base.areas.push(face_val.0);
        }
        return Ok(());
    }

    let mut new_faces: Vec<TopoShape> = Vec::new();
    let mut hole_faces: Vec<TopoShape> = Vec::new();
    let mut mhe: HashSet<usize> = HashSet::new();

    for (_i_loop, wire_val) in fb.base.loops.clone().into_iter().enumerate() {
        let wire = Wire(wire_val.clone());
        let face_val = make_face_on(fb, &wire)?;
        let mut b_is_growth = is_growth_wire(&wire_val, &mhe);
        if !b_is_growth {
            if let Ok(c) = FClass2d::new(&face_val, tol) {
                b_is_growth = !c.is_hole();
            }
        }
        if b_is_growth {
            new_faces.push(face_val.0);
        } else {
            hole_faces.push(face_val.0);
            for e in edges_of(&wire_val) {
                mhe.insert(sk(&e.0));
            }
        }
    }

    if hole_faces.is_empty() {
        fb.base.areas.extend(new_faces);
        return Ok(());
    }

    let nb_h = hole_faces.len();
    let mut box_tree = Box2dTree::new();
    box_tree.set_size(nb_h);
    for (i, h_face) in hole_faces.iter().enumerate() {
        let mut bnd_box = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&Face(h_face.clone()), &mut bnd_box);
        box_tree.add(i as i32, bnd_box);
    }
    box_tree.build();

    let mut hole_face_map: HashMap<usize, TopoShape> = HashMap::new();
    let mut selector = Box2dTreeSelector::new();
    selector.set_bvh_set(&box_tree);

    for face_val in &new_faces {
        let af = Face(face_val.clone());
        let mut bnd_box = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&af, &mut bnd_box);
        selector.clear();
        selector.set_box(bnd_box);
        let li = selector.select();
        for k in li {
            let idx = k as usize;
            if idx >= hole_faces.len() {
                continue;
            }
            let hole_val = &hole_faces[idx];
            if !face_perform_is_inside(hole_val, &af, &mut ctx) {
                continue;
            }
            let hk = sk(hole_val);
            if let Some(was) = hole_face_map.get(&hk).cloned() {
                if face_perform_is_inside(face_val, &Face(was), &mut ctx) {
                    hole_face_map.insert(hk, face_val.clone());
                }
            } else {
                hole_face_map.insert(hk, face_val.clone());
            }
        }
    }

    let mut face_holes: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for (hk, face_val) in &hole_face_map {
        let Some(hole) = hole_faces.iter().find(|h| sk(h) == *hk) else {
            continue;
        };
        face_holes
            .entry(sk(face_val))
            .or_default()
            .push(hole.clone());
    }

    if hole_faces.len() != hole_face_map.len() {
        let box_f = shape_bbox(&gf.0);
        if box_f.is_open_xmin()
            || box_f.is_open_xmax()
            || box_f.is_open_ymin()
            || box_f.is_open_ymax()
            || box_f.is_open_zmin()
            || box_f.is_open_zmax()
        {
            let face_val = make_empty_face(fb, false)?;
            let mut unused: Vec<TopoShape> = Vec::new();
            for hole_val in &hole_faces {
                if !hole_face_map.contains_key(&sk(hole_val)) {
                    unused.push(hole_val.clone());
                }
            }
            face_holes.insert(sk(&face_val.0), unused);
            new_faces.push(face_val.0);
        }
    }

    let b = TopoBuilder::new();
    for face_val in new_faces {
        let holes = face_holes.get(&sk(&face_val)).cloned().unwrap_or_default();
        let mut face = Face(face_val);
        if !holes.is_empty() {
            for f_hole in &holes {
                if let Some(w) = first_wire_of(f_hole) {
                    b.add_wire(&mut face, &w);
                }
            }
            ctx.fclass2d_invalidate(&face);
            let _ = ctx.fclass2d(&face);
            let _ = tol;
        }
        fb.base.areas.push(face.0);
    }
    Ok(())
}

/// `MakeInternalWires`.
pub fn make_internal_wires(me: &[TopoShape]) -> Vec<TopoShape> {
    let mut wires: Vec<TopoShape> = Vec::new();
    let mve = map_vertex_edges(me);
    let mut added: HashSet<usize> = HashSet::new();
    for ee in me {
        if !added.insert(sk(ee)) {
            continue;
        }
        let mut edges: Vec<Edge> = Vec::new();
        let mut e0 = ee.clone();
        e0.set_orientation(Orientation::Internal);
        edges.push(Edge(e0));
        let mut growing = true;
        while growing {
            growing = false;
            let current = edges.clone();
            for e in &current {
                for v in vertices_of(&e.0) {
                    if let Some(le) = mve.get(&sk(&v.0)) {
                        for el_val in le {
                            if added.insert(sk(el_val)) {
                                let mut el = el_val.clone();
                                el.set_orientation(Orientation::Internal);
                                edges.push(Edge(el));
                                growing = true;
                            }
                        }
                    }
                }
            }
        }
        wires.push(TopoBuilder::new().make_wire(&edges).0);
    }
    wires
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

    let mut box_tree = Box2dTree::new();
    let mut an_edges_map: Vec<TopoShape> = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();
    for w in &fb.base.internal {
        for e in edges_of(w) {
            if !seen.insert(sk(&e.0)) {
                continue;
            }
            let mut box_e = occt_core::bnd::BndBox2d::new();
            add_uv_bounds_edge(&gf, &e, &mut box_e);
            box_tree.add(an_edges_map.len() as i32, box_e);
            an_edges_map.push(e.0);
        }
    }
    box_tree.build();

    let mut me_done: HashSet<i32> = HashSet::new();
    let mut ctx = IntToolsContext::new();
    let areas = fb.base.areas.clone();
    for f in &areas {
        let mut face = Face(f.clone());
        let mut box_f = occt_core::bnd::BndBox2d::new();
        add_uv_bounds_face(&face, &mut box_f);
        let mut selector = Box2dTreeSelector::new();
        selector.set_bvh_set(&box_tree);
        selector.set_box(box_f);
        let li = selector.select();
        if li.is_empty() {
            continue;
        }
        let mut an_edges_inside: Vec<TopoShape> = Vec::new();
        for n_e in li {
            if me_done.contains(&n_e) {
                continue;
            }
            let idx = n_e as usize;
            if idx >= an_edges_map.len() {
                continue;
            }
            let e_val = &an_edges_map[idx];
            if face_perform_is_inside(e_val, &face, &mut ctx) {
                an_edges_inside.push(e_val.clone());
                me_done.insert(n_e);
            }
        }
        if an_edges_inside.is_empty() {
            continue;
        }
        let lsi = make_internal_wires(&an_edges_inside);
        let b = TopoBuilder::new();
        for wi in lsi {
            b.add_wire(&mut face, &Wire(wi));
        }
        if let Some(pos) = fb.base.areas.iter().position(|a| a.same_tshape(f)) {
            fb.base.areas[pos] = face.0;
        }
        if me_done.len() == an_edges_map.len() {
            return Ok(());
        }
    }

    let unused: Vec<TopoShape> = an_edges_map
        .iter()
        .enumerate()
        .filter(|(i, _)| !me_done.contains(&(*i as i32)))
        .map(|(_, e)| e.clone())
        .collect();
    if !unused.is_empty() {
        let lsi = make_internal_wires(&unused);
        let b = TopoBuilder::new();
        let mut w_shape = b.make_compound_of(&[]);
        b.add_compound(&mut w_shape, &gf.0);
        if lsi.len() == 1 {
            b.add_compound(&mut w_shape, &lsi[0]);
        } else {
            let mut ce = b.make_compound_of(&[]);
            for w in &lsi {
                b.add_compound(&mut ce, w);
            }
            b.add_compound(&mut w_shape, &ce.0);
        }
        let _ = w_shape;
        fb.base.warnings.push(format!(
            "BOPAlgo_AlertFaceBuilderUnusedEdges: {} unused internal edge group(s)",
            lsi.len()
        ));
    }
    Ok(())
}

/// `BOPAlgo_BuilderFace::CheckData`.
pub fn check_face_data(fb: &FaceBuilder) -> Result<(), String> {
    if fb.face.is_none() {
        return Err("BOPAlgo_AlertNullInputShapes".into());
    }
    Ok(())
}

/// `BOPAlgo_BuilderFace::Perform`.
pub fn perform_face_builder(fb: &mut FaceBuilder) -> Result<(), String> {
    fb.base.errors.clear();
    check_face_data(fb)?;
    perform_shapes_to_avoid(fb)?;
    perform_loops(fb)?;
    perform_areas(fb)?;
    perform_internal_shapes(fb)?;
    Ok(())
}

/// Keep `Vertex` referenced for MapShapesAndAncestors helpers.
pub fn shape_vertex_key(v: &Vertex) -> usize {
    sk(&v.0)
}
