//! `BOPAlgo_BOP` — Fuse / Cut / Common after the General Fuse images exist.
//!
//! Source: `BOPAlgo_BOP.cxx`. The GF builder (`BopBuilder`) fills images and
//! a compound of every split (`BuildResult`). This module replaces that
//! compound with the Boolean result:
//!
//! * [`check_data`] — dimension rules (`CheckData`);
//! * [`build_rc`] — keep object splits that are (Common) / are not (Cut) in
//!   the tool split map, or all unique dim-pieces (Fuse);
//! * [`build_solid`] — Fuse of 3-D: unique outer faces → `BuilderSolid`;
//! * [`build_shape`] — open-solid `BuildBOP` first, else BuildRC + BuildSolid.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools_range::{dimension, dimensions};
use crate::bop_builder2::{BoolOp2, BopBuilder};
use crate::bop_occt_util::{
    explore, images_of, iter_children, shape_key, source_solids, treat_compound,
};
use crate::bop_tools_set::{BopToolsSet, BopToolsSetMap};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::builder_solid::BuilderSolid;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools_full::{edges_of, faces_of};

/// Shape type of the building element for a given dimension
/// (`TypeToExplore` in `BOPAlgo_BOP.cxx`).
pub fn type_to_explore(dim: i32) -> Option<ShapeType> {
    match dim {
        0 => Some(ShapeType::Vertex),
        1 => Some(ShapeType::Edge),
        2 => Some(ShapeType::Face),
        3 => Some(ShapeType::Solid),
        _ => None,
    }
}

/// `BOPTools_AlgoTools3D::IsEmptyShape` (`BOPTools_AlgoTools3D.cxx:732-741`):
/// `!HasGeometry` over the **whole sub-shape tree** (`Add`, `:745-786`).
///
/// The port's previous test ("no vertex **and** no face") counted *sub-shape
/// presence* rather than *geometry*, which is the opposite of OCCT in both
/// directions: OCCT reports a wire whose edge carries a 3D curve as **non-empty**
/// (the vertex rule is `:796-799`, the edge rule `:803-834`), and a face with no
/// surface and no triangulation as **empty** (`:838-851`).
pub fn is_empty_shape(s: &TopoShape) -> bool {
    let mut visited: HashSet<usize> = HashSet::new();
    !tree_has_geometry(s, &mut visited)
}

/// `Add` (`BOPTools_AlgoTools3D.cxx:745-786`): depth-first over
/// `TopoDS_Iterator(aSx, false, false)` — the same sub-shape tree, with the
/// visited map that keeps it linear.
fn tree_has_geometry(s: &TopoShape, visited: &mut HashSet<usize>) -> bool {
    let key = shape_key(s);
    if !visited.insert(key) {
        return false;
    }
    if has_geometry(s) {
        return true;
    }
    iter_children(s)
        .into_iter()
        .any(|c| tree_has_geometry(&c, visited))
}

/// `BOPAlgo_BOP::TreatEmptyShape` (`BOPAlgo_BOP.cxx:214-322`), called only when
/// `CheckData` raised `BOPAlgo_AlertEmptyShape` (`:162-167`).
///
/// `None` ⇒ OCCT returns `false` and the normal pipeline continues; `Some(v)` ⇒
/// OCCT adds `v` to `myShape` and returns `true`, i.e. the result *is* `v` (and is
/// an empty result when `v` is empty).
pub fn treat_empty_shape(
    objects: &[TopoShape],
    tools: &[TopoShape],
    op: BoolOp2,
) -> Option<Vec<TopoShape>> {
    // `:223-236`: find the non-empty objects and tools.
    let valid_obj: Vec<TopoShape> = objects
        .iter()
        .filter(|s| !is_empty_shape(s))
        .cloned()
        .collect();
    let valid_tool: Vec<TopoShape> = tools
        .iter()
        .filter(|s| !is_empty_shape(s))
        .cloned()
        .collect();
    let (has_obj, has_tool) = (!valid_obj.is_empty(), !valid_tool.is_empty());
    // `:240-243`: both groups hold valid shapes ⇒ continue the operation.
    if has_obj && has_tool {
        return None;
    }
    // `:245-250`: every shape is empty ⇒ the result is always an empty shape.
    if !has_obj && !has_tool {
        return Some(Vec::new());
    }
    // `:252-320`: exactly one group is all-empty, so the result can be built at
    // once — unless the operation would first have to split the survivors.
    Some(match op {
        BoolOp2::Fuse => {
            // `:267-272`: more than one valid shape must be split before adding.
            if valid_obj.len() + valid_tool.len() > 1 {
                return None;
            }
            if has_obj {
                valid_obj
            } else {
                valid_tool
            }
        }
        BoolOp2::Cut => {
            // `:280-285`: the objects must be split before adding.
            if valid_obj.len() > 1 {
                return None;
            }
            valid_obj
        }
        // `:305-307`: a Common with an empty group is always empty.
        BoolOp2::Common => Vec::new(),
    })
}

/// `HasGeometry` (`BOPTools_AlgoTools3D.cxx:790-854`).
fn has_geometry(s: &TopoShape) -> bool {
    let reg = crate::tgeometry::GeometryRegistry::global();
    match s.shape_type() {
        // `:796-799`: `TopAbs_VERTEX` always has geometry (its point).
        ShapeType::Vertex => true,
        // `:803-834`: any curve representation counts — `IsCurve3D()` with a
        // non-null curve, `IsCurveOnSurface()`, `IsRegularity()`, a non-null
        // `Polygon3D()`, `IsPolygonOnTriangulation()`, `IsPolygonOnSurface()`.
        // UNPORTED: the three polygon arms — the port's `EdgeGeom` keeps no edge
        // polygons (`tgeometry.rs`, `EdgeGeom`), so only curves/pcurves are seen.
        ShapeType::Edge => {
            let has_pcurve = reg
                .edge_geom(s)
                .map(|g| !g.pcurves.is_empty())
                .unwrap_or(false);
            has_pcurve || BRepTool::edge_curve(&Edge(s.clone())).is_some()
        }
        // `:838-851`: `Surface()` or `Triangulation()`.
        // UNPORTED: the triangulation arm — the port has no shape→triangulation
        // store reachable here, so only the surface is seen.
        ShapeType::Face => reg.face_geom(s).is_some(),
        // `:853`: wires / shells / solids / compounds carry no geometry of their
        // own — they are decided by their children through `Add`.
        _ => false,
    }
}

/// Min/max dimension of a group (`CheckData` loop).
fn group_dims(group: &[TopoShape]) -> (i32, i32, bool) {
    let mut d_min = 3;
    let mut d_max = 0;
    let mut has_valid = false;
    for s in group {
        if is_empty_shape(s) {
            continue;
        }
        let (a, b) = dimensions(s);
        if a < 0 {
            continue;
        }
        if a < d_min {
            d_min = a;
        }
        if b > d_max {
            d_max = b;
        }
        has_valid = true;
    }
    (d_min, d_max, has_valid)
}

/// `BOPAlgo_BOP::CheckData` dimension rules. Returns `(dim_objects, dim_tools)`.
pub fn check_data(objects: &[TopoShape], tools: &[TopoShape], op: BoolOp2) -> Result<(i32, i32), String> {
    if objects.is_empty() || tools.is_empty() {
        return Err("BOPAlgo_AlertTooFewArguments".into());
    }
    let (omin, omax, o_ok) = group_dims(objects);
    let (tmin, tmax, t_ok) = group_dims(tools);
    if op == BoolOp2::Fuse {
        if o_ok && omin != omax {
            return Err("BOPAlgo_AlertBOPNotAllowed".into());
        }
        if t_ok && tmin != tmax {
            return Err("BOPAlgo_AlertBOPNotAllowed".into());
        }
    }
    if o_ok && t_ok {
        let bad = match op {
            BoolOp2::Fuse => omax != tmax,
            BoolOp2::Cut => omax > tmin,
            BoolOp2::Common => false,
        };
        if bad {
            return Err("BOPAlgo_AlertBOPNotAllowed".into());
        }
    }
    let d0 = if o_ok { omin } else { tmin };
    let d1 = if t_ok { tmin } else { omin };
    Ok((d0, d1))
}

fn splits_of(b: &BopBuilder, s: &TopoShape) -> Vec<TopoShape> {
    match b.history().image(s) {
        Some(list) if !list.is_empty() => list.to_vec(),
        _ => vec![s.clone()],
    }
}

fn collect_building(group: &[TopoShape], _dim_filter: i32) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut fence: HashSet<usize> = HashSet::new();
    for s in group {
        let mut flat = Vec::new();
        treat_compound(s, &mut flat, &mut HashSet::new());
        for ss in flat {
            let d = dimension(&ss);
            let Some(ty) = type_to_explore(d) else {
                continue;
            };
            for sub in explore(&ss, ty) {
                if fence.insert(shape_key(&sub)) {
                    out.push(sub);
                }
            }
        }
    }
    out
}

fn map_splits(b: &BopBuilder, sources: &[TopoShape]) -> (HashSet<usize>, HashMap<BopToolsSet, TopoShape>) {
    let mut keys: HashSet<usize> = HashSet::new();
    let mut sets: HashMap<BopToolsSet, TopoShape> = HashMap::new();
    for s in sources {
        if s.shape_type() == ShapeType::Edge {
            let e = crate::shape::Edge(s.clone());
            if crate::brep_tool::BRepTool::is_degenerated(&e) {
                continue;
            }
        }
        for im in splits_of(b, s) {
            keys.insert(shape_key(&im));
            if im.shape_type() == ShapeType::Solid {
                let st = BopToolsSet::from_shape(&im, ShapeType::Face);
                sets.entry(st).or_insert_with(|| im.clone());
            }
        }
    }
    (keys, sets)
}

fn expand_lower(parts: &[TopoShape], dim_min: i32) -> HashSet<usize> {
    let mut out: HashSet<usize> = HashSet::new();
    for s in parts {
        let d_max = dimension(s);
        if d_max < 0 {
            continue;
        }
        let mut d = dim_min;
        while d < d_max {
            if let Some(ty) = type_to_explore(d) {
                for sub in explore(s, ty) {
                    out.insert(shape_key(&sub));
                }
            }
            d += 1;
        }
        out.insert(shape_key(s));
    }
    out
}

/// `BOPAlgo_BOP::BuildRC`.
pub fn build_rc(b: &BopBuilder, op: BoolOp2, dim0: i32, dim1: i32) -> TopoShape {
    let bb = TopoBuilder::new();
    if op == BoolOp2::Fuse {
        let Some(ty) = type_to_explore(dim0) else {
            return bb.make_compound_of(&[]).0;
        };
        let mut fence: HashSet<usize> = HashSet::new();
        let mut kept: Vec<TopoShape> = Vec::new();
        for s in explore(b.result(), ty) {
            if fence.insert(shape_key(&s)) {
                kept.push(s);
            }
        }
        return bb.make_compound_of(&kept).0;
    }

    let obj_src = collect_building(b.objects(), dim0);
    let tool_src = collect_building(b.tools(), dim1);
    let (obj_im, _obj_sets) = map_splits(b, &obj_src);
    let (_tool_im, tool_sets) = map_splits(b, &tool_src);

    let dim_min = dim0.min(dim1);
    let common = op == BoolOp2::Common;
    let it_keys: Vec<usize> = obj_im.iter().copied().collect();
    let it_shapes: Vec<TopoShape> = obj_src
        .iter()
        .flat_map(|s| splits_of(b, s))
        .collect();
    let check_keys = expand_lower(
        &tool_src.iter().flat_map(|s| splits_of(b, s)).collect::<Vec<_>>(),
        dim_min,
    );
    let it_exp = if common {
        expand_lower(&it_shapes, dim_min)
    } else {
        it_keys.into_iter().collect()
    };

    let mut kept: Vec<TopoShape> = Vec::new();
    let mut fence: HashSet<usize> = HashSet::new();
    for s in &it_shapes {
        let k = shape_key(s);
        if common && !it_exp.contains(&k) {
            continue;
        }
        let mut contains = check_keys.contains(&k);
        if !contains && s.shape_type() == ShapeType::Solid {
            let st = BopToolsSet::from_shape(s, ShapeType::Face);
            contains = tool_sets.contains_key(&st);
        }
        let take = if common { contains } else { !contains };
        if take && fence.insert(k) {
            kept.push(s.clone());
        }
    }
    bb.make_compound_of(&kept).0
}

/// `TopExp::MapShapesAndAncestors(sol, FACE, SOLID)` for the aMTSols scan.
fn map_faces_and_ancestors(sol: &TopoShape, mfs: &mut HashMap<usize, Vec<TopoShape>>) {
    for f in faces_of(sol) {
        mfs.entry(shape_key(&f.0)).or_default().push(sol.clone());
    }
}

/// One `theMFS` entry: first-inserted face (key orientation) and ancestor solids.
struct FaceToSolids {
    face: TopoShape,
    first_ori: Orientation,
    sols: Vec<TopoShape>,
}

/// `MapFacesToBuildSolids` (`BOPAlgo_BOP.cxx:1754`): skip INTERNAL; same-TShape
/// with opposite orientation appends a second solid; same orientation does not.
fn map_faces_to_build_solids(sol: &TopoShape, mfs: &mut HashMap<usize, FaceToSolids>) {
    for f in faces_of(sol) {
        if f.0.orientation() == Orientation::Internal {
            continue;
        }
        let k = shape_key(&f.0);
        match mfs.get_mut(&k) {
            None => {
                mfs.insert(
                    k,
                    FaceToSolids {
                        face: f.0.clone(),
                        first_ori: f.0.orientation(),
                        sols: vec![sol.clone()],
                    },
                );
            }
            Some(entry) => {
                if entry.first_ori != f.0.orientation() {
                    entry.sols.push(sol.clone());
                }
            }
        }
    }
}

/// `BOPAlgo_BOP::BuildSolid` (`BOPAlgo_BOP.cxx:1097`).
///
/// `CollectContainers` (`:1587`) gathers WIRE/SHELL/COMPSOLID arguments.
/// When that list is empty, the result is the compound of rebuilt solids
/// (`:1267`). CompSolid arguments rebuild connexity blocks of new solids
/// (`:1267-1377`). `BuilderSolid::SetContext` is unported.
pub fn build_solid(b: &BopBuilder, rc: &TopoShape) -> Result<TopoShape, String> {
    let bb = TopoBuilder::new();
    let mut arg_solids: HashSet<usize> = HashSet::new();
    let mut ancestors: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let mut a_lsc: Vec<TopoShape> = Vec::new();
    for group in [b.objects(), b.tools()] {
        for s in group {
            collect_containers(s, &mut a_lsc);
            for sol in explore(s, ShapeType::Solid) {
                arg_solids.insert(shape_key(&sol));
                map_faces_and_ancestors(&sol, &mut ancestors);
            }
        }
    }
    let mut shared: HashSet<usize> = HashSet::new();
    for sols in ancestors.values() {
        if sols.len() > 1 {
            for s in sols {
                shared.insert(shape_key(s));
            }
        }
    }

    let mut mfs: HashMap<usize, FaceToSolids> = HashMap::new();
    let mut mu_sols: Vec<TopoShape> = Vec::new();
    for sx in iter_children(rc) {
        if arg_solids.contains(&shape_key(&sx)) && !shared.contains(&shape_key(&sx)) {
            mu_sols.push(sx);
            continue;
        }
        map_faces_to_build_solids(&sx, &mut mfs);
    }

    let mut dmsts = BopToolsSetMap::new();
    for sx in mu_sols {
        let mapped = faces_of(&sx)
            .iter()
            .any(|f| mfs.contains_key(&shape_key(&f.0)));
        if mapped {
            map_faces_to_build_solids(&sx, &mut mfs);
        } else {
            let st = BopToolsSet::from_shape(&sx, ShapeType::Face);
            if !dmsts.contains(&st) {
                dmsts.add(st);
            }
        }
    }

    let mut sfs: Vec<TopoShape> = Vec::new();
    // `BOPAlgo_BOP::BuildSolid` walks the R/C faces in the data-structure order
    // and looks `theMFS` up per face, so the builder always receives them in
    // that order. Iterating `mfs.values()` instead hands `BuilderSolid` a
    // per-process `HashMap` order, and its growth/hole walk is order sensitive
    // (`fuse_box_cylinder_is_closed_solid` produced 0 areas in ~40% of runs).
    for f in faces_of(rc) {
        if let Some(entry) = mfs.get(&shape_key(&f.0)) {
            if entry.sols.len() == 1 {
                sfs.push(entry.face.clone());
            }
        }
    }
    let mut out: Vec<TopoShape> = Vec::new();
    if !sfs.is_empty() {
        let mut bs = BuilderSolid::new();
        bs.set_shapes(sfs);
        bs.set_fuzzy(b.fuzzy_value());
        bs.set_avoid_internal_shapes(true);
        bs.perform()?;
        out.extend(bs.areas().iter().cloned());
    }
    for set in dmsts.sets() {
        out.push(set.shape().clone());
    }
    let rc = bb.make_compound_of(&out).0;
    if a_lsc.is_empty() {
        return Ok(rc);
    }
    Ok(rebuild_compsolids(b, &rc, &a_lsc))
}

/// `CollectContainers` (`BOPAlgo_BOP.cxx:1587`).
fn collect_containers(s: &TopoShape, out: &mut Vec<TopoShape>) {
    match s.shape_type() {
        ShapeType::Wire | ShapeType::Shell | ShapeType::CompSolid => out.push(s.clone()),
        ShapeType::Compound => {
            for c in iter_children(s) {
                collect_containers(&c, out);
            }
        }
        _ => {}
    }
}

/// Compsolid rebuild after `BuildSolid` (`BOPAlgo_BOP.cxx:1267-1377`).
fn rebuild_compsolids(b: &BopBuilder, rc: &TopoShape, lsc: &[TopoShape]) -> TopoShape {
    let bb = TopoBuilder::new();
    let solids: Vec<TopoShape> = iter_children(rc)
        .into_iter()
        .filter(|s| s.shape_type() == ShapeType::Solid)
        .collect();
    if solids.is_empty() {
        return rc.clone();
    }
    if solids.len() == 1 {
        let mut cs = TopoShape::new(ShapeType::CompSolid);
        bb.add(&mut cs, &solids[0]);
        return bb.make_compound_of(&[cs]).0;
    }
    let mut mfcs: HashSet<usize> = HashSet::new();
    for cs in lsc {
        if cs.shape_type() != ShapeType::CompSolid {
            continue;
        }
        for f in explore(cs, ShapeType::Face) {
            let imgs = images_of(b.history(), &f);
            if imgs.is_empty() {
                mfcs.insert(shape_key(&f));
            } else {
                for im in imgs {
                    mfcs.insert(shape_key(im));
                }
            }
        }
    }
    let mut result_parts: Vec<TopoShape> = Vec::new();
    for block in connexity_solids_by_face(&solids) {
        let has_cs_face = block.iter().any(|sol| {
            faces_of(sol)
                .iter()
                .any(|f| mfcs.contains(&shape_key(&f.0)))
        });
        if !has_cs_face {
            result_parts.extend(block);
            continue;
        }
        let mut cs = TopoShape::new(ShapeType::CompSolid);
        for sol in &block {
            bb.add(&mut cs, sol);
        }
        result_parts.push(cs);
    }
    bb.make_compound_of(&result_parts).0
}

fn connexity_solids_by_face(solids: &[TopoShape]) -> Vec<Vec<TopoShape>> {
    let n = solids.len();
    let face_sets: Vec<HashSet<usize>> = solids
        .iter()
        .map(|s| faces_of(s).iter().map(|f| shape_key(&f.0)).collect())
        .collect();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] != x {
            p[x] = find(p, p[x]);
        }
        p[x]
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let (ra, rb) = (find(p, a), find(p, b));
        if ra != rb {
            p[ra] = rb;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if face_sets[i].intersection(&face_sets[j]).next().is_some() {
                union(&mut parent, i, j);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(solids[i].clone());
    }
    groups.into_values().collect()
}

fn edge_closed_on_face(edge: &TopoShape, face: &TopoShape) -> bool {
    edges_of(face)
        .into_iter()
        .filter(|e| e.0.same_tshape(edge))
        .count()
        > 1
}

fn edge_internal_on_face(edge: &TopoShape, face: &TopoShape) -> bool {
    edges_of(face).into_iter().any(|e| {
        e.0.same_tshape(edge) && e.0.orientation() == Orientation::Internal
    })
}

/// `BOPAlgo_BOP::CheckArgsForOpenSolid` (`BOPAlgo_BOP.cxx:1382`).
///
/// FaceBuilder unused-internal-edge warnings are not this check. Only a
/// `BOPAlgo_AlertSolidBuilderUnusedFaces` analogue plus a free-edge test on
/// the source solids (or new INTERNAL faces on their images) may divert to
/// the open-solid `BuildBOP` rebuild.
fn check_args_for_open_solid(b: &BopBuilder) -> bool {
    let failed_from_warning = b.warnings().iter().any(|w| {
        w.contains("AlertSolidBuilderUnusedFaces")
            || (w.contains("BuilderSolid:") && w.contains("fit no region"))
    });
    for a_solid in source_solids(b.ds()) {
        let mut mef: HashMap<usize, (TopoShape, Vec<TopoShape>)> = HashMap::new();
        let mut mf_internal: HashSet<usize> = HashSet::new();
        for a_sh in iter_children(&a_solid) {
            if a_sh.shape_type() != ShapeType::Shell {
                continue;
            }
            for a_f in iter_children(&a_sh) {
                if a_f.shape_type() != ShapeType::Face {
                    continue;
                }
                if a_f.orientation() == Orientation::Internal {
                    for im in images_of(b.history(), &a_f) {
                        mf_internal.insert(shape_key(im));
                    }
                    if images_of(b.history(), &a_f).is_empty() {
                        mf_internal.insert(shape_key(&a_f));
                    }
                    continue;
                }
                for e in edges_of(&a_f) {
                    let ent = mef.entry(shape_key(&e.0)).or_insert_with(|| (e.0.clone(), Vec::new()));
                    if !ent.1.iter().any(|f| f.same_tshape(&a_f)) {
                        ent.1.push(a_f.clone());
                    }
                }
            }
        }
        let mut is_closed = true;
        for (edge, faces) in mef.values() {
            if BRepTool::is_degenerated(&Edge(edge.clone())) {
                continue;
            }
            if faces.len() > 1 {
                continue;
            }
            let Some(face) = faces.first() else {
                continue;
            };
            is_closed = edge_closed_on_face(edge, face) || edge_internal_on_face(edge, face);
            if !is_closed {
                break;
            }
        }
        if is_closed {
            continue;
        }
        if failed_from_warning {
            return true;
        }
        let images = images_of(b.history(), &a_solid);
        if images.is_empty() {
            continue;
        }
        for a_sim in images {
            for a_sh in iter_children(a_sim) {
                if a_sh.shape_type() != ShapeType::Shell {
                    continue;
                }
                for a_f in iter_children(&a_sh) {
                    if a_f.shape_type() == ShapeType::Face
                        && a_f.orientation() == Orientation::Internal
                        && !mf_internal.contains(&shape_key(&a_f))
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// `BOPAlgo_BOP::BuildShape` — the only Fuse/Cut/Common filter after GF
/// FillImages + BuildResult (`BOP.cxx:871`).
pub fn build_shape(b: &mut BopBuilder) -> Result<(), String> {
    let op = match (b.obj_state(), b.tools_state()) {
        (FaceState::Out, FaceState::Out) => BoolOp2::Fuse,
        (FaceState::Out, FaceState::In) => BoolOp2::Cut,
        (FaceState::In, FaceState::In) => BoolOp2::Common,
        _ => BoolOp2::Cut,
    };
    let (dim0, dim1) = check_data(b.objects(), b.tools(), op)?;
    let open = dim0 == 3 && dim1 == 3 && check_args_for_open_solid(b);
    if open {
        if crate::bop_build_bop::build_bop(b).is_ok() {
            return Ok(());
        }
    }
    let rc = build_rc(b, op, dim0, dim1);
    if op == BoolOp2::Fuse && dim0 == 3 {
        let shape = build_solid(b, &rc)?;
        b.set_result_shape(shape);
        return Ok(());
    }
    b.set_result_shape(rc);
    Ok(())
}
