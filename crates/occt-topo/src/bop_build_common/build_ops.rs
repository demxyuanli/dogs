use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Host contract
// ---------------------------------------------------------------------------

/// The minimal surface the building-common operations need from the builder
/// main class (`BOPAlgo_Builder` / `crate::bop_builder2::BopBuilder`).
///
/// `history` holds the images naming table (old shape → its splits); `origins`
/// is the back-map (split → the originals it came from). `ds` is the
/// intersection data structure. `arguments` feeds [`fill_internal_shapes`];
/// `fuzzy_value` is the additional tolerance for the 2-D vertex classification
/// in [`fill_internal_vertices`].

pub trait BopBuildOps {
    /// The data structure of the algorithm.
    fn ds(&self) -> &BopdsDS;
    /// The images history (naming side tables).
    fn history(&self) -> &BopHistory;
    /// Mutable images history.
    fn history_mut(&mut self) -> &mut BopHistory;
    /// The additional tolerance of the operation.
    fn fuzzy_value(&self) -> f64;
    /// The arguments of the operation.
    fn arguments(&self) -> &[TopoShape];
    /// The origins back-map (split shape → original shapes).
    fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>>;
}

// ---------------------------------------------------------------------------
// Shape-set helpers
// ---------------------------------------------------------------------------

/// Stable identity key of a shape (the address of its shared `TShape`).
pub(super) fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Add `s` to the ordered set `set`; returns true when it was newly inserted.
pub(super) fn set_add(set: &mut Vec<TopoShape>, s: TopoShape) -> bool {
    if set.iter().any(|x| x.same_tshape(&s)) {
        false
    } else {
        set.push(s);
        true
    }
}

/// True when `set` already holds a shape identical to `s`.
pub(super) fn set_contains(set: &[TopoShape], s: &TopoShape) -> bool {
    set.iter().any(|x| x.same_tshape(s))
}

/// Direct structural children of `shape`.
pub(super) fn direct_children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .clone()
}

/// Reverse the orientation of a shape view.
pub(super) fn reverse_orientation(s: &mut TopoShape) {
    let o = s.orientation();
    s.set_orientation(o.reversed());
}

// ---------------------------------------------------------------------------
// Orientation of a split relative to its original
// ---------------------------------------------------------------------------

/// Whether `split` must be reversed to match the direction of `original`.
///
/// Mirrors `BOPTools_AlgoTools::IsSplitToReverse`: when both shapes share the
/// same underlying curve/surface the orientations are compared directly;
/// otherwise the tangent/normal directions of the two geometries are compared
/// at a common point. `false` is returned for shapes without registered
/// geometry or of a non face/edge type.
pub(crate) fn is_split_to_reverse(split: &TopoShape, original: &TopoShape) -> bool {
    if split.shape_type() != original.shape_type() {
        return false;
    }
    match split.shape_type() {
        ShapeType::Face => face_split_to_reverse(&Face(split.clone()), &Face(original.clone())),
        ShapeType::Edge => edge_split_to_reverse(&Edge(split.clone()), &Edge(original.clone())),
        _ => false,
    }
}

/// A 3-D point on `face`: the surface centre when its UV range is bounded, else
/// the midpoint of the first boundary edge.
pub(super) fn face_sample_point(face: &Face) -> Option<GpPnt> {
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && v1.is_finite() {
        let s = BRepTool::face_surface(face)?;
        return Some(s.d0(0.5 * (u1 + u2), 0.5 * (v1 + v2)));
    }
    let e = edges_of(&face.0).into_iter().next()?;
    let (a, b) = BRepTool::edge_parameters(&e);
    if a.is_finite() && b.is_finite() {
        let c = BRepTool::edge_curve(&e)?;
        Some(c.d0(0.5 * (a + b)))
    } else {
        None
    }
}

/// Face variant: compare the surface normals at a point of the split face.
pub(super) fn face_split_to_reverse(split: &Face, orig: &Face) -> bool {
    let (Some(s), Some(o)) = (BRepTool::face_surface(split), BRepTool::face_surface(orig)) else {
        return false;
    };
    if Arc::ptr_eq(&s, &o) {
        return split.orientation() != orig.orientation();
    }
    let Some(p) = face_sample_point(split) else { return false };
    let (uo, vo) = surface_closest_params(o.as_ref(), &p, 32, 32);
    if !uo.is_finite() || !vo.is_finite() {
        return false;
    }
    // Normal of the split face at `p` (project back onto the split surface).
    let (u, v) = surface_closest_params(s.as_ref(), &p, 32, 32);
    let mut ns = surface_normal_checked(s.as_ref(), u, v);
    let mut no = surface_normal_checked(o.as_ref(), uo, vo);
    if ns.square_magnitude() < 1e-30 || no.square_magnitude() < 1e-30 {
        return false;
    }
    if split.orientation() == Orientation::Reversed {
        ns = ns.reversed();
    }
    if orig.orientation() == Orientation::Reversed {
        no = no.reversed();
    }
    ns.dot(&no) < 0.0
}

/// Edge variant: compare the tangent vectors at a point of the split edge.
pub(super) fn edge_split_to_reverse(split: &Edge, orig: &Edge) -> bool {
    let (Some(cs), Some(co)) = (BRepTool::edge_curve(split), BRepTool::edge_curve(orig)) else {
        return false;
    };
    if Arc::ptr_eq(&cs, &co) {
        return split.orientation() != orig.orientation();
    }
    let (a, b) = BRepTool::edge_parameters(split);
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    let t = 0.5 * (a + b);
    let p = cs.d0(t);
    let (to, _) = closest_point_on_edge(orig, &p, 32);
    let ts = cs.d1(t).1;
    let to = co.d1(to).1;
    if ts.square_magnitude() < 1e-30 || to.square_magnitude() < 1e-30 {
        return false;
    }
    ts.dot(&to) < 0.0
}

/// Surface normal, falling back to the zero vector on error.
pub(super) fn surface_normal_checked(s: &dyn occt_geom::Surface, u: f64, v: f64) -> occt_core::gp::GpVec {
    let n = surface_closest_normal(s, u, v);
    if n.square_magnitude() < 1e-30 {
        occt_core::gp::GpVec::zero()
    } else {
        n
    }
}

/// Normal of a surface at `(u, v)` via the finite-difference helper.
pub(super) fn surface_closest_normal(s: &dyn occt_geom::Surface, u: f64, v: f64) -> occt_core::gp::GpVec {
    let du = 1e-6;
    let p0 = s.d0(u, v);
    let pu = s.d0(u + du, v);
    let pv = s.d0(u, v + du);
    let a = occt_core::gp::GpVec::from_pnts(&p0, &pu);
    let b = occt_core::gp::GpVec::from_pnts(&p0, &pv);
    a.crossed(&b)
}

// ---------------------------------------------------------------------------
// Containers (wires / shells / compsolids)
// ---------------------------------------------------------------------------

/// Fills the images of the container shapes of the data structure.
///
/// Mirrors `BOPAlgo_Builder::FillImagesContainers(theType)`: every source
/// container of `container_type` (wire, shell or comp-solid) is rebuilt from
/// the splits of its direct sub-shapes. A container whose sub-shapes were not
/// modified keeps no image (it is returned as-is by the caller). The caller
/// invokes this once per container type at the matching stage of the build
/// sequence (WIRE after the edges, SHELL after the faces, COMPSOLID after the
/// solids), so each type's direct sub-shape images are already filled.
pub fn fill_images_containers<B: BopBuildOps>(
    f: &mut B,
    container_type: ShapeType,
) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() == container_type {
            let c = si.shape().clone();
            fill_images_container(f, &c, container_type)?;
        }
    }
    Ok(())
}

/// Builds the image of a single container (`BOPAlgo_Builder::FillImagesContainer`).
///
/// A wire is rebuilt with [`TopoBuilder::make_wire`] over the edge splits; a
/// shell is reassembled with [`ShellSplitter`] so the produced image shells are
/// closed. When the face splits do not form any closed shell (an open source
/// shell), a plain shell holding all splits is produced instead, matching the
/// OCCT behaviour of always emitting a container image.
pub(super) fn fill_images_container<B: BopBuildOps>(
    f: &mut B,
    container: &TopoShape,
    container_type: ShapeType,
) -> Result<(), String> {
    let children = direct_children(container);
    // Check if any direct sub-shape carries a non-identity image.
    let mut modified = false;
    for ss in &children {
        if let Some(im) = f.history().image(ss) {
            if im.len() != 1 || !im[0].same_tshape(ss) {
                modified = true;
                break;
            }
        }
    }
    if !modified {
        return Ok(());
    }

    let bld = TopoBuilder::new();
    match container_type {
        ShapeType::Wire => {
            let mut edges: Vec<Edge> = Vec::new();
            for ss in &children {
                if ss.shape_type() != ShapeType::Edge {
                    continue;
                }
                collect_splits(f, ss, &mut |im| {
                    edges.push(Edge(im));
                });
            }
            let wire = bld.make_wire(&edges);
            wire.0.set_closed(crate::topo_tools_full::wire_is_closed(&wire));
            f.history_mut().add_image(container, wire.0);
        }
        ShapeType::Shell | ShapeType::CompSolid => {
            let mut faces: Vec<TopoShape> = Vec::new();
            for ss in &children {
                if ss.shape_type() != ShapeType::Face {
                    continue;
                }
                collect_splits(f, ss, &mut |im| faces.push(im));
            }
            if faces.is_empty() {
                return Ok(());
            }
            // Reassemble the faces into closed shells.
            let mut splitter = ShellSplitter::new();
            for fc in &faces {
                splitter.add_start_element(fc.clone());
            }
            splitter.perform()?;
            let shells = splitter.shells().to_vec();
            if shells.is_empty() {
                // No closed shell (open source shell): emit a plain shell.
                let shell_faces: Vec<Face> = faces.iter().map(|s| Face(s.clone())).collect();
                let shell = bld.make_shell(&shell_faces);
                shell.0.set_closed(!AlgoTools::is_open_shell(&shell.0));
                f.history_mut().add_image(container, shell.0);
            } else {
                for s in shells {
                    f.history_mut().add_image(container, s);
                }
            }
        }
        _ => {}
    }
    Ok(())
}

/// Pushes the image splits of `sub` into `out`, reversing those whose
/// orientation is inverted relative to `sub`. When `sub` has no image, it is
/// pushed itself (`BOPAlgo_Builder::FillImagesContainer` inner loop).
pub(super) fn collect_splits<B: BopBuildOps, F: FnMut(TopoShape)>(f: &B, sub: &TopoShape, out: &mut F) {
    if let Some(im) = f.history().image(sub) {
        for ims in im {
            let mut s = ims.clone();
            if is_split_to_reverse(&s, sub) {
                reverse_orientation(&mut s);
            }
            out(s);
        }
    } else {
        out(sub.clone());
    }
}

// ---------------------------------------------------------------------------
// Compounds
// ---------------------------------------------------------------------------

/// Fills the images of the compound shapes of the data structure.
///
/// Mirrors `BOPAlgo_Builder::FillImagesCompounds`: every source compound is
/// recursively rebuilt from the images of its sub-shapes.
pub fn fill_images_compounds<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    let mut fence: Vec<TopoShape> = Vec::new();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() == ShapeType::Compound {
            let c = si.shape().clone();
            fill_images_compound(f, &c, &mut fence)?;
        }
    }
    Ok(())
}

/// Builds the image of a single compound (`BOPAlgo_Builder::FillImagesCompound`).
///
/// The compound is rebuilt only when at least one of its sub-shapes carries an
/// image (directly or through a nested compound). Each image is oriented as its
/// original sub-shape.
pub(super) fn fill_images_compound<B: BopBuildOps>(
    f: &mut B,
    compound: &TopoShape,
    fence: &mut Vec<TopoShape>,
) -> Result<(), String> {
    if !set_add(fence, compound.clone()) {
        return Ok(());
    }
    let children = direct_children(compound);
    let mut interfered = false;
    for sx in &children {
        if sx.shape_type() == ShapeType::Compound {
            fill_images_compound(f, sx, fence)?;
        }
        if f.history().has_image(sx) {
            interfered = true;
        }
    }
    if !interfered {
        return Ok(());
    }

    let bld = TopoBuilder::new();
    let mut c_im = bld.make_shape(ShapeType::Compound);
    for sx in &children {
        let or = sx.orientation();
        if f.history().has_image(sx) {
            let im = f.history().image(sx).unwrap();
            for ims in im {
                let mut s = ims.clone();
                s.set_orientation(or);
                bld.add(&mut c_im, &s);
            }
        } else {
            bld.add(&mut c_im, sx);
        }
    }
    f.history_mut().add_image(compound, c_im);
    Ok(())
}

// ---------------------------------------------------------------------------
// Internal vertices on faces
// ---------------------------------------------------------------------------

/// Classifies the alone vertices of every split face and adds those falling
/// inside a face image as `INTERNAL` vertices of the image.
///
/// Mirrors `BOPAlgo_Builder::FillInternalVertices`. An *alone* vertex of a face
/// is a vertex belonging to the face (its sub-shape set) that does not belong
/// to any boundary edge of the face.
pub fn fill_internal_vertices<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    // Collect (vertex oriented INTERNAL, face image) pairs to classify.
    let mut tasks: Vec<(Vertex, Face)> = Vec::new();
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Face {
            continue;
        }
        let face = si.shape().clone();
        let Some(images) = f.history().image(&face) else { continue };
        if images.is_empty() {
            continue;
        }
        let alone = alone_vertices(f.ds(), i);
        if alone.is_empty() {
            continue;
        }
        for v_idx in alone {
            let Some(v_shape) = f.ds().shape(v_idx).cloned() else { continue };
            let mut v = Vertex(v_shape);
            v.0.set_orientation(Orientation::Internal);
            for im in images {
                tasks.push((v.clone(), Face(im.clone())));
            }
        }
    }

    let fuzzy = f.fuzzy_value();
    let mut ctx = IntToolsContext::new();
    let bld = TopoBuilder::new();
    for (v, face_im) in tasks {
        // `ComputeVF` returns 0 when the vertex is strictly inside the face.
        if ctx.compute_vf(&v, &face_im, fuzzy) == 0 {
            let mut f_im = face_im.0;
            bld.add(&mut f_im, &v.0);
        }
    }
    Ok(())
}

/// DS indices of the alone vertices of the face with index `face_idx`
/// (`BOPDS_DS::AloneVertices`): vertices of the face that belong to no edge of
/// the face.
pub(super) fn alone_vertices(ds: &BopdsDS, face_idx: usize) -> Vec<usize> {
    let Some(si) = ds.shape_info(face_idx) else { return Vec::new() };
    let mut edge_vertices: HashSet<usize> = HashSet::new();
    let mut vertex_indices: Vec<usize> = Vec::new();
    for &sub in &si.sub_indices {
        let Some(sub_info) = ds.shape_info(sub) else { continue };
        match sub_info.shape_type() {
            ShapeType::Edge => {
                for &v in &sub_info.sub_indices {
                    edge_vertices.insert(v);
                }
            }
            ShapeType::Vertex => vertex_indices.push(sub),
            _ => {}
        }
    }
    vertex_indices
        .into_iter()
        .filter(|v| !edge_vertices.contains(v))
        .collect()
}

// ---------------------------------------------------------------------------
// Internal shapes on solids
// ---------------------------------------------------------------------------

/// Settles the vertices/edges/wires of the arguments (and those inside the
/// source solids) that are located inside a split solid as `INTERNAL` children
/// of that solid.
///
/// Mirrors `BOPAlgo_Builder::FillInternalShapes`:
/// 1. collect the candidate shapes (argument vertices/edges/wires and the
///    internal sub-shapes of the source solids), mapped through their images;
/// 2. build the vertex→edge / vertex→face / edge→face ancestor map of the split
///    solids and drop the candidates already tied to a face;
/// 3. for each candidate, classify one of its points against each split solid
///    and add it as `INTERNAL`; a candidate settled into an *original* (not
///    yet split) solid first copies the solid so the input shape is preserved.
pub fn fill_internal_shapes<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    // 1.1 — shapes from the pure arguments.
    let mut a_lsc: Vec<TopoShape> = Vec::new();
    let mut a_fence: Vec<TopoShape> = Vec::new();
    for a in f.arguments() {
        treat_compound(a, &mut a_lsc, &mut a_fence);
    }
    let mut a_largs: Vec<TopoShape> = Vec::new();
    a_fence.clear();
    for s in &a_lsc {
        match s.shape_type() {
            ShapeType::Wire => {
                for e in edges_of_wire(&Wire(s.clone())) {
                    if set_add(&mut a_fence, e.0.clone()) {
                        a_largs.push(e.0);
                    }
                }
            }
            ShapeType::Vertex | ShapeType::Edge => a_largs.push(s.clone()),
            _ => {}
        }
    }
    a_fence.clear();
    let mut a_msi: Vec<TopoShape> = Vec::new();
    for s in &a_largs {
        if !set_add(&mut a_fence, s.clone()) {
            continue;
        }
        let t = s.shape_type();
        if t == ShapeType::Vertex || t == ShapeType::Edge || t == ShapeType::Wire {
            if f.history().has_image(s) {
                for im in f.history().image(s).unwrap() {
                    set_add(&mut a_msi, im.clone());
                }
            } else {
                set_add(&mut a_msi, s.clone());
            }
        }
    }

    // 2. — internal vertices/edges from the source solids + ancestor map.
    a_fence.clear();
    let mut a_lsd: Vec<TopoShape> = Vec::new();
    let mut a_msx: HashMap<usize, Vec<TopoShape>> = HashMap::new();
    let mut a_msor: Vec<TopoShape> = Vec::new();
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let solid = si.shape().clone();
        // Own internal (non-shell) sub-shapes.
        for s in own_internal_shapes(&solid) {
            if f.history().has_image(&s) {
                for im in f.history().image(&s).unwrap() {
                    set_add(&mut a_msi, im.clone());
                }
            } else {
                set_add(&mut a_msi, s);
            }
        }
        // Ancestors of the splits (or of the solid itself when unsplit).
        if f.history().has_image(&solid) {
            for sp in f.history().image(&solid).unwrap() {
                if set_add(&mut a_fence, sp.clone()) {
                    map_ancestors_into(sp, &mut a_msx);
                    a_lsd.push(sp.clone());
                }
            }
        } else if set_add(&mut a_fence, solid.clone()) {
            map_ancestors_into(&solid, &mut a_msx);
            a_lsd.push(solid.clone());
            set_add(&mut a_msor, solid);
        }
    }

    // 3. — keep only the candidates not already tied to a face of the solids.
    let mut a_lsi: Vec<TopoShape> = Vec::new();
    for s in &a_msi {
        let tied = a_msx.contains_key(&shape_key(s));
        if !tied {
            a_lsi.push(s.clone());
        }
    }

    // 4. — nothing to settle.
    if a_lsi.is_empty() {
        return Ok(());
    }

    // 5. — settle the candidates into the split solids.
    let bld = TopoBuilder::new();
    for sd in &mut a_lsd {
        let mut i = 0;
        while i < a_lsi.len() {
            let mut si_shape = a_lsi[i].clone();
            si_shape.set_orientation(Orientation::Internal);
            let state =
                compute_state_by_one_point(&si_shape, sd, 1e-11).unwrap_or(FaceState::Unknown);
            if state != FaceState::In {
                i += 1;
                continue;
            }
            if set_contains(&a_msor, sd) {
                // Make a copy of the original solid so the input stays intact.
                let mut sdx = Solid::new();
                for sh in direct_children(sd) {
                    bld.add(&mut sdx.0, &sh);
                }
                bld.add(&mut sdx.0, &si_shape);
                f.history_mut().add_image(sd, sdx.0.clone());
                f.origins_mut()
                    .entry(shape_key(&sdx.0))
                    .or_default()
                    .push(sd.clone());
                a_msor.retain(|x| !x.same_tshape(sd));
                *sd = sdx.0;
            } else {
                bld.add(sd, &si_shape);
            }
            a_lsi.remove(i);
        }
    }
    Ok(())
}

/// Flattens `s` into `out`: a compound is expanded into its children, any other
/// shape is appended (once, through the fence) (`BOPTools_AlgoTools::TreatCompound`).
pub(super) fn treat_compound(s: &TopoShape, out: &mut Vec<TopoShape>, fence: &mut Vec<TopoShape>) {
    if s.shape_type() != ShapeType::Compound {
        if set_add(fence, s.clone()) {
            out.push(s.clone());
        }
        return;
    }
    for h in s.tshape.read().unwrap().children.clone() {
        treat_compound(&h, out, fence);
    }
}

/// The non-shell direct children of a solid: its internal vertices/edges/wires
/// (`BOPAlgo_Builder::OwnInternalShapes`).
pub(super) fn own_internal_shapes(solid: &TopoShape) -> Vec<TopoShape> {
    direct_children(solid)
        .into_iter()
        .filter(|c| c.shape_type() != ShapeType::Shell)
        .collect()
}

/// Merges the vertex→edge, vertex→face and edge→face ancestor relations of
/// `shape` into `out` (`TopExp::MapShapesAndAncestors`).
pub(super) fn map_ancestors_into(shape: &TopoShape, out: &mut HashMap<usize, Vec<TopoShape>>) {
    // vertex -> edges
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Edge, out);
    // vertex -> faces
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Face, out);
    // edge -> faces
    merge_ancestors(shape, ShapeType::Edge, ShapeType::Face, out);
}

/// For every ancestor of type `ancestor_type` in `shape`, records each of its
/// boundary sub-shapes of type `child_type` → the ancestor.
pub(super) fn merge_ancestors(
    shape: &TopoShape,
    child_type: ShapeType,
    ancestor_type: ShapeType,
    out: &mut HashMap<usize, Vec<TopoShape>>,
) {
    let ancestors = crate::topo_tools_full::shapes_of(shape, ancestor_type);
    for a in &ancestors {
        let kids: Vec<TopoShape> = match (child_type, ancestor_type) {
            (ShapeType::Vertex, ShapeType::Edge) | (ShapeType::Vertex, ShapeType::Face) => {
                vertices_of(a).into_iter().map(|v| v.0).collect()
            }
            (ShapeType::Edge, ShapeType::Face) => edges_of(a).into_iter().map(|e| e.0).collect(),
            _ => Vec::new(),
        };
        for k in kids {
            let list = out.entry(shape_key(&k)).or_default();
            if !list.iter().any(|x| x.same_tshape(a)) {
                list.push(a.clone());
            }
        }
    }
}
