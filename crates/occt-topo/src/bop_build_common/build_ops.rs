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
    /// The origins back-map (split shape → original shapes). Keyed by
    /// `TopTools_ShapeMapHasher` like `BOPAlgo_Builder::myOrigins`
    /// (`BOPAlgo_Builder.hxx:501`).
    fn origins_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>>;
}

// ---------------------------------------------------------------------------
// Shape-set helpers
// ---------------------------------------------------------------------------

/// `NCollection_Map` / `NCollection_IndexedMap` insert, keyed by
/// `TopTools_ShapeMapHasher` (`aMFence`, `aMSOr`, `aMSI`,
/// `BOPAlgo_Builder_3.cxx:637-643`) — [`same_shape`] = `TopoDS_Shape::IsSame`.
pub(super) fn set_add(set: &mut Vec<TopoShape>, s: TopoShape) -> bool {
    if set.iter().any(|x| same_shape(x, &s)) {
        false
    } else {
        set.push(s);
        true
    }
}

/// True when `set` already holds a shape identical to `s` under
/// [`same_shape`] (`TopTools_ShapeMapHasher`).
pub(super) fn set_contains(set: &[TopoShape], s: &TopoShape) -> bool {
    set.iter().any(|x| same_shape(x, s))
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
/// Thin wrapper over the single faithful
/// [`crate::algo_tools_face::is_split_to_reverse`]
/// (`BOPTools_AlgoTools::IsSplitToReverse(Sp, Sr, ctx, err)`,
/// `BOPTools_AlgoTools.cxx:1255-1291`): edges go through
/// `BOPTools_AlgoTools.cxx:1432-1500`, faces through `:1316-1427`.
pub(crate) fn is_split_to_reverse_ctxless(split: &TopoShape, original: &TopoShape) -> bool {
    let mut ctx = IntToolsContext::new();
    crate::algo_tools_face::is_split_to_reverse(split, original, &mut ctx).unwrap_or(false)
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

/// Builds the image of a single container (`BOPAlgo_Builder::FillImagesContainer`,
/// `BOPAlgo_Builder_1.cxx:221-274`).
///
/// The container is rebuilt as **one** container of the same type: every direct
/// sub-shape contributes its splits (a split whose orientation is inverted
/// relative to the original is reversed first, `cxx:262-270`), or itself when it
/// has no image (`cxx:253-258`). The assembled container is flagged with
/// `BRep_Tool::IsClosed` (`cxx:274`).
///
/// `BOPTools_AlgoTools::MakeContainer` (`BOPTools_AlgoTools.cxx:1600-1645`) is
/// an empty `MakeWire` / `MakeShell` / `MakeCompSolid`; the port's
/// `TopoBuilder::make_wire` / `make_shell` and `TopoShape::new` are that empty
/// constructor, so no splitting of the container into several closed shells
/// happens here (that happens later, in `FillImagesSolids`).
pub(super) fn fill_images_container<B: BopBuildOps>(
    f: &mut B,
    container: &TopoShape,
    container_type: ShapeType,
) -> Result<(), String> {
    let children = direct_children(container);
    // `cxx:225-233`: check if any of the sub-shapes of the container has been
    // modified. A sub-shape whose only image is itself is not a modification.
    let mut modified = false;
    for ss in &children {
        if let Some(im) = f.history().image(ss) {
            if im.len() != 1 || !same_shape(&im[0], ss) {
                modified = true;
                break;
            }
        }
    }
    if !modified {
        return Ok(());
    }

    // `cxx:247-271`: add every sub-shape or its splits into the new container.
    let mut parts: Vec<TopoShape> = Vec::new();
    for ss in &children {
        collect_splits(f, ss, &mut |im| parts.push(im));
    }
    let bld = TopoBuilder::new();
    let Some(image) = make_container_image(&bld, &parts, container_type) else {
        return Ok(());
    };
    // `cxx:274`: `aCIm.Closed(BRep_Tool::IsClosed(aCIm))`.
    image.set_closed(BRepTool::is_closed(&image));
    f.history_mut().add_image(container, image);
    Ok(())
}

/// `BOPTools_AlgoTools::MakeContainer(theType, C)` + the `BRep_Builder::Add`
/// loop of `FillImagesContainer`: an empty container of `container_type` filled
/// with `parts`.
fn make_container_image(
    bld: &TopoBuilder,
    parts: &[TopoShape],
    container_type: ShapeType,
) -> Option<TopoShape> {
    match container_type {
        ShapeType::Wire => {
            let edges: Vec<Edge> = parts
                .iter()
                .filter(|s| s.shape_type() == ShapeType::Edge)
                .map(|s| Edge(s.clone()))
                .collect();
            Some(bld.make_wire(&edges).0)
        }
        ShapeType::Shell => {
            let faces: Vec<Face> = parts
                .iter()
                .filter(|s| s.shape_type() == ShapeType::Face)
                .map(|s| Face(s.clone()))
                .collect();
            Some(bld.make_shell(&faces).0)
        }
        ShapeType::CompSolid | ShapeType::Compound => {
            let mut container = TopoShape::new(container_type);
            for p in parts {
                bld.add(&mut container, p);
            }
            Some(container)
        }
        _ => None,
    }
}

/// Pushes the image splits of `sub` into `out`, reversing those whose
/// orientation is inverted relative to `sub` (`BOPAlgo_Builder_1.cxx:264-269`:
/// `!aSSIm.IsEqual(aSS) && IsSplitToReverse(aSSIm, aSS)`). When `sub` has no
/// image, it is pushed itself (`cxx:253-258`).
pub(super) fn collect_splits<B: BopBuildOps, F: FnMut(TopoShape)>(f: &B, sub: &TopoShape, out: &mut F) {
    if let Some(im) = f.history().image(sub) {
        for ims in im {
            let mut s = ims.clone();
            if !crate::topo_tools::is_equal(&s, sub) && is_split_to_reverse_ctxless(&s, sub) {
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
    let mut lsc: Vec<TopoShape> = Vec::new();
    let mut fence: Vec<TopoShape> = Vec::new();
    for a in f.arguments() {
        treat_compound(a, &mut lsc, &mut fence);
    }
    let mut largs: Vec<TopoShape> = Vec::new();
    fence.clear();
    for s in &lsc {
        match s.shape_type() {
            ShapeType::Wire => {
                for e in edges_of_wire(&Wire(s.clone())) {
                    if set_add(&mut fence, e.0.clone()) {
                        largs.push(e.0);
                    }
                }
            }
            ShapeType::Vertex | ShapeType::Edge => largs.push(s.clone()),
            _ => {}
        }
    }
    fence.clear();
    let mut msi: Vec<TopoShape> = Vec::new();
    for s in &largs {
        if !set_add(&mut fence, s.clone()) {
            continue;
        }
        let t = s.shape_type();
        if t == ShapeType::Vertex || t == ShapeType::Edge || t == ShapeType::Wire {
            if f.history().has_image(s) {
                for im in f.history().image(s).unwrap() {
                    set_add(&mut msi, im.clone());
                }
            } else {
                set_add(&mut msi, s.clone());
            }
        }
    }

    // 2. — internal vertices/edges from the source solids + ancestor map.
    fence.clear();
    let mut lsd: Vec<TopoShape> = Vec::new();
    let mut msx: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut msor: Vec<TopoShape> = Vec::new();
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
                    set_add(&mut msi, im.clone());
                }
            } else {
                set_add(&mut msi, s);
            }
        }
        // Ancestors of the splits (or of the solid itself when unsplit).
        if f.history().has_image(&solid) {
            for sp in f.history().image(&solid).unwrap() {
                if set_add(&mut fence, sp.clone()) {
                    map_ancestors_into(sp, &mut msx);
                    lsd.push(sp.clone());
                }
            }
        } else if set_add(&mut fence, solid.clone()) {
            map_ancestors_into(&solid, &mut msx);
            lsd.push(solid.clone());
            set_add(&mut msor, solid);
        }
    }

    // 3. — keep only the candidates not already tied to a face of the solids.
    let mut lsi: Vec<TopoShape> = Vec::new();
    for s in &msi {
        let tied = msx.contains_key(&shape_key(s));
        if !tied {
            lsi.push(s.clone());
        }
    }

    // 4. — nothing to settle.
    if lsi.is_empty() {
        return Ok(());
    }

    // 5. — settle the candidates into the split solids.
    let bld = TopoBuilder::new();
    for sd in &mut lsd {
        let mut i = 0;
        while i < lsi.len() {
            let mut si_shape = lsi[i].clone();
            si_shape.set_orientation(Orientation::Internal);
            let state =
                compute_state_by_one_point(&si_shape, sd, 1e-11).unwrap_or(FaceState::Unknown);
            if state != FaceState::In {
                i += 1;
                continue;
            }
            if set_contains(&msor, sd) {
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
                msor.retain(|x| !same_shape(x, sd));
                *sd = sdx.0;
            } else {
                bld.add(sd, &si_shape);
            }
            lsi.remove(i);
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
/// `shape` into `out` (`TopExp::MapShapesAndAncestors`,
/// `BOPAlgo_Builder_3.cxx:770-783`). `out` is the `aMSx`
/// `NCollection_IndexedDataMap<TopoDS_Shape, List, TopTools_ShapeMapHasher>`
/// (`BOPAlgo_Builder_3.cxx:635-636`), so its key is [`shape_key`].
pub(super) fn map_ancestors_into(shape: &TopoShape, out: &mut HashMap<ShapeMapKey, Vec<TopoShape>>) {
    // vertex -> edges
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Edge, out);
    // vertex -> faces
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Face, out);
    // edge -> faces
    merge_ancestors(shape, ShapeType::Edge, ShapeType::Face, out);
}

/// For every ancestor occurrence of type `ancestor_type` in `shape`, records
/// each of its boundary sub-shape occurrences of type `child_type` → the
/// ancestor (`TopExp::MapShapesAndAncestors`, `TopExp.cxx:90-107`). Both
/// explorers are occurrence walks (`TopExp_Explorer`), not the uniquifying
/// `TopExp::MapShapes`, so a seam edge stored twice under one face contributes
/// the face twice and vice versa.
pub(super) fn merge_ancestors(
    shape: &TopoShape,
    child_type: ShapeType,
    ancestor_type: ShapeType,
    out: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) {
    let ancestors = crate::topo_tools::explorer_shapes(shape, ancestor_type);
    for a in &ancestors {
        let kids: Vec<TopoShape> = match (child_type, ancestor_type) {
            (ShapeType::Vertex, ShapeType::Edge) | (ShapeType::Vertex, ShapeType::Face) => {
                crate::topo_tools::explorer_shapes(a, ShapeType::Vertex)
            }
            (ShapeType::Edge, ShapeType::Face) => {
                crate::topo_tools::explorer_shapes(a, ShapeType::Edge)
            }
            _ => Vec::new(),
        };
        for k in kids {
            // `TopExp.cxx:96-105` appends unconditionally: one entry per
            // `(ancestor, sub-shape)` occurrence pair, duplicates included.
            out.entry(shape_key(&k)).or_default().push(a.clone());
        }
    }
}
