//! Leftover GF+BOP mix: [`build_split_solids_full`].
//!
//! Not called by `BopBuilder`. Live GF solids are [`crate::bop_split_solids_occt`];
//! live Fuse/Cut/Common filter is [`crate::bop_bop::build_shape`]. Do not fold
//! obj/tool states back into FillImagesSolids. Private helpers stay on the
//! parent module.
//!
//! UNPORTED: a `merge_sharing_faces` pass used to run at the end of
//! [`build_split_solids_full`] and re-assembled pieces that share a face into
//! one solid. OCCT has no such stage (`BOPAlgo_Builder_3.cxx:579-616`); see the
//! note at the image-recording step below.

use super::*;

/// Splits every interfered source solid into its closed-shell pieces, adds the
/// internal faces of the other solids that lie inside it, and records the
/// selected pieces as solid images — the deepened [`build_split_solids`].
///
/// Mirrors the `FillIn3DParts` → `BuildSplitSolids` flow with the 3-D face
/// classification ported as a per-piece solid classification:
///
/// 1. for each source solid, collect the split faces of its own faces
///    ([`collect_solid_split_faces`]) plus the *internal faces* — faces of the
///    other argument solids lying strictly inside it (each added FORWARD and
///    REVERSED so both adjacent pieces close);
/// 2. group the whole set into closed shells ([`crate::shell_splitter::ShellSplitter`])
///    and wrap each shell into a solid;
/// 3. classify each piece with [`classify_solid_state`] against every *other*
///    argument solid: a piece lying strictly inside another solid is an
///    interior region — it contributes to the General-Fuse result only once, so
///    duplicate interior regions (the same geometric volume produced from
///    several source solids) collapse to their first occurrence;
/// 4. the surviving pieces are recorded as images of the source solid with the
///    origins back-map.
///
/// A source solid that was neither split nor hosts any internal face keeps no
/// image (it is returned as-is by the caller).
pub fn build_split_solids_full<B: BopBuildOps>(
    f: &mut B,
    objects: &[TopoShape],
    tools: &[TopoShape],
    obj_state: FaceState,
    tools_state: FaceState,
) -> Result<(), String> {
    let op = BopOp::from_states(obj_state, tools_state)?;
    let n = f.ds().nb_source_shapes();

    // The argument solids: used both as the source of internal faces and as the
    // reference set for the per-piece classification.
    let solids: Vec<TopoShape> = (0..n)
        .filter_map(|i| {
            let si = f.ds().shape_info(i)?;
            if si.shape_type() == ShapeType::Solid {
                Some(si.shape().clone())
            } else {
                None
            }
        })
        .collect();
    let tol = f.fuzzy_value().max(1e-7);

    // FillIn3DParts input shared by every solid's classification: the global
    // candidate faces (all source faces + their splits) and the solids' boxes.
    let candidates = collect_all_candidate_faces(f);
    let solid_boxes: HashMap<usize, BndBox> =
        solids.iter().map(|s| (shape_key(s), shape_bbox(s))).collect();

    // Interior regions already claimed by a previous source solid: an interior
    // piece with a bounding box already seen is a duplicate overlap region.
    let mut seen: Vec<BBoxKey> = Vec::new();
    // Every selected piece, tagged with its source solid, for the cross-solid
    // merge pass at the end.
    let mut collected: Vec<(TopoShape, TopoShape)> = Vec::new();

    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let solid = si.shape().clone();
        let is_object = objects.iter().any(|o| o.same_tshape(&solid));
        let is_tool = tools.iter().any(|t| t.same_tshape(&solid));
        // A Cut operation removes the tools entirely; their pieces never enter
        // the result (the tool faces inside the object are added as internal
        // faces of the object's pieces, so the cut stays closed).
        if op == BopOp::Cut && is_tool && !is_object {
            continue;
        }

        // 1. Split faces of the solid itself.
        let mut faces = collect_solid_split_faces(f, &solid);

        // 2. Internal faces — the FillIn3DParts classification of the global
        //    candidates (every source face and its split images) against this
        //    solid. A face strictly inside the solid is internal; a face on the
        //    solid's boundary whose edges all lie on the solid's own split
        //    faces covers an open boundary — e.g. the cylinder base over the
        //    box-top hole — and is included the same way so the hole closes.
        let own_edges: HashSet<EKey> = faces.iter().flat_map(|f| face_edges(&Face(f.clone()))).collect();
        let own_faces: HashSet<usize> = faces_of(&solid).into_iter().flat_map(|fc| {
            let k = shape_key(&fc.0);
            match f.history().image(&fc.0) {
                Some(imgs) => imgs.iter().map(|im| shape_key(im)).chain([k]).collect(),
                None => vec![k],
            }
        }).collect();
        let in_faces = classify_faces_in_solid(
            &candidates,
            &solid,
            &own_faces,
            &faces,
            &own_edges,
            &solid_boxes[&shape_key(&solid)],
            tol,
        );

        // A solid that neither splits nor hosts internal faces is unchanged.
        if !solid_interfered(f, &solid) && in_faces.is_empty() {
            continue;
        }
        faces.extend(in_faces);

        // 3. Group into closed shells and assemble one solid per shell.
        let pieces = split_faces_into_solids(&faces)?;
        if pieces.is_empty() {
            continue;
        }

        // 4. Select the pieces per the operation. A zero-volume shell (e.g. the
        //    two coincident orientations of a section face closing on
        //    themselves) cannot bound a solid and is dropped.
        let mut kept: Vec<TopoShape> = Vec::new();
        for p in pieces {
            let c = piece_center(&p);
            if piece_bbox_volume(&p) < 1e-9 {
                continue;
            }
            let interior = solids.iter().any(|other| {
                !other.same_tshape(&solid)
                    && classify_solid_state(other, &piece_center(&p), tol) == FaceState::In
            });
            let pass = match op {
                // Union: every distinct region survives; a duplicate interior
                // region collapses to its first occurrence.
                BopOp::Fuse => true,
                // Intersection: keep only the regions inside the other solid.
                BopOp::Common => interior,
                // Difference: keep the object's regions outside the tools.
                BopOp::Cut => !interior,
            };
            if !pass {
                continue;
            }
            if interior {
                let key = piece_bbox_key(&p);
                if seen.contains(&key) {
                    continue;
                }
                seen.push(key);
            }
            kept.push(p);
        }

        for p in kept {
            collected.push((solid.clone(), p));
        }
    }

    // 5. Record the images and the origins back-map. Mirrors the tail of
    //    `BOPAlgo_Builder::BuildSplitSolids` (`BOPAlgo_Builder_3.cxx:579-616`):
    //    every area built for a source solid is interned through the
    //    same-domain face-set map (`aMST`) and recorded as an image (with the
    //    source appended to `myOrigins`) of that source. The areas of one
    //    source stay distinct images.
    //
    //    UNPORTED: a `merge_sharing_faces` pass used to run here. It grouped
    //    the pieces that share a face — same `TShape`, or a coplanar pair with
    //    an area ratio under a constant — and rebuilt each group with
    //    `ShellSplitter`, pushing the merged solid once per contributing
    //    source. OCCT has no such stage: `BuildSplitSolids` keeps the
    //    `BOPAlgo_SplitSolid` areas as built, and the same-domain collapse is
    //    the `aMST` face-set intern above, which maps an area to an existing
    //    representative instead of re-closing shells (the live
    //    [`crate::bop_split_solids_occt`] does exactly that and also has no
    //    merge). The pass broke the piece partition: three union pieces of
    //    volume 0.5 were re-closed into one 1.5 solid recorded three times.
    for (solid, p) in collected {
        f.history_mut().add_image(&solid, p.clone());
        f.origins_mut().entry(shape_key(&p)).or_default().push(solid.clone());
    }

    // 7. Settle the internal vertices/edges/wires of the arguments (and those
    //    inside the source solids) into the split solids — the third
    //    `FillImagesSolids` stage, `BOPAlgo_Builder::FillInternalShapes`
    //    (`BOPAlgo_Builder_3.cxx`), which this port keeps in
    //    `crate::bop_build_common::fill_internal_shapes`. It reads the solid
    //    images just recorded above.
    crate::bop_build_common::fill_internal_shapes(f)?;
    Ok(())
}

/// The global candidate face list of the FillIn3DParts classification
/// (`BOPAlgo_Builder::FillIn3DParts`): every source FACE shape of the data
/// structure replaced by its image splits, an un-split face kept as-is,
/// deduplicated by `TShape` identity (the OCCT `aMFence` fence map).
fn collect_all_candidate_faces<B: BopBuildOps>(f: &B) -> Vec<TopoShape> {
    let n = f.ds().nb_source_shapes();
    let mut fence: Vec<TopoShape> = Vec::new();
    let mut out: Vec<TopoShape> = Vec::new();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Face {
            continue;
        }
        let face = si.shape().clone();
        match f.history().image(&face) {
            Some(imgs) => {
                for im in imgs {
                    if !fence.iter().any(|x| x.same_tshape(im)) {
                        fence.push(im.clone());
                        out.push(im.clone());
                    }
                }
            }
            None => {
                if !fence.iter().any(|x| x.same_tshape(&face)) {
                    fence.push(face.clone());
                    out.push(face);
                }
            }
        }
    }
    out
}

/// Classifies the global candidate faces against `solid` and returns those
/// lying inside it or covering an open boundary of its split faces, each as a
/// FORWARD + REVERSED pair (the shell of a split piece needs both orientations
/// to close against its neighbours).
///
/// Mirrors one `BOPAlgo_FillIn3DParts::Perform` task of
/// `BOPAlgo_Tools::ClassifyFaces`:
/// 1. candidates already part of the solid — its own faces and their split
///    images (`aMSF`) — are skipped;
/// 2. the pairwise box cull rejects faces whose bounding box does not reach
///    the solid's box (the `BOPTools_BoxTree` BVH selector, ported as a
///    [`shape_bbox_with_margin`] overlap check; `shape_bbox` samples a face's
///    surface on a 16×16 UV grid, so the cull over-approximates the face);
/// 3. the survivors are grouped into connexity blocks through edges that are
///    not on the solid and not degenerated ([`connexity_blocks`], the
///    `BOPAlgo_FillIn3DParts::MakeConnexityBlock` BFS), and one representative
///    of each block — the first face carrying a solid/degenerated edge, else
///    the block start — is classified with [`face_state_in_solid`] (the
///    `BOPTools_AlgoTools::ComputeState(Face, Solid)` fallback path of
///    `IsInternalFace`). Live FillIn3DParts uses [`crate::algo_tools_face::is_internal_face`]
///    (`GetFaceOff`). This leftover mix is not called by `BopBuilder`.
///    The representative's verdict applies to the whole block.
fn classify_faces_in_solid(
    candidates: &[TopoShape],
    solid: &TopoShape,
    own_faces: &HashSet<usize>,
    faces: &[TopoShape],
    own_edges: &HashSet<EKey>,
    solid_box: &BndBox,
    tol: f64,
) -> Vec<TopoShape> {
    // 1. Box cull + own-face filter (the `BOPTools_BoxTree` selector and the
    //    `aMSF` fence): candidates already part of the solid or whose box does
    //    not reach the solid's box are dropped before grouping.
    let mut sel: Vec<TopoShape> = Vec::new();
    for im in candidates {
        if own_faces.contains(&shape_key(im)) {
            continue;
        }
        if shape_bbox_with_margin(im, tol).is_out_box(solid_box) {
            continue;
        }
        sel.push(im.clone());
    }
    // 2. Connexity-block grouping of the survivors, classified through one
    //    representative face per block (OCCT `BOPAlgo_FillIn3DParts::Perform`).
    let mut in_faces: Vec<TopoShape> = Vec::new();
    for (block, rep) in connexity_blocks(&sel, own_edges) {
        let state = face_state_in_solid(&Face(rep.clone()), solid, own_edges, tol);
        let covering = state == FaceState::On && is_covering_face(&rep, faces, own_edges);
        if state != FaceState::In && !covering {
            continue;
        }
        for im in block {
            let mut fwd = im.clone();
            fwd.set_orientation(Orientation::Forward);
            in_faces.push(fwd);
            let mut rev = im.clone();
            rev.set_orientation(Orientation::Reversed);
            in_faces.push(rev);
        }
    }
    in_faces
}

/// Connexity blocks of `faces`, grouping faces that connect through edges
/// *not* on the solid and not degenerated, with the block's classification
/// representative.
///
/// Mirrors `BOPAlgo_FillIn3DParts::MakeConnexityBlock`:
/// - an edge on the solid (`solid_edges`, OCCT `aMSE`) or a degenerated edge
///   is a barrier: traversal does not cross it (so candidates glued along a
///   solid boundary are separate blocks), and the first face of the block
///   carrying such an edge becomes the representative `theFaceToClassify`;
/// - a block with no barrier edge is classified by its start face.
///
/// Returns `(block faces, representative face)` per block.
fn connexity_blocks(faces: &[TopoShape], solid_edges: &HashSet<EKey>) -> Vec<(Vec<TopoShape>, TopoShape)> {
    let n = faces.len();
    if n == 0 {
        return Vec::new();
    }
    // Per-face edges with the degeneracy flag (OCCT `BRep_Tool::Degenerated`)
    // and their undirected keys.
    let face_eds: Vec<Vec<(Edge, bool)>> = faces
        .iter()
        .map(|f| {
            edges_of(f)
                .into_iter()
                .map(|e| (e.clone(), BRepTool::is_degenerated(&e)))
                .collect()
        })
        .collect();
    let face_keys: Vec<Vec<EKey>> = face_eds
        .iter()
        .map(|es| es.iter().map(|(e, _)| edge_key(e)).collect())
        .collect();
    // Edge -> faces containing it (the candidate EF map `aMEFP`).
    let mut edge_faces: HashMap<EKey, Vec<usize>> = HashMap::new();
    for (i, keys) in face_keys.iter().enumerate() {
        for &k in keys {
            edge_faces.entry(k).or_default().push(i);
        }
    }
    let mut visited = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        // Breadth-first traversal over non-barrier shared edges, mirroring the
        // growing-list iteration of `MakeConnexityBlock`.
        let mut block: Vec<usize> = Vec::new();
        let mut rep: Option<usize> = None;
        let mut queue = std::collections::VecDeque::from([start]);
        visited[start] = true;
        while let Some(i) = queue.pop_front() {
            block.push(i);
            for ((_, deg), &k) in face_eds[i].iter().zip(&face_keys[i]) {
                if solid_edges.contains(&k) || *deg {
                    // Barrier edge: does not connect the block, and the first
                    // face carrying one is the classification representative.
                    if rep.is_none() {
                        rep = Some(i);
                    }
                    continue;
                }
                if let Some(neigh) = edge_faces.get(&k) {
                    for &j in neigh {
                        if !visited[j] {
                            visited[j] = true;
                            queue.push_back(j);
                        }
                    }
                }
            }
        }
        out.push((
            block.iter().map(|&i| faces[i].clone()).collect(),
            faces[rep.unwrap_or(start)].clone(),
        ));
    }
    out
}

/// Whether `im` is a covering face of the solid's own split `faces`: it is
/// On the solid's boundary (checked by the caller), at least one of its
/// boundary edges lies on the solid's own split faces, and it does not
/// duplicate an existing face (it fills a hole rather than coinciding with a
/// face already there).
///
/// Not every edge needs to be on the solid: a triangulated fan face that
/// covers a hole carries internal spoke edges (hub→rim) that are not edges of
/// the solid — OCCT's `IsInternalFace` likewise needs only the shared section
/// edge, not all of the face's edges.
/// On-boundary covering used only by leftover `build_split_solids_full`
/// (not `BopBuilder`). Live FillIn3DParts uses `GetFaceOff` via
/// [`crate::algo_tools_face::is_internal_face`].
fn is_covering_face(im: &TopoShape, own_faces: &[TopoShape], own_edges: &HashSet<EKey>) -> bool {
    let fk: HashSet<EKey> = face_edges(&Face(im.clone())).into_iter().collect();
    if fk.is_empty() || !fk.iter().any(|k| own_edges.contains(k)) {
        return false;
    }
    !own_faces.iter().any(|g| {
        let gk: HashSet<EKey> = face_edges(&Face(g.clone())).into_iter().collect();
        gk.len() == fk.len() && gk.iter().all(|k| fk.contains(k))
    })
}

/// The boolean operation derived from the object/tool face states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BopOp {
    /// A ∪ B — `(Out, Out)`.
    Fuse,
    /// A − B — `(Out, In)`.
    Cut,
    /// A ∩ B — `(In, In)`.
    Common,
}

impl BopOp {
    /// Derives the operation from the `BOPAlgo_BOP` states.
    pub fn from_states(obj: FaceState, tools: FaceState) -> Result<BopOp, String> {
        match (obj, tools) {
            (FaceState::Out, FaceState::Out) => Ok(BopOp::Fuse),
            (FaceState::Out, FaceState::In) => Ok(BopOp::Cut),
            (FaceState::In, FaceState::In) => Ok(BopOp::Common),
            _ => Err("BOPAlgo: unsupported object/tool states".to_string()),
        }
    }
}

/// Quantized axis-aligned bounding box of a piece, used as the identity key of
/// a geometric region for the interior-duplicate collapse.
type BBoxKey = (i64, i64, i64, i64, i64, i64);

/// The axis-aligned bounding box `(xmin, ymin, zmin, xmax, ymax, zmax)` of
/// `solid`'s vertices.
fn piece_bbox(solid: &TopoShape) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut mx = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in vertices_of(solid) {
        let p = BRepTool::vertex_point(&v);
        mn.0 = mn.0.min(p.x());
        mn.1 = mn.1.min(p.y());
        mn.2 = mn.2.min(p.z());
        mx.0 = mx.0.max(p.x());
        mx.1 = mx.1.max(p.y());
        mx.2 = mx.2.max(p.z());
    }
    (mn.0, mn.1, mn.2, mx.0, mx.1, mx.2)
}

/// The quantized bounding box of `solid`'s vertices.
fn piece_bbox_key(solid: &TopoShape) -> BBoxKey {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    let q = |x: f64| (x / 1e-6).round() as i64;
    (q(mnx), q(mny), q(mnz), q(mxx), q(mxy), q(mxz))
}

/// Volume of the bounding box of `solid`'s vertices (the axis-aligned region
/// extent) — zero for a degenerate shell.
fn piece_bbox_volume(solid: &TopoShape) -> f64 {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    (mxx - mnx) * (mxy - mny) * (mxz - mnz)
}

/// The centre of the bounding box of `solid`'s vertices — a representative
/// interior point for a closed piece.
fn piece_center(solid: &TopoShape) -> GpPnt {
    let (mnx, mny, mnz, mxx, mxy, mxz) = piece_bbox(solid);
    GpPnt::new(0.5 * (mnx + mxx), 0.5 * (mny + mxy), 0.5 * (mnz + mxz))
}

/// Whether `face` lies inside `solid` or on its boundary. Mirrors OCCT
/// `BOPTools_AlgoTools::ComputeState(Face, Solid)`: an edge of the face that is
/// not on the solid determines the state — its midpoint is classified; a face
/// all of whose edges lie on the solid is classified by an interior point.
/// `On` is not internal (a face on the solid's boundary is a shared outer face,
/// not a section face).
///
/// `solid_edges` is the solid's own split faces' edge set (the draft-solid
/// boundary). The boundary-vertex centroid the previous test used is ambiguous
/// for a face that straddles the solid (the box top's square-minus-circle has
/// its centroid inside the cylinder base, yet is mostly outside it); an edge
/// midpoint is not.
fn face_state_in_solid(
    face: &Face,
    solid: &TopoShape,
    solid_edges: &HashSet<EKey>,
    tol: f64,
) -> FaceState {
    for e in edges_of(&face.0) {
        if solid_edges.contains(&edge_key(&e)) {
            continue;
        }
        // An edge of the face not on the solid: classify its midpoint.
        let (first, last) = BRepTool::edge_parameters(&e);
        if first.is_finite() && last.is_finite() {
            if let Some(c) = BRepTool::edge_curve(&e) {
                return classify_solid_state(solid, &c.d0(0.5 * (first + last)), tol);
            }
        }
        let (Some(a), Some(b)) = crate::topo_tools_full::edge_vertices(&e) else { break };
        let pa = BRepTool::vertex_point(&a);
        let pb = BRepTool::vertex_point(&b);
        let mid = GpPnt::new(0.5 * (pa.x() + pb.x()), 0.5 * (pa.y() + pb.y()), 0.5 * (pa.z() + pb.z()));
        return classify_solid_state(solid, &mid, tol);
    }
    // All edges of the face lie on the solid: classify an interior point.
    if let Some(p) = face_sample_point(face) {
        return classify_solid_state(solid, &p, tol);
    }
    let vs = vertices_of(&face.0);
    if vs.is_empty() {
        return FaceState::Unknown;
    }
    let mut acc = GpPnt::new(0.0, 0.0, 0.0);
    let n = vs.len() as f64;
    for v in &vs {
        let p = BRepTool::vertex_point(v);
        acc = GpPnt::new(acc.x() + p.x(), acc.y() + p.y(), acc.z() + p.z());
    }
    classify_solid_state(solid, &GpPnt::new(acc.x() / n, acc.y() / n, acc.z() / n), tol)
}

#[cfg(test)]
#[path = "bop_build_solids_leftover_tests.rs"]
mod tests;
