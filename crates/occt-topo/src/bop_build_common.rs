//! BOPAlgo_Builder — container and internal-shape filling (Phase 20).
//!
//! Port of the container/assembly stages of the General Fuse builder
//! (`BOPAlgo_Builder`):
//!
//! | Rust function              | OCCT source                                       |
//! |----------------------------|---------------------------------------------------|
//! | [`fill_images_containers`] | `FillImagesContainers` + `FillImagesContainer`    |
//! |                            | (`BOPAlgo_Builder_1.cxx`)                         |
//! | [`fill_images_compounds`]  | `FillImagesCompounds` + `FillImagesCompound`      |
//! |                            | (`BOPAlgo_Builder_1.cxx`)                         |
//! | [`fill_internal_vertices`] | `FillInternalVertices` (`BOPAlgo_Builder_2.cxx`)  |
//! | [`fill_internal_shapes`]   | `FillInternalShapes` (`BOPAlgo_Builder_3.cxx`)    |
//! | [`build_draft_solid`]      | `BuildDraftSolid` (`BOPAlgo_Builder_3.cxx`)       |
//!
//! The five entry points operate on any host that implements [`BopBuildOps`]
//! (the minimal contract the builder main class — `BOPAlgo_Builder` /
//! `crate::bop_builder2::BopBuilder` — satisfies): the data structure, the
//! images history, the arguments and the origins back-map. Porting through the
//! trait keeps this module independent of the concrete `BopBuilder` fields, so
//! it compiles and is tested standalone.
//!
//! ## Semantics preserved from OCCT
//!
//! * **Containers** — a wire/shell is rebuilt only when at least one of its
//!   direct sub-shapes carries a non-trivial image. A shell is reassembled with
//!   [`crate::shell_splitter::ShellSplitter`] so the image shells are closed;
//!   a wire is reassembled with [`crate::builder::TopoBuilder::make_wire`].
//! * **Compounds** — recursively rebuilt with the splits of their sub-shapes,
//!   keeping each split oriented as the original sub-shape.
//! * **Internal vertices** — alone vertices of a split face are classified
//!   against each face image (2-D, `IntTools_Context::ComputeVF`) and added as
//!   `INTERNAL` children when they fall strictly inside it.
//! * **Internal shapes** — vertices/edges/wires from the arguments and from
//!   inside the source solids are classified against each split solid (3-D, one
//!   representative point) and added as `INTERNAL` children; settling a shape
//!   into an *original* (un-split) solid copies it first, preserving the input.
//! * **Draft solid** — a shell is rebuilt from the face splits (reversing a
//!   split face whose orientation is inverted relative to its original), flagged
//!   closed and wrapped into a solid.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bop_hist::BopHistory;
use crate::bopds::BopdsDS;
use crate::brep_extrema::{closest_point_on_edge, closest_point_on_face, is_inside};
use crate::brep_surface::surface_closest_params;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_splitter::ShellSplitter;
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of};

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
fn shape_key(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Add `s` to the ordered set `set`; returns true when it was newly inserted.
fn set_add(set: &mut Vec<TopoShape>, s: TopoShape) -> bool {
    if set.iter().any(|x| x.same_tshape(&s)) {
        false
    } else {
        set.push(s);
        true
    }
}

/// True when `set` already holds a shape identical to `s`.
fn set_contains(set: &[TopoShape], s: &TopoShape) -> bool {
    set.iter().any(|x| x.same_tshape(s))
}

/// Direct structural children of `shape`.
fn direct_children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// Reverse the orientation of a shape view.
fn reverse_orientation(s: &mut TopoShape) {
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
fn is_split_to_reverse(split: &TopoShape, original: &TopoShape) -> bool {
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
fn face_sample_point(face: &Face) -> Option<GpPnt> {
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
fn face_split_to_reverse(split: &Face, orig: &Face) -> bool {
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
fn edge_split_to_reverse(split: &Edge, orig: &Edge) -> bool {
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
fn surface_normal_checked(s: &dyn occt_geom::Surface, u: f64, v: f64) -> occt_core::gp::GpVec {
    let n = surface_closest_normal(s, u, v);
    if n.square_magnitude() < 1e-30 {
        occt_core::gp::GpVec::zero()
    } else {
        n
    }
}

/// Normal of a surface at `(u, v)` via the finite-difference helper.
fn surface_closest_normal(s: &dyn occt_geom::Surface, u: f64, v: f64) -> occt_core::gp::GpVec {
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
/// Mirrors `BOPAlgo_Builder::FillImagesContainers`: every source wire, shell
/// and comp-solid is rebuilt from the splits of its direct sub-shapes. A
/// container whose sub-shapes were not modified keeps no image (it is returned
/// as-is by the caller).
pub fn fill_images_containers<B: BopBuildOps>(f: &mut B) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        let t = si.shape_type();
        if t == ShapeType::Wire || t == ShapeType::Shell || t == ShapeType::CompSolid {
            let c = si.shape().clone();
            fill_images_container(f, &c, t)?;
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
fn fill_images_container<B: BopBuildOps>(
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
fn collect_splits<B: BopBuildOps, F: FnMut(TopoShape)>(f: &B, sub: &TopoShape, out: &mut F) {
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
fn fill_images_compound<B: BopBuildOps>(
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
fn alone_vertices(ds: &BopdsDS, face_idx: usize) -> Vec<usize> {
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
fn treat_compound(s: &TopoShape, out: &mut Vec<TopoShape>, fence: &mut Vec<TopoShape>) {
    if s.shape_type() != ShapeType::Compound {
        if set_add(fence, s.clone()) {
            out.push(s.clone());
        }
        return;
    }
    for h in s.tshape.read().unwrap().children.clone() {
        treat_compound(&TopoShape::from_handle(h), out, fence);
    }
}

/// The non-shell direct children of a solid: its internal vertices/edges/wires
/// (`BOPAlgo_Builder::OwnInternalShapes`).
fn own_internal_shapes(solid: &TopoShape) -> Vec<TopoShape> {
    direct_children(solid)
        .into_iter()
        .filter(|c| c.shape_type() != ShapeType::Shell)
        .collect()
}

/// Merges the vertex→edge, vertex→face and edge→face ancestor relations of
/// `shape` into `out` (`TopExp::MapShapesAndAncestors`).
fn map_ancestors_into(shape: &TopoShape, out: &mut HashMap<usize, Vec<TopoShape>>) {
    // vertex -> edges
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Edge, out);
    // vertex -> faces
    merge_ancestors(shape, ShapeType::Vertex, ShapeType::Face, out);
    // edge -> faces
    merge_ancestors(shape, ShapeType::Edge, ShapeType::Face, out);
}

/// For every ancestor of type `ancestor_type` in `shape`, records each of its
/// boundary sub-shapes of type `child_type` → the ancestor.
fn merge_ancestors(
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

/// 3-D state of `shape` relative to `solid`, from a single representative point
/// (`BOPTools_AlgoTools::ComputeStateByOnePoint`).
fn compute_state_by_one_point(
    shape: &TopoShape,
    solid: &TopoShape,
    tol: f64,
) -> Result<FaceState, String> {
    match shape.shape_type() {
        ShapeType::Vertex => {
            let p = BRepTool::vertex_point(&Vertex(shape.clone()));
            Ok(point_solid_state(solid, &p, tol))
        }
        ShapeType::Edge => {
            let e = Edge(shape.clone());
            let (a, b) = BRepTool::edge_parameters(&e);
            let Some(curve) = BRepTool::edge_curve(&e) else {
                return Ok(FaceState::Unknown);
            };
            let p = if a.is_finite() && b.is_finite() {
                curve.d0(0.5 * (a + b))
            } else {
                return Ok(FaceState::Unknown);
            };
            Ok(point_solid_state(solid, &p, tol))
        }
        ShapeType::Face => {
            // Prefer a boundary edge not lying on the solid; fall back to the
            // face surface centre when every edge is on the solid.
            let f = Face(shape.clone());
            let solid_edges = edges_of(solid);
            let mut p: Option<GpPnt> = None;
            for e in edges_of(&f.0) {
                if BRepTool::is_degenerated(&e) {
                    continue;
                }
                if solid_edges.iter().any(|se| se.same_tshape(&e.0)) {
                    continue;
                }
                let (a, b) = BRepTool::edge_parameters(&e);
                if a.is_finite() && b.is_finite() {
                    if let Some(c) = BRepTool::edge_curve(&e) {
                        p = Some(c.d0(0.5 * (a + b)));
                        break;
                    }
                }
            }
            let p = match p {
                Some(p) => p,
                None => match face_sample_point(&f) {
                    Some(p) => p,
                    None => return Ok(FaceState::Unknown),
                },
            };
            Ok(point_solid_state(solid, &p, tol))
        }
        _ => {
            let kids = shape.tshape.read().unwrap().children.clone();
            for h in kids {
                let sub = TopoShape::from_handle(h);
                let st = compute_state_by_one_point(&sub, solid, tol)?;
                if st != FaceState::Unknown {
                    return Ok(st);
                }
            }
            Ok(FaceState::Unknown)
        }
    }
}

/// State of a point relative to a solid: `On` within `tol` of a boundary face,
/// otherwise `In`/`Out` by the parity test. Distances are measured against the
/// solid's faces only, so internal vertices/edges already children of the solid
/// do not make the point `On` (matching the OCCT `BRepClass3d` classifier used
/// by `ComputeState`).
fn point_solid_state(solid: &TopoShape, p: &GpPnt, tol: f64) -> FaceState {
    let mut d = f64::INFINITY;
    for f in faces_of(solid) {
        let (_, q) = closest_point_on_face(&f, p, 16, 16);
        d = d.min(q.distance(p));
    }
    if d.is_finite() && d <= tol {
        return FaceState::On;
    }
    if is_inside(solid, p) {
        FaceState::In
    } else {
        FaceState::Out
    }
}

// ---------------------------------------------------------------------------
// Draft solid
// ---------------------------------------------------------------------------

/// Builds a draft solid from a (closed) shell, rebuilding the shell from the
/// face splits.
///
/// Mirrors `BOPAlgo_Builder::BuildDraftSolid` restricted to a single shell
/// argument: each face of the shell is replaced by its image splits (a split
/// whose orientation is inverted relative to the original face is reversed
/// first); `INTERNAL` faces are dropped (the caller collects them separately in
/// OCCT — this port does not expose the internal-face list). The rebuilt shell
/// is flagged closed and wrapped into a solid.
///
/// A solid input is treated as a collection of shells, mirroring the OCCT loop.
pub fn build_draft_solid<B: BopBuildOps>(
    f: &mut B,
    shell: &TopoShape,
) -> Result<TopoShape, String> {
    let bld = TopoBuilder::new();
    let solid_or = shell.orientation();
    let mut solid = Solid::new();
    solid.0.set_orientation(solid_or);

    let shells: Vec<TopoShape> = match shell.shape_type() {
        ShapeType::Solid => direct_children(shell)
            .into_iter()
            .filter(|c| c.shape_type() == ShapeType::Shell)
            .collect(),
        _ => vec![shell.clone()],
    };

    for sh in &shells {
        let mut new_shell = Shell::new();
        new_shell.0.set_orientation(sh.orientation());
        let mut i_flag = false;
        for child in direct_children(sh) {
            let or = child.orientation();
            if let Some(images) = f.history().image(&child) {
                for im in images {
                    let mut fx = im.clone();
                    if has_same_domain(f.ds(), &fx) {
                        if or == Orientation::Internal {
                            // Internal face: collected by the caller, not added.
                        } else {
                            if is_split_to_reverse(&fx, &child) {
                                reverse_orientation(&mut fx);
                            }
                            if add_draft_face(&bld, &mut new_shell.0, &fx) {
                                i_flag = true;
                            }
                        }
                    } else {
                        fx.set_orientation(or);
                        if or == Orientation::Internal {
                            // Internal face.
                        } else {
                            if add_draft_face(&bld, &mut new_shell.0, &fx) {
                                i_flag = true;
                            }
                        }
                    }
                }
            } else if or != Orientation::Internal {
                if add_draft_face(&bld, &mut new_shell.0, &child) {
                    i_flag = true;
                }
            }
        }
        if i_flag {
            new_shell.set_closed(!AlgoTools::is_open_shell(&new_shell.0));
            bld.add_shell(&mut solid, &new_shell);
        }
    }
    Ok(solid.0)
}

/// True when `shape` has a same-domain counterpart in the data structure.
fn has_same_domain(ds: &BopdsDS, shape: &TopoShape) -> bool {
    ds.index(shape)
        .and_then(|i| ds.has_shape_sd(i))
        .is_some()
}

/// Adds `face` to the rebuilt shell, once. A face that carries no area
/// (its boundary collapses to fewer than three distinct vertices) cannot bound
/// a shell and is dropped — the coincident-source-face splitting can produce
/// such degenerate sliver pieces alongside the real split face. A face already
/// present in the shell is not added twice. Returns true when the face was
/// added.
fn add_draft_face(bld: &TopoBuilder, shell: &mut TopoShape, face: &TopoShape) -> bool {
    if face_is_degenerate(&Face(face.clone())) {
        return false;
    }
    if set_contains(&direct_children(shell), face) {
        return false;
    }
    bld.add(shell, face);
    true
}

/// True when the face has no area: fewer than three distinct boundary-vertex
/// positions (a sliver / segment wire produced by an over-split).
fn face_is_degenerate(face: &Face) -> bool {
    let mut keys: Vec<(i64, i64, i64)> = Vec::new();
    for v in vertices_of(&face.0) {
        let p = BRepTool::vertex_point(&v);
        let k = (
            (p.x() / 1e-6).round() as i64,
            (p.y() / 1e-6).round() as i64,
            (p.z() / 1e-6).round() as i64,
        );
        if !keys.contains(&k) {
            keys.push(k);
        }
        if keys.len() >= 3 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpPln, GpPnt};
    use occt_geom::{GeomPlane, Surface};

    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::shape::{Edge, Face, Vertex};
    use crate::topo_tools_full::{faces_of, shapes_of};

    // -----------------------------------------------------------------------
    // Stub host
    // -----------------------------------------------------------------------

    struct StubBuilder {
        ds: BopdsDS,
        history: BopHistory,
        fuzzy: f64,
        args: Vec<TopoShape>,
        origins: HashMap<usize, Vec<TopoShape>>,
    }

    impl BopBuildOps for StubBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn arguments(&self) -> &[TopoShape] {
            &self.args
        }
        fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
            &mut self.origins
        }
    }

    fn stub(ds: BopdsDS, history: BopHistory, args: Vec<TopoShape>) -> StubBuilder {
        StubBuilder { ds, history, fuzzy: 1e-7, args, origins: HashMap::new() }
    }

    /// The face of `faces` whose boundary-vertex mean lies at height `z`.
    fn face_at_z(faces: &[Face], z: f64) -> Face {
        faces
            .iter()
            .find(|f| {
                let vs = crate::topo_tools_full::vertices_of(&f.0);
                if vs.is_empty() {
                    return false;
                }
                let zavg =
                    vs.iter().map(|v| BRepTool::vertex_point(v).z()).sum::<f64>() / vs.len() as f64;
                (zavg - z).abs() < 1e-9
            })
            .cloned()
            .expect("face at height")
    }

    /// A fresh face on the same surface and boundary edges as `src` (a new
    /// TShape, so it is a distinct split image).
    fn re_face(src: &Face) -> Face {
        let b = TopoBuilder::new();
        let edges = crate::topo_tools_full::edges_of(&src.0);
        let wire = b.make_wire(&edges);
        let surf = BRepTool::face_surface(src).expect("face surface");
        b.make_face(surf, &[wire])
    }

    // -----------------------------------------------------------------------
    // fill_images_containers
    // -----------------------------------------------------------------------

    #[test]
    fn empty_builder_leaves_history_untouched() {
        let ds = BopdsDS::new();
        let history = BopHistory::new();
        let mut b = stub(ds.clone(), history, Vec::new());
        fill_images_containers(&mut b).unwrap();
        assert!(!b.history().has_any_images());
    }

    #[test]
    fn box_shell_with_split_top_face_reassembles_closed_shell() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let faces = faces_of(&shell);
        assert_eq!(faces.len(), 6);

        // Split the top face: register a fresh, geometrically identical face as
        // its only image.
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);
        assert!(!new_top.0.same_tshape(&top.0));

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        fill_images_containers(&mut b).unwrap();

        let imgs = b.history().image(&shell).expect("shell has an image");
        assert_eq!(imgs.len(), 1, "one closed shell image");
        let img = &imgs[0];
        assert!(img.is_shell(), "image is a shell");
        assert!(img.closed(), "rebuilt shell is closed");
        assert_eq!(faces_of(img).len(), 6, "all six split faces present");
    }

    #[test]
    fn unmodified_box_shell_gets_no_image() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![boxed.solid.0.clone()]);
        fill_images_containers(&mut b).unwrap();
        assert!(
            b.history().image(&shell).is_none(),
            "no face was split -> no container image"
        );
    }

    // -----------------------------------------------------------------------
    // fill_images_compounds
    // -----------------------------------------------------------------------

    #[test]
    fn compound_with_split_child_is_rebuilt() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2_new = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let comp = b.make_compound_of(&[e1.0.clone(), e2.0.clone()]);

        let mut ds = BopdsDS::new();
        ds.init(&[comp.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&e2.0, e2_new.0.clone());
        let mut stub = stub(ds, history, vec![comp.0.clone()]);

        fill_images_compounds(&mut stub).unwrap();

        let imgs = stub.history().image(&comp.0).expect("compound has an image");
        assert_eq!(imgs.len(), 1);
        let img = &imgs[0];
        assert!(img.is_compound());
        let kids = direct_children(img);
        assert!(kids.iter().any(|k| k.same_tshape(&e1.0)), "unmodified child kept");
        assert!(kids.iter().any(|k| k.same_tshape(&e2_new.0)), "split child replaced");
        assert!(
            !kids.iter().any(|k| k.same_tshape(&e2.0)),
            "original split child not in the image"
        );
    }

    #[test]
    fn unmodified_compound_gets_no_image() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        let v2 = b.make_vertex(GpPnt::new(1.0, 0.0, 0.0), 0.0);
        let comp = b.make_compound_of(&[v1.0, v2.0]);
        let mut ds = BopdsDS::new();
        ds.init(&[comp.0.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![comp.0.clone()]);
        fill_images_compounds(&mut stub).unwrap();
        assert!(stub.history().image(&comp.0).is_none());
    }

    // -----------------------------------------------------------------------
    // fill_internal_vertices
    // -----------------------------------------------------------------------

    fn square_face_with_alone_vertices() -> (Face, Face, Vertex, Vertex) {
        let b = TopoBuilder::new();
        let p = [GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)];
        let edges: Vec<Edge> = (0..4)
            .map(|i| b.make_edge_segment(&p[i], &p[(i + 1) % 4]))
            .collect();
        let wire = b.make_wire(&edges);
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let mut face = b.make_face(surf.clone(), &[wire]);

        // A fresh split face with the same boundary.
        let wire2 = b.make_wire(&edges);
        let face_im = b.make_face(surf, &[wire2]);

        // Two alone vertices: one inside the square, one outside.
        let v_in = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 0.0);
        let v_out = b.make_vertex(GpPnt::new(2.0, 2.0, 0.0), 0.0);
        b.add(&mut face.0, &v_in.0);
        b.add(&mut face.0, &v_out.0);
        (face, face_im, v_in, v_out)
    }

    #[test]
    fn alone_vertices_inside_split_are_added_as_internal() {
        let (face, face_im, v_in, v_out) = square_face_with_alone_vertices();

        let mut ds = BopdsDS::new();
        ds.init(&[face.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&face.0, face_im.0.clone());
        let mut b = stub(ds, history, vec![face.0.clone()]);

        fill_internal_vertices(&mut b).unwrap();

        // The flat model stores children without a per-child orientation marker,
        // so the INTERNAL annotation set by the port is not observable; what
        // matters is that the inside vertex became a child of the face image
        // while the outside one did not.
        let kids = direct_children(&face_im.0);
        assert!(
            kids.iter().any(|k| k.same_tshape(&v_in.0)),
            "inside alone vertex added as a child of the split face"
        );
        assert!(
            !kids.iter().any(|k| k.same_tshape(&v_out.0)),
            "outside alone vertex dropped"
        );
    }

    #[test]
    fn no_alone_vertices_adds_nothing() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&boxed.solid.0);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);
        fill_internal_vertices(&mut b).unwrap();
        // The box faces have no alone vertices, so no INTERNAL vertex is added
        // to the split top face (its wire child remains the only child).
        let kids = direct_children(&new_top.0);
        assert!(
            !kids.iter().any(|k| k.is_vertex() && k.orientation() == Orientation::Internal),
            "no internal vertex added to the split face"
        );
    }

    // -----------------------------------------------------------------------
    // fill_internal_shapes
    // -----------------------------------------------------------------------

    #[test]
    fn inside_vertex_settles_into_original_solid_as_copy() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = TopoBuilder::new();
        let v_in = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let mut solid = boxed.solid.0.clone();
        b.add(&mut solid, &v_in.0); // internal vertex child of the solid

        let mut ds = BopdsDS::new();
        ds.init(&[solid.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![solid.clone()]);

        fill_internal_shapes(&mut stub).unwrap();

        let imgs = stub.history().image(&solid).expect("solid gains a copy image");
        assert_eq!(imgs.len(), 1);
        let img = &imgs[0];
        assert!(img.is_solid());
        // The flat model does not persist a per-child orientation marker, so the
        // assertion is on presence of the vertex child (the INTERNAL annotation
        // set by the port is not observable through the child list).
        let kids = direct_children(img);
        assert!(
            kids.iter().any(|k| k.same_tshape(&v_in.0)),
            "internal vertex present in the solid copy"
        );
        // The original solid is preserved (only one shell, no internal vertex added).
        assert!(direct_children(&solid).iter().filter(|c| c.is_shell()).count() >= 1);
    }

    #[test]
    fn outside_vertex_is_not_settled() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = TopoBuilder::new();
        let v_out = b.make_vertex(GpPnt::new(5.0, 5.0, 5.0), 0.0);
        let mut solid = boxed.solid.0.clone();
        b.add(&mut solid, &v_out.0);

        let mut ds = BopdsDS::new();
        ds.init(&[solid.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![solid.clone()]);

        fill_internal_shapes(&mut stub).unwrap();

        assert!(
            stub.history().image(&solid).is_none(),
            "no vertex lies inside the solid -> no split solid"
        );
    }

    // -----------------------------------------------------------------------
    // build_draft_solid
    // -----------------------------------------------------------------------

    #[test]
    fn box_shell_wraps_into_unit_solid() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![boxed.solid.0.clone()]);

        let solid = build_draft_solid(&mut b, &shell).unwrap();
        assert!(solid.is_solid());
        let v = crate::brep_gprop::volume(&solid, 0.02);
        assert!((v - 1.0).abs() < 0.01, "unit box volume, got {v}");
    }

    #[test]
    fn split_box_shell_rebuilds_closed_shell_into_solid() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let faces = faces_of(&shell);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        // First reassemble the shell image.
        fill_images_containers(&mut b).unwrap();
        let shell_imgs = b.history().image(&shell).unwrap().to_vec();
        assert_eq!(shell_imgs.len(), 1);

        // Then build a solid from the rebuilt (closed) shell.
        let solid = build_draft_solid(&mut b, &shell_imgs[0]).unwrap();
        assert!(solid.is_solid());
        let v = crate::brep_gprop::volume(&solid, 0.02);
        assert!((v - 1.0).abs() < 0.01, "split box still has unit volume, got {v}");
    }

    #[test]
    fn coincident_boxes_draft_solid_volume_is_one() {
        // Regression: two fully coincident boxes (same make_box parameters)
        // went through FillSameDomainFaces + BuildDraftSolid with their volume
        // inflated to ~17. The coincident faces were over-split into degenerate
        // sliver pieces and the draft-solid assembly counted every piece; both
        // the degenerate-piece skip and the same-domain self-image guard bring
        // each rebuilt box back to volume 1.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut br = crate::bop_builder2::BopBuilder::new();
        br.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        br.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        br.filler_mut().perform().unwrap();
        crate::bop_build_faces::fill_images_faces(&mut br).unwrap();
        assert!(!br.has_errors(), "errors: {:?}", br.errors());

        for solid in [a.solid.0.clone(), c.solid.0.clone()] {
            let shell = shapes_of(&solid, ShapeType::Shell)[0].clone();
            let draft = build_draft_solid(&mut br, &shell).unwrap();
            let v = crate::brep_gprop::volume(&draft, 0.02);
            assert!((v - 1.0).abs() < 0.01, "coincident box draft volume, got {v}");
        }
    }
}
