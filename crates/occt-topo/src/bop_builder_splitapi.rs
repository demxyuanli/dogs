//! Face split / edge class / tolerance heal. Split from `bop_builder.rs`.

use std::collections::{HashMap};
use std::sync::Arc;


use occt_core::gp::{GpPnt, GpPnt2d, GpVec};


use crate::abs::ShapeType;
use crate::brep_extrema::is_inside;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;

use crate::shape::{Face, Shell, TopoShape};
use crate::shell_check::{shell_is_closed};

use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};

use crate::bop_builder_core::{BoolOp};
use crate::bop_builder_dispatch::expand_compound;
use crate::bop_builder_repair::face_is_degenerate;
use crate::bop_builder_planar::*;

// ---------------------------------------------------------------------------
// Face splitting along surface intersections
// ---------------------------------------------------------------------------

/// Split every face of `shape` along its intersection curves with the faces in
/// `pairs`, producing clean sub-faces with complete boundary wires.
///
/// `pairs` lists pairs of face indices (into [`faces_of`]) whose intersection
/// should be cut. For each pair the plane–plane intersection segment is computed
/// (via [`face_face_segments_local`]) and both faces are split along it with the
/// same polygon splitter and weld/edge-map the exact boolean uses, so the split
/// edges coincide and sub-faces share boundary vertices. A face that is not
/// cut (the segment misses it, or it is non-planar) is kept whole. The rebuilt
/// boundary is returned as a solid when closed, otherwise as a shell.
pub fn split_faces_along_intersections(shape: &TopoShape, pairs: &[(usize, usize)], tol: f64) -> Result<TopoShape, String> {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    for &(i, j) in pairs {
        if i >= faces.len() || j >= faces.len() {
            return Err(format!(
                "split_faces_along_intersections: face index ({i}, {j}) out of range for {} faces",
                faces.len()
            ));
        }
    }
    let mut segs: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); faces.len()];
    for &(i, j) in pairs {
        let s = face_face_segments_local(&faces[i], &faces[j], tol);
        for seg in &s {
            segs[i].push(*seg);
            segs[j].push(*seg);
        }
    }
    let bld = TopoBuilder::new();
    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let mut result_faces: Vec<Face> = Vec::new();
    for (i, f) in faces.iter().enumerate() {
        if segs[i].is_empty() {
            result_faces.push(f.clone());
            continue;
        }
        let Some(pln) = face_plane_local(f) else {
            result_faces.push(f.clone());
            continue;
        };
        let Some(poly2d) = face_polygon_local(f, &pln) else {
            result_faces.push(f.clone());
            continue;
        };
        let segs2d: Vec<(GpPnt2d, GpPnt2d)> = segs[i]
            .iter()
            .filter_map(|(a, b)| {
                let a2 = project_point_to_plane(&pln, a);
                let b2 = project_point_to_plane(&pln, b);
                if a2.distance(&b2) < 1e-12 {
                    None
                } else {
                    Some((a2, b2))
                }
            })
            .collect();
        // The same BuilderFace-style planar arrangement the boolean uses:
        // open section lines split the face into regions, closed loops carve
        // holes (a split face may keep a hole loop as a multi-wire sub-face).
        let regions = trace_planar_regions(&poly2d, &segs2d);
        if regions.len() >= 2 {
            let subs: Vec<SubFace> = regions
                .iter()
                .filter_map(|r| region_to_subface(&bld, &pln, r, &mut weld, &mut edge_map))
                .collect();
            if subs.len() >= 2 {
                result_faces.extend(subs.into_iter().map(|sf| sf.face));
                continue;
            }
        }
        result_faces.push(f.clone());
    }
    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    let out: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };
    Ok(out)
}

// ---------------------------------------------------------------------------
// Boolean-result edge classification (TopOpeBRep_BuildTool)
// ---------------------------------------------------------------------------

/// Classification of an edge of a boolean result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClass {
    /// The edge is shared by several boundary patches (referenced by three or
    /// more result faces) — the coincidence of two bodies' edges in the result.
    Shared,
    /// The edge lies in the interior of a face's surface: its two adjacent
    /// faces are coplanar, so the edge is a seam/split line on a flat region.
    OnFace,
    /// The edge is not part of a closed two-face boundary: a free edge or a
    /// wire edge interior to the result topology.
    Internal,
    /// A regular outer-boundary edge (exactly two non-coplanar faces meet).
    External,
}

/// Per-class counts of a [`classify_boolean_edges`] result.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeClassCounts {
    /// Number of `Shared` edges.
    pub shared: usize,
    /// Number of `OnFace` edges.
    pub on_face: usize,
    /// Number of `Internal` edges.
    pub internal: usize,
    /// Number of `External` edges.
    pub external: usize,
}

/// Classify every edge of a boolean result (`TopOpeBRep_BuildTool` style).
///
/// For each distinct edge of `result` the faces referencing it are counted:
///
/// * 3+ faces → [`EdgeClass::Shared`] (several boundary patches — the two
///   bodies' coincident edges — meet on one edge);
/// * exactly 2 coplanar faces → [`EdgeClass::OnFace`] (the edge is a seam lying
///   on a face's surface interior);
/// * exactly 2 non-coplanar faces → [`EdgeClass::External`] (regular outer
///   boundary edge);
/// * 0 or 1 face → [`EdgeClass::Internal`] (free / interior wire edge).
///
/// The operation is accepted for API symmetry with TopOpeBRep; the geometric
/// classification is operation-independent. Use [`edge_class_counts`] to
/// aggregate, and [`classify_edges_with_operands`] for the operand-aware
/// classification that reports edges lying on either input body's surface.
pub fn classify_boolean_edges(result: &TopoShape, _op: BoolOp) -> Vec<EdgeClass> {
    let faces = faces_of(result);
    let mut face_refs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        for e in edges_of(&f.0) {
            face_refs.entry(Arc::as_ptr(&e.0.tshape) as usize).or_default().push(fi);
        }
    }
    let edges = edges_of(result);
    edges
        .iter()
        .map(|e| {
            let ptr = Arc::as_ptr(&e.0.tshape) as usize;
            let fids = face_refs.get(&ptr).map(|v| v.as_slice()).unwrap_or(&[]);
            classify_result_edge(&faces, fids)
        })
        .collect()
}

/// Classify one edge from the indices of the faces referencing it.
fn classify_result_edge(faces: &[Face], fids: &[usize]) -> EdgeClass {
    match fids.len() {
        0 | 1 => EdgeClass::Internal,
        2 => {
            let (f1, f2) = (&faces[fids[0]], &faces[fids[1]]);
            match (face_plane_local(f1), face_plane_local(f2)) {
                (Some(p1), Some(p2)) if planes_coincident(&p1, &p2, 1e-6) => EdgeClass::OnFace,
                _ => EdgeClass::External,
            }
        }
        _ => EdgeClass::Shared,
    }
}

/// Aggregate [`EdgeClass`]es into [`EdgeClassCounts`].
pub fn edge_class_counts(classes: &[EdgeClass]) -> EdgeClassCounts {
    let mut c = EdgeClassCounts::default();
    for cl in classes {
        match cl {
            EdgeClass::Shared => c.shared += 1,
            EdgeClass::OnFace => c.on_face += 1,
            EdgeClass::Internal => c.internal += 1,
            EdgeClass::External => c.external += 1,
        }
    }
    c
}

/// One-line summary of [`EdgeClassCounts`].
pub fn edge_class_counts_summary(c: &EdgeClassCounts) -> String {
    format!(
        "{} edge(s): {} shared, {} on-face, {} internal, {} external",
        c.shared + c.on_face + c.internal + c.external,
        c.shared,
        c.on_face,
        c.internal,
        c.external
    )
}

/// Operand-aware edge classification for a boolean result.
///
/// For each edge of `result`, its midpoint is probed against both operand
/// boundaries:
///
/// * on **both** operands' surfaces → [`EdgeClass::Shared`] (the two bodies
///   meet along this edge);
/// * on exactly one operand's surface → [`EdgeClass::OnFace`];
/// * inside both operands (a seam hidden inside the union/intersection) →
///   [`EdgeClass::Internal`];
/// * otherwise → [`EdgeClass::External`].
///
/// This is the full `TopOpeBRep` classification; [`classify_boolean_edges`] is
/// its result-only form.
pub fn classify_edges_with_operands(a: &TopoShape, b: &TopoShape, result: &TopoShape, _op: BoolOp, tol: f64) -> Vec<EdgeClass> {
    let tol = tol.max(1e-9);
    let faces_a = faces_of(a);
    let faces_b = faces_of(b);
    edges_of(result)
        .iter()
        .map(|e| {
            let m = match BRepTool::edge_vertices(e) {
                Some((p1, p2)) => {
                    let v = GpVec::from_pnts(&p1, &p2);
                    p1.translated_vec(&v.multiplied_scalar(0.5))
                }
                None => GpPnt::zero(),
            };
            let on_a = point_on_surface(&faces_a, &m, tol);
            let on_b = point_on_surface(&faces_b, &m, tol);
            if on_a && on_b {
                EdgeClass::Shared
            } else if on_a || on_b {
                EdgeClass::OnFace
            } else if is_inside(a, &m) && is_inside(b, &m) {
                EdgeClass::Internal
            } else {
                EdgeClass::External
            }
        })
        .collect()
}

/// Is `p` (within `tol`) on the surface of any of `faces`?
fn point_on_surface(faces: &[Face], p: &GpPnt, tol: f64) -> bool {
    for f in faces {
        if let Some(pln) = face_plane_local(f) {
            if point_in_face_polygon(f, &pln, p, tol) {
                return true;
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Tolerance healing
// ---------------------------------------------------------------------------

/// Report of [`heal_tolerance`] / [`heal_tolerance_report`].
#[derive(Debug, Clone)]
pub struct HealToleranceReport {
    /// The healed shape.
    pub healed: TopoShape,
    /// Number of vertex pairs merged (weld coincident vertices within `tol`).
    pub welded_vertices: usize,
    /// Number of edges removed (shorter than the heal tolerance).
    pub removed_edges: usize,
    /// Number of degenerate faces removed.
    pub removed_faces: usize,
    /// Non-fatal diagnostics (open shells, dropped faces, …).
    pub warnings: Vec<String>,
}

/// Heal a shape at a tolerance: weld near-coincident vertices, remove small
/// edges and degenerate faces (`ShapeFix_Shape`-style).
///
/// This is the tolerance-welding wrapper used to tidy a boolean result before
/// downstream processing. See [`heal_tolerance_report`] for the detailed
/// statistics; this entry point returns just the healed shape.
pub fn heal_tolerance(shape: &TopoShape, tol: f64) -> Result<TopoShape, String> {
    Ok(heal_tolerance_report(shape, tol)?.healed)
}

/// Heal a shape at a tolerance and report what was fixed.
///
/// Welds vertices closer than `tol` to one canonical instance, removes edges
/// shorter than `8·tol` (the small-edge threshold), and drops degenerate faces
/// (fewer than 3 boundary edges, or a planar face of near-zero area). A
/// `Compound` is healed component-wise. The healed shape is rebuilt through the
/// existing healing machinery, so shared edges/vertices are preserved.
pub fn heal_tolerance_report(shape: &TopoShape, tol: f64) -> Result<HealToleranceReport, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let children = expand_compound(shape);
        let mut healed_children: Vec<TopoShape> = Vec::with_capacity(children.len());
        let mut welded = 0usize;
        let mut removed_edges = 0usize;
        let mut removed_faces = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        for c in &children {
            let r = heal_tolerance_report(c, tol)?;
            welded += r.welded_vertices;
            removed_edges += r.removed_edges;
            removed_faces += r.removed_faces;
            warnings.extend(r.warnings);
            healed_children.push(r.healed);
        }
        let bld = TopoBuilder::new();
        return Ok(HealToleranceReport {
            healed: bld.make_compound_of(&healed_children).0,
            welded_vertices: welded,
            removed_edges,
            removed_faces,
            warnings,
        });
    }

    let (s1, welded) = crate::shhealing::weld_coincident_vertices(shape, tol);
    let (s2, removed_edges) = crate::shhealing::remove_small_edges(&s1, tol * 8.0);
    let (s3, removed_faces) = remove_degenerate_faces(&s2, tol);
    let mut warnings: Vec<String> = Vec::new();
    for sh in shapes_of(&s3, ShapeType::Shell) {
        if !shell_is_closed(&Shell(sh)) {
            warnings.push("heal_tolerance: an open shell remains after healing".into());
        }
    }
    Ok(HealToleranceReport { healed: s3, welded_vertices: welded, removed_edges, removed_faces, warnings })
}

/// Is a face degenerate — fewer than 3 boundary edges, or a (planar) face with
/// near-zero area?
fn face_is_degenerate_tol(face: &Face, tol: f64) -> bool {
    let mut ecount = 0usize;
    for w in wires_of_face(face) {
        ecount += edges_of_wire(&w).len();
    }
    ecount < 3 || face_is_degenerate(face, tol)
}

/// Rebuild a shape without its degenerate faces. Shapes without shell structure
/// (bare wires/faces) are returned unchanged.
fn remove_degenerate_faces(shape: &TopoShape, tol: f64) -> (TopoShape, usize) {
    let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
    if shells.is_empty() {
        return (shape.clone(), 0);
    }
    let bld = TopoBuilder::new();
    let mut removed = 0usize;
    let mut out: Vec<Shell> = Vec::new();
    for sh in &shells {
        let faces = faces_of(&sh.0);
        let keep: Vec<Face> = faces.iter().filter(|f| !face_is_degenerate_tol(f, tol)).cloned().collect();
        removed += faces.len() - keep.len();
        if !keep.is_empty() {
            out.push(bld.make_shell(&keep));
        }
    }
    if out.is_empty() {
        return (bld.make_compound_of(&[]).0, removed);
    }
    if out.len() == 1 {
        let s = out.pop().unwrap();
        if shape.is_solid() && shell_is_closed(&s) {
            return (bld.make_solid(&[s]).0, removed);
        }
        return (s.0, removed);
    }
    let shapes: Vec<TopoShape> = out.into_iter().map(|s| s.0).collect();
    (bld.make_compound_of(&shapes).0, removed)
}
