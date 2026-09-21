use super::prelude::*;
use super::*;


/// Align the pcurve endpoints with the surface projection of the edge's 3D
/// endpoints (`BRepLib::SameParameter`-style endpoint correction).
///
/// When the pcurve spans the full edge range `[a, b]`, each endpoint is
/// compared against the projection of the edge's 3D curve point at `a` / `b`
/// onto the face surface. If the 2D deviation exceeds `tol`, the pcurve is
/// re-fitted as a degree-1 B-spline over `[a, b]` whose endpoints are the
/// projected values. A pcurve trimmed to a proper sub-range already ends on
/// the face boundary and is returned unchanged; periodic (circle) pcurves are
/// also left alone.
pub(super) fn align_pcurve_endpoints(
    pc: &Arc<dyn Curve2d>,
    edge: &Edge,
    face: &Face,
    tol: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    if pc.is_periodic() {
        return Ok(pc.clone());
    }
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Ok(pc.clone());
    }
    let (p0, p1) = (pc.first_parameter(), pc.last_parameter());
    // A bounded pcurve that no longer spans the edge range was clipped by the
    // trim step — its endpoints are the intended face-boundary points.
    if p0.is_finite() && p1.is_finite() && ((a - p0).abs() > 1e-9 || (b - p1).abs() > 1e-9) {
        return Ok(pc.clone());
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else { return Ok(pc.clone()) };
    let q0 = pc.d0(a);
    let q1 = pc.d0(b);
    let Some((u0, v0, _)) = project_point_on_face(face, &c3d.d0(a)) else { return Ok(pc.clone()) };
    let Some((u1, v1, _)) = project_point_on_face(face, &c3d.d0(b)) else { return Ok(pc.clone()) };
    let t0t = GpPnt2d::new(u0, v0);
    let t1t = GpPnt2d::new(u1, v1);
    if q0.distance(&t0t) <= tol && q1.distance(&t1t) <= tol {
        return Ok(pc.clone());
    }
    // Re-fit a degree-1 B-spline over [a, b] with the corrected endpoints.
    let n = 33;
    let mut pts: Vec<GpPnt2d> = (0..n)
        .map(|i| pc.d0(a + (b - a) * i as f64 / (n - 1) as f64))
        .collect();
    pts[0] = t0t;
    pts[n - 1] = t1t;
    let bs = bspline_from_pts_2d(&pts, a, b, 1)?;
    Ok(Arc::new(bs))
}
/// Fit a degree-`degree` clamped B-spline through `pts` over `[a, b]` — a
/// polyline for degree 1, mirroring the fitter used by
/// [`crate::pcurve_full`].
pub(super) fn bspline_from_pts_2d(
    pts: &[GpPnt2d],
    a: f64,
    b: f64,
    degree: usize,
) -> Result<Geom2dBSplineCurve, String> {
    if pts.len() < 2 {
        return Err("bspline_from_pts_2d: need at least 2 points".into());
    }
    if degree >= pts.len() {
        return Err("bspline_from_pts_2d: degree below pole count".into());
    }
    let mut knots = occt_core::bspl::knots::build_uniform_knots(pts.len(), degree);
    for k in knots.iter_mut() {
        *k = a + (b - a) * *k;
    }
    let xs: Vec<f64> = pts.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, knots, degree).map_err(|e| e.to_string())
}

/// Update the endpoint-vertex tolerances of `edge` after a pcurve was built on
/// `face` (OCCT `UpdateVertices`): grow each vertex to cover the deviation
/// between the 3D curve point and the surface point mapped by the pcurve at
/// the boundary parameter, and floor it at the face 2D tolerance.
pub(super) fn update_vertices_full(edge: &Edge, face: &Face, pc: &Arc<dyn Curve2d>, face_tol: f64) {
    let Some(c3d) = GeometryRegistry::global().edge_curve(&edge.0) else { return };
    let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else { return };
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() {
        return;
    }
    let verts = edge_vertex_shapes(edge);
    for (j, t) in [(0usize, a), (1, b)] {
        if let Some(v) = verts.get(j) {
            let vtx = Vertex(v.clone());
            let tol = BRepTool::vertex_tolerance(&vtx);
            let p3d = c3d.d0(t);
            let q = pc.d0(t);
            let p3dx = surf.d0(q.x(), q.y());
            let d = p3d.distance(&p3dx);
            let mut want = tol;
            if d > tol {
                want = d + crate::algo_tools::D_TOLERANCE;
            }
            if face_tol > want {
                want = face_tol;
            }
            if want > tol {
                vtx.set_tolerance(want);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ProcessDE helpers
// ---------------------------------------------------------------------------

/// Find the pave blocks of the edges of the face `face` that pass through the
/// vertex with DS index `v`.
///
/// Port of `BOPAlgo_PaveFiller::FindPaveBlocks`. The face's pave list (edge
/// index + parameter range) is scanned and every block of those edges whose
/// first or last bound references `v` is appended to `out`.
pub fn find_pave_blocks<F: PaveFillerLike>(f: &F, v: usize, face: usize, out: &mut Vec<BopdsPaveBlock>) {
    let Some(fi) = f.ds().face_info(face) else {
        return;
    };
    let tuples: Vec<(usize, f64, f64)> = fi
        .paves_in()
        .iter()
        .chain(fi.paves_on().iter())
        .chain(fi.paves().iter())
        .copied()
        .collect();
    for (edge, t1, t2) in tuples {
        let blocks = f.ds().pave_blocks(edge);
        for pb in blocks {
            let (n1, n2) = pb.indices();
            if n1 == v || n2 == v {
                out.push(pb.clone());
            }
        }
        let _ = (t1, t2);
    }
}

/// Add split points to the pave block `pbd` of a degenerated edge: intersect
/// the 2D curve of the degenerated edge with the 2D curves of every passing
/// pave block and record the intersection parameters as extra paves.
///
/// Port of `BOPAlgo_PaveFiller::FillPaves`. When the 2D intersection produces
/// no point, the endpoint of the passing curve that corresponds to `v` is
/// projected onto the degenerated curve instead (the OCCT fallback).
pub fn fill_paves<F: PaveFillerLike>(
    f: &mut F,
    v: usize,
    e: usize,
    face: usize,
    lpb: &[BopdsPaveBlock],
    pbd: &mut BopdsPaveBlock,
) {
    let Some(de_shape) = f.ds().shape(e).cloned() else { return };
    let Some(face_shape) = f.ds().shape(face).cloned() else { return };
    let de = Edge(de_shape);
    let fa = Face(face_shape);
    let Ok(c2d_de) = boptools_2d::make_2d(&de, &fa) else { return };
    let tol = PCONFUSION;
    for pb in lpb {
        let n_e = pb.edge();
        if n_e >= f.ds().nb_shapes() {
            continue;
        }
        let Some(pb_shape) = f.ds().shape(n_e).cloned() else { continue };
        let pbe = Edge(pb_shape);
        let Ok(c2d) = boptools_2d::make_2d(&pbe, &fa) else { continue };
        let hits = intersect_curves(c2d_de.as_ref(), c2d.as_ref(), tol);
        if hits.is_empty() {
            let t = if v == pb.pave1().index() {
                pb.pave1().parameter()
            } else {
                pb.pave2().parameter()
            };
            let p2d = c2d.d0(t);
            if let Some(proj) = project_point_on_curve(c2d_de.as_ref(), &p2d, tol) {
                add_split_point(pbd, BopdsPave::new(v, proj.parameter), tol);
            }
        } else {
            for hit in hits {
                add_split_point(pbd, BopdsPave::new(v, hit.u1), tol);
            }
        }
    }
}

/// Validate and add `pave` as an extra pave of `pbd`: the parameter must be
/// strictly inside the block range and not collide with an existing pave.
/// Returns true when the pave was added (OCCT `AddSplitPoint`).
pub(super) fn add_split_point(pbd: &mut BopdsPaveBlock, pave: BopdsPave, tol: f64) -> bool {
    let (td1, td2) = pbd.range();
    let t = pave.parameter();
    if t - td1 < tol || td2 - t < tol {
        return false;
    }
    if pbd.contains_parameter(t, tol).is_some() {
        return false;
    }
    pbd.append_ext_pave1(pave);
    true
}

/// Split the degenerated edge `de` on the face `df`, creating a new
/// (degenerated) sub-edge for each of its pave blocks delimited by a new
/// vertex. Port of `BOPAlgo_PaveFiller::MakeSplitEdge` (ProcessDE branch).
pub fn make_split_edge_de<F: PaveFillerLike>(f: &mut F, de: usize, df: usize) -> Result<(), String> {
    let Some(de_shape) = f.ds().shape(de).cloned() else { return Ok(()) };
    let blocks = f.ds().pave_blocks(de).to_vec();
    if blocks.is_empty() {
        return Ok(());
    }
    let multi = blocks.len() > 1;
    let mut split_indices: Vec<Option<usize>> = Vec::with_capacity(blocks.len());
    for pb in &blocks {
        let (n_v1, n_v2) = pb.indices();
        let (a_t1, a_t2) = pb.range();
        if f.ds().is_new_shape(n_v1) || multi {
            let v1 = f.ds().shape(n_v1).cloned().ok_or_else(|| format!("make_split_edge_de: vertex {n_v1} not found"))?;
            let v2 = f.ds().shape(n_v2).cloned().ok_or_else(|| format!("make_split_edge_de: vertex {n_v2} not found"))?;
            let sp = make_degenerate_split_edge(&Edge(de_shape.clone()), &v1, a_t1, &v2, a_t2)?;
            let mut si = BopdsShapeInfo::new(sp.0);
            si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
            let n_sp = f.ds_mut().append_info(si);
            split_indices.push(Some(n_sp));
        } else {
            // No split needed: drop the block (the whole list is cleared below).
            split_indices.push(None);
        }
    }
    // Clear the block list of the degenerated edge and re-point surviving blocks.
    if split_indices.iter().any(|s| s.is_some()) {
        let blocks = f.ds_mut().change_pave_blocks_mut(de);
        for (pb, si) in blocks.iter_mut().zip(split_indices.iter()) {
            if let Some(n_sp) = si {
                pb.set_edge(*n_sp);
            }
        }
    } else {
        {
            let si = f.ds_mut().change_shape_info(de).expect("de shape info");
            si.pb_reference = -1;
        }
        f.ds_mut().change_pave_blocks_mut(de).clear();
    }
    let _ = df;
    Ok(())
}

/// Build a degenerated sub-edge from the base edge `orig` bounded by `v1`/`v2`
/// (OCCT `MakeSplitEdge1`): a collapsed copy carrying the two vertices and the
/// degenerated flag.
pub(super) fn make_degenerate_split_edge(
    orig: &Edge,
    v1: &TopoShape,
    t1: f64,
    v2: &TopoShape,
    t2: f64,
) -> Result<Edge, String> {
    let b = TopoBuilder::new();
    let p = if let Some(c) = BRepTool::edge_curve(orig) {
        c.d0(0.5 * (t1 + t2))
    } else {
        GpPnt::zero()
    };
    let mut e = b.make_edge_segment(&p, &p);
    if t1 < t2 {
        b.add(&mut e.0, &v1.oriented(Orientation::Forward));
        b.add(&mut e.0, &v2.oriented(Orientation::Reversed));
    } else {
        b.add(&mut e.0, &v1.oriented(Orientation::Reversed));
        b.add(&mut e.0, &v2.oriented(Orientation::Forward));
    }
    if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
        g.degenerated = true;
        g.tolerance = CONFUSION;
        GeometryRegistry::global().set_edge(&e.0, g);
    }
    Ok(e)
}

/// Port of `BOPAlgo_PaveFiller::ProcessDE` (`BOPAlgo_PaveFiller_8.cxx`).
/// Degenerated edges are those whose `ShapeInfo` carries `HasFlag(nF)`.
pub fn process_de<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    for an_edge_index in 0..n {
        let Some(an_edge_info) = f.ds().shape_info(an_edge_index) else {
            continue;
        };
        if an_edge_info.shape_type() != ShapeType::Edge {
            continue;
        }
        let Some(n_f) = an_edge_info.flag() else {
            continue;
        };
        let Some(a_sif) = f.ds().shape_info(n_f) else {
            continue;
        };
        let n_v0 = an_edge_info.sub_shapes().first().copied().unwrap_or(0);
        let n_v = f.ds().get_same_domain_index(n_v0);
        match a_sif.shape_type() {
            ShapeType::Face => {
                let mut a_lpb_out: Vec<BopdsPaveBlock> = Vec::new();
                find_pave_blocks(f, n_v, n_f, &mut a_lpb_out);
                if !a_lpb_out.is_empty() {
                    let pbd = {
                        let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                        let Some(first) = blocks.first().cloned() else {
                            continue;
                        };
                        first
                    };
                    let mut pbd = pbd;
                    fill_paves(f, n_v, an_edge_index, n_f, &a_lpb_out, &mut pbd);
                    let mut out = Vec::new();
                    pbd.update(&mut out, true);
                    let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                    blocks.clear();
                    blocks.extend(out);
                }
                make_split_edge_de(f, an_edge_index, n_f)?;
            }
            ShapeType::Edge => {
                let Some(a_de) = f.ds().shape(an_edge_index).cloned() else {
                    continue;
                };
                let Some(a_vn) = f.ds().shape(n_v).cloned() else {
                    continue;
                };
                let sp = make_degenerate_split_edge(&Edge(a_de), &a_vn, 0.0, &a_vn, 0.0)?;
                let a_si = BopdsShapeInfo::new(sp.0);
                let n_en = f.ds_mut().append_info(a_si);
                let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                if let Some(a_pbd) = blocks.first_mut() {
                    a_pbd.set_edge(n_en);
                }
            }
            _ => {}
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// SplitEdge (single-edge split helper)
// ---------------------------------------------------------------------------

/// Create a new sub-edge of `edge` bounded by the vertices `v1`/`v2` at the
/// parameters `t1`/`t2` and append it to the data structure.
///
/// Port of `BOPAlgo_PaveFiller::SplitEdge`. Returns the DS index of the new
/// edge.
pub fn make_split_edge<F: PaveFillerLike>(
    f: &mut F,
    edge: usize,
    v1: usize,
    t1: f64,
    v2: usize,
    t2: f64,
) -> Result<usize, String> {
    let e_shape = f
        .ds()
        .shape(edge)
        .ok_or_else(|| format!("make_split_edge: edge {edge} not found"))?
        .clone();
    let v1_shape = f
        .ds()
        .shape(v1)
        .ok_or_else(|| format!("make_split_edge: vertex {v1} not found"))?
        .clone();
    let v2_shape = f
        .ds()
        .shape(v2)
        .ok_or_else(|| format!("make_split_edge: vertex {v2} not found"))?
        .clone();
    let sp = AlgoTools::make_split_edge(&Edge(e_shape), Some(&v1_shape), t1, Some(&v2_shape), t2)?;
    let mut si = BopdsShapeInfo::new(sp.0);
    si.change_sub_shapes().extend_from_slice(&[v1, v2]);
    Ok(f.ds_mut().append_info(si))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// The vertex child shapes of an edge (its two boundary vertices).
pub(super) fn edge_vertex_shapes(edge: &Edge) -> Vec<TopoShape> {
    edge.0
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == ShapeType::Vertex)
        .cloned()
        .collect()
}
