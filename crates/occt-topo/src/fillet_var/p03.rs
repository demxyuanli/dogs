use super::prelude::*;
use super::*;

pub(super) fn find_end_face(solid: &TopoShape, f1: &Face, f2: &Face, p: &GpPnt) -> Result<Face, String> {
    let ends: Vec<Face> = faces_of(solid)
        .into_iter()
        .filter(|f| face_contains_point(f, p) && !is_same(&f.0, &f1.0) && !is_same(&f.0, &f2.0))
        .collect();
    if ends.len() != 1 {
        return Err(format!(
            "fillet_corner: expected one end face at {:?}, got {}",
            p,
            ends.len()
        ));
    }
    Ok(ends[0].clone())
}

/// The tangency trim of `edge` on `face`: `(edge_p0, u_dir, off)` such that the
/// retained side is `dot(from_p0, u_dir) >= off`.
pub(super) fn edge_tangency_trim(e: &EdgeInfo, face: &Face, radius: f64) -> Result<(GpPnt, GpVec, f64), String> {
    let spec = VarFilletSpec::new(radius, radius);
    let (_, u1, u2, cot_half) = sample_blend(&e.p0, &e.p1, &e.n1, &e.n2, &spec, 2)?;
    let is_f1 = is_same(&face.0, &e.f1.0);
    let is_f2 = is_same(&face.0, &e.f2.0);
    if is_f1 {
        Ok((e.p0, u1, radius * cot_half))
    } else if is_f2 {
        Ok((e.p0, u2, radius * cot_half))
    } else {
        Err("fillet_corner: edge is not adjacent to this face".to_string())
    }
}

/// The faces meeting at the corner vertex.
pub(super) fn corner_faces_at(solid: &TopoShape, p: &GpPnt) -> Vec<Face> {
    faces_of(solid)
        .into_iter()
        .filter(|f| face_contains_point(f, p))
        .collect()
}

/// A corner face to be rebuilt (a face at one or more shared vertices).
pub(super) struct CornerFaceRebuild {
    pub(super) face: Face,
    /// Tangency trims `(edge_p0, u_dir, off)` from adjacent filleted edges.
    pub(super) trims: Vec<(GpPnt, GpVec, f64)>,
    /// Corner arcs that live on this face (their circle centre/radius also
    /// define the corner disk that is cut out of the face).
    pub(super) arcs: Vec<(GpPnt, GpVec, f64, GpPnt, GpPnt)>,
    /// Sampled tangency polylines `(first, last, points)` of the adjacent
    /// blend faces, shared with them.
    pub(super) polylines: Vec<(GpPnt, GpPnt, Vec<GpPnt>)>,
}

/// A far end face to be rebuilt (perpendicular to an edge, not at a corner).
pub(super) struct EndFaceRebuild {
    pub(super) face: Face,
    pub(super) replacements: Vec<EndArcReplacement>,
}

/// Build the combined fillet of a run of consecutive edges that share vertices.
/// Each edge gets a constant-radius blend face trimmed to `[2R, len - 2R]` at
/// its corner ends, each shared vertex gets a spherical corner patch, and the
/// corner faces / far end faces are rebuilt with shared edges so the resulting
/// shell is closed.
pub(super) fn build_combined_run(
    solid: &TopoShape,
    run_edges: &[Edge],
    run_specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    let n = run_edges.len();
    if n < 2 {
        return Err("build_combined_run: needs at least 2 consecutive edges".to_string());
    }
    if n != run_specs.len() {
        return Err("build_combined_run: specs length mismatch".to_string());
    }

    // Per-edge geometry.
    let mut infos: Vec<EdgeInfo> = Vec::with_capacity(n);
    for (k, edge) in run_edges.iter().enumerate() {
        let spec = run_specs[k];
        spec.check()?;
        let adjacent = faces_touching_edge(solid, edge);
        if adjacent.len() != 2 {
            return Err(format!(
                "fillet_chain_corner: edge {k} touches {} faces (expected 2)",
                adjacent.len()
            ));
        }
        let f1 = adjacent[0].clone();
        let f2 = adjacent[1].clone();
        let n1 = face_outward_normal(&f1)?;
        let n2 = face_outward_normal(&f2)?;
        let (p0, p1) = BRepTool::edge_vertices(edge).ok_or("fillet_chain_corner: edge has no curve")?;
        let axis = GpVec::from_pnts(&p0, &p1).normalized();
        let len = p0.distance(&p1);
        infos.push(EdgeInfo { edge: edge.clone(), p0, p1, f1, f2, n1, n2, axis, len, spec });
    }

    // Shared corners between consecutive edges.
    pub(super) struct CornerInfo {
        pub(super) vertex: GpPnt,
        pub(super) radius: f64,
        pub(super) geom: CornerGeom,
        pub(super) k: usize,
    }
    let mut corners: Vec<CornerInfo> = Vec::new();
    for k in 0..n - 1 {
        let (pa0, pa1) = (infos[k].p0, infos[k].p1);
        let (pb0, pb1) = (infos[k + 1].p0, infos[k + 1].p1);
        let shared = find_shared_vertex(&pa0, &pa1, &pb0, &pb1)
            .ok_or("fillet_chain_corner: consecutive edges do not share a vertex")?;
        let ra = radius_at_endpoint(&infos[k].spec, &pa0, &pa1, &shared);
        let rb = radius_at_endpoint(&infos[k + 1].spec, &pb0, &pb1, &shared);
        let radius = ra.max(rb);
        let vshape = find_vertex_at(solid, &shared)
            .ok_or("fillet_chain_corner: shared vertex not found in the solid")?;
        let geom = CornerGeom::from_edges(solid, &vshape, &infos[k].edge, &infos[k + 1].edge, radius)?;
        corners.push(CornerInfo { vertex: shared, radius, geom, k });
    }

    // Corner ends per edge.
    let mut corner_at_start: Vec<Option<usize>> = vec![None; n];
    let mut corner_at_end: Vec<Option<usize>> = vec![None; n];
    for (ci, c) in corners.iter().enumerate() {
        if infos[c.k].p0.distance(&c.vertex) < 1e-9 {
            corner_at_start[c.k] = Some(ci);
        } else {
            corner_at_end[c.k] = Some(ci);
        }
        if infos[c.k + 1].p0.distance(&c.vertex) < 1e-9 {
            corner_at_start[c.k + 1] = Some(ci);
        } else {
            corner_at_end[c.k + 1] = Some(ci);
        }
    }

    // Per-edge corner radius (all incident corners of an edge must agree).
    let mut edge_radius = vec![0.0f64; n];
    for k in 0..n {
        let mut r: Option<f64> = None;
        for ci in [corner_at_start[k], corner_at_end[k]].iter().flatten() {
            let cr = corners[*ci].radius;
            if let Some(prev) = r {
                if (prev - cr).abs() > tol.max(1e-9) {
                    return Err(
                        "fillet_chain_corner: incident corner radii differ along an edge".to_string(),
                    );
                }
            }
            r = Some(cr);
        }
        edge_radius[k] = r.unwrap_or_else(|| infos[k].spec.r_start.max(infos[k].spec.r_end));
    }

    // Per-edge trimmed blend samples (used for the tangency polylines, the
    // blend faces, and the far-end tangency points).
    let mut edge_samples: Vec<Option<Vec<BlendSample>>> = vec![None; n];
    for k in 0..n {
        let e = &infos[k];
        let r = edge_radius[k];
        let t_lo = if corner_at_start[k].is_some() { 2.0 * r / e.len } else { 0.0 };
        let t_hi = if corner_at_end[k].is_some() { 1.0 - 2.0 * r / e.len } else { 1.0 };
        if t_hi - t_lo < 1e-9 {
            return Err("fillet_chain_corner: corner radius is too large for the edge".to_string());
        }
        let sub_p0 = e.p0.translated_vec(&e.axis.multiplied_scalar(t_lo * e.len));
        let sub_p1 = e.p0.translated_vec(&e.axis.multiplied_scalar(t_hi * e.len));
        let const_spec = VarFilletSpec::new(r, r);
        let (samples, _, _, _) = sample_blend(&sub_p0, &sub_p1, &e.n1, &e.n2, &const_spec, DEFAULT_VAR_SAMPLES)?;
        edge_samples[k] = Some(samples);
    }

    // Collect the corner faces to rebuild.
    let mut face_map: HashMap<usize, CornerFaceRebuild> = HashMap::new();
    let mut corner_face_keys: Vec<usize> = Vec::new();
    for c in &corners {
        for f in corner_faces_at(solid, &c.vertex) {
            let n_out = match face_outward_normal(&f) {
                Ok(n) => n,
                Err(_) => continue,
            };
            let Some(arc) = c.geom.arc_on_face(&n_out) else { continue };
            let key = Arc::as_ptr(&f.0.tshape) as usize;
            let entry = face_map.entry(key).or_insert_with(|| {
                corner_face_keys.push(key);
                CornerFaceRebuild {
                    face: f.clone(),
                    trims: Vec::new(),
                    arcs: Vec::new(),
                    polylines: Vec::new(),
                }
            });
            entry.arcs.push(arc);
            for &eidx in &[c.k, c.k + 1] {
                let e = &infos[eidx];
                if is_same(&e.f1.0, &f.0) {
                    let (p0, u, off) = edge_tangency_trim(e, &f, c.radius)?;
                    entry.trims.push((p0, u, off));
                    let poly: Vec<GpPnt> =
                        edge_samples[eidx].as_ref().unwrap().iter().map(|s| s.t1).collect();
                    entry.polylines.push((poly[0], *poly.last().unwrap(), poly));
                } else if is_same(&e.f2.0, &f.0) {
                    let (p0, u, off) = edge_tangency_trim(e, &f, c.radius)?;
                    entry.trims.push((p0, u, off));
                    let poly: Vec<GpPnt> =
                        edge_samples[eidx].as_ref().unwrap().iter().map(|s| s.t2).collect();
                    entry.polylines.push((poly[0], *poly.last().unwrap(), poly));
                }
            }
        }
    }

    // Collect the far end faces to rebuild.
    let mut end_map: HashMap<usize, EndFaceRebuild> = HashMap::new();
    let mut end_face_keys: Vec<usize> = Vec::new();
    for k in 0..n {
        let e = &infos[k];
        let r = edge_radius[k];
        if corner_at_start[k].is_none() {
            let end_face = find_end_face(solid, &e.f1, &e.f2, &e.p0)?;
            let s = sample_endpoint_blend(e, 0.0, r)?;
            let key = Arc::as_ptr(&end_face.0.tshape) as usize;
            let entry = end_map.entry(key).or_insert_with(|| {
                end_face_keys.push(key);
                EndFaceRebuild { face: end_face.clone(), replacements: Vec::new() }
            });
            entry.replacements.push(EndArcReplacement {
                corner: e.p0,
                t1: s.t1,
                t2: s.t2,
                center: s.center,
                axis: e.axis,
                radius: r,
                n1: e.n1,
            });
        }
        if corner_at_end[k].is_none() {
            let end_face = find_end_face(solid, &e.f1, &e.f2, &e.p1)?;
            let s = sample_endpoint_blend(e, 1.0, r)?;
            let key = Arc::as_ptr(&end_face.0.tshape) as usize;
            let entry = end_map.entry(key).or_insert_with(|| {
                end_face_keys.push(key);
                EndFaceRebuild { face: end_face.clone(), replacements: Vec::new() }
            });
            entry.replacements.push(EndArcReplacement {
                corner: e.p1,
                t1: s.t1,
                t2: s.t2,
                center: s.center,
                axis: e.axis,
                radius: r,
                n1: e.n1,
            });
        }
    }

    // Build everything with a single shared cache.
    let mut cache = EdgeCache::new();
    let mut new_faces: Vec<Face> = Vec::new();

    // Rebuilt corner faces.
    for &key in &corner_face_keys {
        let fr = &face_map[&key];
        let rebuilt = rebuild_corner_face(&fr.face, &fr.trims, &fr.arcs, &fr.polylines, tol, &mut cache)?;
        new_faces.push(rebuilt);
    }

    // Rebuilt far end faces.
    for &key in &end_face_keys {
        let er = &end_map[&key];
        let rebuilt = rebuild_face_with_end_arcs(&er.face, &er.replacements, &mut cache)?;
        new_faces.push(rebuilt);
    }

    // Blend faces (trimmed to the corner radius span).
    for k in 0..n {
        let e = &infos[k];
        let samples = edge_samples[k].as_ref().unwrap();
        let blend = build_var_blend_face(samples, &e.axis, &mut cache)?;
        new_faces.push(blend);
    }

    // Corner patch faces.
    for c in &corners {
        let patch = build_corner_patch_face(&c.geom, &mut cache)?;
        new_faces.push(patch);
    }

    // Assemble: keep every face except the corner faces and far end faces.
    let mut kept: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        let key = Arc::as_ptr(&f.0.tshape) as usize;
        if corner_face_keys.contains(&key) || end_face_keys.contains(&key) {
            continue;
        }
        kept.push(f);
    }
    kept.extend(new_faces);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&kept);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Sample the blend cross-section at an edge endpoint fraction (0 or 1) and
/// return the tangency points and arc centre at that cross-section.
pub(super) fn sample_endpoint_blend(
    e: &EdgeInfo,
    t: f64,
    radius: f64,
) -> Result<BlendSample, String> {
    let spec = VarFilletSpec::new(radius, radius);
    let (samples, _, _, _) = sample_blend(&e.p0, &e.p1, &e.n1, &e.n2, &spec, 2)?;
    Ok(if t < 0.5 { samples[0].clone() } else { samples[1].clone() })
}

/// Split `edge_indices` into maximal runs of edges that share a vertex.
pub(super) fn split_consecutive_runs(solid: &TopoShape, edge_indices: &[usize]) -> Vec<Vec<usize>> {
    let es = edges_of(solid);
    let mut runs: Vec<Vec<usize>> = Vec::new();
    for &i in edge_indices {
        if let Some(last) = runs.last_mut() {
            if let Some(&prev) = last.last() {
                if let (Some(pe), Some(ce)) = (es.get(prev), es.get(i)) {
                    if edges_share_vertex(pe, ce) {
                        last.push(i);
                        continue;
                    }
                }
            }
        }
        runs.push(vec![i]);
    }
    runs
}

pub(super) fn edges_share_vertex(a: &Edge, b: &Edge) -> bool {
    let (a0, a1) = edge_vertices(a);
    let (b0, b1) = edge_vertices(b);
    let (Some(a0), Some(a1)) = (a0, a1) else { return false };
    let (Some(b0), Some(b1)) = (b0, b1) else { return false };
    let (pa0, pa1) = (vertex_position(&a0), vertex_position(&a1));
    let (pb0, pb1) = (vertex_position(&b0), vertex_position(&b1));
    pa0.distance(&pb0) < 1e-9
        || pa0.distance(&pb1) < 1e-9
        || pa1.distance(&pb0) < 1e-9
        || pa1.distance(&pb1) < 1e-9
}

/// Blend a corner where the given edges meet: build the spherical corner patch
/// and trim the three corner faces. The edges themselves must already be (or
/// subsequently be) filleted; this inserts the patch face and re-closes the
/// shell. `edge_specs` maps edge indices (into `edges_of(solid)`) to specs.
pub fn fillet_corner_blend(
    solid: &TopoShape,
    corner_vertex: &Vertex,
    edge_specs: &[(usize, VarFilletSpec)],
    tol: f64,
) -> Result<TopoShape, String> {
    let p = vertex_position(corner_vertex);
    let es = edges_of(solid);
    let mut incident: Vec<(Edge, VarFilletSpec)> = Vec::new();
    for (i, spec) in edge_specs {
        let e = es.get(*i).ok_or("fillet_corner_blend: edge index out of range")?;
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            if vertex_position(&va).distance(&p) < 1e-9
                || vertex_position(&vb).distance(&p) < 1e-9
            {
                incident.push((e.clone(), *spec));
            }
        }
    }
    if incident.len() != 2 {
        return Err(format!(
            "fillet_corner_blend: expected 2 incident edge specs, got {}",
            incident.len()
        ));
    }
    let radius = incident
        .iter()
        .map(|(_, s)| s.r_end.max(s.r_start))
        .fold(0.0, f64::max);
    let geom = CornerGeom::from_edges(solid, corner_vertex, &incident[0].0, &incident[1].0, radius)?;

    let mut cache = EdgeCache::new();
    let all_faces = faces_of(solid);
    let cf = corner_faces_at(solid, &p);
    let mut rebuilt = Vec::new();
    for f in &cf {
        let n_out = match face_outward_normal(f) {
            Ok(n) => n,
            Err(_) => continue,
        };
        let Some(arc) = geom.arc_on_face(&n_out) else { continue };
        rebuilt.push(rebuild_corner_face(f, &[], &[arc], &[], tol, &mut cache)?);
    }
    let patch = build_corner_patch_face(&geom, &mut cache)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if cf.iter().any(|c| is_same(&f.0, &c.0)) {
            continue;
        }
        faces.push(f);
    }
    faces.extend(rebuilt);
    faces.push(patch);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Fillet a chain of edges, adding a spherical corner patch at every vertex
/// shared by two consecutive edges so the result stays a closed shell.
///
/// Non-adjacent edges are filleted sequentially with `fillet_edge_var` (as in
/// `fillet_edge_var_chain`). A run of consecutive edges sharing vertices is
/// built in one pass: each edge gets a constant-radius blend face trimmed at
/// the shared corners, and each shared vertex gets a rolling-ball corner patch
/// (`fillet_corner_blend`). The legacy `fillet_edge_var_chain` is unchanged.
pub fn fillet_edges_chain_with_corner(
    solid: &TopoShape,
    edge_indices: &[usize],
    specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    if edge_indices.len() != specs.len() {
        return Err(
            "fillet_edges_chain_with_corner: specs.len() must equal edge_indices.len()".to_string(),
        );
    }
    let original_edges = edges_of(solid);
    let runs = split_consecutive_runs(solid, edge_indices);
    let mut current = solid.clone();
    let mut spec_offset = 0;
    for run in &runs {
        let run_specs = &specs[spec_offset..spec_offset + run.len()];
        if run.len() == 1 {
            let i = run[0];
            let oe = original_edges
                .get(i)
                .ok_or_else(|| format!("fillet_edges_chain_with_corner: edge index {i} out of range"))?;
            let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("edge has no curve")?;
            let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
                format!(
                    "fillet_edges_chain_with_corner: edge {i} was consumed by an earlier fillet"
                )
            })?;
            current = fillet_edge_var(&current, &e, &run_specs[0], tol)?;
        } else {
            let mut run_edges: Vec<Edge> = Vec::with_capacity(run.len());
            for &i in run {
                let oe = original_edges
                    .get(i)
                    .ok_or_else(|| format!("fillet_edges_chain_with_corner: edge index {i} out of range"))?;
                let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("edge has no curve")?;
                let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
                    format!(
                        "fillet_edges_chain_with_corner: edge {i} was consumed by an earlier fillet"
                    )
                })?;
                run_edges.push(e);
            }
            current = build_combined_run(&current, &run_edges, run_specs, tol)?;
        }
        spec_offset += run.len();
    }
    Ok(current)
}
