use super::prelude::*;
use super::*;

/// Replace a sharp straight edge of `solid` with a variable-radius rolling-ball
/// fillet. `spec` carries the two end radii and the interpolation law; `tol` is
/// the geometric tolerance used for the curved-face clip boundary.
///
/// The two adjacent faces must be planar; the end faces (perpendicular to the
/// edge, containing its endpoints) are rebuilt with the blend arc at the local
/// radius. Returns the rebuilt solid, or an error when the edge is
/// non-manifold, a face is non-planar, or a radius does not fit.
pub fn fillet_edge_var(
    solid: &TopoShape,
    edge: &Edge,
    spec: &VarFilletSpec,
    tol: f64,
) -> Result<TopoShape, String> {
    spec.check()?;
    // Small slack for the curved clip boundary so the tangency points stay
    // robustly inside the retained region. Kept far below `tol` (the chamfer
    // search tolerance is derived from `tol` in `rebuild_var_face`).
    let _ = tol;
    let eps = 1e-9;

    let adjacent = faces_touching_edge(solid, edge);
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet_var: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let n1 = face_outward_normal(f1)?;
    let n2 = face_outward_normal(f2)?;

    let (p0, p1) = BRepTool::edge_vertices(edge).ok_or("fillet_var: edge has no evaluable curve")?;
    let (samples, u1, u2, cot_half) = sample_blend(&p0, &p1, &n1, &n2, spec, DEFAULT_VAR_SAMPLES)?;
    let axis_dir = GpVec::from_pnts(&p0, &p1).normalized();
    let edge_len = p0.distance(&p1);

    // Retained region of each adjacent face: points whose distance from the
    // edge (along the in-wedge direction) reaches the local tangency distance
    // R(t)·cot(θ/2), where t is the point's projection onto the edge axis.
    let inside1 = |q: &GpPnt| {
        let t = (GpVec::from_pnts(&p0, q).dot(&axis_dir) / edge_len).clamp(0.0, 1.0);
        let off = radius_at(spec, t) * cot_half;
        GpVec::from_pnts(&p0, q).dot(&u1) >= off - eps
    };
    let inside2 = |q: &GpPnt| {
        let t = (GpVec::from_pnts(&p0, q).dot(&axis_dir) / edge_len).clamp(0.0, 1.0);
        let off = radius_at(spec, t) * cot_half;
        GpVec::from_pnts(&p0, q).dot(&u2) >= off - eps
    };

    // End faces: faces (other than the two adjacent) containing p0 / p1.
    let all_faces = faces_of(solid);
    let end0: Vec<Face> = all_faces
        .iter()
        .filter(|f| {
            face_contains_point(f, &p0)
                && !is_same(&f.0, &f1.0)
                && !is_same(&f.0, &f2.0)
        })
        .cloned()
        .collect();
    let end1: Vec<Face> = all_faces
        .iter()
        .filter(|f| {
            face_contains_point(f, &p1)
                && !is_same(&f.0, &f1.0)
                && !is_same(&f.0, &f2.0)
        })
        .cloned()
        .collect();
    if end0.len() != 1 || end1.len() != 1 {
        return Err(format!(
            "fillet_var: expected one end face per edge endpoint, got {} and {}",
            end0.len(),
            end1.len()
        ));
    }

    let last = samples.len() - 1;
    let t1_poly: Vec<GpPnt> = samples.iter().map(|s| s.t1).collect();
    let t2_poly: Vec<GpPnt> = samples.iter().map(|s| s.t2).collect();

    let mut cache = EdgeCache::new();
    let trimmed1 = rebuild_var_face(f1, &inside1, &samples[0].t1, &samples[last].t1, &t1_poly, tol, &mut cache)?;
    let trimmed2 = rebuild_var_face(f2, &inside2, &samples[0].t2, &samples[last].t2, &t2_poly, tol, &mut cache)?;
    let end0_face = rebuild_end_face(
        &end0[0],
        &p0,
        &samples[0].t1,
        &samples[0].t2,
        &samples[0].center,
        &axis_dir,
        samples[0].radius,
        &n1,
        &mut cache,
    )?;
    let end1_face = rebuild_end_face(
        &end1[0],
        &p1,
        &samples[last].t1,
        &samples[last].t2,
        &samples[last].center,
        &axis_dir,
        samples[last].radius,
        &n1,
        &mut cache,
    )?;
    let blend_face = build_var_blend_face(&samples, &axis_dir, &mut cache)?;

    // Assemble the new shell: keep every face except the replaced ones.
    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if is_same(&f.0, &f1.0)
            || is_same(&f.0, &f2.0)
            || is_same(&f.0, &end0[0].0)
            || is_same(&f.0, &end1[0].0)
        {
            continue;
        }
        faces.push(f);
    }
    faces.push(trimmed1);
    faces.push(trimmed2);
    faces.push(end0_face);
    faces.push(end1_face);
    faces.push(blend_face);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Apply variable-radius fillets to a chain of edges sequentially, where
/// consecutive edges share a vertex. `specs.len()` must equal
/// `edge_indices.len()`. At a shared vertex the two fillets meet only when
/// `specs[k].r_start` equals `specs[k - 1].r_end` — the caller encodes that
/// continuity; this function does not rewrite a mismatch.
///
/// Each edge is filleted on the result of the previous step. Edges are looked
/// up by their geometric endpoints in the current solid, so a chain whose
/// edges do not interfere (no shared face is rebuilt twice) behaves exactly
/// like `fillet_edge_chain`. When an earlier fillet already replaced an edge
/// (a chain whose shared vertex is consumed by the first fillet), the later
/// edge cannot be re-found and the call errors rather than silently filletting
/// the wrong geometry.
pub fn fillet_edge_var_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    specs: &[VarFilletSpec],
    tol: f64,
) -> Result<TopoShape, String> {
    if edge_indices.len() != specs.len() {
        return Err(
            "fillet_edge_var_chain: specs.len() must equal edge_indices.len()".to_string()
        );
    }
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for (k, &i) in edge_indices.iter().enumerate() {
        let edge = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_edge_var_chain: edge index {i} out of range"))?;
        let (p0, p1) =
            BRepTool::edge_vertices(edge).ok_or("fillet_edge_var_chain: edge has no curve")?;
        let target = find_edge_by_endpoints(&current, &p0, &p1);
        let e = target.ok_or_else(|| {
            format!(
                "fillet_edge_var_chain: edge {i} ({p0:?} → {p1:?}) was consumed by an earlier \
                 fillet and can no longer be re-found; chain a non-adjacent edge or supply a \
                 single-pass chain"
            )
        })?;
        current = fillet_edge_var(&current, &e, &specs[k], tol)?;
    }
    Ok(current)
}

/// Constant-radius fillet along a chain of edges (the ball rolls continuously
/// around a shared vertex). The Phase-6 `fillet_edge::fillet_edge_chain`
/// already implements sequential constant-radius chains; this is a thin wrapper
/// that accepts a tolerance (used only for API symmetry).
pub fn fillet_edge_chain_smooth(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
    _tol: f64,
) -> Result<TopoShape, String> {
    crate::fillet_edge::fillet_edge_chain(solid, edge_indices, radius)
}

/// The `(t, R(t))` pairs actually used by a variable-radius fillet of `edge`,
/// sampled uniformly at `samples` points on `[0, 1]`.
pub fn fillet_radius_profile(
    solid: &TopoShape,
    edge: &Edge,
    spec: &VarFilletSpec,
    samples: usize,
) -> Vec<(f64, f64)> {
    let _ = (solid, edge);
    let n = samples.max(2);
    (0..n).map(|i| {
        let t = i as f64 / (n - 1) as f64;
        (t, radius_at(spec, t))
    }).collect()
}

/// Find an edge of `shape` whose endpoints coincide with `p0`/`p1` (either
/// order). Used by the chain builder to re-locate an edge after a previous
/// fillet rebuilt the solid.
pub(super) fn find_edge_by_endpoints(shape: &TopoShape, p0: &GpPnt, p1: &GpPnt) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        (a.distance(p0) < 1e-6 && b.distance(p1) < 1e-6)
            || (a.distance(p1) < 1e-6 && b.distance(p0) < 1e-6)
    })
}

// ===========================================================================
// Rolling-ball corner patch (shared-vertex blend)
//
// Two consecutive edge fillets meeting at a box corner leave a gap between the
// two blend surfaces and the three corner faces. The gap is closed by a
// spherical patch centred at the rolling-ball corner centre
// `C = v + R·(d0 + d1 + d2)` — the point at distance `R` from each of the
// three faces. A sphere of radius `R·√2` centred at `C` passes through the
// three fillet tangency points, and its boundary is five circular arcs:
//
//   · one arc in each of the three corner-face planes,
//   · one arc shared with each of the two trimmed blend faces.
//
// Each edge fillet is trimmed to start at distance `2R` from the corner, where
// its end arc coincides with the sphere-boundary arc. The result is a closed
// shell with one added face per shared vertex. The blend faces of a corner
// chain are built at the corner radius `R` (the maximum incident radius), so
// the shared tangency geometry is exactly consistent.
// ===========================================================================

/// The three unit directions of the edges meeting at `vertex`, pointing away
/// from the vertex (into the material for a convex corner).
pub(super) fn corner_in_dirs(solid: &TopoShape, vertex: &Vertex) -> Result<[GpVec; 3], String> {
    let p = vertex_position(vertex);
    let es = edges_of(solid);
    let mut dirs: Vec<GpVec> = Vec::new();
    for e in &es {
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            let pa = vertex_position(&va);
            let pb = vertex_position(&vb);
            if pa.distance(&p) < 1e-9 || pb.distance(&p) < 1e-9 {
                let other = if pa.distance(&p) < 1e-9 { pb } else { pa };
                let v = GpVec::from_pnts(&p, &other);
                if v.magnitude() > 1e-12 {
                    dirs.push(v.normalized());
                }
            }
        }
    }
    if dirs.len() != 3 {
        return Err(format!(
            "fillet_corner: vertex has {} incident edges (expected 3)",
            dirs.len()
        ));
    }
    Ok([dirs[0], dirs[1], dirs[2]])
}

/// The rolling-ball corner centre: the point at distance `radius` from each of
/// the three faces meeting at `vertex` (for a box corner, `v + R·Σdirᵢ`).
pub fn corner_center(solid: &TopoShape, vertex: &Vertex, radius: f64, _tol: f64) -> Option<GpPnt> {
    if radius <= 0.0 {
        return None;
    }
    let dirs = corner_in_dirs(solid, vertex).ok()?;
    let p = vertex_position(vertex);
    let s = dirs[0].added(&dirs[1]).added(&dirs[2]);
    Some(p.translated_vec(&s.multiplied_scalar(radius)))
}

/// The effective corner radius at `vertex`: the maximum of the incident-edge
/// fillet radii at that vertex (`r_start` when the vertex is the edge's start,
/// `r_end` otherwise). `specs` is aligned with `edges_of(solid)`; edges without
/// a matching spec are ignored.
pub fn corner_patch_radius_at_vertex(
    solid: &TopoShape,
    vertex: &Vertex,
    specs: &[VarFilletSpec],
    tol: f64,
) -> f64 {
    let p = vertex_position(vertex);
    let es = edges_of(solid);
    let mut rmax = 0.0f64;
    for (i, e) in es.iter().enumerate() {
        let spec = match specs.get(i) {
            Some(s) => s,
            None => continue,
        };
        let (a, b) = edge_vertices(e);
        let (Some(va), Some(vb)) = (a, b) else { continue };
        let pa = vertex_position(&va);
        let pb = vertex_position(&vb);
        let r = if pa.distance(&p) < tol {
            spec.r_start
        } else if pb.distance(&p) < tol {
            spec.r_end
        } else {
            continue;
        };
        rmax = rmax.max(r);
    }
    rmax
}

/// The unit direction of `edge` away from point `p` (one of its endpoints).
pub(super) fn in_edge_dir(edge: &Edge, p: &GpPnt) -> Result<GpVec, String> {
    let (a, b) = edge_vertices(edge);
    let pa = vertex_position(&a.ok_or("fillet_corner: edge has no start vertex")?);
    let pb = vertex_position(&b.ok_or("fillet_corner: edge has no end vertex")?);
    let other = if pa.distance(p) < 1e-9 { pb } else { pa };
    let v = GpVec::from_pnts(p, &other);
    if v.magnitude() < 1e-12 {
        return Err("fillet_corner: degenerate incident edge".to_string());
    }
    Ok(v.normalized())
}

/// Geometry of a shared-vertex corner where two consecutive filleted edges meet.
pub(super) struct CornerGeom {
    pub(super) vertex: GpPnt,
    pub(super) radius: f64,
    /// Rolling-ball corner centre (the corner-patch sphere centre).
    pub(super) center: GpPnt,
    /// Direction of the first filleted edge (e_a), from the corner into material.
    pub(super) d0: GpVec,
    /// Direction of the second filleted edge (e_b).
    pub(super) d1: GpVec,
    /// Direction of the third (unfilleted) edge.
    pub(super) d2: GpVec,
}

impl CornerGeom {
    pub(super) fn from_edges(
        solid: &TopoShape,
        vertex: &Vertex,
        edge_a: &Edge,
        edge_b: &Edge,
        radius: f64,
    ) -> Result<Self, String> {
        let dirs = corner_in_dirs(solid, vertex)?;
        let p = vertex_position(vertex);
        let da = in_edge_dir(edge_a, &p)?;
        let db = in_edge_dir(edge_b, &p)?;
        let mut d0 = GpVec::zero();
        let mut d1 = GpVec::zero();
        let mut d2 = GpVec::zero();
        let mut found = [false; 3];
        for d in &dirs {
            if d.xyz().crossed(&da.xyz()).modulus() < 1e-6 {
                d0 = *d;
                found[0] = true;
            } else if d.xyz().crossed(&db.xyz()).modulus() < 1e-6 {
                d1 = *d;
                found[1] = true;
            } else {
                d2 = *d;
                found[2] = true;
            }
        }
        if !found.iter().all(|f| *f) {
            return Err("fillet_corner: cannot classify the corner edges".to_string());
        }
        let center = p.translated_vec(&d0.added(&d1).added(&d2).multiplied_scalar(radius));
        Ok(CornerGeom { vertex: p, radius, center, d0, d1, d2 })
    }

    pub(super) fn sphere_radius(&self) -> f64 {
        SQRT_2 * self.radius
    }

    /// A point of the corner: `v + R·(a·d0 + b·d1 + c·d2)`.
    pub(super) fn pt(&self, a: f64, b: f64, c: f64) -> GpPnt {
        let v0 = self.d0.multiplied_scalar(a * self.radius);
        let v1 = self.d1.multiplied_scalar(b * self.radius);
        let v2 = self.d2.multiplied_scalar(c * self.radius);
        self.vertex.translated_vec(&v0.added(&v1).added(&v2))
    }

    /// The five boundary arcs of the corner patch, each `(center, normal,
    /// radius, a, b)`. In order they form a closed loop: face ⊥d2, edge b,
    /// face ⊥d0, face ⊥d1, edge a.
    pub(super) fn arcs(&self) -> Vec<(GpPnt, GpVec, f64, GpPnt, GpPnt)> {
        let r = self.radius;
        vec![
            // Arc in the face perpendicular to d2 (the face shared by both edges).
            (self.pt(1.0, 1.0, 0.0), self.d2, r, self.pt(2.0, 1.0, 0.0), self.pt(1.0, 2.0, 0.0)),
            // Arc with edge b (the second filleted edge).
            (self.pt(1.0, 2.0, 1.0), self.d1, r, self.pt(1.0, 2.0, 0.0), self.pt(0.0, 2.0, 1.0)),
            // Arc in the face perpendicular to d0 (the face edge a is not on).
            (self.pt(0.0, 1.0, 1.0), self.d0.multiplied_scalar(-1.0), r, self.pt(0.0, 2.0, 1.0), self.pt(0.0, 0.0, 1.0)),
            // Arc in the face perpendicular to d1 (the face edge b is not on).
            (self.pt(1.0, 0.0, 1.0), self.d1, r, self.pt(0.0, 0.0, 1.0), self.pt(2.0, 0.0, 1.0)),
            // Arc with edge a (the first filleted edge).
            (self.pt(2.0, 1.0, 1.0), self.d0, r, self.pt(2.0, 0.0, 1.0), self.pt(2.0, 1.0, 0.0)),
        ]
    }

    /// The corner arc that lives on the face whose outward normal is `n_out`.
    pub(super) fn arc_on_face(&self, n_out: &GpVec) -> Option<(GpPnt, GpVec, f64, GpPnt, GpPnt)> {
        let arcs = self.arcs();
        let idx = if n_out.xyz().crossed(&self.d2.xyz()).modulus() < 1e-6 {
            Some(0)
        } else if n_out.xyz().crossed(&self.d0.xyz()).modulus() < 1e-6 {
            Some(2)
        } else if n_out.xyz().crossed(&self.d1.xyz()).modulus() < 1e-6 {
            Some(3)
        } else {
            None
        };
        idx.map(|i| arcs[i].clone())
    }
}

/// Build the spherical corner-patch face for `corner`, with its five boundary
/// arcs shared with the corner faces and the two trimmed blend faces.
pub(super) fn build_corner_patch_face(corner: &CornerGeom, cache: &mut EdgeCache) -> Result<Face, String> {
    let arcs = corner.arcs();
    let mut edges = Vec::with_capacity(arcs.len());
    for (center, normal, r, a, b) in arcs {
        let nd = GpDir::from_vec(&normal).map_err(|e| e.to_string())?;
        edges.push(cache.arc(&center, &nd, r, &a, &b)?);
    }
    let b = TopoBuilder::new();
    let wire = b.make_wire(&edges);
    let mut sph = GpSphere::new(GpAx3::standard(), corner.sphere_radius()).map_err(|e| e.to_string())?;
    sph.set_location(corner.center);
    Ok(b.make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// Parameters `t ∈ [0,1]` where the segment `p→q` crosses the circle centred at
/// `c` of radius `r` (the disk boundary).
pub(super) fn line_circle_hits(p: &GpPnt, q: &GpPnt, c: &GpPnt, r: f64) -> Vec<f64> {
    let d = GpVec::from_pnts(p, q);
    let f = GpVec::from_pnts(p, c); // f = c − p
    let a = d.dot(&d);
    if a < 1e-30 {
        return Vec::new();
    }
    // |p + t·d − c|² = r²  ⇒  a·t² − 2(d·f)·t + (|f|² − r²) = 0.
    let b = -2.0 * f.dot(&d);
    let cc = f.dot(&f) - r * r;
    let disc = b * b - 4.0 * a * cc;
    if disc < 0.0 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    let mut ts: Vec<f64> = Vec::new();
    for t in [(-b + sq) / (2.0 * a), (-b - sq) / (2.0 * a)] {
        if t > -1e-9 && t < 1.0 + 1e-9 {
            ts.push(t.clamp(0.0, 1.0));
        }
    }
    ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ts
}

/// Remove the disk (centre `c`, radius `r`) from the polygon, keeping the
/// outside portion. The disk cap is replaced by a straight chord; the returned
/// polygon has that chord as one boundary edge, to be replaced by the shared
/// corner arc.
pub(super) fn cut_disk_poly(poly: &[GpPnt], c: &GpPnt, r: f64) -> Result<Vec<GpPnt>, String> {
    let n = poly.len();
    if n < 3 {
        return Err("cut_disk: degenerate polygon".to_string());
    }
    let mut out: Vec<GpPnt> = Vec::new();
    for i in 0..n {
        let p = &poly[i];
        let q = &poly[(i + 1) % n];
        let hits = line_circle_hits(p, q, c, r);
        let mut ts: Vec<f64> = vec![0.0];
        for t in &hits {
            if *t > 1e-9 && *t < 1.0 - 1e-9 {
                ts.push(*t);
            }
        }
        ts.push(1.0);
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for w in ts.windows(2) {
            let (t0, t1) = (w[0], w[1]);
            if t1 - t0 < 1e-12 {
                continue;
            }
            let tm = 0.5 * (t0 + t1);
            let mid = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(tm));
            if mid.distance(c) >= r - 1e-9 {
                let a0 = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(t0));
                let a1 = p.translated_vec(&GpVec::from_pnts(p, q).multiplied_scalar(t1));
                if out.last().map_or(true, |last| last.distance(&a0) > 1e-6) {
                    out.push(a0);
                }
                if a0.distance(&a1) > 1e-6 {
                    out.push(a1);
                }
            }
        }
    }
    // Close the loop: drop a trailing point equal to the first.
    while out.len() >= 2 && out[0].distance(out.last().unwrap()) < 1e-9 {
        out.pop();
    }
    if out.len() < 3 {
        return Err("cut_disk: the disk removed the whole polygon".to_string());
    }
    Ok(out)
}

/// Rebuild a corner face: clip its polygon by the half-plane tangency trims,
/// then cut out each corner disk and replace the resulting chord with the
/// shared circular-arc edge. `arcs` holds the corner arcs that lie on this
/// face, one per adjacent corner. `polylines` holds the sampled tangency
/// polylines of the adjacent blend faces, which replace the straight tangency
/// edges so the corner face shares those edges with the blend faces.
pub(super) fn rebuild_corner_face(
    face: &Face,
    trims: &[(GpPnt, GpVec, f64)],
    arcs: &[(GpPnt, GpVec, f64, GpPnt, GpPnt)],
    polylines: &[(GpPnt, GpPnt, Vec<GpPnt>)],
    tol: f64,
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let poly = register_face_edges(face, cache)?;
    let eps = 1e-9;
    let mut clipped = poly;
    for (p0, u, off) in trims {
        let inside = |q: &GpPnt| GpVec::from_pnts(p0, q).dot(u) >= off - eps;
        clipped = clip_polygon(&clipped, &inside);
        if clipped.len() < 3 {
            return Err("fillet_corner: clipping a corner face left no retained region".to_string());
        }
    }
    let mut arc_edges: Vec<(GpPnt, GpPnt, Edge)> = Vec::new();
    for (center, normal, radius, a, b) in arcs {
        let cut = cut_disk_poly(&clipped, center, *radius)?;
        if cut.len() < 3 {
            return Err("fillet_corner: cutting the corner disk left no region".to_string());
        }
        let nd = GpDir::from_vec(normal).map_err(|e| e.to_string())?;
        let arc_edge = cache.arc(center, &nd, *radius, a, b)?;
        arc_edges.push((*a, *b, arc_edge));
        clipped = cut;
    }
    let match_tol = (10.0 * tol).max(1e-4).min(1e-3);
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    let m = clipped.len();
    for i in 0..m {
        let a = &clipped[i];
        let c = &clipped[(i + 1) % m];
        if let Some((_, _, e)) = arc_edges.iter().find(|(ta, tb, _)| {
            (a.distance(ta) < match_tol && c.distance(tb) < match_tol)
                || (a.distance(tb) < match_tol && c.distance(ta) < match_tol)
        }) {
            edges.push(e.clone());
        } else if let Some((_, _, poly)) = polylines.iter().find(|(f, l, _)| {
            (a.distance(f) < match_tol && c.distance(l) < match_tol)
                || (a.distance(l) < match_tol && c.distance(f) < match_tol)
        }) {
            if a.distance(&poly[0]) < match_tol {
                for w in poly.windows(2) {
                    edges.push(cache.seg(&w[0], &w[1]));
                }
            } else {
                for w in poly.windows(2).rev() {
                    edges.push(cache.seg(&w[1], &w[0]));
                }
            }
        } else if a.distance(c) > 1e-6 {
            edges.push(cache.seg(a, c));
        }
    }
    let wire = b.make_wire(&edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_corner: corner face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// A single end-arc replacement: the corner vertex `corner` of an end face is
/// replaced by the blend arc `t1 → t2` (shared with the blend face).
#[derive(Clone)]
pub(super) struct EndArcReplacement {
    pub(super) corner: GpPnt,
    pub(super) t1: GpPnt,
    pub(super) t2: GpPnt,
    pub(super) center: GpPnt,
    pub(super) axis: GpVec,
    pub(super) radius: f64,
    /// Outward normal of adjacent face 1 (used to classify the incident edges).
    pub(super) n1: GpVec,
}

/// Rebuild a face by replacing several corner vertices with their end arcs.
/// Used for far end faces (perpendicular to an edge) that may carry arcs from
/// more than one edge (e.g. the far +X face of two opposite X edges).
pub(super) fn rebuild_face_with_end_arcs(
    face: &Face,
    reps: &[EndArcReplacement],
    cache: &mut EdgeCache,
) -> Result<Face, String> {
    let mut poly = register_face_edges(face, cache)?;
    let mut arc_pairs: Vec<(GpPnt, GpPnt, Edge)> = Vec::new();
    for rep in reps {
        let n = poly.len();
        let idx = poly
            .iter()
            .position(|q| q.distance(&rep.corner) < 1e-9)
            .ok_or("fillet_corner: end face does not contain the corner vertex")?;
        let prev = poly[(idx + n - 1) % n];
        let next = poly[(idx + 1) % n];
        let dir_prev = GpVec::from_pnts(&prev, &rep.corner).normalized();
        let dir_next = GpVec::from_pnts(&rep.corner, &next).normalized();
        let prev_is_face1 = dir_prev.dot(&rep.n1).abs() < 1e-6;
        let next_is_face1 = dir_next.dot(&rep.n1).abs() < 1e-6;
        let (t_first, t_last) = if prev_is_face1 && !next_is_face1 {
            (rep.t1, rep.t2)
        } else if next_is_face1 && !prev_is_face1 {
            (rep.t2, rep.t1)
        } else {
            return Err("fillet_corner: cannot classify end-face incident edges".to_string());
        };
        let zd = GpDir::from_vec(&rep.axis).map_err(|e| e.to_string())?;
        let arc_edge = cache.arc(&rep.center, &zd, rep.radius, &rep.t1, &rep.t2)?;
        arc_pairs.push((rep.t1, rep.t2, arc_edge));
        poly[idx] = t_first;
        poly.insert(idx + 1, t_last);
    }

    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    let m = poly.len();
    for i in 0..m {
        let a = &poly[i];
        let c = &poly[(i + 1) % m];
        if let Some((_, _, e)) = arc_pairs.iter().find(|(ta, tb, _)| {
            (a.distance(ta) < 1e-6 && c.distance(tb) < 1e-6)
                || (a.distance(tb) < 1e-6 && c.distance(ta) < 1e-6)
        }) {
            edges.push(e.clone());
        } else if a.distance(c) > 1e-6 {
            edges.push(cache.seg(a, c));
        }
    }
    let wire = b.make_wire(&edges);
    let surf = BRepTool::face_surface(face).ok_or("fillet_corner: end face has no surface")?;
    Ok(b.make_face(surf, &[wire]))
}

/// Per-edge data for a corner-chain run.
pub(super) struct EdgeInfo {
    pub(super) edge: Edge,
    pub(super) p0: GpPnt,
    pub(super) p1: GpPnt,
    pub(super) f1: Face,
    pub(super) f2: Face,
    pub(super) n1: GpVec,
    pub(super) n2: GpVec,
    pub(super) axis: GpVec,
    pub(super) len: f64,
    pub(super) spec: VarFilletSpec,
}

pub(super) fn find_shared_vertex(a0: &GpPnt, a1: &GpPnt, b0: &GpPnt, b1: &GpPnt) -> Option<GpPnt> {
    for p in [a0, a1] {
        if p.distance(b0) < 1e-9 || p.distance(b1) < 1e-9 {
            return Some(*p);
        }
    }
    None
}

pub(super) fn radius_at_endpoint(spec: &VarFilletSpec, p0: &GpPnt, p1: &GpPnt, v: &GpPnt) -> f64 {
    if p0.distance(v) < 1e-9 {
        spec.r_start
    } else if p1.distance(v) < 1e-9 {
        spec.r_end
    } else {
        0.5 * (spec.r_start + spec.r_end)
    }
}

pub(super) fn find_vertex_at(solid: &TopoShape, p: &GpPnt) -> Option<Vertex> {
    crate::topo_tools_full::vertices_of(solid)
        .into_iter()
        .find(|v| BRepTool::vertex_point(v).distance(p) < 1e-6)
}
