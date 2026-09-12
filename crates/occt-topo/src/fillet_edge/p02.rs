use super::prelude::*;
use super::*;

/// Fillet several edges of `solid` sequentially (each step feeds the next).
/// `edge_indices` indexes into `edges_of(solid)`. Edges are filleted in the
/// given order; the first failure aborts with an error.
pub fn fillet_edge_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
) -> Result<TopoShape, String> {
    let mut current = solid.clone();
    for &i in edge_indices {
        let es = edges_of(&current);
        let e = es
            .get(i)
            .ok_or_else(|| format!("fillet_edge_chain: edge index {i} out of range"))?;
        current = fillet_edge(&current, e, radius)?;
    }
    Ok(current)
}

/// Build a spherical-octant blend face for a box corner. `p` is the corner
/// point, `dirs` the three unit edge directions from `p` into the material,
/// `radius` the blend radius. The sphere is centred at `p` (an octant of a
/// sphere centred on the corner), which is the smallest single patch that
/// closes the trimmed faces.
///
/// ponytail: a rolling-ball corner patch that is tangent to all three faces
/// (a sphere centred at p + R·Σdirᵢ) cannot be closed by a single spherical
/// face — it meets each plane at a single point and needs the three edge
/// fillet cylinders. The corner here is therefore a spherical octant cut,
/// which forms a closed shell with exactly one added face.
///
/// The octant boundary is three quarter-circle arcs, each shared with one of
/// the three trimmed planar faces: `t_i = p + R·dirs[i]` and the arc between
/// `t_j`/`t_k` lies in the plane of the face perpendicular to `dirs[i]`.
pub(super) fn build_corner_blend(p: &GpPnt, dirs: &[GpVec; 3], radius: f64, cache: &mut EdgeCache) -> Result<Face, String> {
    let d = [
        dirs[0].normalized(),
        dirs[1].normalized(),
        dirs[2].normalized(),
    ];
    let t = [
        p.translated_vec(&d[0].multiplied_scalar(radius)),
        p.translated_vec(&d[1].multiplied_scalar(radius)),
        p.translated_vec(&d[2].multiplied_scalar(radius)),
    ];
    let nd0 = GpDir::from_xyz(&d[0].xyz()).map_err(|e| e.to_string())?;
    let nd1 = GpDir::from_xyz(&d[1].xyz()).map_err(|e| e.to_string())?;
    let nd2 = GpDir::from_xyz(&d[2].xyz()).map_err(|e| e.to_string())?;
    let a01 = cache.arc(p, &nd2, radius, &t[0], &t[1])?; // in plane of face ⊥ d2
    let a12 = cache.arc(p, &nd0, radius, &t[1], &t[2])?; // in plane of face ⊥ d0
    let a20 = cache.arc(p, &nd1, radius, &t[2], &t[0])?; // in plane of face ⊥ d1
    let wire = TopoBuilder::new().make_wire(&[a01, a12, a20]);
    let mut sph = GpSphere::new(GpAx3::standard(), radius).map_err(|e| e.to_string())?;
    sph.set_location(*p);
    Ok(TopoBuilder::new().make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// Fillet the three edges meeting at a box corner: replaces the corner with a
/// spherical-quadrant blend and trims the three adjacent planar faces.
/// `corner_vertex` must be a vertex of `solid` with exactly three incident
/// edges (a box corner). Returns a closed solid with one additional face.
pub fn fillet_corner_solid(
    solid: &TopoShape,
    corner_vertex: &Vertex,
    radius: f64,
) -> Result<TopoShape, String> {
    FilletSpec { radius }.check()?;
    let p = BRepTool::vertex_point(corner_vertex);

    // The three incident edges and their unit directions into the material.
    let es = edges_of(solid);
    let mut incident: Vec<Edge> = Vec::new();
    for e in &es {
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            let pa = vertex_position(&va);
            let pb = vertex_position(&vb);
            if pa.distance(&p) < 1e-9 || pb.distance(&p) < 1e-9 {
                incident.push(e.clone());
            }
        }
    }
    if incident.len() != 3 {
        return Err(format!(
            "fillet_corner_solid: vertex has {} incident edges (expected 3)",
            incident.len()
        ));
    }
    let mut dirs: Vec<GpVec> = Vec::with_capacity(3);
    for e in &incident {
        let (a, b) = edge_vertices(e);
        let (pa, pb) = (
            vertex_position(&a.ok_or("fillet: edge has no start vertex")?),
            vertex_position(&b.ok_or("fillet: edge has no end vertex")?),
        );
        let other = if pa.distance(&p) < 1e-9 { pb } else { pa };
        dirs.push(GpVec::from_pnts(&p, &other).normalized());
    }
    let dirs = [dirs[0], dirs[1], dirs[2]];

    // The three faces containing the corner.
    let all_faces = faces_of(solid);
    let corner_faces: Vec<Face> = all_faces
        .iter()
        .filter(|f| face_contains_point(f, &p))
        .cloned()
        .collect();
    if corner_faces.len() != 3 {
        return Err(format!(
            "fillet_corner_solid: corner touches {} faces (expected 3)",
            corner_faces.len()
        ));
    }

    // For each corner face, find its outward normal and its polygon. The two
    // in-plane directions from the corner are the two edges on the face.
    let mut cache = EdgeCache::new();
    let mut rebuilt: Vec<Face> = Vec::new();
    let mut idx_by_normal: Vec<(usize, GpVec)> = Vec::new();
    for f in &corner_faces {
        let n = face_outward_normal(f)?;
        // Match this face to the axis direction it is perpendicular to.
        let mut match_idx = usize::MAX;
        for (k, d) in dirs.iter().enumerate() {
            if n.xyz().crossed(&d.xyz()).modulus() < 1e-6 {
                match_idx = k;
                break;
            }
        }
        if match_idx == usize::MAX {
            return Err("fillet_corner_solid: cannot match corner face to axis".to_string());
        }
        idx_by_normal.push((match_idx, n));
    }
    // Build the retained polygon for each face: outside the sphere centred at p
    // with radius R (in the two in-plane coordinates).
    for (f, (k, n)) in corner_faces.iter().zip(idx_by_normal.iter()) {
        let poly = register_face_edges(f, &mut cache)?;
        let d_j = dirs[(k + 1) % 3];
        let d_k = dirs[(k + 2) % 3];
        // Retained: the point is outside the radius-R disk in the (d_j, d_k)
        // plane (i.e. its squared distance from p along the two in-plane axes
        // is at least R²).
        let inside = |q: &GpPnt| {
            let v = GpVec::from_pnts(&p, q);
            let a = v.dot(&d_j);
            let b = v.dot(&d_k);
            a * a + b * b >= radius * radius - 1e-9
        };
        let clipped = clip_polygon(&poly, &inside);
        if clipped.len() < 3 {
            return Err("fillet_corner_solid: clipping left no face".to_string());
        }
        // The clipped polygon has a chamfer corner at the two tangency points;
        // replace it with the spherical arc shared with the blend face.
        let wire = build_corner_face_wire(&p, &dirs, radius, &clipped, &mut cache, &n)?;
        let surf = BRepTool::face_surface(f).ok_or("fillet: face has no surface")?;
        rebuilt.push(cache.b.make_face(surf, &[wire]));
    }

    let blend_face = build_corner_blend(&p, &dirs, radius, &mut cache)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in all_faces {
        if corner_faces.iter().any(|cf| is_same(&f.0, &cf.0)) {
            continue;
        }
        faces.push(f);
    }
    faces.extend(rebuilt);
    faces.push(blend_face);

    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Build the wire of a corner-trimmed planar face: the clipped polygon with its
/// straight chamfer edge (between the two tangency points) replaced by the
/// spherical arc shared with the corner blend face.
///
/// `n` is the face's outward normal; the face is the one perpendicular to the
/// axis direction `dirs[k]` that is parallel to `n`. `dirs` are the three
/// in-wedge directions from the corner.
pub(super) fn build_corner_face_wire(
    p: &GpPnt,
    dirs: &[GpVec; 3],
    radius: f64,
    clipped: &[GpPnt],
    cache: &mut EdgeCache,
    n: &GpVec,
) -> Result<Wire, String> {
    let k = (0..3)
        .find(|&i| dirs[i].xyz().crossed(&n.xyz()).modulus() < 1e-6)
        .ok_or("fillet: cannot classify corner face")?;
    let d = [
        dirs[0].normalized(),
        dirs[1].normalized(),
        dirs[2].normalized(),
    ];
    let t_a = p.translated_vec(&d[(k + 1) % 3].multiplied_scalar(radius));
    let t_b = p.translated_vec(&d[(k + 2) % 3].multiplied_scalar(radius));
    // The spherical arc shared with the blend face: in this face's plane
    // (normal d[k]), from t_a to t_b.
    let nd = GpDir::from_xyz(&d[k].xyz()).map_err(|e| e.to_string())?;
    let arc = cache.arc(p, &nd, radius, &t_a, &t_b)?;

    // Find the chamfer edge of the clipped polygon: the straight segment
    // connecting t_a and t_b (either orientation).
    let m = clipped.len();
    let mut chamfer = None;
    for i in 0..m {
        let a = &clipped[i];
        let c = &clipped[(i + 1) % m];
        let is_ab = (a.distance(&t_a) < 1e-6 && c.distance(&t_b) < 1e-6)
            || (a.distance(&t_b) < 1e-6 && c.distance(&t_a) < 1e-6);
        if is_ab {
            chamfer = Some(i);
            break;
        }
    }
    let ci = chamfer.ok_or("fillet: corner polygon has no chamfer edge")?;

    // Replace the chamfer edge with the arc (which shares its endpoint pair,
    // so the cache returns the same instance the blend face uses).
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    for i in 0..m {
        if i == ci {
            edges.push(arc.clone());
        } else {
            edges.push(cache.seg(&clipped[i], &clipped[(i + 1) % m]));
        }
    }
    Ok(b.make_wire(&edges))
}

/// Classify a surface, adding cylinder detection on top of
/// `brep_surface::classify_surface` (which only distinguishes Plane / Sphere /
/// Other). A cylinder is recognised by constant distance from an axis line.
pub fn classify_surface_full(s: &dyn Surface) -> SurfaceKind {
    let k = classify_surface(s);
    if k != SurfaceKind::Other {
        return k;
    }
    if cylinder_radius(s).is_some() {
        SurfaceKind::Cylinder
    } else {
        SurfaceKind::Other
    }
}

/// If `s` is a cylinder, return its radius. All sampled surface normals must be
/// perpendicular to one axis direction, and all sampled points equidistant
/// from that axis line.
pub fn cylinder_radius(s: &dyn Surface) -> Option<f64> {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let nu = 6;
    let nv = 6;
    let hu = ((u1 - u0) * 1e-4).max(1e-7);
    let hv = ((v1 - v0) * 1e-4).max(1e-7);
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut nrm: Vec<GpVec> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            pts.push(s.d0(u, v));
            // Local finite-difference normal: `brep_surface::surface_normal`
            // falls back to an infinite step on unbounded parameter ranges
            // (e.g. the cylinder's infinite v range), so compute it here with
            // the clamped sampling bounds.
            let p0 = s.d0(u, v);
            let du = GpVec::from_pnts(&p0, &s.d0(u + hu, v));
            let dv = GpVec::from_pnts(&p0, &s.d0(u, v + hv));
            let nn = du.xyz().crossed(&dv.xyz());
            let m = nn.modulus();
            nrm.push(if m > 1e-30 { GpVec::from_xyz(&nn.divided(m)) } else { GpVec::zero() });
        }
    }
    // Axis direction: perpendicular to all normals. Pick two non-parallel
    // normals and cross them.
    let mut axis = None;
    for a in 0..nrm.len() {
        for b in (a + 1)..nrm.len() {
            let c = nrm[a].xyz().crossed(&nrm[b].xyz());
            if c.modulus() > 1e-6 {
                axis = Some(GpVec::from_xyz(&c).normalized());
                break;
            }
        }
        if axis.is_some() {
            break;
        }
    }
    let axis = axis?;
    // All normals must be perpendicular to the axis (dot product ≈ 0).
    for nn in &nrm {
        if nn.dot(&axis).abs() > 1e-3 {
            return None;
        }
    }
    // Project points onto the plane perpendicular to the axis; the axis line
    // passes through the 2D circumcenter of three non-collinear projections.
    let (x2, y2) = project_basis(&axis);
    let proj: Vec<(f64, f64)> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&pts[0], q);
            (v.dot(&x2), v.dot(&y2))
        })
        .collect();
    let (mut cx, mut cy) = (f64::NAN, f64::NAN);
    'outer: for a in 0..proj.len() {
        for b in (a + 1)..proj.len() {
            for c in (b + 1)..proj.len() {
                let (pa, pb, pc) = (proj[a], proj[b], proj[c]);
                let d = 2.0 * (pa.0 * (pb.1 - pc.1) + pb.0 * (pc.1 - pa.1) + pc.0 * (pa.1 - pb.1));
                if d.abs() < 1e-12 {
                    continue;
                }
                let ux = ((pa.0 * pa.0 + pa.1 * pa.1) * (pb.1 - pc.1)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pc.1 - pa.1)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pa.1 - pb.1))
                    / d;
                let uy = ((pa.0 * pa.0 + pa.1 * pa.1) * (pc.0 - pb.0)
                    + (pb.0 * pb.0 + pb.1 * pb.1) * (pa.0 - pc.0)
                    + (pc.0 * pc.0 + pc.1 * pc.1) * (pb.0 - pa.0))
                    / d;
                cx = ux;
                cy = uy;
                break 'outer;
            }
        }
    }
    if !cx.is_finite() {
        return None;
    }
    let axis_pt = pts[0].translated_vec(&x2.multiplied_scalar(cx).added(&y2.multiplied_scalar(cy)));
    let mut radii: Vec<f64> = pts
        .iter()
        .map(|q| {
            let v = GpVec::from_pnts(&axis_pt, q);
            let along = v.dot(&axis);
            let perp = v.subtracted(&axis.multiplied_scalar(along));
            perp.magnitude()
        })
        .collect();
    if radii.is_empty() {
        return None;
    }
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = radii[radii.len() / 2];
    if med < 1e-9 {
        return None;
    }
    for r in &radii {
        if (r - med).abs() > 1e-3 * med.abs().max(1.0) {
            return None;
        }
    }
    Some(med)
}

/// Two orthonormal vectors spanning the plane perpendicular to `axis`.
pub(super) fn project_basis(axis: &GpVec) -> (GpVec, GpVec) {
    let ref_v = if axis.x().abs() < 0.9 {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    let x2 = axis.crossed(&ref_v).normalized();
    let y2 = axis.crossed(&x2).normalized();
    (x2, y2)
}

/// Finite, sane sampling bounds for a surface (unbounded ranges clamp to ±1).
pub(super) fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}
