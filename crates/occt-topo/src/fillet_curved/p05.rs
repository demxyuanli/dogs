use super::prelude::*;
use super::*;

/// Fillet one sphere+sphere arc edge of the tri-sphere corner into a torus
/// sector spanning the azimuth interval between the two rolling-ball corner
/// centres `c0` (near the edge's first endpoint) and `c1` (near its last). The
/// two adjacent sphere caps are rebuilt preserving their other boundaries, and
/// the two sector cross-section arcs (at `c0` and `c1`) are returned so the
/// caller can wire the spherical corner patches.
#[allow(clippy::too_many_arguments)]
pub(super) fn fillet_sphere_sector(
    solid: &TopoShape,
    edge: &Edge,
    f1: &Face,
    f2: &Face,
    s1: &SphereInfo,
    s2: &SphereInfo,
    radius: f64,
    c0: &GpPnt,
    c1: &GpPnt,
) -> Result<(TopoShape, Edge, Edge), String> {
    let b = TopoBuilder::new();
    let g = sphere_sphere_blend(s1, s2, radius)?;
    let dhat = g.axis.normalized();
    let (e1, e2) = project_basis(&dhat);

    let az = |p: &GpPnt| point_azimuth(p, &s1.center, &dhat, &e1, &e2);
    let a0 = az(c0);
    let a1 = az(c1);

    // Determine the edge's forward azimuth direction (increasing or decreasing
    // parameter) from the midpoint, then express c0/c1 azimuths in that frame.
    let (p0, _p1) = BRepTool::edge_vertices(edge).ok_or("fillet_curved: edge has no vertices")?;
    let az_p0 = az(&p0);
    let curve = BRepTool::edge_curve(edge).ok_or("fillet_curved: edge has no curve")?;
    let (first, last) = BRepTool::edge_parameters(edge);
    let p_mid = curve.d0(0.5 * (first + last));
    let rel_mid = (az(&p_mid) - az_p0).rem_euclid(2.0 * PI);
    let increases = rel_mid < PI;
    let forward = |a: f64| -> f64 {
        if increases {
            (a - az_p0).rem_euclid(2.0 * PI)
        } else {
            (az_p0 - a).rem_euclid(2.0 * PI)
        }
    };
    // Sector starts at the corner near p0 and ends at the corner near p1.
    let start = az_p0 + if increases { forward(a0) } else { -forward(a0) };
    let end = az_p0 + if increases { forward(a1) } else { -forward(a1) };

    let tang_pt = |ct: &ContactCircleGeom, th: f64| {
        ct.center.translated_vec(
            &e1.multiplied_scalar(ct.radius * th.cos()).added(&e2.multiplied_scalar(ct.radius * th.sin())),
        )
    };
    let pt_a0 = tang_pt(&g.contact_a, start);
    let pt_a1 = tang_pt(&g.contact_a, end);
    let pt_b0 = tang_pt(&g.contact_b, start);
    let pt_b1 = tang_pt(&g.contact_b, end);

    let tang_a = build_contact_arc(&b, &g.contact_a, &dhat, &e1, start, end, &pt_a0, &pt_a1)?;
    let tang_b = build_contact_arc(&b, &g.contact_b, &dhat, &e1, start, end, &pt_b0, &pt_b1)?;
    let cross_start = build_sector_cross_arc(&b, &g, radius, start, &e1, &e2, &dhat)?;
    let cross_end = build_sector_cross_arc(&b, &g, radius, end, &e1, &e2, &dhat)?;

    // Torus sector face.
    let surf = torus_surface(&g)?;
    let sector_wire = b.make_wire(&[tang_a.clone(), cross_end.clone(), tang_b.clone(), cross_start.clone()]);
    let sector_face = b.make_face(surf, &[sector_wire]);

    // Rebuilt adjacent sphere caps.
    let w_a = b.make_wire(&[tang_a.clone()]);
    let rebuilt1 = rebuild_sphere_keep(f1, edge, &w_a, &b)?;
    let w_b = b.make_wire(&[tang_b.clone()]);
    let rebuilt2 = rebuild_sphere_keep(f2, edge, &w_b, &b)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0) || is_same(&f.0, &f2.0) {
            continue;
        }
        faces.push(f);
    }
    faces.push(rebuilt1);
    faces.push(rebuilt2);
    faces.push(sector_face);
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok((solid_out.0, cross_start, cross_end))
}

/// Rolling-ball corner centre for the tri-sphere corner: the point at distance
/// `radius` from every face incident at `pv`, found by Newton iteration. `hint`
/// seeds the iteration toward the exterior of the corner (there are two such
/// points — one near each of the two common vertices — and the hint picks the
/// right one).
pub(super) fn corner_center_ball(solid: &TopoShape, pv: &GpPnt, radius: f64, hint: &GpVec) -> Result<GpPnt, String> {
    let faces: Vec<Face> = faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| {
                    let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
                    a.distance(pv) < 1e-6 || b.distance(pv) < 1e-6
                })
            })
        })
        .collect();
    if faces.len() < 3 {
        return Err(format!(
            "fillet_curved: corner vertex touches {} faces (expected >= 3)",
            faces.len()
        ));
    }
    let mut c = pv.translated_vec(&hint.multiplied_scalar(radius));
    for _ in 0..300 {
        let mut gx = 0.0f64;
        let mut gy = 0.0f64;
        let mut gz = 0.0f64;
        let mut grad = [[0.0; 3]; 3];
        let mut err_sum = 0.0;
        for f in &faces {
            let s = BRepTool::face_surface(f).ok_or("fillet_curved: corner face has no surface")?;
            let (q, n) = closest_point_on_surface(s.as_ref(), &c).ok_or("fillet_curved: cannot project onto a corner face")?;
            let residual = q.distance(&c) - radius;
            err_sum += residual.abs();
            let nxyz = n.xyz();
            let (nx, ny, nz) = (nxyz.x, nxyz.y, nxyz.z);
            gx += residual * nx;
            gy += residual * ny;
            gz += residual * nz;
            let comps = [nx, ny, nz];
            for a in 0..3 {
                for b in 0..3 {
                    grad[a][b] += comps[a] * comps[b];
                }
            }
        }
        if err_sum / (faces.len() as f64) < 1e-7 {
            return Ok(c);
        }
        if let Some(delta) = solve3x3(&grad, &[gx, gy, gz]) {
            let step = GpVec::new(delta[0], delta[1], delta[2]);
            if step.magnitude() < 1e-12 {
                return Ok(c);
            }
            // Gauss-Newton minimises Σ residual², so the step is subtracted
            // (otherwise the iteration converges to the *internal* tangent ball
            // instead of the exterior rolling ball of the fillet).
            c = c.translated_vec(&step.multiplied_scalar(-1.0));
        } else {
            break;
        }
    }
    Ok(c)
}

/// Deterministic corner patch for 3+ curved edges meeting at `corner_vertex`.
///
/// Every edge is filleted as a torus sector (sphere+sphere arc edges of the
/// tri-sphere corner), the adjacent sphere caps become bands bounded by their
/// two tangency arcs, and each shared vertex is closed by a spherical patch of
/// radius `radius` centred at the rolling-ball corner centre, bounded by the
/// three sector cross-section arcs. Unlike the best-effort
/// `fillet_curved_multi_normal`, the result is a closed shell: every edge of
/// the rebuilt solid is referenced by exactly two faces.
pub fn fillet_curved_corner_patch(
    solid: &TopoShape,
    corner_vertex: &crate::shape::Vertex,
    edges: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_curved_corner_patch: radius must be a positive finite value".to_string());
    }
    if edges.len() < 3 {
        return Err(format!(
            "fillet_curved_corner_patch: needs at least 3 edges, got {}",
            edges.len()
        ));
    }
    let original_edges = edges_of(solid);
    let mut e_objs: Vec<Edge> = Vec::new();
    for &i in edges {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_curved_corner_patch: edge index {i} out of range"))?;
        // Use the edge object by index (several tri-sphere edges share the same
        // geometric endpoints, so `find_edge_by_endpoints` would be ambiguous).
        e_objs.push(oe.clone());
    }

    let pv = BRepTool::vertex_point(corner_vertex);

    // The vertices shared by >= 2 of the edges (a tri-sphere corner has two).
    let mut shared: Vec<GpPnt> = vec![pv];
    for e in &e_objs {
        let (a, b) = BRepTool::edge_vertices(e).ok_or("fillet_curved_corner_patch: edge has no vertices")?;
        for p in [a, b] {
            if !shared.iter().any(|q| q.distance(&p) < 1e-6) {
                shared.push(p);
            }
        }
    }
    shared.retain(|q| {
        let count = e_objs
            .iter()
            .filter(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.distance(q) < 1e-6 || b.distance(q) < 1e-6
            })
            .count();
        count >= 2
    });
    if shared.is_empty() {
        return Err("fillet_curved_corner_patch: the given edges share no vertex".to_string());
    }

    // Rolling-ball corner centre for each shared vertex, seeded toward the
    // exterior (away from the other shared vertices).
    let mut centers: Vec<GpPnt> = Vec::new();
    for q in &shared {
        let mut hint = GpVec::zero();
        for r in &shared {
            if r.distance(q) > 1e-6 {
                let d = GpVec::from_pnts(r, q);
                let m = d.magnitude();
                if m > 1e-12 {
                    hint = hint.added(&d.multiplied_scalar(1.0 / m));
                }
            }
        }
        if hint.magnitude() < 1e-12 {
            hint = GpVec::new(1.0, 0.0, 0.0);
        }
        let c = corner_center_ball(solid, q, radius, &hint.normalized())?;
        centers.push(c);
    }

    // Fillet each edge as a torus sector, collecting the cross-section arcs per
    // corner centre for the spherical patches.
    let mut current = solid.clone();
    let mut cross_by_center: Vec<Vec<Edge>> = vec![Vec::new(); shared.len()];
    for (ei, e) in e_objs.iter().enumerate() {
        let adj = adjacent_faces(&current, e);
        if adj.len() != 2 {
            return Err(format!(
                "fillet_curved_corner_patch: edge {ei} is adjacent to {} faces (expected 2)",
                adj.len()
            ));
        }
        let (f1, f2) = (&adj[0], &adj[1]);
        let s1 = BRepTool::face_surface(f1).ok_or("fillet_curved_corner_patch: face 1 has no surface")?;
        let s2 = BRepTool::face_surface(f2).ok_or("fillet_curved_corner_patch: face 2 has no surface")?;
        let (k1, k2) = (classify_surface_analytic(s1.as_ref()), classify_surface_analytic(s2.as_ref()));
        if k1 != SurfaceKind::Sphere || k2 != SurfaceKind::Sphere {
            return Err(format!(
                "fillet_curved_corner_patch: edge {ei} is not a sphere+sphere pair ({k1:?}, {k2:?})"
            ));
        }
        let sp1 = sphere_from_surface(s1.as_ref()).ok_or("fillet_curved_corner_patch: cannot extract sphere 1")?;
        let sp2 = sphere_from_surface(s2.as_ref()).ok_or("fillet_curved_corner_patch: cannot extract sphere 2")?;

        let (p0, p1) = BRepTool::edge_vertices(e).unwrap();
        let idx0 = shared
            .iter()
            .position(|q| q.distance(&p0) < 1e-6)
            .ok_or("fillet_curved_corner_patch: edge endpoint 0 is not a shared vertex")?;
        let idx1 = shared
            .iter()
            .position(|q| q.distance(&p1) < 1e-6)
            .ok_or("fillet_curved_corner_patch: edge endpoint 1 is not a shared vertex")?;
        let (new_solid, cross_start, cross_end) =
            fillet_sphere_sector(&current, e, f1, f2, &sp1, &sp2, radius, &centers[idx0], &centers[idx1])?;
        current = new_solid;
        cross_by_center[idx0].push(cross_start);
        cross_by_center[idx1].push(cross_end);
    }

    // Spherical corner patches.
    let b = TopoBuilder::new();
    let mut faces = faces_of(&current);
    for (i, c) in centers.iter().enumerate() {
        if cross_by_center[i].len() < 3 {
            return Err(format!(
                "fillet_curved_corner_patch: corner centre {i} has only {} cross-section arcs",
                cross_by_center[i].len()
            ));
        }
        let wire = b.make_wire(&cross_by_center[i]);
        let mut sph = GpSphere::new(GpAx3::standard(), radius).map_err(|e| e.to_string())?;
        sph.set_location(*c);
        faces.push(b.make_face(Arc::new(GeomSphere::new(sph)), &[wire]));
    }

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    let _ = tol;
    Ok(solid_out.0)
}
