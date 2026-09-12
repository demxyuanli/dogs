use super::prelude::*;
use super::*;

/// Closest point and outward normal on an analytic surface to a given point.
/// Handles plane / sphere / cylinder / cone exactly; other surfaces fall back
/// to a coarse grid projection (which may return `None`).
pub fn closest_point_on_surface(s: &dyn Surface, p: &GpPnt) -> Option<(GpPnt, GpVec)> {
    match classify_surface_analytic(s) {
        SurfaceKind::Plane => {
            let pln = plane_from_surface(s)?;
            let n = pln.normal.normalized();
            let d = GpVec::from_pnts(&pln.origin, p).dot(&n);
            let q = p.translated_vec(&n.multiplied_scalar(-d));
            Some((q, n))
        }
        SurfaceKind::Sphere => {
            let sph = sphere_from_surface(s)?;
            let v = GpVec::from_pnts(&sph.center, p);
            let m = v.magnitude();
            if m < 1e-30 {
                return None;
            }
            let n = v.normalized();
            let q = sph.center.translated_vec(&n.multiplied_scalar(sph.radius));
            Some((q, n))
        }
        SurfaceKind::Cylinder => {
            let cyl = cylinder_from_surface(s)?;
            let v = GpVec::from_pnts(&cyl.axis_origin, p);
            let along = v.dot(&cyl.axis);
            let radial = v.subtracted(&cyl.axis.multiplied_scalar(along));
            let m = radial.magnitude();
            if m < 1e-12 {
                return None;
            }
            let n = radial.normalized();
            let q = cyl
                .axis_origin
                .translated_vec(&cyl.axis.multiplied_scalar(along))
                .translated_vec(&n.multiplied_scalar(cyl.radius));
            Some((q, n))
        }
        SurfaceKind::Cone => {
            let cone = cone_from_surface(s)?;
            let v = GpVec::from_pnts(&cone.apex, p);
            let h = v.dot(&cone.axis);
            if h <= 1e-9 {
                // Near or behind the apex.
                let radial = v.subtracted(&cone.axis.multiplied_scalar(h));
                let m = radial.magnitude();
                if m < 1e-12 {
                    return Some((cone.apex, cone.axis.multiplied_scalar(-1.0)));
                }
                return Some((cone.apex, radial.normalized().multiplied_scalar(-1.0)));
            }
            let radial = v.subtracted(&cone.axis.multiplied_scalar(h));
            let m = radial.magnitude();
            if m < 1e-12 {
                return Some((cone.apex, cone.axis.multiplied_scalar(-1.0)));
            }
            // Meridional closest point on the cone curve rho(z) = z·tan α.
            let tan_alpha = cone.semi_angle.tan();
            let rho = |z: f64| z * tan_alpha;
            let (rho_t, z_t) = closest_point_meridional(&rho, m, h, 0.0, h * 4.0 + 1.0)?;
            let q = cone
                .apex
                .translated_vec(&cone.axis.multiplied_scalar(z_t))
                .translated_vec(&radial.normalized().multiplied_scalar(rho_t));
            // Outward normal: cos α·ρ̂ + sin α·(−axis) (points away from the axis).
            let n = radial
                .normalized()
                .multiplied_scalar(cone.semi_angle.cos())
                .added(&cone.axis.multiplied_scalar(-cone.semi_angle.sin()));
            Some((q, n.normalized()))
        }
        _ => {
            // Fallback: grid closest point.
            let (u, v) = surface_closest_params(s, p, 24, 24);
            let q = s.d0(u, v);
            let n = surface_normal(s, u, v);
            if n.xyz().square_modulus() < 1e-30 {
                None
            } else {
                Some((q, n))
            }
        }
    }
}

/// Sample the rolling-ball centre *locus* for two surfaces as one or more 3D
/// curves. For the analytic surface-of-revolution pairs the locus is a circle
/// (returned as one closed polyline); for two parallel cylinders it is one or
/// two straight lines (each returned as a two-point polyline).
pub fn rolling_ball_center_surface(
    f1: &dyn Surface,
    f2: &dyn Surface,
    radius: f64,
    tol: f64,
) -> Result<Vec<Vec<GpPnt>>, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("rolling_ball_center_surface: radius must be positive".to_string());
    }
    let k1 = classify_surface_analytic(f1);
    let k2 = classify_surface_analytic(f2);
    let circle = |g: &TorusBlendGeom| sample_circle(&CircleGeom { center: g.center, normal: g.axis, radius: g.major }, CENTERLINE_SAMPLES);
    match (k1, k2) {
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let s1 = sphere_from_surface(f1).ok_or("cannot extract sphere 1")?;
            let s2 = sphere_from_surface(f2).ok_or("cannot extract sphere 2")?;
            let c = sphere_sphere_centerline(&s1, &s2, radius)?;
            Ok(vec![sample_circle(&c, CENTERLINE_SAMPLES)])
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let (pln, sph) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, sphere_from_surface(f2).ok_or("cannot extract sphere")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, sphere_from_surface(f1).ok_or("cannot extract sphere")?)
            };
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(vec![sample_circle(&c, CENTERLINE_SAMPLES)])
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let (pln, cyl) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, cylinder_from_surface(f2).ok_or("cannot extract cylinder")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, cylinder_from_surface(f1).ok_or("cannot extract cylinder")?)
            };
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Plane, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Plane) => {
            let (pln, cone) = if k1 == SurfaceKind::Plane {
                (plane_from_surface(f1).ok_or("cannot extract plane")?, cone_from_surface(f2).ok_or("cannot extract cone")?)
            } else {
                (plane_from_surface(f2).ok_or("cannot extract plane")?, cone_from_surface(f1).ok_or("cannot extract cone")?)
            };
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cone) | (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
            let (cyl, cone) = if k1 == SurfaceKind::Cylinder {
                (cylinder_from_surface(f1).ok_or("cannot extract cylinder")?, cone_from_surface(f2).ok_or("cannot extract cone")?)
            } else {
                (cylinder_from_surface(f2).ok_or("cannot extract cylinder")?, cone_from_surface(f1).ok_or("cannot extract cone")?)
            };
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cone, SurfaceKind::Sphere) | (SurfaceKind::Sphere, SurfaceKind::Cone) => {
            let (cone, sph) = if k1 == SurfaceKind::Cone {
                (cone_from_surface(f1).ok_or("cannot extract cone")?, sphere_from_surface(f2).ok_or("cannot extract sphere")?)
            } else {
                (cone_from_surface(f2).ok_or("cannot extract cone")?, sphere_from_surface(f1).ok_or("cannot extract sphere")?)
            };
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cone, SurfaceKind::Cone) => {
            let c1 = cone_from_surface(f1).ok_or("cannot extract cone 1")?;
            let c2 = cone_from_surface(f2).ok_or("cannot extract cone 2")?;
            let g = cone_cone_blend(&c1, &c2, radius, tol)?;
            Ok(vec![circle(&g)])
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => {
            let c1 = cylinder_from_surface(f1).ok_or("cannot extract cylinder 1")?;
            let c2 = cylinder_from_surface(f2).ok_or("cannot extract cylinder 2")?;
            let lines = parallel_cylinder_centerlines(&c1, &c2, radius, tol)?;
            Ok(lines)
        }
        _ => Err(format!(
            "rolling_ball_center_surface: unsupported pair ({k1:?}, {k2:?})"
        )),
    }
}

/// The straight rolling-ball centrelines for two parallel cylinders: the lines
/// at distance `r` from both cylinders (the intersections of the two offset
/// cylinders in the cross-section plane). Returns each line as two points.
pub(super) fn parallel_cylinder_centerlines(
    c1: &CylinderInfo,
    c2: &CylinderInfo,
    r: f64,
    tol: f64,
) -> Result<Vec<Vec<GpPnt>>, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("parallel_cylinder_centerlines: cylinders are not parallel".to_string());
    }
    let ax = a1;
    let (xh, yh) = project_basis(&ax);
    let ref_pt = c1.axis_origin;
    let o1 = Xy::from_proj(&c1.axis_origin, &xh, &yh, &ref_pt);
    let o2 = Xy::from_proj(&c2.axis_origin, &xh, &yh, &ref_pt);
    let hits = circle_circle_intersections(&o1, c1.radius + r, &o2, c2.radius + r);
    if hits.is_empty() {
        return Err(format!(
            "parallel_cylinder_centerlines: a ball of radius {r} cannot roll between these cylinders"
        ));
    }
    Ok(hits
        .iter()
        .map(|c| vec![c.to_pnt(-1.0, &xh, &yh, &ax, &ref_pt), c.to_pnt(1.0, &xh, &yh, &ax, &ref_pt)])
        .collect())
}

/// Verify that a blend surface is tangent to both adjacent faces at their
/// contact regions: sample the blend, and wherever a sample lies on either
/// face (within tolerance), require the blend normal to be parallel to the
/// face normal.
pub fn blend_tangency_ok(f1: &dyn Surface, f2: &dyn Surface, blend: &dyn Surface, tol: f64) -> bool {
    let (u0, u1) = blend.u_range();
    let (v0, v1) = blend.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (0.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    let eps = tol.max(1e-3);
    let (nu, nv) = (20, 20);
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = blend.d0(u, v);
            let nb = surface_normal(blend, u, v);
            if nb.xyz().square_modulus() < 1e-30 {
                continue;
            }
            for face in [f1, f2] {
                if let Some((q, nf)) = closest_point_on_surface(face, &p) {
                    if p.distance(&q) < eps {
                        if nf.xyz().crossed(&nb.xyz()).modulus() > 0.1 {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

/// Fillet a multi-normal corner: several curved edges meeting at `vertex`.
///
/// Each edge is filleted with `fillet_edge_curved_general`, then a spherical
/// corner patch centred on the rolling-ball locus closes the vertex. The patch
/// is a sphere of radius `radius` whose boundary is wired to the tangency
/// circles of the three incident blends (best-effort — the patch is added and
/// the shell re-closed; see the test for the supported configuration).
pub fn fillet_curved_multi_normal(
    solid: &TopoShape,
    vertex: &crate::shape::Vertex,
    edges: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_curved_multi_normal: radius must be a positive finite value".to_string());
    }
    if edges.len() < 3 {
        return Err(format!(
            "fillet_curved_multi_normal: needs at least 3 edges, got {}",
            edges.len()
        ));
    }
    // When the incident edges form the supported tri-sphere corner (three
    // sphere+sphere arcs), delegate to the deterministic sector + spherical
    // patch construction, which returns a closed shell.
    if let Ok(out) = fillet_curved_corner_patch(solid, vertex, edges, radius, tol) {
        return Ok(out);
    }
    // Best-effort fallback: fillet every edge sequentially, then add a
    // spherical corner patch centred on the rolling-ball locus.
    let original_edges = edges_of(solid);
    let mut current = solid.clone();
    for &i in edges {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_curved_multi_normal: edge index {i} out of range"))?;
        let (p0, p1) = BRepTool::edge_vertices(oe).ok_or("fillet_curved_multi_normal: edge has no curve")?;
        let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
            format!(
                "fillet_curved_multi_normal: edge {i} was consumed by an earlier fillet; \
                 the corner patch cannot re-find it"
            )
        })?;
        current = fillet_edge_curved_general(&current, &e, radius, tol)?;
    }
    // Add the spherical corner patch.
    let pv = BRepTool::vertex_point(vertex);
    let center = multi_normal_corner_center(&current, &pv, radius, tol)?;
    let patch = build_multi_normal_patch(&current, &pv, &center, radius)?;
    let mut faces = faces_of(&current);
    faces.push(patch);
    let b = TopoBuilder::new();
    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// The rolling-ball corner centre: the point at distance `radius` from each
/// face meeting at `pv`, found by least-squares on the signed distance error.
pub(super) fn multi_normal_corner_center(
    solid: &TopoShape,
    pv: &GpPnt,
    radius: f64,
    _tol: f64,
) -> Result<GpPnt, String> {
    // Collect the faces incident at the vertex.
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
            "multi_normal_corner_center: vertex touches {} faces (expected >= 3)",
            faces.len()
        ));
    }
    // Signed distance function per face: outward positive.
    let mut c = pv.translated_vec(&occt_core::gp::GpVec::new(radius, radius, radius));
    for _ in 0..200 {
        let mut gx = 0.0f64;
        let mut gy = 0.0f64;
        let mut gz = 0.0f64;
        let mut grad = [[0.0; 3]; 3];
        let mut err_sum = 0.0;
        for f in &faces {
            let s = BRepTool::face_surface(f).ok_or("corner face has no surface")?;
            let (q, n) = closest_point_on_surface(s.as_ref(), &c).ok_or("cannot project onto a corner face")?;
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
        if err_sum / (faces.len() as f64) < 1e-6 {
            return Ok(c);
        }
        // Newton step: grad · Δ = g.
        if let Some(delta) = solve3x3(&grad, &[gx, gy, gz]) {
            let step = GpVec::new(delta[0], delta[1], delta[2]);
            if step.magnitude() < 1e-12 {
                return Ok(c);
            }
            c = c.translated_vec(&step);
        } else {
            break;
        }
    }
    Ok(c)
}

/// Build the spherical corner patch face for a multi-normal corner: a sphere of
/// radius `radius` centred at `center`, bounded by a wire of three circular
/// arcs that pass through points on the three incident edges at distance
/// `2·radius` from the vertex.
pub(super) fn build_multi_normal_patch(
    solid: &TopoShape,
    pv: &GpPnt,
    center: &GpPnt,
    radius: f64,
) -> Result<Face, String> {
    // Direction of the three incident edges (from the vertex into the solid).
    let mut dirs: Vec<GpVec> = Vec::new();
    for e in edges_of(solid) {
        let (a, b) = BRepTool::edge_vertices(&e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        let other = if a.distance(pv) < 1e-6 { b } else if b.distance(pv) < 1e-6 { a } else { continue };
        let v = GpVec::from_pnts(pv, &other);
        if v.magnitude() > 1e-12 {
            dirs.push(v.normalized());
        }
    }
    if dirs.len() < 3 {
        return Err(format!(
            "build_multi_normal_patch: vertex has {} incident edges (expected >= 3)",
            dirs.len()
        ));
    }
    let dirs = [dirs[0], dirs[1], dirs[2]];
    // Points on the three edges at distance 2R from the vertex.
    let p = [
        pv.translated_vec(&dirs[0].multiplied_scalar(2.0 * radius)),
        pv.translated_vec(&dirs[1].multiplied_scalar(2.0 * radius)),
        pv.translated_vec(&dirs[2].multiplied_scalar(2.0 * radius)),
    ];
    // Build the three arcs on the sphere surface.
    let b = TopoBuilder::new();
    let mut edges: Vec<Edge> = Vec::new();
    for k in 0..3 {
        let (a, bb) = (p[k], p[(k + 1) % 3]);
        let mid_dir = dirs[k]
            .crossed(&dirs[(k + 1) % 3])
            .normalized();
        let nd = GpDir::from_xyz(&mid_dir.xyz()).map_err(|e| e.to_string())?;
        let edge = build_sphere_arc(&b, center, radius, &nd, &a, &bb)?;
        edges.push(edge);
    }
    let wire = b.make_wire(&edges);
    let mut sph = GpSphere::new(GpAx3::standard(), radius).map_err(|e| e.to_string())?;
    sph.set_location(*center);
    Ok(b.make_face(Arc::new(GeomSphere::new(sph)), &[wire]))
}

/// A circular arc on a sphere of radius `radius` centred at `center`, in the
/// plane through `center` with normal `nd`, from `a` to `b` (both on the
/// sphere).
pub(super) fn build_sphere_arc(
    b: &TopoBuilder,
    center: &GpPnt,
    radius: f64,
    nd: &GpDir,
    a: &GpPnt,
    bb: &GpPnt,
) -> Result<Edge, String> {
    let xd = GpDir::from_vec(&GpVec::from_pnts(center, a))
        .map_err(|_| "fillet_curved: sphere arc start coincides with its centre".to_string())?;
    let ax2 = GpAx2::new(*center, *nd, xd).map_err(|e| format!("fillet_curved: sphere arc frame: {e}"))?;
    let va = GpVec::from_pnts(center, a);
    let vb = GpVec::from_pnts(center, bb);
    let xdir = *ax2.x_direction();
    let ydir = *ax2.y_direction();
    let a1 = va.coord.dot(&ydir.xyz()).atan2(va.coord.dot(&xdir.xyz()));
    let a2 = vb.coord.dot(&ydir.xyz()).atan2(vb.coord.dot(&xdir.xyz()));
    let mut e = b.make_edge_circle(&ax2, radius, a1, a2);
    b.add(&mut e.0, &b.make_vertex(*a, 0.0).0);
    b.add(&mut e.0, &b.make_vertex(*bb, 0.0).0);
    Ok(e)
}

/// Fillet several edges of `solid` sequentially with `fillet_edge_curved`.
///
/// Edges are looked up by their geometric endpoints in the current solid, so a
/// chain whose edges do not interfere (no adjacent face is rebuilt twice)
/// behaves like repeated single-edge fillets. When an earlier fillet consumed
/// an edge, the call errors rather than fillet a wrong edge.
pub fn fillet_edge_curved_chain(
    solid: &TopoShape,
    edge_indices: &[usize],
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    let original_edges = edges_of(solid);
    // When >= 3 of the chain's edges meet at a common vertex, close that
    // corner with the deterministic sector + spherical patch construction
    // (the multi-edge chain start/end patch). Otherwise fall back to repeated
    // single-edge fillets.
    let e_objs: Vec<Edge> = edge_indices
        .iter()
        .map(|&i| original_edges.get(i).cloned().ok_or_else(|| format!("fillet_edge_curved_chain: edge index {i} out of range")))
        .collect::<Result<_, _>>()?;
    let mut candidate: Option<crate::shape::Vertex> = None;
    if e_objs.len() >= 3 {
        let mut pts: Vec<GpPnt> = Vec::new();
        for e in &e_objs {
            if let Some((a, b)) = BRepTool::edge_vertices(e) {
                pts.push(a);
                pts.push(b);
            }
        }
        for p in &pts {
            let count = pts.iter().filter(|q| q.distance(p) < 1e-6).count();
            if count >= 3 {
                if let Some(v) = vertices_of(solid).into_iter().find(|v| BRepTool::vertex_point(v).distance(p) < 1e-6) {
                    candidate = Some(v);
                    break;
                }
            }
        }
    }
    if let Some(v) = candidate {
        return fillet_curved_corner_patch(solid, &v, edge_indices, radius, tol);
    }
    let mut current = solid.clone();
    for &i in edge_indices {
        let oe = original_edges
            .get(i)
            .ok_or_else(|| format!("fillet_edge_curved_chain: edge index {i} out of range"))?;
        let (p0, p1) =
            BRepTool::edge_vertices(oe).ok_or("fillet_edge_curved_chain: edge has no curve")?;
        let e = find_edge_by_endpoints(&current, &p0, &p1).ok_or_else(|| {
            format!(
                "fillet_edge_curved_chain: edge {i} was consumed by an earlier fillet and can no \
                 longer be re-found"
            )
        })?;
        current = fillet_edge_curved(&current, &e, radius, tol)?;
    }
    Ok(current)
}

/// Find an edge of `shape` whose endpoints coincide with `p0`/`p1` (either
/// order). For a closed circle both endpoints are the seam point.
pub(super) fn find_edge_by_endpoints(shape: &TopoShape, p0: &GpPnt, p1: &GpPnt) -> Option<Edge> {
    edges_of(shape).into_iter().find(|e| {
        let (a, b) = BRepTool::edge_vertices(e).unwrap_or((GpPnt::zero(), GpPnt::zero()));
        (a.distance(p0) < 1e-6 && b.distance(p1) < 1e-6)
            || (a.distance(p1) < 1e-6 && b.distance(p0) < 1e-6)
    })
}

// ===========================================================================
// Deterministic multi-edge corner patch (Phase 12)
// ===========================================================================
//
// The tri-sphere corner: three spheres, pairwise intersecting in circular arcs
// that meet at two common vertices, bound a closed solid of three spherical
// caps. Filleting one sphere+sphere arc edge independently with a full annulus
// would overlap the adjacent fillets near the vertices, so each edge is
// filleted as a torus SECTOR spanning the azimuth interval between the two
// rolling-ball corner centres. The trimmed caps become spherical bands bounded
// by two tangency arcs, and each vertex is closed by a spherical patch (the
// rolling ball tangent to all three spheres) bounded by the three sector
// cross-section arcs. Every edge is then shared by exactly two faces, so the
// rebuilt shell is closed deterministically.

/// Azimuth (angle in the plane perpendicular to `axis`) of `p` around `axis`,
/// measured from the basis `(e1, e2)` with reference point `ref_pt` on the
/// axis. All circles of a sphere+sphere blend (the edge arc, the centreline,
/// the two tangency circles) are coaxial, so the same azimuth identifies the
/// matching points across them.
pub(super) fn point_azimuth(p: &GpPnt, ref_pt: &GpPnt, axis: &GpVec, e1: &GpVec, e2: &GpVec) -> f64 {
    let v = GpVec::from_pnts(ref_pt, p);
    let along = v.dot(axis);
    let perp = v.subtracted(&axis.multiplied_scalar(along));
    perp.dot(e2).atan2(perp.dot(e1))
}

/// The faces of `solid` that reference `edge` in one of their wires.
pub(super) fn adjacent_faces(solid: &TopoShape, edge: &Edge) -> Vec<Face> {
    faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f)
                .iter()
                .any(|w| edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0)))
        })
        .collect()
}

/// Rebuild a spherical adjacent face: drop the wire holding `fillet_edge`, keep
/// every other wire, and add `new_wire`. This is what lets a sphere that is
/// adjacent to several fillet edges keep its other tangency arcs (the
/// tri-sphere corner needs each cap to become a band bounded by two arcs).
pub(super) fn rebuild_sphere_keep(face: &Face, fillet_edge: &Edge, new_wire: &Wire, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: sphere face has no surface")?;
    let mut wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(face) {
        let edges = edges_of_wire(&w);
        let orig_len = edges.len();
        let kept: Vec<Edge> = edges
            .into_iter()
            .filter(|e| !is_same(&e.0, &fillet_edge.0))
            .collect();
        if kept.is_empty() {
            continue;
        }
        wires.push(if kept.len() < orig_len { b.make_wire(&kept) } else { w });
    }
    wires.push(new_wire.clone());
    Ok(b.make_face(surf, &wires))
}

/// Build a circular-arc edge on a tangency circle `ct` (normal `dhat`, basis
/// `e1` at azimuth 0) over the azimuth range `[start, end]`.
pub(super) fn build_contact_arc(
    b: &TopoBuilder,
    ct: &ContactCircleGeom,
    dhat: &GpVec,
    e1: &GpVec,
    start: f64,
    end: f64,
    p_start: &GpPnt,
    p_end: &GpPnt,
) -> Result<Edge, String> {
    let nd = GpDir::from_xyz(&dhat.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&e1.xyz()).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(ct.center, nd, xd).map_err(|e| format!("fillet_curved: tangency arc frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, ct.radius, start, end);
    b.add(&mut e.0, &b.make_vertex(*p_start, 0.0).0);
    b.add(&mut e.0, &b.make_vertex(*p_end, 0.0).0);
    Ok(e)
}

/// The tube (minor) angle of a tangency circle in a torus band: the circle is
/// at `g.center + dhat·r·sin(v)`, radius `g.major + r·cos(v)`.
pub(super) fn tube_angle(g: &TorusBlendGeom, ct: &ContactCircleGeom, radius: f64) -> f64 {
    let sin_v = GpVec::from_pnts(&g.center, &ct.center).dot(&g.axis) / radius;
    let cos_v = (ct.radius - g.major) / radius;
    sin_v.atan2(cos_v)
}

/// The cross-section arc of a torus sector at azimuth `theta`: the minor circle
/// (radius `radius`) through the two tangency points at that azimuth, in the
/// plane perpendicular to the centreline tangent. This is the arc shared with
/// the spherical corner patch.
pub(super) fn build_sector_cross_arc(
    b: &TopoBuilder,
    g: &TorusBlendGeom,
    radius: f64,
    theta: f64,
    e1: &GpVec,
    e2: &GpVec,
    dhat: &GpVec,
) -> Result<Edge, String> {
    let rho = e1
        .multiplied_scalar(theta.cos())
        .added(&e2.multiplied_scalar(theta.sin()));
    let center = g.center.translated_vec(&rho.multiplied_scalar(g.major));
    let v_a = tube_angle(g, &g.contact_a, radius);
    let v_b = tube_angle(g, &g.contact_b, radius);
    // Minor circle plane normal: −tangent at `theta` (the +v direction is +dhat).
    let t = dhat.crossed(&rho).normalized();
    let nd = GpDir::from_vec(&t.multiplied_scalar(-1.0)).map_err(|e| e.to_string())?;
    let xd = GpDir::from_vec(&rho).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(center, nd, xd).map_err(|e| format!("fillet_curved: cross arc frame: {e}"))?;
    let mut arc = b.make_edge_circle(&ax2, radius, v_a, v_b);
    let p_a = center
        .translated_vec(&rho.multiplied_scalar(radius * v_a.cos()).added(&dhat.multiplied_scalar(radius * v_a.sin())));
    let p_b = center
        .translated_vec(&rho.multiplied_scalar(radius * v_b.cos()).added(&dhat.multiplied_scalar(radius * v_b.sin())));
    b.add(&mut arc.0, &b.make_vertex(p_a, 0.0).0);
    b.add(&mut arc.0, &b.make_vertex(p_b, 0.0).0);
    Ok(arc)
}
