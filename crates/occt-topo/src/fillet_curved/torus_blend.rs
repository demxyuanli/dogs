use super::prelude::*;
use super::*;

/// The point on a cone-tangency circle that shares a generator with the cone's
/// base-circle seam, so the rebuilt lateral face's seam stays a straight
/// generator from the contact circle to the apex.
pub(super) fn cone_contact_seam(_face: &Face, fillet_edge: &Edge, contact: &ContactCircleGeom) -> Result<GpPnt, String> {
    let (se, _) = BRepTool::edge_vertices(fillet_edge)
        .ok_or("fillet_curved: cone base circle has no vertices")?;
    let base_seam = se;
    let v = GpVec::from_pnts(&contact.center, &base_seam);
    let along = v.dot(&contact.normal);
    let radial = v.subtracted(&contact.normal.multiplied_scalar(along));
    let m = radial.magnitude();
    if m < 1e-9 {
        return Err("fillet_curved: cone base seam lies on the axis".to_string());
    }
    Ok(contact.center.translated_vec(&radial.multiplied_scalar(contact.radius / m)))
}

/// Rebuild a cone lateral face: the base (filleted) circle is replaced by the
/// tangency circle and the seam is re-generated from the contact circle up to
/// the apex. A *truncated* cone face (bounded by a second closed rim circle, as
/// when the cone carries a top cap) keeps that rim and is rebuilt as a band
/// from the tangency circle up to it, mirroring the cylinder rebuild.
pub(super) fn rebuild_cone_face(
    face: &Face,
    fillet_edge: &Edge,
    contact_edge: &Edge,
    contact_seam: &GpPnt,
    b: &TopoBuilder,
) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: cone face has no surface")?;
    if let Ok(top) = find_top_circle(face, fillet_edge) {
        let (s_top, _) = BRepTool::edge_vertices(&top)
            .ok_or("fillet_curved: cone rim has no vertices")?;
        let seam_up = b.make_edge_segment(contact_seam, &s_top);
        let wire = b.make_wire(&[contact_edge.clone(), seam_up.clone(), top, seam_up]);
        return Ok(b.make_face(surf, &[wire]));
    }
    let apex = find_cone_apex(face, fillet_edge)?;
    let new_seam = b.make_edge_segment(contact_seam, &apex);
    let wire = b.make_wire(&[contact_edge.clone(), new_seam.clone(), new_seam]);
    Ok(b.make_face(surf, &[wire]))
}

/// The `GeomTorus` surface of a torus band blend.
pub(super) fn torus_surface(geom: &TorusBlendGeom) -> Result<Arc<dyn Surface>, String> {
    let zd = GpDir::from_xyz(&geom.axis.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&perp_dir(&geom.axis).xyz()).map_err(|e| e.to_string())?;
    let ax3 = GpAx3::new(geom.center, zd, &xd).map_err(|e| format!("fillet_curved: torus frame: {e}"))?;
    let t = GpTorus::new(ax3, geom.major, geom.minor).map_err(|e| e.to_string())?;
    Ok(Arc::new(GeomTorus::new(t)))
}

/// Build the blend face: the torus band bounded by the two tangency circles
/// (each shared with one trimmed adjacent face).
pub(super) fn build_torus_blend_face(geom: &TorusBlendGeom, edge_a: &Edge, edge_b: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = torus_surface(geom)?;
    let wa = b.make_wire(&[edge_a.clone()]);
    let wb = b.make_wire(&[edge_b.clone()]);
    // The larger tangency circle is the outer boundary.
    let (w_outer, w_inner) = if geom.contact_a.radius >= geom.contact_b.radius {
        (wa, wb)
    } else {
        (wb, wa)
    };
    Ok(b.make_face(surf, &[w_outer, w_inner]))
}

/// Build a trimmed adjacent face, returning the rebuilt face and the (shared)
/// tangency circle edge it uses.
#[allow(clippy::type_complexity)]
pub(super) fn build_adjacent_face(
    face: &Face,
    kind: AdjacentKind,
    fillet_edge: &Edge,
    contact: &ContactCircleGeom,
    b: &TopoBuilder,
) -> Result<(Face, Edge), String> {
    match kind {
        AdjacentKind::Plane => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_plane_face(face, fillet_edge, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Sphere => {
            let e = build_circle_edge(b, contact, &default_seam(contact))?;
            let f = rebuild_single_circle_face(face, &e, b)?;
            Ok((f, e))
        }
        AdjacentKind::Cylinder => {
            let seam = cylinder_contact_seam(face, fillet_edge, contact)?;
            let e = build_circle_edge(b, contact, &seam)?;
            let f = rebuild_cylinder_face(face, fillet_edge, &e, &seam, b)?;
            Ok((f, e))
        }
        AdjacentKind::Cone => {
            let seam = cone_contact_seam(face, fillet_edge, contact)?;
            let e = build_circle_edge(b, contact, &seam)?;
            let f = rebuild_cone_face(face, fillet_edge, &e, &seam, b)?;
            Ok((f, e))
        }
    }
}

// ===========================================================================
// Public API
// ===========================================================================

/// Sample the rolling-ball centre locus for two surfaces. For the supported
/// analytic pairs the locus is a circle (plane+sphere, sphere+sphere,
/// plane+cylinder); the returned points lie on it.
///
/// Each returned centre point is at distance `radius` from the first surface
/// and `radius` from the second surface (measuring against the *offset*
/// surfaces — the plane pushed out by `r`, the sphere/cylinder enlarged by
/// `R + r`). Unsupported pairs return an error.
pub fn rolling_ball_center_line(
    f1: &dyn Surface,
    f2: &dyn Surface,
    radius: f64,
    tol: f64,
) -> Result<Vec<GpPnt>, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("rolling_ball_center_line: radius must be positive".to_string());
    }
    let k1 = classify_surface_full(f1);
    let k2 = classify_surface_full(f2);
    match (k1, k2) {
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let s1 = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere 1")?;
            let s2 = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere 2")?;
            let c = sphere_sphere_centerline(&s1, &s2, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pln = plane_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let sph = sphere_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let sph = sphere_from_surface(f1).ok_or("rolling_ball_center_line: cannot extract sphere")?;
            let pln = plane_from_surface(f2).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let c = plane_sphere_centerline(&pln, &sph, radius)?;
            Ok(sample_circle(&c, CENTERLINE_SAMPLES))
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) | (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let (pln_s, cyl_s) = if k1 == SurfaceKind::Plane { (f1, f2) } else { (f2, f1) };
            let pln = plane_from_surface(pln_s).ok_or("rolling_ball_center_line: cannot extract plane")?;
            let cyl = cylinder_from_surface(cyl_s).ok_or("rolling_ball_center_line: cannot extract cylinder")?;
            let blend = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            Ok(sample_circle(
                &CircleGeom { center: blend.center, normal: blend.axis, radius: blend.major },
                CENTERLINE_SAMPLES,
            ))
        }
        _ => Err(format!(
            "rolling_ball_center_line: unsupported surface pair ({k1:?}, {k2:?})"
        )),
    }
}

/// Replace a sharp edge of `solid` with a curved-face rolling-ball fillet.
///
/// This is the original supported-pairs entry point (plane+sphere,
/// plane+cylinder, sphere+sphere, and — through `fillet_edge_curved_general` —
/// the cone/cylinder pairs). The blend surface is a torus band whose boundary —
/// the two tangency circles — is shared with the trimmed adjacent faces, so the
/// rebuilt shell stays closed.
pub fn fillet_edge_curved(solid: &TopoShape, edge: &Edge, radius: f64, tol: f64) -> Result<TopoShape, String> {
    fillet_edge_curved_general(solid, edge, radius, tol)
}

/// Assemble a solid from a torus-band blend and its two trimmed adjacent faces.
#[allow(clippy::too_many_arguments)]
pub(super) fn assemble_torus_blend(
    solid: &TopoShape,
    edge: &Edge,
    f1: &Face,
    f2: &Face,
    blend: &TorusBlendGeom,
    kind_a: AdjacentKind,
    kind_b: AdjacentKind,
    contact_a: ContactCircleGeom,
    contact_b: ContactCircleGeom,
) -> Result<TopoShape, String> {
    let b = TopoBuilder::new();
    let (trimmed1, e1) = build_adjacent_face(f1, kind_a, edge, &contact_a, &b)?;
    let (trimmed2, e2) = build_adjacent_face(f2, kind_b, edge, &contact_b, &b)?;
    let blend_face = build_torus_blend_face(blend, &e1, &e2, &b)?;

    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0) || is_same(&f.0, &f2.0) {
            continue;
        }
        faces.push(f);
    }
    faces.push(trimmed1);
    faces.push(trimmed2);
    faces.push(blend_face);

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// The general curved-face rolling-ball fillet. Supports the plane/sphere/
/// cylinder/cone pairs whose offset surfaces intersect in a circle (torus
/// blend), two parallel cylinders (cylindrical band blend), and returns a
/// documented error for torus/torus or other non-analytic pairs.
pub fn fillet_edge_curved_general(
    solid: &TopoShape,
    edge: &Edge,
    radius: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    if radius <= 0.0 || !radius.is_finite() {
        return Err("fillet_edge_curved_general: radius must be a positive finite value".to_string());
    }
    let adjacent = faces_of(solid)
        .into_iter()
        .filter(|f| {
            wires_of_face(f).iter().any(|w| {
                edges_of_wire(w).iter().any(|e| is_same(&e.0, &edge.0))
            })
        })
        .collect::<Vec<Face>>();
    if adjacent.len() != 2 {
        return Err(format!(
            "fillet_edge_curved_general: edge is adjacent to {} faces (expected 2)",
            adjacent.len()
        ));
    }
    let f1 = &adjacent[0];
    let f2 = &adjacent[1];
    let s1 = BRepTool::face_surface(f1).ok_or("fillet_edge_curved_general: face 1 has no surface")?;
    let s2 = BRepTool::face_surface(f2).ok_or("fillet_edge_curved_general: face 2 has no surface")?;
    let k1 = classify_surface_analytic(s1.as_ref());
    let k2 = classify_surface_analytic(s2.as_ref());

    match (k1, k2) {
        (SurfaceKind::Plane, SurfaceKind::Sphere) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let sph = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere")?;
            let g = plane_sphere_blend(&pln, &sph, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Sphere, SurfaceKind::Plane) => {
            let sph = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_sphere_blend(&pln, &sph, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Plane, SurfaceKind::Cylinder) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let cyl = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder")?;
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Cylinder, ca, cb)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Plane) => {
            let cyl = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_cylinder_blend(&pln, &cyl, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cylinder, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Sphere, SurfaceKind::Sphere) => {
            let sp1 = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere 1")?;
            let sp2 = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere 2")?;
            let g = sphere_sphere_blend(&sp1, &sp2, radius)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Plane, SurfaceKind::Cone) => {
            let pln = plane_from_surface(s1.as_ref()).ok_or("cannot extract plane")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Plane, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cone, SurfaceKind::Plane) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let pln = plane_from_surface(s2.as_ref()).ok_or("cannot extract plane")?;
            let g = plane_cone_blend(&pln, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Plane, cb, ca)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cone) => {
            let cyl = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cylinder, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cone, SurfaceKind::Cylinder) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let cyl = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder")?;
            let g = cylinder_cone_blend(&cyl, &cone, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Cylinder, cb, ca)
        }
        (SurfaceKind::Cone, SurfaceKind::Sphere) => {
            let cone = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone")?;
            let sph = sphere_from_surface(s2.as_ref()).ok_or("cannot extract sphere")?;
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Sphere, ca, cb)
        }
        (SurfaceKind::Sphere, SurfaceKind::Cone) => {
            let sph = sphere_from_surface(s1.as_ref()).ok_or("cannot extract sphere")?;
            let cone = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone")?;
            let g = cone_sphere_blend(&cone, &sph, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Sphere, AdjacentKind::Cone, cb, ca)
        }
        (SurfaceKind::Cone, SurfaceKind::Cone) => {
            let c1 = cone_from_surface(s1.as_ref()).ok_or("cannot extract cone 1")?;
            let c2 = cone_from_surface(s2.as_ref()).ok_or("cannot extract cone 2")?;
            let g = cone_cone_blend(&c1, &c2, radius, tol)?;
            let (ca, cb) = (g.contact_a, g.contact_b);
            assemble_torus_blend(solid, edge, f1, f2, &g, AdjacentKind::Cone, AdjacentKind::Cone, ca, cb)
        }
        (SurfaceKind::Cylinder, SurfaceKind::Cylinder) => {
            let c1 = cylinder_from_surface(s1.as_ref()).ok_or("cannot extract cylinder 1")?;
            let c2 = cylinder_from_surface(s2.as_ref()).ok_or("cannot extract cylinder 2")?;
            fillet_edge_parallel_cylinders(solid, edge, f1, f2, &c1, &c2, radius, tol)
        }
        _ => Err(format!(
            "fillet_edge_curved_general: unsupported face pair ({k1:?}, {k2:?}); supported pairs \
             are plane+sphere, plane+cylinder, plane+cone, sphere+sphere, coaxial cylinder+cone, \
             sphere+cone, cone+cone and parallel cylinder+cylinder"
        )),
    }
}

/// A 2D cross-section point (in the plane perpendicular to the cylinder axes).
#[derive(Debug, Clone, Copy)]
pub(super) struct Xy {
    pub(super) x: f64,
    pub(super) y: f64,
}

impl Xy {
    pub(super) fn from_proj(p: &GpPnt, xh: &GpVec, yh: &GpVec, ref_pt: &GpPnt) -> Xy {
        let v = GpVec::from_pnts(ref_pt, p);
        Xy { x: v.dot(xh), y: v.dot(yh) }
    }
    pub(super) fn to_pnt(&self, z: f64, xh: &GpVec, yh: &GpVec, ax: &GpVec, ref_pt: &GpPnt) -> GpPnt {
        ref_pt.translated_vec(&xh.multiplied_scalar(self.x).added(&yh.multiplied_scalar(self.y)).added(&ax.multiplied_scalar(z)))
    }
    pub(super) fn sub(&self, o: &Xy) -> Xy { Xy { x: self.x - o.x, y: self.y - o.y } }
    pub(super) fn add(&self, o: &Xy) -> Xy { Xy { x: self.x + o.x, y: self.y + o.y } }
    pub(super) fn scale(&self, s: f64) -> Xy { Xy { x: self.x * s, y: self.y * s } }
    pub(super) fn dot(&self, o: &Xy) -> f64 { self.x * o.x + self.y * o.y }
    pub(super) fn cross(&self, o: &Xy) -> f64 { self.x * o.y - self.y * o.x }
    pub(super) fn norm(&self) -> f64 { self.dot(self).sqrt() }
    #[allow(dead_code)]
    pub(super) fn normalized(&self) -> Xy {
        let n = self.norm();
        if n > 1e-30 { self.scale(1.0 / n) } else { Xy { x: 1.0, y: 0.0 } }
    }
}

/// Intersections of two circles (centres `c1`/`c2`, radii `r1`/`r2`) in the
/// cross-section plane. Returns 0, 1 or 2 points.
pub(super) fn circle_circle_intersections(c1: &Xy, r1: f64, c2: &Xy, r2: f64) -> Vec<Xy> {
    let d = c2.sub(c1).norm();
    if d < 1e-12 {
        return Vec::new();
    }
    if d > r1 + r2 + 1e-12 || d < (r1 - r2).abs() - 1e-12 {
        return Vec::new();
    }
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let h2 = r1 * r1 - a * a;
    let h = if h2 < 0.0 { 0.0 } else { h2.sqrt() };
    let dir = c2.sub(c1).scale(1.0 / d);
    let n = Xy { x: -dir.y, y: dir.x };
    let base = c1.add(&dir.scale(a));
    if h < 1e-12 {
        vec![base]
    } else {
        vec![base.add(&n.scale(h)), base.sub(&n.scale(h))]
    }
}

/// Build a circular-arc edge in the plane perpendicular to `ax` at height `z`
/// (relative to `ref_pt`), from cross-section angles `a1` to `a2` (radians).
pub(super) fn build_xy_arc(
    b: &TopoBuilder,
    center_xy: &Xy,
    z: f64,
    r: f64,
    a1: f64,
    a2: f64,
    xh: &GpVec,
    yh: &GpVec,
    ax: &GpVec,
    ref_pt: &GpPnt,
) -> Result<Edge, String> {
    let center = center_xy.to_pnt(z, xh, yh, ax, ref_pt);
    let nd = GpDir::from_xyz(&ax.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_xyz(&xh.xyz()).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(center, nd, xd).map_err(|e| format!("fillet_curved: arc frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, r, a1, a2);
    let v1 = b.make_vertex(center_xy.add(&Xy { x: r * a1.cos(), y: r * a1.sin() }).to_pnt(z, xh, yh, ax, ref_pt), 0.0);
    let v2 = b.make_vertex(center_xy.add(&Xy { x: r * a2.cos(), y: r * a2.sin() }).to_pnt(z, xh, yh, ax, ref_pt), 0.0);
    b.add_edge_vertices(&mut e, &v1, &v2);
    Ok(e)
}

/// Rolling-ball fillet between two parallel (non-coaxial) cylinder lateral
/// faces sharing a generator line. The blend is a cylindrical band of radius
/// `r` whose axis is the rolling-ball centreline (the line at distance `r`
/// from both cylinders, on the material side of the shared edge).
pub(super) fn fillet_edge_parallel_cylinders(
    solid: &TopoShape,
    edge: &Edge,
    f1: &Face,
    f2: &Face,
    c1: &CylinderInfo,
    c2: &CylinderInfo,
    r: f64,
    tol: f64,
) -> Result<TopoShape, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_edge_parallel_cylinders: the two cylinders are not parallel".to_string());
    }
    let ax = a1;
    let (xh, yh) = project_basis(&ax);
    let ref_pt = c1.axis_origin;
    let o1 = Xy::from_proj(&c1.axis_origin, &xh, &yh, &ref_pt);
    let o2 = Xy::from_proj(&c2.axis_origin, &xh, &yh, &ref_pt);

    // Axial extent of the fillet edge.
    let (ep0, ep1) = BRepTool::edge_vertices(edge).ok_or("fillet_edge_parallel_cylinders: edge has no vertices")?;
    let z0 = GpVec::from_pnts(&ref_pt, &ep0).dot(&ax);
    let z1 = GpVec::from_pnts(&ref_pt, &ep1).dot(&ax);
    let (z_lo, z_hi) = (z0.min(z1), z0.max(z1));
    if (z_hi - z_lo).abs() < 1e-9 {
        return Err("fillet_edge_parallel_cylinders: degenerate edge".to_string());
    }
    let edge_mid = GpPnt::new(
        0.5 * (ep0.x() + ep1.x()),
        0.5 * (ep0.y() + ep1.y()),
        0.5 * (ep0.z() + ep1.z()),
    );
    let p_edge = Xy::from_proj(&edge_mid, &xh, &yh, &ref_pt);

    // Inner offset circles: the ball centre is at distance r inside both.
    let r1o = c1.radius - r;
    let r2o = c2.radius - r;
    if r1o <= 1e-9 || r2o <= 1e-9 {
        return Err(format!(
            "fillet_edge_parallel_cylinders: radius {r} does not fit between the cylinders"
        ));
    }
    let hits = circle_circle_intersections(&o1, r1o, &o2, r2o);
    if hits.is_empty() {
        return Err("fillet_edge_parallel_cylinders: no rolling-ball centreline".to_string());
    }
    let side_edge = p_edge.sub(&o1).cross(&o2.sub(&o1)).signum();
    let c_xy = hits
        .iter()
        .find(|p| p.sub(&o1).cross(&o2.sub(&o1)).signum() == side_edge)
        .copied()
        .unwrap_or(hits[0]);
    let t1 = o1.add(&c_xy.sub(&o1).scale(c1.radius / r1o));
    let t2 = o2.add(&c_xy.sub(&o2).scale(c2.radius / r2o));

    // The other intersection line of the two cylinder surfaces.
    let raw = circle_circle_intersections(&o1, c1.radius, &o2, c2.radius);
    if raw.len() < 2 {
        return Err("fillet_edge_parallel_cylinders: cylinders do not cross in a lens".to_string());
    }
    let p2 = raw
        .iter()
        .find(|p| p.sub(&o1).cross(&o2.sub(&o1)).signum() != side_edge)
        .copied()
        .unwrap_or(raw[0]);

    let b = TopoBuilder::new();

    let ang = |p: &Xy, o: &Xy| (p.y - o.y).atan2(p.x - o.x);
    let a_t1 = ang(&t1, &o1);
    let a_p2_1 = ang(&p2, &o1);
    let a_t2 = ang(&t2, &o2);
    let a_p2_2 = ang(&p2, &o2);

    let outer1 = o1.sub(&o2.sub(&o1).normalized().scale(c1.radius));
    let a_outer1 = ang(&outer1, &o1);
    let a_t1_hi = sweep_through(a_t1, a_p2_1, a_outer1);
    let outer2 = o2.sub(&o1.sub(&o2).normalized().scale(c2.radius));
    let a_outer2 = ang(&outer2, &o2);
    let a_t2_hi = sweep_through(a_t2, a_p2_2, a_outer2);

    let t1_bottom = t1.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let t1_top = t1.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);
    let t2_bottom = t2.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let t2_top = t2.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);
    let p2_bottom = p2.to_pnt(z_lo, &xh, &yh, &ax, &ref_pt);
    let p2_top = p2.to_pnt(z_hi, &xh, &yh, &ax, &ref_pt);

    let tangency1 = b.make_edge_segment(&t1_bottom, &t1_top);
    let tangency2 = b.make_edge_segment(&t2_bottom, &t2_top);
    let p2_line = b.make_edge_segment(&p2_bottom, &p2_top);

    // Blend arc in the two cap planes (the SHORT arc from T1 to T2).
    let blend_center = c_xy;
    let mut blend_sweep = a_t2 - a_t1;
    while blend_sweep > PI {
        blend_sweep -= 2.0 * PI;
    }
    while blend_sweep < -PI {
        blend_sweep += 2.0 * PI;
    }
    let a_blend_end = a_t1 + blend_sweep;
    let blend_bottom_arc = build_xy_arc(&b, &blend_center, z_lo, r, a_t1, a_blend_end, &xh, &yh, &ax, &ref_pt)?;
    let blend_top_arc = build_xy_arc(&b, &blend_center, z_hi, r, a_t1, a_blend_end, &xh, &yh, &ax, &ref_pt)?;

    // Retained arcs on each cylinder lateral face.
    let arc1_bottom = build_xy_arc(&b, &o1, z_lo, c1.radius, a_t1, a_t1_hi, &xh, &yh, &ax, &ref_pt)?;
    let arc1_top = build_xy_arc(&b, &o1, z_hi, c1.radius, a_t1_hi, a_t1 + 2.0 * PI, &xh, &yh, &ax, &ref_pt)?;
    let arc2_bottom = build_xy_arc(&b, &o2, z_lo, c2.radius, a_t2, a_t2_hi, &xh, &yh, &ax, &ref_pt)?;
    let arc2_top = build_xy_arc(&b, &o2, z_hi, c2.radius, a_t2_hi, a_t2 + 2.0 * PI, &xh, &yh, &ax, &ref_pt)?;

    // Lateral faces.
    let surf1 = BRepTool::face_surface(f1).ok_or("cylinder face 1 has no surface")?;
    let lat1_wire = b.make_wire(&[arc1_bottom.clone(), p2_line.clone(), arc1_top.clone(), tangency1.clone()]);
    let lat1 = b.make_face(surf1, &[lat1_wire]);
    let surf2 = BRepTool::face_surface(f2).ok_or("cylinder face 2 has no surface")?;
    let lat2_wire = b.make_wire(&[arc2_bottom.clone(), p2_line.clone(), arc2_top.clone(), tangency2.clone()]);
    let lat2 = b.make_face(surf2, &[lat2_wire]);

    // Blend face: a cylinder of radius r around the centreline.
    let blend_axis = c_xy.to_pnt(0.0, &xh, &yh, &ax, &ref_pt);
    let ax3 = GpAx3::new(blend_axis, GpDir::from_xyz(&ax.xyz()).map_err(|e| e.to_string())?, &GpDir::from_xyz(&xh.xyz()).map_err(|e| e.to_string())?)
        .map_err(|e| format!("fillet_curved: blend cylinder frame: {e}"))?;
    let gcyl = GpCylinder::new(ax3, r).map_err(|e| e.to_string())?;
    let blend_surf: Arc<dyn Surface> = Arc::new(GeomCylinder::new(gcyl));
    let blend_wire = b.make_wire(&[tangency1.clone(), blend_top_arc.clone(), tangency2.clone(), blend_bottom_arc.clone()]);
    let blend_face = b.make_face(blend_surf, &[blend_wire]);

    // End caps: rebuild the bottom and top cap faces.
    let (cap_bottom, cap_top, orig_bottom, orig_top) = rebuild_cylinder_lens_caps(
        solid, &arc1_bottom, &arc2_bottom, &blend_bottom_arc, &arc1_top, &arc2_top, &blend_top_arc, &b,
    )?;

    // Assemble.
    let mut faces: Vec<Face> = Vec::new();
    for f in faces_of(solid) {
        if is_same(&f.0, &f1.0)
            || is_same(&f.0, &f2.0)
            || is_same(&f.0, &orig_bottom.0)
            || is_same(&f.0, &orig_top.0)
        {
            continue;
        }
        faces.push(f);
    }
    faces.push(lat1);
    faces.push(lat2);
    faces.push(blend_face);
    faces.push(cap_bottom);
    faces.push(cap_top);

    let shell = b.make_shell(&faces);
    let solid_out = b.make_solid(&[shell]);
    Ok(solid_out.0)
}

/// Extend the start angle `a0` by a positive sweep that lands on `a1` after
/// passing through `a_through` (all modulo 2π). Returns the end angle `> a0`.
pub(super) fn sweep_through(a0: f64, a1: f64, a_through: f64) -> f64 {
    let two_pi = 2.0 * PI;
    let mut t_through = a_through - a0;
    while t_through < 0.0 {
        t_through += two_pi;
    }
    let mut t1 = a1 - a0;
    while t1 < 0.0 {
        t1 += two_pi;
    }
    if t_through > t1 {
        t1 += two_pi;
    }
    a0 + t1
}

/// Whether a face is planar (used to identify the cap faces of the lens solid).
pub(super) fn face_is_planar_face(f: &Face) -> bool {
    match BRepTool::face_surface(f) {
        Some(s) => is_planar(s.as_ref(), 8, 8, 1e-6),
        None => false,
    }
}

/// Rebuild the two planar lens caps of the parallel-cylinder solid with the
/// blend arcs replacing the filleted edge corner. Returns the rebuilt bottom
/// and top caps together with the original cap faces.
#[allow(clippy::too_many_arguments)]
pub(super) fn rebuild_cylinder_lens_caps(
    solid: &TopoShape,
    arc1_bottom: &Edge,
    arc2_bottom: &Edge,
    blend_bottom: &Edge,
    arc1_top: &Edge,
    arc2_top: &Edge,
    blend_top: &Edge,
    b: &TopoBuilder,
) -> Result<(Face, Face, Face, Face), String> {
    let planars: Vec<Face> = faces_of(solid).into_iter().filter(|f| face_is_planar_face(f)).collect();
    let mut bottom: Option<(Face, f64)> = None;
    let mut top: Option<(Face, f64)> = None;
    for f in &planars {
        let s = BRepTool::face_surface(f).unwrap();
        let p = s.d0(0.0, 0.0);
        if bottom.is_none() || p.z() < bottom.as_ref().unwrap().1 {
            bottom = Some((f.clone(), p.z()));
        }
        if top.is_none() || p.z() > top.as_ref().unwrap().1 {
            top = Some((f.clone(), p.z()));
        }
    }
    let (bf, _) = bottom.ok_or("lens: no bottom cap")?;
    let (tf, _) = top.ok_or("lens: no top cap")?;
    let bs = BRepTool::face_surface(&bf).ok_or("lens: bottom cap has no surface")?;
    let ts = BRepTool::face_surface(&tf).ok_or("lens: top cap has no surface")?;
    let bwire = b.make_wire(&[arc1_bottom.clone(), blend_bottom.clone(), arc2_bottom.clone()]);
    let twire = b.make_wire(&[arc1_top.clone(), blend_top.clone(), arc2_top.clone()]);
    Ok((b.make_face(bs, &[bwire]), b.make_face(ts, &[twire]), bf, tf))
}
