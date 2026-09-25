use super::prelude::*;
use super::*;

pub(super) fn disjoint_result(a: &TopoShape, b: &TopoShape, op: BoolOp) -> BooleanResult {
    let bld = TopoBuilder::new();
    match op {
        BoolOp::Fuse => {
            let comp = bld.make_compound_of(&[a.clone(), b.clone()]);
            BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            }
        }
        BoolOp::Cut => {
            let solid = Solid::wrap(a.clone());
            let shells: Vec<Shell> = shapes_of(a, ShapeType::Shell).into_iter().map(Shell).collect();
            BooleanResult {
                shape: a.clone(),
                solid,
                shells,
                faces: faces_of(a),
                warnings: vec![],
            }
        }
        BoolOp::Common => {
            let comp = bld.make_compound_of(&[]);
            BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Core mesh boolean
// ---------------------------------------------------------------------------

/// Triangulate the retained region of a whole (non-crossing) face.
pub(super) fn mesh_whole_face(face: &Face) -> TriMesh {
    face_full_mesh(face, 24, 24)
}

/// Triangulate the retained region of a crossing face.
///
/// Spherical faces crossing another sphere are trimmed to the analytic cap; the
/// cap's rim is a shared discretization of the intersection circle, so the caps
/// from both solids close the mesh. Other surface types fall back to keeping
/// whole grid cells whose centroid lies on the retained side.
pub(super) fn mesh_crossing_face(
    face: &Face,
    other: &TopoShape,
    keep_inside: bool,
    rim_shared: Option<&IntersectionRim>,
    tol: f64,
) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };

    // Spherical cap trimming.
    if classify_surface(surf.as_ref()) == SurfaceKind::Sphere {
        if let Some((c1, r1)) = face_sphere(face) {
            for (c2, r2) in solid_spheres(other) {
                let d = c1.distance(&c2);
                if d <= 1e-12 {
                    continue;
                }
                let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
                if a > r1 + 1e-9 || a < -r1 - 1e-9 {
                    continue;
                }
                let dir = GpVec::from_pnts(&c1, &c2).divided(d);
                let (pole_vec, vb) = if keep_inside {
                    (dir, (a / r1).clamp(-1.0, 1.0).acos())
                } else {
                    (dir.reversed(), (-a / r1).clamp(-1.0, 1.0).acos())
                };
                let v_lo = FRAC_PI_2 - vb;
                let rim = match rim_shared {
                    Some(r) => r.clone(),
                    None => {
                        let (cc, rc, n) = sphere_sphere_circle(c1, r1, c2, r2).unwrap_or((c1, 0.0, dir));
                        rim_points(cc, rc, &n, 48)
                    }
                };
                return sphere_cap_mesh(c1, r1, pole_vec, v_lo, &rim, 48, 24);
            }
        }
    }

    // General fallback: keep whole cells whose centroid is on the retained side.
    let mesh = surface_grid_mesh(surf.as_ref(), 24, 24);
    let nu = 24;
    let nv = 24;
    let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(surf.as_ref());
    let stride = nv + 1;
    let mut kept: Vec<(usize, usize, usize)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            let uc = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let vc = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let centroid = surf.d0(uc, vc);
            let inside = point_in_solid(centroid, other, tol);
            if inside == keep_inside {
                kept.push((a, b, c));
                kept.push((b, d, c));
            }
        }
    }
    TriMesh { verts: mesh.verts, tris: kept }
}

// ---------------------------------------------------------------------------
// General (non-analytic) curved-face meshing
// ---------------------------------------------------------------------------

/// Sample 3D points along a face's boundary edges (endpoints plus a few
/// interior curve samples). Used to infer a finite UV window when the surface's
/// own range is unbounded (planes, cylinders, …).
pub(super) fn face_boundary_samples(face: &Face) -> Vec<GpPnt> {
    let mut out: Vec<GpPnt> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for wire in face_kids {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for e in wire_kids {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(e) = Edge::wrap(e) {
                if let Some((p0, p1)) = BRepTool::edge_vertices(&e) {
                    out.push(p0);
                    out.push(p1);
                }
                if let Some(c) = BRepTool::edge_curve(&e) {
                    let (t0, t1) = BRepTool::edge_parameters(&e);
                    for k in 0..4usize {
                        let t = if t0.is_finite() && t1.is_finite() && t1 > t0 {
                            t0 + (t1 - t0) * k as f64 / 4.0
                        } else {
                            k as f64 - 1.0
                        };
                        out.push(c.d0(t));
                    }
                }
            }
        }
    }
    out
}

/// Finite parametric window `(u_min, u_max, v_min, v_max)` for a face.
///
/// Uses the surface's own range when finite; otherwise projects boundary edge
/// samples onto the surface and takes the bounding box of the `(u, v)` points.
pub(super) fn face_uv_window_local(face: &Face) -> (f64, f64, f64, f64) {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return (-1.0, 1.0, -1.0, 1.0),
    };
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() {
        return (u1, u2, v1, v2);
    }
    let pts = face_boundary_samples(face);
    if pts.is_empty() {
        return crate::intpatch::sample_bounds(surf.as_ref());
    }
    let mut us = Vec::with_capacity(pts.len());
    let mut vs = Vec::with_capacity(pts.len());
    for p in &pts {
        let (u, v) = crate::intpatch::project_params(surf.as_ref(), p);
        us.push(u);
        vs.push(v);
    }
    let ua = us.iter().copied().fold(f64::INFINITY, f64::min);
    let ub = us.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let va = vs.iter().copied().fold(f64::INFINITY, f64::min);
    let vb = vs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (ua, ub) = if ub - ua < 1e-12 { (ua - 1.0, ub + 1.0) } else { (ua, ub) };
    let (va, vb) = if vb - va < 1e-12 { (va - 1.0, vb + 1.0) } else { (va, vb) };
    (ua, ub, va, vb)
}

/// Mesh a whole (non-crossing) face over its finite trimmed UV window.
pub(super) fn mesh_face_window(face: &Face, nu: usize, nv: usize) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };
    let (u0, u1, v0, v1) = face_uv_window_local(face);
    surface_grid_mesh_uv(surf.as_ref(), u0, u1, v0, v1, nu, nv)
}

/// Axis-aligned bounding box of a face's surface samples (bbox pre-filter for
/// the face-pair intersection loop).
pub(super) fn face_bbox(face: &Face) -> BndBox {
    let mut b = BndBox::new();
    if let Some(s) = BRepTool::face_surface(face) {
        let (u0, u1, v0, v1) = crate::intpatch::sample_bounds(s.as_ref());
        for i in 0..=8usize {
            for j in 0..=8usize {
                let u = u0 + (u1 - u0) * i as f64 / 8.0;
                let v = v0 + (v1 - v0) * j as f64 / 8.0;
                b.add_point(&s.d0(u, v));
            }
        }
    }
    b
}

/// Nearest point of `pts` to `p`.
pub(super) fn nearest_curve_point(pts: &[GpPnt], p: GpPnt) -> Option<GpPnt> {
    let mut best: Option<(f64, GpPnt)> = None;
    for q in pts {
        let d = q.distance(&p);
        if best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, *q));
        }
    }
    best.map(|(_, q)| q)
}

/// Triangulate the retained region of a crossing face using the intersection
/// polylines to pull grid vertices onto the shared boundary curve.
///
/// Sphere faces keep the analytic cap path ([`mesh_crossing_face`]); other
/// surfaces are grid-cell classified (a cell is kept when its centroid lies on
/// the retained side) and grid vertices within `snap_tol` of an intersection
/// curve are snapped onto the nearest curve point so both solids' retained
/// meshes close along the shared intersection.
pub(super) fn mesh_crossing_face_curves(
    face: &Face,
    other: &TopoShape,
    keep_inside: bool,
    curves_on_face: &[Vec<GpPnt>],
    rim_shared: Option<&IntersectionRim>,
    tol: f64,
) -> TriMesh {
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return TriMesh::default(),
    };
    // Spherical caps: analytic trimming only when both solids have the matching
    // sphere face (sphere–sphere). A sphere crossing a non-sphere is meshed by
    // the general grid-cell path below so its boundary snaps to the traced
    // intersection curve, closing the mesh against the other solid's wall.
    let is_sphere = classify_surface(surf.as_ref()) == SurfaceKind::Sphere;
    if is_sphere && !solid_spheres(other).is_empty() {
        return mesh_crossing_face(face, other, keep_inside, rim_shared, tol);
    }

    let (u0, u1, v0, v1) = face_uv_window_local(face);
    let (nu, nv) = (32usize, 32usize);
    let stride = nv + 1;
    let mut verts: Vec<GpPnt> = Vec::with_capacity((nu + 1) * (nv + 1));
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            verts.push(surf.d0(u, v));
        }
    }
    // Snap grid vertices near an intersection curve onto the curve.
    let curve_pts: Vec<GpPnt> = curves_on_face.iter().flatten().copied().collect();
    if !curve_pts.is_empty() {
        let mut lo = verts[0];
        let mut hi = verts[0];
        for p in &verts {
            lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
            hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
        }
        let snap_tol = tol.max(lo.distance(&hi) * 0.03);
        for v in &mut verts {
            if let Some(near) = nearest_curve_point(&curve_pts, *v) {
                if near.distance(v) <= snap_tol {
                    *v = near;
                }
            }
        }
    }
    let mut tris: Vec<(usize, usize, usize)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let a = i * stride + j;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            let uc = u0 + (u1 - u0) * (i as f64 + 0.5) / nu as f64;
            let vc = v0 + (v1 - v0) * (j as f64 + 0.5) / nv as f64;
            let centroid = surf.d0(uc, vc);
            let inside = point_in_solid_curved(centroid, other, tol);
            if inside == keep_inside {
                tris.push((a, b, c));
                tris.push((b, d, c));
            }
        }
    }
    TriMesh { verts, tris }
}

/// Compute the retained triangle mesh of the boolean, using the per-face
/// intersection polylines to close the mesh along the shared intersection
/// curves of the general path.
pub(super) fn boolean_mesh_curves(
    a: &TopoShape,
    b: &TopoShape,
    op: BoolOp,
    tol: f64,
    pair_curves: &[Vec<Vec<Vec<GpPnt>>>],
) -> Result<TriMesh, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face_general(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face_general(f, a, tol)).collect();

    // Shared sphere-sphere intersection rim (reused by both caps).
    let mut rim_shared: Option<IntersectionRim> = None;
    'outer: for fa_i in &fa {
        let Some(surf) = BRepTool::face_surface(fa_i) else { continue };
        if classify_surface(surf.as_ref()) != SurfaceKind::Sphere {
            continue;
        }
        let (c1, r1) = match face_sphere(fa_i) {
            Some(x) => x,
            None => continue,
        };
        for fb_j in &fb {
            let (c2, r2) = match face_sphere(fb_j) {
                Some(x) => x,
                None => continue,
            };
            if let Some((cc, rc, n)) = sphere_sphere_circle(c1, r1, c2, r2) {
                rim_shared = Some(rim_points(cc, rc, &n, 48));
                break 'outer;
            }
        }
    }

    let centroid_a = solid_centroid(a);
    let centroid_b = solid_centroid(b);
    let flip_b = op == BoolOp::Cut;

    let mut parts: Vec<TriMesh> = Vec::new();
    for (idx, (f, region)) in fa.iter().zip(&regions_a).enumerate() {
        if select_keep(*region, true, op) {
            let mut m = match region {
                FaceRegion::On => {
                    let curves: Vec<Vec<GpPnt>> = pair_curves[idx].iter().flatten().cloned().collect();
                    mesh_crossing_face_curves(f, b, keep_inside(true, op), &curves, rim_shared.as_ref(), tol)
                }
                _ => mesh_face_window(f, 24, 24),
            };
            orient_outward(&mut m, &centroid_a);
            parts.push(m);
        }
    }
    for (idx, (f, region)) in fb.iter().zip(&regions_b).enumerate() {
        if select_keep(*region, false, op) {
            let mut m = match region {
                FaceRegion::On => {
                    let curves: Vec<Vec<GpPnt>> = pair_curves
                        .iter()
                        .map(|row| row[idx].clone())
                        .flatten()
                        .collect();
                    mesh_crossing_face_curves(f, a, keep_inside(false, op), &curves, rim_shared.as_ref(), tol)
                }
                _ => mesh_face_window(f, 24, 24),
            };
            orient_outward(&mut m, &centroid_b);
            if flip_b {
                for t in &mut m.tris {
                    let (i, j, k) = *t;
                    *t = (i, k, j);
                }
            }
            parts.push(m);
        }
    }

    let mut all = TriMesh::default();
    for p in &parts {
        let offset = all.verts.len();
        all.verts.extend(p.verts.iter().copied());
        for &(i, j, k) in &p.tris {
            all.tris.push((offset + i, offset + j, offset + k));
        }
    }
    let weld_tol = tol.max(1e-9);
    Ok(weld_mesh(&all, weld_tol))
}

/// Compute the retained triangle mesh of the boolean.
pub(super) fn boolean_mesh(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TriMesh, String> {
    let fa = faces_of(a);
    let fb = faces_of(b);
    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face(f, a, tol)).collect();

    // Pre-compute a shared intersection-circle rim when a sphere face of A
    // crosses a sphere face of B (used by both caps).
    let mut rim_shared: Option<IntersectionRim> = None;
    'outer: for fa_i in &fa {
        let Some(surf) = BRepTool::face_surface(fa_i) else { continue };
        if classify_surface(surf.as_ref()) != SurfaceKind::Sphere {
            continue;
        }
        let (c1, r1) = match face_sphere(fa_i) {
            Some(x) => x,
            None => continue,
        };
        for fb_j in &fb {
            let (c2, r2) = match face_sphere(fb_j) {
                Some(x) => x,
                None => continue,
            };
            if let Some((cc, rc, n)) = sphere_sphere_circle(c1, r1, c2, r2) {
                rim_shared = Some(rim_points(cc, rc, &n, 48));
                break 'outer;
            }
        }
    }

    let mut parts: Vec<TriMesh> = Vec::new();
    for (f, region) in fa.iter().zip(&regions_a) {
        if select_keep(*region, true, op) {
            match region {
                FaceRegion::On => {
                    let keep_in = keep_inside(true, op);
                    parts.push(mesh_crossing_face(f, b, keep_in, rim_shared.as_ref(), tol));
                }
                _ => parts.push(mesh_whole_face(f)),
            }
        }
    }
    for (f, region) in fb.iter().zip(&regions_b) {
        if select_keep(*region, false, op) {
            match region {
                FaceRegion::On => {
                    let keep_in = keep_inside(false, op);
                    parts.push(mesh_crossing_face(f, a, keep_in, rim_shared.as_ref(), tol));
                }
                _ => parts.push(mesh_whole_face(f)),
            }
        }
    }

    // Concatenate all fragments into one mesh and weld shared vertices.
    let mut all = TriMesh::default();
    for p in &parts {
        let offset = all.verts.len();
        all.verts.extend(p.verts.iter().copied());
        for &(i, j, k) in &p.tris {
            all.tris.push((offset + i, offset + j, offset + k));
        }
    }
    let weld_tol = tol.max(1e-9);
    Ok(weld_mesh(&all, weld_tol))
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Boolean on solids whose faces include curved surfaces.
///
/// Delegates to `crate::bop_builder::boolean` when both inputs are planar.
/// For curved inputs it classifies each face, keeps whole faces that do not
/// cross the boundary, trims crossing faces to their retained region (analytic
/// spherical caps for sphere–sphere), welds the retained triangles into a
/// closed mesh, and rebuilds a solid via `mesh_to_brep`.
pub fn curved_boolean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    // Planar inputs → exact polygon boolean.
    if all_faces_planar(a) && all_faces_planar(b) {
        return crate::bop_builder::boolean(a, b, op, tol);
    }

    // Disjoint shortcut.
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return Ok(disjoint_result(a, b, op));
    }

    let fa = faces_of(a);
    let fb = faces_of(b);
    if fa.is_empty() || fb.is_empty() {
        return Err("curved_boolean: input has no faces".into());
    }

    let regions_a: Vec<FaceRegion> = fa.iter().map(|f| classify_face(f, b, tol)).collect();
    let regions_b: Vec<FaceRegion> = fb.iter().map(|f| classify_face(f, a, tol)).collect();
    let has_on = regions_a.contains(&FaceRegion::On) || regions_b.contains(&FaceRegion::On);

    // When no face crosses the boundary, the result is assembled from the whole
    // (original) kept faces — curved faces stay analytic, so `shape_volume` is
    // reliable for the curved result.
    if !has_on {
        let mut faces: Vec<Face> = Vec::new();
        for (f, region) in fa.iter().zip(&regions_a) {
            if select_keep(*region, true, op) {
                faces.push(f.clone());
            }
        }
        for (f, region) in fb.iter().zip(&regions_b) {
            if select_keep(*region, false, op) {
                faces.push(f.clone());
            }
        }
        if faces.is_empty() {
            let bld = TopoBuilder::new();
            let comp = bld.make_compound_of(&[]);
            return Ok(BooleanResult {
                shape: comp.0,
                solid: None,
                shells: vec![],
                faces: vec![],
                warnings: vec![],
            });
        }
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&faces);
        let solid = bld.make_solid(&[shell.clone()]);
        return Ok(BooleanResult {
            shape: solid.0.clone(),
            solid: Some(solid),
            shells: vec![shell],
            faces,
            warnings: vec![],
        });
    }

    // Crossing faces: build the closed mesh and rebuild a faceted BRep solid.
    let mesh = boolean_mesh(a, b, op, tol)?;
    let vol = mesh_volume(&mesh);
    if mesh.tris.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return Ok(BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![format!("empty curved boolean (mesh volume {vol:.4})")],
        });
    }

    let smesh = crate::mesh::ShapeMesh {
        vertices: mesh.verts.clone(),
        triangles: mesh
            .tris
            .iter()
            .map(|&(a, b, c)| occt_core::poly::triangulation::Triangle::new(a, b, c))
            .collect(),
        source_shape: ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&smesh);
    let shape = brep.solid.clone().map(|s| s.0).unwrap_or_else(|| brep.shell.0.clone());
    Ok(BooleanResult {
        shape,
        solid: brep.solid,
        shells: vec![brep.shell],
        faces: brep.faces,
        warnings: vec![format!("curved boolean mesh volume ≈ {vol:.4}")],
    })
}

/// The enclosed volume of the curved boolean result (exact for the closed
/// welded mesh; for planar inputs delegates to the planar boolean and meshes
/// its result).
pub fn curved_boolean_volume(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<f64, String> {
    if all_faces_planar(a) && all_faces_planar(b) {
        let r = crate::bop_builder::boolean(a, b, op, tol)?;
        return Ok(crate::shape_mesh::shape_volume(&r.shape, 0.02));
    }
    let bbox_a = crate::bbox_from_geometry::shape_bbox(a);
    let bbox_b = crate::bbox_from_geometry::shape_bbox(b);
    if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) {
        return match op {
            BoolOp::Fuse => Ok(crate::shape_mesh::shape_volume(a, 0.02) + crate::shape_mesh::shape_volume(b, 0.02)),
            BoolOp::Cut => Ok(crate::shape_mesh::shape_volume(a, 0.02)),
            BoolOp::Common => Ok(0.0),
        };
    }
    let mesh = boolean_mesh(a, b, op, tol)?;
    Ok(mesh_volume(&mesh))
}

// ---------------------------------------------------------------------------
// General (non-analytic) curved boolean
// ---------------------------------------------------------------------------

/// Build a `BooleanResult` from the retained whole faces when no face crosses
/// the boundary (the analytic faces are preserved verbatim).
pub(super) fn result_from_kept_faces(
    fa: &[Face],
    fb: &[Face],
    regions_a: &[FaceRegion],
    regions_b: &[FaceRegion],
    op: BoolOp,
) -> BooleanResult {
    let mut faces: Vec<Face> = Vec::new();
    for (f, region) in fa.iter().zip(regions_a) {
        if select_keep(*region, true, op) {
            faces.push(f.clone());
        }
    }
    for (f, region) in fb.iter().zip(regions_b) {
        if select_keep(*region, false, op) {
            faces.push(f.clone());
        }
    }
    if faces.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![],
        };
    }
    let bld = TopoBuilder::new();
    let shell = bld.make_shell(&faces);
    let solid = bld.make_solid(&[shell.clone()]);
    BooleanResult {
        shape: solid.0.clone(),
        solid: Some(solid),
        shells: vec![shell],
        faces,
        warnings: vec![],
    }
}

/// Rebuild a faceted `BooleanResult` solid from a closed retained mesh.
pub(super) fn result_from_mesh(mesh: TriMesh) -> BooleanResult {
    let vol = mesh_volume(&mesh);
    if mesh.tris.is_empty() {
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[]);
        return BooleanResult {
            shape: comp.0,
            solid: None,
            shells: vec![],
            faces: vec![],
            warnings: vec![format!("empty curved boolean (mesh volume {vol:.4})")],
        };
    }
    let smesh = crate::mesh::ShapeMesh {
        vertices: mesh.verts.clone(),
        triangles: mesh
            .tris
            .iter()
            .map(|&(a, b, c)| occt_core::poly::triangulation::Triangle::new(a, b, c))
            .collect(),
        source_shape: ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&smesh);
    let shape = brep.solid.clone().map(|s| s.0).unwrap_or_else(|| brep.shell.0.clone());
    BooleanResult {
        shape,
        solid: brep.solid,
        shells: vec![brep.shell],
        faces: brep.faces,
        warnings: vec![format!("curved boolean mesh volume ≈ {vol:.4}")],
    }
}

/// Face–face intersection polylines for the general path (pairs whose bounding
/// boxes do not overlap are skipped).
pub(super) fn general_pair_curves(fa: &[Face], fb: &[Face], tol: f64) -> Vec<Vec<Vec<Vec<GpPnt>>>> {
    let mut pair_curves: Vec<Vec<Vec<Vec<GpPnt>>>> = vec![vec![Vec::new(); fb.len()]; fa.len()];
    for (i, f_i) in fa.iter().enumerate() {
        let Some(sa) = BRepTool::face_surface(f_i) else { continue };
        let ba = face_bbox(f_i);
        for (j, f_j) in fb.iter().enumerate() {
            let Some(sb) = BRepTool::face_surface(f_j) else { continue };
            let bb = face_bbox(f_j);
            if !ba.is_void() && !bb.is_void() && ba.is_out_box(&bb) {
                continue;
            }
            let curves = crate::intpatch::intersection_curve_points(sa.as_ref(), sb.as_ref(), tol, 40);
            if !curves.is_empty() {
                pair_curves[i][j] = curves;
            }
        }
    }
    pair_curves
}
