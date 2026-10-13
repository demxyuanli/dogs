use super::prelude::*;


/// Result of a feature operation: the modified solid and its volume.
#[derive(Debug, Clone)]

pub struct FeatResult {
    pub shape: TopoShape,
    pub volume: f64,
}

/// Draft one planar face of `solid`: rotate the face's supporting plane about
/// the hinge (the intersection line of the face plane and the pivot plane) by
/// `angle_deg`, replacing the face's surface geometry in a copy of the solid.
///
/// The face is located by `face_pts[0]` (the planar face whose surface passes
/// nearest that point). `angle_deg = 0` returns the solid unchanged. The
/// boundary wire is kept as-is, so the result is the "drafted" solid only in
/// the sense that the face now leans by the draft angle.
pub fn draft(solid: &Solid, face_pts: &[GpPnt], angle_deg: f64, pivot: &GpPln, tol: f64) -> Result<FeatResult, String> {
    if face_pts.is_empty() {
        return Err("draft: face_pts must contain the point on the face to draft".into());
    }
    let angle = angle_deg.to_radians();
    let identity = GpTrsf::identity();
    let copy = crate::shape_ops::transformed_copy(&solid.0, &identity)?;
    if angle.abs() < 1e-12 {
        let volume = crate::shape_mesh::shape_volume(&copy, tol.max(1e-4));
        return Ok(FeatResult { shape: copy, volume });
    }
    let p0 = face_pts[0];
    let face = face_containing_point(&Solid::wrap(copy.clone()).ok_or("draft: result is not a solid")?, &p0, tol)
        .map_err(|e| format!("draft: {e}"))?;
    let fpln = crate::brep_surface::face_plane(&face).ok_or("draft: the target face is not planar")?;
    let hinge = draft_hinge_axis(&fpln, pivot)?;
    let rotated = fpln.rotated(&hinge, angle);
    GeometryRegistry::global().set_face(
        &face.0,
        FaceGeom { surface: Arc::new(GeomPlane::new(rotated)), tolerance: 0.0, natural_restriction: true },
    );
    let volume = crate::shape_mesh::shape_volume(&copy, tol.max(1e-4));
    Ok(FeatResult { shape: copy, volume })
}

/// Groove: revolve the 2D `profile` (x = radius ≥ 0, y = height along the
/// axis) around `axis` into an annular tool and SUBTRACT it from `solid`.
pub fn groove(solid: &Solid, profile: &[GpPnt2d], axis: &GpAx1, steps: usize, tol: f64) -> Result<FeatResult, String> {
    require_min_profile(profile, "groove")?;
    let tool = revolve_profile_about(profile, steps, axis)?;
    boolean_feature(solid, &tool, BoolOp::Cut, tol)
}

/// Neck: revolve the 2D `profile` around `axis` into a solid tool and FUSE it
/// onto `solid`.
pub fn neck(solid: &Solid, profile: &[GpPnt2d], axis: &GpAx1, steps: usize, tol: f64) -> Result<FeatResult, String> {
    require_min_profile(profile, "neck")?;
    let tool = revolve_profile_about(profile, steps, axis)?;
    boolean_feature(solid, &tool, BoolOp::Fuse, tol)
}

/// Rib: a thin plate feature standing `height` above `plane`. The profile's
/// in-plane bounding rectangle is the footprint; a dimension thinner than
/// `thickness` is padded (centered) so a line-like profile becomes a plate of
/// the rib's wall thickness. The plate is fused onto `solid`.
pub fn rib(solid: &Solid, profile: &[GpPnt2d], plane: &GpPln, thickness: f64, height: f64, tol: f64) -> Result<FeatResult, String> {
    require_min_profile(profile, "rib")?;
    require_positive(thickness, "thickness", "rib")?;
    require_positive(height, "height", "rib")?;
    let (origin, xd, yd, nd) = plane_frame(plane);
    let (umin, umax, vmin, vmax) =
        profile_bounds(profile).ok_or("rib: degenerate profile footprint")?;
    // Pad a footprint that is thinner than the rib wall thickness.
    let (du, dv) = (umax - umin, vmax - vmin);
    let (du, dv) = (du.max(thickness), dv.max(thickness));
    let u0 = if du > umax - umin { (umin + umax - du) / 2.0 } else { umin };
    let v0 = if dv > vmax - vmin { (vmin + vmax - dv) / 2.0 } else { vmin };

    let box_s = crate::primitives::BRepPrimBox::make_box(du, dv, height);
    let mut t = GpTrsf::identity();
    t.matrix = GpMat::from_cols(xd.xyz(), yd.xyz(), nd.xyz());
    t.shape = TrsfForm::CompoundTrsf;
    t.loc = origin
        .translated_vec(&xd.multiplied_scalar(u0))
        .translated_vec(&yd.multiplied_scalar(v0))
        .coord;
    let tool = crate::shape_ops::transformed_copy(&box_s.solid.0, &t)?;
    boolean_feature(solid, &tool, BoolOp::Fuse, tol)
}

/// Boss through all: fuse a cylinder of `radius` whose axis passes through
/// `center` (axis +Z) and whose height spans the solid's bounding-box height,
/// so it pierces the solid from below to above.
pub fn boss_thru_all(solid: &Solid, center: &GpPnt, radius: f64, tol: f64) -> Result<FeatResult, String> {
    require_positive(radius, "radius", "boss_thru_all")?;
    let (z0, z1) = solid_z_span(&solid.0).map_err(|e| format!("boss_thru_all: {e}"))?;
    let margin = (z1 - z0).max(radius).max(1.0);
    let height = (z1 - z0) + 2.0 * margin;
    let tool = translated_mesh_cylinder(center, radius, height, z0 - margin)?;
    boolean_feature(solid, &tool, BoolOp::Fuse, tol)
}

/// Boss: fuse a cylindrical boss of `radius` and `height` onto `solid`. The
/// cylinder's axis is +Z through `center`, with its base at the solid's top
/// face (`center`'s height). Source: `BRepFeat_MakePrism` boss family.
pub fn boss(solid: &Solid, center: &GpPnt, radius: f64, height: f64, tol: f64) -> Result<FeatResult, String> {
    require_positive(radius, "radius", "boss")?;
    require_positive(height, "height", "boss")?;
    // Base the boss on the solid's top face (the cylinder axis passes through
    // `center`'s x/y, the base sits at the top z). Placing it at the bottom z
    // hid the cylinder inside the solid, so a fuse swallowed it entirely.
    let (_, z1) = solid_z_span(&solid.0).map_err(|e| format!("boss: {e}"))?;
    let tool = translated_mesh_cylinder(center, radius, height, z1)?;
    boolean_feature(solid, &tool, BoolOp::Fuse, tol)
}

/// Analytic volume of the rib plate a [`rib`] call would build from `profile`
/// on a plane with the given `thickness` and `height`: the padded in-plane
/// footprint times the standing height.
pub fn rib_volume(profile: &[GpPnt2d], thickness: f64, height: f64) -> Result<f64, String> {
    require_min_profile(profile, "rib_volume")?;
    require_positive(thickness, "thickness", "rib_volume")?;
    require_positive(height, "height", "rib_volume")?;
    let (umin, umax, vmin, vmax) = profile_bounds(profile).ok_or("rib_volume: degenerate profile")?;
    let du = (umax - umin).max(thickness);
    let dv = (vmax - vmin).max(thickness);
    Ok(du * dv * height)
}

/// `|after.volume − volume(before)|`, the material added or removed by the
/// feature. Both volumes use [`solid_volume`] (`BRepGProp::VolumeProperties`),
/// the divergence-theorem integral over the shape's actual faces.
pub fn feature_before_after_delta(before: &Solid, after: &FeatResult) -> f64 {
    let v0 = solid_volume(&before.0);
    (after.volume - v0).abs()
}

/// Whether a feature adds material to the target solid (fuse) or removes it
/// (cut). Mirrors the two boolean families behind `BRepFeat_MakeRevol`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatKind {
    /// Material is added: neck, rib, boss.
    Additive,
    /// Material is removed: groove, pocket-like cuts.
    Subtractive,
}

impl FeatResult {
    /// The volume this feature added (`> 0`) or removed (`< 0`), relative to
    /// the solid it was applied to. Negative means a cut (groove).
    pub fn volume_delta(&self, before: &Solid) -> f64 {
        self.volume - solid_volume(&before.0)
    }

    /// `after / before` volume ratio: `> 1` for an additive feature,
    /// `< 1` for a subtractive one.
    pub fn volume_ratio(&self, before: &Solid) -> f64 {
        let v0 = solid_volume(&before.0);
        if v0.abs() < 1e-12 { f64::NAN } else { self.volume / v0 }
    }

    /// Human-readable classification of the feature's effect.
    pub fn kind(&self, before: &Solid) -> FeatKind {
        if self.volume > solid_volume(&before.0) {
            FeatKind::Additive
        } else {
            FeatKind::Subtractive
        }
    }
}

/// The draft hinge line: the intersection of a face's supporting plane and the
/// pivot plane. The drafted face rotates about this line (an axis through the
/// intersection with the intersection direction). Source: the pivot logic of
/// `BRepFeat_MakeDPrism`.
pub fn draft_hinge_axis(face_plane: &GpPln, pivot: &GpPln) -> Result<GpAx1, String> {
    let (p, d) = crate::face_face::plane_plane_intersection(face_plane, pivot)
        .ok_or("draft: face plane and pivot plane are parallel (no hinge line)")?;
    let d = d.normalized();
    let dir = GpDir::from_vec(&d).map_err(|_| "draft: degenerate hinge direction")?;
    Ok(GpAx1::new(p, dir))
}

/// The planar face of `solid` whose surface passes nearest to `p`, if that
/// distance is within `tol`. Used by [`draft`] to select the face to lean.
pub fn face_containing_point(solid: &Solid, p: &GpPnt, tol: f64) -> Result<crate::shape::Face, String> {
    faces_of(&solid.0)
        .into_iter()
        .find(|f| face_point_distance(f, p) < tol.max(1e-6))
        .ok_or_else(|| format!("no face contains point {p:?}"))
}

/// Analytic volume of the solid swept by revolving a closed 2D `profile`
/// around the axis, by Pappus's centroid theorem (`sweep_revolve`). Useful for
/// predicting a groove's material removal or a neck's added volume.
pub fn revolved_profile_volume(profile: &[GpPnt2d], steps: usize) -> Result<f64, String> {
    require_min_profile(profile, "revolved_profile_volume")?;
    let rev = crate::sweep_revolve::revolve_polyline_around_z(profile, steps)?;
    Ok(crate::sweep_revolve::revolved_volume(&rev))
}

/// Combine `solid` with an already-built `tool` solid by the given boolean
/// operation, returning the feature result. The public face of [`boolean_feature`]:
/// useful for feature types built from a custom tool (a `BRepFeat_MakePrism`
/// profile, for example).
pub fn feature_apply(solid: &Solid, tool: &TopoShape, op: BoolOp, tol: f64) -> Result<FeatResult, String> {
    boolean_feature(solid, tool, op, tol)
}

/// Validate a revolve `profile` for `groove`/`neck`: enough points, non-negative
/// radii, and at least four tessellation steps (mirrors `sweep_revolve`'s own
/// checks, surfaced before any geometry is built).
pub fn validate_revolve_profile(profile: &[GpPnt2d], steps: usize) -> Result<(), String> {
    require_min_profile(profile, "revolve")?;
    if steps < 4 {
        return Err("revolve: steps must be >= 4".into());
    }
    if profile.iter().any(|p| p.x() < 0.0) {
        return Err("revolve: profile radii must be >= 0".into());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Internals
// ---------------------------------------------------------------------------

/// The `(z_min, z_max)` extent of a shape's bounding box.
pub(super) fn solid_z_span(shape: &TopoShape) -> Result<(f64, f64), String> {
    let (_, _, _, _, z0, z1) = crate::bbox_from_geometry::shape_bbox(shape)
        .get()
        .ok_or("shape has no bounding box")?;
    Ok((z0, z1))
}

/// Validate that a revolve/rib profile has at least the three points needed to
/// define a closed contour.
pub(super) fn require_min_profile(profile: &[GpPnt2d], feature: &str) -> Result<(), String> {
    if profile.len() < 3 {
        Err(format!("{feature}: profile needs at least 3 points"))
    } else {
        Ok(())
    }
}

/// Validate that a feature dimension is strictly positive.
pub(super) fn require_positive(value: f64, name: &str, feature: &str) -> Result<(), String> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(format!("{feature}: {name} must be positive"))
    }
}

/// The right-handed frame of a plane: `(origin, x_dir, y_dir, normal)`, each
/// direction as a unit vector.
pub(super) fn plane_frame(plane: &GpPln) -> (GpPnt, GpVec, GpVec, GpVec) {
    let pos = plane.position();
    (
        plane.location(),
        GpVec::from_xyz(pos.x_direction().xyz()),
        GpVec::from_xyz(pos.y_direction().xyz()),
        GpVec::from_xyz(pos.direction().xyz()),
    )
}

/// Axis-aligned bounds of a 2D profile: `(u_min, u_max, v_min, v_max)`.
pub(super) fn profile_bounds(profile: &[GpPnt2d]) -> Option<(f64, f64, f64, f64)> {
    let mut umin = f64::INFINITY;
    let mut umax = f64::NEG_INFINITY;
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for p in profile {
        umin = umin.min(p.x());
        umax = umax.max(p.x());
        vmin = vmin.min(p.y());
        vmax = vmax.max(p.y());
    }
    if umin.is_finite() && umax > umin && vmax > vmin {
        Some((umin, umax, vmin, vmax))
    } else {
        None
    }
}

/// Whether every face of `shape` carries a planar surface (eligible for the
/// exact planar boolean).
pub(super) fn all_faces_planar(shape: &TopoShape) -> bool {
    let fa = faces_of(shape);
    if fa.is_empty() {
        return false;
    }
    fa.iter().all(|f| match crate::brep_tool::BRepTool::face_surface(f) {
        Some(s) => crate::brep_surface::classify_surface(s.as_ref())
            == crate::brep_surface::SurfaceKind::Plane,
        None => false,
    })
}

/// Combine `solid` with the `tool` solid.
///
/// T-41: OCCT's only boolean is `BRepAlgoAPI_*` / `BOPAlgo_BOP` (audit A5: the
/// port's `bop_curved` voxel/mesh boolean has no counterpart). Both the planar
/// and the curved case run the exact `BOPAlgo_BOP` path.
pub(super) fn boolean_feature(solid: &Solid, tool: &TopoShape, op: BoolOp, tol: f64) -> Result<FeatResult, String> {
    let r = crate::bop_builder::boolean(&solid.0, tool, op, tol.max(1e-9))
        .map_err(|e| format!("feature boolean: {e}"))?;
    let shape = r.solid.map(|s| s.0).unwrap_or(r.shape);
    // `BRepGProp::VolumeProperties` over the boolean's own faces.
    let volume = solid_volume(&shape);
    Ok(FeatResult { shape, volume })
}

/// Revolve `profile` around the Z axis and align the tool so its axis is the
/// given `axis` (a rotation mapping Z onto the axis direction, then a
/// translation to the axis location).
///
/// UNPORTED workaround (`BRepFeat_MakeRevol` has no such step: OCCT takes the
/// tool `TopoDS_Shape` as given and hands it to `BOPAlgo_BOP`): the analytic
/// revolve solid `sweep_revolve::revolve_polyline_around_z` returns cannot yet
/// be combined by `crate::bop_builder::boolean`, so the tool is rebuilt as a
/// planar-faced (lathe) mesh from the revolution's ring vertices. Measured
/// 2026-10-13 on the `brepfeat` baselines:
///
/// * box ∪ analytic neck (touching at z = 1) -> one solid, but the tool is
///   dropped: volume stays the box's 4.000000 (expected 4.489835);
/// * box ∪ analytic through-boss -> `res_solid = false` and the whole
///   cylinder is added untrimmed: 2.706858 (expected 2.565487);
/// * faceted cylinder − analytic groove ring (`Cut`) -> `res_solid = false`,
///   and the cut *adds* material: 9.337869 > 9.317486.
///
/// With a faceted tool all three agree with the baseline. The gap is in the
/// port's `BOPAlgo_BOP` path (audit A5 / T-79, `bop_curved`), not in the
/// mesher: `GeomCylinder::d1` does return the analytic partials
/// (`occt-geom/src/cylinder.rs:20-22` -> `elib::surface_eval::cylinder_d1`).
/// Consecutive rings form quad strips; a ring on the axis collapses to a fan,
/// so closed profiles produce closed solids.
pub(super) fn revolve_profile_about(profile: &[GpPnt2d], steps: usize, axis: &GpAx1) -> Result<TopoShape, String> {
    let rev = crate::sweep_revolve::revolve_polyline_around_z(profile, steps)?;
    let mut verts: Vec<GpPnt> = Vec::new();
    let mut tris: Vec<occt_core::poly::triangulation::Triangle> = Vec::new();
    for w in rev.rings.windows(2) {
        let (r0, r1) = (&w[0], &w[1]);
        let n = r0.len().min(r1.len());
        if n < 2 {
            continue;
        }
        let base = verts.len();
        for v in r0 {
            verts.push(crate::brep_tool::BRepTool::vertex_point(v));
        }
        for v in r1 {
            verts.push(crate::brep_tool::BRepTool::vertex_point(v));
        }
        for j in 0..n - 1 {
            let (a, b) = (base + j, base + j + 1);
            let (c, d) = (base + n + j, base + n + j + 1);
            // Single shared diagonal (a, d) so adjacent strip quads share the
            // vertical edges (r0[j+1], r1[j+1]) and (r0[j], r1[j]) with their
            // neighbours — a closed quad strip.
            tris.push(occt_core::poly::triangulation::Triangle::new(a, b, d));
            tris.push(occt_core::poly::triangulation::Triangle::new(a, d, c));
        }
    }
    // Drop degenerate triangles (axis rings collapse to fans).
    tris.retain(|t| {
        verts[t.n0].distance(&verts[t.n1]) > 1e-9
            && verts[t.n1].distance(&verts[t.n2]) > 1e-9
            && verts[t.n0].distance(&verts[t.n2]) > 1e-9
    });
    if tris.is_empty() {
        return Err("revolve: profile produces no faces".into());
    }
    let mesh = crate::mesh::ShapeMesh {
        vertices: verts,
        triangles: tris,
        source_shape: crate::abs::ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&mesh);
    let tool = brep
        .solid
        .ok_or("revolve: profile mesh is not a closed solid (profile must close)")?
        .0;
    let t = align_z_to_axis(axis)?;
    crate::shape_ops::transformed_copy(&tool, &t)
}

/// Build a faceted cylinder of `radius` × `height` (axis +Z) translated so its
/// base sits at `z_base` and its axis passes through `center`.
pub(super) fn translated_mesh_cylinder(
    center: &GpPnt,
    radius: f64,
    height: f64,
    z_base: f64,
) -> Result<TopoShape, String> {
    let cyl = mesh_cylinder(radius, height, 24);
    crate::shape_ops::translated_copy(&cyl, &GpVec::new(center.x(), center.y(), z_base))
}

/// Divergence-theorem volume of a triangle mesh (`|Σ a·(b×c)| / 6`).
pub(super) fn mesh_volume_signed(mesh: &crate::mesh::ShapeMesh) -> f64 {
    let mut vol = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t.n0].coord;
        vol += a.dot_cross(&mesh.vertices[t.n1].coord, &mesh.vertices[t.n2].coord);
    }
    (vol / 6.0).abs()
}

/// Faceted cylinder solid (axis +Z) from an exact triangle mesh, so all faces
/// are planar and mesh reliably. UNPORTED workaround for the same
/// `bop_builder::boolean` gap as [`revolve_profile_about`]: OCCT's
/// `BRepPrimAPI_MakeCylinder` tool is analytic, and its fuse with a box is not
/// handled by the port yet (`BRepPrimCylinder::make_cylinder` + `Fuse` here
/// yields the untrimmed union, 2.706858 instead of 2.565487).
pub(super) fn mesh_cylinder(radius: f64, height: f64, slices: usize) -> TopoShape {
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..=slices {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let (x, y) = (radius * theta.cos(), radius * theta.sin());
        vertices.push(GpPnt::new(x, y, 0.0));
        vertices.push(GpPnt::new(x, y, height));
    }
    for i in 0..slices {
        let (a, b, c, d) = (2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3);
        triangles.push(occt_core::poly::triangulation::Triangle::new(a, b, c));
        triangles.push(occt_core::poly::triangulation::Triangle::new(b, d, c));
    }
    let (top_idx, bot_idx) = (vertices.len(), vertices.len() + 1);
    vertices.push(GpPnt::new(0.0, 0.0, height));
    vertices.push(GpPnt::new(0.0, 0.0, 0.0));
    for i in 0..slices {
        let (tb, tt) = (2 * i + 1, 2 * i + 3);
        let (bb, bt) = (2 * i, 2 * i + 2);
        triangles.push(occt_core::poly::triangulation::Triangle::new(top_idx, tt, tb));
        triangles.push(occt_core::poly::triangulation::Triangle::new(bot_idx, bb, bt));
    }
    let mesh = crate::mesh::ShapeMesh {
        vertices,
        triangles,
        source_shape: crate::abs::ShapeType::Solid,
    };
    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&mesh);
    brep.solid.expect("mesh_cylinder: mesh is closed").0
}

/// Transform mapping the +Z axis through the origin onto the line `axis`
/// (direction and location).
///
/// The rotation is computed with Rodrigues' formula about the axis
/// `Z × direction` (the perpendicular that tilts Z onto the axis direction);
/// when Z and the direction are anti-parallel a 180° rotation about X is used,
/// and when parallel the rotation is the identity. The translation then carries
/// the origin to the axis location. The composed transform is exactly what
/// `TopoDS_Shape::Move` would apply to relocate a Z-revolved tool onto an
/// arbitrary feature axis.
pub(super) fn align_z_to_axis(axis: &GpAx1) -> Result<GpTrsf, String> {
    let z = GpVec::new(0.0, 0.0, 1.0);
    let d = GpVec::from_xyz(axis.direction().xyz()).normalized();
    let cos = z.dot(&d).clamp(-1.0, 1.0);
    let mut rot = GpTrsf::identity();
    let cross = z.crossed(&d);
    if cross.xyz().square_modulus() > 1e-30 {
        let ax = GpAx1::new(GpPnt::zero(), GpDir::from_vec(&cross).map_err(|_| "degenerate feature axis")?);
        rot.set_rotation_ax1(&ax, cos.acos()).map_err(|e| e.to_string())?;
    } else if cos < 0.0 {
        let ax = GpAx1::new(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).expect("X axis"));
        rot.set_rotation_ax1(&ax, std::f64::consts::PI).map_err(|e| e.to_string())?;
    }
    let mut tr = GpTrsf::identity();
    tr.set_translation_vec(&GpVec::from_xyz(&axis.location().coord));
    Ok(tr.multiplied(&rot))
}

/// Voxel resolution from the combined bounding-box span so a cell is about
/// `tol` wide (clamped for speed and robustness).
///
/// The target cell width is `tol`, but the resolution is clamped to
/// `[16, 64]` cells across the largest span: fine enough to resolve the tool
/// features, coarse enough to keep the voxel classification tractable.
pub(super) fn resolution_for(solid: &Solid, tool: &TopoShape, tol: f64) -> usize {
    let mut bb = crate::bbox_from_geometry::shape_bbox(&solid.0);
    bb.add_box(&crate::bbox_from_geometry::shape_bbox(tool));
    let (x0, x1, y0, y1, z0, z1) = bb.get().unwrap_or((0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
    let span = (x1 - x0).max(y1 - y0).max(z1 - z0).max(1e-9);
    let tol = if tol.is_finite() && tol > 0.0 { tol } else { span / 40.0 };
    let r = (span / tol.max(span / 64.0)).ceil() as usize;
    r.clamp(16, 64)
}

/// Distance from `p` to the nearest point of the face's surface. For planar
/// faces the perpendicular plane distance is used (independent of the unbounded
/// UV window); other surfaces use OCCT's point-on-surface projection
/// (`Extrema_ExtPS`), the same primitive `BRepExtrema_DistShapeShape` uses.
pub(super) fn face_point_distance(face: &crate::shape::Face, p: &GpPnt) -> f64 {
    let Some(surf) = crate::brep_tool::BRepTool::face_surface(face) else {
        return f64::INFINITY;
    };
    if crate::brep_surface::classify_surface(surf.as_ref()) == crate::brep_surface::SurfaceKind::Plane {
        if let Some(pln) = crate::brep_surface::face_plane(face) {
            let n = GpVec::from_xyz(pln.axis().direction().xyz()).normalized();
            return GpVec::from_pnts(&pln.location(), p).dot(&n).abs();
        }
    }
    occt_geom::geom_api::project_point_on_surface(
        surf.as_ref(),
        p,
        occt_core::precision::CONFUSION,
    )
    .map_or(f64::INFINITY, |ps| ps.distance)
}

/// Enclosed volume of a solid — `BRepGProp::VolumeProperties(S, Props)`
/// (`BRepGProp.cxx:555-586`), i.e. the divergence-theorem surface integral of
/// `BRepGProp_Vinert` over `BRepGProp_Gauss::Vinert` (`BRepGProp_Vinert.cxx:
/// 197-204`). That is the same routine OCCT's feature code uses, and it is
/// exact for the planar faces a feature boolean produces: the previous
/// hand-rolled per-face boundary-polygon fan was only valid for polyhedra and
/// under-counted every curved face (a 24-slice cylinder boss fused onto a box
/// measured 4.000000 instead of 4.155291).
pub(super) fn solid_volume(shape: &TopoShape) -> f64 {
    // `GProp_GProps` starts empty and accumulates whatever faces the Gauss
    // integrator could build; a face the port cannot represent is the
    // `BRepGProp_Gauss::Compute` failure OCCT swallows (`BRepGProp_Vinert.cxx:
    // 192`, the returned error is discarded).
    crate::brep_gprop_full::volume_properties(shape)
        .map(|p| p.mass().abs())
        .unwrap_or(0.0)
}
