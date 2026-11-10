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
        let volume = crate::brep_gprop::volume(&copy, tol.max(1e-4));
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
    let volume = crate::brep_gprop::volume(&copy, tol.max(1e-4));
    Ok(FeatResult { shape: copy, volume })
}

/// How far a `BRepFeat_MakeRevol` feature travels: one variant per `Perform`
/// overload of the reference (`hxx:54-102`).
///
/// The five arms share the class state and only differ in the sweep angle and
/// in whether a limit shape of the basis stops the sweep, so the feature layer
/// takes the limit instead of exposing five near-identical entry points.
#[derive(Debug, Clone)]
pub enum RevolLimit {
    /// `Perform(Angle)` (`BRepFeat_MakeRevol.cxx:200-336`).
    Angle(f64),
    /// `PerformThruAll()` (`cxx:655-667`), i.e. `Perform(2*M_PI)`.
    ThruAll,
    /// `Perform(Until)` (`cxx:340-467`): sweep an almost complete turn and stop
    /// the result at `Until`, a shape of the basis solid.
    Until(TopoShape),
    /// `Perform(From, Until)` (`cxx:473-648`): keep the part of the revolution
    /// between the two limit shapes of the basis.
    FromUntil(TopoShape, TopoShape),
    /// `PerformUntilAngle(Until, Angle)` (`cxx:670-792`): sweep `Angle` and stop
    /// the result at `Until`.
    UntilAngle(TopoShape, f64),
}

/// Groove: `BRepFeat_MakeRevol` in `Mode == 0` (cut) over the 2D `profile`
/// (`x` = radius >= 0, `y` = height along the axis) revolved around `axis`,
/// as a full turn (`PerformThruAll()`, `cxx:655-667`).
///
/// `tol` is kept for the caller's signature only: the reference's Boolean arm
/// is a plain `BRepAlgoAPI_Cut` with no fuzzy value (`cxx:308-327`), so the
/// port passes none either. See [`revolve_feature`].
pub fn groove(solid: &Solid, profile: &[GpPnt2d], axis: &GpAx1, tol: f64) -> Result<FeatResult, String> {
    groove_limited(solid, profile, axis, &RevolLimit::ThruAll, tol)
}

/// Neck: `BRepFeat_MakeRevol` in `Mode == 1` (fuse) over the same profile, as a
/// full turn (`PerformThruAll()`, `cxx:655-667`); `tol` is unused, see
/// [`groove`].
pub fn neck(solid: &Solid, profile: &[GpPnt2d], axis: &GpAx1, tol: f64) -> Result<FeatResult, String> {
    neck_limited(solid, profile, axis, &RevolLimit::ThruAll, tol)
}

/// [`groove`] with an explicit [`RevolLimit`]: the same cut, run through the
/// `Perform` overload the limit names (`BRepFeat_MakeRevol.cxx:200-792`).
///
/// The `Until` / `FromUntil` / `UntilAngle` arms take real shapes of the basis
/// solid, exactly as the reference does: the sweep is cut by the tool solid
/// built from the limit face (`BRepFeat::Tool`, `BRepFeat.cxx:643-712`) and the
/// piece that still carries the generating face is kept (`cxx:438-447`,
/// `:756-770`).
pub fn groove_limited(
    solid: &Solid,
    profile: &[GpPnt2d],
    axis: &GpAx1,
    limit: &RevolLimit,
    _tol: f64,
) -> Result<FeatResult, String> {
    require_min_profile(profile, "groove")?;
    revolve_feature(solid, profile, axis, false, limit, "groove")
}

/// [`neck`] with an explicit [`RevolLimit`]; `tol` is unused, see [`groove`].
pub fn neck_limited(
    solid: &Solid,
    profile: &[GpPnt2d],
    axis: &GpAx1,
    limit: &RevolLimit,
    _tol: f64,
) -> Result<FeatResult, String> {
    require_min_profile(profile, "neck")?;
    revolve_feature(solid, profile, axis, true, limit, "neck")
}

/// The planar profile face of a `groove`/`neck` cross-section `profile`
/// (`x` = radius, `y` = height along the axis), built in the XZ half-plane so
/// its normal is `-Y` and the revolution about `+Z` sweeps counter-clockwise.
///
/// This is `BRepFeat_MakeRevol`'s generating face `myPbase`: OCCT passes a
/// `TopoDS_Face` (a proper `TopoDS_Wire` of `TopoDS_Edge`s sharing their
/// `TopoDS_Vertex`es) to `BRepSweep_Revol::Perform` (`BRepFeat_MakeRevol.cxx:235`).
/// The shared vertices are load-bearing, not cosmetic: `BRepSweep_Tool` maps the
/// generating shape's sub-shapes with `TopExp::MapShapes`, so two vertices that
/// are merely coincident (distinct `TShape`s) are two entries and each spawns its
/// own directing edge -- the two caps then no longer share their circles with
/// the lateral face, and `BOPAlgo_BOP` sees duplicated coincident boundaries.
fn profile_face(profile: &[GpPnt2d]) -> Result<TopoShape, String> {
    use occt_core::gp::GpAx3;
    let b = crate::builder::TopoBuilder::new();
    let pln = GpPln::new(
        GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, -1.0, 0.0).map_err(|e| e.to_string())?,
            &GpDir::new(1.0, 0.0, 0.0).map_err(|e| e.to_string())?,
        )
        .map_err(|e| e.to_string())?,
    );
    let surf: std::sync::Arc<dyn occt_geom::surface::Surface> =
        std::sync::Arc::new(occt_geom::GeomPlane::new(pln));
    // A closed contour may repeat its first point as its last one (the
    // "already closed" form the `brepfeat` baselines use): the wrap edge built
    // below *is* that segment, so a repeated point would only make it degenerate.
    let mut n = profile.len();
    if n > 1 && profile[n - 1].distance(&profile[0]) <= occt_core::precision::CONFUSION {
        n -= 1;
    }
    let profile = &profile[..n];
    let verts: Vec<crate::shape::Vertex> = profile
        .iter()
        .map(|p| b.make_vertex(GpPnt::new(p.x(), 0.0, p.y()), occt_core::precision::CONFUSION))
        .collect();
    let mut edges: Vec<crate::shape::Edge> = Vec::new();
    for i in 0..n {
        let p = profile[i];
        let q = profile[(i + 1) % n];
        edges.push(b.make_edge_segment_with_vertices(
            &GpPnt::new(p.x(), 0.0, p.y()),
            &GpPnt::new(q.x(), 0.0, q.y()),
            &verts[i],
            &verts[(i + 1) % n],
        ));
    }
    let wire = b.make_wire(&edges);
    let face = b.make_face(surf, &[wire]);
    for e in &edges {
        let _ = crate::boptools_2d::make_2d(e, &face);
    }
    Ok(face.0)
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
    let tool = translated_cylinder(center, radius, height, z0 - margin)?;
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
    let tool = translated_cylinder(center, radius, height, z1)?;
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

/// The generating face `myPbase` of a [`revolve_feature`] call, placed on
/// `axis`: [`profile_face`] builds it in the XZ half-plane (revolution about
/// `+Z` through the origin), and [`align_z_to_axis`] carries that frame onto
/// the feature axis.
///
/// History (P042/T-171), still load-bearing for the generating face: the
/// analytic solid used to be rejected -- the `box + neck` fuse dropped the tool
/// (4.000000, expected 4.489835) and the `cylinder - groove` cut *added*
/// material. P039/T-169 aligned the sweep's wire/seam/pcurve conventions
/// (`BRepPrim_OneAxis.cxx:388-439`, `:666-679`; `BRepPrim_Builder.cxx:129-132`)
/// and P041/T-170 aligned `map_vertex_edges` to
/// `TopExp::MapShapesAndAncestors` occurrence counting (`TopExp.cxx:80-120`;
/// `BOPAlgo_BuilderFace.cxx:211-225`); both were necessary but not sufficient.
/// The remaining difference was the generating face itself: a profile whose
/// consecutive edges carry *coincident but distinct* vertices makes
/// `BRepSweep_Tool`'s `TopExp::MapShapes` see duplicate vertices, so each spawns
/// its own directing edge and the caps stop sharing their circles with the
/// lateral face -- `BOPAlgo_PaveFiller` then paves the coplanar target face
/// twice. With `profile_face` sharing the vertices, the tool of a `neck` is the
/// 3-edge/2-face topology `BRepPrimCylinder` has (occurrences=6, unique=3) and
/// every `brepfeat` baseline is exact (`box + neck` = 4.502654825 = 4 + pi*0.4^2).
fn profile_face_on(profile: &[GpPnt2d], axis: &GpAx1) -> Result<TopoShape, String> {
    let face = profile_face(profile)?;
    let t = align_z_to_axis(axis)?;
    crate::shape_ops::transformed_copy(&face, &t)
}

/// [`groove`] / [`neck`]: `BRepFeat_MakeRevol`'s `Perform` overloads
/// (`BRepFeat_MakeRevol.cxx:200-336`, `:340-467`, `:473-648`, `:655-667`,
/// `:670-792`).
///
/// `profile_face_on` is the generating face `myPbase` (`cxx:239-243`), swept by
/// `LocOpe_Revol` (`cxx:241-242`), recorded in `myMap` by `MajMap`
/// (`cxx:245-246`) and stored as `myGShape` (`cxx:248`).
///
/// `Skface` is a write-only slot on this path, so it is filled with `myPbase`
/// rather than a fabricated face: `Perform(Angle)` nullifies `mySkface` for a
/// complete turn (`cxx:236-238`), which skips both the `cxx:258-301` traversal
/// and `VerifGluedFaces` (`cxx:303-307`); with no `Add(E, F)` call `myGluedF`
/// stays empty, so the answer is the plain `BRepAlgoAPI_Fuse|Cut` of
/// `cxx:308-327`. (`TopoDS_Face()` would be the reference's value; this tree's
/// `TopoShape::is_null` is constantly `false` — `shape.rs:69-73` — so an
/// `Option` would be the equivalent and the class signature takes `&Face`.)
///
/// The `Until` arms consult it when they test `Until.IsSame(mySkface)`
/// (`cxx:359`, `:690`) and when `transform_shape_fu` / `GlobalPerform` read the
/// sketch-face flag; the generating face is what the reference's `Init` would
/// have received for a feature whose profile lies on the basis.
///
/// `Modify` (`cxx:109-117`) only reaches `GlobalPerform`, which the `Until` arms
/// do run, so the reference's `Standard_False` default is passed: the port has
/// no `Modify`-toggle above this layer either.
fn revolve_feature(
    solid: &Solid,
    profile: &[GpPnt2d],
    axis: &GpAx1,
    fuse: bool,
    limit: &RevolLimit,
    feature: &str,
) -> Result<FeatResult, String> {
    let pbase = profile_face_on(profile, axis).map_err(|e| format!("{feature}: {e}"))?;
    let skface = crate::shape::Face::wrap(pbase.clone())
        .ok_or_else(|| format!("{feature}: the generating face is not a face"))?;
    let mut feat = crate::brep_feat::MakeRevol::with_shapes(
        &solid.0,
        &pbase,
        &skface,
        axis,
        if fuse { 1 } else { 0 },
        false,
    );
    // `cxx:200-792`: the five `Perform` overloads, one arm per variant.
    match limit {
        RevolLimit::Angle(angle) => feat.perform_angle(*angle),
        RevolLimit::ThruAll => feat.perform_thru_all(),
        RevolLimit::Until(until) => feat.perform_until(until),
        RevolLimit::FromUntil(from, until) => feat.perform_from_until(from, until),
        RevolLimit::UntilAngle(until, angle) => feat.perform_until_angle(until, *angle),
    }
    .map_err(|e| format!("{feature}: {e}"))?;
    let shape = feat
        .form
        .shape()
        .cloned()
        .ok_or_else(|| format!("{feature}: NotDone (no myShape)"))?;
    // `BRepGProp::VolumeProperties` over the feature's own faces.
    let volume = solid_volume(&shape);
    Ok(FeatResult { shape, volume })
}

/// Analytic cylinder of `radius` × `height` (axis +Z) translated so its base
/// sits at `z_base` and its axis passes through `center`.
/// Source: `BRepPrimAPI_MakeCylinder`.
pub(super) fn translated_cylinder(
    center: &GpPnt,
    radius: f64,
    height: f64,
    z_base: f64,
) -> Result<TopoShape, String> {
    let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(radius, height);
    crate::shape_ops::translated_copy(&cyl.solid.0, &GpVec::new(center.x(), center.y(), z_base))
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
