use super::prelude::*;


/// Result of a feature operation: the modified solid and its volume.
#[derive(Debug, Clone)]

pub struct FeatResult {
    pub shape: TopoShape,
    pub volume: f64,
}

/// Draft one planar face of `solid` through `BRepOffsetAPI_DraftAngle`
/// (`BRepOffsetAPI_DraftAngle.cxx:91-105` `Add`, `:254-268` `Build`).
///
/// The face is located by `face_pts[0]` (nearest planar face). `angle_deg = 0`
/// returns a copy unchanged (and matches `Add`'s `|Angle| <= 1e-4` no-op).
/// `NeutralPlane` is `pivot`; `Direction` is the neutral-plane normal (the
/// mold pull of `Draft_Modification::Add` / `FindRotation`, `_1.cxx:2403-2453`)
/// — not the face normal, which collapses `FindRotation` toward a bad tilt.
pub fn draft(
    solid: &Solid,
    face_pts: &[GpPnt],
    angle_deg: f64,
    pivot: &GpPln,
    tol: f64,
) -> Result<FeatResult, String> {
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
    let solid_copy =
        Solid::wrap(copy.clone()).ok_or_else(|| "draft: result is not a solid".to_string())?;
    let face = face_containing_point(&solid_copy, &p0, tol).map_err(|e| format!("draft: {e}"))?;
    let _ = crate::brep_surface::face_plane(&face)
        .ok_or_else(|| "draft: the target face is not planar".to_string())?;
    // `Direction` of `Add` = pull = normal of the neutral plane.
    let direction = *pivot.axis().direction();

    let mut da = crate::brep_offset_api_draft_angle::BRepOffsetApiDraftAngle::with_shape(&copy);
    // `Flag = true` (`hxx` default of `Add`).
    da.add(&face, &direction, angle, pivot, true)
        .map_err(|e| format!("draft: {e}"))?;
    if !da
        .add_done()
        .map_err(|e| format!("draft: {e}"))?
    {
        let st = da
            .status()
            .map(|s| format!("{s:?}"))
            .unwrap_or_else(|_| "unknown".into());
        return Err(format!("draft: DraftAngle::Add not done ({st})"));
    }
    da.build().map_err(|e| format!("draft: {e}"))?;
    if !da.is_done() {
        let st = da
            .status()
            .map(|s| format!("{s:?}"))
            .unwrap_or_else(|_| "unknown".into());
        return Err(format!("draft: DraftAngle::Build not done ({st})"));
    }
    let shape = da
        .shape()
        .cloned()
        .ok_or_else(|| "draft: DraftAngle returned no shape".to_string())?;
    let volume = crate::brep_gprop::volume(&shape, tol.max(1e-4));
    Ok(FeatResult { shape, volume })
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

/// Through hole: `BRepFeat_MakeCylindricalHole::Perform(Radius)` followed by
/// `Build()` (`BRepFeat_MakeCylindricalHole.cxx:56-100`, `:630-648`) - a
/// `BOPAlgo_BOP` cut with the cylinder of `radius` about `axis`, made long
/// enough to cross the whole solid (`cxx:76-81`: an infinite cylinder cannot be
/// used for a topological operation, so the tool spans twice the bounding-box
/// extent along the axis).
///
/// `WithControl` is the reference's `true` default (`hxx:66-68`), so `Build`
/// runs `Validate` (`cxx:652-688`) and a result that still carries one of the
/// tool's caps is reported as an error instead of being accepted.
pub fn through_hole(solid: &Solid, axis: &GpAx1, radius: f64) -> Result<FeatResult, String> {
    let mut hole = crate::brep_feat::MakeCylindricalHole::new();
    hole.init(&solid.0, axis);
    hole.perform(radius)
        .map_err(|e| format!("through_hole: {e}"))?;
    finish_hole(hole, "through_hole")
}

/// Blind hole: `BRepFeat_MakeCylindricalHole::PerformBlind(Radius, Length, Cont)`
/// plus `Build()` (`BRepFeat_MakeCylindricalHole.cxx:474-626`, `:630-648`).
///
/// `length` is measured **from the origin of `axis`**, exactly as the reference
/// documents it (`hxx:80-87`) and enforces it (`cxx:516-528`: the next hit's
/// parameter must exceed `length`, otherwise `BRepFeat_HoleTooLong`). The tool
/// is the finite cylinder of `cxx:539-556`, which reaches from below the
/// bounding box up to `length`.
pub fn blind_hole(
    solid: &Solid,
    axis: &GpAx1,
    radius: f64,
    length: f64,
) -> Result<FeatResult, String> {
    let mut hole = crate::brep_feat::MakeCylindricalHole::new();
    hole.init(&solid.0, axis);
    hole.perform_blind(radius, length, true)
        .map_err(|e| format!("blind_hole: {e}"))?;
    finish_hole(hole, "blind_hole")
}

/// `Build()` (`cxx:630-648`) for the two hole entry points, then
/// [`solid_volume`] over the result.
///
/// A non-`NoError` status is the reference's "`myShape` is never assigned"
/// (`cxx:636-643`), so it is answered as an error carrying the status name.
fn finish_hole(
    mut hole: crate::brep_feat::MakeCylindricalHole,
    feature: &str,
) -> Result<FeatResult, String> {
    hole.build();
    let status = hole.status();
    if status != crate::brep_feat::status::Status::NoError {
        return Err(format!("{feature}: {status:?}"));
    }
    let shape = hole.shape().clone();
    let volume = solid_volume(&shape);
    Ok(FeatResult { shape, volume })
}

/// Rib: `BRepFeat_MakeLinearForm` (`Init` `cxx:84-857` + `Perform` `cxx:905-1028`).
///
/// Product params are the OCCT `Init` arguments (T-511 note (2) / T-536 note (2)):
/// - `wire`: sketch wire on `plane`; its endpoints must meet extreme faces of
///   `solid` (`ExtremeFaces`, `BRepFeat_RibSlot.cxx`);
/// - `plane`: sketch plane of the wire;
/// - `direc` / `direc1`: thickness vectors for the prism of the planar profile
///   (`Perform` sweeps `myPbase` along `|direc|+|direc1|`). They must leave the
///   sketch plane (OCCT demo uses the plane normal); a null `direc1` centres
///   the rib (`cxx:181-196`);
/// - `sliding`: `Modify` — try sliding faces on planar / cylindrical supports
///   (`cxx:280-353`). Fuse mode is fixed (`Mode = 1`).
pub fn rib(
    solid: &Solid,
    wire: &crate::shape::Wire,
    plane: &GpPln,
    direc: &GpVec,
    direc1: &GpVec,
    sliding: bool,
) -> Result<FeatResult, String> {
    let nulldir = GpVec::zero();
    let tol = occt_core::precision::CONFUSION;
    if direc.is_equal_tol(&nulldir, tol, tol) && direc1.is_equal_tol(&nulldir, tol, tol) {
        return Err("rib: Direc and Direc1 are both null".into());
    }
    let mut feat = crate::brep_feat::MakeLinearForm::with_shapes(
        &solid.0, wire, plane, direc, direc1, 1, sliding,
    )
    .map_err(|e| format!("rib: {e}"))?;
    if !feat.is_done() {
        return Err(format!("rib: {:?}", feat.current_status_error()));
    }
    feat.perform().map_err(|e| format!("rib: {e}"))?;
    if !feat.is_done() {
        return Err(format!("rib: {:?}", feat.current_status_error()));
    }
    let shape = feat
        .shape()
        .cloned()
        .ok_or_else(|| "rib: NotDone (no myShape)".to_string())?;
    let volume = solid_volume(&shape);
    Ok(FeatResult { shape, volume })
}

/// Boss through all: `BRepFeat_MakePrism::PerformThruAll` (`cxx:849-897`) over
/// a circular generating face on the solid's mid-height plane, fused into the
/// basis. The prism spans the bounding box, so the boss pierces the solid.
pub fn boss_thru_all(solid: &Solid, center: &GpPnt, radius: f64, _tol: f64) -> Result<FeatResult, String> {
    require_positive(radius, "radius", "boss_thru_all")?;
    let (z0, z1) = solid_z_span(&solid.0).map_err(|e| format!("boss_thru_all: {e}"))?;
    let z_mid = 0.5 * (z0 + z1);
    let pbase = disk_face_at(center, z_mid, radius).map_err(|e| format!("boss_thru_all: {e}"))?;
    prism_feature(solid, &pbase, true, PrismLimit::ThruAll, "boss_thru_all")
}

/// Boss: `BRepFeat_MakePrism::Perform(Length)` (`cxx:210-319`) over a circular
/// generating face on the solid's top plane. The cylinder's axis is +Z through
/// `center`; the base sits at the top z so the fuse adds material outside the
/// solid instead of swallowing an interior tool.
pub fn boss(solid: &Solid, center: &GpPnt, radius: f64, height: f64, _tol: f64) -> Result<FeatResult, String> {
    require_positive(radius, "radius", "boss")?;
    require_positive(height, "height", "boss")?;
    let (_, z1) = solid_z_span(&solid.0).map_err(|e| format!("boss: {e}"))?;
    let pbase = disk_face_at(center, z1, radius).map_err(|e| format!("boss: {e}"))?;
    prism_feature(solid, &pbase, true, PrismLimit::Length(height), "boss")
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
pub fn feature_apply(solid: &Solid, tool: &TopoShape, op: BoolOp, _tol: f64) -> Result<FeatResult, String> {
    boolean_feature(solid, tool, op)
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

/// How far a [`prism_feature`] boss travels: one variant per the
/// `BRepFeat_MakePrism` overload used by the product entry points.
enum PrismLimit {
    /// `Perform(Length)` (`BRepFeat_MakePrism.cxx:210-319`).
    Length(f64),
    /// `PerformThruAll()` (`cxx:849-897`).
    ThruAll,
}

/// Circular generating face of a boss, built as the top cap of
/// `BRepPrimAPI_MakeCylinder` (shared seam vertex, cap pcurve) and translated
/// so its plane sits at `z`. Reusing the primitive avoids a hand-rolled disk
/// that would re-open the FirstShape / WireSplitter issues T-511 hit on Revol.
fn disk_face_at(center: &GpPnt, z: f64, radius: f64) -> Result<TopoShape, String> {
    let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(radius, 1.0);
    let top = faces_of(&cyl.solid.0)
        .into_iter()
        .max_by(|a, b| {
            face_avg_z(a)
                .partial_cmp(&face_avg_z(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .ok_or_else(|| "disk_face_at: cylinder has no faces".to_string())?;
    let dz = z - face_avg_z(&top);
    crate::shape_ops::translated_copy(
        &top.0,
        &GpVec::new(center.x(), center.y(), dz),
    )
}

fn face_avg_z(face: &crate::shape::Face) -> f64 {
    let vs = crate::topo_tools::vertices_of(&face.0);
    if vs.is_empty() {
        return 0.0;
    }
    vs.iter()
        .map(|v| crate::topo_tools::vertex_position(v).z())
        .sum::<f64>()
        / vs.len() as f64
}

/// [`boss`] / [`boss_thru_all`]: `BRepFeat_MakePrism` over a circular
/// `myPbase`, fuse mode, `Modify = false` (the length / thru-all arms that do
/// not hit `GlobalPerform` still keep the Boolean history via
/// [`crate::brep_feat::algo_api`]; thru-all with fuse does call
/// `GlobalPerform`).
fn prism_feature(
    solid: &Solid,
    pbase: &TopoShape,
    fuse: bool,
    limit: PrismLimit,
    feature: &str,
) -> Result<FeatResult, String> {
    let skface = crate::shape::Face::wrap(pbase.clone())
        .ok_or_else(|| format!("{feature}: the generating face is not a face"))?;
    let dir = GpDir::new(0.0, 0.0, 1.0).map_err(|e| format!("{feature}: {e}"))?;
    let mut feat = crate::brep_feat::MakePrism::with_shapes(
        &solid.0,
        pbase,
        &skface,
        &dir,
        if fuse { 1 } else { 0 },
        false,
    );
    match limit {
        PrismLimit::Length(length) => feat.perform_length(length),
        PrismLimit::ThruAll => feat.perform_thru_all(),
    }
    .map_err(|e| format!("{feature}: {e}"))?;
    let shape = feat
        .form
        .shape()
        .cloned()
        .ok_or_else(|| format!("{feature}: NotDone (no myShape)"))?;
    let volume = solid_volume(&shape);
    Ok(FeatResult { shape, volume })
}

/// Combine `solid` with the `tool` solid through the history-preserving
/// `BRepAlgoAPI_*` facade (`brep_feat::algo_api`), the same door
/// `BRepFeat_Form::UpdateDescendants` uses after `GlobalPerform`.
pub(super) fn boolean_feature(solid: &Solid, tool: &TopoShape, op: BoolOp) -> Result<FeatResult, String> {
    let mut bop = match op {
        BoolOp::Fuse => crate::brep_feat::algo_api::fuse(&solid.0, tool),
        BoolOp::Cut => crate::brep_feat::algo_api::cut(&solid.0, tool),
        BoolOp::Common => crate::brep_feat::algo_api::common(&solid.0, tool),
    };
    bop.build()
        .map_err(|e| format!("feature boolean: {e}"))?;
    let shape = bop
        .shape()
        .cloned()
        .ok_or_else(|| "feature boolean: NotDone (no Shape)".to_string())?;
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
