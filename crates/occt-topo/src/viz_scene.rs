//! Phase 6 module: viz_scene — scene graph, camera, and raster pipeline.
//!
//! A lightweight port of OCCT's `AIS_Shape` (scene item), `V3d_View`
//! (camera) and `V3d_Viewer` (projection) for offline rendering. A
//! [`VizScene`] holds [`SceneShape`]s — each a `TopoShape` plus a world
//! transform and an optional color. A [`Camera`] provides a right-handed
//! look-at view with perspective or orthographic projection; world points are
//! mapped to screen pixels by [`project_point`].
//!
//! Two renderers are provided:
//! - [`render_scene_svg`] — triangles are projected, depth-sorted
//!   (painter's algorithm) and emitted as filled `<polygon>`s; a wireframe
//!   variant ([`render_scene_svg_wireframe`]) emits the projected edges.
//! - [`render_scene_ppm`] / [`render_scene_raster`] — a ray is cast per
//!   pixel through a BVH over the scene's triangles; hits are shaded with a
//!   Lambert term against a front light and written as binary PPM.
//! - [`render_scene_ppm_shaded`] / [`render_scene_raster_zbuffer`] — a
//!   hardware-style software renderer: per-shape [`Material`]s, point
//!   [`Light`]s, an ambient term and flat / Gouraud / Phong shading, produced
//!   either by per-pixel ray casting (the nearest BVH hit acts as the depth
//!   buffer) or by a screen-space z-buffer triangle rasterizer.
//!
//! Camera interaction ports a subset of `V3d_View`: orbit, pan, zoom, fit-all
//! and a world-space picking ray ([`Camera::camera_orbit`],
//! [`Camera::camera_pan`], [`Camera::camera_zoom`], [`Camera::camera_fit`],
//! [`Camera::camera_ray`]).

use std::cmp::Ordering;

use occt_core::bnd::BndBox;
use occt_core::bvh::bvh_ops::bvh_ray_cast;
use occt_core::bvh::builder_tri::build_tri_bvh;
use occt_core::gp::{GpMat, GpPnt, GpTrsf, GpVec, GpXyz};
use occt_core::poly::triangulation::Triangle;

use crate::mesh::ShapeMesh;
use crate::shape::TopoShape;

// ---------------------------------------------------------------------------
// Materials and lighting
// ---------------------------------------------------------------------------

/// A surface material: diffuse/specular/emissive reflectance, a shininess
/// exponent and an opacity.
///
/// The channels are normalized RGB in `[0, 1]`. Mirrors the important aspects
/// of OCCT's `Graphic3d_MaterialAspect` (DIFFUSE, SPECULAR, EMISSIVE,
/// SHININESS and TRANSPARENCY).
#[derive(Debug, Clone, Copy)]
pub struct Material {
    /// Diffuse reflectance — the base surface color.
    pub diffuse: (f64, f64, f64),
    /// Specular reflectance — the color of the highlight.
    pub specular: (f64, f64, f64),
    /// Emissive self-illumination (added unconditionally).
    pub emissive: (f64, f64, f64),
    /// Phong specular exponent; higher values give a tighter highlight.
    pub shininess: f64,
    /// Opacity in `[0, 1]`; `1.0` is fully opaque.
    pub opacity: f64,
}

impl Default for Material {
    /// Light steel gray diffuse, white specular, black emissive, 32 shininess,
    /// fully opaque.
    fn default() -> Self {
        Self {
            diffuse: (0.72, 0.72, 0.78),
            specular: (1.0, 1.0, 1.0),
            emissive: (0.0, 0.0, 0.0),
            shininess: 32.0,
            opacity: 1.0,
        }
    }
}

impl Material {
    /// A material with a plain diffuse color and default everything else.
    pub fn from_diffuse(diffuse: (f64, f64, f64)) -> Self {
        Self { diffuse, ..Default::default() }
    }

    /// `true` when the material lets some light pass through (`opacity < 1`).
    pub fn is_transparent(&self) -> bool {
        self.opacity < 1.0
    }
}

/// A point light: a world-space position, an RGB color and an intensity.
///
/// The light direction at a surface point is `position − point`; intensity
/// scales both the diffuse and the specular contribution.
#[derive(Debug, Clone, Copy)]
pub struct Light {
    /// Position of the point light in world coordinates.
    pub position: GpPnt,
    /// Light color (normalized RGB).
    pub color: (f64, f64, f64),
    /// Radiant intensity multiplier.
    pub intensity: f64,
}

impl Default for Light {
    /// A warm white key light above, to the right of and in front of the
    /// origin.
    fn default() -> Self {
        Self {
            position: GpPnt::new(3.0, 4.0, 6.0),
            color: (1.0, 1.0, 1.0),
            intensity: 1.0,
        }
    }
}

/// How the color varies across a triangle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadingMode {
    /// One color per triangle, from the face normal.
    Flat,
    /// Vertex normals are shaded and the colors are interpolated across the
    /// triangle (per-vertex lighting).
    Gouraud,
    /// The vertex normals are interpolated and every pixel is shaded with its
    /// own interpolated normal (per-pixel lighting).
    Phong,
}

/// Global renderer settings shared by the shaded renderers.
///
/// Holds the light list, the background color, a scalar ambient term and the
/// [`ShadingMode`]. Mirrors the light / background / shading state of a
/// `V3d_Viewer`.
#[derive(Debug, Clone)]
pub struct RenderSettings {
    /// Active point lights.
    pub lights: Vec<Light>,
    /// Background color (normalized RGB).
    pub background: (f64, f64, f64),
    /// Ambient light level in `[0, 1]` added to every lit surface.
    pub ambient: f64,
    /// Shading model used by the rasterizer.
    pub shading: ShadingMode,
}

impl Default for RenderSettings {
    /// One key light, a dark blue-gray background, 15% ambient and Phong
    /// shading.
    fn default() -> Self {
        Self {
            lights: vec![Light::default()],
            background: (0.03, 0.03, 0.05),
            ambient: 0.15,
            shading: ShadingMode::Phong,
        }
    }
}

// ---------------------------------------------------------------------------
// Scene graph
// ---------------------------------------------------------------------------

/// One item in a [`VizScene`]: a shape, its world transform, an optional
/// per-shape color and a [`Material`].
///
/// Mirrors OCCT's `AIS_Shape` which attaches a `TopoDS_Shape` and a display
/// transform to the interactive context. The legacy [`SceneShape::color`]
/// field remains a shorthand for the diffuse channel: when set, it overrides
/// `material.diffuse` (see [`SceneShape::material`]).
#[derive(Debug, Clone)]
pub struct SceneShape {
    pub shape: TopoShape,
    pub transform: GpTrsf,
    pub color: Option<(f64, f64, f64)>,
    /// Surface material used by the shaded renderers.
    pub material: Material,
}

impl SceneShape {
    /// A shape at the identity transform with the default material and no
    /// explicit color.
    pub fn new(shape: TopoShape) -> Self {
        Self { shape, transform: GpTrsf::identity(), color: None, material: Material::default() }
    }

    /// Copy with an explicit color. The color also becomes the material's
    /// diffuse channel so both views of the appearance stay in sync.
    pub fn with_color(mut self, color: (f64, f64, f64)) -> Self {
        self.color = Some(color);
        self.material.diffuse = color;
        self
    }

    /// Copy with an explicit material. A legacy `color`, when set, still wins
    /// for the diffuse channel.
    pub fn with_material(mut self, material: Material) -> Self {
        self.material = material;
        if let Some(c) = self.color {
            self.material.diffuse = c;
        }
        self
    }

    /// Copy with a world transform.
    pub fn with_transform(mut self, transform: GpTrsf) -> Self {
        self.transform = transform;
        self
    }

    /// The effective material: the stored material with its diffuse channel
    /// replaced by `color` when an explicit color is set.
    pub fn material(&self) -> Material {
        let mut m = self.material;
        if let Some(c) = self.color {
            m.diffuse = c;
        }
        m
    }

    /// Mesh the shape and transform every vertex into world space.
    pub fn mesh_world(&self, deflection: f64) -> ShapeMesh {
        let mut mesh = crate::shape_mesh::mesh_shape(&self.shape, deflection);
        for v in &mut mesh.vertices {
            v.transform(&self.transform);
        }
        mesh
    }

    /// World-space mesh vertices of this shape.
    pub fn world_vertices(&self, deflection: f64) -> Vec<GpPnt> {
        self.mesh_world(deflection).vertices
    }
}

/// An ordered collection of [`SceneShape`]s, the analogue of OCCT's
/// `AIS_InteractiveContext` display list.
#[derive(Debug, Clone, Default)]
pub struct VizScene {
    pub shapes: Vec<SceneShape>,
}

impl VizScene {
    /// An empty scene.
    pub fn new() -> Self {
        Self::default()
    }

    /// A scene built from an existing list of shapes.
    pub fn from_shapes(shapes: Vec<SceneShape>) -> Self {
        Self { shapes }
    }

    /// Append a shape to the scene.
    pub fn add(&mut self, shape: SceneShape) {
        self.shapes.push(shape);
    }

    /// Remove every shape.
    pub fn clear(&mut self) {
        self.shapes.clear();
    }

    /// `true` when the scene holds no shapes.
    pub fn is_empty(&self) -> bool {
        self.shapes.is_empty()
    }

    /// Number of shapes in the scene.
    pub fn len(&self) -> usize {
        self.shapes.len()
    }

    /// Borrow the shape at `idx`.
    pub fn get(&self, idx: usize) -> Option<&SceneShape> {
        self.shapes.get(idx)
    }

    /// Iterator over the scene's shapes.
    pub fn iter(&self) -> std::slice::Iter<'_, SceneShape> {
        self.shapes.iter()
    }

    /// Remove and return the shape at `idx`, if any.
    pub fn remove(&mut self, idx: usize) -> Option<SceneShape> {
        if idx < self.shapes.len() {
            Some(self.shapes.remove(idx))
        } else {
            None
        }
    }

    /// Total number of triangles across all meshed shapes.
    pub fn triangle_count(&self, deflection: f64) -> usize {
        self.shapes
            .iter()
            .map(|ss| crate::shape_mesh::mesh_shape(&ss.shape, deflection).triangles.len())
            .sum()
    }

    /// Replace the material of the shape at `idx`; returns `false` when `idx`
    /// is out of range.
    pub fn set_material(&mut self, idx: usize, material: Material) -> bool {
        match self.shapes.get_mut(idx) {
            Some(ss) => {
                ss.material = material;
                true
            }
            None => false,
        }
    }

    /// The effective material of the shape at `idx`, or `None` when out of
    /// range.
    pub fn shape_material(&self, idx: usize) -> Option<Material> {
        self.shapes.get(idx).map(|ss| ss.material())
    }
}

// ---------------------------------------------------------------------------
// Camera
// ---------------------------------------------------------------------------

/// Projection type of a [`Camera`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CameraProjection {
    /// Perspective projection: points are divided by their view depth, so
    /// far geometry appears smaller.
    Perspective,
    /// Orthographic projection: parallel projection with no depth scaling.
    Orthographic,
}

impl CameraProjection {
    pub fn is_perspective(&self) -> bool {
        matches!(self, CameraProjection::Perspective)
    }

    pub fn is_orthographic(&self) -> bool {
        matches!(self, CameraProjection::Orthographic)
    }
}

/// A look-at camera: an eye position, a target, an up hint, a vertical field
/// of view and near/far clip planes.
///
/// The view frame is built from `forward = normalize(target − eye)`,
/// `right = normalize(forward × up)` and `up2 = right × forward`, a
/// right-handed basis in which `+z` points in front of the camera. Ports the
/// camera state of OCCT's `V3d_View`.
#[derive(Debug, Clone, Copy)]
pub struct Camera {
    /// Camera position in world coordinates.
    pub eye: GpPnt,
    /// Point the camera looks at (defines the forward axis).
    pub target: GpPnt,
    /// Up hint; need not be unit, must not be parallel to `target − eye`.
    pub up: GpVec,
    /// Vertical field of view in degrees.
    pub fov_deg: f64,
    /// Perspective or orthographic projection.
    pub projection: CameraProjection,
    /// Near clip distance along the forward axis.
    pub near: f64,
    /// Far clip distance along the forward axis (informational).
    pub far: f64,
}

impl Default for Camera {
    /// Eye at `(0, 0, 5)` looking at the origin with `+Y` up, 45° vertical
    /// FOV, perspective projection.
    fn default() -> Self {
        Self {
            eye: GpPnt::new(0.0, 0.0, 5.0),
            target: GpPnt::zero(),
            up: GpVec::new(0.0, 1.0, 0.0),
            fov_deg: 45.0,
            projection: CameraProjection::Perspective,
            near: 0.1,
            far: 100.0,
        }
    }
}

impl Camera {
    /// Build a perspective camera from an explicit eye/target/up.
    pub fn look_at(eye: GpPnt, target: GpPnt, up: GpVec) -> Self {
        Self {
            eye,
            target,
            up,
            fov_deg: 45.0,
            projection: CameraProjection::Perspective,
            near: 0.1,
            far: 100.0,
        }
    }

    /// Copy with a different vertical field of view.
    pub fn with_fov(mut self, fov_deg: f64) -> Self {
        self.fov_deg = fov_deg;
        self
    }

    /// Copy forced to orthographic projection.
    pub fn orthographic(mut self) -> Self {
        self.projection = CameraProjection::Orthographic;
        self
    }

    /// Copy forced to perspective projection.
    pub fn perspective(mut self) -> Self {
        self.projection = CameraProjection::Perspective;
        self
    }

    /// Unit forward axis (from the eye toward the target).
    pub fn forward(&self) -> GpVec {
        camera_basis(self).2
    }

    /// Unit right axis of the view frame.
    pub fn right(&self) -> GpVec {
        camera_basis(self).0
    }

    /// Unit up axis of the view frame (perpendicular to forward).
    pub fn up(&self) -> GpVec {
        camera_basis(self).1
    }

    /// The look-at view matrix (see [`camera_view_matrix`]).
    pub fn view_matrix(&self) -> GpMat {
        camera_view_matrix(self)
    }

    /// The projection scale matrix (see [`camera_projection_matrix`]).
    pub fn projection_matrix(&self) -> GpMat {
        camera_projection_matrix(self)
    }

    /// Convenience wrapper over [`project_point`].
    pub fn project(&self, p: GpPnt, width: usize, height: usize) -> Option<(f64, f64, f64)> {
        project_point(self, p, width, height)
    }

    /// Map a camera-space point `(x, y, z)` (with `+z` in front of the
    /// camera) back to world coordinates.
    pub fn camera_to_world(&self, cam_x: f64, cam_y: f64, cam_z: f64) -> GpPnt {
        let world = camera_view_matrix(self).transpose().multiplied(&GpXyz::new(cam_x, cam_y, cam_z));
        GpPnt::new(self.eye.x() + world.x, self.eye.y() + world.y, self.eye.z() + world.z)
    }

    /// Ray through the centre of pixel `(x, y)` (perspective: from the eye;
    /// orthographic: parallel to the forward axis).
    pub fn ray_for_pixel(&self, width: usize, height: usize, x: usize, y: usize) -> (GpPnt, GpVec) {
        let (nx, ny) = screen_to_ndc(width, height, x as f64 + 0.5, y as f64 + 0.5);
        self.ray_through_ndc(nx, ny)
    }

    /// Ray through a point given in NDC `[-1, 1]²` (y up).
    pub fn ray_through_ndc(&self, x_ndc: f64, y_ndc: f64) -> (GpPnt, GpVec) {
        match self.projection {
            CameraProjection::Perspective => {
                let f = perspective_scale(self);
                // A camera-space point one unit in front of the eye that maps
                // to (x_ndc, y_ndc), transformed back to a world direction.
                let cam_pt = GpXyz::new(x_ndc / f, y_ndc / f, 1.0);
                let world = camera_view_matrix(self).transpose().multiplied(&cam_pt).normalized();
                (self.eye, GpVec::from_xyz(&world))
            }
            CameraProjection::Orthographic => {
                let os = ortho_scale(self);
                let off = self
                    .right()
                    .multiplied_scalar(x_ndc * os)
                    .add(&self.up().multiplied_scalar(y_ndc * os));
                let origin = GpPnt::new(self.eye.x() + off.x(), self.eye.y() + off.y(), self.eye.z() + off.z());
                (origin, self.forward())
            }
        }
    }

    /// Orbit the eye around the target.
    ///
    /// Yaws by `dx_deg` around the world-up axis through the target, then
    /// pitches by `dy_deg` around the camera's right axis. The pitch is
    /// clamped so the eye never crosses the poles (which would flip the view
    /// upside down). The target itself does not move, so the framed object
    /// stays centered. Ports `V3d_View::Rotate` for the orbit case.
    pub fn camera_orbit(&mut self, dx_deg: f64, dy_deg: f64) {
        let d = self.eye.distance(&self.target);
        if d < 1e-12 {
            return;
        }
        // Yaw around the world-up axis through the target.
        if dx_deg != 0.0 {
            let up_world = GpVec::new(0.0, 1.0, 0.0);
            let offset = rotate_vec(&GpVec::from_pnts(&self.target, &self.eye), &up_world, dx_deg.to_radians());
            self.eye = GpPnt::new(
                self.target.x() + offset.x(),
                self.target.y() + offset.y(),
                self.target.z() + offset.z(),
            );
        }
        // Pitch around the camera-right axis through the target, clamped to
        // keep the up hint sane.
        if dy_deg != 0.0 {
            let (right, _up2, _fwd) = camera_basis(self);
            let v = GpVec::from_pnts(&self.target, &self.eye);
            let sin_cur = (v.y() / d).clamp(-1.0, 1.0);
            let cur = sin_cur.asin();
            let max_p = (0.98_f64).asin();
            let total = (cur + dy_deg.to_radians()).clamp(-max_p, max_p);
            let delta = total - cur;
            if delta.abs() > 1e-12 {
                // Positive pitch raises the eye; the right axis is +X for the
                // default view, so rotate by the negated angle.
                let offset = rotate_vec(&v, &right, -delta);
                self.eye = GpPnt::new(
                    self.target.x() + offset.x(),
                    self.target.y() + offset.y(),
                    self.target.z() + offset.z(),
                );
            }
        }
    }

    /// Pan the view: move the eye and the target together along the camera's
    /// right/up plane.
    ///
    /// `dx`/`dy` are in screen pixels at a nominal 600-px viewport height;
    /// positive `dx` moves the scene right on screen, positive `dy` moves it
    /// up. Ports `V3d_View::Pan`.
    pub fn camera_pan(&mut self, dx: f64, dy: f64) {
        let (right, up2, _fwd) = camera_basis(self);
        let dist = self.eye.distance(&self.target);
        // World units per screen pixel, matching the frustum half-height at
        // the target distance.
        let scale = dist * (self.fov_deg.to_radians() * 0.5).tan() * 2.0 / 600.0;
        let off = right
            .multiplied_scalar(-dx * scale)
            .add(&up2.multiplied_scalar(-dy * scale));
        self.eye = GpPnt::new(self.eye.x() + off.x(), self.eye.y() + off.y(), self.eye.z() + off.z());
        self.target = GpPnt::new(
            self.target.x() + off.x(),
            self.target.y() + off.y(),
            self.target.z() + off.z(),
        );
    }

    /// Zoom by moving the eye toward (or away from) the target.
    ///
    /// `factor > 1` zooms in: the eye→target distance is divided by `factor`,
    /// so `camera_zoom(2.0)` halves the distance. The target is fixed. Ports
    /// `V3d_View::SetZoom`.
    pub fn camera_zoom(&mut self, factor: f64) {
        let d = self.eye.distance(&self.target);
        if factor <= 1e-9 || d < 1e-12 {
            return;
        }
        let new_d = d / factor;
        let fwd = GpVec::from_pnts(&self.eye, &self.target).normalized();
        self.eye = GpPnt::new(
            self.target.x() - fwd.x() * new_d,
            self.target.y() - fwd.y() * new_d,
            self.target.z() - fwd.z() * new_d,
        );
    }

    /// Frame the whole scene: replace this camera with [`fit_all_camera`].
    ///
    /// Ports `V3d_View::FitAll`.
    pub fn camera_fit(&mut self, scene: &VizScene) {
        *self = fit_all_camera(scene);
    }

    /// A world-space picking ray through the centre of pixel `(x, y)`.
    ///
    /// Perspective cameras cast from the eye through the pixel; orthographic
    /// cameras cast parallel to the forward axis. Wraps
    /// [`Camera::ray_for_pixel`]. Ports `V3d_View::Convert` + `Probe`.
    pub fn camera_ray(&self, width: usize, height: usize, px: usize, py: usize) -> (GpPnt, GpVec) {
        self.ray_for_pixel(width, height, px, py)
    }
}

/// Right-handed look-at view matrix.
///
/// Returns a 3×3 matrix whose rows are the camera basis
/// `(right, up2, forward)` built from the look-at vectors. Applying it to a
/// world *direction* yields camera coordinates `(x, y, z)` with `+z` in front
/// of the camera; the translation to the eye is applied separately by
/// [`project_point`] (a 3×3 matrix cannot carry the translation).
pub fn camera_view_matrix(cam: &Camera) -> GpMat {
    let (right, up2, fwd) = camera_basis(cam);
    GpMat::new(
        right.x(),
        right.y(),
        right.z(),
        up2.x(),
        up2.y(),
        up2.z(),
        fwd.x(),
        fwd.y(),
        fwd.z(),
    )
}

/// Projection scale matrix.
///
/// For a camera-space point `(x, y, z)` the matrix maps:
/// - perspective: `(f·x, f·y, z)` where `f = 1/tan(fov/2)` — the perspective
///   divide by `z` is applied by [`project_point`];
/// - orthographic: `(x/s, y/s, z)` where `s` is the ortho viewport half-height
///   (`dist(eye, target)·tan(fov/2)`), matching the perspective frustum at the
///   target distance.
pub fn camera_projection_matrix(cam: &Camera) -> GpMat {
    match cam.projection {
        CameraProjection::Perspective => {
            let f = perspective_scale(cam);
            GpMat::new(f, 0.0, 0.0, 0.0, f, 0.0, 0.0, 0.0, 1.0)
        }
        CameraProjection::Orthographic => {
            let inv = 1.0 / ortho_scale(cam);
            GpMat::new(inv, 0.0, 0.0, 0.0, inv, 0.0, 0.0, 0.0, 1.0)
        }
    }
}

/// Project a world point to screen coordinates.
///
/// Returns `(sx, sy, depth)` where `(sx, sy)` are pixel coordinates with the
/// origin at the top-left and `depth` is the positive view depth (distance
/// along the camera forward axis). Returns `None` when the point is behind
/// the near plane.
pub fn project_point(cam: &Camera, p: GpPnt, width: usize, height: usize) -> Option<(f64, f64, f64)> {
    let view = camera_view_matrix(cam);
    let d = GpXyz::new(p.x() - cam.eye.x(), p.y() - cam.eye.y(), p.z() - cam.eye.z());
    let cam_x = view.row(1).dot(&d);
    let cam_y = view.row(2).dot(&d);
    let cam_z = view.row(3).dot(&d);
    if cam_z < cam.near {
        return None;
    }
    let proj = camera_projection_matrix(cam);
    let px = proj.m[0][0] * cam_x;
    let py = proj.m[1][1] * cam_y;
    let (x_ndc, y_ndc) = match cam.projection {
        CameraProjection::Perspective => {
            if cam_z <= 0.0 {
                return None;
            }
            (px / cam_z, py / cam_z)
        }
        CameraProjection::Orthographic => (px, py),
    };
    let sx = (x_ndc + 1.0) * 0.5 * width as f64;
    let sy = (1.0 - y_ndc) * 0.5 * height as f64;
    Some((sx, sy, cam_z))
}

/// Frame a camera so the whole scene fits on screen.
///
/// Positions a default-FOV perspective camera along `+Z` looking at the scene
/// centre, at a distance that encloses the scene's bounding sphere. Ports the
/// behaviour of `V3d_View::FitAll`.
pub fn fit_all_camera(scene: &VizScene) -> Camera {
    let (lo, hi) = scene_bounds(scene)
        .unwrap_or_else(|| (GpPnt::new(-1.0, -1.0, -1.0), GpPnt::new(1.0, 1.0, 1.0)));
    let center = GpPnt::new((lo.x() + hi.x()) * 0.5, (lo.y() + hi.y()) * 0.5, (lo.z() + hi.z()) * 0.5);
    let radius = center.distance(&hi).max(1e-6);
    let cam = Camera::default();
    let dist = radius * 1.3 / (cam.fov_deg.to_radians() * 0.5).tan();
    Camera::look_at(
        GpPnt::new(center.x(), center.y(), center.z() + dist),
        center,
        GpVec::new(0.0, 1.0, 0.0),
    )
}

// ---------------------------------------------------------------------------
// Raster
// ---------------------------------------------------------------------------

/// A `width`×`height` color raster.
///
/// `pixels` holds three normalized RGB components per pixel in `[0, 1]`,
/// row-major, so `pixels.len() == width * height * 3`. Convert to binary PPM
/// bytes with [`Raster::to_ppm_bytes`].
#[derive(Debug, Clone)]
pub struct Raster {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<f64>,
}

impl Raster {
    /// A raster filled with black.
    pub fn new(width: usize, height: usize) -> Self {
        Self { width, height, pixels: vec![0.0; width * height * 3] }
    }

    /// Set a pixel's color (components are clamped to `[0, 1]`).
    pub fn set_pixel(&mut self, x: usize, y: usize, rgb: (f64, f64, f64)) {
        let i = (y * self.width + x) * 3;
        if i + 2 < self.pixels.len() {
            self.pixels[i] = rgb.0.clamp(0.0, 1.0);
            self.pixels[i + 1] = rgb.1.clamp(0.0, 1.0);
            self.pixels[i + 2] = rgb.2.clamp(0.0, 1.0);
        }
    }

    /// Read a pixel's color.
    pub fn get_pixel(&self, x: usize, y: usize) -> (f64, f64, f64) {
        let i = (y * self.width + x) * 3;
        (self.pixels[i], self.pixels[i + 1], self.pixels[i + 2])
    }

    /// Fill every pixel with one color.
    pub fn clear(&mut self, rgb: (f64, f64, f64)) {
        for chunk in self.pixels.chunks_exact_mut(3) {
            chunk[0] = rgb.0.clamp(0.0, 1.0);
            chunk[1] = rgb.1.clamp(0.0, 1.0);
            chunk[2] = rgb.2.clamp(0.0, 1.0);
        }
    }

    /// Fill the raster by evaluating `f(x, y)` at every pixel.
    pub fn fill(&mut self, mut f: impl FnMut(usize, usize) -> (f64, f64, f64)) {
        for y in 0..self.height {
            for x in 0..self.width {
                self.set_pixel(x, y, f(x, y));
            }
        }
    }

    /// Binary PPM bytes: `P6\n<width> <height>\n255\n` followed by RGB bytes.
    pub fn to_ppm_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(13 + self.pixels.len());
        out.extend_from_slice(format!("P6\n{} {}\n255\n", self.width, self.height).as_bytes());
        for chunk in self.pixels.chunks_exact(3) {
            out.push((chunk[0].clamp(0.0, 1.0) * 255.0) as u8);
            out.push((chunk[1].clamp(0.0, 1.0) * 255.0) as u8);
            out.push((chunk[2].clamp(0.0, 1.0) * 255.0) as u8);
        }
        out
    }

    /// One `[r, g, b]` byte triple per pixel, row-major.
    pub fn as_rgb8(&self) -> Vec<[u8; 3]> {
        self.pixels
            .chunks_exact(3)
            .map(|c| {
                [
                    (c[0].clamp(0.0, 1.0) * 255.0) as u8,
                    (c[1].clamp(0.0, 1.0) * 255.0) as u8,
                    (c[2].clamp(0.0, 1.0) * 255.0) as u8,
                ]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Triangle soup
// ---------------------------------------------------------------------------

/// A world-space triangle with its shape's color. The basic primitive shared
/// by the SVG and PPM renderers.
#[derive(Debug, Clone, Copy)]
pub struct SceneTriangle {
    pub a: GpPnt,
    pub b: GpPnt,
    pub c: GpPnt,
    pub color: (f64, f64, f64),
}

impl SceneTriangle {
    pub fn new(a: GpPnt, b: GpPnt, c: GpPnt, color: (f64, f64, f64)) -> Self {
        Self { a, b, c, color }
    }

    /// Unit (unnormalized fallback) face normal via `(b − a) × (c − a)`.
    pub fn normal(&self) -> GpVec {
        GpVec::from_pnts(&self.a, &self.b).crossed(&GpVec::from_pnts(&self.a, &self.c)).normalized()
    }

    /// Centroid of the triangle.
    pub fn centroid(&self) -> GpPnt {
        GpPnt::new(
            (self.a.x() + self.b.x() + self.c.x()) / 3.0,
            (self.a.y() + self.b.y() + self.c.y()) / 3.0,
            (self.a.z() + self.b.z() + self.c.z()) / 3.0,
        )
    }

    /// Surface area of the triangle.
    pub fn area(&self) -> f64 {
        0.5 * GpVec::from_pnts(&self.a, &self.b).crossed(&GpVec::from_pnts(&self.a, &self.c)).magnitude()
    }

    /// The three edges `(a,b)`, `(b,c)`, `(c,a)`.
    pub fn edges(&self) -> [(GpPnt, GpPnt); 3] {
        [(self.a, self.b), (self.b, self.c), (self.c, self.a)]
    }

    /// `true` when the triangle is degenerate (zero area within `tol`).
    pub fn is_degenerate(&self, tol: f64) -> bool {
        self.area() <= tol
    }

    /// `true` when the face normal points along `view` (i.e. the triangle
    /// faces toward `view`).
    pub fn faces_toward(&self, view: &GpVec) -> bool {
        self.normal().dot(view) > 0.0
    }
}

/// Default color (light steel) for shapes without an explicit color.
fn default_color() -> (f64, f64, f64) {
    (0.72, 0.72, 0.78)
}

/// Gather every triangle of every scene shape, transformed to world space.
///
/// Each shape is tessellated at the given deflection, its vertices are pushed
/// through the shape transform, and every mesh triangle is emitted with the
/// shape's color (or the default).
pub fn scene_triangles(scene: &VizScene, deflection: f64) -> Vec<SceneTriangle> {
    let mut out = Vec::new();
    for ss in &scene.shapes {
        let mesh = crate::shape_mesh::mesh_shape(&ss.shape, deflection);
        let color = ss.color.unwrap_or_else(default_color);
        for t in &mesh.triangles {
            out.push(SceneTriangle::new(
                mesh.vertices[t.n0].transformed(&ss.transform),
                mesh.vertices[t.n1].transformed(&ss.transform),
                mesh.vertices[t.n2].transformed(&ss.transform),
                color,
            ));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Bounds and vertex queries
// ---------------------------------------------------------------------------

/// Combined axis-aligned bounds of every transformed shape in the scene.
///
/// Returns `(min_corner, max_corner)`, or `None` for an empty scene. The
/// per-shape box is taken from the registered geometry and then transformed
/// into world space.
pub fn scene_bounds(scene: &VizScene) -> Option<(GpPnt, GpPnt)> {
    let mut b = BndBox::new();
    for ss in &scene.shapes {
        let bb = crate::bbox_from_geometry::shape_bbox(&ss.shape).transformed(&ss.transform);
        b.add_box(&bb);
    }
    b.get()
        .map(|(x0, x1, y0, y1, z0, z1)| (GpPnt::new(x0, y0, z0), GpPnt::new(x1, y1, z1)))
}

/// World-space mesh vertices of the shape at `idx` in `scene`.
pub fn shape_world_vertices(scene: &VizScene, idx: usize, deflection: f64) -> Vec<GpPnt> {
    scene.shapes[idx].world_vertices(deflection)
}

/// Centre of the scene's combined bounds (falls back to the origin when the
/// scene is empty).
pub fn scene_center(scene: &VizScene) -> GpPnt {
    let (lo, hi) = scene_bounds(scene)
        .unwrap_or_else(|| (GpPnt::new(-1.0, -1.0, -1.0), GpPnt::new(1.0, 1.0, 1.0)));
    GpPnt::new((lo.x() + hi.x()) * 0.5, (lo.y() + hi.y()) * 0.5, (lo.z() + hi.z()) * 0.5)
}

/// Möller–Trumbore ray/triangle intersection over a [`SceneTriangle`].
///
/// Returns the ray parameter `t` at the hit, or `None` when the ray is
/// parallel to the triangle, misses it, or the hit lies at or behind the
/// origin. Wraps [`occt_core::bvh::bvh_ops::ray_triangle_t`].
pub fn ray_hits_triangle(origin: &GpPnt, dir: &GpVec, tri: &SceneTriangle) -> Option<f64> {
    occt_core::bvh::bvh_ops::ray_triangle_t(origin, dir, &tri.a, &tri.b, &tri.c)
}

// ---------------------------------------------------------------------------
// Rendering — SVG
// ---------------------------------------------------------------------------

/// A triangle projected to screen space, ready for painter sorting.
struct ProjectedTri {
    depth: f64,
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    color: (f64, f64, f64),
}

/// Project every scene triangle and return them sorted far-to-near.
///
/// Triangles behind the camera or with a degenerate (zero-area) projection
/// are dropped.
fn project_scene(scene: &VizScene, cam: &Camera, width: usize, height: usize, deflection: f64) -> Vec<ProjectedTri> {
    let mut out = Vec::new();
    for t in scene_triangles(scene, deflection) {
        let (Some(pa), Some(pb), Some(pc)) = (
            project_point(cam, t.a, width, height),
            project_point(cam, t.b, width, height),
            project_point(cam, t.c, width, height),
        ) else {
            continue;
        };
        let area = (pb.0 - pa.0) * (pc.1 - pa.1) - (pb.1 - pa.1) * (pc.0 - pa.0);
        if area.abs() < 1e-9 {
            continue;
        }
        out.push(ProjectedTri {
            depth: (pa.2 + pb.2 + pc.2) / 3.0,
            a: (pa.0, pa.1),
            b: (pb.0, pb.1),
            c: (pc.0, pc.1),
            color: t.color,
        });
    }
    out.sort_by(|x, y| y.depth.partial_cmp(&x.depth).unwrap_or(Ordering::Equal));
    out
}

/// Normalized RGB in `[0, 1]` → `#rrggbb`.
fn color_hex(c: (f64, f64, f64)) -> String {
    let r = (c.0.clamp(0.0, 1.0) * 255.0).round() as u8;
    let g = (c.1.clamp(0.0, 1.0) * 255.0).round() as u8;
    let b = (c.2.clamp(0.0, 1.0) * 255.0).round() as u8;
    format!("#{r:02x}{g:02x}{b:02x}")
}

/// Render the scene as a filled SVG document.
///
/// Each shape is meshed, projected through `cam`, and the visible triangles
/// are drawn as `<polygon>`s in painter order (far to near) so nearer
/// geometry overpaints farther geometry. The SVG uses a
/// `viewBox="0 0 <width> <height>"` canvas.
pub fn render_scene_svg(scene: &VizScene, cam: &Camera, width: usize, height: usize, deflection: f64) -> String {
    let drawn = project_scene(scene, cam, width, height, deflection);
    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width} {height}\">\n"
    ));
    for t in &drawn {
        out.push_str(&format!(
            "<polygon points=\"{:.3},{:.3} {:.3},{:.3} {:.3},{:.3}\" fill=\"{}\" stroke=\"black\" stroke-width=\"0.5\" />\n",
            t.a.0, t.a.1, t.b.0, t.b.1, t.c.0, t.c.1, color_hex(t.color)
        ));
    }
    out.push_str("</svg>\n");
    out
}

/// Render the scene as an SVG wireframe.
///
/// Every triangle edge that projects inside the viewport is emitted as a
/// `<polyline>`. This is a full (transparent) wireframe rather than an HLR
/// view; use [`crate::hlr::wireframe_projection`] for hidden-line removal.
pub fn render_scene_svg_wireframe(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
) -> String {
    let tris = scene_triangles(scene, deflection);
    let mut edges: Vec<((f64, f64), (f64, f64))> = Vec::new();
    for t in &tris {
        for (p, q) in t.edges() {
            if let (Some(pa), Some(pb)) = (project_point(cam, p, width, height), project_point(cam, q, width, height)) {
                edges.push(((pa.0, pa.1), (pb.0, pb.1)));
            }
        }
    }
    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width} {height}\">\n"
    ));
    for ((x1, y1), (x2, y2)) in edges {
        out.push_str(&format!(
            "<polyline points=\"{:.3},{:.3} {:.3},{:.3}\" fill=\"none\" stroke=\"black\" stroke-width=\"0.5\" />\n",
            x1, y1, x2, y2
        ));
    }
    out.push_str("</svg>\n");
    out
}

// ---------------------------------------------------------------------------
// Rendering — raster / PPM
// ---------------------------------------------------------------------------

/// Lambert shading of `tri` under a light direction.
///
/// `shade = 0.35 + 0.65 · max(0, n·l)` — an ambient term keeps back-angled
/// faces visible while `n·l` darkens faces turned away from the light.
fn shade_lambert(tri: &SceneTriangle, light: &GpVec) -> (f64, f64, f64) {
    let lambert = tri.normal().dot(light).max(0.0);
    let s = 0.35 + 0.65 * lambert;
    (tri.color.0 * s, tri.color.1 * s, tri.color.2 * s)
}

/// Rasterize the scene by casting one ray per pixel through a BVH.
///
/// The returned [`Raster`] is a `width`×`height` image whose pixels hold the
/// shaded color of the nearest triangle hit (front-lit Lambert shading), or a
/// dark background where no triangle is hit.
pub fn render_scene_raster(scene: &VizScene, cam: &Camera, width: usize, height: usize, deflection: f64) -> Raster {
    let tris = scene_triangles(scene, deflection);
    let tri_pts: Vec<(GpPnt, GpPnt, GpPnt)> = tris.iter().map(|t| (t.a, t.b, t.c)).collect();
    let bvh = build_tri_bvh(&tri_pts, 8);
    let light = GpVec::from_pnts(&cam.eye, &cam.target).normalized();
    let mut raster = Raster::new(width, height);
    raster.fill(|x, y| {
        let (nx, ny) = screen_to_ndc(width, height, x as f64 + 0.5, y as f64 + 0.5);
        let (origin, dir) = cam.ray_through_ndc(nx, ny);
        match bvh_ray_cast(&bvh, &tri_pts, &origin, &dir) {
            Some((idx, _t)) => shade_lambert(&tris[idx], &light),
            None => (0.03, 0.03, 0.05),
        }
    });
    raster
}

/// Rasterize the scene with a screen-space z-buffer (alternative to the
/// ray-cast [`render_scene_raster`]).
///
/// Triangles are projected, then every pixel inside each triangle's bounding
/// box is filled when the barycentrically interpolated view depth is nearer
/// than the current z-buffer value. Cheaper than ray casting for small scenes
/// but uses an affine depth approximation (fine for near-planar triangles).
pub fn render_scene_zbuffer_raster(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
) -> Raster {
    let tris = scene_triangles(scene, deflection);
    let light = GpVec::from_pnts(&cam.eye, &cam.target).normalized();
    let mut raster = Raster::new(width, height);
    raster.clear((0.03, 0.03, 0.05));
    let mut depth = vec![f64::INFINITY; width * height];
    for t in &tris {
        let (Some(pa), Some(pb), Some(pc)) = (
            project_point(cam, t.a, width, height),
            project_point(cam, t.b, width, height),
            project_point(cam, t.c, width, height),
        ) else {
            continue;
        };
        let area2 = (pb.0 - pa.0) * (pc.1 - pa.1) - (pb.1 - pa.1) * (pc.0 - pa.0);
        if area2.abs() < 1e-9 {
            continue; // edge-on triangle
        }
        let xmin = pa.0.min(pb.0).min(pc.0).max(0.0) as usize;
        let xmax = pa.0.max(pb.0).max(pc.0).min(width as f64 - 1.0) as usize;
        let ymin = pa.1.min(pb.1).min(pc.1).max(0.0) as usize;
        let ymax = pa.1.max(pb.1).max(pc.1).min(height as f64 - 1.0) as usize;
        let shade = shade_lambert(t, &light);
        for y in ymin..=ymax {
            for x in xmin..=xmax {
                let (px, py) = (x as f64 + 0.5, y as f64 + 0.5);
                let w0 = (pb.0 - pa.0) * (py - pa.1) - (pb.1 - pa.1) * (px - pa.0);
                let w1 = (pc.0 - pb.0) * (py - pb.1) - (pc.1 - pb.1) * (px - pb.0);
                let w2 = (pa.0 - pc.0) * (py - pc.1) - (pa.1 - pc.1) * (px - pc.0);
                let inside = (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0)
                    || (w0 <= 0.0 && w1 <= 0.0 && w2 <= 0.0);
                if !inside {
                    continue;
                }
                let d = (w0 * pa.2 + w1 * pb.2 + w2 * pc.2) / (w0 + w1 + w2);
                let i = y * width + x;
                if d < depth[i] {
                    depth[i] = d;
                    raster.set_pixel(x, y, shade);
                }
            }
        }
    }
    raster
}

/// Render the scene to binary PPM bytes.
///
/// Equivalent to rasterizing with [`render_scene_raster`] and serializing to
/// the P6 PPM format (`P6\n<width> <height>\n255\n` + RGB bytes). This is the
/// analogue of `V3d_View::Dump` writing a `Write_PPM` image.
pub fn render_scene_ppm(scene: &VizScene, cam: &Camera, width: usize, height: usize, deflection: f64) -> Vec<u8> {
    render_scene_raster(scene, cam, width, height, deflection).to_ppm_bytes()
}

/// Write PPM bytes to `path`.
pub fn write_ppm(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    std::fs::write(path, bytes)
}

// ---------------------------------------------------------------------------
// Rendering — shaded ray cast / PPM
// ---------------------------------------------------------------------------

/// The reflection of `light_dir` about `normal`, `r = 2(n·l)n − l`.
///
/// Both vectors are normalized internally; the returned vector is unit (or
/// zero when the vectors are degenerate). This is the `R` used by the Phong
/// specular term `(r·v)^shininess`.
pub fn phong_reflect(normal: GpVec, light_dir: GpVec) -> GpVec {
    let n = normal.normalized();
    let l = light_dir.normalized();
    n.multiplied_scalar(2.0 * n.dot(&l)).subtracted(&l).normalized()
}

/// The unit direction from `point` toward a light.
///
/// Returns a zero vector when the light coincides with `point`.
pub fn light_dir_at(light: &Light, point: GpPnt) -> GpVec {
    GpVec::from_pnts(&point, &light.position).normalized()
}

/// Phong-shade a surface point under one light.
///
/// `normal`, `view_dir` and `light_dir` need not be unit; each is normalized
/// internally. Returns the light's contribution as ambient + diffuse +
/// specular:
///
/// - ambient: `diffuse · ambient`;
/// - diffuse: `diffuse · light.color · intensity · max(0, n·l)`;
/// - specular: `specular · light.color · intensity · max(0, r·v)^shininess`
///   with `r = 2(n·l)n − l`, only when `n·l > 0` (a light behind the surface
///   produces no highlight).
///
/// The components are not clamped here; the rasterizers clamp on write.
pub fn shade_phong(
    material: &Material,
    light: &Light,
    normal: GpVec,
    view_dir: GpVec,
    light_dir: GpVec,
    ambient: f64,
) -> (f64, f64, f64) {
    let n = normal.normalized();
    let l = light_dir.normalized();
    let v = view_dir.normalized();
    let ndotl = n.dot(&l).max(0.0);
    let ambient_rgb = (
        material.diffuse.0 * ambient,
        material.diffuse.1 * ambient,
        material.diffuse.2 * ambient,
    );
    let diffuse_rgb = (
        material.diffuse.0 * light.color.0 * light.intensity * ndotl,
        material.diffuse.1 * light.color.1 * light.intensity * ndotl,
        material.diffuse.2 * light.color.2 * light.intensity * ndotl,
    );
    let spec_rgb = if ndotl > 0.0 {
        let r = phong_reflect(n, l);
        let rdotv = r.dot(&v).max(0.0).powf(material.shininess);
        (
            material.specular.0 * light.color.0 * light.intensity * rdotv,
            material.specular.1 * light.color.1 * light.intensity * rdotv,
            material.specular.2 * light.color.2 * light.intensity * rdotv,
        )
    } else {
        (0.0, 0.0, 0.0)
    };
    (
        ambient_rgb.0 + diffuse_rgb.0 + spec_rgb.0,
        ambient_rgb.1 + diffuse_rgb.1 + spec_rgb.1,
        ambient_rgb.2 + diffuse_rgb.2 + spec_rgb.2,
    )
}

/// Sum the contribution of every light in `settings` at a surface point, then
/// add the material's emissive term.
///
/// `hit` is the world-space position of the shaded point — it is used to
/// compute each light direction. Opacity is applied as a blend toward
/// `settings.background`.
fn shade_point(
    material: &Material,
    settings: &RenderSettings,
    normal: GpVec,
    view_dir: GpVec,
    hit: GpPnt,
) -> (f64, f64, f64) {
    let mut rgb = material.emissive;
    for light in &settings.lights {
        let c = shade_phong(
            material,
            light,
            normal,
            view_dir,
            light_dir_at(light, hit),
            settings.ambient,
        );
        rgb.0 += c.0;
        rgb.1 += c.1;
        rgb.2 += c.2;
    }
    if material.opacity < 1.0 {
        let a = material.opacity;
        rgb.0 = rgb.0 * a + settings.background.0 * (1.0 - a);
        rgb.1 = rgb.1 * a + settings.background.1 * (1.0 - a);
        rgb.2 = rgb.2 * a + settings.background.2 * (1.0 - a);
    }
    rgb
}

/// Mesh every shape of the scene into world-space [`ShadeMesh`]es.
///
/// The returned vector is parallel to `scene.shapes`, so shape index `i` maps
/// to `meshes[i]` and its effective material is `meshes[i].material`.
fn scene_shade_meshes(scene: &VizScene, deflection: f64) -> Vec<ShadeMesh> {
    scene
        .shapes
        .iter()
        .map(|ss| ShadeMesh::from_scene_shape(ss, deflection))
        .collect()
}

/// Flatten every mesh of `meshes` into one triangle buffer with `(shape,
/// triangle)` tags, ready for [`build_tri_bvh`] + [`bvh_ray_cast`].
fn flatten_shade_meshes(meshes: &[ShadeMesh]) -> (Vec<(GpPnt, GpPnt, GpPnt)>, Vec<(usize, usize)>) {
    let mut tri_pts = Vec::new();
    let mut tags = Vec::new();
    for (si, m) in meshes.iter().enumerate() {
        for (ti, &(n0, n1, n2)) in m.triangles.iter().enumerate() {
            tri_pts.push((m.verts[n0], m.verts[n1], m.verts[n2]));
            tags.push((si, ti));
        }
    }
    (tri_pts, tags)
}

/// Barycentric coordinates of `p` in triangle `(a, b, c)`, via signed areas.
///
/// Returns `(la, lb, lc)` with `la + lb + lc = 1`. A degenerate triangle falls
/// back to the centroid weights `(1/3, 1/3, 1/3)`.
fn triangle_barycentric(a: &GpPnt, b: &GpPnt, c: &GpPnt, p: &GpPnt) -> (f64, f64, f64) {
    let n = GpVec::from_pnts(a, b).crossed(&GpVec::from_pnts(a, c));
    let n2 = n.square_magnitude();
    if n2 < 1e-24 {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }
    let la = GpVec::from_pnts(b, c).crossed(&GpVec::from_pnts(b, p)).dot(&n) / n2;
    let lb = GpVec::from_pnts(c, a).crossed(&GpVec::from_pnts(c, p)).dot(&n) / n2;
    let lc = GpVec::from_pnts(a, b).crossed(&GpVec::from_pnts(a, p)).dot(&n) / n2;
    (la, lb, lc)
}

/// Ray-cast the scene into a [`Raster`] with per-shape materials.
///
/// A single [`TriBvh`] is built over every shape's world triangles; each
/// triangle keeps its owning shape and triangle index so the correct material
/// and vertex normals are used for shading. One ray is cast per pixel; the
/// nearest BVH hit acts as the depth buffer.
///
/// The shading mode is honored through the hit's world-space barycentric
/// coordinates:
///
/// - `Flat` — the face normal shades the hit;
/// - `Gouraud` / `Phong` — the per-vertex normals are interpolated at the hit
///   and the point is shaded with the interpolated normal.
///
/// Missed pixels get `settings.background`.
pub fn render_scene_raster_shaded(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    settings: &RenderSettings,
) -> Raster {
    let meshes = scene_shade_meshes(scene, deflection);
    let (tri_pts, tags) = flatten_shade_meshes(&meshes);
    let bvh = build_tri_bvh(&tri_pts, 8);
    let mut raster = Raster::new(width, height);
    raster.fill(|x, y| {
        let (nx, ny) = screen_to_ndc(width, height, x as f64 + 0.5, y as f64 + 0.5);
        let (origin, dir) = cam.ray_through_ndc(nx, ny);
        match bvh_ray_cast(&bvh, &tri_pts, &origin, &dir) {
            Some((idx, t)) => {
                let (si, ti) = tags[idx];
                let m = &meshes[si];
                let (n0, n1, n2) = m.triangles[ti];
                let (a, b, c) = (m.verts[n0], m.verts[n1], m.verts[n2]);
                let hit = GpPnt::new(
                    origin.x() + dir.x() * t,
                    origin.y() + dir.y() * t,
                    origin.z() + dir.z() * t,
                );
                let normal = match settings.shading {
                    ShadingMode::Flat => {
                        GpVec::from_pnts(&a, &b).crossed(&GpVec::from_pnts(&a, &c)).normalized()
                    }
                    ShadingMode::Gouraud | ShadingMode::Phong => {
                        let (la, lb, lc) = triangle_barycentric(&a, &b, &c, &hit);
                        m.normals[n0]
                            .multiplied_scalar(la)
                            .add(&m.normals[n1].multiplied_scalar(lb))
                            .add(&m.normals[n2].multiplied_scalar(lc))
                            .normalized()
                    }
                };
                let view_dir = GpVec::from_pnts(&hit, &cam.eye);
                let (r, g, b) = shade_point(&m.material, settings, normal, view_dir, hit);
                (r, g, b)
            }
            None => settings.background,
        }
    });
    raster
}

/// Pick the nearest shape under a screen pixel.
///
/// Casts a world-space ray through pixel `(px, py)` (see [`Camera::camera_ray`])
/// against a BVH over every shape's world triangles and returns the owning
/// shape index, the world-space hit point and the ray parameter `t`. This is
/// the picking analogue of `V3d_View::Pick` / `Select`; the returned shape
/// index can be fed back into [`VizScene::get`].
pub fn pick_shape(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    px: usize,
    py: usize,
    deflection: f64,
) -> Option<(usize, GpPnt, f64)> {
    let meshes = scene_shade_meshes(scene, deflection);
    let (tri_pts, tags) = flatten_shade_meshes(&meshes);
    let bvh = build_tri_bvh(&tri_pts, 8);
    let (origin, dir) = cam.camera_ray(width, height, px, py);
    bvh_ray_cast(&bvh, &tri_pts, &origin, &dir).map(|(idx, t)| {
        let hit = GpPnt::new(
            origin.x() + dir.x() * t,
            origin.y() + dir.y() * t,
            origin.z() + dir.z() * t,
        );
        (tags[idx].0, hit, t)
    })
}

/// Ray-cast the scene and serialize to shaded binary PPM bytes.
///
/// Equivalent to [`render_scene_raster_shaded`] followed by
/// [`Raster::to_ppm_bytes`]; the PPM header is `P6\n<width> <height>\n255\n`
/// followed by one RGB byte per pixel. This is the analogue of
/// `V3d_View::Dump` with a shaded (non-wireframe) renderer.
pub fn render_scene_ppm_shaded(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    settings: &RenderSettings,
) -> Vec<u8> {
    render_scene_raster_shaded(scene, cam, width, height, deflection, settings).to_ppm_bytes()
}

// ---------------------------------------------------------------------------
// Rendering — shaded z-buffer rasterizer
// ---------------------------------------------------------------------------

/// A meshed scene shape in world space: vertices, per-vertex unit normals,
/// triangles and the effective material.
#[derive(Debug, Clone)]
struct ShadeMesh {
    material: Material,
    verts: Vec<GpPnt>,
    normals: Vec<GpVec>,
    triangles: Vec<(usize, usize, usize)>,
}

impl ShadeMesh {
    /// Mesh `ss`, transform every vertex into world space and compute the
    /// per-vertex normals.
    fn from_scene_shape(ss: &SceneShape, deflection: f64) -> Self {
        let mesh = crate::shape_mesh::mesh_shape(&ss.shape, deflection);
        let verts: Vec<GpPnt> = mesh
            .vertices
            .iter()
            .map(|v| v.transformed(&ss.transform))
            .collect();
        let normals = mesh_vertex_normals(&verts, &mesh.triangles);
        let triangles: Vec<(usize, usize, usize)> = mesh
            .triangles
            .iter()
            .map(|t| (t.n0, t.n1, t.n2))
            .collect();
        Self { material: ss.material(), verts, normals, triangles }
    }
}

/// Per-vertex unit normals from a triangle soup.
///
/// Each triangle's face normal is accumulated (unweighted) into its three
/// vertices and the sums are normalized, giving the smooth normals needed by
/// Gouraud/Phong shading. Exposed publicly for callers that want to reuse the
/// smoothing pass on a raw [`crate::mesh::ShapeMesh`].
pub fn mesh_vertex_normals(verts: &[GpPnt], triangles: &[Triangle]) -> Vec<GpVec> {
    let mut acc = vec![GpVec::zero(); verts.len()];
    for t in triangles {
        let n = GpVec::from_pnts(&verts[t.n0], &verts[t.n1])
            .crossed(&GpVec::from_pnts(&verts[t.n0], &verts[t.n2]))
            .normalized();
        acc[t.n0] = acc[t.n0].add(&n);
        acc[t.n1] = acc[t.n1].add(&n);
        acc[t.n2] = acc[t.n2].add(&n);
    }
    acc.iter().map(|v| v.normalized()).collect()
}

/// A screen-space z-buffer rasterizer with per-shape materials and lights.
///
/// Every shape is meshed and projected; for each projected triangle the pixels
/// inside its bounding box are filled when the barycentrically interpolated
/// view depth is nearer than the current z-buffer value. The per-pixel color
/// depends on [`ShadingMode`]:
///
/// - `Flat` — the face normal shades the whole triangle once;
/// - `Gouraud` — the vertex normals shade each vertex and the colors are
///   interpolated across the triangle;
/// - `Phong` — the vertex normals are interpolated and every pixel is shaded
///   with its own interpolated normal.
///
/// The z-buffer makes occlusion exact for the small triangles produced by the
/// shape tessellator. Ports the rasterization stage of a `V3d_View`-style
/// software pipeline.
pub fn render_scene_raster_zbuffer(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    settings: &RenderSettings,
) -> Raster {
    let mut raster = Raster::new(width, height);
    raster.clear(settings.background);
    let mut depth = vec![f64::INFINITY; width * height];

    for ss in &scene.shapes {
        let m = ShadeMesh::from_scene_shape(ss, deflection);
        for &(n0, n1, n2) in &m.triangles {
            let (Some(pa), Some(pb), Some(pc)) = (
                project_point(cam, m.verts[n0], width, height),
                project_point(cam, m.verts[n1], width, height),
                project_point(cam, m.verts[n2], width, height),
            ) else {
                continue;
            };
            let area2 = (pb.0 - pa.0) * (pc.1 - pa.1) - (pb.1 - pa.1) * (pc.0 - pa.0);
            if area2.abs() < 1e-9 {
                continue; // edge-on triangle
            }
            let xmin = pa.0.min(pb.0).min(pc.0).max(0.0) as usize;
            let xmax = pa.0.max(pb.0).max(pc.0).min(width as f64 - 1.0) as usize;
            let ymin = pa.1.min(pb.1).min(pc.1).max(0.0) as usize;
            let ymax = pa.1.max(pb.1).max(pc.1).min(height as f64 - 1.0) as usize;

            // Flat color: face normal at the first vertex.
            let face_normal = GpVec::from_pnts(&m.verts[n0], &m.verts[n1])
                .crossed(&GpVec::from_pnts(&m.verts[n0], &m.verts[n2]))
                .normalized();
            let flat_color = shade_point(
                &m.material,
                settings,
                face_normal,
                GpVec::from_pnts(&m.verts[n0], &cam.eye),
                m.verts[n0],
            );
            // Per-vertex colors for Gouraud (and normal sources for Phong).
            let smooth = settings.shading != ShadingMode::Flat;
            let vcolors = if smooth {
                [
                    shade_point(
                        &m.material,
                        settings,
                        m.normals[n0],
                        GpVec::from_pnts(&m.verts[n0], &cam.eye),
                        m.verts[n0],
                    ),
                    shade_point(
                        &m.material,
                        settings,
                        m.normals[n1],
                        GpVec::from_pnts(&m.verts[n1], &cam.eye),
                        m.verts[n1],
                    ),
                    shade_point(
                        &m.material,
                        settings,
                        m.normals[n2],
                        GpVec::from_pnts(&m.verts[n2], &cam.eye),
                        m.verts[n2],
                    ),
                ]
            } else {
                [(0.0, 0.0, 0.0); 3]
            };

            for y in ymin..=ymax {
                for x in xmin..=xmax {
                    let (px, py) = (x as f64 + 0.5, y as f64 + 0.5);
                    let w0 = (pb.0 - pa.0) * (py - pa.1) - (pb.1 - pa.1) * (px - pa.0);
                    let w1 = (pc.0 - pb.0) * (py - pb.1) - (pc.1 - pb.1) * (px - pb.0);
                    let w2 = (pa.0 - pc.0) * (py - pc.1) - (pa.1 - pc.1) * (px - pc.0);
                    let inside = (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0)
                        || (w0 <= 0.0 && w1 <= 0.0 && w2 <= 0.0);
                    if !inside {
                        continue;
                    }
                    let sum = w0 + w1 + w2;
                    if sum.abs() < 1e-12 {
                        continue;
                    }
                    // Barycentric weights (w1→A, w2→B, w0→C).
                    let la = w1 / sum;
                    let lb = w2 / sum;
                    let lc = w0 / sum;
                    let d = la * pa.2 + lb * pb.2 + lc * pc.2;
                    let i = y * width + x;
                    if d >= depth[i] {
                        continue;
                    }
                    depth[i] = d;
                    let rgb = match settings.shading {
                        ShadingMode::Flat => flat_color,
                        ShadingMode::Gouraud => (
                            la * vcolors[0].0 + lb * vcolors[1].0 + lc * vcolors[2].0,
                            la * vcolors[0].1 + lb * vcolors[1].1 + lc * vcolors[2].1,
                            la * vcolors[0].2 + lb * vcolors[1].2 + lc * vcolors[2].2,
                        ),
                        ShadingMode::Phong => {
                            let n = m.normals[n0]
                                .multiplied_scalar(la)
                                .add(&m.normals[n1].multiplied_scalar(lb))
                                .add(&m.normals[n2].multiplied_scalar(lc))
                                .normalized();
                            let hit = GpPnt::new(
                                la * m.verts[n0].x() + lb * m.verts[n1].x() + lc * m.verts[n2].x(),
                                la * m.verts[n0].y() + lb * m.verts[n1].y() + lc * m.verts[n2].y(),
                                la * m.verts[n0].z() + lb * m.verts[n1].z() + lc * m.verts[n2].z(),
                            );
                            let view_dir = GpVec::from_pnts(&hit, &cam.eye);
                            shade_point(&m.material, settings, n, view_dir, hit)
                        }
                    };
                    raster.set_pixel(x, y, rgb);
                }
            }
        }
    }
    raster
}

/// Render a depth map of the scene as a grayscale [`Raster`].
///
/// Every pixel holds the view depth (distance along the camera forward axis)
/// of the nearest surface, mapped linearly between the camera's near plane
/// (white) and far plane (black). Pixels with no geometry stay black. This is
/// the classic "depth buffer" debug view of a z-buffer renderer: nearer
/// geometry appears brighter.
pub fn render_scene_depth(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
) -> Raster {
    let meshes = scene_shade_meshes(scene, deflection);
    let (tri_pts, _tags) = flatten_shade_meshes(&meshes);
    let bvh = build_tri_bvh(&tri_pts, 8);
    let range = (cam.far - cam.near).max(1e-6);
    let fwd = cam.forward();
    let mut raster = Raster::new(width, height);
    raster.fill(|x, y| {
        let (nx, ny) = screen_to_ndc(width, height, x as f64 + 0.5, y as f64 + 0.5);
        let (origin, dir) = cam.ray_through_ndc(nx, ny);
        match bvh_ray_cast(&bvh, &tri_pts, &origin, &dir) {
            Some((_, t)) => {
                let hit = GpPnt::new(
                    origin.x() + dir.x() * t,
                    origin.y() + dir.y() * t,
                    origin.z() + dir.z() * t,
                );
                let depth_z = GpVec::from_pnts(&cam.eye, &hit).dot(&fwd);
                let v = (1.0 - ((depth_z - cam.near) / range).clamp(0.0, 1.0)).clamp(0.0, 1.0);
                (v, v, v)
            }
            None => (0.0, 0.0, 0.0),
        }
    });
    raster
}

/// Render the depth map and serialize it to binary PPM bytes.
///
/// Equivalent to [`render_scene_depth`] followed by [`Raster::to_ppm_bytes`].
pub fn render_scene_ppm_depth(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
) -> Vec<u8> {
    render_scene_depth(scene, cam, width, height, deflection).to_ppm_bytes()
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Camera-space `(x, y)` in `[-1, 1]²` for a screen point `(sx, sy)` with the
/// origin at the top-left of a `width`×`height` image (y flipped).
fn screen_to_ndc(width: usize, height: usize, sx: f64, sy: f64) -> (f64, f64) {
    let x_ndc = 2.0 * sx / width as f64 - 1.0;
    let y_ndc = 1.0 - 2.0 * sy / height as f64;
    (x_ndc, y_ndc)
}

/// Perspective scale `f = 1/tan(fov/2)`.
fn perspective_scale(cam: &Camera) -> f64 {
    (cam.fov_deg.to_radians() * 0.5).tan().recip()
}

/// Orthographic viewport half-height in world units.
///
/// Matches the perspective frustum width at the target distance so switching
/// projection keeps the on-screen size approximately constant.
fn ortho_scale(cam: &Camera) -> f64 {
    let s = cam.eye.distance(&cam.target) * (cam.fov_deg.to_radians() * 0.5).tan();
    if s.abs() < 1e-12 { 1.0 } else { s }
}

/// The orthonormal camera basis `(right, up2, forward)`.
///
/// Handles a degenerate `forward` (eye == target) by defaulting to `−Z`, and
/// an `up` parallel to `forward` by picking a reference axis not parallel to
/// the forward direction.
fn camera_basis(cam: &Camera) -> (GpVec, GpVec, GpVec) {
    let v = GpVec::from_pnts(&cam.eye, &cam.target);
    let fwd = if v.square_magnitude() < 1e-24 {
        GpVec::new(0.0, 0.0, -1.0)
    } else {
        v.normalized()
    };
    let mut right = fwd.crossed(&cam.up);
    if right.square_magnitude() < 1e-24 {
        let reference = if fwd.z().abs() < 0.9 { GpVec::new(0.0, 0.0, 1.0) } else { GpVec::new(1.0, 0.0, 0.0) };
        right = fwd.crossed(&reference);
    }
    let right = right.normalized();
    let up2 = right.crossed(&fwd);
    (right, up2, fwd)
}

/// Rotate `v` about `axis` by `angle` radians (Rodrigues' rotation formula).
///
/// `axis` need not be unit — it is normalized internally. `v` is returned
/// unchanged when `axis` is degenerate (zero).
fn rotate_vec(v: &GpVec, axis: &GpVec, angle: f64) -> GpVec {
    let k = axis.normalized();
    if k.square_magnitude() < 1e-24 {
        return *v;
    }
    let s = angle.sin();
    let c = angle.cos();
    let dot = v.dot(&k);
    let cross = k.crossed(v);
    k.multiplied_scalar(dot * (1.0 - c))
        .add(&v.multiplied_scalar(c))
        .add(&cross.multiplied_scalar(s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::render_svg::svg_polygon_count;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0
    }

    fn translate(x: f64, y: f64, z: f64) -> GpTrsf {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(x, y, z));
        t
    }

    fn box_scene() -> VizScene {
        let mut s = VizScene::new();
        s.add(SceneShape::new(unit_box()));
        s
    }

    /// Screen-space bbox width of shape `idx` under `cam`.
    fn projected_screen_width(
        scene: &VizScene,
        idx: usize,
        cam: &Camera,
        w: usize,
        h: usize,
        deflection: f64,
    ) -> f64 {
        let mut minx = f64::INFINITY;
        let mut maxx = f64::NEG_INFINITY;
        for v in shape_world_vertices(scene, idx, deflection) {
            if let Some((sx, _sy, _d)) = project_point(cam, v, w, h) {
                minx = minx.min(sx);
                maxx = maxx.max(sx);
            }
        }
        maxx - minx
    }

    #[test]
    fn camera_look_at_origin() {
        let cam = Camera::default();
        // The origin maps to the screen centre.
        let (ox, oy, _od) = project_point(&cam, GpPnt::zero(), 200, 150).unwrap();
        assert!((ox - 100.0).abs() < 1.0, "origin x {ox}");
        assert!((oy - 75.0).abs() < 1.0, "origin y {oy}");
        // +X world maps to the right half of the screen.
        let (px, _py, _pd) = project_point(&cam, GpPnt::new(1.0, 0.0, 0.0), 200, 150).unwrap();
        assert!(px > 101.0, "+X point maps right, got x={px}");
        // +Y world maps to the upper half (screen y decreases upward).
        let (_qx, qy, _qd) = project_point(&cam, GpPnt::new(0.0, 1.0, 0.0), 200, 150).unwrap();
        assert!(qy < 74.0, "+Y point maps to upper half, got y={qy}");
    }

    #[test]
    fn camera_projection_behind() {
        let cam = Camera::default();
        // z > eye.z is behind the camera looking along −Z.
        assert!(project_point(&cam, GpPnt::new(0.0, 0.0, 6.0), 200, 150).is_none());
        assert!(project_point(&cam, GpPnt::new(0.0, 0.0, 5.5), 200, 150).is_none());
    }

    #[test]
    fn view_matrix_rotation() {
        let cam = Camera {
            eye: GpPnt::zero(),
            target: GpPnt::new(0.0, 0.0, -1.0),
            up: GpVec::new(0.0, 1.0, 0.0),
            ..Default::default()
        };
        let m = camera_view_matrix(&cam);
        // +X world → +X camera, +Y world → +Y camera, +Z world → −Z camera.
        let vx = m.multiplied(&GpXyz::new(1.0, 0.0, 0.0));
        assert!((vx.x - 1.0).abs() < 1e-12 && vx.y.abs() < 1e-12 && vx.z.abs() < 1e-12, "vx {vx:?}");
        let vy = m.multiplied(&GpXyz::new(0.0, 1.0, 0.0));
        assert!(vy.x.abs() < 1e-12 && (vy.y - 1.0).abs() < 1e-12 && vy.z.abs() < 1e-12, "vy {vy:?}");
        let vz = m.multiplied(&GpXyz::new(0.0, 0.0, 1.0));
        assert!(vz.x.abs() < 1e-12 && vz.y.abs() < 1e-12 && (vz.z + 1.0).abs() < 1e-12, "vz {vz:?}");
    }

    #[test]
    fn scene_add_clear() {
        let mut scene = VizScene::new();
        assert!(scene.is_empty());
        scene.add(SceneShape::new(unit_box()));
        scene.add(SceneShape::new(unit_box()));
        assert_eq!(scene.shapes.len(), 2);
        assert_eq!(scene.len(), 2);
        assert!(!scene.is_empty());
        assert!(scene.get(1).is_some());
        scene.remove(1);
        assert_eq!(scene.len(), 1);
        scene.clear();
        assert!(scene.is_empty());
    }

    #[test]
    fn render_svg_contains_polygons() {
        let scene = box_scene();
        let svg = render_scene_svg(&scene, &Camera::default(), 400, 300, 0.25);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("viewBox=\"0 0 400 300\""));
        assert!(svg.contains("<polygon"));
        assert!(svg_polygon_count(&svg) > 0);
    }

    #[test]
    fn render_ppm_dimensions() {
        let scene = box_scene();
        let ppm = render_scene_ppm(&scene, &Camera::default(), 64, 48, 0.25);
        assert!(ppm.starts_with(b"P6\n64 48\n255\n"));
        // 13-byte header ("P6\n", "64 48\n", "255\n") + one byte per component.
        assert_eq!(ppm.len(), 13 + 64 * 48 * 3);
    }

    #[test]
    fn render_ppm_not_blank() {
        let scene = box_scene();
        let ppm = render_scene_ppm(&scene, &Camera::default(), 64, 48, 0.25);
        let bg = [8u8, 8, 14];
        let has_hit = ppm[13..].chunks_exact(3).any(|p| p != bg);
        assert!(has_hit, "expected at least one non-background pixel");
    }

    #[test]
    fn render_zbuffer_not_blank() {
        let scene = box_scene();
        let raster = render_scene_zbuffer_raster(&scene, &Camera::default(), 64, 48, 0.25);
        let bg = (0.03, 0.03, 0.05);
        let has_hit = raster.pixels.chunks_exact(3).any(|p| p != [bg.0, bg.1, bg.2]);
        assert!(has_hit, "expected at least one shaded pixel");
        // Same header convention as the ray-cast PPM path.
        let ppm = raster.to_ppm_bytes();
        assert!(ppm.starts_with(b"P6\n64 48\n255\n"));
    }

    #[test]
    fn orthographic_parallel() {
        let persp = Camera::default();
        let mut ortho = Camera::default();
        ortho.projection = CameraProjection::Orthographic;
        let a = GpPnt::new(1.0, 1.0, 0.0);
        let b = GpPnt::new(1.0, 1.0, 2.0);
        // Orthographic: same (x, y) regardless of depth.
        let (oa, _oa_y, _oa_d) = project_point(&ortho, a, 200, 150).unwrap();
        let (ob, _ob_y, _ob_d) = project_point(&ortho, b, 200, 150).unwrap();
        assert!((oa - ob).abs() < 1e-9, "ortho x {oa} vs {ob}");
        // Perspective: different depth gives different screen x.
        let (pa, _pa_y, _pa_d) = project_point(&persp, a, 200, 150).unwrap();
        let (pb, _pb_y, _pb_d) = project_point(&persp, b, 200, 150).unwrap();
        assert!((pa - pb).abs() > 1e-3, "perspective x {pa} vs {pb}");
    }

    #[test]
    fn scene_bounds_box() {
        let scene = box_scene();
        let (lo, hi) = scene_bounds(&scene).unwrap();
        assert!((lo.x()).abs() < 1e-9 && (lo.y()).abs() < 1e-9 && (lo.z()).abs() < 1e-9);
        assert!((hi.x() - 1.0).abs() < 1e-9 && (hi.y() - 1.0).abs() < 1e-9 && (hi.z() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn world_vertices_transformed() {
        let mut scene = VizScene::new();
        let mut ss = SceneShape::new(unit_box());
        ss.transform = translate(10.0, 0.0, 0.0);
        scene.add(ss);
        let verts = shape_world_vertices(&scene, 0, 0.25);
        assert!(!verts.is_empty());
        let minx = verts.iter().map(|p| p.x()).fold(f64::INFINITY, f64::min);
        let maxx = verts.iter().map(|p| p.x()).fold(f64::NEG_INFINITY, f64::max);
        assert!((minx - 10.0).abs() < 1e-6, "minx {minx}");
        assert!((maxx - 11.0).abs() < 1e-6, "maxx {maxx}");
        // Every world vertex carries the translation in x.
        assert!(verts.iter().all(|p| p.x() >= 9.999 && p.x() <= 11.001));
    }

    #[test]
    fn perspective_vanishing() {
        let mut scene = VizScene::new();
        scene.add(SceneShape::new(unit_box())); // near: z in [0, 1]
        let mut far = SceneShape::new(unit_box());
        far.transform = translate(0.0, 0.0, -5.0); // far: z in [-5, -4]
        scene.add(far);
        let cam = Camera::default();
        let near_w = projected_screen_width(&scene, 0, &cam, 400, 300, 0.25);
        let far_w = projected_screen_width(&scene, 1, &cam, 400, 300, 0.25);
        assert!(near_w > far_w, "near {near_w} should be wider than far {far_w}");
    }

    fn sphere_scene() -> VizScene {
        let mut s = VizScene::new();
        s.add(SceneShape::new(BRepPrimSphere::make_sphere(1.0).solid.0));
        s
    }

    /// Colors of the pixels in a square neighbourhood around `(cx, cy)`.
    fn sample_region(raster: &Raster, cx: f64, cy: f64, r: i32) -> Vec<(f64, f64, f64)> {
        let mut out = Vec::new();
        for dy in -r..=r {
            for dx in -r..=r {
                let x = cx as i32 + dx;
                let y = cy as i32 + dy;
                if x >= 0 && y >= 0 && x < raster.width as i32 && y < raster.height as i32 {
                    out.push(raster.get_pixel(x as usize, y as usize));
                }
            }
        }
        out
    }

    /// Number of distinct colors after a coarse 5-bit quantization.
    fn count_distinct(raster: &Raster) -> usize {
        let mut set = std::collections::HashSet::new();
        for chunk in raster.pixels.chunks_exact(3) {
            set.insert(((chunk[0] * 31.0) as u8, (chunk[1] * 31.0) as u8, (chunk[2] * 31.0) as u8));
        }
        set.len()
    }

    #[test]
    fn phong_shading_light_hits() {
        let mat = Material::default();
        let light = Light {
            position: GpPnt::new(0.0, 0.0, 2.0),
            color: (1.0, 1.0, 1.0),
            intensity: 1.0,
        };
        // Front-lit: normal, light and view all along +Z.
        let front = shade_phong(
            &mat,
            &light,
            GpVec::new(0.0, 0.0, 1.0),
            GpVec::new(0.0, 0.0, 1.0),
            GpVec::new(0.0, 0.0, 1.0),
            0.15,
        );
        // Back-lit: the normal points away from the light → no diffuse, and
        // (because the light is behind) no specular either.
        let back = shade_phong(
            &mat,
            &light,
            GpVec::new(0.0, 0.0, -1.0),
            GpVec::new(0.0, 0.0, 1.0),
            GpVec::new(0.0, 0.0, 1.0),
            0.15,
        );
        let ambient = mat.diffuse.0 * 0.15;
        assert!(front.0 > ambient, "front-lit {:.3} should exceed ambient {:.3}", front.0, ambient);
        assert!(
            (back.0 - ambient).abs() < 1e-9,
            "back-lit should be ambient only, got {:.3}",
            back.0
        );
        // Straight-on view adds a specular spike on the lit face.
        assert!(front.0 > mat.diffuse.0, "front-lit includes specular, got {:.3}", front.0);
    }

    #[test]
    fn phong_specular_high_shininess() {
        let high = Material { shininess: 128.0, ..Default::default() };
        let low = Material { shininess: 4.0, ..Default::default() };
        let light = Light {
            position: GpPnt::new(0.0, 0.0, 2.0),
            color: (1.0, 1.0, 1.0),
            intensity: 1.0,
        };
        let normal = GpVec::new(0.0, 0.0, 1.0);
        let light_dir = GpVec::new(0.0, 0.0, 1.0);
        let view_mirror = GpVec::new(0.0, 0.0, 1.0);
        let view_off = GpVec::new(1.0, 0.0, 1.0).normalized(); // 45° off the mirror
        let hi_mirror = shade_phong(&high, &light, normal, view_mirror, light_dir, 0.15).0;
        let hi_off = shade_phong(&high, &light, normal, view_off, light_dir, 0.15).0;
        let lo_mirror = shade_phong(&low, &light, normal, view_mirror, light_dir, 0.15).0;
        let lo_off = shade_phong(&low, &light, normal, view_off, light_dir, 0.15).0;
        // The highlight peaks at the mirror direction for both materials.
        assert!(hi_mirror > hi_off, "mirror {hi_mirror:.3} should beat off-axis {hi_off:.3}");
        assert!(lo_mirror > lo_off, "mirror {lo_mirror:.3} should beat off-axis {lo_off:.3}");
        // High shininess collapses 45° off-axis; low shininess keeps a wide glow.
        assert!(lo_off > hi_off, "low shininess keeps a wider highlight {lo_off:.3} vs {hi_off:.3}");
        let hi_drop = hi_mirror - hi_off;
        let lo_drop = lo_mirror - lo_off;
        assert!(hi_drop > lo_drop, "high shininess should fall off faster ({hi_drop:.3} vs {lo_drop:.3})");
    }

    #[test]
    fn render_shaded_not_blank() {
        let scene = box_scene();
        // A slightly off-axis camera exposes the side faces, which are darker
        // than the front face under the key light.
        let cam = Camera::look_at(
            GpPnt::new(1.6, 1.2, 5.0),
            GpPnt::new(0.5, 0.5, 0.5),
            GpVec::new(0.0, 1.0, 0.0),
        );
        let settings = RenderSettings::default();
        let ppm = render_scene_ppm_shaded(&scene, &cam, 64, 48, 0.25, &settings);
        assert!(ppm.starts_with(b"P6\n64 48\n255\n"));
        let bg = [
            (settings.background.0 * 255.0) as u8,
            (settings.background.1 * 255.0) as u8,
            (settings.background.2 * 255.0) as u8,
        ];
        let body = &ppm[13..];
        let has_hit = body.chunks_exact(3).any(|p| p != bg);
        assert!(has_hit, "expected at least one non-background pixel");
        let max_r = body.chunks_exact(3).map(|p| p[0]).max().unwrap_or(0);
        let has_dark = body
            .chunks_exact(3)
            .any(|p| p[0] < 100 && p[0] < max_r.saturating_sub(30));
        assert!(has_dark, "expected a dark side/back pixel (max r {max_r})");
    }

    #[test]
    fn material_per_shape_color() {
        let mut scene = VizScene::new();
        let mut red = SceneShape::new(unit_box());
        red.transform = translate(-1.5, 0.0, 0.0);
        red.material = Material::from_diffuse((0.8, 0.1, 0.1));
        scene.add(red);
        let mut blue = SceneShape::new(unit_box());
        blue.transform = translate(1.5, 0.0, 0.0);
        blue.material = Material::from_diffuse((0.1, 0.1, 0.8));
        scene.add(blue);
        let raster = render_scene_raster_zbuffer(
            &scene,
            &Camera::default(),
            96,
            72,
            0.25,
            &RenderSettings::default(),
        );
        let cam = Camera::default();
        // Centres of the two front faces.
        let (lx, ly, _) = project_point(&cam, GpPnt::new(-1.5, 0.5, 1.0), 96, 72).unwrap();
        let (rx, ry, _) = project_point(&cam, GpPnt::new(1.5, 0.5, 1.0), 96, 72).unwrap();
        let left = sample_region(&raster, lx, ly, 2);
        let right = sample_region(&raster, rx, ry, 2);
        let left_red = left.iter().filter(|&&(r, _g, b)| r > b && r > 0.2).count();
        let right_blue = right.iter().filter(|&&(r, _g, b)| b > r && b > 0.2).count();
        assert!(left_red >= 3, "left half should be red-dominant, got {left_red}/{}", left.len());
        assert!(right_blue >= 3, "right half should be blue-dominant, got {right_blue}/{}", right.len());
    }

    #[test]
    fn depth_buffer_occludes() {
        let mut scene = VizScene::new();
        let mut near = SceneShape::new(unit_box());
        near.material = Material::from_diffuse((0.9, 0.1, 0.1));
        scene.add(near);
        let mut far = SceneShape::new(unit_box());
        far.transform = translate(0.0, 0.0, -4.0);
        far.material = Material::from_diffuse((0.1, 0.1, 0.9));
        scene.add(far);
        let raster = render_scene_raster_zbuffer(
            &scene,
            &Camera::default(),
            96,
            72,
            0.25,
            &RenderSettings::default(),
        );
        // The far box's front face centre projects here.
        let cam = Camera::default();
        let (fx, fy, _) = project_point(&cam, GpPnt::new(0.5, 0.5, -3.0), 96, 72).unwrap();
        let region = sample_region(&raster, fx, fy, 5);
        assert!(!region.is_empty(), "sampled region should not be empty");
        let blueish = region.iter().filter(|&&(r, _g, b)| b > r + 0.1).count();
        assert_eq!(blueish, 0, "far box should be fully occluded; found {blueish} blue pixels");
    }

    #[test]
    fn depth_map_nearer_is_brighter() {
        let mut scene = VizScene::new();
        let mut near = SceneShape::new(unit_box());
        near.transform = translate(-2.0, 0.0, 0.0); // z in [0, 1], close
        scene.add(near);
        let mut far = SceneShape::new(unit_box());
        far.transform = translate(2.0, 0.0, -3.0); // z in [-3, -2], far
        scene.add(far);
        let mut cam = Camera::default();
        cam.far = 10.0; // tighter range so the depth gradient is visible
        let raster = render_scene_depth(&scene, &cam, 96, 72, 0.25);
        // Centres of the two visible front faces.
        let (nx, ny, _) = project_point(&cam, GpPnt::new(-1.5, 0.5, 1.0), 96, 72).unwrap();
        let (fx, fy, _) = project_point(&cam, GpPnt::new(2.5, 0.5, -2.0), 96, 72).unwrap();
        let near_max = sample_region(&raster, nx, ny, 2).iter().map(|c| c.0).fold(0.0, f64::max);
        let far_max = sample_region(&raster, fx, fy, 2).iter().map(|c| c.0).fold(0.0, f64::max);
        assert!(
            near_max > far_max,
            "nearer geometry should be brighter: {near_max:.3} vs {far_max:.3}"
        );
    }

    #[test]
    fn camera_orbit_keeps_target() {
        let mut cam = Camera::default();
        let target_before = cam.target;
        cam.camera_orbit(90.0, 0.0);
        // The target never moves.
        assert!(cam.target.distance(&target_before) < 1e-9);
        // The eye keeps its distance from the target.
        assert!((cam.eye.distance(&cam.target) - 5.0).abs() < 1e-9, "distance {}", cam.eye.distance(&cam.target));
        // The target still projects to the screen centre.
        let (sx, sy, _) = project_point(&cam, cam.target, 200, 150).unwrap();
        assert!((sx - 100.0).abs() < 1.0, "target x {sx}");
        assert!((sy - 75.0).abs() < 1.0, "target y {sy}");
    }

    #[test]
    fn camera_zoom_changes_distance() {
        let mut cam = Camera::default();
        let d0 = cam.eye.distance(&cam.target);
        cam.camera_zoom(2.0);
        let d1 = cam.eye.distance(&cam.target);
        assert!((d1 - d0 / 2.0).abs() < 1e-9, "distance {d0} -> {d1}");
        // The origin still projects to the screen centre.
        let (sx, sy, _) = project_point(&cam, GpPnt::zero(), 200, 150).unwrap();
        assert!((sx - 100.0).abs() < 1.0, "origin x {sx}");
        assert!((sy - 75.0).abs() < 1.0, "origin y {sy}");
    }

    #[test]
    fn camera_pan_moves_view() {
        let mut cam = Camera::default();
        let before = project_point(&cam, GpPnt::zero(), 200, 150).unwrap();
        cam.camera_pan(50.0, 0.0);
        let after = project_point(&cam, GpPnt::zero(), 200, 150).unwrap();
        assert!(
            after.0 > before.0,
            "pan right should shift the origin right: {} -> {}",
            before.0,
            after.0
        );
        // The target moves with the eye, so it stays centred on screen.
        let (tx, ty, _) = project_point(&cam, cam.target, 200, 150).unwrap();
        assert!((tx - 100.0).abs() < 1.0 && (ty - 75.0).abs() < 1.0);
    }

    #[test]
    fn camera_ray_through_center() {
        let cam = Camera::default();
        let (origin, dir) = cam.camera_ray(200, 150, 100, 75);
        // Closest point on the ray to the target.
        let t = GpVec::from_pnts(&origin, &cam.target).dot(&dir);
        let closest = GpPnt::new(
            origin.x() + dir.x() * t,
            origin.y() + dir.y() * t,
            origin.z() + dir.z() * t,
        );
        let dist = closest.distance(&cam.target);
        // The centre pixel's sub-pixel ray is offset by half a pixel from the
        // exact optical axis, so the ray passes within a fraction of a world
        // unit of the target rather than exactly through it.
        assert!(dist < 0.1, "centre ray should pass near the target, dist {dist}");
    }

    #[test]
    fn shading_modes_flat_vs_gouraud() {
        let scene = sphere_scene();
        let mut flat_settings = RenderSettings::default();
        flat_settings.shading = ShadingMode::Flat;
        let flat = render_scene_raster_zbuffer(
            &scene,
            &Camera::default(),
            64,
            48,
            0.25,
            &flat_settings,
        );
        let mut gouraud_settings = RenderSettings::default();
        gouraud_settings.shading = ShadingMode::Gouraud;
        let gouraud = render_scene_raster_zbuffer(
            &scene,
            &Camera::default(),
            64,
            48,
            0.25,
            &gouraud_settings,
        );
        let n_flat = count_distinct(&flat);
        let n_gouraud = count_distinct(&gouraud);
        assert!(
            n_gouraud > n_flat,
            "Gouraud should smooth away flat banding: flat {n_flat}, gouraud {n_gouraud}"
        );
    }

    #[test]
    fn pick_center_hits_front_shape() {
        let scene = box_scene();
        let cam = Camera::default();
        // The box lives in the +X/+Y/+Z octant, so aim at the projected
        // centre of its front face rather than the image centre.
        let (sx, sy, _) = project_point(&cam, GpPnt::new(0.5, 0.5, 1.0), 64, 48).unwrap();
        let (idx, hit, _t) = pick_shape(&scene, &cam, 64, 48, sx as usize, sy as usize, 0.25).unwrap();
        assert_eq!(idx, 0, "pick returns the only shape");
        assert!((hit.z() - 1.0).abs() < 1e-6, "front face sits at z=1, got {}", hit.z());
    }

    #[test]
    fn material_scene_accessors() {
        let mut scene = box_scene();
        assert!(scene.shape_material(0).is_some());
        assert!(scene.shape_material(1).is_none());
        let red = Material::from_diffuse((0.9, 0.2, 0.1));
        assert!(scene.set_material(0, red));
        assert!(!scene.set_material(9, red));
        let m = scene.shape_material(0).unwrap();
        assert!((m.diffuse.0 - 0.9).abs() < 1e-9, "material updated, diffuse {}", m.diffuse.0);
    }
}
