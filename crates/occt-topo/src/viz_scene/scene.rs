use super::prelude::*;
use super::*;

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
/// `material.diffuse` (see [`SceneShape::material`]). The optional `uv` field
/// carries per-vertex texture coordinates (in the mesh's local parameter
/// space) so a shape imported from an OBJ/glTF/VRML node — whose [`MeshNode`]
/// stores `uv` — can be rendered with a texture; when `None`, the textured
/// renderer synthesizes a planar projection (see
/// [`textured_mesh_from_scene_shape`]).
#[derive(Debug, Clone)]
pub struct SceneShape {
    pub shape: TopoShape,
    pub transform: GpTrsf,
    pub color: Option<(f64, f64, f64)>,
    /// Surface material used by the shaded renderers.
    pub material: Material,
    /// Optional per-vertex texture coordinates, one [`GpPnt2d`] per mesh
    /// vertex (same length as the shape's tessellation). `None` (the default)
    /// means "no UV — synthesize a planar projection".
    pub uv: Option<Vec<GpPnt2d>>,
}

impl SceneShape {
    /// A shape at the identity transform with the default material and no
    /// explicit color.
    pub fn new(shape: TopoShape) -> Self {
        Self {
            shape,
            transform: GpTrsf::identity(),
            color: None,
            material: Material::default(),
            uv: None,
        }
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

    /// Copy with per-vertex texture coordinates (one [`GpPnt2d`] per mesh
    /// vertex). Pass `None` to fall back to synthesized planar projection UVs.
    pub fn with_uv(mut self, uv: Option<Vec<GpPnt2d>>) -> Self {
        self.uv = uv;
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
