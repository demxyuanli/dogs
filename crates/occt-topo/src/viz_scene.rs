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

use std::cmp::Ordering;

use occt_core::bnd::BndBox;
use occt_core::bvh::bvh_ops::bvh_ray_cast;
use occt_core::bvh::builder_tri::build_tri_bvh;
use occt_core::gp::{GpMat, GpPnt, GpTrsf, GpVec, GpXyz};

use crate::mesh::ShapeMesh;
use crate::shape::TopoShape;

// ---------------------------------------------------------------------------
// Scene graph
// ---------------------------------------------------------------------------

/// One item in a [`VizScene`]: a shape, its world transform and an optional
/// per-shape color (normalized RGB in `[0, 1]`).
///
/// Mirrors OCCT's `AIS_Shape` which attaches a `TopoDS_Shape` and a display
/// transform to the interactive context.
#[derive(Debug, Clone)]
pub struct SceneShape {
    pub shape: TopoShape,
    pub transform: GpTrsf,
    pub color: Option<(f64, f64, f64)>,
}

impl SceneShape {
    /// A shape at the identity transform with no explicit color.
    pub fn new(shape: TopoShape) -> Self {
        Self { shape, transform: GpTrsf::identity(), color: None }
    }

    /// Copy with an explicit color.
    pub fn with_color(mut self, color: (f64, f64, f64)) -> Self {
        self.color = Some(color);
        self
    }

    /// Copy with a world transform.
    pub fn with_transform(mut self, transform: GpTrsf) -> Self {
        self.transform = transform;
        self
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;
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
}
