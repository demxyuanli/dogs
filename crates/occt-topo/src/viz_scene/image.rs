use super::prelude::*;
use super::*;

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
pub(super) fn default_color() -> (f64, f64, f64) {
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
pub(super) struct ProjectedTri {
    pub(super) depth: f64,
    pub(super) a: (f64, f64),
    pub(super) b: (f64, f64),
    pub(super) c: (f64, f64),
    pub(super) color: (f64, f64, f64),
}

/// Project every scene triangle and return them sorted far-to-near.
///
/// Triangles behind the camera or with a degenerate (zero-area) projection
/// are dropped.
pub(super) fn project_scene(scene: &VizScene, cam: &Camera, width: usize, height: usize, deflection: f64) -> Vec<ProjectedTri> {
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
pub(super) fn color_hex(c: (f64, f64, f64)) -> String {
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
pub(super) fn shade_lambert(tri: &SceneTriangle, light: &GpVec) -> (f64, f64, f64) {
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
pub(super) fn shade_point(
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
pub(super) fn scene_shade_meshes(scene: &VizScene, deflection: f64) -> Vec<ShadeMesh> {
    scene
        .shapes
        .iter()
        .map(|ss| ShadeMesh::from_scene_shape(ss, deflection))
        .collect()
}

/// Flatten every mesh of `meshes` into one triangle buffer with `(shape,
/// triangle)` tags, ready for [`build_tri_bvh`] + [`bvh_ray_cast`].
pub(super) fn flatten_shade_meshes(meshes: &[ShadeMesh]) -> (Vec<(GpPnt, GpPnt, GpPnt)>, Vec<(usize, usize)>) {
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
pub(super) fn triangle_barycentric(a: &GpPnt, b: &GpPnt, c: &GpPnt, p: &GpPnt) -> (f64, f64, f64) {
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

/// Pick the world-space hit point under a screen pixel.
///
/// Convenience wrapper over [`pick_shape`] that discards the shape index and
/// the ray parameter and returns only the hit point. Returns `None` when the
/// ray through `(px, py)` misses every shape. Ports the "give me the point"
/// form of `V3d_View::Convert` used by interactive snapping.
pub fn pick_point(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    px: usize,
    py: usize,
    deflection: f64,
) -> Option<GpPnt> {
    pick_shape(scene, cam, width, height, px, py, deflection).map(|(_, hit, _)| hit)
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
