use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Rendering — shaded z-buffer rasterizer
// ---------------------------------------------------------------------------

/// A meshed scene shape in world space: vertices, per-vertex unit normals,
/// triangles and the effective material.
#[derive(Debug, Clone)]
pub(super) struct ShadeMesh {
    pub(super) material: Material,
    pub(super) verts: Vec<GpPnt>,
    pub(super) normals: Vec<GpVec>,
    pub(super) triangles: Vec<(usize, usize, usize)>,
}

impl ShadeMesh {
    /// Mesh `ss`, transform every vertex into world space and compute the
    /// per-vertex normals.
    pub(super) fn from_scene_shape(ss: &SceneShape, deflection: f64) -> Self {
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

/// A selection highlight material: bright yellow diffuse with a matching
/// emissive so the shape reads as "selected" even in deep shadow.
///
/// Preserves the base material's specular, shininess and opacity so the
/// highlighted shape keeps its surface character — only the tint changes.
/// Mirrors the default highlight color OCCT applies to a picked `AIS_Shape`.
pub fn highlight_material(base: &Material) -> Material {
    Material {
        diffuse: (1.0, 0.85, 0.2),
        specular: base.specular,
        emissive: (0.3, 0.22, 0.05),
        shininess: base.shininess,
        opacity: base.opacity,
    }
}

/// Render the scene with the selected shapes highlighted, as binary PPM bytes.
///
/// The scene is rendered with the same z-buffer pipeline as
/// [`render_scene_raster_zbuffer`] (same materials, lights and shading modes),
/// except every shape whose index appears in `selection` is given the
/// [`highlight_material`] before rasterizing. Occlusion is unchanged — a
/// selected shape hidden behind an unselected one stays hidden, exactly as in
/// the unselected render. Ports the selection highlight of
/// `AIS_InteractiveContext::SetSelected`, which recolors the chosen
/// `AIS_Shape`s.
pub fn render_scene_with_selection(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    settings: &RenderSettings,
    selection: &[usize],
) -> Vec<u8> {
    let selected: std::collections::HashSet<usize> = selection.iter().copied().collect();
    let mut highlighted = scene.clone();
    for (i, ss) in highlighted.shapes.iter_mut().enumerate() {
        if selected.contains(&i) {
            ss.material = highlight_material(&ss.material);
        }
    }
    render_scene_raster_zbuffer(&highlighted, cam, width, height, deflection, settings).to_ppm_bytes()
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
// Texture mapping
// ---------------------------------------------------------------------------

/// A width×height RGB texture image.
///
/// `pixels` holds one `[r, g, b]` byte triple per texel in row-major order, so
/// `pixels.len() == width * height`. This is the minimal surface image needed
/// by [`render_textured_ppm`]; it is the analogue of the pixel array behind an
/// OCCT `Image_PixMap` / `Graphic3d_Texture2D`.
#[derive(Debug, Clone)]
pub struct Texture {
    /// Texture width in texels.
    pub width: usize,
    /// Texture height in texels.
    pub height: usize,
    /// Row-major RGB texels, one `[r, g, b]` byte triple per texel.
    pub pixels: Vec<[u8; 3]>,
}

/// Build a procedural checkerboard texture.
///
/// The `width`×`height` image is divided into a `cells`×`cells` grid (each
/// cell is a square block of texels); adjacent cells alternate between `c1`
/// and `c2`. A cell `(cx, cy)` is `c1` when `(cx + cy)` is even and `c2`
/// otherwise. No file I/O is involved — the image is generated in memory.
pub fn checkerboard_texture(
    width: usize,
    height: usize,
    cells: usize,
    c1: (u8, u8, u8),
    c2: (u8, u8, u8),
) -> Texture {
    let cells = cells.max(1);
    let w = width.max(1);
    let h = height.max(1);
    let mut pixels = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            let cx = x * cells / w;
            let cy = y * cells / h;
            let p = if (cx + cy) % 2 == 0 { c1 } else { c2 };
            pixels.push([p.0, p.1, p.2]);
        }
    }
    Texture { width: w, height: h, pixels }
}

/// Build a solid-color texture: every texel is `c`.
pub fn solid_texture(width: usize, height: usize, c: (u8, u8, u8)) -> Texture {
    let pixels = vec![[c.0, c.1, c.2]; width * height];
    Texture { width, height, pixels }
}

/// A world-space textured triangle mesh, the geometry a textured renderer
/// consumes.
///
/// `uv` holds one texture coordinate per vertex (parallel to `vertices`); a
/// triangle's texture coordinate at any interior point is the barycentric
/// interpolation of its three vertex UVs. `texture` is the per-mesh surface
/// image (kept optional so a shape can fall back to its plain `color`), and
/// `color` is the normalized RGB used when no texture is available.
#[derive(Debug, Clone)]
pub struct TexturedMesh {
    /// World-space vertices.
    pub vertices: Vec<GpPnt>,
    /// Triangle indices into `vertices`.
    pub triangles: Vec<(usize, usize, usize)>,
    /// One texture coordinate per vertex (parallel to `vertices`).
    pub uv: Vec<GpPnt2d>,
    /// Optional surface image; `None` means "use `color`".
    pub texture: Option<Texture>,
    /// Normalized RGB base color used when `texture` is `None`.
    pub color: (f64, f64, f64),
}

/// Mesh a [`SceneShape`] into a [`TexturedMesh`] in world space.
///
/// The shape is tessellated at `deflection` and its vertices are pushed
/// through the shape's world transform. Texture coordinates come from
/// `shape.uv` when it is present and matches the vertex count (this is how a
/// node imported by [`crate::rwmesh`] — whose [`crate::rwmesh::MeshNode`]
/// stores `uv` — keeps its authored UVs); otherwise a planar projection is
/// synthesized onto the dominant axis plane so every vertex lands in
/// `[0, 1]²` (see [`synthesize_planar_uv`]).
pub fn textured_mesh_from_scene_shape(shape: &SceneShape, deflection: f64) -> TexturedMesh {
    let mesh = crate::shape_mesh::mesh_shape(&shape.shape, deflection);
    let vertices: Vec<GpPnt> = mesh
        .vertices
        .iter()
        .map(|v| v.transformed(&shape.transform))
        .collect();
    let triangles: Vec<(usize, usize, usize)> = mesh
        .triangles
        .iter()
        .map(|t| (t.n0, t.n1, t.n2))
        .collect();
    let uv = match &shape.uv {
        Some(u) if u.len() == vertices.len() && !u.is_empty() => u.clone(),
        _ => synthesize_planar_uv(&vertices),
    };
    let color = shape.color.unwrap_or_else(default_color);
    TexturedMesh { vertices, triangles, uv, texture: None, color }
}

/// Synthesize planar-projection texture coordinates for `verts`.
///
/// The dominant axis (the world axis with the largest bounding-box extent;
/// ties resolve in favour of `Z`, then `Y`, then `X`) is dropped and the other
/// two coordinates are normalized into `[0, 1]²` against the mesh's bounding
/// box. This is the classic "project onto the largest face" UV mapping used
/// when a mesh carries no authored texture coordinates.
pub(super) fn synthesize_planar_uv(verts: &[GpPnt]) -> Vec<GpPnt2d> {
    if verts.is_empty() {
        return Vec::new();
    }
    let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut zmin, mut zmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in verts {
        xmin = xmin.min(v.x());
        xmax = xmax.max(v.x());
        ymin = ymin.min(v.y());
        ymax = ymax.max(v.y());
        zmin = zmin.min(v.z());
        zmax = zmax.max(v.z());
    }
    let ex = xmax - xmin;
    let ey = ymax - ymin;
    let ez = zmax - zmin;
    // Dominant axis = the largest extent; tie-break Z, then Y, then X.
    // (u_axis, v_axis) are the world coordinate indices kept as u and v.
    let (u_axis, v_axis) = if ez >= ey && ez >= ex {
        (0, 1) // project onto the XY plane
    } else if ey >= ex {
        (0, 2) // project onto the XZ plane
    } else {
        (1, 2) // project onto the YZ plane
    };
    let coord = |p: &GpPnt, axis: usize| -> f64 {
        match axis {
            0 => p.x(),
            1 => p.y(),
            _ => p.z(),
        }
    };
    let (umin, umax) = (coord(&GpPnt::new(xmin, ymin, zmin), u_axis), coord(&GpPnt::new(xmax, ymax, zmax), u_axis));
    let (vmin, vmax) = (coord(&GpPnt::new(xmin, ymin, zmin), v_axis), coord(&GpPnt::new(xmax, ymax, zmax), v_axis));
    let ur = (umax - umin).max(1e-12);
    let vr = (vmax - vmin).max(1e-12);
    verts
        .iter()
        .map(|p| {
            let u = coord(p, u_axis);
            let v = coord(p, v_axis);
            GpPnt2d::new(((u - umin) / ur).clamp(0.0, 1.0), ((v - vmin) / vr).clamp(0.0, 1.0))
        })
        .collect()
}

/// Sample a texture at a normalized `(u, v)` coordinate in `[0, 1]²`.
///
/// The coordinate is clamped to the unit square and the texel at
/// `floor(u·width) × floor(v·height)` is returned (edge-clamped so the last
/// texel covers `u = 1`). A degenerate texture (zero area or empty pixels)
/// yields white.
pub(super) fn sample_texture(tex: &Texture, u: f64, v: f64) -> [u8; 3] {
    if tex.width == 0 || tex.height == 0 || tex.pixels.is_empty() {
        return [255, 255, 255];
    }
    let tx = ((u.clamp(0.0, 1.0) * tex.width as f64).floor() as usize).min(tex.width - 1);
    let ty = ((v.clamp(0.0, 1.0) * tex.height as f64).floor() as usize).min(tex.height - 1);
    tex.pixels[ty * tex.width + tx]
}

/// Flatten every [`TexturedMesh`] into one triangle buffer with `(shape,
/// triangle)` tags, ready for [`build_tri_bvh`] + [`bvh_ray_cast`].
pub(super) fn flatten_textured_meshes(meshes: &[TexturedMesh]) -> (Vec<(GpPnt, GpPnt, GpPnt)>, Vec<(usize, usize)>) {
    let mut tri_pts = Vec::new();
    let mut tags = Vec::new();
    for (si, m) in meshes.iter().enumerate() {
        for (ti, &(n0, n1, n2)) in m.triangles.iter().enumerate() {
            tri_pts.push((m.vertices[n0], m.vertices[n1], m.vertices[n2]));
            tags.push((si, ti));
        }
    }
    (tri_pts, tags)
}

/// Ray-cast the scene into a [`Raster`] with UV texture mapping.
///
/// A single [`occt_core::bvh::builder_tri::TriBvh`] is built over every
/// shape's world triangles; each triangle keeps its owning shape and triangle
/// index so the correct per-vertex UVs are used. One ray is cast per pixel;
/// the nearest BVH hit acts as the depth buffer. At a hit the per-vertex UVs
/// are barycentrically interpolated, the texture is sampled at the resulting
/// `(u, v)`, and the texel color is scaled by a simple Lambert term
/// (`0.35 + 0.65·max(0, n·l)`) with the light pointing from the hit toward the
/// camera — front-facing surfaces appear fully lit. Shapes whose mesh has no
/// UV (or that have no geometry) fall back to the mesh's base `color` under
/// the same Lambert term. Missed pixels get a dark background.
pub fn render_textured_raster(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    texture: &Texture,
) -> Raster {
    let meshes: Vec<TexturedMesh> = scene
        .shapes
        .iter()
        .map(|ss| textured_mesh_from_scene_shape(ss, deflection))
        .collect();
    let (tri_pts, tags) = flatten_textured_meshes(&meshes);
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
                let (a, b, c) = (m.vertices[n0], m.vertices[n1], m.vertices[n2]);
                let hit = GpPnt::new(
                    origin.x() + dir.x() * t,
                    origin.y() + dir.y() * t,
                    origin.z() + dir.z() * t,
                );
                let normal = GpVec::from_pnts(&a, &b).crossed(&GpVec::from_pnts(&a, &c)).normalized();
                let light = GpVec::from_pnts(&hit, &cam.eye).normalized();
                let shade = 0.35 + 0.65 * normal.dot(&light).max(0.0);
                if m.uv.len() == m.vertices.len() && !m.uv.is_empty() {
                    let (la, lb, lc) = triangle_barycentric(&a, &b, &c, &hit);
                    let u = la * m.uv[n0].x() + lb * m.uv[n1].x() + lc * m.uv[n2].x();
                    let v = la * m.uv[n0].y() + lb * m.uv[n1].y() + lc * m.uv[n2].y();
                    let t = sample_texture(texture, u, v);
                    (
                        t[0] as f64 / 255.0 * shade,
                        t[1] as f64 / 255.0 * shade,
                        t[2] as f64 / 255.0 * shade,
                    )
                } else {
                    (m.color.0 * shade, m.color.1 * shade, m.color.2 * shade)
                }
            }
            None => (0.03, 0.03, 0.05),
        }
    });
    raster
}

/// Render the scene with UV texture mapping to binary PPM bytes.
///
/// Equivalent to [`render_textured_raster`] followed by
/// [`Raster::to_ppm_bytes`]. Shapes with per-vertex UV (kept from a
/// [`crate::rwmesh`] import or synthesized as a planar projection) sample
/// `texture`; shapes without UV fall back to their material color. This is
/// the analogue of an `V3d_View::Dump` with a texture-mapped display.
pub fn render_textured_ppm(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    texture: &Texture,
) -> Vec<u8> {
    render_textured_raster(scene, cam, width, height, deflection, texture).to_ppm_bytes()
}

// ---------------------------------------------------------------------------
// Bitmap font labels
// ---------------------------------------------------------------------------

/// A tiny 5×7 bitmap font.
///
/// Each glyph is a 7-row × 5-column pixel pattern: `[[u8; 5]; 7]` with
/// `1` for a foreground pixel and `0` for a transparent pixel. The built-in
/// [`BitmapFont::new`] font covers `A–Z`, `0–9`, space and `':'` — enough for
/// short scene labels and viewport captions. This is the analogue of the
/// 5×7 glyph tables OCCT bundles for its `V3d_Viewer` / `AIS` trihedron text.
#[derive(Debug, Clone)]
pub struct BitmapFont {
    /// Glyph patterns, keyed by character.
    pub glyphs: std::collections::HashMap<char, [[u8; 5]; 7]>,
}
