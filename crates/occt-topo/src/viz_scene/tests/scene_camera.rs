use super::*;



    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::render_svg::svg_polygon_count;

    pub(super) fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0
    }

    pub(super) fn translate(x: f64, y: f64, z: f64) -> GpTrsf {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(x, y, z));
        t
    }

    pub(super) fn box_scene() -> VizScene {
        let mut s = VizScene::new();
        s.add(SceneShape::new(unit_box()));
        s
    }

    /// Screen-space bbox width of shape `idx` under `cam`.
    pub(super) fn projected_screen_width(
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

    pub(super) fn sphere_scene() -> VizScene {
        let mut s = VizScene::new();
        s.add(SceneShape::new(BRepPrimSphere::make_sphere(1.0).solid.0));
        s
    }

    /// Colors of the pixels in a square neighbourhood around `(cx, cy)`.
    pub(super) fn sample_region(raster: &Raster, cx: f64, cy: f64, r: i32) -> Vec<(f64, f64, f64)> {
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
    pub(super) fn count_distinct(raster: &Raster) -> usize {
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

    // -- Texture / font / multi-view tests -----------------------------------

    #[test]
    fn checkerboard_pixels_distinct() {
        let tex = checkerboard_texture(4, 4, 2, (255, 0, 0), (0, 0, 255));
        assert_eq!(tex.width, 4);
        assert_eq!(tex.height, 4);
        let has_red = tex.pixels.contains(&[255, 0, 0]);
        let has_blue = tex.pixels.contains(&[0, 0, 255]);
        assert!(has_red, "expected red cells");
        assert!(has_blue, "expected blue cells");
        assert!(has_red && has_blue, "checkerboard should mix both colors");
    }

    #[test]
    fn texture_dimensions() {
        let tex = checkerboard_texture(8, 8, 4, (10, 20, 30), (40, 50, 60));
        assert_eq!(tex.width, 8);
        assert_eq!(tex.height, 8);
        assert_eq!(tex.pixels.len(), 64);
    }

    #[test]
    fn textured_render_not_blank() {
        let scene = box_scene();
        let cam = Camera::default();
        let tex = checkerboard_texture(16, 16, 4, (220, 40, 40), (40, 40, 220));
        let textured = render_textured_ppm(&scene, &cam, 96, 72, 0.25, &tex);
        let solid = render_scene_ppm_shaded(&scene, &cam, 96, 72, 0.25, &RenderSettings::default());
        assert_ne!(textured, solid, "textured output should differ from the solid-color render");
        let body = &textured[13..];
        let near = |p: &[u8], c: (u8, u8, u8)| {
            (p[0] as i16 - c.0 as i16).abs() <= 6
                && (p[1] as i16 - c.1 as i16).abs() <= 6
                && (p[2] as i16 - c.2 as i16).abs() <= 6
        };
        let has_c1 = body.chunks_exact(3).any(|p| near(p, (220, 40, 40)));
        let has_c2 = body.chunks_exact(3).any(|p| near(p, (40, 40, 220)));
        assert!(has_c1, "expected checkerboard color 1 pixels");
        assert!(has_c2, "expected checkerboard color 2 pixels");
    }

    #[test]
    fn uv_synthesis_dominant_axis() {
        let ss = SceneShape::new(unit_box()); // [0,1]^3 — equal extents, dominant Z
        let tm = textured_mesh_from_scene_shape(&ss, 0.25);
        assert_eq!(tm.uv.len(), tm.vertices.len());
        assert!(!tm.uv.is_empty());
        // Dominant axis is Z, so UV = (x, y) normalized to [0,1]².
        for (p, uv) in tm.vertices.iter().zip(&tm.uv) {
            assert!(uv.x() >= 0.0 && uv.x() <= 1.0 && uv.y() >= 0.0 && uv.y() <= 1.0, "uv {uv:?}");
            assert!((uv.x() - p.x()).abs() < 1e-6, "u should be x: {} vs {}", uv.x(), p.x());
            assert!((uv.y() - p.y()).abs() < 1e-6, "v should be y: {} vs {}", uv.y(), p.y());
        }
    }

    #[test]
    fn bitmap_font_has_glyphs() {
        let font = BitmapFont::default();
        for c in ['A', 'B', 'Z', '0', '9', ' '] {
            assert!(font.has_glyph(c), "missing glyph {c:?}");
        }
        let a = font.glyph('A').expect("A glyph present");
        assert_eq!(a.len(), 7);
        assert!(a.iter().all(|row| row.len() == 5));
        assert!(a[0].contains(&1), "top row of 'A' should have lit pixels");
    }

    #[test]
    fn text_raster_size_scales() {
        let font = BitmapFont::default();
        assert_eq!(text_raster_size("A", &font, 1), (5, 7));
        assert_eq!(text_raster_size("AB", &font, 1), (10, 7));
        assert_eq!(text_raster_size("A", &font, 2), (10, 14));
    }

    #[test]
    fn overlay_text_pixels() {
        let font = BitmapFont::default();
        let (tw, th) = text_raster_size("A", &font, 1);
        assert_eq!((tw, th), (5, 7));
        let mut body = vec![0u8; 40 * 30 * 3]; // raw RGB body, no header
        overlay_text(&mut body, 40, 30, "A", &font, 1, 2, 2, (255, 255, 255));
        let has_white = body.chunks_exact(3).enumerate().any(|(i, p)| {
            let px = i % 40;
            let py = i / 40;
            p == [255, 255, 255] && px >= 2 && px < 2 + 5 && py >= 2 && py < 2 + 7
        });
        assert!(has_white, "expected white text pixels at the overlay position");
    }

    #[test]
    fn render_view_grid_tiles() {
        let mut scene_a = VizScene::new();
        scene_a.add(SceneShape::new(unit_box()));
        let mut scene_b = VizScene::new();
        let mut box_b = SceneShape::new(unit_box());
        box_b.transform = translate(1.5, 0.0, 0.0);
        scene_b.add(box_b);
        let cam = Camera::default();
        let ppm = render_view_grid(&[&scene_a, &scene_b], &cam, 40, 30, ViewLayout { cols: 2, rows: 1 }, 0.25);
        assert!(ppm.starts_with(b"P6\n80 30\n255\n"), "expected 2×1 grid of 40×30 tiles");
        let body = &ppm[13..];
        let bg = [
            (0.03 * 255.0) as u8,
            (0.03 * 255.0) as u8,
            (0.05 * 255.0) as u8,
        ];
        // The full image is 80 wide (2×40 columns); the left tile occupies
        // columns 0–39 and the right tile columns 40–79 of every row.
        let grid_w = 80usize;
        let left_has = body.chunks_exact(3).enumerate().any(|(i, p)| p != bg && (i % grid_w) < 40);
        let right_has = body.chunks_exact(3).enumerate().any(|(i, p)| p != bg && (i % grid_w) >= 40);
        assert!(left_has, "left tile should have content");
        assert!(right_has, "right tile should have content");
    }

    #[test]
    fn label_overlay() {
        let scene = box_scene();
        let font = BitmapFont::default();
        // The default scene produces no exact-white pixels (background is dark
        // blue-gray and the material's specular is faint), so white label text
        // is unambiguous in the top rows.
        let ppm = render_scene_with_label(&scene, &Camera::default(), 64, 48, 0.25, "HI", &font, 2);
        let body = &ppm[13..];
        let w = 64usize;
        let has_label = body.chunks_exact(3).enumerate().any(|(i, p)| p == [255, 255, 255] && i / w < 12);
        assert!(has_label, "expected white label pixels in the top rows");
    }

    #[test]
    fn textured_sphere_looks_textured() {
        let scene = sphere_scene();
        let tex = checkerboard_texture(16, 16, 4, (240, 60, 60), (60, 60, 240));
        let ppm = render_textured_ppm(&scene, &Camera::default(), 64, 48, 0.25, &tex);
        let body = &ppm[13..];
        let has_red = body.chunks_exact(3).any(|p| p[0] as u16 > p[2] as u16 + 40);
        let has_blue = body.chunks_exact(3).any(|p| p[2] as u16 > p[0] as u16 + 40);
        assert!(has_red && has_blue, "sphere should show both checkerboard colors (red {has_red}, blue {has_blue})");
    }

    #[test]
    fn solid_texture_uniform() {
        let tex = solid_texture(4, 5, (7, 8, 9));
        assert_eq!(tex.width, 4);
        assert_eq!(tex.height, 5);
        assert_eq!(tex.pixels.len(), 20);
        assert!(tex.pixels.iter().all(|p| p == &[7, 8, 9]));
    }
