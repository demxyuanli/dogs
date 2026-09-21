use super::*;

use super::*;

    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::render_svg::svg_polygon_count;

    // -- Phase 12: picking / selection highlight / font styles ----------------

    #[test]
    fn pick_two_shape_scene_hits_correct() {
        // A box (index 0) and a sphere (index 1), side by side so a pick
        // unambiguously resolves to one of them.
        let mut scene = VizScene::new();
        let mut box_ss = SceneShape::new(unit_box());
        box_ss.transform = translate(-1.5, 0.0, 0.0); // x in [-1.5, -0.5]
        box_ss.material = Material::from_diffuse((0.8, 0.1, 0.1));
        scene.add(box_ss);
        let mut sphere_ss = SceneShape::new(BRepPrimSphere::make_sphere(0.5).solid.0);
        sphere_ss.transform = translate(1.5, 0.0, 0.0); // x in [1.0, 2.0]
        sphere_ss.material = Material::from_diffuse((0.1, 0.1, 0.8));
        scene.add(sphere_ss);
        let cam = Camera::default();
        // Box front-face centre and the sphere's front-most point toward the
        // camera (a head-on hit; grazing rays can slip between the flat
        // tessellation triangles of the coarse sphere mesh).
        let (bx, by, _) = project_point(&cam, GpPnt::new(-1.0, 0.5, 1.0), 200, 150).unwrap();
        let (sx, sy, _) = project_point(&cam, GpPnt::new(1.3565, 0.0, 0.479), 200, 150).unwrap();
        let (idx_b, hit_b, _t) = pick_shape(&scene, &cam, 200, 150, bx as usize, by as usize, 0.25).unwrap();
        assert_eq!(idx_b, 0, "box pick returns the box, got {idx_b}");
        assert!((hit_b.z() - 1.0).abs() < 1e-6, "box front face at z=1, got {}", hit_b.z());
        let (idx_s, hit_s, _t) = pick_shape(&scene, &cam, 200, 150, sx as usize, sy as usize, 0.25).unwrap();
        assert_eq!(idx_s, 1, "sphere pick returns the sphere, got {idx_s}");
        let dist_c = hit_s.distance(&GpPnt::new(1.5, 0.0, 0.0));
        assert!((dist_c - 0.5).abs() < 0.03, "hit sits on the sphere surface, dist {dist_c:.3}");
        assert!(hit_s.z() > 0.3, "front hemisphere of the sphere, got z={}", hit_s.z());
        // pick_point is the hit-point form of the same query.
        let p = pick_point(&scene, &cam, 200, 150, bx as usize, by as usize, 0.25).unwrap();
        assert!(p.distance(&hit_b) < 1e-9, "pick_point agrees with pick_shape");
    }

    #[test]
    fn render_selection_highlights_picked() {
        let scene = box_scene();
        let cam = Camera::default();
        let settings = RenderSettings::default();
        // Pick a pixel on the box's front face.
        let (sx, sy, _) = project_point(&cam, GpPnt::new(0.5, 0.5, 1.0), 96, 72).unwrap();
        let (px, py) = (sx as usize, sy as usize);
        let (idx, _hit, _t) = pick_shape(&scene, &cam, 96, 72, px, py, 0.25).unwrap();
        assert_eq!(idx, 0);
        // Base and highlighted renders use the same z-buffer pipeline, so the
        // only difference is the selection material.
        let base = render_scene_raster_zbuffer(&scene, &cam, 96, 72, 0.25, &settings);
        let sel = render_scene_with_selection(&scene, &cam, 96, 72, 0.25, &settings, &[idx]);
        assert!(sel.starts_with(b"P6\n96 72\n255\n"), "PPM header");
        let base_px = base.get_pixel(px, py);
        let i = 13 + (py * 96 + px) * 3;
        let sel_px = (sel[i] as f64 / 255.0, sel[i + 1] as f64 / 255.0, sel[i + 2] as f64 / 255.0);
        // The highlight is yellow-dominant while the base front face is a
        // neutral gray, so the red-minus-blue gap grows.
        let gap_base = base_px.0 - base_px.2;
        let gap_sel = sel_px.0 - sel_px.2;
        assert!(
            gap_sel > gap_base + 0.2,
            "highlight gap {gap_sel:.3} should beat base gap {gap_base:.3}"
        );
        let delta =
            (base_px.0 - sel_px.0).abs() + (base_px.1 - sel_px.1).abs() + (base_px.2 - sel_px.2).abs();
        assert!(delta > 0.2, "highlighted pixel should differ from base, delta {delta:.3}");
    }

    #[test]
    fn font_style_sizes() {
        let std = font_for_style(FontStyle::Standard);
        assert_eq!((std.glyph_w, std.glyph_h), (5, 7));
        assert_eq!(font_raster_size("A", &std, 1), (5, 7));
        assert_eq!(font_raster_size("AB", &std, 2), (20, 14));
        let large = font_for_style(FontStyle::Large);
        assert_eq!((large.glyph_w, large.glyph_h), (7, 9));
        assert_eq!(font_raster_size("A", &large, 1), (7, 9), "large font is 7×9");
        // Every glyph from the 5×7 set survives scaling, including symbols.
        for c in ['A', 'Z', '0', '9', ' ', '-', '!'] {
            assert!(large.has_glyph(c), "scaled font missing {c:?}");
        }
        assert!(large.glyph('A').unwrap().contains(&1), "scaled glyph keeps lit pixels");
    }

    #[test]
    fn font_charset_extended() {
        let font = font_for_style(FontStyle::Standard);
        for c in ['-', '.', ',', '!', '?', '+', '=', '/', '_', '(', ')', ';'] {
            assert!(font.has_glyph(c), "missing extended symbol {c:?}");
        }
        let px = font_to_pixels("!", &font, 1);
        assert!(px.contains(&[255, 255, 255]), "exclamation should have lit pixels");
    }

    #[test]
    fn overlay_text_rect_draws_box() {
        let font = font_for_style(FontStyle::Standard);
        let mut raster = Raster::new(40, 30);
        raster.clear((0.0, 0.0, 0.0));
        overlay_text_rect(&mut raster, "A", &font, 1, 10, 10, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        // Padding is 2 and "A" is 5×7, so the white box spans x∈[8,16], y∈[8,18].
        assert_eq!(raster.get_pixel(8, 8), (1.0, 1.0, 1.0), "box top-left corner filled");
        assert_eq!(raster.get_pixel(16, 8), (1.0, 1.0, 1.0), "box top-right corner filled");
        assert_eq!(raster.get_pixel(6, 8), (0.0, 0.0, 0.0), "outside the box stays untouched");
        // 'A' top row has a lit pixel at source column 1 → absolute (11, 10).
        assert_eq!(raster.get_pixel(11, 10), (0.0, 0.0, 0.0), "text pixel stamped black");
    }

    #[test]
    fn set_font_registers_style() {
        let mut registry = std::collections::HashMap::new();
        let custom = Font::new(3, 5);
        set_font(&mut registry, FontStyle::Standard, custom.clone());
        assert_eq!(registry.len(), 1);
        assert_eq!(registry[&FontStyle::Standard].glyph_w, 3);
        assert_eq!(registry[&FontStyle::Standard].glyph_h, 5);
    }

    #[test]
    fn overlay_text_font_large() {
        let font = font_for_style(FontStyle::Large);
        let mut body = vec![0u8; 30 * 20 * 3];
        overlay_text_font(&mut body, 30, 20, "A", &font, 1, 0, 0, (255, 255, 255));
        let has_white = body.chunks_exact(3).any(|p| p == [255, 255, 255]);
        assert!(has_white, "large font should render lit pixels");
    }
