use super::prelude::*;
use super::*;

/// Draw `text` at pixel position `(x, y)` onto a `width`×`height` RGB byte
/// buffer, tinted `color`.
///
/// The buffer holds `width * height * 3` raw RGB bytes (the PPM pixel body
/// without the header). If `ppm` is a complete P6 PPM image — `len ==
/// width * height * 3 + 13` and it starts with `P6\n` — the 13-byte header is
/// skipped automatically, so both raw bodies and full PPM buffers work.
/// Transparent (black) raster pixels are skipped, leaving the underlying image
/// untouched; glyph pixels are written as `color` at the scaled 5×7 positions.
pub fn overlay_text(
    ppm: &mut Vec<u8>,
    width: usize,
    height: usize,
    text: &str,
    font: &BitmapFont,
    scale: usize,
    x: usize,
    y: usize,
    color: (u8, u8, u8),
) {
    let pixels = text_to_pixels(text, font, scale);
    let (tw, th) = text_raster_size(text, font, scale);
    let header_len = if ppm.len() == width * height * 3 + 13 && ppm.starts_with(b"P6\n") {
        13
    } else {
        0
    };
    let body = &mut ppm[header_len..];
    for py in 0..th {
        for px in 0..tw {
            if pixels[py * tw + px] == [0, 0, 0] {
                continue;
            }
            let dx = x + px;
            let dy = y + py;
            if dx < width && dy < height {
                let i = (dy * width + dx) * 3;
                body[i] = color.0;
                body[i + 1] = color.1;
                body[i + 2] = color.2;
            }
        }
    }
}

/// Render the scene with a text label overlaid at the top-left corner.
///
/// The scene is shaded with [`render_scene_ppm_shaded`] using default
/// [`RenderSettings`], then `label` is drawn in white at `(4, 4)` with the
/// given `font` and `scale`. The result is a complete P6 PPM image — the
/// analogue of a `V3d_View::Dump` with the viewer's text caption enabled.
pub fn render_scene_with_label(
    scene: &VizScene,
    cam: &Camera,
    width: usize,
    height: usize,
    deflection: f64,
    label: &str,
    font: &BitmapFont,
    scale: usize,
) -> Vec<u8> {
    let settings = RenderSettings::default();
    let ppm = render_scene_ppm_shaded(scene, cam, width, height, deflection, &settings);
    let mut body = ppm[13..].to_vec();
    overlay_text(&mut body, width, height, label, font, scale, 4, 4, (255, 255, 255));
    let mut out = Vec::with_capacity(body.len() + 13);
    out.extend_from_slice(format!("P6\n{width} {height}\n255\n").as_bytes());
    out.extend_from_slice(&body);
    out
}

// ---------------------------------------------------------------------------
// Font styles (TKService full glyph texture set)
// ---------------------------------------------------------------------------

/// A size-generic bitmap font: each glyph is a fixed `glyph_w`×`glyph_h`
/// pattern of `1`/`0` bits.
///
/// Unlike [`BitmapFont`] (hard-coded 5×7 cells), a [`Font`] carries its own
/// glyph cell size, so one text rasterizer serves every style. This is the
/// analogue of OCCT's `Font_FTFont` text drawing, where one font object owns
/// its glyph cache and metrics and the renderer only asks "how big is a cell"
/// and "what bits are in this glyph".
#[derive(Debug, Clone)]
pub struct Font {
    /// Glyph cell width in pixels.
    pub glyph_w: usize,
    /// Glyph cell height in pixels.
    pub glyph_h: usize,
    /// Row-major glyph bit patterns, `glyph_w * glyph_h` `1`/`0` per character.
    pub glyphs: std::collections::HashMap<char, Vec<u8>>,
}

impl Font {
    /// An empty font with a given glyph cell size.
    pub fn new(glyph_w: usize, glyph_h: usize) -> Self {
        Self {
            glyph_w: glyph_w.max(1),
            glyph_h: glyph_h.max(1),
            glyphs: std::collections::HashMap::new(),
        }
    }

    /// Insert a `glyph_w × glyph_h` glyph pattern for `c`.
    ///
    /// `pattern` holds one `1`/`0` bit per cell, row-major; entries shorter
    /// than `glyph_w * glyph_h` are zero-padded, longer entries truncated.
    pub fn insert(&mut self, c: char, pattern: &[u8]) {
        let mut bits = vec![0u8; self.glyph_w * self.glyph_h];
        let n = pattern.len().min(bits.len());
        bits[..n].copy_from_slice(&pattern[..n]);
        self.glyphs.insert(c, bits);
    }

    /// The row-major bit pattern for `c`, or `None` when the font has no glyph.
    pub fn glyph(&self, c: char) -> Option<&[u8]> {
        self.glyphs.get(&c).map(|g| g.as_slice())
    }

    /// `true` when the font has a glyph for `c`.
    pub fn has_glyph(&self, c: char) -> bool {
        self.glyphs.contains_key(&c)
    }

    /// A 5×7 [`Font`] built from the built-in [`BitmapFont`] glyph set.
    pub fn from_bitmap(bitmap: &BitmapFont) -> Self {
        let mut glyphs = std::collections::HashMap::new();
        for (&c, g) in &bitmap.glyphs {
            let bits: Vec<u8> = g.iter().flatten().copied().collect();
            glyphs.insert(c, bits);
        }
        Self { glyph_w: 5, glyph_h: 7, glyphs }
    }

    /// A copy of this font whose glyphs are scaled to `glyph_w`×`glyph_h`
    /// cells by nearest-neighbour sampling.
    ///
    /// Used to derive larger (or smaller) styles from the same master glyph
    /// set without re-authoring every pattern.
    pub fn scaled(&self, glyph_w: usize, glyph_h: usize) -> Self {
        let gw = glyph_w.max(1);
        let gh = glyph_h.max(1);
        let mut glyphs = std::collections::HashMap::new();
        for (&c, bits) in &self.glyphs {
            let mut out = vec![0u8; gw * gh];
            for ty in 0..gh {
                let sy = ty * self.glyph_h / gh;
                for tx in 0..gw {
                    let sx = tx * self.glyph_w / gw;
                    if bits[sy * self.glyph_w + sx] != 0 {
                        out[ty * gw + tx] = 1;
                    }
                }
            }
            glyphs.insert(c, out);
        }
        Self { glyph_w: gw, glyph_h: gh, glyphs }
    }
}

/// A selectable text style, the analogue of OCCT's `Font_FontAspect`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FontStyle {
    /// The built-in 5×7 bitmap font.
    Standard,
    /// A 7×9 variant scaled from the 5×7 glyph set.
    Large,
}

/// The [`Font`] for a [`FontStyle`].
///
/// [`FontStyle::Standard`] is the 5×7 built-in glyph set (A–Z, 0–9 and the
/// common symbols in [`BITMAP_GLYPHS`]); [`FontStyle::Large`] is the same set
/// scaled to 7×9 cells.
pub fn font_for_style(style: FontStyle) -> Font {
    let standard = Font::from_bitmap(&BitmapFont::default());
    match style {
        FontStyle::Standard => standard,
        FontStyle::Large => standard.scaled(7, 9),
    }
}

/// Register `font` as the glyph source for `style` in a mutable registry.
///
/// A [`std::collections::HashMap`] keyed by [`FontStyle`] lets an application
/// keep one font per style (e.g. a themed viewport caption font) and switch at
/// draw time. Prefer [`font_for_style`] when the built-in styles are enough.
pub fn set_font(
    registry: &mut std::collections::HashMap<FontStyle, Font>,
    style: FontStyle,
    font: Font,
) {
    registry.insert(style, font);
}

/// The raster size of `text` in a generic [`Font`].
///
/// Each glyph is `font.glyph_w` wide and `font.glyph_h` tall, scaled by
/// `scale`; the width is `chars × glyph_w × scale` (no inter-glyph padding).
pub fn font_raster_size(text: &str, font: &Font, scale: usize) -> (usize, usize) {
    let scale = scale.max(1);
    (text.chars().count() * font.glyph_w * scale, font.glyph_h * scale)
}

/// Rasterize `text` in a generic [`Font`] into RGB pixels.
///
/// The returned buffer is `width × height × 3` bytes (row-major) where
/// `(width, height) = font_raster_size(text, font, scale)`. Foreground pixels
/// are white `[255, 255, 255]`; background pixels are black and count as
/// transparent for [`overlay_text_font`] / [`overlay_text_rect`]. Characters
/// without a glyph in `font` are skipped.
pub fn font_to_pixels(text: &str, font: &Font, scale: usize) -> Vec<[u8; 3]> {
    let scale = scale.max(1);
    let (w, h) = font_raster_size(text, font, scale);
    let mut pixels = vec![[0u8, 0, 0]; w * h];
    for (ci, ch) in text.chars().enumerate() {
        let Some(glyph) = font.glyph(ch) else { continue };
        for (ri, row) in glyph.chunks(font.glyph_w).enumerate() {
            for (rx, &on) in row.iter().enumerate() {
                if on == 0 {
                    continue;
                }
                for sy in 0..scale {
                    for sx in 0..scale {
                        let px = ci * font.glyph_w * scale + rx * scale + sx;
                        let py = ri * scale + sy;
                        if px < w && py < h {
                            pixels[py * w + px] = [255, 255, 255];
                        }
                    }
                }
            }
        }
    }
    pixels
}

/// Draw `text` in a generic [`Font`] at pixel position `(x, y)` onto a raw RGB
/// byte buffer, tinted `color`.
///
/// Mirrors [`overlay_text`] but for a size-carrying [`Font`], so the same call
/// serves every glyph cell size. Like [`overlay_text`], the buffer may be a
/// raw RGB body or a complete P6 PPM image (the 13-byte header is skipped
/// automatically); transparent (black) raster pixels are skipped.
pub fn overlay_text_font(
    ppm: &mut Vec<u8>,
    width: usize,
    height: usize,
    text: &str,
    font: &Font,
    scale: usize,
    x: usize,
    y: usize,
    color: (u8, u8, u8),
) {
    let pixels = font_to_pixels(text, font, scale);
    let (tw, th) = font_raster_size(text, font, scale);
    let header_len = if ppm.len() == width * height * 3 + 13 && ppm.starts_with(b"P6\n") {
        13
    } else {
        0
    };
    let body = &mut ppm[header_len..];
    for py in 0..th {
        for px in 0..tw {
            if pixels[py * tw + px] == [0, 0, 0] {
                continue;
            }
            let dx = x + px;
            let dy = y + py;
            if dx < width && dy < height {
                let i = (dy * width + dx) * 3;
                body[i] = color.0;
                body[i + 1] = color.1;
                body[i + 2] = color.2;
            }
        }
    }
}

/// Draw `text` with a filled background box onto a [`Raster`], tinted `color`.
///
/// A `PADDING`-pixel solid box of `background` is stamped first, then the text
/// glyphs (in `font` at `scale`) are stamped in `color` on top; both are
/// clipped to the raster. `(x, y)` is the top-left of the *text*; the box
/// extends `PADDING` pixels around it. Working on a [`Raster`] makes the
/// overlay composable with [`render_view_grid`]'s tile buffers so captions can
/// be stamped per viewport — the classic "caption chip" of a multi-view
/// `V3d_Viewer` window.
pub fn overlay_text_rect(
    raster: &mut Raster,
    text: &str,
    font: &Font,
    scale: usize,
    x: usize,
    y: usize,
    color: (f64, f64, f64),
    background: (f64, f64, f64),
) {
    pub(super) const PADDING: usize = 2;
    let (tw, th) = font_raster_size(text, font, scale);
    for by in 0..(th + PADDING * 2) {
        for bx in 0..(tw + PADDING * 2) {
            let dx = x as i64 + bx as i64 - PADDING as i64;
            let dy = y as i64 + by as i64 - PADDING as i64;
            if dx >= 0 && dy >= 0 && (dx as usize) < raster.width && (dy as usize) < raster.height {
                raster.set_pixel(dx as usize, dy as usize, background);
            }
        }
    }
    let pixels = font_to_pixels(text, font, scale);
    for py in 0..th {
        for px in 0..tw {
            if pixels[py * tw + px] == [0, 0, 0] {
                continue;
            }
            let dx = x as i64 + px as i64;
            let dy = y as i64 + py as i64;
            if dx >= 0 && dy >= 0 && (dx as usize) < raster.width && (dy as usize) < raster.height {
                raster.set_pixel(dx as usize, dy as usize, color);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Multi-view layout
// ---------------------------------------------------------------------------

/// A viewport grid for [`render_view_grid`].
///
/// `cols` × `rows` cells, each cell rendering one scene into a
/// `width` × `height` tile; the full output is `width*cols` × `height*rows`.
/// Mirrors the `V3d_Viewer` multi-view window layout where several `V3d_View`s
/// (front, top, side, isometric, …) share one screen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ViewLayout {
    /// Number of viewport columns.
    pub cols: usize,
    /// Number of viewport rows.
    pub rows: usize,
}

/// Render several scenes into a tiled multi-view grid as PPM bytes.
///
/// Every scene is rendered with the *same* [`Camera`] into a
/// `width` × `height` tile (the classic four-view viewer renders each viewport
/// with a different camera, but a shared camera keeps the grid useful for
/// comparing scenes side by side); the tiles are laid out row-major per
/// [`ViewLayout`] and the full image is returned as a P6 PPM buffer of size
/// `width*cols` × `height*rows`. Cells with no scene (when `scenes` is shorter
/// than the grid) stay at the background color.
pub fn render_view_grid(
    scenes: &[&VizScene],
    cam: &Camera,
    width: usize,
    height: usize,
    layout: ViewLayout,
    deflection: f64,
) -> Vec<u8> {
    let cols = layout.cols.max(1);
    let rows = layout.rows.max(1);
    let settings = RenderSettings::default();
    let mut out = Raster::new(width * cols, height * rows);
    out.clear(settings.background);
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            let Some(scene) = scenes.get(i) else { continue };
            let tile = render_scene_raster_zbuffer(scene, cam, width, height, deflection, &settings);
            for ty in 0..height {
                for tx in 0..width {
                    out.set_pixel(c * width + tx, r * height + ty, tile.get_pixel(tx, ty));
                }
            }
        }
    }
    out.to_ppm_bytes()
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Camera-space `(x, y)` in `[-1, 1]²` for a screen point `(sx, sy)` with the
/// origin at the top-left of a `width`×`height` image (y flipped).
pub(super) fn screen_to_ndc(width: usize, height: usize, sx: f64, sy: f64) -> (f64, f64) {
    let x_ndc = 2.0 * sx / width as f64 - 1.0;
    let y_ndc = 1.0 - 2.0 * sy / height as f64;
    (x_ndc, y_ndc)
}

/// Perspective scale `f = 1/tan(fov/2)`.
pub(super) fn perspective_scale(cam: &Camera) -> f64 {
    (cam.fov_deg.to_radians() * 0.5).tan().recip()
}

/// Orthographic viewport half-height in world units.
///
/// Matches the perspective frustum width at the target distance so switching
/// projection keeps the on-screen size approximately constant.
pub(super) fn ortho_scale(cam: &Camera) -> f64 {
    let s = cam.eye.distance(&cam.target) * (cam.fov_deg.to_radians() * 0.5).tan();
    if s.abs() < 1e-12 { 1.0 } else { s }
}

/// The orthonormal camera basis `(right, up2, forward)`.
///
/// Handles a degenerate `forward` (eye == target) by defaulting to `−Z`, and
/// an `up` parallel to `forward` by picking a reference axis not parallel to
/// the forward direction.
pub(super) fn camera_basis(cam: &Camera) -> (GpVec, GpVec, GpVec) {
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
pub(super) fn rotate_vec(v: &GpVec, axis: &GpVec, angle: f64) -> GpVec {
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
