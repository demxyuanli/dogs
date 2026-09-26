
use super::*;

impl Default for BitmapFont {
    fn default() -> Self {
        Self::new()
    }
}

impl BitmapFont {
    /// The built-in font: `A–Z`, `0–9`, space and `':'` at 5×7.
    pub fn new() -> Self {
        let mut glyphs = std::collections::HashMap::new();
        for &(c, g) in BITMAP_GLYPHS {
            glyphs.insert(c, g);
        }
        Self { glyphs }
    }

    /// The 5×7 pattern for `c`, or `None` when the font has no glyph for it.
    pub fn glyph(&self, c: char) -> Option<&[[u8; 5]; 7]> {
        self.glyphs.get(&c)
    }

    /// `true` when the font has a glyph for `c`.
    pub fn has_glyph(&self, c: char) -> bool {
        self.glyphs.contains_key(&c)
    }
}

/// The raster size of `text` in pixels.
///
/// Each glyph is 5 wide and 7 tall, scaled by `scale`; the width is
/// `chars × 5 × scale` (no inter-glyph padding) and the height `7 × scale`.
pub fn text_raster_size(text: &str, _font: &BitmapFont, scale: usize) -> (usize, usize) {
    let scale = scale.max(1);
    let n = text.chars().count();
    (n * 5 * scale, 7 * scale)
}

/// Rasterize `text` into an RGB pixel array.
///
/// The returned buffer is `width × height × 3` bytes (row-major) where
/// `(width, height) = text_raster_size(text, font, scale)`. Foreground pixels
/// are white `[255, 255, 255]`; background pixels are black `[0, 0, 0]` and
/// count as transparent for [`overlay_text`]. Characters without a glyph in
/// `font` are skipped.
pub fn text_to_pixels(text: &str, font: &BitmapFont, scale: usize) -> Vec<[u8; 3]> {
    let scale = scale.max(1);
    let (w, h) = text_raster_size(text, font, scale);
    let mut pixels = vec![[0u8, 0, 0]; w * h];
    for (ci, ch) in text.chars().enumerate() {
        let Some(glyph) = font.glyph(ch) else { continue };
        for (ry, row) in glyph.iter().enumerate() {
            for (rx, &on) in row.iter().enumerate() {
                if on == 0 {
                    continue;
                }
                for sy in 0..scale {
                    for sx in 0..scale {
                        let px = ci * 5 * scale + rx * scale + sx;
                        let py = ry * scale + sy;
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
