//! Color manipulation utilities. Source: `Quantity_Color` advanced ops.
use super::{Color, ColorRGBA};

impl Color {
    /// Convert to HSV color space. Returns (hue[0,360), sat[0,1], val[0,1]).
    pub fn to_hsv(&self) -> (f64, f64, f64) {
        let r = self.r as f64; let g = self.g as f64; let b = self.b as f64;
        let max = r.max(g).max(b);
        let min = r.min(g).min(b);
        let d = max - min;
        let h = if d == 0.0 { 0.0 }
            else if max == r { 60.0 * ((g - b) / d).rem_euclid(6.0) }
            else if max == g { 60.0 * ((b - r) / d + 2.0) }
            else { 60.0 * ((r - g) / d + 4.0) };
        let s = if max == 0.0 { 0.0 } else { d / max };
        (h, s, max)
    }

    /// Construct from HSV. h in [0,360), s in [0,1], v in [0,1].
    pub fn from_hsv(h: f64, s: f64, v: f64) -> Self {
        let c = v * s;
        let hp = h / 60.0;
        let x = c * (1.0 - (hp.rem_euclid(2.0) - 1.0).abs());
        let m = v - c;
        let (r, g, b) = match hp.floor() as i32 {
            0 => (c, x, 0.0),
            1 => (x, c, 0.0),
            2 => (0.0, c, x),
            3 => (0.0, x, c),
            4 => (x, 0.0, c),
            _ => (c, 0.0, x),
        };
        Self { r: (r + m) as f32, g: (g + m) as f32, b: (b + m) as f32 }
    }

    /// Linear interpolation between two colors.
    pub fn lerp(a: &Color, b: &Color, t: f32) -> Self {
        Self {
            r: a.r + (b.r - a.r) * t,
            g: a.g + (b.g - a.g) * t,
            b: a.b + (b.b - a.b) * t,
        }
    }

    /// Relative luminance (Rec. 709).
    pub fn luminance(&self) -> f32 {
        0.2126 * self.r + 0.7152 * self.g + 0.0722 * self.b
    }

    /// Blend with alpha over background.
    pub fn blend_over(&self, bg: &Color, alpha: f32) -> Self {
        Self {
            r: self.r * alpha + bg.r * (1.0 - alpha),
            g: self.g * alpha + bg.g * (1.0 - alpha),
            b: self.b * alpha + bg.b * (1.0 - alpha),
        }
    }

    /// Convert to RGBA with given alpha.
    pub fn to_rgba(&self, alpha: f32) -> ColorRGBA {
        ColorRGBA { r: self.r, g: self.g, b: self.b, a: alpha }
    }

    /// Parse "r,g,b" (0-255) or "#RRGGBB".
    pub fn parse(s: &str) -> Option<Self> {
        let s = s.trim();
        if let Some(hex) = s.strip_prefix('#') {
            if hex.len() == 6 {
                let r = u8::from_str_radix(&hex[0..2], 16).ok()?;
                let g = u8::from_str_radix(&hex[2..4], 16).ok()?;
                let b = u8::from_str_radix(&hex[4..6], 16).ok()?;
                return Some(Self { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0 });
            }
            return None;
        }
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() == 3 {
            let r: f64 = parts[0].trim().parse().ok()?;
            let g: f64 = parts[1].trim().parse().ok()?;
            let b: f64 = parts[2].trim().parse().ok()?;
            return Some(Self { r: r as f32 / 255.0, g: g as f32 / 255.0, b: b as f32 / 255.0 });
        }
        None
    }

    /// Convert to hex string "#RRGGBB".
    pub fn to_hex_string(&self) -> String {
        let (r, g, b) = self.to_rgb_u8();
        format!("#{r:02x}{g:02x}{b:02x}")
    }
}

impl ColorRGBA {
    /// Premultiplied alpha representation.
    pub fn premultiplied(&self) -> Self {
        Self { r: self.r * self.a, g: self.g * self.a, b: self.b * self.a, a: self.a }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hsv_roundtrip() {
        let c = Color::RED;
        let (h, s, v) = c.to_hsv();
        assert!((h - 0.0).abs() < 1e-6);
        let back = Color::from_hsv(h, s, v);
        assert!((back.r - c.r).abs() < 1e-3);
    }

    #[test]
    fn hsv_green() {
        let c = Color::GREEN;
        let (h, _, _) = c.to_hsv();
        assert!((h - 120.0).abs() < 1e-6);
    }

    #[test]
    fn parse_hex() {
        let c = Color::parse("#ff0000").unwrap();
        assert!((c.r - 1.0).abs() < 1e-6);
        assert_eq!(c.to_hex_string(), "#ff0000");
    }

    #[test]
    fn parse_rgb() {
        let c = Color::parse("255, 0, 0").unwrap();
        assert!((c.r - 1.0).abs() < 1e-6);
    }
}
