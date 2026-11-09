//! Physical quantities and named colors. Source: `Quantity/`
//! Maps OCCT's Quantity_Color, Quantity_Length, Quantity_Angle, etc.

/// Named colors — the subset of `Quantity_NameOfColor`
/// (`Quantity_NameOfColor.hxx:25-530`, ~500 X11 names) that this port uses, so
/// the full enum is UNPORTED. Every variant below exists in OCCT under the same
/// spelling; `Quantity_NOC_SILVER`, `_NEUTRAL`, `_DARKBLUE` and `_LIGHTGREEN`
/// do **not** exist in OCCT 8.0, and the invented variants for them are gone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NameOfColor {
    Black, Blue, Brown, Cyan, Gold, Gray, Green, Magenta, Orange, Pink,
    Red, Violet, White, Yellow, DarkGreen, DarkOrange, LightBlue, LightGray,
    LightYellow,
}

/// RGB color (0.0—1.0). Source: `Quantity_Color.hxx`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Color { pub r: f32, pub g: f32, pub b: f32 }
impl Color {
    pub const BLACK:   Self = Self{r:0.,g:0.,b:0.};
    pub const WHITE:   Self = Self{r:1.,g:1.,b:1.};
    pub const RED:     Self = Self{r:1.,g:0.,b:0.};
    pub const GREEN:   Self = Self{r:0.,g:1.,b:0.};
    pub const BLUE:    Self = Self{r:0.,g:0.,b:1.};
    pub const YELLOW:  Self = Self{r:1.,g:1.,b:0.};
    pub const CYAN:    Self = Self{r:0.,g:1.,b:1.};
    pub const MAGENTA: Self = Self{r:1.,g:0.,b:1.};
    pub fn new(r:f32,g:f32,b:f32) -> Self { Self{r,g,b} }
    /// `Quantity_Color(theName, Quantity_TOC_RGB)` -> `valuesOf(name,
    /// Quantity_TOC_RGB)` (`Quantity_ColorTable.pxx`): the **RGB (linear)**
    /// column of the table entry, not the `sRGB` column. The previous values
    /// were the sRGB column truncated to two decimals, which is a different
    /// color (BROWN is off by ~1.7x); the line references below are the
    /// `RawColor(...)` rows of `Quantity_ColorTable.pxx`.
    pub fn from_name(n: NameOfColor) -> Self {
        match n {
            NameOfColor::Black      => Self{r:0.000000,g:0.000000,b:0.000000}, // `:14`
            NameOfColor::Blue       => Self{r:0.000000,g:0.000000,b:1.000000}, // `:171`
            NameOfColor::Brown      => Self{r:0.376262,g:0.023153,b:0.023153}, // `:185`
            NameOfColor::Cyan       => Self{r:0.000000,g:1.000000,b:1.000000}, // `:435`
            NameOfColor::Gold       => Self{r:1.000000,g:0.679542,b:0.000000}, // `:1009`
            NameOfColor::Gray       => Self{r:0.527115,g:0.527115,b:0.527115}, // `:1063`
            NameOfColor::Green      => Self{r:0.000000,g:1.000000,b:0.000000}, // `:1162`
            NameOfColor::Magenta    => Self{r:1.000000,g:0.000000,b:1.000000}, // `:1866`
            NameOfColor::Orange     => Self{r:1.000000,g:0.376262,b:0.000000}, // `:2287`
            NameOfColor::Pink       => Self{r:1.000000,g:0.527115,b:0.597201}, // `:2610`
            NameOfColor::Red        => Self{r:1.000000,g:0.000000,b:0.000000}, // `:2671`
            NameOfColor::Violet     => Self{r:0.854992,g:0.223228,b:0.854992}, // `:3349`
            NameOfColor::White      => Self{r:1.000000,g:1.000000,b:1.000000}, // `:3456`
            NameOfColor::Yellow     => Self{r:1.000000,g:1.000000,b:0.000000}, // `:3415`
            NameOfColor::DarkGreen  => Self{r:0.000000,g:0.127438,b:0.000000}, // `:489`
            NameOfColor::DarkOrange => Self{r:1.000000,g:0.262251,b:0.000000}, // `:559`
            NameOfColor::LightBlue  => Self{r:0.417885,g:0.686686,b:0.791299}, // `:1425`
            NameOfColor::LightGray  => Self{r:0.651406,g:0.651406,b:0.651406}, // `:1585`
            NameOfColor::LightYellow=> Self{r:1.000000,g:1.000000,b:0.745404}, // `:1815`
        }
    }
    pub fn to_rgb_u8(&self) -> (u8,u8,u8) { ((self.r*255.) as u8, (self.g*255.) as u8, (self.b*255.) as u8) }
}

/// RGBA color with alpha.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColorRGBA { pub r: f32, pub g: f32, pub b: f32, pub a: f32 }
impl Default for ColorRGBA { fn default() -> Self { Self{r:0.,g:0.,b:0.,a:1.} } }

/// Length quantity type. Source: `Quantity_Length`
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Length(pub f64);
impl Length {
    pub fn mm(v: f64)  -> Self { Self(v * 0.001) }
    pub fn cm(v: f64)  -> Self { Self(v * 0.01) }
    pub fn m(v: f64)   -> Self { Self(v) }
    pub fn km(v: f64)  -> Self { Self(v * 1000.0) }
    pub fn inch(v: f64)-> Self { Self(v * 0.0254) }
    pub fn as_meters(&self) -> f64 { self.0 }
    pub fn as_mm(&self) -> f64 { self.0 * 1000.0 }
}

/// Plane angle in radians. Source: `Quantity_PlaneAngle`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneAngle(pub f64);
impl PlaneAngle {
    pub fn rad(v: f64) -> Self { Self(v) }
    pub fn deg(v: f64) -> Self { Self(v * std::f64::consts::PI / 180.0) }
    pub fn as_radians(&self) -> f64 { self.0 }
    pub fn as_degrees(&self) -> f64 { self.0 * 180.0 / std::f64::consts::PI }
}

/// Volume quantity. Source: `Quantity_Volume`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Volume(pub f64);

/// Area quantity. Source: `Quantity_Area`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Area(pub f64);

/// Mass quantity. Source: `Quantity_Mass`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mass(pub f64);

/// Density quantity. Source: `Quantity_Density`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Density(pub f64);

/// Ratio (dimensionless). Source: `Quantity_Ratio`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ratio(pub f64);

/// Parameter on a curve (dimensionless). Source: `Quantity_Parameter`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parameter(pub f64);

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn length_convert() { let l=Length::mm(100.); assert!((l.as_meters()-0.1).abs()<1e-14); }
    #[test] fn angle_convert() { let a=PlaneAngle::deg(180.); assert!((a.as_radians()-std::f64::consts::PI).abs()<1e-14); }
    #[test] fn color_name() { assert_eq!(Color::from_name(NameOfColor::Red), Color::RED); }
    #[test] fn color_to_u8() { let c = Color::RED; assert_eq!(c.to_rgb_u8(), (255,0,0)); }
}
pub mod color_tools;
