//! Physical quantities and named colors. Source: `Quantity/`
//! Maps OCCT's Quantity_Color, Quantity_Length, Quantity_Angle, etc.

/// Named colors from OCCT's Quantity_NameOfColor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NameOfColor {
    Black, Blue, Brown, Cyan, Gold, Gray, Green, Magenta, Orange, Pink,
    Red, Silver, Violet, White, Yellow, DarkBlue, DarkGreen, DarkOrange,
    LightBlue, LightGray, LightGreen, LightYellow, Neutral,
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
    pub fn from_name(n: NameOfColor) -> Self {
        match n {
            NameOfColor::Black=>Self::BLACK, NameOfColor::White=>Self::WHITE,
            NameOfColor::Red=>Self::RED, NameOfColor::Green=>Self::GREEN,
            NameOfColor::Blue=>Self::BLUE, NameOfColor::Yellow=>Self::YELLOW,
            NameOfColor::Cyan=>Self::CYAN, NameOfColor::Magenta=>Self::MAGENTA,
            NameOfColor::Gray=>Self{r:0.5,g:0.5,b:0.5},
            NameOfColor::Silver=>Self{r:0.75,g:0.75,b:0.75},
            NameOfColor::Brown=>Self{r:0.65,g:0.16,b:0.16},
            NameOfColor::Orange=>Self{r:1.,g:0.65,b:0.},
            NameOfColor::Gold=>Self{r:1.,g:0.84,b:0.},
            NameOfColor::Pink=>Self{r:1.,g:0.75,b:0.8},
            NameOfColor::Violet=>Self{r:0.93,g:0.51,b:0.93},
            NameOfColor::DarkBlue=>Self{r:0.,g:0.,b:0.55},
            NameOfColor::DarkGreen=>Self{r:0.,g:0.39,b:0.},
            NameOfColor::DarkOrange=>Self{r:1.,g:0.55,b:0.},
            NameOfColor::LightBlue=>Self{r:0.68,g:0.85,b:0.9},
            NameOfColor::LightGray=>Self{r:0.83,g:0.83,b:0.83},
            NameOfColor::LightGreen=>Self{r:0.56,g:0.93,b:0.56},
            NameOfColor::LightYellow=>Self{r:1.,g:1.,b:0.88},
            NameOfColor::Neutral=>Self{r:0.7,g:0.7,b:0.7},
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
