//! Unit conversion system. Source: `Units/` + `UnitsAPI/`
//! Maps OCCT's unit system to simple conversion factors.

/// Length unit with conversion to meters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LengthUnit { Millimeter, Centimeter, Meter, Kilometer, Inch, Foot, Mile }
impl LengthUnit {
    pub fn to_meters(&self) -> f64 { match self { Self::Millimeter=>0.001, Self::Centimeter=>0.01, Self::Meter=>1.0, Self::Kilometer=>1000.0, Self::Inch=>0.0254, Self::Foot=>0.3048, Self::Mile=>1609.344 } }
    pub fn convert(&self, value: f64, target: Self) -> f64 { value * self.to_meters() / target.to_meters() }
}

/// Angle unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AngleUnit { Radian, Degree, Gradian }
impl AngleUnit {
    pub fn to_radians(&self) -> f64 { match self { Self::Radian=>1.0, Self::Degree=>std::f64::consts::PI/180.0, Self::Gradian=>std::f64::consts::PI/200.0 } }
    pub fn convert(&self, value: f64, target: Self) -> f64 { value * self.to_radians() / target.to_radians() }
}

/// Mass unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MassUnit { Gram, Kilogram, Tonne, Pound, Ounce }
impl MassUnit {
    pub fn to_kg(&self) -> f64 { match self { Self::Gram=>0.001, Self::Kilogram=>1.0, Self::Tonne=>1000.0, Self::Pound=>0.45359237, Self::Ounce=>0.028349523125 } }
    pub fn convert(&self, value: f64, target: Self) -> f64 { value * self.to_kg() / target.to_kg() }
}

/// Time unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TimeUnit { Second, Minute, Hour }
impl TimeUnit {
    pub fn to_seconds(&self) -> f64 { match self { Self::Second=>1.0, Self::Minute=>60.0, Self::Hour=>3600.0 } }
    pub fn convert(&self, value: f64, target: Self) -> f64 { value * self.to_seconds() / target.to_seconds() }
}

/// Convert between OCCT's internal unit system (mm) and user units.
/// OCCT internally uses millimeters for all geometric computations.
pub struct UnitConverter { pub length_unit: LengthUnit }
impl UnitConverter {
    pub fn new(unit: LengthUnit) -> Self { Self { length_unit: unit } }
    pub fn to_occt(&self, value: f64) -> f64 { self.length_unit.convert(value, LengthUnit::Millimeter) }
    pub fn from_occt(&self, value: f64) -> f64 { LengthUnit::Millimeter.convert(value, self.length_unit) }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn length_mm_to_m() { assert!((LengthUnit::Millimeter.convert(1000.0, LengthUnit::Meter)-1.0).abs()<1e-14); }
    #[test] fn angle_deg_to_rad() { assert!((AngleUnit::Degree.convert(180.0, AngleUnit::Radian)-std::f64::consts::PI).abs()<1e-14); }
    #[test] fn mass_kg_to_g() { assert!((MassUnit::Kilogram.convert(1.0, MassUnit::Gram)-1000.0).abs()<1e-14); }
}
