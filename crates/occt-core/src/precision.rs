//! Precision constants — exact OCCT values.
//! Source: `src/FoundationClasses/TKernel/Precision/Precision.hxx`

pub const ANGULAR: f64 = 1e-12;
pub const CONFUSION: f64 = 1e-7;
pub const SQUARE_CONFUSION: f64 = CONFUSION * CONFUSION;
pub const COMPUTATIONAL: f64 = f64::EPSILON;
pub const SQUARE_COMPUTATIONAL: f64 = COMPUTATIONAL * COMPUTATIONAL;
pub const INTERSECTION: f64 = CONFUSION * 0.01;
pub const APPROXIMATION: f64 = CONFUSION * 10.0;
pub const INFINITE: f64 = 2e100;
pub const PCONFUSION: f64 = CONFUSION * 0.01;
pub const PINTERSECTION: f64 = INTERSECTION * 0.01;
pub const PAPPROXIMATION: f64 = APPROXIMATION * 0.01;
pub const RESOLUTION: f64 = 1e-12;

pub struct Precision;

impl Precision {
    pub const ANGULAR: f64 = ANGULAR;
    pub const CONFUSION: f64 = CONFUSION;
    pub const SQUARE_CONFUSION: f64 = SQUARE_CONFUSION;
    pub const COMPUTATIONAL: f64 = COMPUTATIONAL;
    pub const SQUARE_COMPUTATIONAL: f64 = SQUARE_COMPUTATIONAL;
    pub const INTERSECTION: f64 = INTERSECTION;
    pub const APPROXIMATION: f64 = APPROXIMATION;
    pub const INFINITE: f64 = INFINITE;
    pub const PCONFUSION: f64 = PCONFUSION;
    pub const PINTERSECTION: f64 = PINTERSECTION;
    pub const PAPPROXIMATION: f64 = PAPPROXIMATION;
    #[inline] pub const fn parametric(p: f64, t: f64) -> f64 { p / t }
    #[inline] pub fn is_infinite(r: f64) -> bool { r.abs() >= (0.5 * INFINITE) }
    #[inline] pub const fn is_positive_infinite(r: f64) -> bool { r >= (0.5 * INFINITE) }
    #[inline] pub const fn is_negative_infinite(r: f64) -> bool { r <= -(0.5 * INFINITE) }
}
