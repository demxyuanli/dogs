//! `TopOpeBRepDS_CurveData` — a `TopOpeBRepDS_Curve` plus its interference
//! list.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_CurveData.hxx:22-40`,
//! `TopOpeBRepDS_CurveData.cxx:22-31`.

use super::curve::Curve;
use super::geometry_data::GeometryData;

/// Source: `TopOpeBRepDS_CurveData.hxx:36-39`
#[derive(Debug, Clone, Default)]
pub struct CurveData {
    /// The `TopOpeBRepDS_GeometryData` base.
    pub(super) geometry: GeometryData,
    /// `myCurve`
    pub(super) curve: Curve,
}

impl CurveData {
    /// `TopOpeBRepDS_CurveData()` (`cxx:24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_CurveData(const TopOpeBRepDS_Curve& C)` (`cxx:28-31`).
    pub fn with_curve(c: Curve) -> Self {
        Self {
            geometry: GeometryData::new(),
            curve: c,
        }
    }

    /// The inherited `TopOpeBRepDS_GeometryData` part.
    pub fn geometry(&self) -> &GeometryData {
        &self.geometry
    }

    /// The inherited `TopOpeBRepDS_GeometryData` part, mutably.
    pub fn geometry_mut(&mut self) -> &mut GeometryData {
        &mut self.geometry
    }

    /// The stored `TopOpeBRepDS_Curve`.
    pub fn curve(&self) -> &Curve {
        &self.curve
    }

    /// The stored `TopOpeBRepDS_Curve`, mutably.
    pub fn curve_mut(&mut self) -> &mut Curve {
        &mut self.curve
    }
}
