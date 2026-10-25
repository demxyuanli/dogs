//! `TopOpeBRepDS_SurfaceData` — a `TopOpeBRepDS_Surface` plus its interference
//! list.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_SurfaceData.hxx:22-40`,
//! `TopOpeBRepDS_SurfaceData.cxx:22-31`.

use super::geometry_data::GeometryData;
use super::surface::Surface;

/// Source: `TopOpeBRepDS_SurfaceData.hxx:36-39`
#[derive(Debug, Clone, Default)]
pub struct SurfaceData {
    /// The `TopOpeBRepDS_GeometryData` base.
    pub(super) geometry: GeometryData,
    /// `mySurface`
    pub(super) surface: Surface,
}

impl SurfaceData {
    /// `TopOpeBRepDS_SurfaceData()` (`cxx:24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_SurfaceData(const TopOpeBRepDS_Surface& S)` (`cxx:28-31`).
    pub fn with_surface(s: Surface) -> Self {
        Self {
            geometry: GeometryData::new(),
            surface: s,
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

    /// The stored `TopOpeBRepDS_Surface`.
    pub fn surface(&self) -> &Surface {
        &self.surface
    }

    /// The stored `TopOpeBRepDS_Surface`, mutably.
    pub fn surface_mut(&mut self) -> &mut Surface {
        &mut self.surface
    }
}
