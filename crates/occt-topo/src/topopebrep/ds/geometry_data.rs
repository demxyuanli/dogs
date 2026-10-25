//! `TopOpeBRepDS_GeometryData` — mother class of `SurfaceData`, `CurveData`
//! and `PointData`: the list of interferences attached to one geometry.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_GeometryData.hxx:26-56`,
//! `TopOpeBRepDS_GeometryData.cxx:23-69`.

use super::interference::InterferenceHandle;

/// Source: `TopOpeBRepDS_GeometryData.hxx:29-56`
#[derive(Debug, Clone, Default)]
pub struct GeometryData {
    /// `myInterferences`
    interferences: Vec<InterferenceHandle>,
}

impl GeometryData {
    /// `TopOpeBRepDS_GeometryData()` (`cxx:24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Assign(Other)` (`cxx:36-47`): clears the list then appends Other's
    /// handles in order. The handles are shared, not copied, as in OCCT.
    pub fn assign(&mut self, other: &GeometryData) {
        self.interferences.clear();
        for it in &other.interferences {
            self.interferences.push(it.clone());
        }
    }

    /// `Interferences()` (`cxx:51-55`).
    pub fn interferences(&self) -> &[InterferenceHandle] {
        &self.interferences
    }

    /// `ChangeInterferences()` (`cxx:59-63`).
    pub fn change_interferences(&mut self) -> &mut Vec<InterferenceHandle> {
        &mut self.interferences
    }

    /// `AddInterference(I)` (`cxx:65-68`).
    pub fn add_interference(&mut self, i: InterferenceHandle) {
        self.interferences.push(i);
    }
}
