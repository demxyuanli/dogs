//! Port of OCCT `BRepMesh_DegreeOfFreedom` — Wave 4 BRepMesh.
//!
//! The seven-state enum was already ported in Wave 1 together with
//! `BRepMesh_Deflection` (see [`super::deflection::DegreeOfFreedom`]). This
//! module is the dedicated home of the concept: it re-exports that enum, adds
//! the `from_u8` raw-value constructor and the bidirectional conversions to/from
//! the Delaunay [`VertexState`] (the same seven states, consumed by the Wave 2
//! `delaun_types` value types).
//!
//! Source: `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DegreeOfFreedom.hxx`.

use super::delaun_types::VertexState;

pub use super::deflection::DegreeOfFreedom;

impl DegreeOfFreedom {
    /// Raw enum value as in the OCCT C enum (`Free = 0`, ..., `Deleted = 6`).
    /// Returns `None` for out-of-range values.
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            0 => Some(Self::Free),
            1 => Some(Self::InVolume),
            2 => Some(Self::OnSurface),
            3 => Some(Self::OnCurve),
            4 => Some(Self::Fixed),
            5 => Some(Self::Frontier),
            6 => Some(Self::Deleted),
            _ => None,
        }
    }
}

/// Lossless conversion from the Delaunay vertex state.
impl From<VertexState> for DegreeOfFreedom {
    fn from(state: VertexState) -> Self {
        match state {
            VertexState::Free => Self::Free,
            VertexState::InVolume => Self::InVolume,
            VertexState::OnSurface => Self::OnSurface,
            VertexState::OnCurve => Self::OnCurve,
            VertexState::Fixed => Self::Fixed,
            VertexState::Frontier => Self::Frontier,
            VertexState::Deleted => Self::Deleted,
        }
    }
}

/// Lossless conversion back to the Delaunay vertex state.
impl From<DegreeOfFreedom> for VertexState {
    fn from(dof: DegreeOfFreedom) -> Self {
        match dof {
            DegreeOfFreedom::Free => Self::Free,
            DegreeOfFreedom::InVolume => Self::InVolume,
            DegreeOfFreedom::OnSurface => Self::OnSurface,
            DegreeOfFreedom::OnCurve => Self::OnCurve,
            DegreeOfFreedom::Fixed => Self::Fixed,
            DegreeOfFreedom::Frontier => Self::Frontier,
            DegreeOfFreedom::Deleted => Self::Deleted,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_u8_matches_occt_enum_order() {
        let all = [
            DegreeOfFreedom::Free,
            DegreeOfFreedom::InVolume,
            DegreeOfFreedom::OnSurface,
            DegreeOfFreedom::OnCurve,
            DegreeOfFreedom::Fixed,
            DegreeOfFreedom::Frontier,
            DegreeOfFreedom::Deleted,
        ];
        for (i, &dof) in all.iter().enumerate() {
            assert_eq!(DegreeOfFreedom::from_u8(i as u8), Some(dof));
            assert_eq!(dof.index(), i);
        }
        assert_eq!(DegreeOfFreedom::from_u8(7), None);
        assert_eq!(DegreeOfFreedom::from_u8(255), None);
    }

    #[test]
    fn degree_of_freedom_roundtrips_with_vertex_state() {
        for dof in [
            DegreeOfFreedom::Free,
            DegreeOfFreedom::InVolume,
            DegreeOfFreedom::OnSurface,
            DegreeOfFreedom::OnCurve,
            DegreeOfFreedom::Fixed,
            DegreeOfFreedom::Frontier,
            DegreeOfFreedom::Deleted,
        ] {
            let state: VertexState = dof.into();
            let back: DegreeOfFreedom = state.into();
            assert_eq!(back, dof, "roundtrip {dof:?}");
            assert_eq!(state.index(), dof.index(), "same numeric value");
        }
    }

    #[test]
    fn to_str_names_are_occt_names() {
        assert_eq!(DegreeOfFreedom::Free.to_str(), "Free");
        assert_eq!(DegreeOfFreedom::InVolume.to_str(), "InVolume");
        assert_eq!(DegreeOfFreedom::OnSurface.to_str(), "OnSurface");
        assert_eq!(DegreeOfFreedom::OnCurve.to_str(), "OnCurve");
        assert_eq!(DegreeOfFreedom::Fixed.to_str(), "Fixed");
        assert_eq!(DegreeOfFreedom::Frontier.to_str(), "Frontier");
        assert_eq!(DegreeOfFreedom::Deleted.to_str(), "Deleted");
    }
}
