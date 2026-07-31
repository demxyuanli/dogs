//! Topology location — nested datum transforms. Source: `TopLoc/`
//! A TopLoc_Location is a chain of elementary gp_Trsf transforms.
use crate::gp::{GpTrsf, GpPnt};

/// Nested chain of transformations (like a singly-linked list of gp_Trsf).
/// Each node stores one transform and optionally points to a parent.
#[derive(Debug, Clone)]
pub struct TopLocLocation {
    transform: GpTrsf,
    parent: Option<Box<TopLocLocation>>,
}

impl TopLocLocation {
    /// Identity location (no transform).
    pub fn identity() -> Self { Self { transform: GpTrsf::identity(), parent: None } }

    /// Compose a transform with an existing location.
    pub fn composed(t: &GpTrsf, parent: TopLocLocation) -> Self {
        Self { transform: t.clone(), parent: Some(Box::new(parent)) }
    }

    /// Compute the cumulative transformation (product of chain).
    pub fn transformation(&self) -> GpTrsf {
        let mut result = GpTrsf::identity();
        let mut current = self;
        loop {
            result = (&result).multiplied(&current.transform);
            match &current.parent {
                Some(p) => current = p,
                None => break,
            }
        }
        result
    }

    /// Invert the location.
    pub fn inverted(&self) -> Self {
        let inv = self.transform.inverted().unwrap_or(GpTrsf::identity());
        match &self.parent {
            Some(p) => TopLocLocation::composed(&inv, p.inverted()),
            None => Self { transform: inv, parent: None },
        }
    }

    /// Is identity (no transform)?
    pub fn is_identity(&self) -> bool {
        self.transform.form() == crate::gp::TrsfForm::Identity && self.parent.is_none()
    }

    /// Transform a point through the full chain.
    pub fn transforms_point(&self, p: &GpPnt) -> GpPnt { p.transformed(&self.transformation()) }
}

impl Default for TopLocLocation { fn default() -> Self { Self::identity() } }
