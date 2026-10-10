//! Topology location — nested datum transforms. Source: `TopLoc/`
pub mod datum;
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

    /// `TopLoc_Location::IsEqual` (`TopLoc_Location.hxx:136`): the two locations
    /// hold the same series of elementary transforms. Identity links carry no
    /// datum in OCCT, so they are skipped. Datums compare by value here, as in
    /// `datum.rs`.
    pub fn is_equal(&self, other: &Self) -> bool {
        fn trsf_eq(a: &GpTrsf, b: &GpTrsf) -> bool {
            a.scale == b.scale && a.shape == b.shape && a.matrix == b.matrix && a.loc == b.loc
        }
        fn elementary(loc: &TopLocLocation) -> Vec<&GpTrsf> {
            let mut out = Vec::new();
            let mut cur = Some(loc);
            while let Some(node) = cur {
                if node.transform.form() != crate::gp::TrsfForm::Identity {
                    out.push(&node.transform);
                }
                cur = node.parent.as_deref();
            }
            out
        }
        let a = elementary(self);
        let b = elementary(other);
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| trsf_eq(x, y))
    }

    /// Transform a point through the full chain.
    pub fn transforms_point(&self, p: &GpPnt) -> GpPnt { p.transformed(&self.transformation()) }
}

impl Default for TopLocLocation { fn default() -> Self { Self::identity() } }
