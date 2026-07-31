//! Multi-level datum chain for nested transformations. Source: `TopLoc_Datum3D.hxx`
use crate::gp::{GpTrsf, GpPnt};

/// A 3D datum — wraps a transformation and supports equality by value.
#[derive(Debug, Clone)]
pub struct TopLocDatum3D { pub transform: GpTrsf }

impl TopLocDatum3D {
    pub fn new(t: &GpTrsf) -> Self { Self { transform: t.clone() } }
    pub fn identity() -> Self { Self { transform: GpTrsf::identity() } }
    pub fn transformation(&self) -> &GpTrsf { &self.transform }
    pub fn is_identity(&self) -> bool { self.transform.form() == crate::gp::TrsfForm::Identity }
}

/// Location = chain of datums, stored newest-first (like OCCT's TopLoc_Location).
#[derive(Debug, Clone, Default)]
pub struct TopLocLocationChain {
    datums: Vec<TopLocDatum3D>, // list of datums (bottom-up)
}

impl TopLocLocationChain {
    pub fn identity() -> Self { Self { datums: Vec::new() } }
    pub fn from_trsf(t: &GpTrsf) -> Self { Self { datums: vec![TopLocDatum3D::new(t)] } }
    pub fn is_identity(&self) -> bool { self.datums.is_empty() }

    /// Add a datum at the bottom (applied last).
    pub fn add_datum(&mut self, d: &TopLocDatum3D) {
        if !d.is_identity() { self.datums.push(d.clone()); }
    }

    /// Compose current location with t: new_loc = t * old_loc.
    pub fn composed(&self, t: &GpTrsf) -> Self {
        let mut result = Self::from_trsf(t);
        for d in &self.datums { result.datums.push(d.clone()); }
        result
    }

    /// Compute the cumulative transform (bottom-up product).
    pub fn transformation(&self) -> GpTrsf {
        let mut result = GpTrsf::identity();
        for d in &self.datums {
            result = result.multiplied(&d.transform);
        }
        result
    }

    /// Inverse location.
    pub fn inverted(&self) -> Self {
        let mut result = Self::identity();
        for d in self.datums.iter().rev() {
            let inv = d.transform.inverted().unwrap_or(GpTrsf::identity());
            result.datums.push(TopLocDatum3D::new(&inv));
        }
        result
    }

    /// Number of datums in chain.
    pub fn nb_datums(&self) -> usize { self.datums.len() }

    /// Transform a point through the full chain.
    pub fn transforms_point(&self, p: &GpPnt) -> GpPnt {
        let t = self.transformation();
        p.transformed(&t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{GpVec, GpPnt};

    #[test]
    fn chain_translation() {
        let mut v = GpVec::new(1., 0., 0.);
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&v);
        let loc = TopLocLocationChain::from_trsf(&t);
        let p = loc.transforms_point(&GpPnt::new(1., 1., 1.));
        assert!((p.x() - 2.).abs() < 1e-14);
    }

    #[test]
    fn invert_roundtrip() {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(3., 4., 5.));
        let loc = TopLocLocationChain::from_trsf(&t);
        let inv = loc.inverted();
        let p = inv.transforms_point(&GpPnt::new(0., 0., 0.));
        assert!((p.x() + 3.).abs() < 1e-12);
    }
}
