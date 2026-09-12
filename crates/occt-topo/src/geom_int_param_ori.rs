//! `GeomInt_ParameterAndOrientation`. Source: the .cxx / .hxx of the same name.

use crate::abs::Orientation;

/// Parameter on a restriction line plus orientations on each surface domain.
#[derive(Debug, Clone, Copy)]
pub struct ParameterAndOrientation {
    prm: f64,
    or1: Orientation,
    or2: Orientation,
}

impl ParameterAndOrientation {
    pub fn new() -> Self {
        Self {
            prm: 0.0,
            or1: Orientation::Forward,
            or2: Orientation::Forward,
        }
    }

    pub fn with(p: f64, or1: Orientation, or2: Orientation) -> Self {
        Self { prm: p, or1, or2 }
    }

    pub fn set_orientation1(&mut self, or1: Orientation) {
        self.or1 = or1;
    }

    pub fn set_orientation2(&mut self, or2: Orientation) {
        self.or2 = or2;
    }

    pub fn parameter(&self) -> f64 {
        self.prm
    }

    pub fn orientation1(&self) -> Orientation {
        self.or1
    }

    pub fn orientation2(&self) -> Orientation {
        self.or2
    }
}

impl Default for ParameterAndOrientation {
    fn default() -> Self {
        Self::new()
    }
}
