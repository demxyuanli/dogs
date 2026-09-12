//! `BRepClass_FaceClassifier` — 2D classification of a point on a face.
//!
//! Source: `BRepClass_FaceClassifier.cxx`. The OCCT classifier walks the face
//! edges with `BRepClass_FClassifier`; this port uses the already-translated
//! [`crate::fclass2d::FClass2d`] (`IntTools_FClass2d` / CSLib_Class2d) as the
//! same UV IN/ON/OUT test.

use occt_core::gp::GpPnt2d;

use crate::fclass2d::{FaceState, FClass2d};
use crate::shape::Face;

/// `BRepClass_FaceClassifier`.
pub struct FaceClassifier {
    inner: Option<FClass2d>,
    state: FaceState,
}

impl FaceClassifier {
    pub fn new() -> Self {
        Self {
            inner: None,
            state: FaceState::Unknown,
        }
    }

    /// `Perform(F, Puv, Tol)`.
    pub fn perform(&mut self, face: &Face, puv: GpPnt2d, tol: f64) {
        self.state = FaceState::Unknown;
        match FClass2d::new(face, tol) {
            Ok(cl) => {
                self.state = cl.perform(puv);
                self.inner = Some(cl);
            }
            Err(_) => self.inner = None,
        }
    }

    pub fn state(&self) -> FaceState {
        self.state
    }

    pub fn perform_infinite_point(&self) -> FaceState {
        self.inner
            .as_ref()
            .map(|c| c.perform_infinite_point())
            .unwrap_or(FaceState::Unknown)
    }
}

impl Default for FaceClassifier {
    fn default() -> Self {
        Self::new()
    }
}
