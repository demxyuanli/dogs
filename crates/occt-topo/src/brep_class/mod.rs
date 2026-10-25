//! `BRepClass_FaceClassifier` — 2D classification of a point on a face.
//!
//! Two ports live here:
//!
//! * [`FClassifier`] is the faithful `BRepClass_FClassifier`
//!   (`BRepClass_FClassifier.cxx:25-100`) over
//!   `TopClass_FaceClassifier.pxx:28-176`, driving [`FaceExplorer`]
//!   (`BRepClass_FaceExplorer.cxx`), [`BRepClassFacePassiveClassifier`] and
//!   [`BRepClassIntersector`]. It exposes `Edge()` / `EdgeParameter()`, which
//!   `TopOpeBRep_VPointInterClassifier` reads.
//! * [`FaceClassifier`] is the older `BRepClass_FaceClassifier` entry point
//!   still used by the shape-analysis and shape-fix callers. It answers the
//!   same UV IN/ON/OUT question through [`crate::fclass2d::FClass2d`]
//!   (`IntTools_FClass2d` / CSLib_Class2d) and therefore carries no edge
//!   information. Switching it to [`FClassifier`] is a separate step.

use occt_core::gp::GpPnt2d;

use crate::fclass2d::{FaceState, FClass2d};
use crate::shape::Face;

mod edge;
mod f_classifier;
mod face_explorer;
mod face_passive_classifier;
mod intersector;
mod intersector_checks;
mod top_class_classifier2d;

pub use edge::BRepClassEdge;
pub use f_classifier::FClassifier;
pub use face_explorer::FaceExplorer;
pub use face_passive_classifier::BRepClassFacePassiveClassifier;
pub use intersector::BRepClassIntersector;
pub use top_class_classifier2d::TopClassClassifier2d;

/// `BRepClass_FaceClassifier` (see the module note: this is the
/// `IntTools_FClass2d`-based port, without the `FClassifier` edge walk).
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
