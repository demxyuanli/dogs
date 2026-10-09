//! Port of `BRepClass_FacePassiveClassifier` (`BRepClass_FacePassiveClassifier.cxx`,
//! `.hxx`). It composes the `TopClass_Classifier2d` state, the complex
//! curve transition and the `BRepClass_Intersector`, as the OCCT class does
//! through its members. Not wired into the face classifier yet.

use occt_core::gp::GpLin2d;

use crate::abs::Orientation;
use crate::fclass2d::FaceState;
use crate::intpatch::CurveTransition;

use super::edge::BRepClassEdge;
use super::intersector::BRepClassIntersector;
use super::top_class_classifier2d::TopClassClassifier2d;

/// `BRepClass_FacePassiveClassifier`.
pub struct BRepClassFacePassiveClassifier {
    /// The `TopClass_Classifier2d` fields (`myIsSet`, `myFirstCompare`,
    /// `myFirstTrans`, `myLin`, `myParam`, `myTolerance`, `myClosest`,
    /// `myState`, `myIsHeadOrEnd`).
    top: TopClassClassifier2d,
    /// `myTrans` (`TopTrans_CurveTransition`).
    trans: CurveTransition,
    /// `myIntersector` (`BRepClass_Intersector`).
    intersector: BRepClassIntersector,
}

impl Default for BRepClassFacePassiveClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepClassFacePassiveClassifier {
    /// `BRepClass_FacePassiveClassifier()`.
    pub fn new() -> Self {
        Self {
            top: TopClassClassifier2d::new(),
            trans: CurveTransition::new(),
            intersector: BRepClassIntersector::new(),
        }
    }

    /// `Reset(L, P, Tol)` (`BRepClass_FacePassiveClassifier.cxx:40-53`).
    pub fn reset(&mut self, l: &GpLin2d, p: f64, tol: f64) {
        self.top.reset(l, p, tol);
    }

    /// `Compare(E, Or)` (`BRepClass_FacePassiveClassifier.cxx:59-71`).
    pub fn compare(&mut self, e: &BRepClassEdge, orientation: Orientation) {
        self.top
            .compare(&mut self.intersector, &mut self.trans, e, orientation);
    }

    /// `Parameter()`.
    pub fn parameter(&self) -> f64 {
        self.top.param()
    }

    /// `Intersector()`.
    pub fn intersector(&mut self) -> &mut BRepClassIntersector {
        &mut self.intersector
    }

    /// `ClosestIntersection()`.
    pub fn closest_intersection(&self) -> i32 {
        self.top.closest()
    }

    /// `State()`.
    pub fn state(&self) -> FaceState {
        self.top.state()
    }

    /// `IsHeadOrEnd()`.
    pub fn is_head_or_end(&self) -> bool {
        self.top.is_head_or_end()
    }

    /// `myIsSet`.
    pub fn is_set(&self) -> bool {
        self.top.is_set()
    }
}
