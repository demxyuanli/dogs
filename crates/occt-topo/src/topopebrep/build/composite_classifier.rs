//! `TopOpeBRepBuild_CompositeClassifier` — classifies composite loops, i.e.
//! loops that are either a shape or a block of elements.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_CompositeClassifier.hxx:28-64`,
//! `.cxx:36-132`.
//!
//! `TopOpeBRepBuild_CompositeClassifier` is an abstract base: `Compare`
//! (`cxx:46-131`) is implemented in terms of the pure virtuals `CompareShapes`,
//! `CompareElementToShape`, `ResetShape`, `ResetElement`, `CompareElement` and
//! `State` (`hxx:38-58`). The port keeps that as a trait whose `compare_loops`
//! has the body of `Compare`; [`LoopClassifier`] is then provided for every
//! implementor.

use crate::abs::State;
use crate::shape::TopoShape;

use super::block_builder::BlockBuilder;
use super::loop_classifier::LoopClassifier;
use super::loop_shape::LoopHandle;

/// `class TopOpeBRepBuild_CompositeClassifier` (`hxx:33`).
///
/// `myBlockBuilder` (`hxx:62`) is the `const BlockBuilder&` the constructor
/// stores (`cxx:36-40`); the port keeps an owned handle-free reference through
/// [`CompositeClassifier::block_builder`].
pub trait CompositeClassifier {
    /// The `TopOpeBRepBuild_BlockBuilder` the constructor was given.
    fn block_builder(&self) -> &BlockBuilder;

    /// `CompareShapes(B1, B2)` (`hxx:38-40`).
    fn compare_shapes(&mut self, b1: &TopoShape, b2: &TopoShape) -> State;

    /// `CompareElementToShape(E, B)` (`hxx:42-44`).
    fn compare_element_to_shape(&mut self, e: &TopoShape, b: &TopoShape) -> State;

    /// `ResetShape(B)` (`hxx:46-48`).
    fn reset_shape(&mut self, b: &TopoShape);

    /// `ResetElement(E)` (`hxx:50-52`).
    fn reset_element(&mut self, e: &TopoShape);

    /// `CompareElement(E)` (`hxx:54-56`).
    fn compare_element(&mut self, e: &TopoShape) -> bool;

    /// `State()` (`hxx:58-60`).
    fn state(&mut self) -> State;

    /// The virtual `Compare(L1, L2)` (`hxx:37-38`).
    ///
    /// The base body is `cxx:46-131` and is the default; the
    /// `WireEdgeClassifier` override (`WireEdgeClassifier.cxx:62-163`) replaces
    /// it.
    fn compare(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State {
        self.compare_loops(l1, l2)
    }

    /// `Compare(L1, L2)` (`cxx:46-131`), the base body.
    fn compare_loops(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State {
        let (isshape1, isshape2, shape1, shape2, bit1, bit2) = {
            let a = l1.read().expect("poisoned Loop lock");
            let b = l2.read().expect("poisoned Loop lock");
            (
                a.is_shape(),
                b.is_shape(),
                a.shape().cloned(),
                b.shape().cloned(),
                *a.block_iterator(),
                *b.block_iterator(),
            )
        };

        let mut state = State::Unknown;

        if isshape2 && isshape1 {
            // L1 is Shape, L2 is Shape
            if let (Some(s1), Some(s2)) = (shape1.as_ref(), shape2.as_ref()) {
                state = self.compare_shapes(s1, s2);
            }
        } else if isshape2 && !isshape1 {
            // L1 is Block, L2 is Shape
            let mut bit1 = bit1;
            bit1.initialize();
            let mut yena1 = bit1.more();
            while yena1 {
                let s1 = match self.block_builder().element_of_iterator(&bit1) {
                    Some(s) => s.clone(),
                    None => break,
                };
                let s2 = match shape2.as_ref() {
                    Some(s) => s.clone(),
                    None => break,
                };
                state = self.compare_element_to_shape(&s1, &s2);
                yena1 = false;
                if state == State::Unknown {
                    if bit1.more() {
                        bit1.next();
                    }
                    yena1 = bit1.more();
                }
            }
        } else if !isshape2 && isshape1 {
            // L1 is Shape, L2 is Block
            let s1 = match shape1.as_ref() {
                Some(s) => s.clone(),
                None => return state,
            };
            self.reset_shape(&s1);
            let mut bit2 = bit2;
            bit2.initialize();
            while bit2.more() {
                let s2 = match self.block_builder().element_of_iterator(&bit2) {
                    Some(s) => s.clone(),
                    None => break,
                };
                if !self.compare_element(&s2) {
                    break;
                }
                bit2.next();
            }
            state = self.state();
        } else if !isshape2 && !isshape1 {
            // L1 is Block, L2 is Block
            let mut bit1 = bit1;
            bit1.initialize();
            let mut yena1 = bit1.more();
            while yena1 {
                let s1 = match self.block_builder().element_of_iterator(&bit1) {
                    Some(s) => s.clone(),
                    None => break,
                };
                self.reset_element(&s1);
                let mut bit2 = bit2;
                bit2.initialize();
                while bit2.more() {
                    let s2 = match self.block_builder().element_of_iterator(&bit2) {
                        Some(s) => s.clone(),
                        None => break,
                    };
                    self.compare_element(&s2);
                    bit2.next();
                }
                state = self.state();
                yena1 = false;
                if state == State::Unknown {
                    if bit1.more() {
                        bit1.next();
                    }
                    yena1 = bit1.more();
                }
            }
        }

        state
    }
}

impl<T: CompositeClassifier + ?Sized> LoopClassifier for T {
    /// `TopOpeBRepBuild_CompositeClassifier::Compare` (`cxx:46-131`), i.e. the
    /// virtual `Compare` of the concrete subclass.
    fn compare(&mut self, l1: &LoopHandle, l2: &LoopHandle) -> State {
        CompositeClassifier::compare(self, l1, l2)
    }
}
