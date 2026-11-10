//! `TopOpeBRepBuild_ShellFaceClassifier` — classifies faces against shells by
//! a 3D point taken on the first vertex of the face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ShellFaceClassifier.hxx:32-80`,
//! `.cxx:31-205` (OCCT 8.0.0).
//!
//! Shapes are shells, elements are faces (`hxx:32-33`). `CompareElement`
//! (`cxx:128-168`) accumulates the faces of one block into a `TopoDS_Shell`
//! and remembers it in `myFaceShellMap`, so that a face met again reuses the
//! shell of its block instead of starting a new one. `State`
//! (`cxx:179-205`) classifies the 3D point of `ResetElement` inside that shell
//! with `TopOpeBRepTool_SolidClassifier::Classify(SHE, P, Tol)`.

use std::collections::HashMap;

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;

use crate::abs::{ShapeType, State};
use crate::builder::TopoBuilder;
use crate::shape::{Face, Shell, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key, vertex_position, ShapeMapKey};
use crate::topopebrep::tool::{SolidClassifier, SurfaceAdaptor};

use super::block_builder::BlockBuilder;
use super::composite_classifier::CompositeClassifier;

/// `class TopOpeBRepBuild_ShellFaceClassifier` (`hxx:34`).
pub struct ShellFaceClassifier<'a> {
    /// `myBlockBuilder` (`hxx:67` of `CompositeClassifier`), the
    /// `const BlockBuilder&` of the constructor (`cxx:33-37`).
    my_block_builder: &'a BlockBuilder,
    /// `myFirstCompare` (`hxx:73`).
    my_first_compare: bool,
    /// `myPoint3d` (`hxx:74`).
    my_point3d: GpPnt,
    /// `myShell` (`hxx:75`).
    my_shell: Shell,
    /// `myBuilder` (`hxx:76`).
    my_builder: TopoBuilder,
    /// `mySolidClassifier` (`hxx:77`).
    my_solid_classifier: SolidClassifier,
    /// `myFaceShellMap` (`hxx:78`), a
    /// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`.
    my_face_shell_map: HashMap<ShapeMapKey, TopoShape>,
    // `myShape` (`hxx:79`) is declared but never read nor written by `.cxx`,
    // so it has no member here.
}

impl<'a> ShellFaceClassifier<'a> {
    /// `TopOpeBRepBuild_ShellFaceClassifier(BB)` (`cxx:33-37`).
    pub fn new(bb: &'a BlockBuilder) -> Self {
        Self {
            my_block_builder: bb,
            my_first_compare: false,
            my_point3d: GpPnt::new(0.0, 0.0, 0.0),
            my_shell: Shell::new(),
            my_builder: TopoBuilder::new(),
            my_solid_classifier: SolidClassifier::new(),
            my_face_shell_map: HashMap::new(),
        }
    }

    /// `Clear()` (`cxx:41-45`).
    pub fn clear(&mut self) {
        self.my_solid_classifier.clear();
        self.my_face_shell_map.clear();
    }
}

impl<'a> CompositeClassifier for ShellFaceClassifier<'a> {
    fn block_builder(&self) -> &BlockBuilder {
        self.my_block_builder
    }

    /// `CompareShapes(B1, B2)` (`cxx:49-63`).
    fn compare_shapes(&mut self, b1: &TopoShape, b2: &TopoShape) -> State {
        self.reset_shape(b1);
        if let Some(s) = Shell::wrap(b2.clone()) {
            self.my_shell = s;
        }
        self.my_solid_classifier.load_shell(&self.my_shell);
        self.state()
    }

    /// `CompareElementToShape(F, SHE)` (`cxx:67-81`).
    fn compare_element_to_shape(&mut self, f: &TopoShape, she: &TopoShape) -> State {
        self.reset_element(f);
        if let Some(s) = Shell::wrap(she.clone()) {
            self.my_shell = s;
        }
        self.my_solid_classifier.load_shell(&self.my_shell);
        self.state()
    }

    /// `ResetShape(SHE)` (`cxx:85-94`): `ResetElement` on the first face.
    ///
    /// `cxx:91-92` builds a `TopExp_Explorer` on `TopAbs_FACE` and dereferences
    /// `ex.Current()` without testing `ex.More()`; a shell without faces is not
    /// expected there.
    fn reset_shape(&mut self, she: &TopoShape) {
        let faces = explorer_shapes(she, ShapeType::Face);
        if let Some(first) = faces.first() {
            self.reset_element(first);
        }
    }

    /// `ResetElement(F)` (`cxx:98-124`).
    fn reset_element(&mut self, f: &TopoShape) {
        let t = f.shape_type();

        // initialize myPoint3d with first vertex of face <E>
        self.my_first_compare = true;

        let vertices = explorer_shapes(f, ShapeType::Vertex);
        if let Some(v) = vertices.first() {
            if let Some(vertex) = crate::shape::Vertex::wrap(v.clone()) {
                self.my_point3d = vertex_position(&vertex);
            }
        } else if t == ShapeType::Face {
            if let Some(face) = Face::wrap(f.clone()) {
                let bas = SurfaceAdaptor::of_face(&face);
                let u = 0.5 * (bas.first_u_parameter() + bas.last_u_parameter());
                let v = 0.5 * (bas.first_v_parameter() + bas.last_v_parameter());
                if let Some(surface) = bas.surface() {
                    self.my_point3d = surface.value(u, v);
                }
            }
        } else {
            self.my_point3d = GpPnt::new(0.0, 0.0, 0.0);
        }
    }

    /// `CompareElement(F)` (`cxx:128-168`).
    fn compare_element(&mut self, f: &TopoShape) -> bool {
        let mut b_ret = true;
        if self.my_first_compare {
            let found = self.my_face_shell_map.contains_key(&shape_map_key(f));
            if !found {
                // face F is the first one of a block of faces:
                // build a shell out of it and bind F to that shell
                let mut shell = Shell::new();
                self.my_builder.add(&mut shell.0, f);
                self.my_face_shell_map
                    .insert(shape_map_key(f), shell.0.clone());
                self.my_shell = shell;
            } else {
                // face F is the first one of an already shell-ized block of
                // faces: take back the shell it belongs to
                if let Some(s) = self.my_face_shell_map.get(&shape_map_key(f)) {
                    if let Some(shell) = Shell::wrap(s.clone()) {
                        self.my_shell = shell;
                    }
                }
                b_ret = !b_ret;
            }
            self.my_first_compare = false;
        } else {
            // F is not the first face of the block: myShell is defined and is
            // the shell of the block F belongs to
            self.my_builder.add(&mut self.my_shell.0, f);
        }
        b_ret
    }

    /// `State()` (`cxx:179-205`).
    fn state(&mut self) -> State {
        let tol3d = CONFUSION;
        self.my_solid_classifier
            .classify_shell(&self.my_shell, &self.my_point3d, tol3d)
    }
}
