//! `TopOpeBRepTool_SolidClassifier`.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_SolidClassifier.cxx:21-228`
//! (OCCT 8.0.0).
//!
//! A cache of `BRepClass3d_SolidClassifier` objects keyed by the classified
//! solid or shell — `myShapeClassifierMap` (`cxx:47-60`, `cxx:63-79`). The two
//! `Classify` overloads (`cxx:99-160`, `cxx:162-179`) forward to
//! `BRepClass3d_SolidClassifier::Perform` and adjust the state for an
//! `EXTERNAL` / `INTERNAL` deciding face (`cxx:119-158`).
//!
//! `LoadShell` builds `mySolid` out of the shell (`cxx:69-71`) and caches the
//! classifier under the *shell* key (`cxx:76`), so the solid is rebuilt only on
//! the first `Classify` of that shell.

use std::collections::HashMap;

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, State};
use crate::fclass2d::FaceState;
use crate::brep_class3d::SolidClassifier as BRepClass3dSolidClassifier;
use crate::builder::TopoBuilder;
use crate::shape::{Shell, Solid};
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `TopOpeBRepTool_SolidClassifier` (`hxx:26-87`).
pub struct SolidClassifier {
    /// `myPClassifier` — the entry of `myShapeClassifierMap` in use, by key.
    current: Option<ShapeMapKey>,
    /// `myShapeClassifierMap`: `NCollection_DataMap<TopoDS_Shape, void*>`,
    /// storing the classifiers themselves instead of raw pointers.
    shape_classifier_map: HashMap<ShapeMapKey, BRepClass3dSolidClassifier>,
    /// `myState`
    state: State,
    /// `myShell`
    shell: Option<Shell>,
    /// `mySolid`
    solid: Option<Solid>,
}

impl Default for SolidClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl SolidClassifier {
    /// `TopOpeBRepTool_SolidClassifier()` (`cxx:23-26`): `Clear()`.
    pub fn new() -> Self {
        Self {
            current: None,
            shape_classifier_map: HashMap::new(),
            state: State::Unknown,
            shell: None,
            solid: None,
        }
    }

    /// `LoadSolid(SOL)` (`cxx:47-60`).
    pub fn load_solid(&mut self, sol: &Solid) {
        let key = shape_map_key(&sol.0);
        if !self.shape_classifier_map.contains_key(&key) {
            let mut classifier = BRepClass3dSolidClassifier::new();
            classifier.load(sol.0.clone());
            self.shape_classifier_map.insert(key.clone(), classifier);
        }
        self.current = Some(key);
    }

    /// `LoadShell(SHE)` (`cxx:63-79`).
    pub fn load_shell(&mut self, she: &Shell) {
        let key = shape_map_key(&she.0);
        if !self.shape_classifier_map.contains_key(&key) {
            let builder = TopoBuilder::new();
            let solid = builder.make_solid(std::slice::from_ref(she));
            self.solid = Some(solid.clone());
            let mut classifier = BRepClass3dSolidClassifier::new();
            classifier.load(solid.0.clone());
            self.shape_classifier_map.insert(key.clone(), classifier);
        }
        self.current = Some(key);
    }

    /// `Clear()` (`cxx:87-95`): `myPClassifier`, the map, `myState`, `myShell`
    /// and `mySolid`.
    pub fn clear(&mut self) {
        self.current = None;
        self.shape_classifier_map.clear();
        self.state = State::Unknown;
        self.shell = None;
        self.solid = None;
    }

    /// `Classify(SOL, P, Tol)` (`cxx:99-160`).
    pub fn classify_solid(&mut self, sol: &Solid, p: &GpPnt, tol: f64) -> State {
        self.current = None;
        self.state = State::Unknown;
        self.load_solid(sol);
        let Some(key) = self.current.clone() else {
            return self.state;
        };
        let Some(classifier) = self.shape_classifier_map.get_mut(&key) else {
            return self.state;
        };
        classifier.perform(p, tol);
        self.state = to_state(classifier.state());
        // `const TopoDS_Shape& fres = myPClassifier->Face();` (`cxx:113`).
        let Some(fres) = classifier.face() else {
            return self.state;
        };
        let ofres = fres.0.orientation();
        if ofres == Orientation::External {
            self.state = match self.state {
                State::In => State::Out,
                State::Out => State::Out,
                State::On => State::On,
                State::Unknown => State::Out,
            };
        } else if ofres == Orientation::Internal {
            self.state = match self.state {
                State::In => State::In,
                State::Out => State::In,
                State::On => State::On,
                State::Unknown => State::In,
            };
        }
        self.state
    }

    /// `Classify(SHE, P, Tol)` (`cxx:162-179`): no `EXTERNAL` / `INTERNAL`
    /// face adjustment in this overload.
    pub fn classify_shell(&mut self, she: &Shell, p: &GpPnt, tol: f64) -> State {
        self.current = None;
        self.state = State::Unknown;
        self.load_shell(she);
        let Some(key) = self.current.clone() else {
            return self.state;
        };
        let Some(classifier) = self.shape_classifier_map.get_mut(&key) else {
            return self.state;
        };
        classifier.perform(p, tol);
        self.state = to_state(classifier.state());
        self.state
    }

    /// `State()` (`cxx:181-184`).
    pub fn state(&self) -> State {
        self.state
    }
}

impl std::fmt::Debug for SolidClassifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SolidClassifier")
            .field("state", &self.state)
            .field("cached", &self.shape_classifier_map.len())
            .finish()
    }
}

/// `BRepClass3d_SClassifier` states as `TopAbs_State`.
fn to_state(s: FaceState) -> State {
    match s {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}
