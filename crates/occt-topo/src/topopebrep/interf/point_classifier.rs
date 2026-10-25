//! `TopOpeBRep_PointClassifier` — classifies a 2D point of a face into
//! `TopAbs_State`, caching one topological tool per face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_PointClassifier.hxx`,
//! `.cxx:26-73`.
//!
//! Translation notes:
//! - `myHSurface` is a [`SurfaceAdaptor`] (`BRepAdaptor_Surface`);
//!   `myHSurface->Initialize(F)` is [`SurfaceAdaptor::initialize`].
//! - `BRepTopAdaptor_TopolTool` is [`crate::geom_int::TopolTool`];
//!   `Classify(P2d, Tol)` is `TopolTool::classify` and returns
//!   [`FaceState`], which maps one-to-one onto `TopAbs_State`.
//! - `myTopolToolMap` is a `TopTools_ShapeMapHasher` map, i.e. keyed by
//!   `IsSame`; that is [`crate::topo_tools::shape_map_key`].

use std::collections::HashMap;

use occt_core::gp::GpPnt2d;

use crate::abs::State;
use crate::fclass2d::FaceState;
use crate::geom_int::TopolTool;
use crate::shape::Face;
use crate::topo_tools::{shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `TopAbs_State` of a [`FaceState`] — the two enums are the same four cases
/// (`TopOpeBRep_PointClassifier.cxx:63`).
fn state_of(face_state: FaceState) -> State {
    match face_state {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}

/// `class TopOpeBRep_PointClassifier` (`hxx`).
#[derive(Clone)]
pub struct PointClassifier {
    /// `myTopolToolMap` — one tool per face.
    topol_tool_map: HashMap<ShapeMapKey, TopolTool>,
    /// `myHSurface` (`BRepAdaptor_Surface`).
    h_surface: SurfaceAdaptor,
    /// `myTopolTool`.
    topol_tool: Option<TopolTool>,
    /// `myState`.
    state: State,
}

impl Default for PointClassifier {
    /// `TopOpeBRep_PointClassifier()` (`cxx:26-30`): the adaptor is created
    /// empty and `Init()` runs.
    fn default() -> Self {
        Self::new()
    }
}

impl PointClassifier {
    /// `TopOpeBRep_PointClassifier()` (`cxx:26-30`).
    pub fn new() -> Self {
        let mut c = Self {
            topol_tool_map: HashMap::new(),
            h_surface: empty_surface_adaptor(),
            topol_tool: None,
            state: State::Unknown,
        };
        c.init();
        c
    }

    /// `Init()` (`cxx:32-36`).
    pub fn init(&mut self) {
        self.topol_tool_map.clear();
        self.state = State::Unknown;
    }

    /// `Load(F)` (`cxx:40-55`).
    pub fn load(&mut self, f: &Face) {
        let key = shape_map_key(&f.0);
        if let Some(tool) = self.topol_tool_map.get(&key) {
            self.topol_tool = Some(tool.clone());
            return;
        }
        self.h_surface.initialize(f);
        // OCCT builds `BRepTopAdaptor_TopolTool(myHSurface)` unconditionally;
        // with no surface under the face the port has nothing to classify in,
        // so no tool is cached and `Classify` answers `TopAbs_UNKNOWN`.
        let Some(s) = self.h_surface.surface() else {
            self.topol_tool = None;
            return;
        };
        let tool = TopolTool::from_face(f, s.as_ref());
        self.topol_tool_map.insert(key, tool.clone());
        self.topol_tool = Some(tool);
    }

    /// `Classify(F, P2d, Tol)` (`cxx:59-66`).
    pub fn classify(&mut self, f: &Face, p2d: GpPnt2d, tol: f64) -> State {
        self.state = State::Unknown;
        self.load(f);
        self.state = match self.topol_tool.as_ref() {
            Some(tool) => state_of(tool.classify(p2d, tol)),
            None => State::Unknown,
        };
        self.state
    }

    /// `State()` (`cxx:70-73`).
    pub fn state(&self) -> State {
        self.state
    }
}

/// `new BRepAdaptor_Surface()` (`cxx:28`): an adaptor with no face and no
/// surface, which `Load(F)` fills through [`SurfaceAdaptor::initialize`].
fn empty_surface_adaptor() -> SurfaceAdaptor {
    SurfaceAdaptor::of_face(&Face::new())
}
