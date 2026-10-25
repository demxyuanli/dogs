//! `TopOpeBRepTool_ShapeClassifier`.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_ShapeClassifier.hxx:38-149`,
//! `TopOpeBRepTool_ShapeClassifier.cxx:53-717` (OCCT 8.0.0).
//!
//! Ported:
//! - constructor `:53-58`, `SetReference` / `MapRef` (`:110-139`),
//!   `ClearCurrent` (`:78-92`), `ClearAll` (`:70-75`), `SameDomain`
//!   (`:96-107`), `State` (`:685-687`), `P2D` / `P3D` (`:692-709`),
//!   `HasAvLS` (`:714-717`), `ChangeSolidClassifier` (`:215-218`);
//! - `StateP2DReference` (`:585-637`), the `BRepTopAdaptor_FClass2d(F, 1e-8)`
//!   classification the 2d transitions of `TopOpeBRepDS_EXPORT` read through
//!   `FSC_StatePonFace`;
//! - `Perform` (`:318-482`), `FindEdge` (`:222-303`), `FindFace` (`:289-316`),
//!   `StateEdgeReference` (`:487-581`), `StateP3DReference` (`:642-681`) and the
//!   three `StateShapeShape` / `StateShapeReference` entry points
//!   (`:143-213`), on top of [`SolidClassifier`].
//!
//! UNPORTED:
//! - `StateShapeReference(S, AvLS)` (`:200-211`) — the `myPAvLS` overloads are
//!   ported, but nothing outside `Perform` sets `myPAvLS`.
//!
//! The C++ constructor initializes only `myP3Ddef` / `myP2Ddef`
//! (`ShapeClassifier.cxx:53-58`), so `mySameDomain` starts undefined and
//! `MapRef` reads it (`:122`). The port starts it at `-1`, the value
//! `ClearCurrent` writes (`:80`).

use occt_core::gp::{GpPnt, GpPnt2d};

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;

use super::closed::vertices;
use super::solid_classifier::SolidClassifier;

/// `TopOpeBRepTool_ShapeClassifier` (`hxx:38-149`).
#[derive(Debug)]
pub struct ShapeClassifier {
    /// `myS`
    s: Option<TopoShape>,
    /// `myRef`
    reference: Option<TopoShape>,
    /// `mymre` — `MapRef`'s edges of a same-domain reference face
    /// (`cxx:124-125`).
    map_ref_edges: Vec<TopoShape>,
    /// `mymren`
    map_ref_n: usize,
    /// `mymredone`
    map_ref_done: bool,
    /// `myState`
    state: State,
    /// `myP3Ddef`
    p3d_def: bool,
    /// `myP3D`
    p3d: GpPnt,
    /// `myP2Ddef`
    p2d_def: bool,
    /// `myP2D`
    p2d: GpPnt2d,
    /// `mySameDomain`
    same_domain: i32,
    /// `myAvS`
    av_s: Option<TopoShape>,
    /// `myPAvLS`
    av_ls: Option<Vec<TopoShape>>,
    /// `myMapAvS`
    map_av_s: Vec<TopoShape>,
    /// `myEdge`
    edge: Option<TopoShape>,
    /// `myFace`
    face: Option<TopoShape>,
    /// `mySolidClassifier`
    solid_classifier: SolidClassifier,
}

impl Default for ShapeClassifier {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeClassifier {
    /// `TopOpeBRepTool_ShapeClassifier()` (`cxx:53-58`): `myP3Ddef =
    /// myP2Ddef = false`; see the module note on `mySameDomain`.
    pub fn new() -> Self {
        Self {
            s: None,
            reference: None,
            map_ref_edges: Vec::new(),
            map_ref_n: 0,
            map_ref_done: false,
            state: State::Unknown,
            p3d_def: false,
            p3d: GpPnt::new(0.0, 0.0, 0.0),
            p2d_def: false,
            p2d: GpPnt2d::new(0.0, 0.0),
            same_domain: -1,
            av_s: None,
            av_ls: None,
            map_av_s: Vec::new(),
            edge: None,
            face: None,
            solid_classifier: SolidClassifier::new(),
        }
    }

    /// `TopOpeBRepTool_ShapeClassifier(SRef)` (`cxx:61-65`).
    pub fn with_reference(reference: &TopoShape) -> Self {
        let mut c = Self::new();
        c.reference = Some(reference.clone());
        c
    }

    /// `ClearAll()` (`cxx:70-75`).
    pub fn clear_all(&mut self) {
        self.clear_current();
        self.solid_classifier.clear();
    }

    /// `ClearCurrent()` (`cxx:78-92`). `myPAvLS` is *not* cleared — it is a
    /// loaned pointer in OCCT.
    pub fn clear_current(&mut self) {
        self.same_domain = -1;
        self.s = None;
        self.reference = None;
        self.av_s = None;
        self.map_av_s.clear();
        self.map_ref_edges.clear();
        self.map_ref_n = 0;
        self.map_ref_done = false;
        self.state = State::Unknown;
        self.edge = None;
        self.face = None;
        self.p3d_def = false;
        self.p2d_def = false;
    }

    /// `SameDomain()` (`cxx:96-98`).
    pub fn same_domain(&self) -> i32 {
        self.same_domain
    }

    /// `SameDomain(sam)` (`cxx:103-105`).
    pub fn set_same_domain(&mut self, sam: i32) {
        self.same_domain = sam;
    }

    /// `SetReference(SRef)` (`cxx:110-113`).
    pub fn set_reference(&mut self, reference: &TopoShape) {
        self.reference = Some(reference.clone());
        self.map_ref();
    }

    /// `MapRef()` (`cxx:118-139`): fills `mymre` when the reference is a
    /// same-domain face (`mySameDomain == 1`), and drops the single-edge case
    /// whose two vertices coincide.
    fn map_ref(&mut self) {
        self.map_ref_edges.clear();
        self.map_ref_n = 0;
        let Some(reference) = self.reference.clone() else {
            self.map_ref_done = true;
            return;
        };
        if reference.shape_type() == ShapeType::Face && self.same_domain == 1 {
            let mut ex = Explorer::new(&reference, ShapeType::Edge);
            while ex.more() {
                self.map_ref_edges.push(ex.current().clone());
                ex.next();
            }
            self.map_ref_n = self.map_ref_edges.len();
            if self.map_ref_n == 1 {
                let e = self.map_ref_edges[0].clone();
                let (v1, v2) = vertices(&e);
                let same = match (v1, v2) {
                    (Some(a), Some(b)) => crate::topo_tools::is_same(&a, &b),
                    _ => false,
                };
                if same {
                    self.map_ref_n = 0;
                }
            }
        }
        self.map_ref_done = true;
    }

    /// `StateP2DReference(P2D)` (`cxx:585-637`).
    ///
    /// The `mymren == 1` arm (`cxx:596-621`) answers `OUT` on an `EXTERNAL`
    /// boundary edge and `IN` on an `INTERNAL` one, stopping at the first other
    /// orientation with `myState` still `UNKNOWN`; the other arm classifies the
    /// point with `BRepTopAdaptor_FClass2d(F, 1e-8)` (`cxx:622-632`).
    pub fn state_p2d_reference(&mut self, p2d: &GpPnt2d) -> Result<(), String> {
        self.state = State::Unknown;
        let Some(reference) = self.reference.clone() else {
            return Ok(());
        };
        let t_r = reference.shape_type();
        if t_r != ShapeType::Face {
            return Err("StateShapeShape : bad operands".to_string());
        }
        if self.map_ref_n == 1 {
            let mut ex = Explorer::new(&reference, ShapeType::Edge);
            while ex.more() {
                let o = ex.current().orientation();
                if o == Orientation::External {
                    self.state = State::Out;
                } else if o == Orientation::Internal {
                    self.state = State::In;
                } else {
                    break;
                }
                ex.next();
            }
        } else {
            self.p2d = *p2d;
            self.p2d_def = true;
            let f = reference.oriented(Orientation::Forward);
            let tol_class = 1e-8;
            let classifier = FClass2d::new(&crate::shape::Face(f), tol_class)?;
            self.state = match classifier.perform_tab_orien(*p2d) {
                FaceState::In => State::In,
                FaceState::Out => State::Out,
                FaceState::On => State::On,
                FaceState::Unknown => State::Unknown,
            };
        }
        Ok(())
    }

    /// `State()` (`cxx:685-687`).
    pub fn state(&self) -> State {
        self.state
    }

    /// `P2D()` (`cxx:703-710`): `Err` when `myP2Ddef` is false, which raises
    /// `Standard_ProgramError` in OCCT.
    pub fn p2d(&self) -> Result<GpPnt2d, String> {
        if self.p2d_def {
            Ok(self.p2d)
        } else {
            Err("ShapeClassifier::P2D undefined".to_string())
        }
    }

    /// `P3D()` (`cxx:692-701`): `Err` when `myP3Ddef` is false.
    pub fn p3d(&self) -> Result<GpPnt, String> {
        if self.p3d_def {
            Ok(self.p3d)
        } else {
            Err("ShapeClassifier::P3D undefined".to_string())
        }
    }

    /// `HasAvLS()` (`cxx:714-717`).
    pub fn has_av_ls(&self) -> bool {
        self.av_ls.is_some()
    }

    /// `myS` getter.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.s.as_ref()
    }

    /// `ChangeSolidClassifier()` (`cxx:215-218`).
    pub fn change_solid_classifier(&mut self) -> &mut SolidClassifier {
        &mut self.solid_classifier
    }

    /// `StateShapeShape(S, SRef, samedomain)` (`cxx:143-155`).
    pub fn state_shape_shape(
        &mut self,
        s: &TopoShape,
        s_ref: &TopoShape,
        samedomain: i32,
    ) -> Result<State, String> {
        self.clear_current();
        self.same_domain = samedomain;
        self.s = Some(s.clone());
        self.av_s = None;
        self.av_ls = None;
        self.reference = Some(s_ref.clone());
        self.perform()
    }

    /// `StateShapeShape(S, AvS, SRef)` (`cxx:159-170`).
    pub fn state_shape_shape_av_s(
        &mut self,
        s: &TopoShape,
        av_s: &TopoShape,
        s_ref: &TopoShape,
    ) -> Result<State, String> {
        self.clear_current();
        self.s = Some(s.clone());
        self.av_s = Some(av_s.clone());
        self.av_ls = None;
        self.reference = Some(s_ref.clone());
        self.perform()
    }

    /// `StateShapeShape(S, AvS, SRef)` (`cxx:159-170`) with a null `AvS`:
    /// `Perform` takes the `myAvS`-null branch and the state is decided by the
    /// point chosen for `S` (`cxx:380-386`). `KPclasSS(S1, Snull, S2)`
    /// (`KPart.cxx:1605-1610`) is the caller.
    pub fn state_shape_shape_null_avoid(
        &mut self,
        s: &TopoShape,
        s_ref: &TopoShape,
    ) -> Result<State, String> {
        self.clear_current();
        self.s = Some(s.clone());
        self.av_s = None;
        self.av_ls = None;
        self.reference = Some(s_ref.clone());
        self.perform()
    }

    /// `StateShapeShape(S, AvLS, SRef)` (`cxx:173-186`).
    pub fn state_shape_shape_av_ls(
        &mut self,
        s: &TopoShape,
        av_ls: &[TopoShape],
        s_ref: &TopoShape,
    ) -> Result<State, String> {
        self.clear_current();
        self.s = Some(s.clone());
        self.av_s = None;
        self.av_ls = Some(av_ls.to_vec());
        self.reference = Some(s_ref.clone());
        self.perform()
    }

    /// `StateShapeReference(S, AvS)` (`cxx:190-197`): no `ClearCurrent`.
    pub fn state_shape_reference(
        &mut self,
        s: &TopoShape,
        av_s: &TopoShape,
    ) -> Result<State, String> {
        self.s = Some(s.clone());
        self.av_s = Some(av_s.clone());
        self.av_ls = None;
        self.perform()
    }

    /// `StateShapeReference(S, AvS)` (`cxx:190-198`) with a null `AvS`:
    /// `myAvS` stays null, so `Perform` skips the `myAvS`-branch and does not
    /// touch `myMapAvS` (`cxx:335-351`). `TopOpeBRepBuild_Tools::
    /// FindStateThroughVertex` (`TopOpeBRepBuild_Tools.cxx:378`) is the caller.
    pub fn state_shape_reference_null_avoid(&mut self, s: &TopoShape) -> Result<State, String> {
        self.s = Some(s.clone());
        self.av_s = None;
        self.av_ls = None;
        self.perform()
    }

    /// `StateShapeReference(S, AvLS)` (`cxx:202-210`): no `ClearCurrent`.
    pub fn state_shape_reference_av_ls(
        &mut self,
        s: &TopoShape,
        av_ls: &[TopoShape],
    ) -> Result<State, String> {
        self.s = Some(s.clone());
        self.av_s = None;
        self.av_ls = Some(av_ls.to_vec());
        self.perform()
    }

    /// `FindEdge()` (`cxx:222-237`).
    fn find_edge(&mut self) {
        self.edge = None;
        self.face = None;
        let Some(s) = self.s.clone() else {
            return;
        };
        let t = s.shape_type();
        let below_face = matches!(
            t,
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Solid | ShapeType::Shell
        );
        if below_face {
            self.find_face(&s);
            let face = self.face.clone();
            if let Some(face) = face {
                self.find_edge_of(&face);
            }
        } else {
            self.find_edge_of(&s);
        }
    }

    /// `FindEdge(S)` (`cxx:241-287`).
    fn find_edge_of(&mut self, s: &TopoShape) {
        self.edge = None;
        let is_av_ls = self.has_av_ls();
        let is_av_s = self.av_s.is_some();
        if s.is_null() {
            return;
        }
        let t_s = s.shape_type();
        let from = match self.face.clone() {
            Some(face) => face,
            None => s.clone(),
        };
        let mut ex = Explorer::new(&from, ShapeType::Edge);
        while ex.more() {
            let e = ex.current().clone();
            let mut to_avoid = false;
            if is_av_ls || is_av_s {
                to_avoid = to_avoid
                    || self
                        .map_av_s
                        .iter()
                        .any(|s| crate::topo_tools::is_same(s, &e));
                if let Some(av_s) = self.av_s.clone() {
                    to_avoid = to_avoid || crate::topo_tools::is_same(&av_s, &e);
                }
            } else if let Some(edge) = Edge::wrap(e.clone()) {
                if BRepTool::is_degenerated(&edge) {
                    to_avoid = t_s != ShapeType::Edge;
                }
            }
            if to_avoid {
                ex.next();
                continue;
            }
            self.edge = Some(e);
            break;
        }
    }

    /// `FindFace(S)` (`cxx:289-316`).
    fn find_face(&mut self, s: &TopoShape) {
        self.face = None;
        let is_av_ls = self.has_av_ls();
        let is_av_s = self.av_s.is_some();
        let mut ex = Explorer::new(s, ShapeType::Face);
        while ex.more() {
            let f = ex.current().clone();
            let mut to_avoid = false;
            if is_av_ls || is_av_s {
                to_avoid = to_avoid
                    || self
                        .map_av_s
                        .iter()
                        .any(|s| crate::topo_tools::is_same(s, &f));
                if let Some(av_s) = self.av_s.clone() {
                    to_avoid = to_avoid || crate::topo_tools::is_same(&av_s, &f);
                }
            }
            if to_avoid {
                ex.next();
                continue;
            }
            self.face = Some(f);
            break;
        }
    }

    /// `Perform()` (`cxx:318-482`).
    fn perform(&mut self) -> Result<State, String> {
        self.state = State::Unknown;
        let Some(s) = self.s.clone() else {
            return Ok(self.state);
        };
        let Some(reference) = self.reference.clone() else {
            return Ok(self.state);
        };

        if !self.map_ref_done {
            self.map_ref();
        }

        if let Some(av_s) = self.av_s.clone() {
            // `cxx:336-351`
            self.map_av_s.clear();
            let t_av_s = av_s.shape_type();
            if t_av_s == ShapeType::Face {
                self.map_av_s.push(av_s.clone());
                for e in edges_of(&av_s) {
                    self.map_av_s.push(e);
                }
            } else if t_av_s == ShapeType::Edge {
                for e in edges_of(&av_s) {
                    self.map_av_s.push(e);
                }
            }
        } else if let Some(av_ls) = self.av_ls.clone() {
            // `cxx:352-378`
            self.map_av_s.clear();
            let Some(first) = av_ls.first() else {
                return Err("StateShapeShape : bad operands".to_string());
            };
            let t_av_s = first.shape_type();
            if t_av_s == ShapeType::Face {
                for elem in &av_ls {
                    self.map_av_s.push(elem.clone());
                    for e in edges_of(elem) {
                        self.map_av_s.push(e);
                    }
                }
            } else if t_av_s == ShapeType::Edge {
                for elem in &av_ls {
                    for e in edges_of(elem) {
                        self.map_av_s.push(e);
                    }
                }
            }
        } else if s.shape_type() == ShapeType::Face {
            // `cxx:380-386`
            let mut p = GpPnt::new(0.0, 0.0, 0.0);
            let found = match Face::wrap(s.clone()) {
                Some(f) => crate::brep_class3d::find_a_point_in_the_face(&f, &mut p),
                None => false,
            };
            self.p3d_def = found;
            if found {
                self.p3d = p;
            }
        }

        let t_s = s.shape_type();
        let t_r = reference.shape_type();
        let ref_le_solid = matches!(
            t_r,
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Solid
        );

        match t_s {
            ShapeType::Vertex => {
                if ref_le_solid {
                    // `cxx:394-397`
                    let p3d = super::shape_tool::pnt_of(&s)?;
                    self.state_p3d_reference(&p3d)?;
                }
            }
            ShapeType::Edge | ShapeType::Wire => {
                // `cxx:399-413`
                if t_r == ShapeType::Face || ref_le_solid {
                    self.find_edge();
                    self.state_edge_reference()?;
                }
            }
            ShapeType::Face => {
                if t_r == ShapeType::Face {
                    // `cxx:417-439`
                    self.find_edge();
                    if self.same_domain == 1 {
                        self.state_edge_reference()?;
                    } else {
                        if !self.p3d_def {
                            let mut p = GpPnt::new(0.0, 0.0, 0.0);
                            let found = match Face::wrap(s.clone()) {
                                Some(f) => {
                                    crate::brep_class3d::find_a_point_in_the_face(&f, &mut p)
                                }
                                None => false,
                            };
                            self.p3d_def = found;
                            if found {
                                self.p3d = p;
                            }
                        }
                        if self.p3d_def {
                            let p = self.p3d;
                            self.state_p3d_reference(&p)?;
                        } else {
                            self.state = State::Unknown;
                            return Err("TopOpeBRepTool_ShapeClassifier !P3Ddef".to_string());
                        }
                    }
                } else if ref_le_solid {
                    // `cxx:441-452`
                    self.find_edge();
                    if self.p3d_def {
                        let p = self.p3d;
                        self.state_p3d_reference(&p)?;
                    } else {
                        self.state_edge_reference()?;
                    }
                }
            }
            ShapeType::Shell | ShapeType::Solid => {
                // `cxx:453-468`
                if ref_le_solid {
                    self.find_edge();
                    self.state_edge_reference()?;
                }
            }
            _ => {
                return Err("StateShapeShape : bad operands".to_string());
            }
        }

        // `cxx:473-483`: take the reference shape's orientation into account.
        let ori_ref = reference.orientation();
        if (ori_ref == Orientation::External || ori_ref == Orientation::Internal)
            && self.state == State::In
        {
            self.state = State::Out;
        }
        Ok(self.state)
    }

    /// `StateEdgeReference()` (`cxx:487-581`).
    fn state_edge_reference(&mut self) -> Result<(), String> {
        self.state = State::Unknown;
        let Some(edge) = self.edge.clone() else {
            return Ok(());
        };
        let Some(reference) = self.reference.clone() else {
            return Ok(());
        };
        let t_r = reference.shape_type();
        // `myEdge` is an edge of `myS`, not of `myRef` (`cxx:509`).
        if t_r == ShapeType::Face {
            let face = Face::wrap(reference).ok_or("StateShapeShape : bad face")?;
            if self.same_domain != 0 {
                // `cxx:514-535`
                let e = Edge::wrap(edge).ok_or("StateShapeShape : bad edge")?;
                let p_curve = crate::boptools_2d::make_2d(&e, &face)?;
                let (f2d, l2d) =
                    match crate::boptools_2d::curve_on_surface_range(&e, &face) {
                        Some((_, f, l)) => (f, l),
                        None => return Err("StateShapeShape : no 2d curve".to_string()),
                    };
                let t = 0.127956477;
                let p = (1.0 - t) * f2d + t * l2d;
                let p2d = p_curve.d0(p);
                self.state_p2d_reference(&p2d)?;
                return Ok(());
            }
            // `cxx:536-549`
            let e = Edge::wrap(edge).ok_or("StateShapeShape : bad edge")?;
            let p3d = point_on_curve(&e)?;
            self.state_p3d_reference(&p3d)?;
            return Ok(());
        }
        if matches!(
            t_r,
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Solid
        ) {
            // `cxx:551-577`
            let e = Edge::wrap(edge).ok_or("StateShapeShape : bad edge")?;
            if BRepTool::is_degenerated(&e) {
                let (v1, _) = vertices(&e);
                let v = v1.ok_or("StateShapeShape : no first vertex")?;
                let p3d = super::shape_tool::pnt_of(&v)?;
                self.state_p3d_reference(&p3d)?;
                return Ok(());
            }
            let p3d = point_on_curve(&e)?;
            self.state_p3d_reference(&p3d)?;
            return Ok(());
        }
        Err("StateShapeShape : bad operands".to_string())
    }

    /// `StateP3DReference(P3D)` (`cxx:642-681`).
    pub fn state_p3d_reference(&mut self, p3d: &GpPnt) -> Result<(), String> {
        self.state = State::Unknown;
        let Some(reference) = self.reference.clone() else {
            return Ok(());
        };
        let t_r = reference.shape_type();
        if t_r == ShapeType::Solid {
            // `cxx:650-658`
            self.p3d = *p3d;
            self.p3d_def = true;
            let solid = crate::shape::Solid::wrap(reference).ok_or("StateShapeShape : bad solid")?;
            let tol3d = occt_core::precision::CONFUSION;
            self.solid_classifier.classify_solid(&solid, p3d, tol3d);
            self.state = self.solid_classifier.state();
        } else if matches!(t_r, ShapeType::Compound | ShapeType::CompSolid) {
            // `cxx:659-678`
            self.p3d = *p3d;
            self.p3d_def = true;
            let mut ex = Explorer::new(&reference, ShapeType::Solid);
            while ex.more() {
                let so = ex.current().clone();
                if let Some(solid) = crate::shape::Solid::wrap(so) {
                    let tol3d = occt_core::precision::CONFUSION;
                    self.solid_classifier.classify_solid(&solid, p3d, tol3d);
                }
                self.state = self.solid_classifier.state();
                if self.state == State::In || self.state == State::On {
                    break;
                }
                ex.next();
            }
        } else {
            return Err("StateShapeShape : bad operands".to_string());
        }
        Ok(())
    }
}

/// `TopExp_Explorer(S, TopAbs_EDGE)` list.
fn edges_of(s: &TopoShape) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(s, ShapeType::Edge);
    while ex.more() {
        out.push(ex.current().clone());
        ex.next();
    }
    out
}

/// `BRep_Tool::Curve(myEdge, f3d, l3d)` then `C3D->Value(p)` with
/// `p = (1 - 0.127956477)·f3d + 0.127956477·l3d` (`cxx:538-547`, `:564-574`).
fn point_on_curve(e: &Edge) -> Result<GpPnt, String> {
    let adaptor = super::adaptor::CurveAdaptor::of_edge(e);
    if adaptor.curve().is_none() {
        return Err("StateShapeShape : no 3d curve".to_string());
    }
    let f3d = adaptor.first_parameter();
    let l3d = adaptor.last_parameter();
    let t = 0.127956477;
    let p = (1.0 - t) * f3d + t * l3d;
    adaptor
        .d0(p)
        .ok_or_else(|| "StateShapeShape : no 3d curve".to_string())
}
