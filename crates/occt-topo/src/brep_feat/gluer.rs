//! `BRepFeat_Gluer` — glueing a new shape onto a basis shape.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Gluer.hxx:1-89`,
//! `BRepFeat_Gluer.cxx:24-62`, `BRepFeat_Gluer.lxx:19-66`.
//!
//! The class is the `BRepBuilderAPI_MakeShape` facade of `LocOpe_Gluer`; the
//! two `Bind` overloads invert their argument pair on the way in
//! (`lxx:25-27`, `:31-34`): `BRepFeat_Gluer(Snew, Sbase)` builds
//! `LocOpe_Gluer(Sbase, Snew)`.

use crate::locope::{Gluer, LocOpeOperation};
use crate::shape::{Edge, Face, TopoShape};

/// `BRepFeat_Gluer` (`hxx:49-103`).
#[derive(Debug, Clone)]
pub struct FeatureGluer {
    my_gluer: Gluer,
    my_shape: Option<TopoShape>,
    my_done: bool,
}

impl Default for FeatureGluer {
    fn default() -> Self {
        Self::new()
    }
}

impl FeatureGluer {
    /// `BRepFeat_Gluer()` (`lxx:19`).
    pub fn new() -> Self {
        Self {
            my_gluer: Gluer::new(),
            my_shape: None,
            my_done: false,
        }
    }

    /// `BRepFeat_Gluer(Snew, Sbase)` (`lxx:23-27`) — argument order inverted
    /// into `LocOpe_Gluer(Sbase, Snew)`.
    pub fn with_shapes(snew: &TopoShape, sbase: &TopoShape) -> Self {
        Self {
            my_gluer: Gluer::with_shapes(sbase, snew),
            my_shape: None,
            my_done: false,
        }
    }

    /// `Init(Snew, Sbase)` (`lxx:31-34`), same inversion.
    pub fn init(&mut self, snew: &TopoShape, sbase: &TopoShape) {
        self.my_gluer.init(sbase, snew);
        self.my_shape = None;
        self.my_done = false;
    }

    /// `Bind(Fnew, Fbase)` (`lxx:38-41`).
    pub fn bind_face(&mut self, fnew: &Face, fbase: &Face) -> Result<(), String> {
        self.my_gluer.bind_face(fnew, fbase)
    }

    /// `Bind(Enew, Ebase)` (`lxx:45-48`).
    pub fn bind_edge(&mut self, enew: &Edge, ebase: &Edge) -> Result<(), String> {
        self.my_gluer.bind_edge(enew, ebase)
    }

    /// `OpeType()` (`lxx:52-55`).
    pub fn ope_type(&self) -> LocOpeOperation {
        self.my_gluer.ope_type()
    }

    /// `BasisShape()` (`lxx:59-62`).
    pub fn basis_shape(&self) -> Option<&TopoShape> {
        self.my_gluer.basis_shape()
    }

    /// `GluedShape()` (`lxx:66-69`).
    pub fn glued_shape(&self) -> Option<&TopoShape> {
        self.my_gluer.glued_shape()
    }

    /// `Build()` (`cxx:26-35`): `Perform`, then `Done()` + `myShape` on
    /// success and `NotDone()` otherwise.
    pub fn build(&mut self) -> Result<(), String> {
        self.my_gluer.perform()?;
        if self.my_gluer.is_done() {
            self.my_done = true;
            self.my_shape = self.my_gluer.resulting_shape().ok().cloned();
            Ok(())
        } else {
            self.my_done = false;
            self.my_shape = None;
            Ok(())
        }
    }

    /// `IsDone()` of the inherited `BRepBuilderAPI_MakeShape`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// The inherited `Shape()`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `IsDeleted(F)` (`cxx:42-45`): `DescendantFaces(F).IsEmpty()`.
    ///
    /// `DescendantFaces` raises `StdFail_NotDone` before `Perform`
    /// (`LocOpe_Gluer.cxx:339-342`); the port answers the empty list, so a
    /// call before `Build` reads as deleted, exactly as an empty result would.
    pub fn is_deleted(&self, f: &Face) -> bool {
        self.descendants(f).is_empty()
    }

    /// `Modified(F)` (`cxx:47-62`): the descendants when the first one is not
    /// `F` itself, otherwise the empty list.
    pub fn modified(&self, f: &TopoShape) -> Vec<TopoShape> {
        if f.shape_type() != crate::abs::ShapeType::Face {
            return Vec::new();
        }
        let ls = self.descendants(&Face(f.clone()));
        if !ls.is_empty() {
            if !crate::topo_tools::is_same(&ls[0], f) {
                return ls;
            }
        }
        Vec::new()
    }

    fn descendants(&self, f: &Face) -> Vec<TopoShape> {
        self.my_gluer
            .descendant_faces(f)
            .map(|l| l.to_vec())
            .unwrap_or_default()
    }
}
