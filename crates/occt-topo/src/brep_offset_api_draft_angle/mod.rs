//! `BRepOffsetAPI_DraftAngle` - the `BRepBuilderAPI_MakeShape` facade of
//! `Draft_Modification` (`ModelingAlgorithms/TKOffset/BRepOffsetAPI/`,
//! `BRepOffsetAPI_DraftAngle.cxx:52-268` + `.hxx:38-221`, OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepOffsetAPI_DraftAngle()` | `cxx:52` | [`BRepOffsetApiDraftAngle::new`] |
//! | `BRepOffsetAPI_DraftAngle(S)` | `cxx:56-61` | [`BRepOffsetApiDraftAngle::with_shape`] |
//! | `Clear()` | `cxx:64-70` | [`BRepOffsetApiDraftAngle::clear`] |
//! | `Init(S)` | `cxx:74-86` | [`BRepOffsetApiDraftAngle::init`] |
//! | `Add(F, D, Angle, Plane, Flag)` | `cxx:90-104` | [`BRepOffsetApiDraftAngle::add`] |
//! | `AddDone()` | `cxx:108-113` | [`BRepOffsetApiDraftAngle::add_done`] |
//! | `Remove(F)` | `cxx:117-122` | [`BRepOffsetApiDraftAngle::remove`] |
//! | `ProblematicShape()` | `cxx:126-133` | [`BRepOffsetApiDraftAngle::problematic_shape`] |
//! | `Status()` | `cxx:137-142` | [`BRepOffsetApiDraftAngle::status`] |
//! | `ConnectedFaces(F)` | `cxx:145-152` | [`BRepOffsetApiDraftAngle::connected_faces`] |
//! | `ModifiedFaces()` | `cxx:156-163` | [`BRepOffsetApiDraftAngle::modified_faces`] |
//! | `Generated(S)` | `cxx:167-192` | [`BRepOffsetApiDraftAngle::generated`] |
//! | `Modified(S)` | `cxx:196-227` | [`BRepOffsetApiDraftAngle::modified`] |
//! | `ModifiedShape(S)` | `cxx:232-250` | [`BRepOffsetApiDraftAngle::modified_shape`] |
//! | `Build(Range)` | `cxx:254-268` | [`BRepOffsetApiDraftAngle::build`] |
//! | `CorrectWires()` | `cxx:271-909` | `correct_wires` |
//! | `CorrectVertexTol()` | `cxx:913-993` | `correct_vertex_tol` |
//!
//! # Translation notes
//!
//! * `myModification` is an `occ::handle<Draft_Modification>` that the default
//!   constructor leaves null (`hxx:214`); the port keeps it as
//!   `Option<DraftModification>`. `Add` / `Remove` / `Build` dereference that
//!   handle without a null test in the reference, which would crash; the port
//!   answers the `Standard_NullObject` string instead (no OCCT rule added).
//! * `myInitialShape` is a `TopoDS_Shape` that can be null, so
//!   `IsNull()` is `Option::is_none()`. Every `Standard_NullObject_Raise_if`
//!   on it becomes the `Err` of the same message.
//! * `myModifier` is `BRepTools_Modifier` (`hxx:215`); `DoModif(S)` (from
//!   `BRepBuilderAPI_ModifyShape`) is `Modifier::on_shape(S)` + `perform` +
//!   `modified_shape(S)`.
//! * `mySubs` is a `BRepTools_ReShape` (`hxx:220`), modelled by
//!   [`MapReShape`]: `Value(S)` is [`MapReShape::value`] with the "not
//!   recorded answers `S`" fallback of the reference.
//! * `TopTools_DataMapOfShapeShape::IsEmpty` / `IsBound` / `Find` on
//!   `myVtxToReplace` (`hxx:219`) are an insertion-ordered
//!   `(key, value)` list keyed by `IsSame`, matching the map's iteration order
//!   that `CorrectVertexTol` relies on.
//! * `NotDone()` is `BRepBuilderAPI_Command::NotDone` (`hxx:54-55`): it only
//!   clears `myDone`.

mod correct_vertex_tol;
mod correct_wires;

use occt_core::gp::{GpDir, GpPln};

use crate::brep_tools_modifier::Modifier;
use crate::draft::{DraftErrorStatus, DraftModification};
use crate::shape::{Face, TopoShape};
use crate::shape_fix_compose_shell::MapReShape;
use crate::topo_tools::{is_same, shape_map_key};

/// `class BRepOffsetAPI_DraftAngle : public BRepBuilderAPI_ModifyShape`
/// (`hxx:38-221`).
pub struct BRepOffsetApiDraftAngle {
    /// `myInitialShape` (`hxx:212`): null until the shape-taking constructor or
    /// `Init`.
    pub(super) my_initial_shape: Option<TopoShape>,
    /// `myShape` (`BRepBuilderAPI_MakeShape`).
    pub(super) my_shape: Option<TopoShape>,
    /// `myGenerated` (`BRepBuilderAPI_MakeShape`).
    pub(super) my_generated: Vec<TopoShape>,
    /// `myDone` (`BRepBuilderAPI_Command`).
    my_done: bool,
    /// `myModification` (`hxx:214`).
    pub(super) my_modification: Option<DraftModification>,
    /// `myModifier` (`hxx:215`).
    pub(super) my_modifier: Modifier,
    /// `myVtxToReplace` (`hxx:219`), insertion ordered and keyed by `IsSame`.
    pub(super) my_vtx_to_replace: Vec<(TopoShape, TopoShape)>,
    /// `mySubs` (`hxx:220`).
    pub(super) my_subs: MapReShape,
}

impl Default for BRepOffsetApiDraftAngle {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiDraftAngle {
    /// `BRepOffsetAPI_DraftAngle()` (`cxx:52`): every member default
    /// constructed, so `myModification` and `myInitialShape` are null.
    pub fn new() -> Self {
        Self {
            my_initial_shape: None,
            my_shape: None,
            my_generated: Vec::new(),
            my_done: false,
            my_modification: None,
            my_modifier: Modifier::new(false),
            my_vtx_to_replace: Vec::new(),
            my_subs: MapReShape::new(),
        }
    }

    /// `BRepOffsetAPI_DraftAngle(S)` (`cxx:56-61`).
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            my_initial_shape: Some(s.clone()),
            my_modification: Some(DraftModification::new(s)),
            ..Self::new()
        }
    }

    /// `Clear()` (`cxx:63-70`).
    pub fn clear(&mut self) {
        if let Some(dmod) = self.my_modification.as_mut() {
            dmod.clear();
        }
    }

    /// `Init(S)` (`cxx:74-86`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_initial_shape = Some(s.clone());
        // `NotDone()` (`cxx:78`).
        self.not_done();
        match self.my_modification.as_mut() {
            Some(dmod) => dmod.init(s),
            None => self.my_modification = Some(DraftModification::new(s)),
        }
    }

    /// `Add(F, D, Angle, Plane, Flag = true)` (`cxx:91-105`).
    ///
    /// `Err` carries the reference's `Standard_NullObject` / the
    /// `Standard_ConstructionError` of `Draft_Modification::Add`.
    pub fn add(
        &mut self,
        f: &Face,
        d: &GpDir,
        angle: f64,
        plane: &GpPln,
        flag: bool,
    ) -> Result<(), String> {
        // `cxx:96-100`: POP-DPF protection, a silent no-op.
        if angle.abs() <= 1.0e-04 {
            return Ok(());
        }
        // `cxx:101-102`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::Add - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:104`
        let dmod = self.my_modification.as_mut().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::Add - modification is null".to_string()
        })?;
        dmod.add(f, d, angle, plane, flag).map(|_| ())
    }

    /// `AddDone()` (`cxx:109-113`): `ProblematicShape().IsNull()`.
    pub fn add_done(&self) -> Result<bool, String> {
        // `cxx:110-111`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::AddDone - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:112`
        let dmod = self.my_modification.as_ref().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::AddDone - modification is null"
                .to_string()
        })?;
        Ok(dmod.problematic_shape().is_none())
    }

    /// `Remove(F)` (`cxx:117-122`).
    pub fn remove(&mut self, f: &Face) -> Result<(), String> {
        // `cxx:118-120`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::Remove - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:121`
        let dmod = self.my_modification.as_mut().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::Remove - modification is null"
                .to_string()
        })?;
        dmod.remove(f)
    }

    /// `ProblematicShape()` (`cxx:126-133`).
    pub fn problematic_shape(&self) -> Result<Option<&TopoShape>, String> {
        // `cxx:128-130`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::ProblematicShape - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:131`
        let dmod = self.my_modification.as_ref().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::ProblematicShape - modification is null"
                .to_string()
        })?;
        Ok(dmod.problematic_shape())
    }

    /// `Status()` (`cxx:137-142`): the `Draft_ErrorStatus` of the modification.
    pub fn status(&self) -> Result<DraftErrorStatus, String> {
        // `cxx:138-140`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::Status - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:141`
        let dmod = self.my_modification.as_ref().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::Status - modification is null".to_string()
        })?;
        Ok(dmod.error())
    }

    /// `ConnectedFaces(F)` (`cxx:145-152`).
    pub fn connected_faces(&mut self, f: &Face) -> Result<&[TopoShape], String> {
        // `cxx:150-151`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::ConnectedFaces - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:152`
        let dmod = self.my_modification.as_mut().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::ConnectedFaces - modification is null"
                .to_string()
        })?;
        dmod.connected_faces(f)
    }

    /// `ModifiedFaces()` (`cxx:156-163`).
    pub fn modified_faces(&mut self) -> Result<&[TopoShape], String> {
        // `cxx:160-161`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::ModifiedFaces - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:162`
        let dmod = self.my_modification.as_mut().ok_or_else(|| {
            "Standard_NullObject: BRepOffsetAPI_DraftAngle::ModifiedFaces - modification is null"
                .to_string()
        })?;
        dmod.modified_faces()
    }

    /// `Generated(S)` (`cxx:167-192`).
    pub fn generated(&mut self, s: &TopoShape) -> Result<&[TopoShape], String> {
        // `cxx:168`
        self.my_generated.clear();
        // `cxx:169-171`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::Generated - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:172-174`
        if s.shape_type() == crate::abs::ShapeType::Face {
            // `cxx:175-179`
            if self.is_new_surface(s) {
                // `cxx:180-189`
                let m = self.modified_shape(s)?;
                if self.my_vtx_to_replace.is_empty() {
                    self.my_generated.push(m);
                } else {
                    // `mySubs.Value(...)`
                    let v = self.my_subs.value(&m).unwrap_or(m);
                    self.my_generated.push(v);
                }
            }
        }
        Ok(&self.my_generated)
    }

    /// `Modified(S)` (`cxx:196-227`).
    pub fn modified(&mut self, s: &TopoShape) -> Result<&[TopoShape], String> {
        // `cxx:198`
        self.my_generated.clear();
        // `cxx:199-201`
        if self.my_initial_shape.is_none() {
            return Err(
                "Standard_NullObject: BRepOffsetAPI_DraftAngle::Modified - initial shape is not set"
                    .to_string(),
            );
        }
        // `cxx:204`
        if s.shape_type() == crate::abs::ShapeType::Face {
            // `cxx:210`: not a generation, so possibly a modification.
            if !self.is_new_surface(s) {
                let m = self.modified_shape(s)?;
                if self.my_vtx_to_replace.is_empty() {
                    self.my_generated.push(m);
                } else {
                    let v = self.my_subs.value(&m).unwrap_or(m);
                    self.my_generated.push(v);
                }
                // `cxx:222-225`
                if self.my_generated.len() == 1 && is_same(&self.my_generated[0], s) {
                    self.my_generated.clear();
                }
            }
        }
        Ok(&self.my_generated)
    }

    /// `ModifiedShape(S)` (`cxx:231-250`).
    pub fn modified_shape(&self, s: &TopoShape) -> Result<TopoShape, String> {
        // `cxx:235-241`
        if s.shape_type() == crate::abs::ShapeType::Vertex {
            if let Some(v) = self
                .my_vtx_to_replace
                .iter()
                .find(|(k, _)| is_same(k, s))
                .map(|(_, v)| v.clone())
            {
                return Ok(v);
            }
        }
        // `cxx:242-245`
        if self.my_vtx_to_replace.is_empty() {
            return self.my_modifier.modified_shape(s);
        }
        // `cxx:247-249`
        let a_ns = self.my_modifier.modified_shape(s)?;
        Ok(self.my_subs.value(&a_ns).unwrap_or(a_ns))
    }

    /// `Build(Range)` (`cxx:254-268`).
    pub fn build(&mut self) -> Result<(), String> {
        // `cxx:256-257`
        if let Some(dmod) = self.my_modification.as_mut() {
            // `Perform()` is `void` in the reference; the port's `Err` records
            // additive-algorithm failures that `Status()` reports instead.
            let _ = dmod.perform();
        }
        // `cxx:257-259`
        let done = self
            .my_modification
            .as_ref()
            .map(|dmod| dmod.is_done())
            .unwrap_or(false);
        if !done {
            self.not_done();
            return Ok(());
        }
        // `cxx:262-264`
        self.do_modif()?;
        self.correct_wires();
        self.correct_vertex_tol();
        Ok(())
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`): the built shape.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `IsDone()` (`BRepBuilderAPI_Command`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `NotDone()` (`BRepBuilderAPI_Command.hxx:54-55`).
    fn not_done(&mut self) {
        self.my_done = false;
    }

    /// `BRepBuilderAPI_ModifyShape::DoModif(S)` (`BRepBuilderAPI_ModifyShape.cxx`):
    /// `if (!S.IsEqual(myInitialShape) || !IsDone()) { myInitialShape = S;
    /// myModifier.Init(S); DoModif(); }`, where `DoModif()` runs
    /// `myModifier.Perform(myModification)`, and on `myModifier.IsDone()` calls
    /// `Done()` and stores `myShape = myModifier.ModifiedShape(myInitialShape)`.
    ///
    /// `S` is always the facade's own `myInitialShape`, so `IsEqual` holds and
    /// the guard reduces to `!IsDone()`.
    fn do_modif(&mut self) -> Result<(), String> {
        if self.my_done {
            return Ok(());
        }
        let Some(initial) = self.my_initial_shape.clone() else {
            return Ok(());
        };
        let Some(dmod) = self.my_modification.as_ref() else {
            return Ok(());
        };
        let mut modifier = Modifier::on_shape(&initial);
        modifier.perform(dmod)?;
        if modifier.is_done() {
            // `Done()`
            self.my_done = true;
            self.my_shape = Some(modifier.modified_shape(&initial)?);
            self.my_modifier = modifier;
        }
        Ok(())
    }

    /// `DMod->NewSurface(F, Surf, L, Tol, RW, RF)` as a boolean test
    /// (`cxx:179`, `cxx:210`).
    ///
    /// The reference ignores the out parameters here; an `Err` (the
    /// `Standard_DomainError` of `IsDone()`) is not a "false" in OCCT, so it is
    /// propagated as `false` only after `IsDone()` - `Generated` / `Modified`
    /// are queried after `Build`, where the throw cannot trigger.
    fn is_new_surface(&self, s: &TopoShape) -> bool {
        let Some(dmod) = self.my_modification.as_ref() else {
            return false;
        };
        let Some(f) = Face::wrap(s.clone()) else {
            return false;
        };
        let mut surface = None;
        let mut tol = 0.0;
        let (mut rw, mut rf) = (false, false);
        matches!(
            dmod.new_surface_of(&f, &mut surface, &mut tol, &mut rw, &mut rf),
            Ok(true)
        )
    }

    /// `myVtxToReplace.IsBound(V)` / `myVtxToReplace(V)` (`cxx:237-239`).
    pub(super) fn vtx_replacement(&self, v: &TopoShape) -> Option<TopoShape> {
        self.my_vtx_to_replace
            .iter()
            .find(|(k, _)| is_same(k, v))
            .map(|(_, val)| val.clone())
    }

    /// `myVtxToReplace.Bind(V, NewV)` (`cxx:967`).
    pub(super) fn bind_vtx_replacement(&mut self, v: &TopoShape, new_v: &TopoShape) {
        let key = shape_map_key(v);
        for entry in self.my_vtx_to_replace.iter_mut() {
            if shape_map_key(&entry.0) == key {
                entry.1 = new_v.clone();
                return;
            }
        }
        self.my_vtx_to_replace.push((v.clone(), new_v.clone()));
    }
}
