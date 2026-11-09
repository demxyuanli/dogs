//! `Draft_FaceInfo` (`Draft_FaceInfo.cxx:24-117`): the new surface of a face the
//! draft modification rebuilds, plus the two faces that share it.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Draft_FaceInfo()` | `cxx:24` | [`DraftFaceInfo::new`] |
//! | `Draft_FaceInfo(S, HasNewGeometry)` | `cxx:28-41` | [`DraftFaceInfo::with_surface`] |
//! | `RootFace(F)` | `cxx:44-47` | [`DraftFaceInfo::root_face`] |
//! | `Add(F)` | `cxx:52-62` | [`DraftFaceInfo::add`] |
//! | `FirstFace()` / `SecondFace()` | `cxx:66-76` | [`DraftFaceInfo::first_face`] / [`DraftFaceInfo::second_face`] |
//! | `NewGeometry()` | `cxx:81-84` | [`DraftFaceInfo::new_geometry`] |
//! | `Geometry()` / `ChangeGeometry()` | `cxx:87-96` | [`DraftFaceInfo::geometry`] / [`DraftFaceInfo::change_geometry`] |
//! | `Curve()` / `ChangeCurve()` | `cxx:100-110` | [`DraftFaceInfo::curve`] / [`DraftFaceInfo::change_curve`] |
//! | `RootFace()` | `cxx:114-117` | [`DraftFaceInfo::root_face_value`] |
//!
//! # Translation notes
//!
//! * The constructor unwraps a `Geom_RectangularTrimmedSurface` to its basis
//!   (`cxx:31-39`); [`occt_geom::Surface::rectangular_trimmed_basis`] is that
//!   `DynamicType()` test plus `BasisSurface()`. The unwrapped surface is stored
//!   untrimmed, which is what `NewSurface` answers afterwards.
//! * `Add(F)` has no `IsSame` test (`cxx:52-62`): the second slot simply takes
//!   the first face offered when the first slot is already used.
//! * `myCurv` (`hxx:65`) is the circle `InternalAdd` keeps for a postponed
//!   cylindrical / extruded face; it is filled by `ChangeCurve`.

use std::sync::Arc;

use occt_geom::{Curve, Surface};

use crate::shape::Face;

/// `class Draft_FaceInfo` (`Draft_FaceInfo.hxx:29-67`).
#[derive(Clone)]
pub struct DraftFaceInfo {
    /// `myNewGeom` (`hxx:59`).
    my_new_geom: bool,
    /// `myGeom` (`hxx:60`).
    my_geom: Option<Arc<dyn Surface>>,
    /// `myRootFace` (`hxx:61`).
    my_root_face: Option<Face>,
    /// `myF1` (`hxx:62`).
    my_f1: Option<Face>,
    /// `myF2` (`hxx:63`).
    my_f2: Option<Face>,
    /// `myCurv` (`hxx:65`).
    my_curv: Option<Arc<dyn Curve>>,
}

impl Default for DraftFaceInfo {
    fn default() -> Self {
        Self::new()
    }
}

impl DraftFaceInfo {
    /// `Draft_FaceInfo()` (`cxx:24`): the default constructor leaves every
    /// handle null and `myNewGeom` uninitialized in OCCT; the port answers
    /// `false`, the value `Draft_Modification` would see for a face it never
    /// rebuilt.
    pub fn new() -> Self {
        Self {
            my_new_geom: false,
            my_geom: None,
            my_root_face: None,
            my_f1: None,
            my_f2: None,
            my_curv: None,
        }
    }

    /// `Draft_FaceInfo(S, HasNewGeometry)` (`cxx:28-41`).
    ///
    /// `s` is the `occ::handle<Geom_Surface>&` argument: `None` is the null
    /// handle `InternalAdd` stores for a postponed (cylinder / extrusion)
    /// face (`Draft_Modification_1.cxx:266`), which leaves `myGeom` null.
    pub fn with_surface(s: Option<&Arc<dyn Surface>>, has_new_geometry: bool) -> Self {
        let mut me = Self::new();
        me.my_new_geom = has_new_geometry;
        // `cxx:31-39`: `occ::down_cast<Geom_RectangularTrimmedSurface>(S)`.
        me.my_geom = s.map(|s| match s.rectangular_trimmed_basis() {
            Some(basis) => basis,
            None => s.clone(),
        });
        me
    }

    /// `RootFace(F)` (`cxx:44-47`).
    pub fn root_face(&mut self, f: &Face) {
        self.my_root_face = Some(f.clone());
    }

    /// `Add(F)` (`cxx:52-62`).
    pub fn add(&mut self, f: &Face) {
        if self.my_f1.is_none() {
            self.my_f1 = Some(f.clone());
        } else if self.my_f2.is_none() {
            self.my_f2 = Some(f.clone());
        }
    }

    /// `FirstFace()` (`cxx:66-69`).
    pub fn first_face(&self) -> Option<&Face> {
        self.my_f1.as_ref()
    }

    /// `SecondFace()` (`cxx:72-76`).
    pub fn second_face(&self) -> Option<&Face> {
        self.my_f2.as_ref()
    }

    /// `NewGeometry()` (`cxx:81-84`).
    pub fn new_geometry(&self) -> bool {
        self.my_new_geom
    }

    /// `Geometry()` (`cxx:87-90`).
    pub fn geometry(&self) -> Option<&Arc<dyn Surface>> {
        self.my_geom.as_ref()
    }

    /// `ChangeGeometry()` (`cxx:93-96`).
    pub fn change_geometry(&mut self) -> &mut Option<Arc<dyn Surface>> {
        &mut self.my_geom
    }

    /// `Curve()` (`cxx:100-103`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.my_curv.as_ref()
    }

    /// `ChangeCurve()` (`cxx:107-110`).
    pub fn change_curve(&mut self) -> &mut Option<Arc<dyn Curve>> {
        &mut self.my_curv
    }

    /// `RootFace()` (`cxx:114-117`).
    pub fn root_face_value(&self) -> Option<&Face> {
        self.my_root_face.as_ref()
    }
}
