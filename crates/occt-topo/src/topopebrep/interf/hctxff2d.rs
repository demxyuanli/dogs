//! `TopOpeBRep_Hctxff2d` — the face/face intersection context that
//! `TopOpeBRep_FacesIntersector` hands to `TopOpeBRep_Point2d`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_Hctxff2d.hxx:27-68`,
//! `.cxx:28-245`.
//!
//! Translation notes:
//! - `occ::handle<BRepAdaptor_Surface>` becomes [`SurfaceAdaptor`] (the port's
//!   `BRepAdaptor_Surface` slice); `mySurface1.IsNull()` is `Option::is_none`.
//! - `SetFaces`'s `mySurface1 = new BRepAdaptor_Surface; ... Initialize(F1,
//!   computerestriction)` (`cxx:55-65`) is exactly `SurfaceAdaptor::of_face`
//!   with `Restriction == false`; the existing adaptor is overwritten instead
//!   of re-initialized, which reads the same state afterwards.
//! - `TopoDS_Face::IsEqual` is `crate::topo_tools::is_equal`.
//! - `Face(Index)` / `HSurface(Index)` raise `Standard_Failure` for an index
//!   other than 1 or 2 (`cxx:191-195`, `:205-209`); the port returns `Err`.

use crate::abs::Orientation;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::Face;
use crate::topo_tools::is_equal;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::shape_tool::{faces_same_oriented, surfaces_same_oriented};

/// `class TopOpeBRep_Hctxff2d` (`hxx:27-68`).
#[derive(Clone)]
pub struct Hctxff2d {
    /// `myFace1`
    face1: Option<Face>,
    /// `mySurface1`
    surface1: Option<SurfaceAdaptor>,
    /// `mySurfaceType1` — assigned by `SetHSurfacesPrivate` (`cxx:104`). OCCT
    /// exposes no accessor for it, so it is write-only here as well.
    #[allow(dead_code)]
    surface_type1: GeomAbsSurfaceType,
    /// `myf1surf1F_sameoriented`
    f1surf1f_same_oriented: bool,
    /// `myFace2`
    face2: Option<Face>,
    /// `mySurface2`
    surface2: Option<SurfaceAdaptor>,
    /// `mySurfaceType2` — assigned by `SetHSurfacesPrivate` (`cxx:108`);
    /// write-only, as in OCCT.
    #[allow(dead_code)]
    surface_type2: GeomAbsSurfaceType,
    /// `myf2surf1F_sameoriented`
    f2surf1f_same_oriented: bool,
    /// `mySurfacesSameOriented`
    surfaces_same_oriented: bool,
    /// `myFacesSameOriented`
    faces_same_oriented: bool,
    /// `myTol1`
    tol1: f64,
    /// `myTol2`
    tol2: f64,
}

impl Default for Hctxff2d {
    fn default() -> Self {
        Self::new()
    }
}

impl Hctxff2d {
    /// `TopOpeBRep_Hctxff2d()` (`cxx:28-37`).
    pub fn new() -> Self {
        Self {
            face1: None,
            surface1: None,
            surface_type1: GeomAbsSurfaceType::OtherSurface,
            f1surf1f_same_oriented: true,
            face2: None,
            surface2: None,
            surface_type2: GeomAbsSurfaceType::OtherSurface,
            f2surf1f_same_oriented: true,
            surfaces_same_oriented: false,
            faces_same_oriented: false,
            tol1: 0.0,
            tol2: 0.0,
        }
    }

    /// `SetFaces(F1, F2)` (`cxx:41-70`).
    pub fn set_faces(&mut self, f1: &Face, f2: &Face) {
        let newf1 = match &self.face1 {
            Some(face) => !is_equal(face, f1),
            None => true,
        };
        let newf2 = match &self.face2 {
            Some(face) => !is_equal(face, f2),
            None => true,
        };
        let yaduneuf = newf1 || newf2;
        if !yaduneuf {
            return;
        }

        // `computerestriction` is `false` in this method (`cxx:51-54`).
        if newf1 {
            self.surface1 = Some(SurfaceAdaptor::of_face(f1));
        }
        if newf2 {
            self.surface2 = Some(SurfaceAdaptor::of_face(f2));
        }
        self.set_h_surfaces_private();
    }

    /// `SetHSurfaces(HS1, HS2)` (`cxx:74-95`).
    pub fn set_h_surfaces(&mut self, s1: Option<&SurfaceAdaptor>, s2: Option<&SurfaceAdaptor>) {
        let mut newf1 = false;
        let mut newf2 = false;
        if let Some(adaptor) = s1 {
            newf1 = match (&self.face1, adaptor.face()) {
                (Some(face), Some(other)) => !is_equal(face, other),
                (None, Some(_)) => true,
                _ => false,
            };
        }
        if let Some(adaptor) = s2 {
            newf2 = match (&self.face2, adaptor.face()) {
                (Some(face), Some(other)) => !is_equal(face, other),
                (None, Some(_)) => true,
                _ => false,
            };
        }
        let yaduneuf = newf1 || newf2;
        if !yaduneuf {
            return;
        }

        self.surface1 = s1.cloned();
        self.surface2 = s2.cloned();
        self.set_h_surfaces_private();
    }

    /// `SetHSurfacesPrivate()` (`cxx:99-150`).
    fn set_h_surfaces_private(&mut self) {
        // OCCT dereferences both handles here; `set_faces` / `set_h_surfaces`
        // only reach this point after filling them for the arms that changed.
        if let Some(adaptor) = self.surface1.as_ref() {
            self.face1 = adaptor.face().cloned();
            self.surface_type1 = adaptor.get_type();
        }
        if let Some(adaptor) = self.surface2.as_ref() {
            self.face2 = adaptor.face().cloned();
            self.surface_type2 = adaptor.get_type();
        }

        self.surfaces_same_oriented = true;
        self.faces_same_oriented = true;
        self.f1surf1f_same_oriented = true;
        self.f2surf1f_same_oriented = true;

        let face1forward = match self.face1.as_ref() {
            Some(face) => Some(face.oriented(Orientation::Forward)),
            None => None,
        };
        let so11 = match (face1forward.as_ref(), self.face1.as_ref()) {
            (Some(f1f), Some(f1)) => faces_same_oriented(f1f, f1),
            _ => true,
        };
        self.f1surf1f_same_oriented = so11;
        let so21 = match (face1forward.as_ref(), self.face2.as_ref()) {
            (Some(f1f), Some(f2)) => faces_same_oriented(f1f, f2),
            _ => true,
        };
        self.f2surf1f_same_oriented = so21;

        if let (Some(s1), Some(s2)) = (self.surface1.as_ref(), self.surface2.as_ref()) {
            self.surfaces_same_oriented = surfaces_same_oriented(s1, s2);
        }
        if let (Some(f1), Some(f2)) = (self.face1.as_ref(), self.face2.as_ref()) {
            self.faces_same_oriented = faces_same_oriented(f1, f2);
        }
    }

    /// `SetTolerances(Tol1, Tol2)` (`cxx:153-157`).
    pub fn set_tolerances(&mut self, tol1: f64, tol2: f64) {
        self.tol1 = tol1;
        self.tol2 = tol2;
    }

    /// `GetTolerances(Tol1, Tol2)` (`cxx:161-165`).
    pub fn get_tolerances(&self) -> (f64, f64) {
        (self.tol1, self.tol2)
    }

    /// `GetMaxTolerance()` (`cxx:170-174`).
    pub fn get_max_tolerance(&self) -> f64 {
        self.tol1.max(self.tol2)
    }

    /// `Face(Index)` (`cxx:177-195`).
    pub fn face(&self, index: i32) -> Result<&Face, String> {
        match index {
            1 => self.face1.as_ref().ok_or_else(|| "TopOpeBRep_Hctxff2d::Face".into()),
            2 => self.face2.as_ref().ok_or_else(|| "TopOpeBRep_Hctxff2d::Face".into()),
            _ => Err("TopOpeBRep_Hctxff2d::Face".into()),
        }
    }

    /// `HSurface(Index)` (`cxx:197-212`).
    pub fn h_surface(&self, index: i32) -> Result<&SurfaceAdaptor, String> {
        match index {
            1 => self
                .surface1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxff2d::HSurface".into()),
            2 => self
                .surface2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxff2d::HSurface".into()),
            _ => Err("TopOpeBRep_Hctxff2d::HSurface".into()),
        }
    }

    /// `SurfacesSameOriented()` (`cxx:216-219`).
    pub fn surfaces_same_oriented(&self) -> bool {
        self.surfaces_same_oriented
    }

    /// `FacesSameOriented()` (`cxx:222-225`).
    pub fn faces_same_oriented(&self) -> bool {
        self.faces_same_oriented
    }

    /// `FaceSameOrientedWithRef(Index)` (`cxx:228-245`).
    pub fn face_same_oriented_with_ref(&self, index: i32) -> Result<bool, String> {
        match index {
            1 => Ok(self.f1surf1f_same_oriented),
            2 => Ok(self.f2surf1f_same_oriented),
            _ => Err("TopOpeBRep_Hctxff2d::FSO".into()),
        }
    }
}
