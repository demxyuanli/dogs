//! `BRepFill_MultiLine` — the 3D curve and the two p-curves of the
//! intersection between the ruled surface of a bisector and two faces.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_MultiLine.cxx:1-787`
//! (OCCT 8.0.0), `.hxx:31-120`.
//!
//! The class is an `AppCont_Function` handed to the approximation of the sweep
//! (`BRepFill_Sweep`). Every evaluation is done in the parametric space of
//! `Face1` / `Face2`: for a parameter `U` of the bisector it answers
//!
//! * `Value(U)` — the 3D point of `Face1` surface at the (u, v) found on the
//!   face (`cxx:566-580`),
//! * `ValueOnF1(U)` / `ValueOnF2(U)` — the (u, v) of that point on each face
//!   (`cxx:728-739`), solved by [`value_on_face`],
//! * `Curves(...)` — the analytic 3D curve + p-curves when the configuration
//!   is a particular case (`cxx:405-548`), i.e. when the projection of the
//!   bisector on the faces stays on an iso line.
//!
//! Ported layout:
//! * this file — the type, its default constructor, the queries and the
//!   `AppCont_Function` implementation,
//! * [`build`] — the seven-argument constructor (`cxx:78-391`),
//! * [`value_on_face`] — the static `ValueOnFace` solver (`cxx:584-724`),
//! * [`particular_case`] — `Curves` (`cxx:405-548`).
//!
//! UNPORTED: the `myKPart == 2` branch of `Curves` (`cxx:509-547`) calls
//! `GeomProjLib::ProjectOnPlane` (`ProjLib_ProjectOnPlane.cxx:1567` lines),
//! which has no port yet; [`MultiLine::curves`] reports
//! [`Error::UnportedProjectOnPlane`] for it instead of re-deriving the
//! projection with a different control flow.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom::curve::Curve;

use crate::shape::{Edge, Face};

mod build;
mod particular_case;
pub mod value_on_face;

/// Failure modes of [`MultiLine`].
///
/// OCCT signals these with `Standard_ConstructionError` / a null handle; the
/// port has no exceptions, so the same conditions are reported as values. The
/// `.cxx` line of each OCCT throw is named on the variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// `BRepFill_MultiLine : Edge without PCurve` (`cxx:60-62`).
    EdgeWithoutPCurve,
    /// `BRep_Tool::Surface` answered null for a face (`cxx:146`, `cxx:263`).
    NullSurface,
    /// `GeomProjLib::Curve2d` answered a null handle (`GeomProjLib.cxx:108-111`,
    /// reached from `cxx:227`, `cxx:234`, `cxx:346`, `cxx:350`).
    NullProjection,
    /// `Standard_ConstructionError` of `gp_Ax3(P, DZ, gp::DZ())` when the iso
    /// tangent is parallel to the global Z (`cxx:229-231`, `cxx:348-350`).
    ConstructionError,
    /// `Geom_TrimmedCurve` refused the trim of `cxx:433` / `cxx:457`
    /// (`Geom_TrimmedCurve.cxx:99-151`).
    TrimError,
    /// `BRepFill_MultiLine: ValueOnFace` (`cxx:629`).
    ValueOnFace,
    /// `Curves` is only defined when `IsParticularCase()` is true
    /// (`hxx:72`).
    NotParticularCase,
    /// `Geom_TrimmedCurve` / `Geom2d_TrimmedCurve` are built on a null curve
    /// (`cxx:428`, `cxx:455`), or a stored adaptor was never loaded.
    NullCurve,
    /// `Curves` on the `myKPart == 2` branch (`cxx:509-547`) needs
    /// `GeomProjLib::ProjectOnPlane`, which is UNPORTED (see the module doc).
    UnportedProjectOnPlane,
}

/// `class BRepFill_MultiLine : public AppCont_Function` (`hxx:46-120`).
#[derive(Clone)]
pub struct MultiLine {
    /// `myFace1` (`hxx:109`).
    face1: Face,
    /// `myFace2` (`hxx:110`). Only the `myKPart == 2` branch of `Curves`
    /// (`cxx:536-544`) reads it, and that branch is UNPORTED, so the port keeps
    /// the member without a reader yet.
    #[allow(dead_code)]
    face2: Face,
    /// `myU1` (`hxx:111`): the projected u-iso of `myFace1`.
    u1: Geom2dAdaptorCurve,
    /// `myV1` (`hxx:112`): the projected v-iso of `myFace1`.
    v1: Geom2dAdaptorCurve,
    /// `myU2` (`hxx:113`).
    u2: Geom2dAdaptorCurve,
    /// `myV2` (`hxx:114`).
    v2: Geom2dAdaptorCurve,
    /// `myIsoU1` (`hxx:115`): `Edge1` is a `GeomAbs_IsoU` edge of `myFace1`.
    iso_u1: bool,
    /// `myIsoU2` (`hxx:116`).
    iso_u2: bool,
    /// `myBis` (`hxx:117`): the bisector.
    bis: Geom2dAdaptorCurve,
    /// `myKPart` (`hxx:118`): `0` = no particular case, `1` = iso/iso,
    /// `2` = two planes.
    k_part: i32,
    /// `myCont` (`hxx:119`): continuity between the two faces along the
    /// bisector.
    cont: Shape,
}

impl Default for MultiLine {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiLine {
    /// `BRepFill_MultiLine()` (`cxx:71-75`): an unloaded function whose
    /// `AppCont_Function` counts are set (`myNbPnt2d = 2`, `myNbPnt = 1`) while
    /// the faces and adaptors stay default.
    pub fn new() -> Self {
        Self {
            face1: Face::new(),
            face2: Face::new(),
            u1: Geom2dAdaptorCurve::new(),
            v1: Geom2dAdaptorCurve::new(),
            u2: Geom2dAdaptorCurve::new(),
            v2: Geom2dAdaptorCurve::new(),
            iso_u1: false,
            iso_u2: false,
            bis: Geom2dAdaptorCurve::new(),
            k_part: 0,
            cont: Shape::C0,
        }
    }

    /// `BRepFill_MultiLine(Face1, Face2, Edge1, Edge2, Inv1, Inv2, Bissec)`
    /// (`cxx:78-391`).
    ///
    /// `inv1` / `inv2` are the `Inv1` / `Inv2` flags of the constructor: they
    /// say that the first / second iso of the face must be taken at its upper
    /// bound.
    pub fn from_faces(
        face1: Face,
        face2: Face,
        edge1: &Edge,
        edge2: &Edge,
        inv1: bool,
        inv2: bool,
        bissec: Arc<dyn Curve2d>,
    ) -> Result<Self, Error> {
        build::build(face1, face2, edge1, edge2, inv1, inv2, bissec)
    }

    /// `IsParticularCase()` (`cxx:395-397`).
    pub fn is_particular_case(&self) -> bool {
        self.k_part != 0
    }

    /// `Continuity()` (`cxx:761-764`).
    pub fn continuity(&self) -> Shape {
        self.cont
    }

    /// `FirstParameter()` (`cxx:552-555`): the bisector's first parameter.
    pub fn first_parameter(&self) -> f64 {
        self.bis.first_parameter()
    }

    /// `LastParameter()` (`cxx:558-561`).
    pub fn last_parameter(&self) -> f64 {
        self.bis.last_parameter()
    }

    /// `Value(U)` (`cxx:566-580`): the 3D point of `myFace1` at the (u, v)
    /// [`Self::value_on_f1`] answers.
    pub fn value(&self, u: f64) -> Result<GpPnt, Error> {
        let p2d = self.value_on_f1(u)?;
        let s = crate::brep_tool::BRepTool::face_surface_world(&self.face1)
            .ok_or(Error::NullSurface)?;
        Ok(s.d0(p2d.x(), p2d.y()))
    }

    /// `ValueOnF1(U)` (`cxx:728-731`).
    pub fn value_on_f1(&self, u: f64) -> Result<GpPnt2d, Error> {
        value_on_face::value_on_face(u, &self.bis, &self.u1, &self.v1, self.iso_u1)
    }

    /// `ValueOnF2(U)` (`cxx:734-737`).
    pub fn value_on_f2(&self, u: f64) -> Result<GpPnt2d, Error> {
        value_on_face::value_on_face(u, &self.bis, &self.u2, &self.v2, self.iso_u2)
    }

    /// `Value3dOnF1OnF2(U, P3d, PF1, PF2)` (`cxx:745-757`).
    pub fn value_3d_on_f1_on_f2(
        &self,
        u: f64,
    ) -> Result<(GpPnt, GpPnt2d, GpPnt2d), Error> {
        let pf1 = self.value_on_f1(u)?;
        let pf2 = self.value_on_f2(u)?;
        let s = crate::brep_tool::BRepTool::face_surface_world(&self.face1)
            .ok_or(Error::NullSurface)?;
        Ok((s.d0(pf1.x(), pf1.y()), pf1, pf2))
    }

    /// `Curves(Curve, PCurve1, PCurve2)` (`cxx:405-548`): the analytic 3D
    /// curve and its two p-curves, defined only for a particular case
    /// (`hxx:72`).
    pub fn curves(&self) -> Result<(Arc<dyn Curve>, Arc<dyn Curve2d>, Arc<dyn Curve2d>), Error> {
        particular_case::curves(self)
    }

    /// `myKPart` (`hxx:118`), exposed for the sweep's diagnostics.
    pub fn k_part(&self) -> i32 {
        self.k_part
    }

    /// `myIsoU1` (`hxx:115`).
    pub fn is_iso_u1(&self) -> bool {
        self.iso_u1
    }

    /// `myIsoU2` (`hxx:116`).
    pub fn is_iso_u2(&self) -> bool {
        self.iso_u2
    }

    /// `myBis` (`hxx:117`).
    pub fn bissec(&self) -> &Geom2dAdaptorCurve {
        &self.bis
    }

    /// The `myNbPnt` / `myNbPnt2d` pair the constructor sets (`cxx:98-99`,
    /// `cxx:287-288`): one 3D point, two 2D points.
    pub fn number_of_points(&self) -> (i32, i32) {
        (1, 2)
    }
}

impl occt_geom::app_cont::Function for MultiLine {
    /// `GetNumberOfPoints` (`AppCont_Function.hxx:37-41`) on the pair set by
    /// the constructor (`cxx:98-99`).
    fn get_number_of_points(&self) -> (i32, i32) {
        self.number_of_points()
    }

    fn first_parameter(&self) -> f64 {
        Self::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        Self::last_parameter(self)
    }

    /// `Value(theT, thePnt2d, thePnt)` (`cxx:768-777`): `thePnt(1) = Value`,
    /// `thePnt2d(1) = ValueOnF1`, `thePnt2d(2) = ValueOnF2`.
    ///
    /// OCCT answers `true` and lets the geometry calls throw; the port answers
    /// `false` when one of them fails, which is what a caller of
    /// `AppCont_Function` can act on.
    fn value(&self, the_t: f64, the_pnt2d: &mut [GpPnt2d], the_pnt: &mut [GpPnt]) -> bool {
        let (Ok(p3d), Ok(pf1), Ok(pf2)) = (
            self.value(the_t),
            self.value_on_f1(the_t),
            self.value_on_f2(the_t),
        ) else {
            return false;
        };
        if let Some(slot) = the_pnt.first_mut() {
            *slot = p3d;
        }
        if let Some(slot) = the_pnt2d.first_mut() {
            *slot = pf1;
        }
        if let Some(slot) = the_pnt2d.get_mut(1) {
            *slot = pf2;
        }
        true
    }

    /// `D1(theT, theVec2d, theVec)` (`cxx:780-786`): `return false;`.
    fn d1(&self, _the_t: f64, _the_vec2d: &mut [GpVec2d], _the_vec: &mut [GpVec]) -> bool {
        false
    }
}
