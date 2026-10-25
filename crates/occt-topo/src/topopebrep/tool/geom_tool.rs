//! `TopOpeBRepTool_GeomTool` — the flags and tolerances that drive
//! `TopOpeBRepTool_CurveTool::MakeCurves`.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_GeomTool.hxx:27-89`,
//! `.cxx:22-141`.

use occt_core::precision;

use super::out_curve_type::OutCurveType;

/// `class TopOpeBRepTool_GeomTool` (`hxx:27-89`).
///
/// `CompC3D` / `CompPC1` / `CompPC2` say whether the result curves of
/// `MakeCurves` must be computed from the intersection line; `TypeC3D` is the
/// kind of 3d curve built when the line is a walking one (`hxx:34-43`).
#[derive(Debug, Clone, Copy)]
pub struct GeomTool {
    /// `myTypeC3D`
    type_c3d: OutCurveType,
    /// `myCompC3D`
    comp_c3d: bool,
    /// `myCompPC1`
    comp_pc1: bool,
    /// `myCompPC2`
    comp_pc2: bool,
    /// `myTol3d`
    tol3d: f64,
    /// `myTol2d`
    tol2d: f64,
    /// `myNbPntMax`
    nb_pnt_max: i32,
}

impl Default for GeomTool {
    /// `TopOpeBRepTool_GeomTool(TypeC3D = BSPLINE1, CompC3D = true,
    /// CompPC1 = true, CompPC2 = true)` (`hxx:44-47`, `.cxx:26-35`).
    fn default() -> Self {
        Self {
            type_c3d: OutCurveType::BSpline1,
            comp_c3d: true,
            comp_pc1: true,
            comp_pc2: true,
            tol3d: precision::APPROXIMATION,
            tol2d: precision::PAPPROXIMATION,
            nb_pnt_max: 30,
        }
    }
}

impl GeomTool {
    /// `TopOpeBRepTool_GeomTool(TypeC3D, CompC3D, CompPC1, CompPC2)`
    /// (`hxx:44-47`, `.cxx:26-35`).
    pub fn new(type_c3d: OutCurveType, comp_c3d: bool, comp_pc1: bool, comp_pc2: bool) -> Self {
        Self {
            type_c3d,
            comp_c3d,
            comp_pc1,
            comp_pc2,
            ..Self::default()
        }
    }

    /// `Define(TypeC3D, CompC3D, CompPC1, CompPC2)` (`cxx:38-48`).
    pub fn define(
        &mut self,
        type_c3d: OutCurveType,
        comp_c3d: bool,
        comp_pc1: bool,
        comp_pc2: bool,
    ) {
        self.type_c3d = type_c3d;
        self.comp_c3d = comp_c3d;
        self.comp_pc1 = comp_pc1;
        self.comp_pc2 = comp_pc2;
    }

    /// `Define(TypeC3D)` (`cxx:52-55`): only the curve type is written.
    pub fn define_type(&mut self, type_c3d: OutCurveType) {
        self.type_c3d = type_c3d;
    }

    /// `DefineCurves(CompC3D)` (`cxx:59-62`).
    pub fn define_curves(&mut self, comp_c3d: bool) {
        self.comp_c3d = comp_c3d;
    }

    /// `DefinePCurves1(CompPC1)` (`cxx:66-69`).
    pub fn define_pcurves1(&mut self, comp_pc1: bool) {
        self.comp_pc1 = comp_pc1;
    }

    /// `DefinePCurves2(CompPC2)` (`cxx:73-76`).
    pub fn define_pcurves2(&mut self, comp_pc2: bool) {
        self.comp_pc2 = comp_pc2;
    }

    /// `Define(const TopOpeBRepTool_GeomTool& GT)` (`cxx:79-82`): the OCCT body
    /// is `*this = GT`, i.e. every field, tolerances and `NbPntMax` included.
    pub fn define_copy(&mut self, gt: &GeomTool) {
        *self = *gt;
    }

    /// `GetTolerances(tol3d, tol2d)` (`cxx:85-91`).
    pub fn get_tolerances(&self) -> (f64, f64) {
        (self.tol3d, self.tol2d)
    }

    /// `SetTolerances(tol3d, tol2d)` (`cxx:94-98`).
    pub fn set_tolerances(&mut self, tol3d: f64, tol2d: f64) {
        self.tol3d = tol3d;
        self.tol2d = tol2d;
    }

    /// `NbPntMax()` (`cxx:102-105`).
    pub fn nb_pnt_max(&self) -> i32 {
        self.nb_pnt_max
    }

    /// `SetNbPntMax(NbPntMax)` (`cxx:109-112`).
    pub fn set_nb_pnt_max(&mut self, nb_pnt_max: i32) {
        self.nb_pnt_max = nb_pnt_max;
    }

    /// `CompC3D()` (`cxx:117-120`).
    pub fn comp_c3d(&self) -> bool {
        self.comp_c3d
    }

    /// `TypeC3D()` (`cxx:124-127`).
    pub fn type_c3d(&self) -> OutCurveType {
        self.type_c3d
    }

    /// `CompPC1()` (`cxx:131-134`).
    pub fn comp_pc1(&self) -> bool {
        self.comp_pc1
    }

    /// `CompPC2()` (`cxx:138-141`).
    pub fn comp_pc2(&self) -> bool {
        self.comp_pc2
    }
}
