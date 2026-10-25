//! `AppParCurves_MultiCurve` — a set of parallel Bezier curves sharing one
//! pole count, addressed by curve index and evaluated through `BSplCLib`.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_MultiCurve.cxx:27-327`
//! and `AppParCurves_MultiCurve.hxx:44-190` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<NCollection_HArray1<AppParCurves_MultiPoint>> tabPoint`
//!   (`hxx:181`) may be null: the default constructor leaves it so and
//!   `NbCurves` / `NbPoles` report 0 (`cxx:55-71`). The port keeps
//!   `Option<Vec<MultiPoint>>` with `None` as that state; the methods that
//!   dereference the handle unconditionally in OCCT (`Degree`, `SetValue`,
//!   `Dimension`, ...) panic here, OCCT's null dereference counterpart.
//! - `BSplCLib::D0 / D1 / D2` on a pole array without knots (`cxx:198`,
//!   `cxx:231`, `cxx:256`, `cxx:276`, `cxx:297`) is the Bezier evaluation
//!   `BSplCLib.cxx:843-998` with implicit flat knots on `[0, 1]` and
//!   `Degree = Poles.Length() - 1`; [`flat_knots`] builds that knot vector and
//!   [`curve_d0`] / [`curve_d1`] / [`curve_d2`] are the faithful evaluators
//!   (`BSplCLib.cxx` ported in `occt_core::bspl::curve_eval`).
//! - The 2d overloads (`cxx:201-215`, `cxx:258-299`) run the same de Boor
//!   recursion on two coordinates (`BSplCLib_2.cxx` `gp_Pnt2d` arm); the port
//!   lifts each `gp_Pnt2d` to `(x, y, 0)` and reads the first two coordinates
//!   back, which is coordinate-wise identical.
//! - `Pole` / `Pole2d` reject only when **both** indices are `<= 0`
//!   (`cxx:121`, `cxx:131`: `&&`, not `||`); kept verbatim.

use occt_core::bspl::curve_eval::{curve_d0, curve_d1, curve_d2};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec};

use crate::app_par_curves::multi_point::MultiPoint;

/// The implicit flat knot vector of a degree `poles.len() - 1` Bezier curve on
/// `[0, 1]` (`BSplCLib.cxx` `D0` with `Knots == nullptr`, `Mults == nullptr`).
fn flat_knots(degree: usize) -> Vec<f64> {
    let mut knots = vec![0.0; degree + 1];
    knots.extend(std::iter::repeat(1.0).take(degree + 1));
    knots
}

/// `gp_Pnt2d` -> `gp_Pnt` lift used by the 2d evaluation arms.
fn lift2d(p: &GpPnt2d) -> GpPnt {
    GpPnt::new(p.x(), p.y(), 0.0)
}

/// `AppParCurves_MultiCurve` (`AppParCurves_MultiCurve.hxx:44-190`).
#[derive(Debug, Clone)]
pub struct MultiCurve {
    /// `tabPoint` (`hxx:181`): `None` is the null handle.
    tab_point: Option<Vec<MultiPoint>>,
}

impl Default for MultiCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiCurve {
    /// `AppParCurves_MultiCurve()` (`cxx:27`): the handle stays null.
    pub fn new() -> Self {
        Self { tab_point: None }
    }

    /// `AppParCurves_MultiCurve(NbPol)` (`cxx:29-32`): `NbPol` default
    /// `MultiPoint`s.
    pub fn with_nb_poles(nb_pol: usize) -> Self {
        Self {
            tab_point: Some(vec![MultiPoint::new(); nb_pol]),
        }
    }

    /// `AppParCurves_MultiCurve(tabMU)` (`cxx:34-41`).
    pub fn from_tab(tab: &[MultiPoint]) -> Self {
        Self {
            tab_point: Some(tab.to_vec()),
        }
    }

    /// The pole array; `None` reproduces the null OCCT handle.
    pub fn tab_point(&self) -> Option<&[MultiPoint]> {
        self.tab_point.as_deref()
    }

    /// The 1-based `tabPoint->Value(Index)` (`cxx:140`).
    fn value_at(&self, index: i32) -> &MultiPoint {
        let Some(tab) = &self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::Value");
        };
        if index <= 0 || index > tab.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Value");
        }
        &tab[(index - 1) as usize]
    }

    /// `Dimension(Index)` (`cxx:45-53`): the counts of the **first** pole.
    pub fn dimension(&self, index: i32) -> usize {
        let first = self.value_at(1);
        let nb = (first.nb_points() + first.nb_points2d()) as i32;
        if index <= 0 || index > nb {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Dimension");
        }
        first.dimension(index)
    }

    /// `NbCurves` (`cxx:55-63`): `nbP + nbP2d` of the first pole, 0 when the
    /// array is null.
    pub fn nb_curves(&self) -> usize {
        match &self.tab_point {
            None => 0,
            Some(tab) => match tab.first() {
                None => 0,
                Some(mp) => mp.nb_points() + mp.nb_points2d(),
            },
        }
    }

    /// `NbPoles` (`cxx:65-71`).
    pub fn nb_poles(&self) -> usize {
        match &self.tab_point {
            None => 0,
            Some(tab) => tab.len(),
        }
    }

    /// `Degree` (`cxx:73-76`): `tabPoint->Length() - 1`.
    pub fn degree(&self) -> usize {
        match &self.tab_point {
            None => panic!("Standard_NullObject: AppParCurves_MultiCurve::Degree"),
            Some(tab) => tab.len().saturating_sub(1),
        }
    }

    /// `SetNbPoles` (`cxx:80-83`): a fresh array of default `MultiPoint`s.
    pub fn set_nb_poles(&mut self, nb_poles: usize) {
        self.tab_point = Some(vec![MultiPoint::new(); nb_poles]);
    }

    /// `SetValue(Index, MPoint)` (`cxx:85-94`).
    pub fn set_value(&mut self, index: i32, mpoint: &MultiPoint) {
        let Some(tab) = &mut self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::SetValue");
        };
        if index <= 0 || index > tab.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::SetValue");
        }
        tab[(index - 1) as usize] = mpoint.clone();
    }

    /// `Curve(CuIndex, TabPnt)` (`cxx:96-105`): the 3d poles of one curve.
    pub fn curve3d(&self, cu_index: i32) -> Vec<GpPnt> {
        if cu_index <= 0 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Curve");
        }
        let Some(tab) = &self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::Curve");
        };
        (1..=tab.len() as i32)
            .map(|i| *tab[(i - 1) as usize].point(cu_index))
            .collect()
    }

    /// `Curve(CuIndex, TabPnt2d)` (`cxx:107-117`).
    pub fn curve2d(&self, cu_index: i32) -> Vec<GpPnt2d> {
        if cu_index <= 0 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Curve");
        }
        let Some(tab) = &self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::Curve");
        };
        (1..=tab.len() as i32)
            .map(|i| *tab[(i - 1) as usize].point2d(cu_index))
            .collect()
    }

    /// `Pole(CuIndex, Nieme)` (`cxx:119-127`): note the `&&` of `cxx:121`.
    pub fn pole(&self, cu_index: i32, nieme: i32) -> &GpPnt {
        if cu_index <= 0 && nieme <= 0 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Pole");
        }
        self.value_at(nieme).point(cu_index)
    }

    /// `Pole2d(CuIndex, Nieme)` (`cxx:129-137`).
    pub fn pole2d(&self, cu_index: i32, nieme: i32) -> &GpPnt2d {
        if cu_index <= 0 && nieme <= 0 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Pole2d");
        }
        self.value_at(nieme).point2d(cu_index)
    }

    /// `Value(Index)` (`cxx:139-145`).
    pub fn value(&self, index: i32) -> &MultiPoint {
        self.value_at(index)
    }

    /// `Transform(CuIndex, ...)` (`cxx:147-164`).
    pub fn transform(&mut self, cu_index: i32, x: f64, dx: f64, y: f64, dy: f64, z: f64, dz: f64) {
        if self.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Transform");
        }
        let Some(tab) = &mut self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::Transform");
        };
        for mp in tab.iter_mut() {
            mp.transform(cu_index, x, dx, y, dy, z, dz);
        }
    }

    /// `Transform2d(CuIndex, ...)` (`cxx:166-181`).
    pub fn transform2d(&mut self, cu_index: i32, x: f64, dx: f64, y: f64, dy: f64) {
        if self.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Transform2d");
        }
        let Some(tab) = &mut self.tab_point else {
            panic!("Standard_NullObject: AppParCurves_MultiCurve::Transform2d");
        };
        for mp in tab.iter_mut() {
            mp.transform2d(cu_index, x, dx, y, dy);
        }
    }

    /// `Value(CuIndex, U, Pt)` (`cxx:183-199`): the 3d evaluator.
    pub fn value3d(&self, cu_index: i32, u: f64) -> GpPnt {
        if self.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Value");
        }
        let poles = self.curve3d(cu_index);
        let degree = poles.len().saturating_sub(1);
        curve_d0(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None)
    }

    /// `Value(CuIndex, U, Pt2d)` (`cxx:201-215`): the 2d evaluator.
    pub fn value2d(&self, cu_index: i32, u: f64) -> GpPnt2d {
        if self.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::Value");
        }
        let poles: Vec<GpPnt> = self.curve2d(cu_index).iter().map(lift2d).collect();
        let degree = poles.len().saturating_sub(1);
        let p = curve_d0(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None);
        GpPnt2d::new(p.x(), p.y())
    }

    /// `D1(CuIndex, U, Pt, V1)` (`cxx:217-232`).
    pub fn d1_3d(&self, cu_index: i32, u: f64) -> (GpPnt, GpVec) {
        if self.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::D1");
        }
        let poles = self.curve3d(cu_index);
        let degree = poles.len().saturating_sub(1);
        curve_d1(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None)
    }

    /// `D2(CuIndex, U, Pt, V1, V2)` (`cxx:234-256`).
    pub fn d2_3d(&self, cu_index: i32, u: f64) -> (GpPnt, GpVec, GpVec) {
        if self.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::D2");
        }
        let poles = self.curve3d(cu_index);
        let degree = poles.len().saturating_sub(1);
        curve_d2(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None)
    }

    /// `D1(CuIndex, U, Pt2d, V1)` (`cxx:258-277`).
    pub fn d1_2d(&self, cu_index: i32, u: f64) -> (GpPnt2d, occt_core::gp::GpVec2d) {
        if self.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::D1");
        }
        let poles: Vec<GpPnt> = self.curve2d(cu_index).iter().map(lift2d).collect();
        let degree = poles.len().saturating_sub(1);
        let (p, v) = curve_d1(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None);
        (GpPnt2d::new(p.x(), p.y()), occt_core::gp::GpVec2d::new(v.x(), v.y()))
    }

    /// `D2(CuIndex, U, Pt2d, V1, V2)` (`cxx:279-299`).
    pub fn d2_2d(
        &self,
        cu_index: i32,
        u: f64,
    ) -> (GpPnt2d, occt_core::gp::GpVec2d, occt_core::gp::GpVec2d) {
        if self.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiCurve::D2");
        }
        let poles: Vec<GpPnt> = self.curve2d(cu_index).iter().map(lift2d).collect();
        let degree = poles.len().saturating_sub(1);
        let (p, v1, v2) =
            curve_d2(u, 1, degree as i32, false, &poles, None, &flat_knots(degree), None);
        (
            GpPnt2d::new(p.x(), p.y()),
            occt_core::gp::GpVec2d::new(v1.x(), v1.y()),
            occt_core::gp::GpVec2d::new(v2.x(), v2.y()),
        )
    }

    /// `Dump(o)` (`cxx:301-327`): the two header lines OCCT actually emits; the
    /// pole listing is commented out in `cxx:307-325` and stays out.
    pub fn dump(&self) -> String {
        format!(
            "AppParCurves_MultiCurve dump:\n It contains {} Bezier curves of degree {}\n The poles are: \n",
            self.nb_curves(),
            self.nb_poles().saturating_sub(1)
        )
    }
}
