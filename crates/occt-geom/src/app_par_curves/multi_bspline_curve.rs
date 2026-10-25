//! `AppParCurves_MultiBSpCurve` — the B-spline flavour of
//! [`MultiCurve`](crate::app_par_curves::MultiCurve): the same parallel pole
//! arrays, evaluated as one non-periodic B-spline of shared knots / multiplicities.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_MultiBSpCurve.cxx:28-306`
//! and `AppParCurves_MultiBSpCurve.hxx:35-120` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `myknots` / `mymults` (`hxx`) are null until `SetKnots` /
//!   `SetMultiplicities` run; `Option` is that state and the accessors that
//!   would dereference a null handle in OCCT (`cxx:105`, `cxx:110`) panic.
//! - `ComputeDegree` (`cxx:28-36`) is `sum(Mults) - NbPoles - 1`; it runs in
//!   the knot/mult constructors (`cxx:67`, `cxx:80`) and in
//!   `SetMultiplicities` (`cxx:97`), but **not** in `SetKnots` (`cxx:87-90`).
//! - `BSplCLib::D0/D1/D2` with `Index = 0`, `Periodic = false`, distinct
//!   `Knots` and `&Mults` (`cxx:135-145`, `cxx:155-167`, `cxx:181-192`,
//!   `cxx:206-219`, `cxx:236-247`, `cxx:264-276`) is
//!   `occt_core::bspl::curve_eval::{curve_d0, curve_d1, curve_d2}` with
//!   `mults = Some(..)` and `index = 0`.
//! - The 2d arms lift each `gp_Pnt2d` to `(x, y, 0)` and read `x` / `y` back,
//!   coordinate-wise identical to the `gp_Pnt2d` de Boor recursion.

use occt_core::bspl::curve_eval::{curve_d0, curve_d1, curve_d2};
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::multi_point::MultiPoint;

/// `ComputeDegree(Mults, NbPoles)` (`cxx:28-36`).
fn compute_degree(mults: &[i32], nb_poles: usize) -> i32 {
    let sum: i32 = mults.iter().sum();
    sum - nb_poles as i32 - 1
}

/// `gp_Pnt2d` -> `gp_Pnt` lift used by the 2d evaluation arms.
fn lift2d(p: &GpPnt2d) -> GpPnt {
    GpPnt::new(p.x(), p.y(), 0.0)
}

/// `AppParCurves_MultiBSpCurve` (`AppParCurves_MultiBSpCurve.hxx:35-120`).
#[derive(Debug, Clone)]
pub struct MultiBSpCurve {
    /// The `AppParCurves_MultiCurve` base (`cxx:49`, `cxx:59`, `cxx:77`).
    base: MultiCurve,
    /// `myknots` (`hxx`): the distinct knots.
    knots: Option<Vec<f64>>,
    /// `mymults` (`hxx`): the multiplicities of `knots`.
    mults: Option<Vec<i32>>,
    /// `myDegree` (`hxx`).
    degree: i32,
}

impl Default for MultiBSpCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiBSpCurve {
    /// `AppParCurves_MultiBSpCurve()` (`cxx:42-45`).
    pub fn new() -> Self {
        Self {
            base: MultiCurve::new(),
            knots: None,
            mults: None,
            degree: 0,
        }
    }

    /// `AppParCurves_MultiBSpCurve(NbPol)` (`cxx:47-51`).
    pub fn with_nb_poles(nb_pol: usize) -> Self {
        Self {
            base: MultiCurve::with_nb_poles(nb_pol),
            knots: None,
            mults: None,
            degree: 0,
        }
    }

    /// `AppParCurves_MultiBSpCurve(tabMU, Knots, Mults)` (`cxx:53-68`).
    pub fn from_tab(tab: &[MultiPoint], knots: &[f64], mults: &[i32]) -> Self {
        let base = MultiCurve::from_tab(tab);
        let degree = compute_degree(mults, base.nb_poles());
        Self {
            base,
            knots: Some(knots.to_vec()),
            mults: Some(mults.to_vec()),
            degree,
        }
    }

    /// `AppParCurves_MultiBSpCurve(SC, Knots, Mults)` (`cxx:70-81`).
    pub fn from_curve(sc: &MultiCurve, knots: &[f64], mults: &[i32]) -> Self {
        let degree = compute_degree(mults, sc.nb_poles());
        Self {
            base: sc.clone(),
            knots: Some(knots.to_vec()),
            mults: Some(mults.to_vec()),
            degree,
        }
    }

    /// The `AppParCurves_MultiCurve` base.
    pub fn base(&self) -> &MultiCurve {
        &self.base
    }

    /// The mutable base; `SetValue` is inherited from `MultiCurve`.
    pub fn base_mut(&mut self) -> &mut MultiCurve {
        &mut self.base
    }

    /// `SetKnots` (`cxx:85-90`): the degree is **not** recomputed.
    pub fn set_knots(&mut self, knots: &[f64]) {
        self.knots = Some(knots.to_vec());
    }

    /// `SetMultiplicities` (`cxx:94-98`): the degree **is** recomputed.
    pub fn set_multiplicities(&mut self, mults: &[i32]) {
        self.mults = Some(mults.to_vec());
        self.degree = compute_degree(mults, self.base.nb_poles());
    }

    /// `Knots` (`cxx:102-105`).
    pub fn knots(&self) -> &[f64] {
        match &self.knots {
            Some(k) => k,
            None => panic!("Standard_NullObject: AppParCurves_MultiBSpCurve::Knots"),
        }
    }

    /// `Multiplicities` (`cxx:108-111`).
    pub fn multiplicities(&self) -> &[i32] {
        match &self.mults {
            Some(m) => m,
            None => panic!("Standard_NullObject: AppParCurves_MultiBSpCurve::Multiplicities"),
        }
    }

    /// `Degree` (`cxx:114-117`).
    pub fn degree(&self) -> i32 {
        self.degree
    }

    /// `Value(CuIndex, U, Pt)` (`cxx:123-146`).
    pub fn value3d(&self, cu_index: i32, u: f64) -> GpPnt {
        if self.base.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::Value");
        }
        let poles = self.base.curve3d(cu_index);
        curve_d0(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        )
    }

    /// `Value(CuIndex, U, Pt2d)` (`cxx:148-167`).
    pub fn value2d(&self, cu_index: i32, u: f64) -> GpPnt2d {
        if self.base.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::Value");
        }
        let poles: Vec<GpPnt> = self.base.curve2d(cu_index).iter().map(lift2d).collect();
        let p = curve_d0(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        );
        GpPnt2d::new(p.x(), p.y())
    }

    /// `D1(CuIndex, U, Pt, V1)` (`cxx:171-192`).
    pub fn d1_3d(&self, cu_index: i32, u: f64) -> (GpPnt, GpVec) {
        if self.base.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::D1");
        }
        let poles = self.base.curve3d(cu_index);
        curve_d1(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        )
    }

    /// `D2(CuIndex, U, Pt, V1, V2)` (`cxx:196-220`).
    pub fn d2_3d(&self, cu_index: i32, u: f64) -> (GpPnt, GpVec, GpVec) {
        if self.base.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::D2");
        }
        let poles = self.base.curve3d(cu_index);
        curve_d2(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        )
    }

    /// `D1(CuIndex, U, Pt2d, V1)` (`cxx:223-247`).
    pub fn d1_2d(&self, cu_index: i32, u: f64) -> (GpPnt2d, GpVec2d) {
        if self.base.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::D1");
        }
        let poles: Vec<GpPnt> = self.base.curve2d(cu_index).iter().map(lift2d).collect();
        let (p, v) = curve_d1(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        );
        (GpPnt2d::new(p.x(), p.y()), GpVec2d::new(v.x(), v.y()))
    }

    /// `D2(CuIndex, U, Pt2d, V1, V2)` (`cxx:251-276`).
    pub fn d2_2d(&self, cu_index: i32, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        if self.base.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::D2");
        }
        let poles: Vec<GpPnt> = self.base.curve2d(cu_index).iter().map(lift2d).collect();
        let (p, v1, v2) = curve_d2(
            u,
            0,
            self.degree,
            false,
            &poles,
            None,
            self.knots(),
            self.mults.as_deref(),
        );
        (
            GpPnt2d::new(p.x(), p.y()),
            GpVec2d::new(v1.x(), v1.y()),
            GpVec2d::new(v2.x(), v2.y()),
        )
    }

    /// `Dump(o)` (`cxx:280-306`): the emitted lines only; the pole listing is
    /// commented out in OCCT.
    pub fn dump(&self) -> String {
        format!(
            "AppParCurves_MultiBSpCurve dump:\n It contains {} BSpline curves \n The poles are: \n",
            self.base.nb_curves()
        )
    }
}
