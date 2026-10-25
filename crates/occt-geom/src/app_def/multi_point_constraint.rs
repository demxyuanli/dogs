//! `AppDef_MultiPointConstraint` — one constraint column of an
//! `AppDef_MultiLine`, i.e. the 3d / 2d points of a single parameter plus
//! their optional tangency and curvature vectors.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_MultiPointConstraint.hxx:38-196`
//! and `AppDef_MultiPointConstraint.cxx:26-318` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The class derives from `AppParCurves_MultiPoint`; the port holds that
//!   base object and re-exports its pole accessors ([`MultiPointConstraint::point`],
//!   [`MultiPointConstraint::point2d`], ...). `nbP` / `nbP2d` used by the
//!   array allocations below are the base class counts.
//! - The four OCCT handles `tabTang` / `tabCurv` / `tabTang2d` / `tabCurv2d`
//!   are allocated by some constructors and left null by others, and
//!   [`MultiPointConstraint::is_tangency_point`] tests exactly for
//!   "not both null" (`cxx:306-309`). `Option<Vec<..>>` keeps that
//!   distinction; `Some(empty)` is an allocated-but-empty array, `None` is
//!   the null handle.
//! - `Standard_ConstructionError` and `Standard_OutOfRange` become panics,
//!   the same convention as `app_par_curves::multi_point`.
//! - `Tang` / `Curv` / `Tang2d` / `Curv2d` dereference the handle without a
//!   null test in OCCT (`cxx:228-301`); the port panics with
//!   `Standard_NullObject` instead of reading a null handle.
//! - `cxx:84` loops `i <= tabVec.Length()` while filling `tabCurv`; the
//!   preceding check (`cxx:59-63`) already forces `tabCur.Length() ==
//!   tabP.Length() == tabVec.Length()`, so the loop bound is kept verbatim.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::app_par_curves::multi_point::MultiPoint;

/// `AppDef_MultiPointConstraint` (`AppDef_MultiPointConstraint.hxx:49-196`).
#[derive(Debug, Clone)]
pub struct MultiPointConstraint {
    /// The `AppParCurves_MultiPoint` base (`hxx:49`).
    point: MultiPoint,
    /// `tabTang` (`hxx:193`).
    tang: Option<Vec<GpVec>>,
    /// `tabCurv` (`hxx:194`).
    curv: Option<Vec<GpVec>>,
    /// `tabTang2d` (`hxx:195`).
    tang2d: Option<Vec<GpVec2d>>,
    /// `tabCurv2d` (`hxx:196`).
    curv2d: Option<Vec<GpVec2d>>,
}

impl Default for MultiPointConstraint {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiPointConstraint {
    /// `AppDef_MultiPointConstraint()` (`cxx:27`): every array null.
    pub fn new() -> Self {
        Self {
            point: MultiPoint::new(),
            tang: None,
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(NbPoints, NbPoints2d)` (`cxx:29-32`):
    /// the points are default-constructed and filled later.
    pub fn with_size(nb_points: usize, nb_points2d: usize) -> Self {
        Self {
            point: MultiPoint::with_size(nb_points, nb_points2d),
            tang: None,
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP)` (`cxx:34-37`): 3d poles only.
    pub fn from_3d(tab_p: &[GpPnt]) -> Self {
        Self {
            point: MultiPoint::from_3d(tab_p),
            tang: None,
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP)` over `gp_Pnt2d` (`cxx:39-42`):
    /// 2d poles only.
    pub fn from_2d(tab_p2d: &[GpPnt2d]) -> Self {
        Self {
            point: MultiPoint::from_2d(tab_p2d),
            tang: None,
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP, tabP2d)` (`cxx:44-47`).
    pub fn from_3d_2d(tab_p: &[GpPnt], tab_p2d: &[GpPnt2d]) -> Self {
        Self {
            point: MultiPoint::from_3d_2d(tab_p, tab_p2d),
            tang: None,
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP, tabP2d, tabVec, tabVec2d,
    /// tabCur, tabCur2d)` (`cxx:49-92`): tangency and curvature.
    ///
    /// The `Standard_ConstructionError` of `cxx:59-63` is a panic.
    pub fn from_3d_2d_tang_curv(
        tab_p: &[GpPnt],
        tab_p2d: &[GpPnt2d],
        tab_vec: &[GpVec],
        tab_vec2d: &[GpVec2d],
        tab_cur: &[GpVec],
        tab_cur2d: &[GpVec2d],
    ) -> Self {
        if tab_p.len() != tab_vec.len()
            || tab_p2d.len() != tab_vec2d.len()
            || tab_cur.len() != tab_p.len()
            || tab_cur2d.len() != tab_p2d.len()
        {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_3d_2d(tab_p, tab_p2d),
            tang: Some(tab_vec.to_vec()),
            tang2d: Some(tab_vec2d.to_vec()),
            curv: Some(tab_cur.to_vec()),
            curv2d: Some(tab_cur2d.to_vec()),
        }
    }

    /// `AppDef_MultiPointConstraint(tabP, tabP2d, tabVec, tabVec2d)`
    /// (`cxx:94-121`): tangency only.
    pub fn from_3d_2d_tang(
        tab_p: &[GpPnt],
        tab_p2d: &[GpPnt2d],
        tab_vec: &[GpVec],
        tab_vec2d: &[GpVec2d],
    ) -> Self {
        if tab_p.len() != tab_vec.len() || tab_p2d.len() != tab_vec2d.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_3d_2d(tab_p, tab_p2d),
            tang: Some(tab_vec.to_vec()),
            curv: None,
            tang2d: Some(tab_vec2d.to_vec()),
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP, tabVec)` (`cxx:122-139`): 3d
    /// tangency only; `tabTang2d` stays null.
    pub fn from_3d_tang(tab_p: &[GpPnt], tab_vec: &[GpVec]) -> Self {
        if tab_p.len() != tab_vec.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_3d(tab_p),
            tang: Some(tab_vec.to_vec()),
            curv: None,
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP, tabVec, tabCur)` (`cxx:141-166`):
    /// 3d tangency and curvature.
    pub fn from_3d_tang_curv(tab_p: &[GpPnt], tab_vec: &[GpVec], tab_cur: &[GpVec]) -> Self {
        if tab_p.len() != tab_vec.len() || tab_p.len() != tab_cur.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_3d(tab_p),
            tang: Some(tab_vec.to_vec()),
            curv: Some(tab_cur.to_vec()),
            tang2d: None,
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP2d, tabVec2d)` (`cxx:168-187`): 2d
    /// tangency only; `tabTang` stays null.
    pub fn from_2d_tang(tab_p2d: &[GpPnt2d], tab_vec2d: &[GpVec2d]) -> Self {
        if tab_p2d.len() != tab_vec2d.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_2d(tab_p2d),
            tang: None,
            curv: None,
            tang2d: Some(tab_vec2d.to_vec()),
            curv2d: None,
        }
    }

    /// `AppDef_MultiPointConstraint(tabP2d, tabVec2d, tabCur2d)`
    /// (`cxx:189-212`): 2d tangency and curvature.
    pub fn from_2d_tang_curv(
        tab_p2d: &[GpPnt2d],
        tab_vec2d: &[GpVec2d],
        tab_cur2d: &[GpVec2d],
    ) -> Self {
        if tab_p2d.len() != tab_vec2d.len() || tab_cur2d.len() != tab_p2d.len() {
            panic!("Standard_ConstructionError: AppDef_MultiPointConstraint");
        }
        Self {
            point: MultiPoint::from_2d(tab_p2d),
            tang: None,
            curv: None,
            tang2d: Some(tab_vec2d.to_vec()),
            curv2d: Some(tab_cur2d.to_vec()),
        }
    }

    /// The `AppParCurves_MultiPoint` base object.
    pub fn as_multi_point(&self) -> &MultiPoint {
        &self.point
    }

    /// `NbPoints` (`AppParCurves_MultiPoint.lxx:25-30`), inherited.
    pub fn nb_points(&self) -> usize {
        self.point.nb_points()
    }

    /// `NbPoints2d` (`AppParCurves_MultiPoint.lxx:32-35`), inherited.
    pub fn nb_points2d(&self) -> usize {
        self.point.nb_points2d()
    }

    /// `SetPoint(Index, Point)` (`AppParCurves_MultiPoint.cxx:140-146`).
    pub fn set_point(&mut self, index: i32, point: &GpPnt) {
        self.point.set_point(index, point);
    }

    /// `Point(Index)` (`AppParCurves_MultiPoint.cxx:148-153`).
    pub fn point(&self, index: i32) -> &GpPnt {
        self.point.point(index)
    }

    /// `SetPoint2d(Index, Point)` (`AppParCurves_MultiPoint.cxx:155-159`).
    pub fn set_point2d(&mut self, index: i32, point: &GpPnt2d) {
        self.point.set_point2d(index, point);
    }

    /// `Point2d(Index)` (`AppParCurves_MultiPoint.cxx:161-166`).
    pub fn point2d(&self, index: i32) -> &GpPnt2d {
        self.point.point2d(index)
    }

    /// `SetTang(Index, Tang)` (`cxx:214-226`).
    ///
    /// The array is allocated before the range test, as in OCCT.
    pub fn set_tang(&mut self, index: i32, tang: &GpVec) {
        if self.tang.is_none() {
            self.tang = Some(vec![GpVec::zero(); self.point.nb_points()]);
        }
        if index <= 0 || index > self.point.nb_points() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::SetTang");
        }
        self.tang.as_mut().unwrap()[(index - 1) as usize] = *tang;
    }

    /// `Tang(Index)` (`cxx:228-236`).
    pub fn tang(&self, index: i32) -> GpVec {
        if index <= 0 || index > self.point.nb_points() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::Tang");
        }
        match &self.tang {
            Some(t) => t[(index - 1) as usize],
            None => panic!("Standard_NullObject: AppDef_MultiPointConstraint::Tang"),
        }
    }

    /// `SetTang2d(Index, Tang2d)` (`cxx:238-249`): `Index` counts the 3d
    /// points first.
    pub fn set_tang2d(&mut self, index: i32, tang2d: &GpVec2d) {
        if self.tang2d.is_none() {
            self.tang2d = Some(vec![GpVec2d::zero(); self.point.nb_points2d()]);
        }
        let nb_p = self.point.nb_points() as i32;
        if index <= nb_p || index > nb_p + self.point.nb_points2d() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::SetTang2d");
        }
        self.tang2d.as_mut().unwrap()[(index - nb_p - 1) as usize] = *tang2d;
    }

    /// `Tang2d(Index)` (`cxx:251-258`).
    pub fn tang2d(&self, index: i32) -> GpVec2d {
        let nb_p = self.point.nb_points() as i32;
        if index <= nb_p || index > nb_p + self.point.nb_points2d() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::Tang2d");
        }
        match &self.tang2d {
            Some(t) => t[(index - nb_p - 1) as usize],
            None => panic!("Standard_NullObject: AppDef_MultiPointConstraint::Tang2d"),
        }
    }

    /// `SetCurv(Index, Curv)` (`cxx:260-271`).
    pub fn set_curv(&mut self, index: i32, curv: &GpVec) {
        if self.curv.is_none() {
            self.curv = Some(vec![GpVec::zero(); self.point.nb_points()]);
        }
        if index <= 0 || index > self.point.nb_points() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::SetCurv");
        }
        self.curv.as_mut().unwrap()[(index - 1) as usize] = *curv;
    }

    /// `Curv(Index)` (`cxx:273-280`).
    pub fn curv(&self, index: i32) -> GpVec {
        if index <= 0 || index > self.point.nb_points() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::Curv");
        }
        match &self.curv {
            Some(c) => c[(index - 1) as usize],
            None => panic!("Standard_NullObject: AppDef_MultiPointConstraint::Curv"),
        }
    }

    /// `SetCurv2d(Index, Curv2d)` (`cxx:282-293`).
    pub fn set_curv2d(&mut self, index: i32, curv2d: &GpVec2d) {
        if self.curv2d.is_none() {
            self.curv2d = Some(vec![GpVec2d::zero(); self.point.nb_points2d()]);
        }
        let nb_p = self.point.nb_points() as i32;
        if index <= nb_p || index > nb_p + self.point.nb_points2d() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::SetCurv2d");
        }
        self.curv2d.as_mut().unwrap()[(index - nb_p - 1) as usize] = *curv2d;
    }

    /// `Curv2d(Index)` (`cxx:295-302`).
    pub fn curv2d(&self, index: i32) -> GpVec2d {
        let nb_p = self.point.nb_points() as i32;
        if index <= nb_p || index > nb_p + self.point.nb_points2d() as i32 {
            panic!("Standard_OutOfRange: AppDef_MultiPointConstraint::Curv2d");
        }
        match &self.curv2d {
            Some(c) => c[(index - nb_p - 1) as usize],
            None => panic!("Standard_NullObject: AppDef_MultiPointConstraint::Curv2d"),
        }
    }

    /// `IsTangencyPoint` (`cxx:304-308`): not both tangency handles null.
    pub fn is_tangency_point(&self) -> bool {
        !(self.tang.is_none() && self.tang2d.is_none())
    }

    /// `IsCurvaturePoint` (`cxx:310-313`): not both curvature handles null.
    pub fn is_curvature_point(&self) -> bool {
        !(self.curv.is_none() && self.curv2d.is_none())
    }

    /// `Dump(o)` (`cxx:315-318`): only the header line is printed.
    pub fn dump(&self) -> String {
        "AppDef_MultiPointConstraint dump:\n".to_string()
    }
}
