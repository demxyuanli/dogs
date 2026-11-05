//! `GeomPlate_PointConstraint` — a point as a surface-deformation constraint.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_PointConstraint.hxx`,
//! `GeomPlate_PointConstraint.cxx:29-215` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `GeomLProp_SLProps myLProp` is [`SlProps`]; the port's `SlProps` carries no
//!   surface (it takes one per call, see `geom_lprop/sl_props.rs`), so
//!   [`PointConstraint::l_prop_surf`] hands the surface back with the props —
//!   the only consumer is `GeomPlate_BuildPlateSurface.cxx:2775`.
//! - `occ::handle<Geom_Surface> mySurf` is `Option<Arc<dyn Surface>>`: the
//!   point-only constructor (`cxx:31-48`) leaves it null.
//! - The two `Standard_Failure` raises (`cxx:44-47`) become the `Err` of
//!   [`PointConstraint::from_point`].

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};

use crate::geom_lprop::sl_props::SlProps;
use crate::surface::Surface;

/// `class GeomPlate_PointConstraint` (`hxx:32-157`): "defines points as
/// constraints to be used to deform a surface".
///
/// `Clone` stands in for the reference's `occ::handle` sharing: OCCT hands the
/// same constraint object to the builder and keeps its own reference.
#[derive(Clone)]
pub struct PointConstraint {
    /// `myOrder` (`hxx:140`).
    order: i32,
    /// `myLProp` (`hxx:141`).
    l_prop: SlProps,
    /// `myPoint` (`hxx:142`).
    point: GpPnt,
    /// `myPt2d` (`hxx:143`).
    pt2d: GpPnt2d,
    /// `mySurf` (`hxx:144`).
    surf: Option<Arc<dyn Surface>>,
    /// `myD11` (`hxx:145`): the surface `D1U`.
    d11: GpVec,
    /// `myD12` (`hxx:146`): the surface `D1V`.
    d12: GpVec,
    /// `myD21` (`hxx:147`): the surface `D2U`.
    d21: GpVec,
    /// `myD22` (`hxx:148`): the surface `D2V`.
    d22: GpVec,
    /// `myD23` (`hxx:149`): the surface `D2UV`.
    d23: GpVec,
    /// `myU` (`hxx:150`).
    u: f64,
    /// `myV` (`hxx:151`).
    v: f64,
    /// `myTolDist` (`hxx:152`).
    tol_dist: f64,
    /// `myTolAng` (`hxx:153`).
    tol_ang: f64,
    /// `myTolCurv` (`hxx:154`).
    tol_curv: f64,
    /// `hasPnt2dOnSurf` (`hxx:155`).
    has_pnt2d_on_surf: bool,
}

impl PointConstraint {
    /// `GeomPlate_PointConstraint(Pt, Order, TolDist)` (`cxx:31-48`): a
    /// constraint defined by a 3D point alone.
    ///
    /// Raises when `Order > 1 || Order < -1` (`cxx:44-48`).
    pub fn from_point(pt: GpPnt, order: i32, tol_dist: f64) -> Result<Self, String> {
        // `myLProp(2, TolDist)` is the `(N, Resolution)` base constructor
        // (`GeomLProp_SLProps.hxx:93-99`), whose `myCN` starts at 0.
        let l_prop = SlProps::empty(2, tol_dist)?;
        if order > 1 || order < -1 {
            return Err(
                "GeomPlate_PointConstraint : the constraint must 0 or -1 with a point".to_string(),
            );
        }
        Ok(Self {
            order,
            l_prop,
            point: pt,
            pt2d: GpPnt2d::zero(),
            surf: None,
            d11: GpVec::zero(),
            d12: GpVec::zero(),
            d21: GpVec::zero(),
            d22: GpVec::zero(),
            d23: GpVec::zero(),
            u: 0.0,
            v: 0.0,
            tol_dist,
            tol_ang: 0.0,
            tol_curv: 0.0,
            has_pnt2d_on_surf: false,
        })
    }

    /// `GeomPlate_PointConstraint(U, V, Surf, Order, TolDist, TolAng, TolCurv)`
    /// (`cxx:52-73`): a constraint defined by the intersection point of `U` and
    /// `V` on `Surf`.
    pub fn from_surface(
        u: f64,
        v: f64,
        surf: Arc<dyn Surface>,
        order: i32,
        tol_dist: f64,
        tol_ang: f64,
        tol_curv: f64,
    ) -> Result<Self, String> {
        // `Surf->D2(myU, myV, myPoint, myD11, myD12, myD21, myD22, myD23)`
        // (`cxx:70`).
        let (point, d11, d12, d21, d22, d23) = surf.d2(u, v);
        // `myLProp(2, TolDist)` is the `(N, Resolution)` base constructor,
        // then `myLProp.SetSurface(Surf)` (`cxx:71`) raises `myCN` to 4.
        let mut l_prop = SlProps::empty(2, tol_dist)?;
        l_prop.set_surface();
        Ok(Self {
            order,
            l_prop,
            point,
            pt2d: GpPnt2d::zero(),
            surf: Some(surf),
            d11,
            d12,
            d21,
            d22,
            d23,
            u,
            v,
            tol_dist,
            tol_ang,
            tol_curv,
            has_pnt2d_on_surf: false,
        })
    }

    /// `SetOrder(Order)` (`cxx:186-190`).
    pub fn set_order(&mut self, order: i32) {
        self.order = order;
    }

    /// `Order()` (`cxx:178-182`).
    pub fn order(&self) -> i32 {
        self.order
    }

    /// `SetG0Criterion(TolDist)` (`cxx:113-116`).
    pub fn set_g0_criterion(&mut self, tol_dist: f64) {
        self.tol_dist = tol_dist;
    }

    /// `SetG1Criterion(TolAng)` (`cxx:121-124`).
    pub fn set_g1_criterion(&mut self, tol_ang: f64) {
        self.tol_ang = tol_ang;
    }

    /// `SetG2Criterion(TolCurv)` (`cxx:129-132`).
    pub fn set_g2_criterion(&mut self, tol_curv: f64) {
        self.tol_curv = tol_curv;
    }

    /// `G0Criterion()` (`cxx:137-141`).
    pub fn g0_criterion(&self) -> f64 {
        self.tol_dist
    }

    /// `G1Criterion()` (`cxx:144-148`).
    pub fn g1_criterion(&self) -> f64 {
        self.tol_ang
    }

    /// `G2Criterion()` (`cxx:152-156`).
    pub fn g2_criterion(&self) -> f64 {
        self.tol_curv
    }

    /// `D0(P)` (`cxx:77-80`).
    pub fn d0(&self) -> GpPnt {
        self.point
    }

    /// `D1(P, V1, V2)` (`cxx:85-92`).
    pub fn d1(&self) -> (GpPnt, GpVec, GpVec) {
        (self.point, self.d11, self.d12)
    }

    /// `D2(P, V1, V2, V3, V4, V5)` (`cxx:97-108`).
    pub fn d2(&self) -> (GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec) {
        (self.point, self.d11, self.d12, self.d21, self.d22, self.d23)
    }

    /// `HasPnt2dOnSurf()` (`cxx:193-197`).
    pub fn has_pnt2d_on_surf(&self) -> bool {
        self.has_pnt2d_on_surf
    }

    /// `SetPnt2dOnSurf(Pnt2d)` (`cxx:201-205`).
    pub fn set_pnt2d_on_surf(&mut self, pnt2d: GpPnt2d) {
        self.pt2d = pnt2d;
        self.has_pnt2d_on_surf = true;
    }

    /// `Pnt2dOnSurf()` (`cxx:210-214`).
    pub fn pnt2d_on_surf(&self) -> GpPnt2d {
        self.pt2d
    }

    /// `LPropSurf()` (`cxx:168-174`): `myLProp.SetParameters(myU, myV)`, then
    /// the props themselves.
    ///
    /// The port's [`SlProps`] does not store the surface, so it travels with
    /// the props; the second element is `None` for a point-only constraint,
    /// which has no surface to evaluate on.
    pub fn l_prop_surf(&mut self) -> (&mut SlProps, Option<&Arc<dyn Surface>>) {
        if let Some(surf) = self.surf.as_ref() {
            self.l_prop.set_parameters(surf.as_ref(), self.u, self.v);
        }
        (&mut self.l_prop, self.surf.as_ref())
    }
}
