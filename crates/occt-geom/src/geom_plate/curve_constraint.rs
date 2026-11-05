//! `GeomPlate_CurveConstraint` — a curve as a surface-deformation constraint.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_CurveConstraint.hxx`,
//! `GeomPlate_CurveConstraint.cxx:39-434` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<Adaptor3d_Curve> Boundary` plus the
//!   `occ::down_cast<Adaptor3d_CurveOnSurface>` null test (`cxx:80-107`) is the
//!   [`BoundaryCurve`] sum type: the `OnSurface` arm is the `myFrontiere`
//!   (pcurve on a surface), the `Curve3d` arm is `my3dCurve`.
//! - `GeomAdaptor_Surface` in the same `down_cast` chain (`cxx:89-106`) is
//!   [`CurveOnSurface::surface`] — the port's `CurveOnSurface` always carries a
//!   `Geom_Surface`, so the `throw` at `cxx:102` cannot fire; the non
//!   `GeomAdaptor_Surface` adaptor (a `BRepAdaptor_Surface`) is the arm OCCT
//!   itself rejects.
//! - `GeomLProp_SLProps myLProp` is [`SlProps`], which takes the surface per
//!   call; [`CurveConstraint::l_prop_surf`] therefore returns both.
//! - `Curve2dOnSurf()` (`cxx:322-345`) needs `Approx_Curve2d`, which is
//!   **UNPORTED** (no `occt-geom` translation of `Approx/Approx_Curve2d.cxx`);
//!   that arm reports the missing class, the `my2dCurve` arm is ported.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom2d::curve::Curve2d;

use crate::curve::Curve;
use crate::curve_on_surface::CurveOnSurface;
use crate::gcpnts::curve_length;
use crate::geom_lprop::sl_props::SlProps;
use crate::law::function::LawFunction;
use crate::surface::Surface;

/// The `Adaptor3d_Curve` a `GeomPlate_CurveConstraint` is built from, after
/// the `occ::down_cast<Adaptor3d_CurveOnSurface>` test of `cxx:80-86`.
#[derive(Clone)]
pub enum BoundaryCurve {
    /// `myFrontiere` (`hxx:145`): a pcurve on a surface.
    OnSurface(Arc<CurveOnSurface>),
    /// `my3dCurve` (`hxx:149`): a free 3D curve.
    Curve3d(Arc<dyn Curve>),
}

/// `class GeomPlate_CurveConstraint` (`hxx:28-167`): "defines curves as
/// constraints to be used to deform a surface".
pub struct CurveConstraint {
    /// `myFrontiere` (`hxx:145`).
    frontiere: Option<Arc<CurveOnSurface>>,
    /// `myNbPoints` (`hxx:146`).
    nb_points: i32,
    /// `myOrder` (`hxx:147`).
    order: i32,
    /// `my3dCurve` (`hxx:149`).
    curve3d: Option<Arc<dyn Curve>>,
    /// `myTang` (`hxx:150`). OCCT leaves it uninitialised in the boundary
    /// constructor (`cxx:60-70`) and never reads it; the port initialises it so
    /// that no uninitialised value can be observed, and keeps it unread.
    #[allow(dead_code)]
    tang: i32,
    /// `my2dCurve` (`hxx:152`): the pcurve on the surface being built.
    curve2d: Option<Arc<dyn Curve2d>>,
    /// `myHCurve2d` (`hxx:153`): the projected pcurve.
    h_curve2d: Option<Arc<dyn Curve2d>>,
    /// `myG0Crit` (`hxx:154`).
    g0_crit: Option<Arc<dyn LawFunction>>,
    /// `myG1Crit` (`hxx:155`).
    g1_crit: Option<Arc<dyn LawFunction>>,
    /// `myG2Crit` (`hxx:156`).
    g2_crit: Option<Arc<dyn LawFunction>>,
    /// `myConstG0` (`hxx:157`).
    const_g0: bool,
    /// `myConstG1` (`hxx:158`).
    const_g1: bool,
    /// `myConstG2` (`hxx:159`).
    const_g2: bool,
    /// `myLProp` (`hxx:160`).
    l_prop: SlProps,
    /// `myTolDist` (`hxx:161`).
    tol_dist: f64,
    /// `myTolAng` (`hxx:162`).
    tol_ang: f64,
    /// `myTolCurv` (`hxx:163`).
    tol_curv: f64,
    /// `myTolU` (`hxx:164`).
    tol_u: f64,
    /// `myTolV` (`hxx:165`).
    tol_v: f64,
}

impl CurveConstraint {
    /// `GeomPlate_CurveConstraint()` (`cxx:40-53`): an empty constraint.
    pub fn new() -> Result<Self, String> {
        Ok(Self {
            frontiere: None,
            nb_points: 0,
            order: 0,
            curve3d: None,
            tang: 0,
            curve2d: None,
            h_curve2d: None,
            g0_crit: None,
            g1_crit: None,
            g2_crit: None,
            const_g0: false,
            const_g1: false,
            const_g2: false,
            // `myLProp(2, 1.e-4)` (`cxx:47`): the `(N, Resolution)` base
            // constructor, whose `myCN` starts at 0.
            l_prop: SlProps::empty(2, 1.0e-4)?,
            tol_dist: 0.0,
            tol_ang: 0.0,
            tol_curv: 0.0,
            tol_u: 0.0,
            tol_v: 0.0,
        })
    }

    /// `GeomPlate_CurveConstraint(Boundary, Order, NPt, TolDist, TolAng,
    /// TolCurv)` (`cxx:58-117`).
    ///
    /// Raises when `Order < -1 || Order > 2` (`cxx:71-75`).
    pub fn from_boundary(
        boundary: BoundaryCurve,
        order: i32,
        nb_points: i32,
        tol_dist: f64,
        tol_ang: f64,
        tol_curv: f64,
    ) -> Result<Self, String> {
        if !(-1..=2).contains(&order) {
            return Err("GeomPlate : The continuity is not G0 G1 or G2".to_string());
        }
        let mut l_prop = SlProps::empty(2, tol_dist)?;
        let (frontiere, curve3d) = match boundary {
            BoundaryCurve::OnSurface(f) => {
                // `myLProp.SetSurface(Surf)` (`cxx:104`) raises `myCN` to 4.
                l_prop.set_surface();
                (Some(f), None)
            }
            BoundaryCurve::Curve3d(c) => (None, Some(c)),
        };
        Ok(Self {
            frontiere,
            nb_points,
            order,
            curve3d,
            tang: order,
            // `my2dCurve.Nullify()` / `myHCurve2d.Nullify()` (`cxx:108-109`).
            curve2d: None,
            h_curve2d: None,
            g0_crit: None,
            g1_crit: None,
            g2_crit: None,
            const_g0: true,
            const_g1: true,
            const_g2: true,
            l_prop,
            tol_dist,
            tol_ang,
            tol_curv,
            tol_u: 0.0,
            tol_v: 0.0,
        })
    }

    /// `FirstParameter()` (`cxx:120-133`).
    pub fn first_parameter(&self) -> f64 {
        if let Some(c) = self.h_curve2d.as_ref() {
            return c.first_parameter();
        }
        match self.curve3d.as_ref() {
            None => self
                .frontiere
                .as_ref()
                .map_or(f64::NAN, |f| Curve::first_parameter(f.as_ref())),
            Some(c) => c.first_parameter(),
        }
    }

    /// `LastParameter()` (`cxx:138-153`).
    pub fn last_parameter(&self) -> f64 {
        if let Some(c) = self.h_curve2d.as_ref() {
            return c.last_parameter();
        }
        match self.curve3d.as_ref() {
            None => self
                .frontiere
                .as_ref()
                .map_or(f64::NAN, |f| Curve::last_parameter(f.as_ref())),
            Some(c) => c.last_parameter(),
        }
    }

    /// `Length()` (`cxx:158-176`): `GCPnts_AbscissaPoint::Length` on the
    /// frontier, or on the free 3D curve.
    ///
    /// `Err` is the empty constraint, whose both handles are null in OCCT
    /// (a null dereference, not a throw).
    pub fn length(&self) -> Result<f64, String> {
        if let Some(f) = self.frontiere.as_ref() {
            return Ok(curve_length(f.as_ref() as &dyn Curve));
        }
        match self.curve3d.as_ref() {
            Some(c) => Ok(curve_length(c.as_ref())),
            None => Err("GeomPlate_CurveConstraint::Length: no curve".to_string()),
        }
    }

    /// `D0(U, P)` (`cxx:182-193`).
    pub fn d0(&self, u: f64) -> GpPnt {
        if let (None, Some(f)) = (self.curve3d.as_ref(), self.frontiere.as_ref()) {
            let p2d = f.curve().d0(u);
            return f.surface().d0(p2d.x(), p2d.y());
        }
        match self.curve3d.as_ref() {
            Some(c) => c.d0(u),
            None => GpPnt::zero(),
        }
    }

    /// `D1(U, P, V1, V2)` (`cxx:198-209`).
    ///
    /// Raises when the constraint is a free 3D curve (`cxx:202-205`).
    pub fn d1(&self, u: f64) -> Result<(GpPnt, GpVec, GpVec), String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        match self.frontiere.as_ref() {
            Some(f) => {
                let p2d = f.curve().d0(u);
                Ok(f.surface().d1(p2d.x(), p2d.y()))
            }
            None => Err("GeomPlate_CurveConstraint::D1: no curve".to_string()),
        }
    }

    /// `D2(U, P, V1, V2, V3, V4, V5)` (`cxx:213-230`).
    ///
    /// Raises when the constraint is a free 3D curve (`cxx:224-227`).
    pub fn d2(&self, u: f64) -> Result<(GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec), String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        match self.frontiere.as_ref() {
            Some(f) => {
                let p2d = f.curve().d0(u);
                Ok(f.surface().d2(p2d.x(), p2d.y()))
            }
            None => Err("GeomPlate_CurveConstraint::D2: no curve".to_string()),
        }
    }

    /// `SetG0Criterion(G0Crit)` (`cxx:236-239`).
    pub fn set_g0_criterion(&mut self, crit: Arc<dyn LawFunction>) {
        self.g0_crit = Some(crit);
        self.const_g0 = false;
    }

    /// `SetG1Criterion(G1Crit)` (`cxx:244-252`).
    pub fn set_g1_criterion(&mut self, crit: Arc<dyn LawFunction>) -> Result<(), String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        self.g1_crit = Some(crit);
        self.const_g1 = false;
        Ok(())
    }

    /// `SetG2Criterion(G2Crit)` (`cxx:257-265`).
    pub fn set_g2_criterion(&mut self, crit: Arc<dyn LawFunction>) -> Result<(), String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        self.g2_crit = Some(crit);
        self.const_g2 = false;
        Ok(())
    }

    /// `G0Criterion(U)` (`cxx:271-281`).
    pub fn g0_criterion(&self, u: f64) -> Result<f64, String> {
        if self.const_g0 {
            return Ok(self.tol_dist);
        }
        self.g0_crit
            .as_ref()
            .map(|c| c.value(u))
            .ok_or_else(|| "GeomPlate_CurveConstraint::G0Criterion: no law".to_string())
    }

    /// `G1Criterion(U)` (`cxx:285-298`).
    pub fn g1_criterion(&self, u: f64) -> Result<f64, String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        if self.const_g1 {
            return Ok(self.tol_ang);
        }
        self.g1_crit
            .as_ref()
            .map(|c| c.value(u))
            .ok_or_else(|| "GeomPlate_CurveConstraint::G1Criterion: no law".to_string())
    }

    /// `G2Criterion(U)` (`cxx:304-316`).
    pub fn g2_criterion(&self, u: f64) -> Result<f64, String> {
        if self.curve3d.is_some() {
            return Err("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string());
        }
        if self.const_g2 {
            return Ok(self.tol_curv);
        }
        self.g2_crit
            .as_ref()
            .map(|c| c.value(u))
            .ok_or_else(|| "GeomPlate_CurveConstraint::G2Criterion: no law".to_string())
    }

    /// `Curve2dOnSurf()` (`cxx:322-345`).
    ///
    /// The `my2dCurve.IsNull() && !myHCurve2d.IsNull()` arm builds the pcurve
    /// with `Approx_Curve2d`, which has no port; it reports `UNPORTED` and
    /// otherwise returns `my2dCurve`.
    pub fn curve2d_on_surf(&self) -> Result<Option<Arc<dyn Curve2d>>, String> {
        if self.curve2d.is_none() && self.h_curve2d.is_some() {
            return Err(
                "UNPORTED: GeomPlate_CurveConstraint.cxx:329-340 Approx_Curve2d".to_string(),
            );
        }
        Ok(self.curve2d.clone())
    }

    /// `SetCurve2dOnSurf(Curve)` (`cxx:352-355`).
    pub fn set_curve2d_on_surf(&mut self, curve: Arc<dyn Curve2d>) {
        self.curve2d = Some(curve);
    }

    /// `SetCurve2dOnSurf(NullCurve)` (`GeomPlate_BuildPlateSurface.cxx:1898-1901`):
    /// stores a null handle, i.e. drops the p-curve.
    pub fn clear_curve2d_on_surf(&mut self) {
        self.curve2d = None;
    }

    /// `ProjectedCurve()` (`cxx:360-363`).
    pub fn projected_curve(&self) -> Option<Arc<dyn Curve2d>> {
        self.h_curve2d.clone()
    }

    /// `SetProjectedCurve(Curve2d, TolU, TolV)` (`cxx:367-373`).
    pub fn set_projected_curve(&mut self, curve: Arc<dyn Curve2d>, tol_u: f64, tol_v: f64) {
        self.h_curve2d = Some(curve);
        self.tol_u = tol_u;
        self.tol_v = tol_v;
    }

    /// `Curve3d()` (`cxx:379-388`): `my3dCurve` when set, else the frontier
    /// itself as an `Adaptor3d_Curve`.
    pub fn curve3d(&self) -> Option<Arc<dyn Curve>> {
        match self.curve3d.as_ref() {
            Some(c) => Some(c.clone()),
            None => self
                .frontiere
                .as_ref()
                .map(|f| f.clone() as Arc<dyn Curve>),
        }
    }

    /// `NbPoints()` (`cxx:397-400`).
    pub fn nb_points(&self) -> i32 {
        self.nb_points
    }

    /// `SetNbPoints(NewNb)` (`cxx:413-416`).
    pub fn set_nb_points(&mut self, new_nb: i32) {
        self.nb_points = new_nb;
    }

    /// `Order()` (`cxx:405-408`).
    pub fn order(&self) -> i32 {
        self.order
    }

    /// `SetOrder(Order)` (`cxx:421-424`).
    pub fn set_order(&mut self, order: i32) {
        self.order = order;
    }

    /// `LPropSurf(U)` (`cxx:428-434`): set the parameters on the frontier's
    /// pcurve, then the props.
    ///
    /// The port's [`SlProps`] does not store the surface, so it travels with
    /// the props. Raises when the curve is not on a surface (`cxx:429-432`).
    pub fn l_prop_surf(&mut self, u: f64) -> Result<(&mut SlProps, &Arc<dyn Surface>), String> {
        let f = self
            .frontiere
            .as_ref()
            .ok_or_else(|| "GeomPlate_CurveConstraint.cxx : Curve must be on a Surface".to_string())?;
        let p2d = f.curve().d0(u);
        let surf = f.surface();
        self.l_prop.set_parameters(surf.as_ref(), p2d.x(), p2d.y());
        Ok((&mut self.l_prop, surf))
    }
}
