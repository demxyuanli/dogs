//! `LocalAnalysis_SurfaceContinuity` — C0/C1/C2/G1/G2 continuity between two
//! points that sit on two surfaces.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/LocalAnalysis/LocalAnalysis_SurfaceContinuity.cxx:27-822`
//! and `LocalAnalysis_SurfaceContinuity.hxx:33-241`.
//!
//! | OCCT | here |
//! |---|---|
//! | `SurfC0` (`cxx:28-32`) | [`SurfaceContinuity::surf_c0`] |
//! | `SurfC1` (`cxx:36-92`) | [`SurfaceContinuity::surf_c1`] |
//! | `SurfC2` (`cxx:98-194`) | [`SurfaceContinuity::surf_c2`] |
//! | `SurfG1` (`cxx:197-217`) | [`SurfaceContinuity::surf_g1`] |
//! | `SurfG2` (`cxx:222-266`) | [`SurfaceContinuity::surf_g2`] |
//! | ctor `(EpsNul, ..., Maxlen)` (`cxx:268-303`) | [`SurfaceContinuity::new`] |
//! | `ComputeAnalysis` (`cxx:305-340`) | [`SurfaceContinuity::compute_analysis`] |
//! | ctor `(Surf1, u1, v1, Surf2, u2, v2, ...)` (`cxx:344-428`) | [`SurfaceContinuity::from_points`] |
//! | ctor `(curv1, curv2, U, Surf1, Surf2, ...)` (`cxx:433-537`) | [`SurfaceContinuity::from_curves`] |
//! | `IsC0` / `IsC1` / `IsC2` / `IsG1` / `IsG2` (`cxx:540-678`) | `is_c0` / `is_c1` / `is_c2` / `is_g1` / `is_g2` |
//! | `ContinuityStatus` (`cxx:682-688`) | [`SurfaceContinuity::continuity_status`] |
//! | the measured-value accessors (`cxx:692-776`) | [`SurfaceContinuity::c0_value`] and the `c1u_angle` / `c1v_angle` / `c2u_angle` / `c2v_angle` / `g1_angle` / `c1u_ratio` / `c2u_ratio` / `c1v_ratio` / `c2v_ratio` family |
//! | `G2CurvatureGap` (`cxx:802-808`) | [`SurfaceContinuity::g2_curvature_gap`] |
//! | `IsDone` / `StatusError` (`cxx:811-822`) | [`SurfaceContinuity::is_done`] / [`SurfaceContinuity::status_error`] |
//!
//! Translation notes:
//! - `GeomLProp_SLProps` is [`occt_geom::geom_lprop::SlProps`], whose surface is
//!   a per-call argument instead of a stored handle
//!   (`occt-geom/src/geom_lprop/sl_props.rs:419-424`). `ComputeAnalysis` and the
//!   five `Surf*` helpers therefore take `(&mut SlProps, &dyn Surface)` pairs;
//!   the arithmetic is unchanged.
//! - `StdFail_NotDone` (`cxx:543`, `:558`, `:565`, ..., `:806`) becomes
//!   `Err("StdFail_NotDone")`, as in the rest of this tree.
//! - `gp_Dir(X, Y, Z)` raises `Standard_ConstructionError` on a null vector
//!   (`gp_Dir.cxx`); [`GpDir::new`] returns `Err` and that `Err` propagates.
//! - `myErrorStatus` is *not* in any constructor's initialiser list
//!   (`cxx:270-303`, `:360-378`, `:446-464`) and `myTypeCont` is absent from the
//!   parameter-only constructor (`cxx:268-303`). Both are seeded with the enum's
//!   first value (`NullFirstDerivative`, `C0`) here so the port has no
//!   uninitialised read; every path that is observable after `IsDone()` returned
//!   a value writes them first.
//! - `M_PI / 2` and `M_PI` (`cxx:71-78`, `:85-88`, ...) are [`std::f64::consts::FRAC_PI_2`]
//!   and [`std::f64::consts::PI`].

use occt_core::gp::GpDir;
use occt_core::kernel::geomabs::Shape;
use occt_geom2d::Curve2d;

use crate::Surface;
use crate::geom_lprop::SlProps;

use std::f64::consts::{FRAC_PI_2, PI};

/// `LocalAnalysis_StatusErrorType` (`LocalAnalysis_StatusErrorType.hxx:21-29`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusErrorType {
    /// `LocalAnalysis_NullFirstDerivative`
    NullFirstDerivative,
    /// `LocalAnalysis_NullSecondDerivative`
    NullSecondDerivative,
    /// `LocalAnalysis_TangentNotDefined`
    TangentNotDefined,
    /// `LocalAnalysis_NormalNotDefined`
    NormalNotDefined,
    /// `LocalAnalysis_CurvatureNotDefined`
    CurvatureNotDefined,
}

/// `LocalAnalysis_SurfaceContinuity` (`hxx:33-241`).
#[derive(Debug, Clone)]
pub struct SurfaceContinuity {
    /// `myContC0`
    my_cont_c0: f64,
    /// `myContC1U`
    my_cont_c1u: f64,
    /// `myContC1V`
    my_cont_c1v: f64,
    /// `myContC2U`
    my_cont_c2u: f64,
    /// `myContC2V`
    my_cont_c2v: f64,
    /// `myContG1`
    my_cont_g1: f64,
    /// `myLambda1U`
    my_lambda1u: f64,
    /// `myLambda2U`
    my_lambda2u: f64,
    /// `myLambda1V`
    my_lambda1v: f64,
    /// `myLambda2V`
    my_lambda2v: f64,
    /// `myETA1`
    my_eta1: f64,
    /// `myETA2`
    my_eta2: f64,
    /// `myETA`
    my_eta: f64,
    /// `myZETA1`
    my_zeta1: f64,
    /// `myZETA2`
    my_zeta2: f64,
    /// `myZETA`
    my_zeta: f64,
    /// `myAlpha`
    my_alpha: f64,
    /// `myTypeCont`
    my_type_cont: Shape,
    /// `myepsC0`
    my_eps_c0: f64,
    /// `myepsnul`
    my_eps_nul: f64,
    /// `myepsC1`
    my_eps_c1: f64,
    /// `myepsC2`
    my_eps_c2: f64,
    /// `myepsG1`
    my_eps_g1: f64,
    /// `myperce`
    my_perce: f64,
    /// `mymaxlen`
    my_max_len: f64,
    /// `myGap`
    my_gap: f64,
    /// `myIsDone`
    my_is_done: bool,
    /// `myErrorStatus`
    my_error_status: StatusErrorType,
}

impl SurfaceContinuity {
    /// `LocalAnalysis_SurfaceContinuity(EpsNul, EpsC0, EpsC1, EpsC2, EpsG1,
    /// Percent, Maxlen)` (`cxx:268-303`): the study is left for
    /// [`Self::compute_analysis`].
    pub fn new(
        eps_nul: f64,
        eps_c0: f64,
        eps_c1: f64,
        eps_c2: f64,
        eps_g1: f64,
        percent: f64,
        max_len: f64,
    ) -> Self {
        let mut s = Self::blank();
        s.my_eps_nul = eps_nul;
        s.my_eps_c0 = eps_c0;
        s.my_eps_c1 = eps_c1;
        s.my_eps_c2 = eps_c2;
        s.my_eps_g1 = eps_g1;
        s.my_perce = percent;
        s.my_max_len = max_len;
        s.my_is_done = true;
        s
    }

    /// `LocalAnalysis_SurfaceContinuity(Surf1, u1, v1, Surf2, u2, v2, Ordre,
    /// ...)` (`cxx:344-428`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_points(
        surf1: &dyn Surface,
        u1: f64,
        v1: f64,
        surf2: &dyn Surface,
        u2: f64,
        v2: f64,
        order: Shape,
        eps_nul: f64,
        eps_c0: f64,
        eps_c1: f64,
        eps_c2: f64,
        eps_g1: f64,
        percent: f64,
        max_len: f64,
    ) -> Result<Self, String> {
        let mut s = Self::new(eps_nul, eps_c0, eps_c1, eps_c2, eps_g1, percent, max_len);
        s.my_type_cont = order;
        // `GeomLProp_SLProps Surfa1(Surf1, u1, v1, N, myepsnul)` — the derivative
        // order N is the one the switch arm asks for.
        let n = Self::deriv_order(order)?;
        let mut a1 = SlProps::new(surf1, u1, v1, n, s.my_eps_nul)?;
        let mut a2 = SlProps::new(surf2, u2, v2, n, s.my_eps_nul)?;
        s.dispatch(&mut a1, surf1, &mut a2, surf2, order)?;
        Ok(s)
    }

    /// `LocalAnalysis_SurfaceContinuity(curv1, curv2, U, Surf1, Surf2, Ordre,
    /// ...)` (`cxx:433-537`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_curves(
        curv1: &dyn Curve2d,
        curv2: &dyn Curve2d,
        u: f64,
        surf1: &dyn Surface,
        surf2: &dyn Surface,
        order: Shape,
        eps_nul: f64,
        eps_c0: f64,
        eps_c1: f64,
        eps_c2: f64,
        eps_g1: f64,
        percent: f64,
        max_len: f64,
    ) -> Result<Self, String> {
        let mut s = Self::new(eps_nul, eps_c0, eps_c1, eps_c2, eps_g1, percent, max_len);
        s.my_type_cont = order;
        let par_d1 = curv1.first_parameter();
        let par_d2 = curv2.first_parameter();
        let par_f1 = curv1.last_parameter();
        let par_f2 = curv2.last_parameter();
        if u > par_f1 || u < par_d1 || u > par_f2 || u < par_d2 {
            // `cxx:485-487`: `myIsDone = false` with no error status.
            s.my_is_done = false;
            return Ok(s);
        }
        let pt1 = curv1.value(u);
        let pt2 = curv2.value(u);
        let (u1, v1) = (pt1.x(), pt1.y());
        let (u2, v2) = (pt2.x(), pt2.y());
        let n = Self::deriv_order(order)?;
        let mut a1 = SlProps::new(surf1, u1, v1, n, s.my_eps_nul)?;
        let mut a2 = SlProps::new(surf2, u2, v2, n, s.my_eps_nul)?;
        s.dispatch(&mut a1, surf1, &mut a2, surf2, order)?;
        Ok(s)
    }

    /// The `GeomLProp_SLProps(S, U, V, N, Tol)` order each `GeomAbs_Shape` arm
    /// builds (`cxx:392`, `:397`, `:405-406`, `:412-413`, `:419-420`).
    fn deriv_order(order: Shape) -> Result<i32, String> {
        Ok(match order {
            Shape::C0 => 0,
            Shape::C1 | Shape::G1 => 1,
            Shape::C2 | Shape::G2 => 2,
            // No constructor arm exists for C3 / CN; the `default:` arm
            // (`cxx:426-427`) builds nothing and touches no measured value.
            Shape::C3 | Shape::CN => 0,
        })
    }

    /// `ComputeAnalysis(Surf1, Surf2, Order)` (`cxx:305-340`), plus the arm
    /// bodies of the two geometry constructors (`cxx:389-427`, `:494-535`).
    pub fn compute_analysis(
        &mut self,
        surf1: &mut SlProps,
        s1: &dyn Surface,
        surf2: &mut SlProps,
        s2: &dyn Surface,
        order: Shape,
    ) -> Result<(), String> {
        self.my_type_cont = order;
        self.dispatch(surf1, s1, surf2, s2, order)
    }

    fn dispatch(
        &mut self,
        a1: &mut SlProps,
        s1: &dyn Surface,
        a2: &mut SlProps,
        s2: &dyn Surface,
        order: Shape,
    ) -> Result<(), String> {
        match order {
            Shape::C0 => self.surf_c0(a1, a2),
            Shape::C1 => {
                self.surf_c0(a1, a2);
                self.surf_c1(a1, s1, a2, s2);
            }
            Shape::C2 => {
                self.surf_c0(a1, a2);
                self.surf_c1(a1, s1, a2, s2);
                self.surf_c2(a1, s1, a2, s2);
            }
            Shape::G1 => {
                self.surf_c0(a1, a2);
                self.surf_g1(a1, s1, a2, s2);
            }
            Shape::G2 => {
                self.surf_c0(a1, a2);
                self.surf_g1(a1, s1, a2, s2);
                self.surf_g2(a1, s1, a2, s2)?;
            }
            // `default: {}` (`cxx:338-339`).
            Shape::C3 | Shape::CN => {}
        }
        Ok(())
    }

    /// `SurfC0` (`cxx:28-32`): the gap between the two evaluated points.
    pub fn surf_c0(&mut self, surf1: &SlProps, surf2: &SlProps) {
        self.my_cont_c0 = surf1.value().distance(&surf2.value());
    }

    /// `SurfC1` (`cxx:36-92`).
    pub fn surf_c1(
        &mut self,
        surf1: &mut SlProps,
        s1: &dyn Surface,
        surf2: &mut SlProps,
        s2: &dyn Surface,
    ) {
        let v1u = surf1.d1u(s1);
        let v2u = surf2.d1u(s2);
        let v1v = surf1.d1v(s1);
        let v2v = surf2.d1v(s2);

        let norm1u = v1u.magnitude();
        let norm2u = v2u.magnitude();
        let norm1v = v1v.magnitude();
        let norm2v = v2v.magnitude();

        if norm1u > self.my_eps_nul
            && norm2u > self.my_eps_nul
            && norm1v > self.my_eps_nul
            && norm2v > self.my_eps_nul
        {
            self.my_lambda1u = if norm1u >= norm2u {
                norm2u / norm1u
            } else {
                norm1u / norm2u
            };
            self.my_lambda1v = if norm1v >= norm2v {
                norm2v / norm1v
            } else {
                norm1v / norm2v
            };
            let angu = v1u.angle(&v2u);
            self.my_cont_c1u = if angu > FRAC_PI_2 { PI - angu } else { angu };
            let angv = v1v.angle(&v2v);
            self.my_cont_c1v = if angv > FRAC_PI_2 { PI - angv } else { angv };
        } else {
            self.my_is_done = false;
            self.my_error_status = StatusErrorType::NullFirstDerivative;
        }
    }

    /// `SurfC2` (`cxx:98-194`).
    pub fn surf_c2(
        &mut self,
        surf1: &mut SlProps,
        s1: &dyn Surface,
        surf2: &mut SlProps,
        s2: &dyn Surface,
    ) {
        let v11u = surf1.d1u(s1);
        let v12u = surf2.d1u(s2);
        let v21u = surf1.d2u(s1);
        let v22u = surf2.d2u(s2);
        let norm11u = v11u.magnitude();
        let norm12u = v12u.magnitude();
        let norm21u = v21u.magnitude();
        let norm22u = v22u.magnitude();

        if norm11u > self.my_eps_nul && norm12u > self.my_eps_nul {
            if norm21u > self.my_eps_nul && norm22u > self.my_eps_nul {
                if norm11u >= norm12u {
                    self.my_lambda1u = norm12u / norm11u;
                    self.my_lambda2u = norm22u / norm21u;
                } else {
                    self.my_lambda1u = norm11u / norm12u;
                    self.my_lambda2u = norm21u / norm22u;
                }
                let ang = v21u.angle(&v22u);
                self.my_cont_c2u = if ang > FRAC_PI_2 { PI - ang } else { ang };
            } else {
                self.my_is_done = false;
                self.my_error_status = StatusErrorType::NullSecondDerivative;
            }
        } else {
            self.my_is_done = false;
            self.my_error_status = StatusErrorType::NullFirstDerivative;
        }

        let v11v = surf1.d1v(s1);
        let v12v = surf2.d1v(s2);
        let v21v = surf1.d2v(s1);
        let v22v = surf2.d2v(s2);
        let norm11v = v11v.magnitude();
        let norm12v = v12v.magnitude();
        let norm21v = v21v.magnitude();
        let norm22v = v22v.magnitude();

        if norm11v > self.my_eps_nul && norm12v > self.my_eps_nul {
            if norm21v > self.my_eps_nul && norm22v > self.my_eps_nul {
                if norm11v >= norm12v {
                    self.my_lambda1v = norm12v / norm11v;
                    self.my_lambda2v = norm22v / norm21v;
                } else {
                    self.my_lambda1v = norm11v / norm12v;
                    self.my_lambda2v = norm21v / norm22v;
                }
                let ang = v21v.angle(&v22v);
                self.my_cont_c2v = if ang > FRAC_PI_2 { PI - ang } else { ang };
            } else {
                self.my_is_done = false;
                self.my_error_status = StatusErrorType::NullSecondDerivative;
            }
        } else {
            self.my_is_done = false;
            self.my_error_status = StatusErrorType::NullFirstDerivative;
        }
    }

    /// `SurfG1` (`cxx:197-217`).
    pub fn surf_g1(
        &mut self,
        surf1: &mut SlProps,
        s1: &dyn Surface,
        surf2: &mut SlProps,
        s2: &dyn Surface,
    ) {
        if surf1.is_normal_defined(s1) && surf2.is_normal_defined(s2) {
            let d1 = surf1.normal(s1).unwrap_or_else(|_| GpDir::default_dir());
            let d2 = surf2.normal(s2).unwrap_or_else(|_| GpDir::default_dir());
            let ang = d1.angle(&d2);
            self.my_cont_g1 = if ang > FRAC_PI_2 { PI - ang } else { ang };
        } else {
            self.my_is_done = false;
            self.my_error_status = StatusErrorType::NormalNotDefined;
        }
    }

    /// `SurfG2` (`cxx:222-266`).
    pub fn surf_g2(
        &mut self,
        surf1: &mut SlProps,
        s1: &dyn Surface,
        surf2: &mut SlProps,
        s2: &dyn Surface,
    ) -> Result<(), String> {
        if surf1.is_curvature_defined(s1) && surf2.is_curvature_defined(s2) {
            let (dmax1, dmin1) = surf1.curvature_directions(s1)?;
            let (dmax2, dmin2) = surf2.curvature_directions(s2)?;
            let (x1, y1, z1) = (dmin1.x(), dmin1.y(), dmin1.z());
            let (x2, y2, z2) = (dmax1.x(), dmax1.y(), dmax1.z());
            let mcd1 = GpDir::new(
                (x1.abs() + x2.abs()) / 2.0,
                (y1.abs() + y2.abs()) / 2.0,
                (z1.abs() + z2.abs()) / 2.0,
            )
            .map_err(|e| e.to_string())?;
            let (x1, y1, z1) = (dmin2.x(), dmin2.y(), dmin2.z());
            let (x2, y2, z2) = (dmax2.x(), dmax2.y(), dmax2.z());
            let mcd2 = GpDir::new(
                (x1.abs() + x2.abs()) / 2.0,
                (y1.abs() + y2.abs()) / 2.0,
                (z1.abs() + z2.abs()) / 2.0,
            )
            .map_err(|e| e.to_string())?;

            self.my_alpha = mcd1.angle(&mcd2);
            let rmin1 = surf1.min_curvature(s1)?;
            let rmax1 = surf1.max_curvature(s1)?;
            let rmin2 = surf2.min_curvature(s2)?;
            let rmax2 = surf2.max_curvature(s2)?;
            self.my_eta1 = (rmin1 + rmax1) / 2.0;
            self.my_eta2 = (rmin2 + rmax2) / 2.0;
            self.my_eta = (self.my_eta1 + self.my_eta2) / 2.0;
            self.my_zeta1 = (rmax1 - rmin1) / 2.0;
            self.my_zeta2 = (rmax2 - rmin2) / 2.0;
            self.my_zeta = (self.my_zeta1 + self.my_zeta2) / 2.0;
            let deta = (self.my_eta1 - self.my_eta2) / 2.0;
            let dzeta = (self.my_zeta1 - self.my_zeta2) / 2.0;
            self.my_gap = deta.abs()
                + (dzeta * dzeta * self.my_alpha.cos() * self.my_alpha.cos()
                    + self.my_zeta * self.my_zeta * self.my_alpha.sin() * self.my_alpha.sin())
                .sqrt();
        } else {
            self.my_is_done = false;
            self.my_error_status = StatusErrorType::CurvatureNotDefined;
        }
        Ok(())
    }

    /// `IsC0()` (`cxx:540-548`).
    pub fn is_c0(&self) -> Result<bool, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c0 <= self.my_eps_c0)
    }

    /// `IsC1()` (`cxx:552-558`).
    pub fn is_c1(&self) -> Result<bool, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.is_c0()? && self.my_cont_c1u <= self.my_eps_c1 && self.my_cont_c1v <= self.my_eps_c1)
    }

    /// `IsC2()` (`cxx:562-602`).
    pub fn is_c2(&self) -> Result<bool, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        if !self.is_c1()? {
            return Ok(false);
        }
        let eps1u = 0.5 * self.my_eps_c1 * self.my_eps_c1 * self.my_lambda1u;
        let eps1v = 0.5 * self.my_eps_c1 * self.my_eps_c1 * self.my_lambda1v;
        let eps2u = 0.5 * self.my_eps_c2 * self.my_eps_c2 * self.my_lambda2u;
        let eps2v = 0.5 * self.my_eps_c2 * self.my_eps_c2 * self.my_lambda2v;
        if self.my_cont_c2u < self.my_eps_c2 && self.my_cont_c2v < self.my_eps_c2 {
            if (self.my_lambda1u * self.my_lambda1u - self.my_lambda2u).abs()
                <= (eps1u * eps1u + eps2u)
            {
                Ok((self.my_lambda1v * self.my_lambda1v - self.my_lambda2v).abs()
                    <= (eps1v * eps1v + eps2v))
            } else {
                Ok(false)
            }
        } else {
            Ok(false)
        }
    }

    /// `IsG1()` (`cxx:606-614`).
    pub fn is_g1(&self) -> Result<bool, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.is_c0()? && self.my_cont_g1 <= self.my_eps_g1)
    }

    /// `IsG2()` (`cxx:618-678`).
    pub fn is_g2(&self) -> Result<bool, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        let mut itype = 0;
        let epsnl = 8.0 * self.my_eps_c0 / (self.my_max_len * self.my_max_len);
        if !self.is_g1()? {
            return Ok(false);
        }
        if self.my_eta.abs() < epsnl && self.my_zeta.abs() < epsnl {
            return Ok(true);
        }
        if self.my_zeta1.abs() < epsnl && self.my_zeta2.abs() < epsnl {
            itype = 1;
        } else if self.my_eta1.abs() < epsnl && self.my_eta2.abs() < epsnl {
            itype = 1;
        } else if (self.my_zeta.abs() - self.my_eta.abs()).abs() < epsnl {
            itype = 1;
        } else if self.my_eta1 < self.my_zeta1 && self.my_eta2 < self.my_zeta2 {
            itype = 1;
        } else if self.my_eta1 > self.my_zeta1 && self.my_eta2 > self.my_zeta2 {
            itype = 1;
        }
        if itype == 1 {
            if self.my_eta >= 2.0 * self.my_zeta
                && self.my_gap <= self.my_perce * (self.my_eta - self.my_zeta)
            {
                return Ok(true);
            }
            if self.my_zeta >= self.my_eta && self.my_gap <= self.my_perce * self.my_zeta {
                return Ok(true);
            }
            return Ok(self.my_zeta <= self.my_eta
                && self.my_eta <= 2.0 * self.my_zeta
                && self.my_gap <= self.my_perce * self.my_eta);
        }
        Ok(false)
    }

    /// `ContinuityStatus()` (`cxx:682-688`).
    pub fn continuity_status(&self) -> Result<Shape, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_type_cont)
    }

    /// `C0Value()` (`cxx:692-698`).
    pub fn c0_value(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c0)
    }

    /// `C1UAngle()` (`cxx:702-710`).
    pub fn c1u_angle(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c1u)
    }

    /// `C1VAngle()` (`cxx:714-722`).
    pub fn c1v_angle(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c1v)
    }

    /// `C2UAngle()` (`cxx:726-732`).
    pub fn c2u_angle(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c2u)
    }

    /// `C2VAngle()` (`cxx:736-744`).
    pub fn c2v_angle(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_c2v)
    }

    /// `G1Angle()` (`cxx:748-756`).
    pub fn g1_angle(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_cont_g1)
    }

    /// `C1URatio()` (`cxx:760-766`).
    pub fn c1u_ratio(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_lambda1u)
    }

    /// `C2URatio()` (`cxx:770-776`).
    pub fn c2u_ratio(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_lambda2u)
    }

    /// `C1VRatio()` (`cxx:780-786`).
    pub fn c1v_ratio(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_lambda1v)
    }

    /// `C2VRatio()` (`cxx:790-796`).
    pub fn c2v_ratio(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_lambda2v)
    }

    /// `G2CurvatureGap()` (`cxx:802-808`).
    pub fn g2_curvature_gap(&self) -> Result<f64, String> {
        if !self.my_is_done {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.my_gap)
    }

    /// `IsDone()` (`cxx:814-817`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `StatusError()` (`cxx:819-822`).
    pub fn status_error(&self) -> StatusErrorType {
        self.my_error_status
    }

    /// The eighteen measured values, all zero, as in the initialiser lists at
    /// `cxx:276-293`, `:360-377` and `:448-465`.
    fn blank() -> Self {
        Self {
            my_cont_c0: 0.0,
            my_cont_c1u: 0.0,
            my_cont_c1v: 0.0,
            my_cont_c2u: 0.0,
            my_cont_c2v: 0.0,
            my_cont_g1: 0.0,
            my_lambda1u: 0.0,
            my_lambda2u: 0.0,
            my_lambda1v: 0.0,
            my_lambda2v: 0.0,
            my_eta1: 0.0,
            my_eta2: 0.0,
            my_eta: 0.0,
            my_zeta1: 0.0,
            my_zeta2: 0.0,
            my_zeta: 0.0,
            my_alpha: 0.0,
            // Not initialised by any OCCT constructor; see the module notes.
            my_type_cont: Shape::C0,
            my_eps_c0: 0.0,
            my_eps_nul: 0.0,
            my_eps_c1: 0.0,
            my_eps_c2: 0.0,
            my_eps_g1: 0.0,
            my_perce: 0.0,
            my_max_len: 0.0,
            my_gap: 0.0,
            my_is_done: false,
            // Not initialised by any OCCT constructor; see the module notes.
            my_error_status: StatusErrorType::NullFirstDerivative,
        }
    }
}
