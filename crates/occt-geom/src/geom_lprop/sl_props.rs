//! `GeomLProp_SLPropsBase` — local properties of a 3D surface at one
//! `(u, v)`.
//!
//! Source: `ModelingData/TKGeomBase/GeomLProp/GeomLProp_SLProps.hxx:37-315`
//! (the `GeomLProp_SLPropsBase` template, instantiated as
//! `GeomLProp_SLProps` with `occ::handle<Geom_Surface>`) plus the utility
//! namespace `GeomLProp_SurfaceUtils.hxx` (OCCT 8.0.0):
//!
//! | OCCT | here |
//! |---|---|
//! | `SetParameters` (`hxx:488-511`) | [`SlProps::set_parameters`] |
//! | `EvalSurfDerivatives` (`hxx:186-210`) | [`eval_surf_derivatives`] |
//! | `EnsureSurfDeriv` (`hxx:513-534`) | [`SlProps::ensure_surf_deriv`] |
//! | `FindSurfTangentOrder` (`hxx:220-251`) | [`find_surf_tangent_order`] |
//! | `ComputeSurfTangent` (`hxx:263-343`) | [`compute_surf_tangent`] |
//! | `ComputeSurfNormal` (`hxx:346-358`) | [`compute_surf_normal`] |
//! | `ComputeSurfCurvatures` (`hxx:363-481`) | [`compute_surf_curvatures`] |
//! | `IsTangentU/VDefined` (`hxx:536-573`) | [`SlProps::is_tangent_u_defined`] / `_v_` |
//! | `TangentU` / `TangentV` (`hxx:575-607`) | [`SlProps::tangent_u`] / `tangent_v` |
//! | `D1U` / `D1V` / `D2U` / `D2V` (`hxx:145-174`) | [`SlProps::d1u`] / `d1v` / `d2u` / `d2v` |
//! | `IsNormalDefined` / `Normal` (`hxx:609-636`) | [`SlProps::is_normal_defined`] / [`SlProps::normal`] |
//! | `IsCurvatureDefined` (`hxx:638-696`) | [`SlProps::is_curvature_defined`] |
//! | `RequireCurvature` (`hxx:699-704`) | the `max_curvature` / `min_curvature` / `mean_curvature` / `gaussian_curvature` accessors |
//! | `IsUmbilic` (`hxx:706-710`) | [`SlProps::is_umbilic`] |
//! | `CurvatureDirections` (`hxx:712-728`) | [`SlProps::curvature_directions`] |
//!
//! Translation notes:
//! - The `LProp_NotDefined` throw behind `RequireCurvature` / `Normal` /
//!   `TangentU` / `TangentV` / `IsUmbilic` / `CurvatureDirections` becomes
//!   `Err(String)`, as in [`crate::geom_lprop::cl_props`].
//! - The `DirectAccess` policy (`hxx:89-128`) is the only one instantiated in
//!   this tree, so the `Access` template parameter is folded into direct
//!   `Surface::d0`/`d1`/`d2` calls. The `GetSurfBounds(Geom_Surface&)`
//!   overload (`hxx:61-70`) and `GetSurfBounds(Adaptor3d_Surface&)`
//!   (`hxx:73-86`) both read `(U1, U2, V1, V2)`; [`Surface::u_range`] /
//!   [`Surface::v_range`] are that pair.
//! - The `myCN` cache counter (`hxx:290`) starts at `4` in the `(S, U, V, N,
//!   Tol)` and `(S, N, Tol)` constructors and at `0` in the `(N, Tol)` one
//!   (`hxx:63`, `:82`, `:99`); `SetSurface` resets it to `4` (`hxx:113-117`).
//!   The port keeps it as `cn`.
//! - `RealLast()` / `RealFirst()` (`hxx:290` of `ComputeSurfTangent`,
//!   `Standard_Real.hxx:186-200`) are `f64::MAX` / `-f64::MAX`.

use occt_core::cslib::normal::{normal_d1, CSLibDerivativeStatus};
use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;
use occt_core::precision::{epsilon, REAL_EPSILON};

use crate::surface::Surface;

pub use crate::geom_lprop::cl_props::LPropStatus;

/// `RealLast()` (`Standard_Real.hxx:186-189`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx:176-179`).
const REAL_FIRST: f64 = -f64::MAX;
/// `THE_DIVISION_FACTOR` (`GeomLProp_SurfaceUtils.hxx:280`).
const THE_DIVISION_FACTOR: f64 = 1.0e-3;
/// `THE_MIN_STEP` (`GeomLProp_SurfaceUtils.hxx:281`).
const THE_MIN_STEP: f64 = 1.0e-7;

/// `GeomLProp_SLPropsBase` (`GeomLProp_SLProps.hxx:37-315`).
#[derive(Clone, Debug)]
pub struct SlProps {
    u: f64,
    v: f64,
    /// `myDerOrder`.
    der_order: i32,
    /// `myCN`.
    cn: i32,
    /// `myLinTol`.
    lin_tol: f64,
    pnt: GpPnt,
    d1u: GpVec,
    d1v: GpVec,
    d2u: GpVec,
    d2v: GpVec,
    duv: GpVec,
    normal: GpDir,
    min_curv: f64,
    max_curv: f64,
    dir_min_curv: GpDir,
    dir_max_curv: GpDir,
    mean_curv: f64,
    gaus_curv: f64,
    sig_order_u: i32,
    sig_order_v: i32,
    u_tangent_status: LPropStatus,
    v_tangent_status: LPropStatus,
    normal_status: LPropStatus,
    curvature_status: LPropStatus,
}

/// `EvalSurfDerivatives` (`GeomLProp_SurfaceUtils.hxx:186-210`): a `switch`
/// on the order with no `default`, so `theOrder < 0 || > 2` leaves every
/// output untouched — the caller has already rejected that
/// (`hxx:70`, `:86`, `:103`).
fn eval_surf_derivatives(
    surf: &dyn Surface,
    u: f64,
    v: f64,
    order: i32,
    pnt: &mut GpPnt,
    d1u: &mut GpVec,
    d1v: &mut GpVec,
    d2u: &mut GpVec,
    d2v: &mut GpVec,
    duv: &mut GpVec,
) {
    match order {
        0 => {
            *pnt = surf.d0(u, v);
        }
        1 => {
            let (p, u1, v1) = surf.d1(u, v);
            *pnt = p;
            *d1u = u1;
            *d1v = v1;
        }
        2 => {
            let (p, u1, v1, u2, v2, uv) = surf.d2(u, v);
            *pnt = p;
            *d1u = u1;
            *d1v = v1;
            *d2u = u2;
            *d2v = v2;
            *duv = uv;
        }
        _ => {}
    }
}

/// `FindSurfTangentOrder` (`hxx:220-251`).
fn find_surf_tangent_order(
    d1: &GpVec,
    d2: &GpVec,
    cn: i32,
    tol_sq: f64,
    order: &mut i32,
    status: &mut LPropStatus,
) -> bool {
    let derivs = [d1, d2];
    *order = 0;
    while *order < 3 {
        *order += 1;
        if cn >= *order {
            if *order <= 2 && derivs[(*order - 1) as usize].square_magnitude() > tol_sq {
                *status = LPropStatus::Defined;
                return true;
            }
        } else {
            *status = LPropStatus::Undefined;
            return false;
        }
    }
    *status = LPropStatus::Undefined;
    false
}

/// `ComputeSurfTangent` (`hxx:263-343`).
#[allow(clippy::too_many_arguments)]
fn compute_surf_tangent(
    surf: &dyn Surface,
    u: f64,
    v: f64,
    first_deriv: &GpVec,
    sec_deriv: &GpVec,
    sig_order: i32,
    is_u: bool,
    dir: &mut GpDir,
) -> Result<(), String> {
    if sig_order == 1 {
        *dir = GpDir::from_vec(first_deriv)
            .map_err(|e| format!("LProp_NotDefined: {e}"))?;
        return Ok(());
    }

    let (u1, u2) = surf.u_range();
    let (v1, v2) = surf.v_range();

    let mut a_v = *sec_deriv;
    if is_u {
        let a_du = if u2 >= REAL_LAST || u1 <= REAL_FIRST {
            0.0
        } else {
            u2 - u1
        };
        let a_delta = (a_du * THE_DIVISION_FACTOR).max(THE_MIN_STEP);
        let other_u = if u - u1 < a_delta { u + a_delta } else { u - a_delta };
        let p1 = surf.d0(u.min(other_u), v);
        let p2 = surf.d0(u.max(other_u), v);
        let chord = GpVec::from_pnts(&p1, &p2);
        if a_v.dot(&chord) < 0.0 {
            a_v.reverse();
        }
        *dir = GpDir::from_vec(&a_v).map_err(|e| format!("LProp_NotDefined: {e}"))?;
    } else {
        let a_dv = if v2 >= REAL_LAST || v1 <= REAL_FIRST {
            0.0
        } else {
            v2 - v1
        };
        let a_delta = (a_dv * THE_DIVISION_FACTOR).max(THE_MIN_STEP);
        let other_v = if v - v1 < a_delta { v + a_delta } else { v - a_delta };
        let p1 = surf.d0(u, v.min(other_v));
        let p2 = surf.d0(u, v.max(other_v));
        let chord = GpVec::from_pnts(&p1, &p2);
        if a_v.dot(&chord) < 0.0 {
            a_v.reverse();
        }
        *dir = GpDir::from_vec(&a_v).map_err(|e| format!("LProp_NotDefined: {e}"))?;
    }
    Ok(())
}

/// `ComputeSurfNormal` (`hxx:346-358`): `CSLib::Normal(D1U, D1V, LinTol,
/// Status, Normal)` and `Status == CSLib_Done`.
fn compute_surf_normal(
    d1u: &GpVec,
    d1v: &GpVec,
    lin_tol: f64,
    normal: &mut GpDir,
) -> bool {
    let (status, n) = normal_d1(d1u, d1v, lin_tol);
    if let Some(n) = n {
        *normal = n;
    }
    status == CSLibDerivativeStatus::Done
}

/// `ComputeSurfCurvatures` (`hxx:363-481`).
#[allow(clippy::too_many_arguments)]
fn compute_surf_curvatures(
    d1u: &GpVec,
    d1v: &GpVec,
    d2u: &GpVec,
    d2v: &GpVec,
    duv: &GpVec,
    normal: &GpDir,
    min_curv: &mut f64,
    max_curv: &mut f64,
    dir_min: &mut GpDir,
    dir_max: &mut GpDir,
    mean_curv: &mut f64,
    gaus_curv: &mut f64,
) -> bool {
    let a_norm = GpVec::new(normal.x(), normal.y(), normal.z());

    let an_e = d1u.square_magnitude();
    let an_f = d1u.dot(d1v);
    let a_g = d1v.square_magnitude();

    let a_l = a_norm.dot(d2u);
    let a_m = a_norm.dot(duv);
    let a_n = a_norm.dot(d2v);

    let an_a0 = an_e * a_m - an_f * a_l;
    let a_b0 = an_e * a_n - a_g * a_l;
    let a_c0 = an_f * a_n - a_g * a_m;

    let a_max_abc = an_a0.abs().max(a_b0.abs()).max(a_c0.abs());
    if a_max_abc < REAL_EPSILON {
        // Umbilic point.
        if a_g < REAL_EPSILON {
            return false;
        }
        *min_curv = a_n / a_g;
        *max_curv = *min_curv;
        *dir_min = match GpDir::from_vec(d1u) {
            Ok(d) => d,
            Err(_) => return false,
        };
        let cross = d1u.crossed(&a_norm);
        *dir_max = match GpDir::from_vec(&cross) {
            Ok(d) => d,
            Err(_) => return false,
        };
        *mean_curv = *min_curv;
        *gaus_curv = *min_curv * *min_curv;
        return true;
    }

    let an_a = an_a0 / a_max_abc;
    let a_b = a_b0 / a_max_abc;
    let a_c = a_c0 / a_max_abc;

    let (a_curv1, a_curv2, a_vect1, a_vect2);
    if an_a.abs() > REAL_EPSILON {
        let a_root = DirectPolynomialRoots::new3(an_a, a_b, a_c);
        if a_root.nb_solutions() != 2 {
            return false;
        }
        let root1 = a_root.value(1);
        let root2 = a_root.value(2);
        a_curv1 = ((a_l * root1 + 2.0 * a_m) * root1 + a_n) / ((an_e * root1 + 2.0 * an_f) * root1 + a_g);
        a_curv2 = ((a_l * root2 + 2.0 * a_m) * root2 + a_n) / ((an_e * root2 + 2.0 * an_f) * root2 + a_g);
        a_vect1 = d1u.multiplied_scalar(root1).added(d1v);
        a_vect2 = d1u.multiplied_scalar(root2).added(d1v);
    } else if a_c.abs() > REAL_EPSILON {
        let a_root = DirectPolynomialRoots::new3(a_c, a_b, an_a);
        if a_root.nb_solutions() != 2 {
            return false;
        }
        let root1 = a_root.value(1);
        let root2 = a_root.value(2);
        a_curv1 = ((a_n * root1 + 2.0 * a_m) * root1 + a_l) / ((a_g * root1 + 2.0 * an_f) * root1 + an_e);
        a_curv2 = ((a_n * root2 + 2.0 * a_m) * root2 + a_l) / ((a_g * root2 + 2.0 * an_f) * root2 + an_e);
        a_vect1 = d1u.added(&d1v.multiplied_scalar(root1));
        a_vect2 = d1u.added(&d1v.multiplied_scalar(root2));
    } else {
        a_curv1 = a_l / an_e;
        a_curv2 = a_n / a_g;
        a_vect1 = *d1u;
        a_vect2 = *d1v;
    }

    if a_curv1 < a_curv2 {
        *min_curv = a_curv1;
        *max_curv = a_curv2;
        *dir_min = match GpDir::from_vec(&a_vect1) {
            Ok(d) => d,
            Err(_) => return false,
        };
        *dir_max = match GpDir::from_vec(&a_vect2) {
            Ok(d) => d,
            Err(_) => return false,
        };
    } else {
        *min_curv = a_curv2;
        *max_curv = a_curv1;
        *dir_min = match GpDir::from_vec(&a_vect2) {
            Ok(d) => d,
            Err(_) => return false,
        };
        *dir_max = match GpDir::from_vec(&a_vect1) {
            Ok(d) => d,
            Err(_) => return false,
        };
    }

    let an_eg_ff = (an_e * a_g) - (an_f * an_f);
    *mean_curv = ((a_n * an_e) - (2.0 * a_m * an_f) + (a_l * a_g)) / (2.0 * an_eg_ff);
    *gaus_curv = ((a_l * a_n) - (a_m * a_m)) / an_eg_ff;
    true
}

impl SlProps {
    /// `GeomLProp_SLPropsBase(S, U, V, N, Resolution)` (`hxx:54-72`).
    ///
    /// `Standard_OutOfRange_Raise_if(N < 0 || N > 2)` becomes the `N` range
    /// check of the constructor's caller (`Err`).
    pub fn new(
        surf: &dyn Surface,
        u: f64,
        v: f64,
        n: i32,
        resolution: f64,
    ) -> Result<Self, String> {
        if !(0..=2).contains(&n) {
            return Err("GeomLProp_SLPropsBase: N out of range".to_string());
        }
        let mut p = Self::blank(n, resolution);
        p.cn = 4;
        p.set_parameters(surf, u, v);
        Ok(p)
    }

    /// `GeomLProp_SLPropsBase(S, N, Resolution)` (`hxx:75-89`) — no
    /// parameters set, `myU = myV = RealLast()`.
    ///
    /// The port carries no `mySurf` (the surface is passed per call), so `S`
    /// has no stored effect; only `myCN = 4` (`hxx:82`) is observable.
    pub fn with_surface(_surf: &dyn Surface, n: i32, resolution: f64) -> Result<Self, String> {
        if !(0..=2).contains(&n) {
            return Err("GeomLProp_SLPropsBase: N out of range".to_string());
        }
        let mut p = Self::blank(n, resolution);
        p.cn = 4;
        Ok(p)
    }

    /// `GeomLProp_SLPropsBase(N, Resolution)` (`hxx:92-108`): no surface and
    /// no parameters, `myCN = 0`. `SetSurface` raises `myCN` back to `4`.
    pub fn empty(n: i32, resolution: f64) -> Result<Self, String> {
        if !(0..=2).contains(&n) {
            return Err("GeomLProp_SLPropsBase: N out of range (bad level)".to_string());
        }
        Ok(Self::blank(n, resolution))
    }

    fn blank(n: i32, resolution: f64) -> Self {
        Self {
            u: REAL_LAST,
            v: REAL_LAST,
            der_order: n,
            cn: 0,
            lin_tol: resolution,
            pnt: GpPnt::zero(),
            d1u: GpVec::zero(),
            d1v: GpVec::zero(),
            d2u: GpVec::zero(),
            d2v: GpVec::zero(),
            duv: GpVec::zero(),
            normal: GpDir::default_dir(),
            min_curv: 0.0,
            max_curv: 0.0,
            dir_min_curv: GpDir::default_dir(),
            dir_max_curv: GpDir::default_dir(),
            mean_curv: 0.0,
            gaus_curv: 0.0,
            sig_order_u: 0,
            sig_order_v: 0,
            u_tangent_status: LPropStatus::Undecided,
            v_tangent_status: LPropStatus::Undecided,
            normal_status: LPropStatus::Undecided,
            curvature_status: LPropStatus::Undecided,
        }
    }

    /// `SetSurface(S)` (`hxx:111-117`). The port carries no `mySurf` (the
    /// surface is passed per call, like `BRepLProp_SLProps` does), so only the
    /// `myCN = 4` reset (`hxx:115`) has an effect.
    pub fn set_surface(&mut self) {
        self.cn = 4;
    }

    /// `SetParameters(U, V)` (`hxx:120-137` → `hxx:488-511`).
    pub fn set_parameters(&mut self, surf: &dyn Surface, u: f64, v: f64) {
        self.u = u;
        self.v = v;
        eval_surf_derivatives(
            surf,
            u,
            v,
            self.der_order,
            &mut self.pnt,
            &mut self.d1u,
            &mut self.d1v,
            &mut self.d2u,
            &mut self.d2v,
            &mut self.duv,
        );
        self.u_tangent_status = LPropStatus::Undecided;
        self.v_tangent_status = LPropStatus::Undecided;
        self.normal_status = LPropStatus::Undecided;
        self.curvature_status = LPropStatus::Undecided;
    }

    /// `myU` / `myV` as last stored.
    pub fn parameters(&self) -> (f64, f64) {
        (self.u, self.v)
    }

    /// `Value()` (`hxx:141`).
    pub fn value(&self) -> GpPnt {
        self.pnt
    }

    /// `EnsureSurfDeriv(...)` (`hxx:513-534`) for `myD1u`.
    fn ensure_d1u(&mut self, surf: &dyn Surface) -> GpVec {
        if self.der_order < 1 {
            self.der_order = 1;
            self.reeval(surf);
        }
        self.d1u
    }

    fn ensure_d1v(&mut self, surf: &dyn Surface) -> GpVec {
        if self.der_order < 1 {
            self.der_order = 1;
            self.reeval(surf);
        }
        self.d1v
    }

    fn ensure_d2u(&mut self, surf: &dyn Surface) -> GpVec {
        if self.der_order < 2 {
            self.der_order = 2;
            self.reeval(surf);
        }
        self.d2u
    }

    fn ensure_d2v(&mut self, surf: &dyn Surface) -> GpVec {
        if self.der_order < 2 {
            self.der_order = 2;
            self.reeval(surf);
        }
        self.d2v
    }

    fn reeval(&mut self, surf: &dyn Surface) {
        let (u, v, order) = (self.u, self.v, self.der_order);
        eval_surf_derivatives(
            surf,
            u,
            v,
            order,
            &mut self.pnt,
            &mut self.d1u,
            &mut self.d1v,
            &mut self.d2u,
            &mut self.d2v,
            &mut self.duv,
        );
    }

    /// `D1U()` (`hxx:145-150`): `EnsureSurfDeriv(..., 1, myD1u)`, i.e. the
    /// accessor raises `myDerOrder` to 1 and re-evaluates before handing the
    /// derivative out. No definiteness test - that is `TangentU`.
    pub fn d1u(&mut self, surf: &dyn Surface) -> GpVec {
        self.ensure_d1u(surf)
    }

    /// `D1V()` (`hxx:153-158`).
    pub fn d1v(&mut self, surf: &dyn Surface) -> GpVec {
        self.ensure_d1v(surf)
    }

    /// `D2U()` (`hxx:161-166`).
    pub fn d2u(&mut self, surf: &dyn Surface) -> GpVec {
        self.ensure_d2u(surf)
    }

    /// `D2V()` (`hxx:169-174`).
    pub fn d2v(&mut self, surf: &dyn Surface) -> GpVec {
        self.ensure_d2v(surf)
    }

    /// `IsTangentUDefined()` (`hxx:188-195` → `hxx:536-553`).
    pub fn is_tangent_u_defined(&mut self, surf: &dyn Surface) -> bool {
        if self.u_tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.u_tangent_status == LPropStatus::Defined {
            return true;
        }
        let d1 = self.ensure_d1u(surf);
        let d2 = self.ensure_d2u(surf);
        let (cn, tol) = (self.cn, self.lin_tol);
        let mut order = self.sig_order_u;
        let mut status = self.u_tangent_status;
        let ok = find_surf_tangent_order(&d1, &d2, cn, tol * tol, &mut order, &mut status);
        self.sig_order_u = order;
        self.u_tangent_status = status;
        ok
    }

    /// `IsTangentVDefined()` (`hxx:211-218` → `hxx:554-572`).
    pub fn is_tangent_v_defined(&mut self, surf: &dyn Surface) -> bool {
        if self.v_tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.v_tangent_status == LPropStatus::Defined {
            return true;
        }
        let d1 = self.ensure_d1v(surf);
        let d2 = self.ensure_d2v(surf);
        let (cn, tol) = (self.cn, self.lin_tol);
        let mut order = self.sig_order_v;
        let mut status = self.v_tangent_status;
        let ok = find_surf_tangent_order(&d1, &d2, cn, tol * tol, &mut order, &mut status);
        self.sig_order_v = order;
        self.v_tangent_status = status;
        ok
    }

    /// `TangentU(D)` (`hxx:197-206` → `hxx:575-591`).
    pub fn tangent_u(&mut self, surf: &dyn Surface) -> Result<GpDir, String> {
        if !self.is_tangent_u_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        let d1u = self.ensure_d1u(surf);
        let d2u = self.ensure_d2u(surf);
        let (u, v, sig) = (self.u, self.v, self.sig_order_u);
        let mut dir = GpDir::default_dir();
        compute_surf_tangent(surf, u, v, &d1u, &d2u, sig, true, &mut dir)?;
        Ok(dir)
    }

    /// `TangentV(D)` (`hxx:221-230` → `hxx:593-607`).
    pub fn tangent_v(&mut self, surf: &dyn Surface) -> Result<GpDir, String> {
        if !self.is_tangent_v_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        let d1v = self.ensure_d1v(surf);
        let d2v = self.ensure_d2v(surf);
        let (u, v, sig) = (self.u, self.v, self.sig_order_v);
        let mut dir = GpDir::default_dir();
        compute_surf_tangent(surf, u, v, &d1v, &d2v, sig, false, &mut dir)?;
        Ok(dir)
    }

    /// `IsNormalDefined()` (`hxx:233-237` → `hxx:609-627`).
    pub fn is_normal_defined(&mut self, surf: &dyn Surface) -> bool {
        if self.normal_status == LPropStatus::Undefined {
            return false;
        }
        if self.normal_status == LPropStatus::Defined {
            return true;
        }
        let d1u = self.ensure_d1u(surf);
        let d1v = self.ensure_d1v(surf);
        if compute_surf_normal(&d1u, &d1v, self.lin_tol, &mut self.normal) {
            self.normal_status = LPropStatus::Computed;
            return true;
        }
        self.normal_status = LPropStatus::Undefined;
        false
    }

    /// `Normal()` (`hxx:239-240` → `hxx:630-636`): throws
    /// `LProp_NotDefined` when the normal is not defined.
    pub fn normal(&mut self, surf: &dyn Surface) -> Result<GpDir, String> {
        if !self.is_normal_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok(self.normal)
    }

    /// `IsCurvatureDefined()` (`hxx:242-263` → `hxx:638-697`).
    pub fn is_curvature_defined(&mut self, surf: &dyn Surface) -> bool {
        if self.curvature_status == LPropStatus::Undefined {
            return false;
        }
        if self.curvature_status == LPropStatus::Defined {
            return true;
        }
        if self.cn < 2 {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if !self.is_normal_defined(surf) {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if !self.is_tangent_u_defined(surf) || !self.is_tangent_v_defined(surf) {
            self.curvature_status = LPropStatus::Undefined;
            return false;
        }
        if self.der_order < 2 {
            // `theProps.D2U()` — the call is made for its side effect of
            // raising `myDerOrder` and re-evaluating the derivatives.
            let _ = self.ensure_d2u(surf);
        }
        let (d1u, d1v) = (self.d1u, self.d1v);
        let (d2u, d2v, duv) = (self.d2u, self.d2v, self.duv);
        let normal = self.normal;
        let (mut min_curv, mut max_curv) = (self.min_curv, self.max_curv);
        let (mut dir_min, mut dir_max) = (self.dir_min_curv, self.dir_max_curv);
        let (mut mean_curv, mut gaus_curv) = (self.mean_curv, self.gaus_curv);
        let ok = compute_surf_curvatures(
            &d1u,
            &d1v,
            &d2u,
            &d2v,
            &duv,
            &normal,
            &mut min_curv,
            &mut max_curv,
            &mut dir_min,
            &mut dir_max,
            &mut mean_curv,
            &mut gaus_curv,
        );
        self.min_curv = min_curv;
        self.max_curv = max_curv;
        self.dir_min_curv = dir_min;
        self.dir_max_curv = dir_max;
        self.mean_curv = mean_curv;
        self.gaus_curv = gaus_curv;
        if ok {
            self.curvature_status = LPropStatus::Computed;
            return true;
        }
        self.curvature_status = LPropStatus::Undefined;
        false
    }

    /// `IsUmbilic()` (`hxx:266-269` → `hxx:706-710`): `|max - min| <
    /// |Epsilon(max)|`; throws `LProp_NotDefined` through
    /// `IsCurvatureDefined`.
    pub fn is_umbilic(&mut self, surf: &dyn Surface) -> Result<bool, String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok((self.max_curv - self.min_curv).abs() < epsilon(self.max_curv).abs())
    }

    /// `MaxCurvature()` (`hxx:271-272` → `hxx:699-704`).
    pub fn max_curvature(&mut self, surf: &dyn Surface) -> Result<f64, String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok(self.max_curv)
    }

    /// `MinCurvature()` (`hxx:274-275`).
    pub fn min_curvature(&mut self, surf: &dyn Surface) -> Result<f64, String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok(self.min_curv)
    }

    /// `MeanCurvature()` (`hxx:282-283`).
    pub fn mean_curvature(&mut self, surf: &dyn Surface) -> Result<f64, String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok(self.mean_curv)
    }

    /// `GaussianCurvature()` (`hxx:285-286`).
    pub fn gaussian_curvature(&mut self, surf: &dyn Surface) -> Result<f64, String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok(self.gaus_curv)
    }

    /// `CurvatureDirections(MaxD, MinD)` (`hxx:277-280` → `hxx:712-728`).
    pub fn curvature_directions(
        &mut self,
        surf: &dyn Surface,
    ) -> Result<(GpDir, GpDir), String> {
        if !self.is_curvature_defined(surf) {
            return Err("LProp_NotDefined".to_string());
        }
        Ok((self.dir_max_curv, self.dir_min_curv))
    }
}
