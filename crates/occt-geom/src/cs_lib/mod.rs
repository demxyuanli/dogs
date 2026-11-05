//! `CSLib` — surface normal computation and its derivatives.
//!
//! Source: `FoundationClasses/TKMath/CSLib` (OCCT 8.0.0):
//!
//! | OCCT | here |
//! |---|---|
//! | `CSLib.hxx:31-182` | this file |
//! | `CSLib.cxx:46-179` (four `Normal` overloads) | [`normal`], [`normal_taylor`], [`normal_mag`], [`normal_higher`] |
//! | `CSLib.cxx:391-565` (`DNNUV` / `DNNormal`) | [`dnnuv`], [`dnnuv_two`], [`dnnorm`] |
//! | `CSLib_NormalPolyDef.cxx:24-111` | [`normal_poly_def`] |
//! | `CSLib_Class2d.cxx:45-342` | [`class2d`] |
//!
//! Translation notes:
//! - `NCollection_Array2<gp_Vec>` becomes [`Array2Vec`], a dense row-major array
//!   that keeps OCCT's `Lower`/`Upper` bounds on both directions; `Value` panics
//!   on an index OCCT would answer with `Standard_OutOfRange`,
//! - `gp::Resolution()` is `RealSmall()` = `DBL_MIN`
//!   (`occt_core::precision::REAL_SMALL`), see the note on
//!   `occt_core::precision::RESOLUTION`,
//! - a `gp_Dir` whose `gp_Vec` construction would raise
//!   `Standard_ConstructionError` is reported as `None` in the status-carrying
//!   overloads; the two overloads that have no status to report it with
//!   (`normal_taylor`'s `theNormal` and `normal_higher`'s) cannot reach it, as
//!   the source guards every construction with a magnitude test.

pub mod class2d;
pub mod normal_poly_def;

use occt_core::bspl::plib;
use occt_core::gp::{GpDir, GpVec};
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::math_function_roots::FunctionRoots;
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_EPSILON, REAL_SMALL};

pub use class2d::{Class2d, Class2dResult};
pub use normal_poly_def::NormalPolyDef;

/// `CSLib_DerivativeStatus` (`CSLib_DerivativeStatus.hxx:20-35`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivativeStatus {
    /// `CSLib_Done`.
    Done,
    /// `CSLib_D1uIsNull`.
    D1uIsNull,
    /// `CSLib_D1vIsNull`.
    D1vIsNull,
    /// `CSLib_D1IsNull`.
    D1IsNull,
    /// `CSLib_D1uIsParallelD1v`.
    D1uIsParallelD1v,
}

/// `CSLib_NormalStatus` (`CSLib_NormalStatus.hxx:20-37`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalStatus {
    /// `CSLib_Singular`.
    Singular,
    /// `CSLib_Defined`.
    Defined,
    /// `CSLib_InfinityOfSolutions`.
    InfinityOfSolutions,
    /// `CSLib_D1NuIsNull`.
    D1NuIsNull,
    /// `CSLib_D1NvIsNull`.
    D1NvIsNull,
    /// `CSLib_D1NIsNull`.
    D1NIsNull,
    /// `CSLib_D1NuNvRatioIsNull`.
    D1NuNvRatioIsNull,
    /// `CSLib_D1NvNuRatioIsNull`.
    D1NvNuRatioIsNull,
    /// `CSLib_D1NuIsParallelD1Nv`.
    D1NuIsParallelD1Nv,
}

/// Angular tolerance of `CSLib.cxx:34`.
const THE_PARALLEL_ANGULAR_TOL: f64 = 1e-6;
/// `CSLib.cxx:37`.
const THE_MAX_ROOT_ITERATIONS: i32 = 200;
/// `CSLib.cxx:40`.
const THE_ROOT_FINDING_TOL: f64 = 1e-5;

/// `NCollection_Array2<gp_Vec>` with OCCT's two-direction bounds.
#[derive(Debug, Clone)]
pub struct Array2Vec {
    lower_i: i32,
    upper_i: i32,
    lower_j: i32,
    upper_j: i32,
    data: Vec<GpVec>,
}

impl Array2Vec {
    /// Allocates the `[lower_i, upper_i] x [lower_j, upper_j]` box with every
    /// entry zero, the way `NCollection_Array2<gp_Vec>` default-constructs its
    /// elements.
    pub fn new(lower_i: i32, upper_i: i32, lower_j: i32, upper_j: i32) -> Self {
        let rows = (upper_i - lower_i + 1).max(0) as usize;
        let cols = (upper_j - lower_j + 1).max(0) as usize;
        Self {
            lower_i,
            upper_i,
            lower_j,
            upper_j,
            data: vec![GpVec::zero(); rows * cols],
        }
    }

    /// `Value(I, J)` (`NCollection_Array2.hxx`): panics on an out-of-range
    /// index, which OCCT reports as `Standard_OutOfRange`.
    pub fn value(&self, i: i32, j: i32) -> GpVec {
        assert!(
            i >= self.lower_i && i <= self.upper_i && j >= self.lower_j && j <= self.upper_j,
            "Array2Vec::value({i}, {j}) outside [{}, {}] x [{}, {}]",
            self.lower_i,
            self.upper_i,
            self.lower_j,
            self.upper_j
        );
        let cols = (self.upper_j - self.lower_j + 1) as usize;
        self.data[(i - self.lower_i) as usize * cols + (j - self.lower_j) as usize].clone()
    }

    /// `ChangeValue(I, J) = V`.
    pub fn set_value(&mut self, i: i32, j: i32, v: GpVec) {
        assert!(
            i >= self.lower_i && i <= self.upper_i && j >= self.lower_j && j <= self.upper_j,
            "Array2Vec::set_value({i}, {j}) outside [{}, {}] x [{}, {}]",
            self.lower_i,
            self.upper_i,
            self.lower_j,
            self.upper_j
        );
        let cols = (self.upper_j - self.lower_j + 1) as usize;
        self.data[(i - self.lower_i) as usize * cols + (j - self.lower_j) as usize] = v;
    }

    /// `LowerRow()`.
    pub fn lower_row(&self) -> i32 {
        self.lower_i
    }

    /// `UpperRow()`.
    pub fn upper_row(&self) -> i32 {
        self.upper_i
    }
}

/// `CSLib::Normal(theD1U, theD1V, theSinTol, theStatus, theNormal)`
/// (`CSLib.cxx:46-79`).
pub fn normal(d1u: &GpVec, d1v: &GpVec, sin_tol: f64) -> (DerivativeStatus, Option<GpDir>) {
    let d1u_mag = d1u.square_magnitude();
    let d1v_mag = d1v.square_magnitude();
    let d1u_x_d1v = d1u.crossed(d1v);

    if d1u_mag <= REAL_SMALL && d1v_mag <= REAL_SMALL {
        (DerivativeStatus::D1IsNull, None)
    } else if d1u_mag <= REAL_SMALL {
        (DerivativeStatus::D1uIsNull, None)
    } else if d1v_mag <= REAL_SMALL {
        (DerivativeStatus::D1vIsNull, None)
    } else {
        let sin2 = d1u_x_d1v.square_magnitude() / (d1u_mag * d1v_mag);
        if sin2 < sin_tol * sin_tol {
            (DerivativeStatus::D1uIsParallelD1v, None)
        } else {
            (
                DerivativeStatus::Done,
                Some(GpDir::from_vec(&d1u_x_d1v).expect("CSLib::Normal: null cross product")),
            )
        }
    }
}

/// `CSLib::Normal(theD1U, theD1V, theD2U, theD2V, theD2UV, theSinTol,
/// theDone, theStatus, theNormal)` (`CSLib.cxx:83-151`).
pub fn normal_taylor(
    d1u: &GpVec,
    d1v: &GpVec,
    d2u: &GpVec,
    d2v: &GpVec,
    d2uv: &GpVec,
    sin_tol: f64,
) -> (bool, NormalStatus, Option<GpDir>) {
    // `aD1Nu = aD2U ^ aD1V + aD1U ^ aD2UV` / `aD1Nv = aD2UV ^ aD1V + aD1U ^ aD2V`
    // (`cxx:99-104`).
    let d1nu = d2u.crossed(d1v).added(&d1u.crossed(d2uv));
    let d1nv = d2uv.crossed(d1v).added(&d1u.crossed(d2v));

    let l_d1nu = d1nu.square_magnitude();
    let l_d1nv = d1nv.square_magnitude();

    if l_d1nu <= REAL_EPSILON && l_d1nv <= REAL_EPSILON {
        (false, NormalStatus::D1NIsNull, None)
    } else if l_d1nu < REAL_EPSILON {
        (
            true,
            NormalStatus::D1NuIsNull,
            Some(GpDir::from_vec(&d1nv).expect("CSLib::Normal: null d1nv")),
        )
    } else if l_d1nv < REAL_EPSILON {
        (
            true,
            NormalStatus::D1NvIsNull,
            Some(GpDir::from_vec(&d1nu).expect("CSLib::Normal: null d1nu")),
        )
    } else if (l_d1nv / l_d1nu) <= REAL_EPSILON {
        (false, NormalStatus::D1NvNuRatioIsNull, None)
    } else if (l_d1nu / l_d1nv) <= REAL_EPSILON {
        (false, NormalStatus::D1NuNvRatioIsNull, None)
    } else {
        let d1n_cross = d1nu.crossed(&d1nv);
        let sin2 = d1n_cross.square_magnitude() / (l_d1nu * l_d1nv);
        if sin2 < sin_tol * sin_tol {
            (
                true,
                NormalStatus::D1NuIsParallelD1Nv,
                Some(GpDir::from_vec(&d1nu).expect("CSLib::Normal: null d1nu")),
            )
        } else {
            (false, NormalStatus::InfinityOfSolutions, None)
        }
    }
}

/// `CSLib::Normal(theD1U, theD1V, theMagTol, theStatus, theNormal)`
/// (`CSLib.cxx:155-179`).
pub fn normal_mag(d1u: &GpVec, d1v: &GpVec, mag_tol: f64) -> (NormalStatus, Option<GpDir>) {
    let d1u_mag = d1u.magnitude();
    let d1v_mag = d1v.magnitude();
    let d1u_x_d1v = d1u.crossed(d1v);
    let n_mag = d1u_x_d1v.magnitude();

    if n_mag <= mag_tol || d1u_mag <= mag_tol || d1v_mag <= mag_tol {
        (NormalStatus::Singular, None)
    } else {
        // `const gp_Dir aD1UDir(theD1U); const gp_Dir aD1VDir(theD1V);
        //  theNormal = gp_Dir(aD1UDir.Crossed(aD1VDir));` (`cxx:174-176`).
        let d1u_dir = GpDir::from_vec(d1u).expect("CSLib::Normal: null d1u");
        let d1v_dir = GpDir::from_vec(d1v).expect("CSLib::Normal: null d1v");
        match d1u_dir.crossed(&d1v_dir) {
            Ok(n) => (NormalStatus::Defined, Some(n)),
            Err(_) => (NormalStatus::Defined, None),
        }
    }
}

/// `CSLib::Normal(theMaxOrder, theDerNUV, theSinTol, theU, theV, theUmin,
/// theUmax, theVmin, theVmax, theStatus, theNormal, theOrderU, theOrderV)`
/// (`CSLib.cxx:183-387`).
///
/// Answers `(status, normal, order_u, order_v)`; `order_u` / `order_v` are
/// written on every path, exactly like the out-parameters of the source.
pub fn normal_higher(
    max_order: i32,
    der_nuv: &Array2Vec,
    sin_tol: f64,
    u: f64,
    v: f64,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> (NormalStatus, Option<GpDir>, i32, i32) {
    let mut an_order = -1;
    let mut a_found_u_idx = 0;
    let mut a_found = false;
    let mut a_d = GpVec::zero();

    // `while (!aFound && anOrder < theMaxOrder) { ... }` (`cxx:204-216`).
    while !a_found && an_order < max_order {
        an_order += 1;
        a_found_u_idx = an_order;
        while a_found_u_idx >= 0 && !a_found {
            let a_v_idx = an_order - a_found_u_idx;
            a_d = der_nuv.value(a_found_u_idx, a_v_idx);
            a_found = a_d.magnitude() >= sin_tol;
            a_found_u_idx -= 1;
        }
    }

    let the_order_u = a_found_u_idx + 1;
    let the_order_v = an_order - the_order_u;

    if !a_found {
        return (NormalStatus::Singular, None, the_order_u, the_order_v);
    }

    if an_order == 0 {
        return (
            NormalStatus::Defined,
            Some(GpDir::from_vec(&a_d.normalized()).expect("CSLib::Normal: null normal")),
            the_order_u,
            the_order_v,
        );
    }

    let a_vk0 = der_nuv.value(the_order_u, the_order_v);
    let mut a_ratio = vec![0.0f64; (an_order + 1) as usize];

    // `cxx:240-262`.
    let mut a_ratio_idx = 0;
    let mut is_defined = false;
    while a_ratio_idx <= an_order && !is_defined {
        let a_der_vec = der_nuv.value(a_ratio_idx, an_order - a_ratio_idx);
        if a_der_vec.magnitude() <= sin_tol {
            a_ratio[a_ratio_idx as usize] = 0.0;
        } else if a_der_vec.is_parallel_ang(&a_vk0, THE_PARALLEL_ANGULAR_TOL) {
            let mut magnitude_ratio = a_der_vec.magnitude() / a_vk0.magnitude();
            if a_der_vec.is_opposite_ang(&a_vk0, THE_PARALLEL_ANGULAR_TOL) {
                magnitude_ratio = -magnitude_ratio;
            }
            a_ratio[a_ratio_idx as usize] = magnitude_ratio;
        } else {
            is_defined = true;
        }
        a_ratio_idx += 1;
    }

    if is_defined {
        return (
            NormalStatus::Defined,
            Some(GpDir::from_vec(&a_d.normalized()).expect("CSLib::Normal: null normal")),
            the_order_u,
            the_order_v,
        );
    }

    // `cxx:271-315`: the angular domain of the corner / edge / interior case.
    let mut a_inf = -std::f64::consts::PI;
    let mut a_sup = std::f64::consts::PI;

    let is_fu = (u - umin).abs() < PCONFUSION;
    let is_lu = (u - umax).abs() < PCONFUSION;
    let is_fv = (v - vmin).abs() < PCONFUSION;
    let is_lv = (v - vmax).abs() < PCONFUSION;

    if is_lu {
        a_inf = std::f64::consts::FRAC_PI_2;
        a_sup = 3.0 * std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_inf = std::f64::consts::PI;
        }
        if is_fv {
            a_sup = std::f64::consts::PI;
        }
    } else if is_fu {
        a_sup = std::f64::consts::FRAC_PI_2;
        a_inf = -std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_sup = 0.0;
        }
        if is_fv {
            a_inf = 0.0;
        }
    } else if is_lv {
        a_inf = -std::f64::consts::PI;
        a_sup = 0.0;
    } else if is_fv {
        a_inf = 0.0;
        a_sup = std::f64::consts::PI;
    }

    let mut a_changes_sign = false;
    let mut a_vsuiv = 0.0;

    // `CSLib_NormalPolyDef aPoly(anOrder, aRatio); math_FunctionRoots aFindRoots(
    //  aPoly, aInf, aSup, 200, 1.e-5, Precision::Confusion(),
    //  Precision::Confusion());` (`cxx:324-331`).
    let mut a_poly = NormalPolyDef::new(an_order, &a_ratio);
    let a_find_roots = FunctionRoots::new(
        &mut a_poly,
        a_inf,
        a_sup,
        THE_MAX_ROOT_ITERATIONS,
        THE_ROOT_FINDING_TOL,
        CONFUSION,
        CONFUSION,
        0.0,
    );

    if a_find_roots.is_done() && a_find_roots.nb_solutions() > 0 {
        let a_nb_sol = a_find_roots.nb_solutions();
        // `NCollection_Array1<double> aSol(0, aNbSol + 1);` plus the sort of
        // `[1, aNbSol]` (`cxx:337-342`).
        let mut a_sol = vec![0.0f64; (a_nb_sol + 2) as usize];
        for a_root_idx in 1..=a_nb_sol {
            a_sol[a_root_idx as usize] = a_find_roots.value(a_root_idx);
        }
        a_sol[1..=(a_nb_sol as usize)].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        a_sol[0] = a_inf;
        a_sol[(a_nb_sol + 1) as usize] = a_sup;

        // `cxx:349-366`.
        let mut a_first = 0;
        let mut a_vprec = 0.0;
        for a_interval_idx in 0..=a_nb_sol {
            let lo = a_sol[a_interval_idx as usize];
            let hi = a_sol[(a_interval_idx + 1) as usize];
            if (hi - lo).abs() > PCONFUSION {
                let mut f = 0.0;
                let mut d = 0.0;
                a_poly.values((lo + hi) / 2.0, &mut f, &mut d);
                a_vsuiv = f;
                if a_first == 0 {
                    a_first = a_interval_idx;
                    a_changes_sign = false;
                    a_vprec = a_vsuiv;
                } else {
                    a_changes_sign = a_changes_sign || (a_vprec * a_vsuiv) < 0.0;
                    a_vprec = a_vsuiv;
                }
            }
        }
    } else {
        // `aChangesSign = false; aPoly.Value(aInf, aVsuiv);` (`cxx:368-374`).
        a_changes_sign = false;
        let mut f = 0.0;
        a_poly.value(a_inf, &mut f);
        a_vsuiv = f;
    }

    if a_changes_sign {
        (
            NormalStatus::InfinityOfSolutions,
            None,
            the_order_u,
            the_order_v,
        )
    } else {
        let sign = if a_vsuiv > 0.0 { 1.0 } else { -1.0 };
        (
            NormalStatus::Defined,
            Some(
                GpDir::from_vec(&a_vk0.normalized().multiplied_scalar(sign))
                    .expect("CSLib::Normal: null vk0"),
            ),
            the_order_u,
            the_order_v,
        )
    }
}

/// `CSLib::DNNUV(theNu, theNv, theDerSurf)` (`CSLib.cxx:391-408`).
pub fn dnnuv(nu: i32, nv: i32, der_surf: &Array2Vec) -> GpVec {
    let mut result = GpVec::zero();
    for i in 0..=nu {
        for j in 0..=nv {
            let a_vg = der_surf.value(i + 1, j);
            let a_vd = der_surf.value(nu - i, nv + 1 - j);
            let a_bin_coef = plib::bin(nu, i) * plib::bin(nv, j);
            result = result.added(&a_vg.crossed(&a_vd).multiplied_scalar(a_bin_coef));
        }
    }
    result
}

/// `CSLib::DNNUV(theNu, theNv, theDerSurf1, theDerSurf2)` (`CSLib.cxx:411-433`).
pub fn dnnuv_two(nu: i32, nv: i32, der_surf1: &Array2Vec, der_surf2: &Array2Vec) -> GpVec {
    let mut result = GpVec::zero();
    for i in 0..=nu {
        for j in 0..=nv {
            let a_vg = der_surf1.value(i + 1, j);
            let a_vd = der_surf2.value(nu - i, nv + 1 - j);
            let a_bin_coef = plib::bin(nu, i) * plib::bin(nv, j);
            result = result.added(&a_vg.crossed(&a_vd).multiplied_scalar(a_bin_coef));
        }
    }
    result
}

/// `CSLib::DNNormal(theNu, theNv, theDerNUV, theIduref, theIdvref)`
/// (`CSLib.cxx:436-565`).
pub fn dnnorm(nu: i32, nv: i32, der_nuv: &Array2Vec, iduref: i32, idvref: i32) -> GpVec {
    let a_kderiv = nu + nv;

    let mut a_der_vec_nor = Array2Vec::new(0, a_kderiv, 0, a_kderiv);
    let mut a_tab_scal = vec![0.0f64; ((a_kderiv + 1) * (a_kderiv + 1)) as usize];
    let mut a_tab_norm = vec![0.0f64; ((a_kderiv + 1) * (a_kderiv + 1)) as usize];
    let stride = (a_kderiv + 1) as usize;
    let at = |i: i32, j: i32| -> usize { i as usize * stride + j as usize };

    let base = der_nuv.value(iduref, idvref);
    let mut a_der_nor = base.normalized();
    a_der_vec_nor.set_value(0, 0, a_der_nor.clone());

    let a_dnorm0 = base.dot(&a_der_nor);
    a_tab_norm[at(0, 0)] = a_dnorm0;
    a_tab_scal[at(0, 0)] = 0.0;

    for a_mderiv in 1..=a_kderiv {
        for a_pderiv in 0..=a_mderiv {
            let a_qderiv = a_mderiv - a_pderiv;
            if a_pderiv > nu || a_qderiv > nv {
                continue;
            }

            // `cxx:467-516`: `n . d^(p,q) n`.
            let mut a_scal = 0.0;
            if a_pderiv > a_qderiv {
                for a_jderiv in 1..=a_qderiv {
                    a_scal -= plib::bin(a_qderiv, a_jderiv)
                        * a_der_vec_nor
                            .value(0, a_jderiv)
                            .dot(&a_der_vec_nor.value(a_pderiv, a_qderiv - a_jderiv));
                }
                for a_jderiv in 0..a_qderiv {
                    a_scal -= plib::bin(a_qderiv, a_jderiv)
                        * a_der_vec_nor
                            .value(a_pderiv, a_jderiv)
                            .dot(&a_der_vec_nor.value(0, a_qderiv - a_jderiv));
                }
                for a_ideriv in 1..a_pderiv {
                    for a_jderiv in 0..=a_qderiv {
                        a_scal -= plib::bin(a_pderiv, a_ideriv) * plib::bin(a_qderiv, a_jderiv)
                            * a_der_vec_nor.value(a_ideriv, a_jderiv).dot(
                                &a_der_vec_nor
                                    .value(a_pderiv - a_ideriv, a_qderiv - a_jderiv),
                            );
                    }
                }
            } else {
                for a_ideriv in 1..=a_pderiv {
                    a_scal -= plib::bin(a_pderiv, a_ideriv)
                        * a_der_vec_nor
                            .value(a_ideriv, 0)
                            .dot(&a_der_vec_nor.value(a_pderiv - a_ideriv, a_qderiv));
                }
                for a_ideriv in 0..a_pderiv {
                    a_scal -= plib::bin(a_pderiv, a_ideriv)
                        * a_der_vec_nor
                            .value(a_ideriv, a_qderiv)
                            .dot(&a_der_vec_nor.value(a_pderiv - a_ideriv, 0));
                }
                for a_ideriv in 0..=a_pderiv {
                    for a_jderiv in 1..a_qderiv {
                        a_scal -= plib::bin(a_pderiv, a_ideriv) * plib::bin(a_qderiv, a_jderiv)
                            * a_der_vec_nor.value(a_ideriv, a_jderiv).dot(
                                &a_der_vec_nor
                                    .value(a_pderiv - a_ideriv, a_qderiv - a_jderiv),
                            );
                    }
                }
            }
            a_tab_scal[at(a_pderiv, a_qderiv)] = a_scal / 2.0;

            // `cxx:518-537`: the derivative `(p,q)` of `|NUV|`.
            let mut a_dnorm = der_nuv
                .value(a_pderiv + iduref, a_qderiv + idvref)
                .dot(&a_der_vec_nor.value(0, 0));

            for a_jderiv in 0..a_qderiv {
                a_dnorm -= plib::bin(a_qderiv + idvref, a_jderiv + idvref)
                    * a_tab_norm[at(a_pderiv, a_jderiv)]
                    * a_tab_scal[at(0, a_qderiv - a_jderiv)];
            }
            for a_ideriv in 0..a_pderiv {
                for a_jderiv in 0..=a_qderiv {
                    a_dnorm -= plib::bin(a_pderiv + iduref, a_ideriv + iduref)
                        * plib::bin(a_qderiv + idvref, a_jderiv + idvref)
                        * a_tab_norm[at(a_ideriv, a_jderiv)]
                        * a_tab_scal[at(a_pderiv - a_ideriv, a_qderiv - a_jderiv)];
                }
            }
            a_tab_norm[at(a_pderiv, a_qderiv)] = a_dnorm;

            // `cxx:539-563`: the derivative `(p,q)` of `n`.
            a_der_nor = der_nuv.value(a_pderiv + iduref, a_qderiv + idvref);

            for a_jderiv in 1..=a_qderiv {
                a_der_nor = a_der_nor.subtracted(
                    &a_der_vec_nor
                        .value(a_pderiv, a_qderiv - a_jderiv)
                        .multiplied_scalar(
                            plib::bin(a_pderiv + iduref, iduref)
                                * plib::bin(a_qderiv + idvref, a_jderiv + idvref)
                                * a_tab_norm[at(0, a_jderiv)],
                        ),
                );
            }
            for a_ideriv in 1..=a_pderiv {
                for a_jderiv in 0..=a_qderiv {
                    a_der_nor = a_der_nor.subtracted(
                        &a_der_vec_nor
                            .value(a_pderiv - a_ideriv, a_qderiv - a_jderiv)
                            .multiplied_scalar(
                                plib::bin(a_pderiv + iduref, a_ideriv + iduref)
                                    * plib::bin(a_qderiv + idvref, a_jderiv + idvref)
                                    * a_tab_norm[at(a_ideriv, a_jderiv)],
                            ),
                    );
                }
            }

            a_der_nor = a_der_nor.divided(
                plib::bin(a_pderiv + iduref, iduref)
                    * plib::bin(a_qderiv + idvref, idvref)
                    * a_tab_norm[at(0, 0)],
            );
            a_der_vec_nor.set_value(a_pderiv, a_qderiv, a_der_nor.clone());
        }
    }

    a_der_vec_nor.value(nu, nv)
}
