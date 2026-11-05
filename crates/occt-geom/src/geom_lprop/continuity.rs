//! `GeomLProp::Continuity` — the continuity of two curves meeting at a
//! parameter pair. Source: `ModelingData/TKGeomBase/GeomLProp/GeomLProp.cxx:27-213`
//! (OCCT 8.0.0).
//!
//! The order actually compared follows `GeomLProp.cxx:123-146`: a B-spline
//! interior knot downgrades `cont1`/`cont2` from `Geom_Curve::Continuity()` to
//! `Degree - Multiplicity`, and `GeomAbs_CN`/`GeomAbs_C3` ask for order 3,
//! `GeomAbs_C2` for order 2, `GeomAbs_C1` for order 1, everything else 0.
//!
//! Faithfulness notes:
//! * `GeomLProp_CLProps` is built on the **original** curves (`cxx:147-148`),
//!   not on the untrimmed bases the knot lookup used.
//! * `throw Standard_Failure("Courbes non jointives")` (`cxx:151`) becomes
//!   `Err`, the failure channel of this crate.
//! * `GeomAbsToInteger` (`cxx:27-57`) is the identity on the
//!   [`Shape`] discriminants, so the `u8` returned by [`Curve::continuity`] is
//!   already the OCCT integer.

use crate::curve::Curve;
use crate::geom_lprop::cl_props::ClProps;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION};

/// `GeomLProp::Continuity(C1, C2, u1, u2, r1, r2, tl, ta)`
/// (`GeomLProp.cxx:61-203`).
///
/// `tl` is the linear tolerance and `ta` the angular one; `r1`/`r2` reverse the
/// first derivative (and the tangent) of the corresponding curve before the
/// comparison, which is how a closed junction is tested against itself.
pub fn continuity(
    c1: &dyn Curve,
    c2: &dyn Curve,
    u1: f64,
    u2: f64,
    r1: bool,
    r2: bool,
    tl: f64,
    ta: f64,
) -> Result<Shape, String> {
    let mut cont = Shape::C0;
    let mut cont1 = i32::from(c1.continuity());
    let mut cont2 = i32::from(c2.continuity());

    // `cxx:80-92`: one level of `Geom_TrimmedCurve::BasisCurve()`.
    let basis1 = c1.trim_basis();
    let basis2 = c2.trim_basis();
    let a_curve1: &dyn Curve = basis1.as_deref().unwrap_or(c1);
    let a_curve2: &dyn Curve = basis2.as_deref().unwrap_or(c2);

    // `cxx:93-105`.
    if let Some(bs) = a_curve1.as_bspline_curve() {
        let tolerance = bs.resolution(tl);
        let (index1, index2) = bs.locate_u(u1, tolerance);
        let nb = bs.distinct_knots().len() as i32;
        if index1 > 1 && index2 < nb && index1 == index2 {
            cont1 = bs.degree() as i32 - bs.multiplicity(index1);
        } else {
            cont1 = 5;
        }
    }
    // `cxx:106-121`.
    if let Some(bs) = a_curve2.as_bspline_curve() {
        let tolerance = bs.resolution(tl);
        let (index1, index2) = bs.locate_u(u2, tolerance);
        let nb = bs.distinct_knots().len() as i32;
        if index1 > 1 && index2 < nb && index1 == index2 {
            cont2 = bs.degree() as i32 - bs.multiplicity(index1);
        } else {
            cont2 = 5;
        }
    }

    // `cxx:122-146`.
    let mut n1 = 0;
    if cont1 >= 5 {
        n1 = 3;
    } else if cont1 == 4 {
        n1 = 2;
    } else if cont1 == 2 {
        n1 = 1;
    }
    let mut n2 = 0;
    if cont2 >= 5 {
        n2 = 3;
    } else if cont2 == 4 {
        n2 = 2;
    } else if cont2 == 2 {
        n2 = 1;
    }

    let mut clp1 = ClProps::new(c1, u1, n1, tl)?;
    let mut clp2 = ClProps::new(c2, u2, n2, tl)?;
    // `cxx:149-152`.
    if clp1.value().square_distance(&clp2.value()) > tl * tl {
        return Err("GeomLProp::Continuity: Courbes non jointives".into());
    }

    let min = n1.min(n2);
    let mut fini = false;
    if min >= 1 {
        let mut d1 = clp1.d1();
        let mut d2 = clp2.d1();
        if r1 {
            d1.reverse();
        }
        if r2 {
            d2.reverse();
        }
        if d1.is_equal_tol(&d2, tl, ta) {
            cont = Shape::C1;
        } else if clp1.is_tangent_defined() && clp2.is_tangent_defined() {
            let mut dir1 = clp1.tangent()?;
            let mut dir2 = clp2.tangent()?;
            if r1 {
                dir1.reverse();
            }
            if r2 {
                dir2.reverse();
            }
            if dir1.is_equal_tol(&dir2, ta) {
                cont = Shape::G1;
            }
            fini = true;
        } else {
            fini = true;
        }
    }

    // `cxx:193-201`.
    if min >= 2 && !fini {
        let d1 = clp1.d2();
        let d2 = clp2.d2();
        if d1.is_equal_tol(&d2, tl, ta) {
            cont = Shape::C2;
        }
    }
    Ok(cont)
}

/// `GeomLProp::Continuity(C1, C2, u1, u2, r1, r2)` (`GeomLProp.cxx:207-213`):
/// `Precision::Confusion()` / `Precision::Angular()` defaults.
pub fn continuity_default(
    c1: &dyn Curve,
    c2: &dyn Curve,
    u1: f64,
    u2: f64,
    r1: bool,
    r2: bool,
) -> Result<Shape, String> {
    continuity(c1, c2, u1, u2, r1, r2, CONFUSION, ANGULAR)
}
