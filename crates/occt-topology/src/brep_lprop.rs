//! `BRepLProp::Continuity` (`BRepLProp.cxx:28-127`).
//!
//! `BRepLProp_CLProps` is the `LProp` template on a `BRepAdaptor_Curve`.
//! Derivatives come from the curve. `IsTangentDefined` and `Tangent` follow
//! `LProp_CurveUtils.hxx:333-218`. A third derivative that the curve does not
//! override stays the trait default (zero).

use occt_core::gp::{GpDir, GpVec};
use occt_geom::curve::Curve;

use crate::abs::Orientation;

const REAL_LAST: f64 = f64::MAX;
const REAL_FIRST: f64 = -f64::MAX;

/// `GeomAbs_C0` .. `GeomAbs_CN` as the `u8` values already stored by
/// `Curve::continuity`.
pub const C0: u8 = 0;
pub const G1: u8 = 1;
pub const C1: u8 = 2;
pub const C2: u8 = 4;
pub const CN: u8 = 6;

/// Points farther than `tl` throw `Standard_Failure` (`BRepLProp.cxx:67-69`).
pub fn continuity(
    c1: &dyn Curve,
    ori1: Orientation,
    u1: f64,
    range1: (f64, f64),
    c2: &dyn Curve,
    ori2: Orientation,
    u2: f64,
    range2: (f64, f64),
    same_edge: bool,
    tl: f64,
    ta: f64,
) -> Result<u8, ()> {
    let n1 = deriv_order(c1.continuity());
    let n2 = deriv_order(c2.continuity());
    let p1 = c1.value(u1);
    let p2 = c2.value(u2);
    if p1.distance(&p2) > tl {
        return Err(());
    }
    let mut cont = C0;
    let min = n1.min(n2);
    let mut fini = false;
    if min >= 1 {
        let mut d1 = c1.d1(u1).1;
        let mut d2 = c2.d1(u2).1;
        if ori1 == Orientation::Reversed {
            d1 = d1.reversed();
        }
        if ori2 == Orientation::Reversed {
            d2 = d2.reversed();
        }
        if vec_is_equal(&d1, &d2, tl, ta) {
            cont = C1;
        } else if let (Some(dir1), Some(dir2)) = (
            tangent_dir(c1, u1, range1, n1, tl),
            tangent_dir(c2, u2, range2, n2, tl),
        ) {
            let mut dir1 = dir1;
            let mut dir2 = dir2;
            if ori1 == Orientation::Reversed {
                dir1.reverse();
            }
            if ori2 == Orientation::Reversed {
                dir2.reverse();
            }
            if dir1.angle(&dir2) <= ta {
                cont = G1;
            }
            fini = true;
        } else {
            fini = true;
        }
    }
    if min >= 2 && !fini {
        let d1 = c1.d2(u1).2;
        let d2 = c2.d2(u2).2;
        if vec_is_equal(&d1, &d2, tl, ta) {
            cont = C2;
        }
    }
    if same_edge && c1.is_periodic() && cont >= G1 {
        cont = CN;
    }
    Ok(cont)
}

/// `BRepLProp.cxx:41-64`: `C3`/`CN` request 3, `C2` requests 2, `C1` requests 1.
fn deriv_order(cont: u8) -> i32 {
    if cont >= 5 {
        3
    } else if cont == 4 {
        2
    } else if cont == 2 {
        1
    } else {
        0
    }
}

/// `gp_Vec::IsEqual` (`gp_Vec.cxx:32-48`).
fn vec_is_equal(a: &GpVec, b: &GpVec, linear: f64, angular: f64) -> bool {
    let ma = a.magnitude();
    let mb = b.magnitude();
    let diff = (ma - mb).abs();
    if ma <= linear || mb <= linear {
        diff <= linear
    } else {
        diff <= linear && a.angle(b) <= angular
    }
}

/// `IsTangentDefined` (`LProp_CurveUtils.hxx:333-379`) then `ComputeTangent`
/// (`:175-218`). `None` is an undefined tangent (`LProp_Undefined`).
fn tangent_dir(
    curve: &dyn Curve,
    u: f64,
    range: (f64, f64),
    cn: i32,
    lin_tol: f64,
) -> Option<GpDir> {
    let tol_sq = lin_tol * lin_tol;
    let mut sig = 0;
    let mut deriv = GpVec::zero();
    for order in 1..=4 {
        if cn < order {
            return None;
        }
        let av = match order {
            1 => curve.d1(u).1,
            2 => curve.d2(u).2,
            3 => curve.d3(u).3,
            _ => return None,
        };
        if av.square_magnitude() > tol_sq {
            sig = order;
            deriv = av;
            break;
        }
    }
    if sig == 0 {
        return None;
    }
    if sig == 1 {
        return GpDir::from_vec(&deriv).ok();
    }
    let (u_inf, u_sup) = range;
    let du = if u_sup >= REAL_LAST || u_inf <= REAL_FIRST {
        0.0
    } else {
        u_sup - u_inf
    };
    let delta = (du * 1.0e-3).max(1.0e-7);
    let other = if u - u_inf < delta { u + delta } else { u - delta };
    let p1 = curve.value(u.min(other));
    let p2 = curve.value(u.max(other));
    let chord = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
    if deriv.dot(&chord) < 0.0 {
        deriv = deriv.reversed();
    }
    GpDir::from_vec(&deriv).ok()
}
