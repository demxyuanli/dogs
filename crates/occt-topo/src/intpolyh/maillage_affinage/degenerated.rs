//! `IsDegenerated` / `DegeneratedIndex` (`IntPolyh_MaillageAffinage.cxx:3863-3996`):
//! detection of the degenerated iso-lines of a surface, used when filling the
//! arrays of points so that a whole collapsed row/column is flagged.

use occt_geom::Surface;

use super::MY_TOLERANCE;

/// `IsDegenerated(aS, aIndex, aTol2, aDegX)` (`cxx:3920-3996`).
///
/// `aIndex` selects the iso-line: 1/2 -> `V = v1` / `V = v2`, 3/4 ->
/// `U = u1` / `U = u2`. Returns `(bRet, aDegX)`.
///
/// The `.cxx:3954` / `:3979` compute the running parameter as `i * dU`
/// (`i * dV`) instead of `aU1 + i * dU`, i.e. the walked range is anchored at
/// zero rather than at the first parameter. The port keeps that expression.
pub fn is_degenerated(surface: &dyn Surface, index: i32, tol2: f64) -> (bool, f64) {
    let mut ret = false;
    let nb_p: i32 = 3;
    let mut deg_x = 99.0;

    let (u1, u2) = surface.u_range();
    let (v1, v2) = surface.v_range();

    if index < 3 {
        // V = const.
        let v = if index == 2 { v2 } else { v1 };
        let du = (u2 - u1) / f64::from(nb_p - 1);
        let mut u = u1;
        let mut p1 = surface.value(u, v);
        for i in 1..nb_p {
            u = f64::from(i) * du;
            if i == nb_p - 1 {
                u = u2;
            }
            let p2 = surface.value(u, v);
            let d2 = p1.square_distance(&p2);
            if d2 > tol2 {
                return (ret, deg_x);
            }
            p1 = p2;
        }
        deg_x = v;
        ret = !ret;
    } else {
        // U = const.
        let u = if index == 4 { u2 } else { u1 };
        let dv = (v2 - v1) / f64::from(nb_p - 1);
        let mut v = v1;
        let mut p1 = surface.value(u, v);
        for i in 1..nb_p {
            v = f64::from(i) * dv;
            if i == nb_p - 1 {
                v = v2;
            }
            let p2 = surface.value(u, v);
            let d2 = p1.square_distance(&p2);
            if d2 > tol2 {
                return (ret, deg_x);
            }
            p1 = p2;
        }
        ret = !ret;
        deg_x = u;
    }

    (ret, deg_x)
}

/// `DegeneratedIndex(aXpars, aNbX, aS, aIsoDirection, aI1, aI2)`
/// (`cxx:3865-3915`): 1-based indices of the degenerated parameters.
pub fn degenerated_index(
    x_pars: &[f64],
    nb_x: i32,
    surface: &dyn Surface,
    iso_direction: i32,
    i1: &mut i32,
    i2: &mut i32,
) {
    *i1 = 0;
    *i2 = 0;
    let tol2 = MY_TOLERANCE * MY_TOLERANCE;

    let (deg_x1, deg_x2) = if iso_direction == 1 {
        // V = const.
        let (b1, d1) = is_degenerated(surface, 1, tol2);
        let (b2, d2) = is_degenerated(surface, 2, tol2);
        (if b1 { Some(d1) } else { None }, if b2 { Some(d2) } else { None })
    } else if iso_direction == 2 {
        // U = const.
        let (b1, d1) = is_degenerated(surface, 3, tol2);
        let (b2, d2) = is_degenerated(surface, 4, tol2);
        (if b1 { Some(d1) } else { None }, if b2 { Some(d2) } else { None })
    } else {
        return;
    };

    if deg_x1.is_none() && deg_x2.is_none() {
        return;
    }

    for i in 1..=nb_x {
        if let Some(d) = deg_x1 {
            if (x_pars[(i - 1) as usize] - d).abs() < MY_TOLERANCE {
                *i1 = i;
            }
        }
        if let Some(d) = deg_x2 {
            if (x_pars[(i - 1) as usize] - d).abs() < MY_TOLERANCE {
                *i2 = i;
            }
        }
    }
}
