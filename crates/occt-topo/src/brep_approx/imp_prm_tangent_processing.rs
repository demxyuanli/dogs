//! The two tangent post-processors of `ApproxInt_ImpPrmSvSurfaces` — the
//! branches that survive a singular (degenerate) surface point.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_ImpPrmSvSurfaces.gxx`
//! - `IsSingular` (`:32-56`)
//! - `SingularProcessing` (`:80-267`)
//! - `NonSingularProcessing` (`:287-322`)
//!
//! `MyZerImpFunc` uses neither: these are file-static helpers of the `.gxx`,
//! called only from `Compute` (`:640-743`), so they live here as
//! `pub(crate)` functions of the same module namespace.
//!
//! `theAngTol` is used as `sin(theAngTol) ~ theAngTol` throughout the `.gxx`
//! (cross-product comparisons), which the port keeps.

use occt_core::gp::{GpVec, GpVec2d};

/// `IsSingular(theDU, theDV, theSqLinTol, theAngTol)` (`gxx:32-56`): true when
/// either derivative is shorter than `theSqLinTol`, or when the two are
/// parallel within `theAngTol`.
pub(crate) fn is_singular(du: &GpVec, dv: &GpVec, sq_lin_tol: f64, ang_tol: f64) -> bool {
    let a_du = *du;
    let a_dv = *dv;

    let a_sq_magn_du = a_du.square_magnitude();
    let a_sq_magn_dv = a_dv.square_magnitude();

    if a_sq_magn_du < sq_lin_tol {
        return true;
    }
    let a_du = a_du.divided(a_sq_magn_du.sqrt());

    if a_sq_magn_dv < sq_lin_tol {
        return true;
    }
    let a_dv = a_dv.divided(a_sq_magn_dv.sqrt());

    // Here aDU and aDV vectors have magnitude 1.0.
    a_du.crossed(&a_dv).square_magnitude() < ang_tol * ang_tol
}

/// `SingularProcessing(theDU, theDV, theIsTo3DTgCompute, theLinTol, theAngTol,
/// theTg3D, theTg2D)` (`gxx:80-267`).
///
/// The four cases are, in OCCT's order: both derivatives vanish (fail); only
/// `theDU` vanishes; only `theDV` vanishes; both are parallel. With
/// `theIsTo3DTgCompute == false` the 3D tangent is already defined and is only
/// checked to be parallel to the non-vanishing derivative.
pub(crate) fn singular_processing(
    the_du: &GpVec,
    the_dv: &GpVec,
    the_is_to3d_tg_compute: bool,
    the_lin_tol: f64,
    the_ang_tol: f64,
    the_tg3d: &mut GpVec,
    the_tg2d: &mut GpVec2d,
) -> bool {
    // theTg3D has to be normalized when it is already defined.
    let a_sq_tan = the_tg3d.square_magnitude();

    let a_sq_magn_du = the_du.square_magnitude();
    let a_sq_magn_dv = the_dv.square_magnitude();

    // 1.
    if a_sq_magn_du < the_lin_tol && a_sq_magn_dv < the_lin_tol {
        return false;
    }

    // 2.
    if a_sq_magn_du < the_lin_tol {
        if the_is_to3d_tg_compute {
            // theTg3D will be normalized. Its magnitude is
            let a_tg_magn = 1.0;

            let a_norm = a_sq_magn_dv.sqrt();
            *the_tg3d = the_dv.divided(a_norm);
            the_tg2d.set_coord(0.0, a_tg_magn / a_norm);
        } else {
            // theTg3D is already defined; check it is parallel to theDV.
            if the_dv.crossed(the_tg3d).square_magnitude()
                < the_ang_tol * the_ang_tol * a_sq_magn_dv * a_sq_tan
            {
                // Use sign "+" when codirectional, "-" when opposite.
                let a_dp = the_tg3d.dot(the_dv);
                the_tg2d
                    .set_coord(0.0, (a_sq_tan / a_sq_magn_dv).sqrt().copysign(a_dp));
            } else {
                return false;
            }
        }
        return true;
    }

    // 3.
    if a_sq_magn_dv < the_lin_tol {
        if the_is_to3d_tg_compute {
            let a_tg_magn = 1.0;

            let a_norm = a_sq_magn_du.sqrt();
            *the_tg3d = the_du.divided(a_norm);
            the_tg2d.set_coord(a_tg_magn / a_norm, 0.0);
        } else {
            if the_du.crossed(the_tg3d).square_magnitude()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                let a_dp = the_tg3d.dot(the_du);
                the_tg2d
                    .set_coord((a_sq_tan / a_sq_magn_du).sqrt().copysign(a_dp), 0.0);
            } else {
                return false;
            }
        }
        return true;
    }

    // 4. aSqMagnDU > 0.0 && aSqMagnDV > 0.0 but theDV || theDU.
    let a_len_u = a_sq_magn_du.sqrt();
    let a_len_v = a_sq_magn_dv.sqrt();

    // aLenSum > 0.0 definitely
    let a_len_sum = a_len_u + a_len_v;

    if the_dv.dot(the_du) > 0.0 {
        // Vectors theDV and theDU are codirectional.
        if the_is_to3d_tg_compute {
            the_tg2d.set_coord(1.0 / a_len_sum, 1.0 / a_len_sum);
            *the_tg3d = the_du
                .multiplied_scalar(the_tg2d.coord.x)
                .added(&the_dv.multiplied_scalar(the_tg2d.coord.y));
        } else {
            if the_du.crossed(the_tg3d).square_magnitude()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                let a_dp = the_tg3d.dot(the_du);
                let a_len_tg = a_sq_tan.sqrt().copysign(a_dp);
                the_tg2d.set_coord(a_len_tg / a_len_sum, a_len_tg / a_len_sum);
            } else {
                return false;
            }
        }
    } else {
        // Vectors theDV and theDU are opposite.
        if the_is_to3d_tg_compute {
            the_tg2d.set_coord(1.0 / a_len_sum, -1.0 / a_len_sum);
            *the_tg3d = the_du
                .multiplied_scalar(the_tg2d.coord.x)
                .added(&the_dv.multiplied_scalar(the_tg2d.coord.y));
        } else {
            if the_du.crossed(the_tg3d).square_magnitude()
                < the_ang_tol * the_ang_tol * a_sq_magn_du * a_sq_tan
            {
                let a_dp = the_tg3d.dot(the_du);
                let a_len_tg = a_sq_tan.sqrt().copysign(a_dp);
                the_tg2d
                    .set_coord(a_len_tg / a_len_sum, -a_len_tg / a_len_sum);
            } else {
                return false;
            }
        }
    }

    true
}

/// `NonSingularProcessing(theDU, theDV, theTg3D, theLinTol, theAngTol,
/// theTg2D)` (`gxx:287-322`).
pub(crate) fn non_singular_processing(
    the_du: &GpVec,
    the_dv: &GpVec,
    the_tg3d: &GpVec,
    the_lin_tol: f64,
    the_ang_tol: f64,
    the_tg2d: &mut GpVec2d,
) -> bool {
    let a_normal = the_du.crossed(the_dv);
    let a_sq_magn = a_normal.square_magnitude();

    if is_singular(the_du, the_dv, the_lin_tol, the_ang_tol) {
        let mut a_tg = *the_tg3d;
        return singular_processing(
            the_du,
            the_dv,
            false,
            the_lin_tol,
            the_ang_tol,
            &mut a_tg,
            the_tg2d,
        );
    }

    // If T = A*U + B*V then A x T = (A x B)*V and B x T = (B x A)*U.
    let a_tg_u = the_tg3d.crossed(the_du);
    let a_tg_v = the_tg3d.crossed(the_dv);
    let a_delta_u = a_tg_v.square_magnitude() / a_sq_magn;
    let a_delta_v = a_tg_u.square_magnitude() / a_sq_magn;

    the_tg2d.set_coord(
        a_delta_u.sqrt().copysign(a_tg_v.dot(&a_normal)),
        -a_delta_v.sqrt().copysign(a_tg_u.dot(&a_normal)),
    );

    true
}
