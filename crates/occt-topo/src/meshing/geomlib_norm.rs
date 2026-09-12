//! `GeomLib::NormEstim`. Source: `GeomLib.cxx:2562-2689`.

use occt_core::cslib::{normal_d2, CSLibNormalStatus};
use occt_core::gp::GpDir;
use occt_geom::Surface;

const CONE_SINGULARITY_ANGLE_EPS: f64 = 1.0e-4;
const STEP: f64 = 1.0e-5;
const EPS: f64 = 1.0e-16;

/// `GeomLib::NormEstim(theSurf, theUV, theTol, theNorm)`.
/// Returns the C++ status (0 = defined, 1 = quasisingular, 2/3 = fail)
/// and the written normal when status is 0 or 1.
pub fn norm_estim(surface: &dyn Surface, u: f64, v: f64, tol: f64) -> (i32, Option<GpDir>) {
    let tol2 = tol * tol;
    let (_, du, dv) = surface.d1(u, v);
    let mdu = du.square_magnitude();
    let mdv = dv.square_magnitude();
    if mdu >= tol2 && mdv >= tol2 {
        let n = du.crossed(&dv);
        if n.square_magnitude() < tol2 {
            return (3, None);
        }
        return match GpDir::from_vec(&n) {
            Ok(dir) => (0, Some(dir)),
            Err(_) => (3, None),
        };
    }

    let (_, du, dv, d2u, d2v, d2uv) = surface.d2(u, v);
    let (done, status, normal) = normal_d2(&du, &dv, &d2u, &d2v, &d2uv, tol);
    if !done {
        return (
            if status == CSLibNormalStatus::D1NIsNull {
                2
            } else {
                3
            },
            None,
        );
    }
    let mut normal = match normal {
        Some(n) => n,
        None => return (3, None),
    };

    let (umin, umax) = surface.u_range();
    let (vmin, vmax) = surface.v_range();
    let mut sign = -1.0;

    if (v > vmin + STEP) && (v < vmax - STEP) {
        let (_, du_m, dv_m) = surface.d1(u, v - sign * STEP);
        if du_m.xyz().square_modulus() > EPS && dv_m.xyz().square_modulus() > EPS {
            if let Ok(n1) = GpDir::from_vec(&du_m.crossed(&dv_m)) {
                let (_, du_p, dv_p) = surface.d1(u, v + sign * STEP);
                if du_p.xyz().square_modulus() > EPS && dv_p.xyz().square_modulus() > EPS {
                    if let Ok(n2) = GpDir::from_vec(&du_p.crossed(&dv_p)) {
                        if std::f64::consts::PI - n1.angle(&n2) <= CONE_SINGULARITY_ANGLE_EPS {
                            return (2, None);
                        }
                    }
                }
            }
        }
    }

    if mdu < tol2 && mdv >= tol2 {
        if (vmax - v) > (v - vmin) {
            sign = 1.0;
        }
        let (_, du_s, dv_s) = surface.d1(u, v + sign * STEP);
        let mut norm = du_s.crossed(&dv_s);
        if norm.square_magnitude() < EPS {
            let sign1 = if (umax - u) > (u - umin) { 1.0 } else { -1.0 };
            let (_, du_o, dv_o) = surface.d1(u + sign1 * STEP, v + sign * STEP);
            norm = du_o.crossed(&dv_o);
        }
        if norm.square_magnitude() >= EPS && norm.xyz().dot(normal.xyz()) < 0.0 {
            normal.reverse();
        }
    }

    if mdv < tol2 && mdu >= tol2 {
        if (umax - u) > (u - umin) {
            sign = 1.0;
        }
        let (_, du_s, dv_s) = surface.d1(u + sign * STEP, v);
        let mut norm = du_s.crossed(&dv_s);
        if norm.square_magnitude() < EPS {
            let sign1 = if (vmax - v) > (v - vmin) { 1.0 } else { -1.0 };
            let (_, du_o, dv_o) = surface.d1(u + sign * STEP, v + sign1 * STEP);
            norm = du_o.crossed(&dv_o);
        }
        if norm.square_magnitude() >= EPS && norm.xyz().dot(normal.xyz()) < 0.0 {
            normal.reverse();
        }
    }

    if status == CSLibNormalStatus::D1NuIsNull
        || status == CSLibNormalStatus::D1NvIsNull
        || status == CSLibNormalStatus::D1NuIsParallelD1Nv
    {
        return (1, Some(normal));
    }

    (
        if status == CSLibNormalStatus::InfinityOfSolutions {
            2
        } else {
            3
        },
        None,
    )
}
