//! The 2D jets of `ProjLib_CompProjectedCurve`: `d1`, `d2`, `d2CurvOnSurf`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:137-305`
//! (OCCT 8.0.0).
//!
//! These are the implicit derivatives of the normal-projection system. With
//! `E(t, u, v) = |S(u, v) - C(t)|^2 / 2` OCCT writes the 2x2 system
//! `dE/du = 0`, `dE/dv = 0` and differentiates it in `t`:
//!
//! ```text
//! [dE_du_du  dE_du_dv] [du/dt]     [-dE_dt_u]
//! [dE_dv_du  dE_dv_dv] [dv/dt]  =  [-dE_dt_v]
//! ```
//!
//! so `(du/dt, dv/dt) = -M * dE_dt` with `M` the inverse of the Hessian. `d2`
//! differentiates once more (the second derivative of `E` in `t` plus the two
//! chain-rule terms `2 (d2E/dt dX) dX/dt` and `(d2E/dX2)(dX/dt)^2`).
//!
//! `d2CurvOnSurf` returns the same two derivatives carried into 3D
//! (`cxx:240-303`): `V1 = S_u u' + S_v v'` and the second derivative of `S`
//! along the projected curve.
//!
//! Translation notes:
//! - `gp_Resolution()` (`cxx:158`, `cxx:187`, `cxx:257`) is
//!   [`occt_core::precision::RESOLUTION`]; a singular Hessian raises in OCCT
//!   (`Standard_ConstructionError`), so the port answers `Err`.
//! - `gp_Vec2d * gp_Vec2d` in OCCT is the dot product (`gp_Vec2d.hxx`), which
//!   is [`GpVec2d::dot`] here.
//! - every `gp_Vec * gp_Vec` is a dot product, i.e. [`GpVec::dot`].

use occt_core::gp::{GpVec, GpVec2d};
use occt_core::precision::RESOLUTION;

use crate::{Curve, Surface};

/// `d1(t, u, v, V, Curve, Surface)` (`cxx:137-166`): the derivative of the
/// projection parameters in `t`.
pub(crate) fn d1(
    t: f64,
    u: f64,
    v: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
) -> Result<GpVec2d, String> {
    let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv) = surface.d2(u, v);
    let (c, dc1_t) = curve.d1(t);
    // `gp_Vec Ort(C, S)` (`cxx:149`): `S - C`.
    let ort = GpVec::from_pnts(&c, &s);

    let de_dt = GpVec2d::new(-dc1_t.dot(&ds1_u), -dc1_t.dot(&ds1_v));
    let de_du = GpVec2d::new(
        ds1_u.dot(&ds1_u) + ort.dot(&ds2_u),
        ds1_u.dot(&ds1_v) + ort.dot(&ds2_uv),
    );
    let de_dv = GpVec2d::new(
        ds1_v.dot(&ds1_u) + ort.dot(&ds2_uv),
        ds1_v.dot(&ds1_v) + ort.dot(&ds2_v),
    );

    let (m11, m12, m21, m22) = inverse_hessian(&de_du, &de_dv)?;
    Ok(GpVec2d::new(
        -(m11 * de_dt.x() + m12 * de_dt.y()),
        -(m21 * de_dt.x() + m22 * de_dt.y()),
    ))
}

/// `d2(t, u, v, V1, V2, Curve, Surface)` (`cxx:168-230`).
pub(crate) fn d2(
    t: f64,
    u: f64,
    v: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
) -> Result<(GpVec2d, GpVec2d), String> {
    let d3 = surface
        .d3(u, v)
        .ok_or_else(|| "ProjLib_CompProjectedCurve: surface D3 unavailable".to_string())?;
    let (c, dc1_t, dc2_t) = curve.d2(t);
    let ort = GpVec::from_pnts(&c, &d3.point);

    let de_dt = GpVec2d::new(-dc1_t.dot(&d3.d1u), -dc1_t.dot(&d3.d1v));
    let de_du = GpVec2d::new(
        d3.d1u.dot(&d3.d1u) + ort.dot(&d3.d2u),
        d3.d1u.dot(&d3.d1v) + ort.dot(&d3.d2uv),
    );
    let de_dv = GpVec2d::new(
        d3.d1v.dot(&d3.d1u) + ort.dot(&d3.d2uv),
        d3.d1v.dot(&d3.d1v) + ort.dot(&d3.d2v),
    );
    let (m11, m12, m21, m22) = inverse_hessian(&de_du, &de_dv)?;

    // First derivative (`cxx:196`).
    let v1 = GpVec2d::new(
        -(m11 * de_dt.x() + m12 * de_dt.y()),
        -(m21 * de_dt.x() + m22 * de_dt.y()),
    );

    // `d2E_dt` (`cxx:200`).
    let d2e_dt = GpVec2d::new(-dc2_t.dot(&d3.d1u), -dc2_t.dot(&d3.d1v));

    // `S2 = 2 * (d2E/dtdX)(dX/dt)` (`cxx:203-206`).
    let d2e1_dtdx = GpVec2d::new(-dc1_t.dot(&d3.d2u), -dc1_t.dot(&d3.d2uv));
    let d2e2_dtdx = GpVec2d::new(-dc1_t.dot(&d3.d2uv), -dc1_t.dot(&d3.d2v));
    let s2 = GpVec2d::new(
        2.0 * d2e1_dtdx.dot(&v1),
        2.0 * d2e2_dtdx.dot(&v1),
    );

    // `S3 = (d2E/dX2)(dX/dt)^2` (`cxx:208-223`).
    let tmp1 = 2.0 * d3.d1u.dot(&d3.d2uv) + d3.d1v.dot(&d3.d2u) + ort.dot(&d3.d3_uuv);
    let row11 = GpVec2d::new(3.0 * d3.d1u.dot(&d3.d2u) + ort.dot(&d3.d3u), tmp1);
    let row12 = GpVec2d::new(
        tmp1,
        d3.d2v.dot(&d3.d1u) + 2.0 * d3.d1v.dot(&d3.d2uv) + ort.dot(&d3.d3_uvv),
    );
    let tmp2 = 2.0 * d3.d2uv.dot(&d3.d1v) + d3.d1u.dot(&d3.d2v) + ort.dot(&d3.d3_uvv);
    let row21 = GpVec2d::new(
        d3.d2u.dot(&d3.d1v) + 2.0 * d3.d1u.dot(&d3.d2uv) + ort.dot(&d3.d3_uuv),
        tmp2,
    );
    let row22 = GpVec2d::new(tmp2, 3.0 * d3.d1v.dot(&d3.d2v) + ort.dot(&d3.d3v));
    let s3 = GpVec2d::new(
        v1.dot(&GpVec2d::new(row11.dot(&v1), row12.dot(&v1))),
        v1.dot(&GpVec2d::new(row21.dot(&v1), row22.dot(&v1))),
    );

    let sum = GpVec2d::new(
        d2e_dt.x() + s2.x() + s3.x(),
        d2e_dt.y() + s2.y() + s3.y(),
    );
    let v2 = GpVec2d::new(
        -(m11 * sum.x() + m12 * sum.y()),
        -(m21 * sum.x() + m22 * sum.y()),
    );
    Ok((v1, v2))
}

/// `d2CurvOnSurf(t, u, v, V1, V2, Curve, Surface)` (`cxx:235-303`): the first
/// and second 3D derivatives of the projected curve on the surface.
pub(crate) fn d2_curv_on_surf(
    t: f64,
    u: f64,
    v: f64,
    curve: &dyn Curve,
    surface: &dyn Surface,
) -> Result<(GpVec, GpVec), String> {
    let d3 = surface
        .d3(u, v)
        .ok_or_else(|| "ProjLib_CompProjectedCurve: surface D3 unavailable".to_string())?;
    let (c, dc1_t, dc2_t) = curve.d2(t);
    let ort = GpVec::from_pnts(&c, &d3.point);

    let de_dt = GpVec2d::new(-dc1_t.dot(&d3.d1u), -dc1_t.dot(&d3.d1v));
    let de_du = GpVec2d::new(
        d3.d1u.dot(&d3.d1u) + ort.dot(&d3.d2u),
        d3.d1u.dot(&d3.d1v) + ort.dot(&d3.d2uv),
    );
    let de_dv = GpVec2d::new(
        d3.d1v.dot(&d3.d1u) + ort.dot(&d3.d2uv),
        d3.d1v.dot(&d3.d1v) + ort.dot(&d3.d2v),
    );
    let (m11, m12, m21, m22) = inverse_hessian(&de_du, &de_dv)?;

    let v12d = GpVec2d::new(
        -(m11 * de_dt.x() + m12 * de_dt.y()),
        -(m21 * de_dt.x() + m22 * de_dt.y()),
    );

    let d2e_dt = GpVec2d::new(-dc2_t.dot(&d3.d1u), -dc2_t.dot(&d3.d1v));
    let d2e1_dtdx = GpVec2d::new(-dc1_t.dot(&d3.d2u), -dc1_t.dot(&d3.d2uv));
    let d2e2_dtdx = GpVec2d::new(-dc1_t.dot(&d3.d2uv), -dc1_t.dot(&d3.d2v));
    let s2 = GpVec2d::new(
        2.0 * d2e1_dtdx.dot(&v12d),
        2.0 * d2e2_dtdx.dot(&v12d),
    );

    let tmp1 = 2.0 * d3.d1u.dot(&d3.d2uv) + d3.d1v.dot(&d3.d2u) + ort.dot(&d3.d3_uuv);
    let row11 = GpVec2d::new(3.0 * d3.d1u.dot(&d3.d2u) + ort.dot(&d3.d3u), tmp1);
    let row12 = GpVec2d::new(
        tmp1,
        d3.d2v.dot(&d3.d1u) + 2.0 * d3.d1v.dot(&d3.d2uv) + ort.dot(&d3.d3_uvv),
    );
    let tmp2 = 2.0 * d3.d2uv.dot(&d3.d1v) + d3.d1u.dot(&d3.d2v) + ort.dot(&d3.d3_uvv);
    let row21 = GpVec2d::new(
        d3.d2u.dot(&d3.d1v) + 2.0 * d3.d1u.dot(&d3.d2uv) + ort.dot(&d3.d3_uuv),
        tmp2,
    );
    let row22 = GpVec2d::new(tmp2, 3.0 * d3.d1v.dot(&d3.d2v) + ort.dot(&d3.d3v));
    let s3 = GpVec2d::new(
        v12d.dot(&GpVec2d::new(row11.dot(&v12d), row12.dot(&v12d))),
        v12d.dot(&GpVec2d::new(row21.dot(&v12d), row22.dot(&v12d))),
    );

    let sum = GpVec2d::new(
        d2e_dt.x() + s2.x() + s3.x(),
        d2e_dt.y() + s2.y() + s3.y(),
    );
    let v22d = GpVec2d::new(
        -(m11 * sum.x() + m12 * sum.y()),
        -(m21 * sum.x() + m22 * sum.y()),
    );

    // `V1 = DS1_u * V12d.X() + DS1_v * V12d.Y()` (`cxx:299`).
    let v1 = GpVec::zero()
        .set_linear_form_2(v12d.x(), &d3.d1u, v12d.y(), &d3.d1v);
    // `V2 = DS2_u u'^2 + DS1_u u'' + 2 DS2_uv u' v' + DS2_v v'^2 + DS1_v v''`
    // (`cxx:300-302`).
    let v2 = GpVec::zero()
        .set_linear_form_2(v12d.x() * v12d.x(), &d3.d2u, v12d.y() * v12d.y(), &d3.d2v)
        .added(&d3.d2uv.multiplied_scalar(2.0 * v12d.x() * v12d.y()))
        .added(&d3.d1u.multiplied_scalar(v22d.x()))
        .added(&d3.d1v.multiplied_scalar(v22d.y()));
    Ok((v1, v2))
}

/// `gp_Mat2d M(dE_dv.Y()/det, -dE_du.Y()/det, -dE_dv.X()/det, dE_du.X()/det)`
/// (`cxx:161`), read as its two rows. `Err` is the
/// `Standard_ConstructionError` of `|det| < gp::Resolution()`.
fn inverse_hessian(de_du: &GpVec2d, de_dv: &GpVec2d) -> Result<(f64, f64, f64, f64), String> {
    let det = de_du.x() * de_dv.y() - de_du.y() * de_dv.x();
    if det.abs() < RESOLUTION {
        return Err("ProjLib_CompProjectedCurve: singular Hessian (cxx:158)".into());
    }
    Ok((
        de_dv.y() / det,
        -de_du.y() / det,
        -de_dv.x() / det,
        de_du.x() / det,
    ))
}