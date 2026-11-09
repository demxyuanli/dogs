//! `Contap_SurfProps` - the point / normal / normal-derivative evaluations the
//! contour functions are built on.
//!
//! Sources: `Contap_SurfProps.hxx:24-58`, `Contap_SurfProps.cxx:25-351`.
//!
//! Every arm reproduces the analytic formulas of the `.cxx`; the `default` arm
//! is `D1` / `D2` on the adaptor, exactly as `Adaptor3d_HSurfaceTool` forwards.
//! A surface that was never loaded answers `None` where OCCT dereferences a
//! null handle.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::REAL_EPSILON;

use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `gp_Vec` of a `gp_Dir`.
fn dir_vec(d: &occt_core::gp::GpDir) -> GpVec {
    GpVec::from_xyz(d.xyz())
}

/// `Contap_SurfProps::Normale(S, U, V, P, N)` (`cxx:25-127`).
pub fn normale(s: &SurfaceAdaptor, u: f64, v: f64) -> Option<(GpPnt, GpVec)> {
    match s.get_type() {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane()?;
            let mut norm = dir_vec(&pl.axis().direction());
            let p = s.d0(u, v)?;
            if !pl.is_direct() {
                norm.reverse();
            }
            Some((p, norm))
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere()?;
            let p = s.d0(u, v)?;
            let mut norm = GpVec::from_pnts(&sp.location(), &p);
            if sp.is_direct() {
                norm = norm.divide(sp.radius());
            } else {
                norm = norm.divide(-sp.radius());
            }
            Some((p, norm))
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder()?;
            let p = s.d0(u, v)?;
            let mut norm = GpVec::new(0.0, 0.0, 0.0).set_linear_form_2(
                u.cos(),
                &dir_vec(&cy.x_axis().direction()),
                u.sin(),
                &dir_vec(&cy.y_axis().direction()),
            );
            if !cy.is_direct() {
                norm.reverse();
            }
            Some((p, norm))
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone()?;
            let p = s.d0(u, v)?;
            let sina = co.semi_angle().sin();
            let cosa = co.semi_angle().cos();
            let rad = co.radius();
            if (v * sina + rad).abs() <= 1.0e-12 {
                // `cxx:84-94`: at the apex the normal is the zero vector.
                return Some((p, GpVec::new(0.0, 0.0, 0.0)));
            }
            let mut norm = if rad + v * sina < 0.0 {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            } else {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    -sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            };
            if !co.is_direct() {
                norm.reverse();
            }
            Some((p, norm))
        }
        _ => {
            let (p, d1u, d1v) = s.d1(u, v)?;
            Some((p, d1u.crossed(&d1v)))
        }
    }
}

/// `Contap_SurfProps::DerivAndNorm(S, U, V, P, d1u, d1v, N)`
/// (`cxx:133-232`).
pub fn deriv_and_norm(
    s: &SurfaceAdaptor,
    u: f64,
    v: f64,
) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    match s.get_type() {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane()?;
            let mut norm = dir_vec(&pl.axis().direction());
            let (p, d1u, d1v) = s.d1(u, v)?;
            if !pl.is_direct() {
                norm.reverse();
            }
            Some((p, d1u, d1v, norm))
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere()?;
            let (p, d1u, d1v) = s.d1(u, v)?;
            let mut norm = GpVec::from_pnts(&sp.location(), &p);
            if sp.is_direct() {
                norm = norm.divide(sp.radius());
            } else {
                norm = norm.divide(-sp.radius());
            }
            Some((p, d1u, d1v, norm))
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder()?;
            let (p, d1u, d1v) = s.d1(u, v)?;
            let mut norm = GpVec::new(0.0, 0.0, 0.0).set_linear_form_2(
                u.cos(),
                &dir_vec(&cy.x_axis().direction()),
                u.sin(),
                &dir_vec(&cy.y_axis().direction()),
            );
            if !cy.is_direct() {
                norm.reverse();
            }
            Some((p, d1u, d1v, norm))
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone()?;
            let (p, d1u, d1v) = s.d1(u, v)?;
            let sina = co.semi_angle().sin();
            let cosa = co.semi_angle().cos();
            let rad = co.radius();
            let vcalc = cone_vcalc(s, v, sina, rad);
            let mut norm = if rad + vcalc * sina < 0.0 {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            } else {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    -sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            };
            if !co.is_direct() {
                norm.reverse();
            }
            Some((p, d1u, d1v, norm))
        }
        _ => {
            let (p, d1u, d1v) = s.d1(u, v)?;
            let norm = d1u.crossed(&d1v);
            Some((p, d1u, d1v, norm))
        }
    }
}

/// `Contap_SurfProps::NormAndDn(S, U, V, P, N, Dnu, Dnv)` (`cxx:237-351`).
pub fn norm_and_dn(
    s: &SurfaceAdaptor,
    u: f64,
    v: f64,
) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
    match s.get_type() {
        GeomAbsSurfaceType::Plane => {
            let pl = s.plane()?;
            let p = s.d0(u, v)?;
            let mut norm = dir_vec(&pl.axis().direction());
            if !pl.is_direct() {
                norm.reverse();
            }
            Some((
                p,
                norm,
                GpVec::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
            ))
        }
        GeomAbsSurfaceType::Sphere => {
            let sp = s.sphere()?;
            let (p, mut dnu, mut dnv) = s.d1(u, v)?;
            let mut norm = GpVec::from_pnts(&sp.location(), &p);
            let mut rad = sp.radius();
            if !sp.is_direct() {
                rad = -rad;
            }
            norm = norm.divide(rad);
            dnu = dnu.divide(rad);
            dnv = dnv.divide(rad);
            Some((p, norm, dnu, dnv))
        }
        GeomAbsSurfaceType::Cylinder => {
            let cy = s.cylinder()?;
            let p = s.d0(u, v)?;
            let mut norm = GpVec::new(0.0, 0.0, 0.0).set_linear_form_2(
                u.cos(),
                &dir_vec(&cy.x_axis().direction()),
                u.sin(),
                &dir_vec(&cy.y_axis().direction()),
            );
            let mut dnu = GpVec::new(0.0, 0.0, 0.0).set_linear_form_2(
                -u.sin(),
                &dir_vec(&cy.x_axis().direction()),
                u.cos(),
                &dir_vec(&cy.y_axis().direction()),
            );
            if !cy.is_direct() {
                norm.reverse();
                dnu.reverse();
            }
            Some((p, norm, dnu, GpVec::new(0.0, 0.0, 0.0)))
        }
        GeomAbsSurfaceType::Cone => {
            let co = s.cone()?;
            let p = s.d0(u, v)?;
            let sina = co.semi_angle().sin();
            let cosa = co.semi_angle().cos();
            let rad = co.radius();
            let vcalc = cone_vcalc(s, v, sina, rad);
            let mut norm = if rad + vcalc * sina < 0.0 {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            } else {
                GpVec::new(0.0, 0.0, 0.0).set_linear_form_3(
                    -sina,
                    &dir_vec(&co.axis().direction()),
                    cosa * u.cos(),
                    &dir_vec(&co.x_axis().direction()),
                    cosa * u.sin(),
                    &dir_vec(&co.y_axis().direction()),
                )
            };
            let mut dnu = GpVec::new(0.0, 0.0, 0.0).set_linear_form_2(
                -cosa * u.sin(),
                &dir_vec(&co.x_axis().direction()),
                cosa * u.cos(),
                &dir_vec(&co.y_axis().direction()),
            );
            if !co.is_direct() {
                norm.reverse();
                dnu.reverse();
            }
            Some((p, norm, dnu, GpVec::new(0.0, 0.0, 0.0)))
        }
        _ => {
            let (p, d1u, d1v, d2u, d2v, d2uv) = s.d2(u, v)?;
            let norm = d1u.crossed(&d1v);
            let dnu = d2u.crossed(&d1v).add(&d1u.crossed(&d2uv));
            let dnv = d2uv.crossed(&d1v).add(&d1u.crossed(&d2v));
            Some((p, norm, dnu, dnv))
        }
    }
}

/// The `Vcalc` apex adjustment of `DerivAndNorm` (`cxx:187-201`) and
/// `NormAndDn` (`cxx:296-310`): step one unit along V when the requested `V`
/// sits on the apex within `RealEpsilon()`, so the `Rad + Vcalc * Sina` sign
/// picks the valid branch.
fn cone_vcalc(s: &SurfaceAdaptor, v: f64, sina: f64, rad: f64) -> f64 {
    if (v * sina + rad).abs() <= REAL_EPSILON {
        let vfi = s.first_v_parameter();
        if vfi < -rad / sina {
            v - 1.0
        } else {
            v + 1.0
        }
    } else {
        v
    }
}
