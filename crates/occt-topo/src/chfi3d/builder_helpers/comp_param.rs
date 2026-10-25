//! `CompParam` — the static helper of `ChFi3d_Builder_6.cxx` that re-expresses
//! the parameters of a spine point on the fillet surface and on the p-curve of
//! the face.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:209-311`
//! (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - the `Carc` argument is `BRepAdaptor_Curve2d`, which OCCT passes to
//!   `CompParam`'s `const Geom2dAdaptor_Curve&` because it derives from it
//!   (`BRepAdaptor_Curve2d.hxx:41`); the port passes [`Curve2dAdaptor`], whose
//!   `first_parameter` / `last_parameter` are the same adaptor range,
//! - `Geom2dAPI_ProjectPointOnCurve(P, C)` (`cxx:278`) is
//!   [`occt_geom2d::curve_ops::curve2d_closest_point`]; `NbPoints() == 0` is
//!   its `None`,
//! - `Geom2dInt_GInter::Perform(C1, C2, TolConf, Tol)` (`cxx:237-240`) is
//!   [`Geom2dIntGInter::perform`] with the two adaptor domains of
//!   [`Geom2dIntGInter::domain_over_range`],
//! - an unloaded `Carc` raises `Standard_DomainError` from
//!   `Geom2dAdaptor_Curve::Value` in OCCT; the port reports that as `Err`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::{PCONFUSION, PINTERSECTION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::curve_ops::curve2d_closest_point;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;

/// `CompParam(Carc, Ctg, Parc, Ptg, Prefarc, Preftg)` (`cxx:209-311`).
///
/// `Parc` is the parameter on `Carc` (the arc's p-curve on the face), `Ptg` the
/// parameter on `Ctg` (the approximated p-curve of the fillet surface); the two
/// are returned through the `&mut` pairs.
pub fn comp_param(
    carc: &Curve2dAdaptor,
    ctg: &dyn Curve2d,
    parc: &mut f64,
    ptg: &mut f64,
    prefarc: f64,
    preftg: f64,
) -> Result<(), String> {
    let carc_curve = carc
        .curve()
        .ok_or_else(|| "Geom2dAdaptor_Curve::Value: Standard_DomainError".to_string())?;
    let point: GpPnt2d = carc_curve.d0(prefarc);
    let distini = point.distance(&ctg.d0(preftg));
    if distini <= PCONFUSION {
        // (1) The provided parameters are good (`cxx:218-226`).
        *parc = prefarc;
        *ptg = preftg;
        return Ok(());
    }

    // (2) Intersection (`cxx:228-268`).
    let tol = PINTERSECTION;
    let d_ctg = Geom2dIntGInter::domain_over_range(
        ctg,
        ctg.first_parameter(),
        ctg.last_parameter(),
        tol.max(tol),
    );
    let d_carc = Geom2dIntGInter::domain_over_range(
        carc_curve.as_ref(),
        carc.first_parameter(),
        carc.last_parameter(),
        tol.max(tol),
    );
    let mut inter = Geom2dIntGInter::new();
    inter.perform(ctg, &d_ctg, carc_curve.as_ref(), &d_carc, tol, tol);

    let mut found = false;
    let mut dist = f64::INFINITY;
    if inter.is_done() && inter.nb_points() > 0 {
        for i in 1..=inter.nb_points() {
            let int2d = inter.point(i);
            let p1 = int2d.param_on_first();
            let p2 = int2d.param_on_second();
            if (prefarc - p2).abs() < dist {
                *ptg = p1;
                *parc = p2;
                dist = (prefarc - p2).abs();
                found = true;
            }
        }
    }

    if !found {
        // (3) Projection (`cxx:270-310`).
        *parc = prefarc;
        match curve2d_closest_point(ctg, &point, PCONFUSION) {
            // `NbPoints() == 0`: the extremity is taken (`cxx:282-288`).
            None => *ptg = preftg,
            Some((u, q)) => {
                // `LowerDistance() < distini` (`cxx:292-299`).
                if q.distance(&point) < distini {
                    *ptg = u;
                } else {
                    *ptg = preftg;
                }
            }
        }
    }
    Ok(())
}
