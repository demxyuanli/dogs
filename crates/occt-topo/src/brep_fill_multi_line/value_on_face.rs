//! `static gp_Pnt2d ValueOnFace(...)` — the solver that maps one parameter of
//! the bisector to an (u, v) point on a face.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_MultiLine.cxx:584-724`.
//!
//! The bisector lives in 2D; its point at `U` is projected on `TheU` (the
//! projected u-iso) to get the first parameter, and the second parameter is
//! read where `TheV` (the projected v-iso) meets the horizontal line at the
//! projected distance. The four steps of the body are reproduced in order:
//! the projection (`cxx:589-603`), the end-point arbitration (`cxx:604-646`),
//! the `TheV` search (`cxx:648-714`) and the iso transposition (`cxx:716-723`).

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::{elib, precision};
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use super::Error;

/// `ValueOnFace(U, TheBis, TheU, TheV, IsIsoU)` (`cxx:584-724`).
///
/// `the_u` / `the_v` are the two projected isos of one face and `is_iso_u`
/// says whether that face's edge is a `U = const` iso, which decides the order
/// of the two answers.
pub fn value_on_face(
    u: f64,
    the_bis: &Geom2dAdaptorCurve,
    the_u: &Geom2dAdaptorCurve,
    the_v: &Geom2dAdaptorCurve,
    is_iso_u: bool,
) -> Result<GpPnt2d, Error> {
    // `gp_Pnt2d P = TheBis.Value(U);` (`cxx:589`).
    let p = the_bis.value(u).ok_or(Error::NullCurve)?;

    // `Geom2dAPI_ProjectPointOnCurve Ext(P, TheU.Curve(),
    // TheU.FirstParameter(), TheU.LastParameter());` (`cxx:591`).
    let cu = the_u.curve().ok_or(Error::NullCurve)?;
    let ext = ProjectPointOnCurve2d::with_curve_range(
        &p,
        cu,
        the_u.first_parameter(),
        the_u.last_parameter(),
    );

    // `constexpr double mult = 5.; constexpr double eps = mult *
    // Precision::Confusion();` (`cxx:593-594`).
    const MULT: f64 = 5.0;
    let eps = MULT * precision::CONFUSION;

    let mut uu = 0.0;
    let mut dist = precision::INFINITE;
    // `cxx:598-603`.
    if ext.nb_points() != 0 {
        uu = ext.lower_distance_parameter().ok_or(Error::ValueOnFace)?;
        dist = ext.lower_distance().ok_or(Error::ValueOnFace)?;
    }

    // `cxx:605-606`: control with the two ends of `TheU`.
    let first_u = the_u.value(the_u.first_parameter()).ok_or(Error::NullCurve)?;
    let last_u = the_u.value(the_u.last_parameter()).ok_or(Error::NullCurve)?;
    let d1 = p.distance(&first_u);
    let d2 = p.distance(&last_u);

    // `cxx:608-645`.
    if d1 < dist || d2 < dist || (d1 - dist).abs() < eps || (d2 - dist).abs() < eps {
        if (d1 - d2).abs() < eps {
            if the_u.get_type() == Some(CurveType::Circle) {
                // `cxx:611-622`: both ends are equally far, so the parameter is
                // read from the angle of the bisector with the circle's X axis.
                let circle = the_u.circle().ok_or(Error::ValueOnFace)?;
                let axis = GpVec2d::from_dir2d(&circle.x_axis().vdir);
                let d12d = the_bis
                    .dn(the_bis.first_parameter(), 1)
                    .ok_or(Error::NullCurve)?;
                let mut ang = axis.angle(&d12d);
                if !circle.is_direct() {
                    ang = -ang;
                }
                uu = elib::clib::in_period(
                    ang,
                    the_u.first_parameter(),
                    the_u.first_parameter() + 2.0 * std::f64::consts::PI,
                );
                dist = circle.radius();
            } else {
                // `cxx:625-630`: `Standard_ConstructionError`.
                return Err(Error::ValueOnFace);
            }
        } else if d1 < d2 {
            dist = d1;
            uu = the_u.first_parameter();
        } else {
            dist = d2;
            uu = the_u.last_parameter();
        }
    }

    // `cxx:648-714`: the second parameter.
    let tol = precision::CONFUSION;
    let pf = the_v.value(the_v.first_parameter()).ok_or(Error::NullCurve)?;
    let pl = the_v.value(the_v.last_parameter()).ok_or(Error::NullCurve)?;

    let vv;
    if (dist - pf.y().abs()).abs() < tol {
        vv = the_v.first_parameter(); // `cxx:650-653`
    } else if (dist - pl.y().abs()).abs() < tol {
        vv = the_v.last_parameter(); // `cxx:654-657`
    } else {
        // `cxx:659-714`.
        if pf.y().min(pl.y()) < -tol {
            dist = -dist;
        }
        let line = Geom2dLine::from_pnt_dir(
            GpPnt2d::new(0.0, dist),
            GpDir2d::new(1.0, 0.0).expect("gp::DX2d"),
        );

        // `Geom2dInt_GInter Intersector(Cu1, Cu2, TolConf, Tol);`
        // (`cxx:672`): `TolConf = 0.` and both domains are
        // `ComputeDomain(curve, max(TolConf, Tol))` over the adaptor's own
        // range — `TheV`'s trimming range and the infinite range of the line.
        let mut inter = Geom2dIntGInter::new();
        let cu1 = the_v.curve().ok_or(Error::NullCurve)?;
        let d1 = Geom2dIntGInter::domain_over_range(
            cu1.as_ref(),
            the_v.first_parameter(),
            the_v.last_parameter(),
            tol,
        );
        let d2 = Geom2dIntGInter::domain_over_range(
            &line,
            line.first_parameter(),
            line.last_parameter(),
            tol,
        );
        inter.perform(cu1.as_ref(), &d1, &line, &d2, 0.0, tol);

        if !inter.is_done() {
            // `cxx:675-682`: OCCT returns `gp_Pnt2d(0., 0.)`.
            return Ok(GpPnt2d::new(0.0, 0.0));
        }
        if inter.nb_points() > 0 {
            vv = inter.point(1).param_on_first(); // `cxx:687-689`
        } else if inter.nb_segments() > 0 {
            // `cxx:690-697`: interpolate inside the overlapping segment.
            let seg = inter.segment(1);
            let vs1 = seg.first_point().param_on_first();
            let vs2 = seg.last_point().param_on_first();
            let ps1 = the_v.value(vs1).ok_or(Error::NullCurve)?;
            let ps2 = the_v.value(vs2).ok_or(Error::NullCurve)?;
            let alp = (dist - ps1.y()) / (ps2.y() - ps1.y());
            vv = alp * (vs2 - vs1) + vs1;
        } else if (dist - pl.y()).abs() < (dist - pf.y()).abs() {
            vv = the_v.last_parameter(); // `cxx:704-711`
        } else {
            vv = the_v.first_parameter();
        }
    }

    // `cxx:716-723`.
    Ok(if is_iso_u {
        GpPnt2d::new(vv, uu)
    } else {
        GpPnt2d::new(uu, vv)
    })
}
