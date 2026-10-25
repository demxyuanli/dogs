//! `ChFiKPart_ComputeData_Sphere` — the spherical corner of a fillet.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_Sphere.cxx:47-206`
//! (OCCT 8.0.0), `bool ChFiKPart_Sphere(DStr, Data, S1, S2, OrFace1, OrFace2,
//! Or1, Or2, Rad, PS1, P1S2, P2S2)`.
//!
//! The sphere is centred on the line normal to the circle through the three
//! contour points, at the distance `sqrt(Rad^2 - rr^2)` from that circle's
//! centre (`cxx:84-101`); its axis has to point from the centre to the `S1`
//! point, which is the test of `cxx:89-101`. `PS1` is the `S1` parameter pair,
//! `P1S2` / `P2S2` the two `S2` pairs.
//!
//! The p-curve on the spherical surface of the `S2` side is the projection of
//! the `[0, ang]` arc of the corner circle ([`ch_fi_kpart_proj_pc`]); that
//! projection is analytic only when the circle is a meridian or a parallel of
//! the sphere, and otherwise OCCT falls back to `ProjLib_ComputeApprox`
//! (`ProjLib_ProjectedCurve.cxx:717-785`), which is **UNPORTED**. On the
//! generic corner this function therefore answers `Err` with that marker.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{
    GpAx2, GpAx3, GpCirc, GpDir, GpDir2d, GpPnt, GpPnt2d, GpSphere, GpTrsf2d, GpVec, GpVec2d,
};
use occt_core::precision::CONFUSION;
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::gce_make_circ::gce_make_circ_3p;
use occt_geom::sphere::GeomSphere;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::fcts::{
    ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds, ch_fi_kpart_p_curve,
    ch_fi_kpart_proj_pc,
};

/// `gp::DX2d()` (`gp.hxx`) — the 2D X axis direction.
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };

/// `ChFiKPart_Sphere(DStr, Data, S1, S2, OrFace1, OrFace2, Or1, Or2, Rad, PS1,
/// P1S2, P2S2)` (`cxx:47-206`).
///
/// `OrFace2` (`cxx:52`) and the last orientation (`cxx:54`) are unnamed in
/// OCCT and never read; the port keeps them as `_or_face2` / `_or2`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_sphere(
    ds: &mut DataStructure,
    data: &mut SurfData,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    or_face1: Orientation,
    _or_face2: Orientation,
    or1: Orientation,
    _or2: Orientation,
    rad: f64,
    ps1: GpPnt2d,
    p1s2: GpPnt2d,
    p2s2: GpPnt2d,
) -> Result<bool, String> {
    // `double ptol = Precision::Confusion(); gp_Pnt p1, p2, p3; gp_Vec v1, v2;
    // S1->D1(PS1.X(), PS1.Y(), p1, v1, v2);` (`cxx:64-68`).
    let ptol = CONFUSION;
    let (p1, v1, v2) = s1
        .d1(ps1.x(), ps1.y())
        .ok_or_else(|| "ChFiKPart_Sphere: S1->D1".to_string())?;

    // `gp_Dir ds1(v1.Crossed(v2)); gp_Dir df1 = ds1; gp_Dir dnat1 = ds1; if
    // (Or1 == TopAbs_REVERSED) { ds1.Reverse(); } if (OrFace1 ==
    // TopAbs_REVERSED) { df1.Reverse(); }` (`cxx:69-79`).
    let mut ds1 = GpDir::from_vec(&v1.crossed(&v2)).map_err(|e| e.to_string())?;
    let dnat1 = ds1;
    let mut df1 = ds1;
    if or1 == Orientation::Reversed {
        ds1.reverse();
    }
    if or_face1 == Orientation::Reversed {
        df1.reverse();
    }

    // `S2->D0(P1S2.X(), P1S2.Y(), p2); S2->D0(P2S2.X(), P2S2.Y(), p3);`
    // (`cxx:80-81`).
    let p2 = s2
        .d0(p1s2.x(), p1s2.y())
        .ok_or_else(|| "ChFiKPart_Sphere: S2->D0 (P1S2)".to_string())?;
    let p3 = s2
        .d0(p2s2.x(), p2s2.y())
        .ok_or_else(|| "ChFiKPart_Sphere: S2->D0 (P2S2)".to_string())?;

    // `gp_Circ ci = gce_MakeCirc(p1, p2, p3); gp_Dir di = ci.Axis().Direction();
    // gp_Pnt pp = ci.Location(); double rr = ci.Radius();` (`cxx:82-86`).
    let ci = gce_make_circ_3p(&p1, &p2, &p3)?;
    let ci_ax = ci.axis();
    let di = *ci_ax.direction();
    let pp = ci.location();
    let rr = ci.radius();

    // `double delta = sqrt(Rad * Rad - rr * rr); gp_Pnt cen(pp.X() + delta *
    // di.X(), ...); gp_Dir dz(gp_Vec(p1, cen)); if (std::abs(ds1.Dot(dz) - 1.)
    // > ptol) { cen.SetCoord(pp.X() - delta * di.X(), ...); dz =
    // gp_Dir(gp_Vec(p1, cen)); if (std::abs(ds1.Dot(dz) - 1.) > ptol) { return
    // false; } }` (`cxx:87-101`). The `OCCT_DEBUG` print at `cxx:96-98` is not
    // translated.
    let delta = (rad * rad - rr * rr).sqrt();
    let mut cen = GpPnt::new(
        pp.x() + delta * di.x(),
        pp.y() + delta * di.y(),
        pp.z() + delta * di.z(),
    );
    let mut dz = GpDir::from_vec(&GpVec::from_pnts(&p1, &cen)).map_err(|e| e.to_string())?;
    if (ds1.dot(&dz) - 1.0).abs() > ptol {
        cen = GpPnt::new(
            pp.x() - delta * di.x(),
            pp.y() - delta * di.y(),
            pp.z() - delta * di.z(),
        );
        dz = GpDir::from_vec(&GpVec::from_pnts(&p1, &cen)).map_err(|e| e.to_string())?;
        if (ds1.dot(&dz) - 1.0).abs() > ptol {
            return Ok(false);
        }
    }

    // `gp_Dir ddx(gp_Vec(cen, p2)); gp_Dir dddx = ddx; gp_Dir ddy(gp_Vec(cen,
    // p3)); gp_Dir dx = dz.Crossed(ddx.Crossed(dz)); gp_Dir ddz =
    // dz.Reversed();` (`cxx:102-107`).
    let ddx = GpDir::from_vec(&GpVec::from_pnts(&cen, &p2)).map_err(|e| e.to_string())?;
    let mut dddx = ddx;
    let ddy = GpDir::from_vec(&GpVec::from_pnts(&cen, &p3)).map_err(|e| e.to_string())?;
    let dx = dz
        .crossed(&ddx.crossed(&dz).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let mut ddz = dz.reversed();

    // `gp_Ax3 FilAx3(cen, dz, dx); if (FilAx3.YDirection().Dot(ddy) <= 0.) {
    // FilAx3.YReverse(); ddz.Reverse(); dddx.Reverse(); }` (`cxx:108-112`).
    let mut fil_ax3 = GpAx3::from_point_dir_x(cen, dz, &dx).map_err(|e| e.to_string())?;
    if fil_ax3.y_direction().dot(&ddy) <= 0.0 {
        fil_ax3.y_reverse();
        ddz.reverse();
        dddx.reverse();
    }

    // `occ::handle<Geom_SphericalSurface> gsph = new
    // Geom_SphericalSurface(FilAx3, Rad); Data->ChangeSurf(
    // ChFiKPart_IndexSurfaceInDS(gsph, DStr));` (`cxx:113-114`).
    let gsph = GpSphere::new(fil_ax3.clone(), rad).map_err(|e| e.to_string())?;
    let gsph_handle: Arc<dyn GeomSurface> = Arc::new(GeomSphere::new(gsph));
    data.set_surf(ch_fi_kpart_index_surface_in_ds(gsph_handle.clone(), ds));

    // `bool toreverse = (ddz.Dot(df1) <= 0.); if (toreverse) {
    // Data->ChangeOrientation() = TopAbs_REVERSED; } else { ...
    // TopAbs_FORWARD; }` (`cxx:118-126`).
    let toreverse = ddz.dot(&df1) <= 0.0;
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // `double uu1, vv1, uu2, vv2; ElSLib::SphereParameters(FilAx3, Rad, p2, uu1,
    // vv1); uu1 = 0.; ElSLib::SphereParameters(FilAx3, Rad, p3, uu2, vv2);`
    // (`cxx:130-134`).
    let (_u1, vv1) = slib::sphere_parameters(&fil_ax3, &p2);
    let uu1 = 0.0;
    let (uu2, _vv2) = slib::sphere_parameters(&fil_ax3, &p3);

    // Pointed side (`cxx:138-152`).
    // `gp_Pnt2d p2dFil(0., -M_PI / 2.); gp_Lin2d lin2dFil(p2dFil, gp::DX2d());
    // occ::handle<Geom2d_Curve> C2dFil = new Geom2d_Line(lin2dFil);`
    // (`cxx:142-145`).
    let p2dfil = GpPnt2d::new(0.0, -std::f64::consts::FRAC_PI_2);
    let c2dfil_line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2dfil, DX2D));

    // `toreverse = (ddz.Dot(dnat1) <= 0.); TopAbs_Orientation trans =
    // TopAbs_REVERSED; if (toreverse) { trans = TopAbs_FORWARD; }
    // Data->ChangeInterferenceOnS1().SetInterference(
    // ChFiKPart_IndexCurveInDS(C, DStr), trans, C2d, C2dFil);` (`cxx:146-152`).
    // `C` and `C2d` are still null there (`cxx:139-140`).
    let toreverse = ddz.dot(&dnat1) <= 0.0;
    let trans = if toreverse {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let ic = ch_fi_kpart_index_curve_in_ds(None, ds);
    data.interference_on_s1_mut()
        .set_interference(ic, trans, None, Some(c2dfil_line.clone()));

    // The other side (`cxx:156-181`).
    // `double ang = ddx.Angle(ddy); gp_Dir dci = ddx.Crossed(ddy); gp_Ax2
    // axci(cen, dci, ddx); gp_Circ ci2(axci, Rad); C = new Geom_Circle(ci2);`
    // (`cxx:157-161`).
    let ang = ddx.angle_tol(&ddy);
    let dci = ddx.crossed(&ddy).map_err(|e| e.to_string())?;
    let axci = GpAx2::new(cen, dci, ddx).map_err(|e| e.to_string())?;
    let ci2 = GpCirc::new(axci, rad);
    let circle3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(ci2));

    // `GeomAdaptor_Surface AS(gsph); GeomAdaptor_Curve AC(C, 0., ang);
    // ChFiKPart_ProjPC(AC, AS, C2dFil);` (`cxx:162-165`).
    let arc: Arc<dyn GeomCurve> =
        Arc::new(GeomTrimmedCurveBasis::new(circle3d.clone(), 0.0, ang));
    let c2dfil = ch_fi_kpart_proj_pc(gsph_handle, arc)?;

    // `gp_Pnt2d p2dbid = C2dFil->Value(0.); gp_Pnt2d pp2dbid(uu1, vv1); if
    // (!pp2dbid.IsEqual(p2dbid, ptol)) { gp_Vec2d v2dbid(p2dbid, pp2dbid);
    // C2dFil->Translate(v2dbid); }` (`cxx:166-171`).
    let p2dbid = c2dfil.d0(0.0);
    let pp2dbid = GpPnt2d::new(uu1, vv1);
    let c2dfil: Arc<dyn Curve2d> = if !pp2dbid.is_equal(&p2dbid, ptol) {
        let v2dbid = GpVec2d::new(pp2dbid.x() - p2dbid.x(), pp2dbid.y() - p2dbid.y());
        let mut trsf = GpTrsf2d::identity();
        trsf.set_translation_vec(&v2dbid);
        Arc::from(c2dfil.transformed(&trsf))
    } else {
        c2dfil
    };

    // `gp_Vec2d v2d(P1S2, P2S2); gp_Dir2d d2d(v2d); if
    // (std::abs(v2d.Magnitude() - ang) <= ptol) { gp_Lin2d l2d(P1S2, d2d); C2d =
    // new Geom2d_Line(l2d); } else { C2d = ChFiKPart_PCurve(P1S2, P2S2, 0.,
    // ang); }` (`cxx:172-181`).
    let v2d = GpVec2d::new(p2s2.x() - p1s2.x(), p2s2.y() - p1s2.y());
    let c2d: Arc<dyn Curve2d> = if (v2d.magnitude() - ang).abs() <= ptol {
        let d2d = GpDir2d::new(v2d.x(), v2d.y()).map_err(|e| e.to_string())?;
        Arc::new(Geom2dLine::from_pnt_dir(p1s2, d2d))
    } else {
        Arc::new(ch_fi_kpart_p_curve(p1s2, p2s2, 0.0, ang)?)
    };

    // `gp_Pnt pp1; S2->D1(P1S2.X(), P1S2.Y(), pp1, v1, v2); gp_Dir
    // ds2(v1.Crossed(v2)); toreverse = (ds2.Dot(dddx) <= 0.); trans =
    // TopAbs_REVERSED; if (!toreverse) { trans = TopAbs_FORWARD; }
    // Data->ChangeInterferenceOnS2().SetInterference(
    // ChFiKPart_IndexCurveInDS(C, DStr), trans, C2d, C2dFil);` (`cxx:183-193`).
    // `pp1` is not read.
    let (_pp1, v1b, v2b) = s2
        .d1(p1s2.x(), p1s2.y())
        .ok_or_else(|| "ChFiKPart_Sphere: S2->D1".to_string())?;
    let ds2 = GpDir::from_vec(&v1b.crossed(&v2b)).map_err(|e| e.to_string())?;
    let toreverse = ds2.dot(&dddx) <= 0.0;
    let trans = if !toreverse {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circle3d), ds);
    data.interference_on_s2_mut()
        .set_interference(ic, trans, Some(c2d), Some(c2dfil));

    // `Data->ChangeVertexFirstOnS1().SetPoint(p1); ...LastOnS1(p1);
    // ...FirstOnS2(p2); ...LastOnS2(p3);` (`cxx:196-200`).
    data.vertex_first_on_s1_mut().set_point(&p1);
    data.vertex_last_on_s1_mut().set_point(&p1);
    data.vertex_first_on_s2_mut().set_point(&p2);
    data.vertex_last_on_s2_mut().set_point(&p3);

    // `Data->ChangeInterferenceOnS1().SetFirstParameter(0.);
    // SetLastParameter(uu2); ...OnS2().SetFirstParameter(0.);
    // SetLastParameter(ang);` (`cxx:201-204`).
    data.interference_on_s1_mut().set_first_parameter(0.0);
    data.interference_on_s1_mut().set_last_parameter(uu2);
    data.interference_on_s2_mut().set_first_parameter(0.0);
    data.interference_on_s2_mut().set_last_parameter(ang);

    Ok(true)
}
