//! `ChFiKPart_ComputeData_FilPlnCon` — the fillet between a plane and a cone,
//! guideline a circle: a torus, or a sphere when the two radii coincide.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_FilPlnCon.cxx:48-341`
//! (OCCT 8.0.0). This is both cone arms of `ChFiKPart_ComputeData.cxx:177-207`;
//! like OCCT there is only the `GeomAbs_Circle` overload — the plane/cone
//! dispatcher arm does not test `ctyp`.
//!
//! `plandab` is `true` when the plane is `S1` (`ChFiKPart_ComputeData.cxx:181`).

use std::sync::Arc;

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{
    GpAx22d, GpAx3, GpCirc, GpCirc2d, GpCone, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d, GpSphere,
    GpTorus, GpVec, GpVec2d,
};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::intana::{quadric_quadric_plane_cone, QuadricIntersection};
use occt_geom::sphere::GeomSphere;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom::torus::GeomTorus;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChFiKPart_MakeFillet(DStr, Data, Pln, Con, fu, lu, Or1, Or2, Radius, Spine,
/// First, Ofpl, plandab)` (`cxx:48-341`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_fillet_plane_cone(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cone: &GpCone,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> Result<bool, String> {
    let mut c1sphere = false;

    // `gp_Ax3 PosPl = Pln.Position(); gp_Dir Dpnat =
    // PosPl.XDirection().Crossed(PosPl.YDirection()); gp_Dir Dp = Dpnat; gp_Dir
    // Df = Dp; if (Or1 == REVERSED) Dp.Reverse(); if (Ofpl == REVERSED)
    // Df.Reverse();` (`cxx:64-75`).
    let pos_pl = pln.position();
    let dpnat = pos_pl
        .x_direction()
        .crossed(pos_pl.y_direction())
        .map_err(|e| e.to_string())?;
    let mut dp = dpnat;
    let mut df = dp;
    if or1 == Orientation::Reversed {
        dp.reverse();
    }
    if ofpl == Orientation::Reversed {
        df.reverse();
    }

    // `gp_Pnt Or = Con.Location(); double u, v; ElSLib::PlaneParameters(PosPl,
    // Or, u, v); gp_Pnt2d c2dPln(u, v); ElSLib::PlaneD0(u, v, PosPl, Or);
    // gp_Pnt cPln = Or; Or.SetCoord(Or.X() + Radius * Dp.X(), ...);`
    // (`cxx:77-83`).
    let or_base = cone.location();
    let (c2d_u, c2d_v) = slib::plane_parameters(&pos_pl, &or_base);
    let c2d_pln = GpPnt2d::new(c2d_u, c2d_v);
    let c_pln = slib::plane_value(pln, c2d_u, c2d_v);
    let or = shifted(&or_base, radius, &dp);

    // `gp_Pnt PtSp; gp_Vec DSp; ElCLib::D1(First, Spine, PtSp, DSp);` 
    // (`cxx:85-87`).
    let (pt_sp, d_sp) = clib::circle_d1(spine, first);

    // `IntAna_QuadQuadGeo CInt(Pln, Con, Precision::Angular(),
    // Precision::Confusion()); gp_Pnt Pv; if (CInt.IsDone()) { Pv =
    // ElCLib::Value(ElCLib::Parameter(CInt.Circle(1), PtSp),
    // CInt.Circle(1)); } else { return false; }` (`cxx:88-99`).
    let pv = match quadric_quadric_plane_cone(pln, cone, occt_core::precision::ANGULAR, CONFUSION) {
        QuadricIntersection::Circle(c) => {
            let par = clib::circle_parameter(&c.position(), &pt_sp);
            clib::circle_value(&c, par)
        }
        QuadricIntersection::None => return Ok(false),
        // `CInt.Circle(1)` on a non-circle answer raises
        // `Standard_OutOfRange` in OCCT.
        _ => return Err("ChFiKPart_MakeFillet: CInt.Circle(1)".into()),
    };

    // `gp_Dir Dx(gp_Vec(cPln, Pv)); gp_Dir Dy(DSp); ElSLib::Parameters(Con, Pv,
    // u, v); gp_Pnt PtCon; gp_Vec Vu, Vv; ElSLib::D1(u, v, Con, PtCon, Vu, Vv);
    // gp_Dir Dc(Vu.Crossed(Vv)); if (Or2 == REVERSED) { Dc.Reverse(); }`
    // (`cxx:100-110`).
    let dx = GpDir::from_vec(&GpVec::from_pnts(&c_pln, &pv)).map_err(|e| e.to_string())?;
    let dy = GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?;
    let (u, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &pv);
    let (_pt_con, vu, vv) = surface_eval::cone_d1(cone, u, v);
    let mut dc = GpDir::from_vec(&vu.crossed(&vv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dc.reverse();
    }

    // `gp_Dir Dz = Dp; gp_Pnt pp(Pv.X() + Dc.X(), ...); ElSLib::PlaneParameters
    // (PosPl, pp, u, v); ElSLib::PlaneD0(u, v, PosPl, pp); gp_Dir ddp(gp_Vec
    // (Pv, pp));` (`cxx:111-116`). `PlaneD0` recomputes `pp` from the (u, v)
    // just read, i.e. the same point up to round-off — the port round-trips
    // through `plane_parameters`/`plane_value` for the same reason.
    let mut dz = dp;
    let pp_raw = shifted(&pv, 1.0, &dc);
    let (u, v) = slib::plane_parameters(&pos_pl, &pp_raw);
    let pp = slib::plane_value(pln, u, v);
    let ddp = GpDir::from_vec(&GpVec::from_pnts(&pv, &pp)).map_err(|e| e.to_string())?;

    // `ElSLib::Parameters(Con, Pv, u, v); gp_Vec dcu, dcv; ElSLib::D1(u, v,
    // Con, pp, dcu, dcv); gp_Dir ddc(dcv); if (ddc.Dot(Dp) < 0.) {
    // ddc.Reverse(); }` (`cxx:117-124`).
    let (u, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &pv);
    let (_pp, _dcu, dcv) = surface_eval::cone_d1(cone, u, v);
    let mut ddc = GpDir::from_vec(&dcv).map_err(|e| e.to_string())?;
    if ddc.dot(&dp) < 0.0 {
        ddc.reverse();
    }

    // `double Ang = ddp.Angle(ddc); double Rabio = Radius / std::tan(Ang / 2);
    // double Maxrad = cPln.Distance(Pv);` (`cxx:125-127`).
    let ang = ddp.angle(&ddc);
    let rabio = radius / (ang / 2.0).tan();
    let maxrad = c_pln.distance(&pv);

    // `double Rad; bool dedans = Dx.Dot(Dc) <= 0.; if (dedans) { if (!plandab)
    // { Dz.Reverse(); } Rad = Maxrad - Rabio; if (std::abs(Rad) <=
    // Precision::Confusion()) { c1sphere = true; } else if (Rad < 0) { return
    // false; } } else { if (plandab) { Dz.Reverse(); } Rad = Maxrad + Rabio;
    // }` (`cxx:128-156`). The `OCCT_DEBUG` print of the `Rad < 0` arm is not
    // translated.
    let dedans = dx.dot(&dc) <= 0.0;
    let rad;
    if dedans {
        if !plandab {
            dz.reverse();
        }
        rad = maxrad - rabio;
        if rad.abs() <= CONFUSION {
            c1sphere = true;
        } else if rad < 0.0 {
            return Ok(false);
        }
    } else {
        if plandab {
            dz.reverse();
        }
        rad = maxrad + rabio;
    }

    // `gp_Ax3 FilAx3(Or, Dz, Dx); if (FilAx3.YDirection().Dot(Dy) <= 0.) {
    // FilAx3.YReverse(); }` (`cxx:157-161`).
    let mut fil_ax3 = GpAx3::from_point_dir_x(or, dz, &dx).map_err(|e| e.to_string())?;
    if fil_ax3.y_direction().dot(&dy) <= 0.0 {
        fil_ax3.y_reverse();
    }

    // `if (c1sphere) { Data->ChangeSurf(IndexSurfaceInDS(new
    // Geom_SphericalSurface(FilAx3, Radius))); } else {
    // Data->ChangeSurf(IndexSurfaceInDS(new Geom_ToroidalSurface(FilAx3, Rad,
    // Radius))); }` (`cxx:163-172`).
    let fillet_surf: Arc<dyn GeomSurface> = if c1sphere {
        let s = GpSphere::new(fil_ax3.clone(), radius).map_err(|e| e.to_string())?;
        Arc::new(GeomSphere::new(s))
    } else {
        let t = GpTorus::new(fil_ax3.clone(), rad, radius).map_err(|e| e.to_string())?;
        Arc::new(GeomTorus::new(t))
    };
    data.set_surf(ch_fi_kpart_index_surface_in_ds(fillet_surf, ds));

    // `gp_Pnt P, PP; gp_Vec deru, derv; P.SetCoord(cPln.X() + Rad * Dx.X(),
    // ...); if (c1sphere) { ElSLib::SphereParameters(FilAx3, Rad, P, u, v);
    // ElSLib::SphereD1(u, v, FilAx3, Rad, PP, deru, derv); } else {
    // ElSLib::TorusParameters(FilAx3, Rad, Radius, P, u, v); ElSLib::TorusD1(u,
    // v, FilAx3, Rad, Radius, PP, deru, derv); if (!plandab && Ang < M_PI / 2 &&
    // dedans) { v = v + 2 * M_PI; } }` (`cxx:176-192`).
    let p_on_fillet = shifted(&c_pln, rad, &dx);
    let (u_fil, mut v_fil) = if c1sphere {
        slib::sphere_parameters(&fil_ax3, &p_on_fillet)
    } else {
        slib::torus_parameters(&fil_ax3, rad, radius, &p_on_fillet)
    };
    if !c1sphere && !plandab && ang < std::f64::consts::FRAC_PI_2 && dedans {
        v_fil += 2.0 * std::f64::consts::PI;
    }

    // `gp_Pnt2d p2dFil(0., v); gp_Dir norFil(deru.Crossed(derv)); bool
    // toreverse = (norFil.Dot(Df) <= 0.);` (`cxx:193-203`).
    let mut p2d_fil = GpPnt2d::new(0.0, v_fil);
    let nor_fil = fillet_normal(fil_ax3.clone(), rad, radius, c1sphere, u_fil, v_fil)?;
    let toreverse = nor_fil.dot(&df) <= 0.0;
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // Planar face (`cxx:208-251`).
    let mut circ_ax2 = fil_ax3.ax2();
    let mut circ_pln: Option<GpCirc> = None;
    let mut circ2d_pln: Option<GpCirc2d> = None;
    if !c1sphere {
        // `ElSLib::PlaneParameters(PosPl, P, u, v); gp_Pnt2d p2dPln(u, v);
        // gp_Dir2d d2d(DSp.Dot(PosPl.XDirection()), DSp.Dot(PosPl.YDirection()));
        // gp_Ax22d ax2dPln(c2dPln, gp_Dir2d(gp_Vec2d(c2dPln, p2dPln)), d2d);
        // gp_Circ2d circ2dPln(ax2dPln, Rad); circAx2.SetLocation(cPln); gp_Circ
        // circPln(circAx2, Rad);` (`cxx:216-224`).
        let (u, v) = slib::plane_parameters(&pos_pl, &p_on_fillet);
        let p2d_pln = GpPnt2d::new(u, v);
        let d2d = GpDir2d::new(
            d_sp.dot(&GpVec::from_xyz(pos_pl.x_direction().xyz())),
            d_sp.dot(&GpVec::from_xyz(pos_pl.y_direction().xyz())),
        )
        .map_err(|e| e.to_string())?;
        let v2d = GpVec2d::new(p2d_pln.x() - c2d_pln.x(), p2d_pln.y() - c2d_pln.y());
        let vx = GpDir2d::new(v2d.x(), v2d.y()).map_err(|e| e.to_string())?;
        // `gp_Ax22d(P, Vx, Vy)` stores both directions without re-deriving
        // them (`gp_Ax22d.hxx:60-66`); `GpAx22d::new` would validate
        // perpendicularity instead, so the fields are set directly.
        let ax2d_pln = GpAx22d {
            point: c2d_pln,
            vxdir: vx,
            vydir: d2d,
        };
        circ2d_pln = Some(GpCirc2d::new(ax2d_pln, rad));
        circ_ax2.set_location(c_pln);
        circ_pln = Some(GpCirc::new(circ_ax2, rad));
    } else {
        // `ElSLib::PlaneParameters(PosPl, P, u, v); gp_Pnt2d p2dPln(u, v),
        // pbid; if (plandab) { Data->Set2dPoints(p2dPln, p2dPln, pbid, pbid); }
        // else { Data->Set2dPoints(pbid, pbid, p2dPln, p2dPln); }` — the sphere
        // arm of this overload has no `Set2dPoints` call (`cxx:214-225`), so
        // the 2D points stay untouched.
        let _ = slib::plane_parameters(&pos_pl, &p_on_fillet);
    }

    // `gp_Lin2d lin2dFil(p2dFil, gp::DX2d()); GLin2dFil1 = new
    // Geom2d_Line(lin2dFil); toreverse = (norFil.Dot(Dpnat) <= 0.); if
    // ((toreverse && plandab) || (!toreverse && !plandab)) { trans =
    // TopAbs_FORWARD; } else { trans = TopAbs_REVERSED; }` (`cxx:226-237`).
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let lin2d_fil = Geom2dLine::from_pnt_dir(p2d_fil, dx2d);
    let toreverse = nor_fil.dot(&dpnat) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let circ_3d: Option<Arc<dyn GeomCurve>> =
        circ_pln.map(|c| Arc::new(GeomCircle::new(c)) as Arc<dyn GeomCurve>);
    let ic = ch_fi_kpart_index_curve_in_ds(circ_3d, ds);
    let on_face: Option<Arc<dyn Curve2d>> =
        circ2d_pln.map(|c| Arc::new(Geom2dCircle::new(c)) as Arc<dyn Curve2d>);
    let on_fil: Arc<dyn Curve2d> = Arc::new(lin2d_fil);
    if plandab {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    } else {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    }

    // Conic face (`cxx:253-339`).
    // `P.SetCoord(Pv.X() + Rabio * ddc.X(), ...); if (c1sphere) {
    // ElSLib::SphereParameters(FilAx3, Radius, P, u, v); ElSLib::SphereD1(u, v,
    // FilAx3, Radius, PP, deru, derv); } else { ElSLib::TorusParameters
    // (FilAx3, Rad, Radius, P, u, v); ElSLib::TorusD1(u, v, FilAx3, Rad,
    // Radius, PP, deru, derv); if (plandab && Ang < M_PI / 2 && dedans) { v =
    // v + 2 * M_PI; } } norFil = deru.Crossed(derv);`
    // (`cxx:256-271`). The sphere arm here uses `Radius` (`:259-260`), not the
    // `Rad` of the planar face's arm.
    let p_con = shifted(&pv, rabio, &ddc);
    let (u_con, mut v_con) = if c1sphere {
        slib::sphere_parameters(&fil_ax3, &p_con)
    } else {
        slib::torus_parameters(&fil_ax3, rad, radius, &p_con)
    };
    if !c1sphere && plandab && ang < std::f64::consts::FRAC_PI_2 && dedans {
        v_con += 2.0 * std::f64::consts::PI;
    }
    let nor_fil2 = fillet_normal(fil_ax3.clone(), rad, radius, c1sphere, u_con, v_con)?;

    // `p2dFil.SetCoord(0., v); lin2dFil.SetLocation(p2dFil); GLin2dFil2 = new
    // Geom2d_Line(lin2dFil);` (`cxx:272-274`).
    p2d_fil = GpPnt2d::new(0.0, v_con);

    // `ElSLib::Parameters(Con, P, u, v); double tol = Precision::PConfusion();
    // bool careaboutsens = false; if (std::abs(lu - fu - 2 * M_PI) < tol) {
    // careaboutsens = true; } if (u >= fu - tol && u < fu) { u = fu; } if (u <=
    // lu + tol && u > lu) { u = lu; } if (u < fu || u > lu) { u =
    // ElCLib::InPeriod(u, fu, fu + 2 * M_PI); }` (`cxx:275-293`). The period
    // call here is the three-argument `ElCLib::InPeriod`, unlike the
    // `ChFiKPart_InPeriod` used in the cylinder overload.
    let (mut u, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &p_con);
    let tol = PCONFUSION;
    let careaboutsens = (lu - fu - 2.0 * std::f64::consts::PI).abs() < tol;
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    if u < fu || u > lu {
        u = clib::in_period(u, fu, fu + 2.0 * std::f64::consts::PI);
    }

    // `ElSLib::D1(u, v, Con, PP, deru, derv); gp_Dir norCon =
    // deru.Crossed(derv); gp_Dir2d d2dCon = gp::DX2d(); if (deru.Dot(Dy) < 0.)
    // { d2dCon.Reverse(); if (careaboutsens && std::abs(fu - u) < tol) { u =
    // lu; } } else if (careaboutsens && std::abs(lu - u) < tol) { u = fu; }`
    // (`cxx:294-308`).
    let (_pp, deru, derv) = surface_eval::cone_d1(cone, u, v);
    let nor_con = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    let mut d2d_con = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    if deru.dot(&GpVec::from_xyz(dy.xyz())) < 0.0 {
        d2d_con.reverse();
        if careaboutsens && (fu - u).abs() < tol {
            u = lu;
        }
    } else if careaboutsens && (lu - u).abs() < tol {
        u = fu;
    }

    // `double scal = gp_Vec(Dp).Dot(gp_Vec(Pv, P)); PP.SetCoord(cPln.X() + scal
    // * Dp.X(), ...); circAx2.SetLocation(PP); gp_Circ circCon(circAx2,
    // P.Distance(PP));` (`cxx:312-315`).
    let scal = GpVec::from_xyz(dp.xyz()).dot(&GpVec::from_pnts(&pv, &p_con));
    let pp_axis = shifted(&c_pln, scal, &dp);
    circ_ax2.set_location(pp_axis);
    let circ_con = GpCirc::new(circ_ax2, p_con.distance(&pp_axis));

    // `toreverse = (norFil.Dot(norCon) <= 0.); if ((toreverse && plandab) ||
    // (!toreverse && !plandab)) { trans = TopAbs_REVERSED; } else { trans =
    // TopAbs_FORWARD; }` (`cxx:317-325`).
    let toreverse = nor_fil2.dot(&nor_con) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };
    let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_con));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
    let on_face: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(u, v),
        d2d_con,
    ));
    let on_fil: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d_fil, dx2d));
    if plandab {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    } else {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    }

    Ok(true)
}

/// `P.SetCoord(base.X() + d * dir.X(), ...)`, the four occurrences of that
/// pattern in `cxx:83`, `:113`, `:178`, `:256`.
fn shifted(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// `ElSLib::SphereD1` / `ElSLib::TorusD1` followed by `gp_Dir(deru.Crossed
/// (derv))` (`cxx:182-187`, `:260-265`).
fn fillet_normal(
    ax3: GpAx3,
    rad: f64,
    radius: f64,
    c1sphere: bool,
    u: f64,
    v: f64,
) -> Result<GpDir, String> {
    let (deru, derv) = if c1sphere {
        let s = GpSphere::new(ax3, radius).map_err(|e| e.to_string())?;
        let (_, du, dv) = surface_eval::sphere_d1(&s, u, v);
        (du, dv)
    } else {
        let t = GpTorus::new(ax3, rad, radius).map_err(|e| e.to_string())?;
        let (_, du, dv) = surface_eval::torus_d1(&t, u, v);
        (du, dv)
    };
    GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())
}
