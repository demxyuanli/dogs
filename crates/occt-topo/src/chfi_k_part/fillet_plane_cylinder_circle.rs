//! `ChFiKPart_ComputeData_FilPlnCyl` — the fillet between a plane and a
//! cylinder, guideline a circle: a torus, or a sphere when the two radii
//! coincide.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_FilPlnCyl.cxx:283-599`
//! (OCCT 8.0.0). This is the `GeomAbs_Circle` arm of
//! `ChFiKPart_ComputeData.cxx:124-140` / `:157-173`; the `GeomAbs_Line` arm is
//! [`super::fillet_plane_cylinder`].
//!
//! `plandab` is `true` when the plane is `S1` (`ChFiKPart_ComputeData.cxx:137`).

use std::sync::Arc;

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{
    GpAx22d, GpAx3, GpCirc, GpCirc2d, GpCylinder, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d, GpSphere,
    GpTorus, GpVec, GpVec2d,
};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::sphere::GeomSphere;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom::torus::GeomTorus;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_in_period, ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChFiKPart_MakeFillet(DStr, Data, Pln, Cyl, fu, lu, Or1, Or2, Radius, Spine,
/// First, Ofpl, plandab)` (`cxx:283-599`), `Spine` a `gp_Circ`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_fillet_plane_cylinder_circle(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cyl: &GpCylinder,
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
    // Df.Reverse();` (`cxx:300-311`).
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

    // `gp_Pnt Or = Cyl.Location(); double u, v; ElSLib::PlaneParameters(PosPl,
    // Or, u, v); gp_Pnt2d c2dPln(u, v); ElSLib::PlaneD0(u, v, PosPl, Or);
    // gp_Pnt cPln = Or; Or.SetCoord(Or.X() + Radius * Dp.X(), ...);`
    // (`cxx:313-319`). `ElSLib::PlaneD0` rewrites `Or` with the same point it
    // read, so it is a no-op on exact input.
    let or_base = cyl.location();
    let (c2d_u, c2d_v) = slib::plane_parameters(&pos_pl, &or_base);
    let c2d_pln = GpPnt2d::new(c2d_u, c2d_v);
    let c_pln = slib::plane_value(pln, c2d_u, c2d_v);
    let or = shifted(&or_base, radius, &dp);

    // `gp_Pnt PtSp; gp_Vec DSp; gp_Pnt PtSp2; gp_Vec DSp2; double acote =
    // 1e-7; ElCLib::D1(First, Spine, PtSp, DSp); ElSLib::Parameters(Cyl, PtSp,
    // u, v);` (`cxx:320-327`).
    let (mut pt_sp, mut d_sp) = clib::circle_d1(spine, first);
    let (mut u, v) = slib::cylinder_parameters(&cyl.position(), &pt_sp);

    // `if ((std::abs(u) < acote) || (std::abs(u - (2 * M_PI)) < acote)) {
    // ElCLib::D1(First + 0.2, Spine, PtSp2, DSp2); ElSLib::Parameters(Cyl,
    // PtSp2, u2, v2); if (std::abs(u2 - u) > M_PI) { u = (2 * M_PI) - u; PtSp =
    // ElSLib::Value(u, v, Cyl); double PR = ElCLib::Parameter(Spine, PtSp);
    // ElCLib::D1(PR, Spine, PtSp2, DSp); } }` (`cxx:325-341`). `v` is reused
    // from the first `Parameters` call in the `PtSp` reassignment, `DSp` is
    // overwritten by the third `D1`.
    let acote = 1e-7;
    if u.abs() < acote || (u - 2.0 * std::f64::consts::PI).abs() < acote {
        let (pt_sp2, _d_sp2) = clib::circle_d1(spine, first + 0.2);
        let (u2, _v2) = slib::cylinder_parameters(&cyl.position(), &pt_sp2);
        if (u2 - u).abs() > std::f64::consts::PI {
            u = 2.0 * std::f64::consts::PI - u;
            pt_sp = slib::cylinder_value(cyl, u, v);
            let pr = clib::circle_parameter(&spine.position(), &pt_sp);
            let (_pt, d) = clib::circle_d1(spine, pr);
            d_sp = d;
        }
    }

    // `gp_Dir Dx(gp_Vec(Or, PtSp)); Dx = Dp.Crossed(Dx.Crossed(Dp)); gp_Dir
    // Dy(DSp); gp_Pnt PtCyl; gp_Vec Vu, Vv; ElSLib::D1(u, v, Cyl, PtCyl, Vu,
    // Vv); gp_Dir Dc(Vu.Crossed(Vv)); if (Or2 == REVERSED) Dc.Reverse();`
    // (`cxx:343-353`).
    let mut dx = GpDir::from_vec(&GpVec::from_pnts(&or, &pt_sp)).map_err(|e| e.to_string())?;
    let inner = dx.crossed(&dp).map_err(|e| e.to_string())?;
    dx = dp.crossed(&inner).map_err(|e| e.to_string())?;
    let dy = GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?;
    let (_pt_cyl, vu, vv) = surface_eval::cylinder_d1(cyl, u, v);
    let mut dc = GpDir::from_vec(&vu.crossed(&vv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dc.reverse();
    }

    // `gp_Dir Dz = Dp; double Rad, cylrad = Cyl.Radius(); bool dedans =
    // (Dx.Dot(Dc) <= 0.); if (dedans) { if (!plandab) { Dz.Reverse(); } Rad =
    // cylrad - Radius; if (std::abs(Rad) <= Precision::Confusion()) { c1sphere
    // = true; } else if (Rad < 0) { return false; } } else { if (plandab) {
    // Dz.Reverse(); } Rad = cylrad + Radius; }` (`cxx:354-383`). The
    // `OCCT_DEBUG` print of the `Rad < 0` arm is not translated.
    let mut dz = dp;
    let cylrad = cyl.radius();
    let dedans = dx.dot(&dc) <= 0.0;
    let rad;
    if dedans {
        if !plandab {
            dz.reverse();
        }
        rad = cylrad - radius;
        if rad.abs() <= CONFUSION {
            c1sphere = true;
        } else if rad < 0.0 {
            return Ok(false);
        }
    } else {
        if plandab {
            dz.reverse();
        }
        rad = cylrad + radius;
    }

    // `gp_Ax3 FilAx3(Or, Dz, Dx); if (FilAx3.YDirection().Dot(Dy) <= 0.)
    // FilAx3.YReverse();` (`cxx:384-388`).
    let mut fil_ax3 = GpAx3::from_point_dir_x(or, dz, &dx).map_err(|e| e.to_string())?;
    if fil_ax3.y_direction().dot(&dy) <= 0.0 {
        fil_ax3.y_reverse();
    }

    // `if (c1sphere) { Data->ChangeSurf(IndexSurfaceInDS(new
    // Geom_SphericalSurface(FilAx3, Radius))); } else {
    // Data->ChangeSurf(IndexSurfaceInDS(new Geom_ToroidalSurface(FilAx3, Rad,
    // Radius))); }` (`cxx:390-399`).
    let fillet_surf: Arc<dyn GeomSurface> = if c1sphere {
        let s = GpSphere::new(fil_ax3.clone(), radius).map_err(|e| e.to_string())?;
        Arc::new(GeomSphere::new(s))
    } else {
        let t = GpTorus::new(fil_ax3.clone(), rad, radius).map_err(|e| e.to_string())?;
        Arc::new(GeomTorus::new(t))
    };
    data.set_surf(ch_fi_kpart_index_surface_in_ds(fillet_surf, ds));

    // `gp_Pnt P, PP; gp_Vec deru, derv; P.SetCoord(cPln.X() + Rad * Dx.X(),
    // ...); u = 0.; if ((dedans && plandab) || (!dedans && !plandab)) { if
    // (c1sphere) { v = -M_PI / 2; } else { v = 3 * M_PI / 2; } } else { v =
    // M_PI / 2; }` (`cxx:402-420`).
    let p_fil = shifted(&c_pln, rad, &dx);
    let mut u_fil = 0.0;
    let mut v_fil = if (dedans && plandab) || (!dedans && !plandab) {
        if c1sphere {
            -std::f64::consts::FRAC_PI_2
        } else {
            3.0 * std::f64::consts::FRAC_PI_2
        }
    } else {
        std::f64::consts::FRAC_PI_2
    };

    // `if (c1sphere) { ElSLib::SphereD1(u, v, FilAx3, cylrad, PP, deru, derv);
    // norFil = FilAx3.XDirection().Crossed(FilAx3.YDirection()); if (v < 0.)
    // { norFil.Reverse(); } } else { ElSLib::TorusD1(u, v, FilAx3, Rad, Radius,
    // PP, deru, derv); norFil = gp_Dir(deru.Crossed(derv)); }` (`cxx:421-435`).
    // The sphere arm reads the axis normal rather than the `D1` cross product,
    // so its `SphereD1` call has no observable effect and is not translated.
    let nor_fil = if c1sphere {
        let mut n = fil_ax3
            .x_direction()
            .crossed(fil_ax3.y_direction())
            .map_err(|e| e.to_string())?;
        if v_fil < 0.0 {
            n.reverse();
        }
        n
    } else {
        fillet_normal(fil_ax3.clone(), rad, radius, false, cylrad, u_fil, v_fil)?
    };

    // `gp_Pnt2d p2dFil(0., v); bool toreverse = (norFil.Dot(Df) <= 0.);` 
    // (`cxx:436-445`).
    let mut p2d_fil = GpPnt2d::new(0.0, v_fil);
    let toreverse = nor_fil.dot(&df) <= 0.0;
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // Planar face (`cxx:449-503`).
    // `gp_Ax2 circAx2 = FilAx3.Ax2();` (`cxx:452`).
    let mut circ_ax2 = fil_ax3.ax2();
    let mut circ_pln: Option<GpCirc> = None;
    let mut circ2d_pln: Option<GpCirc2d> = None;
    if !c1sphere {
        // `ElSLib::PlaneParameters(PosPl, P, u, v); gp_Pnt2d p2dPln(u, v);
        // gp_Dir2d d2d(DSp.Dot(PosPl.XDirection()), DSp.Dot(PosPl.YDirection()));
        // gp_Ax22d ax2dPln(c2dPln, gp_Dir2d(gp_Vec2d(c2dPln, p2dPln)), d2d);
        // gp_Circ2d circ2dPln(ax2dPln, Rad);` (`cxx:455-459`).
        let (u, v) = slib::plane_parameters(&pos_pl, &p_fil);
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

        // `circAx2.SetLocation(cPln); gp_Circ circPln(circAx2, Rad);`
        // (`cxx:461-462`).
        circ_ax2.set_location(c_pln);
        circ_pln = Some(GpCirc::new(circ_ax2, rad));
    } else {
        // `ElSLib::PlaneParameters(PosPl, P, u, v); gp_Pnt2d p2dPln(u, v),
        // pbid; if (plandab) { Data->Set2dPoints(p2dPln, p2dPln, pbid, pbid); }
        // else { Data->Set2dPoints(pbid, pbid, p2dPln, p2dPln); }`
        // (`cxx:467-477`).
        let (u, v) = slib::plane_parameters(&pos_pl, &p_fil);
        let p2d_pln = GpPnt2d::new(u, v);
        let pbid = GpPnt2d::new(0.0, 0.0);
        if plandab {
            data.set_2d_points(p2d_pln, p2d_pln, pbid, pbid);
        } else {
            data.set_2d_points(pbid, pbid, p2d_pln, p2d_pln);
        }
    }

    // `gp_Lin2d lin2dFil(p2dFil, gp::DX2d()); GLin2dFil1 = new
    // Geom2d_Line(lin2dFil);` (`cxx:478-479`).
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let lin2d_fil = Geom2dLine::from_pnt_dir(p2d_fil, dx2d);

    // `toreverse = (norFil.Dot(Dpnat) <= 0.); if ((toreverse && plandab) ||
    // (!toreverse && !plandab)) { trans = TopAbs_FORWARD; } else { trans =
    // TopAbs_REVERSED; }` (`cxx:480-489`).
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

    // Cylindrical face (`cxx:505-597`).
    // `P.SetCoord(Or.X() + cylrad * Dx.X(), ...); u = 0.; if (dedans) { if
    // (plandab && !c1sphere) { v = 2 * M_PI; } else { v = 0.; } } else { v =
    // M_PI; } p2dFil.SetCoord(u, v);` (`cxx:506-523`).
    let p_cyl = shifted(&or, cylrad, &dx);
    u_fil = 0.0;
    v_fil = if dedans {
        if plandab && !c1sphere {
            2.0 * std::f64::consts::PI
        } else {
            0.0
        }
    } else {
        std::f64::consts::PI
    };
    p2d_fil = GpPnt2d::new(u_fil, v_fil);

    // `if (c1sphere) { ElSLib::SphereD1(u, v, FilAx3, cylrad, PP, deru, derv);
    // } else { ElSLib::TorusD1(u, v, FilAx3, Rad, Radius, PP, deru, derv); }
    // norFil = deru.Crossed(derv);` (`cxx:524-532`).
    let nor_fil2 = fillet_normal(
        fil_ax3.clone(),
        rad,
        radius,
        c1sphere,
        cylrad,
        u_fil,
        v_fil,
    )?;

    // `ElSLib::Parameters(Cyl, P, u, v); double tol = Precision::PConfusion();
    // bool careaboutsens = false; if (std::abs(lu - fu - 2 * M_PI) < tol) {
    // careaboutsens = true; } if (u >= fu - tol && u < fu) { u = fu; } if (u <=
    // lu + tol && u > lu) { u = lu; } if (u < fu || u > lu) { u =
    // ChFiKPart_InPeriod(u, fu, fu + 2 * M_PI, tol); }` (`cxx:535-553`).
    let (mut u, v) = slib::cylinder_parameters(&cyl.position(), &p_cyl);
    let tol = PCONFUSION;
    let careaboutsens = (lu - fu - 2.0 * std::f64::consts::PI).abs() < tol;
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    if u < fu || u > lu {
        u = ch_fi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, tol);
    }

    // `ElSLib::D1(u, v, Cyl, PP, deru, derv); gp_Dir norcyl =
    // deru.Crossed(derv); gp_Dir2d d2dCyl = gp::DX2d(); if (deru.Dot(Dy) < 0.)
    // { d2dCyl.Reverse(); if (careaboutsens && std::abs(fu - u) < tol) { u =
    // lu; } } else if (careaboutsens && std::abs(lu - u) < tol) { u = fu; }`
    // (`cxx:554-568`).
    let (_pp, deru, derv) = surface_eval::cylinder_d1(cyl, u, v);
    let nor_cyl = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    let mut d2d_cyl = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    if deru.dot(&GpVec::from_xyz(dy.xyz())) < 0.0 {
        d2d_cyl.reverse();
        if careaboutsens && (fu - u).abs() < tol {
            u = lu;
        }
    } else if careaboutsens && (lu - u).abs() < tol {
        u = fu;
    }

    // `gp_Pnt2d p2dCyl(u, v); gp_Lin2d lin2dCyl(p2dCyl, d2dCyl);
    // GLin2dCyl = new Geom2d_Line(lin2dCyl); circAx2.SetLocation(Or); gp_Circ
    // circCyl(circAx2, cylrad); GCircCyl = new Geom_Circle(circCyl);`
    // (`cxx:569-574`).
    let lin2d_cyl = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), d2d_cyl);
    circ_ax2.set_location(or);
    let circ_cyl = GpCirc::new(circ_ax2, cylrad);

    // `toreverse = (norFil.Dot(norcyl) <= 0.); if ((toreverse && plandab) ||
    // (!toreverse && !plandab)) { trans = TopAbs_REVERSED; } else { trans =
    // TopAbs_FORWARD; }` (`cxx:575-583`).
    let toreverse = nor_fil2.dot(&nor_cyl) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };
    let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_cyl));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
    let on_face: Arc<dyn Curve2d> = Arc::new(lin2d_cyl);
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

/// `P.SetCoord(base.X() + d * dir.X(), ...)`, the three occurrences of that
/// pattern in `cxx:319`, `:404`, `:506`.
fn shifted(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// `ElSLib::SphereD1` / `ElSLib::TorusD1` followed by `gp_Dir(deru.Crossed
/// (derv))` (`cxx:421-435`, `:524-532`). For the sphere the radius argument is
/// the caller's `cylrad` (`cxx:424`, `:526`), not the fillet radius.
#[allow(clippy::too_many_arguments)]
fn fillet_normal(
    ax3: GpAx3,
    rad: f64,
    radius: f64,
    c1sphere: bool,
    cylrad: f64,
    u: f64,
    v: f64,
) -> Result<GpDir, String> {
    let (deru, derv) = if c1sphere {
        let s = GpSphere::new(ax3, cylrad).map_err(|e| e.to_string())?;
        let (_, du, dv) = surface_eval::sphere_d1(&s, u, v);
        (du, dv)
    } else {
        let t = GpTorus::new(ax3, rad, radius).map_err(|e| e.to_string())?;
        let (_, du, dv) = surface_eval::torus_d1(&t, u, v);
        (du, dv)
    };
    GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())
}
