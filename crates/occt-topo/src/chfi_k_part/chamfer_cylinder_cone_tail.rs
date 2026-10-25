//! The tail shared by the two plane/cylinder chamfer makers whose guideline is
//! a circle: everything from `// compute the parameters of the conical
//! surface` to `return true`. The corner surface there is a cone.
//!
//! `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnCyl.cxx:145-362`
//! and `ChFiKPart_ComputeData_ChAsymPlnCyl.cxx:151-367` (OCCT 8.0.0) are the
//! same statements in that range apart from blank lines — the two makers differ
//! only in how `Dis1`/`Dis2` are derived, which stays in each caller's head.
//! Callers therefore pass the locals alive at that point, notably the already
//! shifted `Or` and the `Vec1`/`Vec2` whose angle is the cone's semi-angle.

use std::sync::Arc;

use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{
    GpAx22d, GpAx3, GpCirc, GpCirc2d, GpCone, GpCylinder, GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec,
    GpVec2d,
};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::cone::GeomCone;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{
    ch_fi_kpart_in_period, ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds,
};

/// `ChPlnCyl.cxx:145-362` / `ChAsymPlnCyl.cxx:151-367`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_load_chamfer_cylinder_cone(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pos_pl: &GpAx3,
    cyl: &GpCylinder,
    dpl: GpDir,
    norf: &GpDir,
    dx: &GpDir,
    dy: &GpDir,
    d_sp: &GpVec,
    or: &GpPnt,
    pt_pl: &GpPnt,
    pt2d_pln: GpPnt2d,
    vec1: &GpDir,
    vec2: &GpDir,
    dedans: bool,
    plandab: bool,
    dis1: f64,
    dis2: f64,
    fu: f64,
    lu: f64,
) -> Result<bool, String> {
    // `if (dedans) { Rad = Cyl.Radius() - Dis1; if (std::abs(Rad) <=
    // Precision::Confusion()) { pointu = true; } if (Rad < 0) { ... return
    // false; } } else { Rad = Cyl.Radius() + Dis1; gp_Dir Dplr =
    // Dpl.Reversed(); Dpl = Dplr; }` (`cxx:146-166`). The `Rad < 0` arm's
    // `std::cout` is `OCCT_DEBUG` only.
    let mut pointu = false;
    let mut dpl = dpl;
    let rad;
    if dedans {
        rad = cyl.radius() - dis1;
        if rad.abs() <= CONFUSION {
            pointu = true;
        }
        if rad < 0.0 {
            return Ok(false);
        }
    } else {
        rad = cyl.radius() + dis1;
        dpl.reverse();
    }

    // `ConRad = Cyl.Radius(); SemiAngl = Vec1.Angle(Vec2); gp_Ax3 ConAx3(Or,
    // Dpl, Dx); occ::handle<Geom_ConicalSurface> gcon = new
    // Geom_ConicalSurface(ConAx3, SemiAngl, ConRad);` (`cxx:167-171`).
    let con_rad = cyl.radius();
    let mut semi_angl = vec1.angle(vec2);
    let mut con_ax3 = GpAx3::from_point_dir_x(*or, dpl, dx).map_err(|e| e.to_string())?;
    let mut gcon = build_cone(&con_ax3, con_rad, semi_angl)?;

    // `if ((dedans && !plandab) || (!dedans && plandab)) { gcon->VReverse();
    // ConAx3 = gcon->Position(); SemiAngl = gcon->SemiAngle(); }`
    // (`cxx:175-181`). `Geom_ConicalSurface::VReverse` negates the semi-angle
    // and reverses the frame's Z axis (`Geom_ConicalSurface.cxx:104-108`).
    if (dedans && !plandab) || (!dedans && plandab) {
        semi_angl = -semi_angl;
        con_ax3.z_reverse();
        gcon = build_cone(&con_ax3, con_rad, semi_angl)?;
    }

    // `if (ConAx3.YDirection().Dot(DSp) <= 0.) { ConAx3.YReverse();
    // gcon->SetPosition(ConAx3); }` (`cxx:185-189`). `SetPosition` only stores
    // the frame (`Geom_ElementarySurface.hxx:74`), so the cone is rebuilt from
    // the same radius and semi-angle.
    if con_ax3
        .y_direction()
        .dot(&GpDir::from_vec(d_sp).map_err(|e| e.to_string())?)
        <= 0.0
    {
        con_ax3.y_reverse();
        gcon = build_cone(&con_ax3, con_rad, semi_angl)?;
    }

    // `Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gcon, DStr));` (`cxx:190`).
    data.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::new(GeomCone::new(gcon.clone())),
        ds,
    ));

    // `gp_Vec deru, derv; ElSLib::ConeD1(0., 0., ConAx3, ConRad, SemiAngl, Pt,
    // deru, derv); gp_Dir norCon(deru.Crossed(derv));` (`cxx:197-200`). The
    // point this call writes into `Pt` — the cone's `(0, 0)` point,
    // `Or + ConRad * Dx` — is the one the plane pcurve is computed at
    // (`cxx:231`), so it is kept.
    let (pt_con00, deru, derv) = surface_eval::cone_d1(&gcon, 0.0, 0.0);
    let nor_con = unit(&deru.crossed(&derv))?;

    // `bool toreverse = (norCon.Dot(norf) <= 0.); if (toreverse) {
    // Data->ChangeOrientation() = TopAbs_REVERSED; } else { ... FORWARD; }`
    // (`cxx:202-207`).
    let toreverse = nor_con.dot(norf) <= 0.0;
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // Planar face (`cxx:213-282`).
    // `occ::handle<Geom2d_Circle> GCir2dPln; occ::handle<Geom_Circle> GCirPln;
    // gp_Ax2 CirAx2 = ConAx3.Ax2(); CirAx2.SetLocation(PtPl);` (`cxx:217-222`).
    let mut circ_ax2 = con_ax3.ax2();
    circ_ax2.set_location(*pt_pl);
    let mut circ_pln: Option<GpCirc> = None;
    let mut circ2d_pln: Option<GpCirc2d> = None;
    if !pointu {
        // `gp_Circ CirPln(CirAx2, Rad); GCirPln = new Geom_Circle(CirPln);`
        // (`cxx:226-228`).
        circ_pln = Some(GpCirc::new(circ_ax2, rad));

        // `ElSLib::PlaneParameters(PosPl, Pt, u, v); gp_Pnt2d p2dPln(u, v);
        // gp_Dir2d d2d(DSp.Dot(PosPl.XDirection()),
        // DSp.Dot(PosPl.YDirection())); gp_Ax22d ax2dPln(pt2dPln,
        // gp_Dir2d(gp_Vec2d(pt2dPln, p2dPln)), d2d); gp_Circ2d
        // cir2dPln(ax2dPln, Rad);` (`cxx:231-239`).
        let (u, v) = slib::plane_parameters(pos_pl, &pt_con00);
        let p2d_pln = GpPnt2d::new(u, v);
        let d2d = GpDir2d::new(
            d_sp.dot(&GpVec::from_xyz(pos_pl.x_direction().xyz())),
            d_sp.dot(&GpVec::from_xyz(pos_pl.y_direction().xyz())),
        )
        .map_err(|e| e.to_string())?;
        let v2d = GpVec2d::new(p2d_pln.x() - pt2d_pln.x(), p2d_pln.y() - pt2d_pln.y());
        let vx = GpDir2d::new(v2d.x(), v2d.y()).map_err(|e| e.to_string())?;
        // `gp_Ax22d(P, Vx, Vy)` stores both directions as given
        // (`gp_Ax22d.hxx:60-66`); `GpAx22d::new` would reject a non
        // perpendicular pair instead, so the fields are set directly.
        let ax2d_pln = GpAx22d {
            point: pt2d_pln,
            vxdir: vx,
            vydir: d2d,
        };
        circ2d_pln = Some(GpCirc2d::new(ax2d_pln, rad));
    }

    // `gp_Pnt2d p2dch; if (plandab) { v = -sqrt(Dis1 * Dis1 + Dis2 * Dis2); }
    // else { v = sqrt(...); } p2dch.SetCoord(0., v); ElSLib::ConeD1(0., v,
    // ConAx3, ConRad, SemiAngl, Pt, deru, derv); gp_Lin2d lin2dch(p2dch,
    // gp::DX2d());` (`cxx:245-252`).
    let mut p2dch = GpPnt2d::new(0.0, 0.0);
    let chord = (dis1 * dis1 + dis2 * dis2).sqrt();
    p2dch.set_coord(0.0, if plandab { -chord } else { chord });
    // `ElSLib::ConeD1(0., v, ConAx3, ConRad, SemiAngl, Pt, deru, derv)`
    // (`cxx:250`) only refreshes `deru`/`derv`, which are overwritten by
    // `ElSLib::ConeD1(0., 0., ...)` at `cxx:295` before being read, so the
    // call has no observable effect and is not reproduced.
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let lin2d_ch = Geom2dLine::from_pnt_dir(p2dch, dx2d);

    // `TopAbs_Orientation trans; gp_Dir norpl = PosPl.XDirection().Crossed
    // (PosPl.YDirection()); toreverse = (norCon.Dot(norpl) <= 0.); if
    // ((toreverse && plandab) || (!toreverse && !plandab)) { trans =
    // TopAbs_FORWARD; } else { trans = TopAbs_REVERSED; }` (`cxx:256-264`).
    let norpl = pos_pl
        .x_direction()
        .crossed(pos_pl.y_direction())
        .map_err(|e| e.to_string())?;
    let toreverse = nor_con.dot(&norpl) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };

    // `if (plandab) { Data->ChangeInterferenceOnS1().SetInterference(
    // IndexCurveInDS(GCirPln, DStr), trans, GCir2dPln, GLin2dCh1); } else {
    // ...OnS2... }` (`cxx:266-281`).
    let circ_3d: Option<Arc<dyn GeomCurve>> =
        circ_pln.map(|c| Arc::new(GeomCircle::new(c)) as Arc<dyn GeomCurve>);
    let ic = ch_fi_kpart_index_curve_in_ds(circ_3d, ds);
    let on_face: Option<Arc<dyn Curve2d>> =
        circ2d_pln.map(|c| Arc::new(Geom2dCircle::new(c)) as Arc<dyn Curve2d>);
    let on_fil: Arc<dyn Curve2d> = Arc::new(lin2d_ch);
    if plandab {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    } else {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    }

    // Cylindrical face (`cxx:284-362`).
    // `CirAx2.SetLocation(Or); gp_Circ CirCyl(CirAx2, ConRad);` (`cxx:288-290`).
    circ_ax2.set_location(*or);
    let circ_cyl = GpCirc::new(circ_ax2, con_rad);

    // `p2dch.SetCoord(0., 0.); ElSLib::ConeD1(0., 0., ConAx3, ConRad,
    // SemiAngl, Pt, deru, derv); lin2dch.SetLocation(p2dch);` (`cxx:294-296`).
    p2dch.set_coord(0.0, 0.0);
    let (_, deru2, derv2) = surface_eval::cone_d1(&gcon, 0.0, 0.0);
    let lin2d_ch2 = Geom2dLine::from_pnt_dir(p2dch, dx2d);

    // `norCon.SetXYZ(deru.Crossed(derv).XYZ());` (`cxx:297`).
    let nor_con = unit(&deru2.crossed(&derv2))?;

    // `Pt.SetCoord(Or.X() + ConRad * Dx.X(), ...); ElSLib::Parameters(Cyl, Pt,
    // u, v); double tol = Precision::PConfusion(); bool careaboutsens = false;
    // if (std::abs(lu - fu - 2 * M_PI) < tol) { careaboutsens = true; }`
    // (`cxx:299-308`).
    let pt_cyl_on_ax = shifted(or, con_rad, dx);
    let (mut u, v) = slib::cylinder_parameters(&cyl.position(), &pt_cyl_on_ax);
    let tol = PCONFUSION;
    let careaboutsens = (lu - fu - 2.0 * std::f64::consts::PI).abs() < tol;
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    // `if (u < fu || u > lu) { u = ChFiKPart_InPeriod(u, fu, fu + 2 * M_PI,
    // tol); }` (`cxx:315-318`).
    if u < fu || u > lu {
        u = ch_fi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, tol);
    }

    // `ElSLib::D1(u, v, Cyl, Pt, deru, derv); gp_Dir norcyl =
    // deru.Crossed(derv); gp_Dir2d d2dCyl = gp::DX2d(); if (deru.Dot(Dy) < 0.)
    // { d2dCyl.Reverse(); if (careaboutsens && std::abs(fu - u) < tol) { u =
    // lu; } } else if (careaboutsens && std::abs(lu - u) < tol) { u = fu; }`
    // (`cxx:319-333`).
    let (_, du, dv) = surface_eval::cylinder_d1(cyl, u, v);
    let nor_cyl = unit(&du.crossed(&dv))?;
    let mut d2d_cyl = dx2d;
    if du.dot(&GpVec::from_xyz(dy.xyz())) < 0.0 {
        d2d_cyl.reverse();
        if careaboutsens && (fu - u).abs() < tol {
            u = lu;
        }
    } else if careaboutsens && (lu - u).abs() < tol {
        u = fu;
    }
    let lin2d_cyl = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), d2d_cyl);

    // `toreverse = (norCon.Dot(norcyl) <= 0.); if ((toreverse && plandab) ||
    // (!toreverse && !plandab)) { trans = TopAbs_REVERSED; } else { ... }`
    // (`cxx:341-348`).
    let toreverse = nor_con.dot(&nor_cyl) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };

    // `if (plandab) { ...OnS2... } else { ...OnS1... }` (`cxx:350-362`).
    let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_cyl));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
    let on_face: Arc<dyn Curve2d> = Arc::new(lin2d_cyl);
    let on_fil: Arc<dyn Curve2d> = Arc::new(lin2d_ch2);
    if plandab {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    } else {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    }

    Ok(true)
}

/// `new Geom_ConicalSurface(ConAx3, SemiAngl, ConRad)` (`cxx:171`), whose
/// constructor raises `Standard_ConstructionError` for a negative radius or an
/// out-of-range semi-angle (`Geom_ConicalSurface.cxx:52-62`).
fn build_cone(ax3: &GpAx3, radius: f64, semi_angle: f64) -> Result<GpCone, String> {
    GpCone::new(ax3.clone(), radius, semi_angle).map_err(|e| e.to_string())
}

/// `P.SetCoord(base.X() + d * dir.X(), ...)`, the pattern of `cxx:299`.
fn shifted(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// `gp_Dir(V)` (`cxx:200`, `:320`): normalizes and rejects a zero vector, which
/// is what `gp_Dir`'s constructor raises `Standard_ConstructionError` for.
fn unit(v: &GpVec) -> Result<GpDir, String> {
    GpDir::from_vec(v).map_err(|e| e.to_string())
}
