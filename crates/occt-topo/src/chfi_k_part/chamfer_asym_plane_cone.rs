//! `ChFiKPart_ComputeData_ChAsymPlnCon` — the distance/angle chamfer between a
//! plane and a cone, guideline a circle.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChAsymPlnCon.cxx:57-705`
//! (OCCT 8.0.0). It serves the `ChFiDS_DistAngle` cone arms of
//! `ChFiKPart_ComputeData.cxx:501-513` / `:618-636`, and also the
//! `ChFiDS_Sym`/`TwoDist` cone arms, which
//! [`super::chamfer_plane_cone`] reduces to one distance plus one angle first.
//!
//! Two ways of reading `Dis`/`Angle` are possible and OCCT picks between them
//! with `(plandab && DisOnP) || (!plandab && !DisOnP)` (`cxx:122`):
//!
//! * the distance is given on the cone, so the chamfer's own half-section gives
//!   `move` along the plane normal and `dis` along the plane's X direction
//!   (`cxx:124-170`);
//! * the distance is given on the plane, so the chamfer is a triangle attached
//!   to the cone's slope (`cxx:172-251`); a zero half-angle there means the
//!   chamfer is a cylinder (`iscylinder`), and a positive one that the chamfer
//!   is parallel to the cone (`isConPar`) — both flag a sign flip later on.
//!
//! Everything from `Or.SetCoord(...)` on is one of two shapes: a cylindrical
//! surface (`cxx:260-464`) or a conical one (`cxx:465-704`).

use std::sync::Arc;

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{
    GpAx22d, GpAx3, GpCirc, GpCirc2d, GpCone, GpCylinder, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d,
    GpVec, GpVec2d,
};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::cone::GeomCone;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::cylinder::GeomCylinder;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChFiKPart_MakeChAsym(DStr, Data, Pln, Con, fu, lu, Or1, Or2, Dis, Angle,
/// Spine, First, Ofpl, plandab, DisOnP)` (`cxx:57-705`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_ch_asym_plane_cone(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cone: &GpCone,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> Result<bool, String> {
    // `gp_Ax3 PosPl = Pln.Position(); gp_Dir Dpl = PosPl.XDirection().Crossed
    // (PosPl.YDirection()); gp_Dir norf = Dpl; if (Ofpl == TopAbs_REVERSED) {
    // norf.Reverse(); } if (Or1 == TopAbs_REVERSED) { Dpl.Reverse(); }`
    // (`cxx:74-84`).
    let pos_pl = pln.position();
    let mut dpl = pos_pl
        .x_direction()
        .crossed(pos_pl.y_direction())
        .map_err(|e| e.to_string())?;
    let mut norf = dpl;
    if ofpl == Orientation::Reversed {
        norf.reverse();
    }
    if or1 == Orientation::Reversed {
        dpl.reverse();
    }

    // `gp_Pnt Or = Con.Location(); double u, v; ElSLib::PlaneParameters(PosPl,
    // Or, u, v); gp_Pnt2d pt2dPln(u, v); ElSLib::PlaneD0(u, v, PosPl, Or);
    // gp_Pnt PtPl = Or;` (`cxx:87-92`).
    let or_base = cone.location();
    let (u0, v0) = slib::plane_parameters(&pos_pl, &or_base);
    let pt2d_pln = GpPnt2d::new(u0, v0);
    let pt_pl = slib::plane_value(pln, u0, v0);

    // `gp_Pnt PtSp; gp_Vec DSp; ElCLib::D1(First, Spine, PtSp, DSp); gp_Dir
    // Dx(gp_Vec(PtPl, PtSp));` (`cxx:94-97`). Note `Dx` runs from the *plane*
    // origin here, unlike the cylinder makers' `gp_Vec(Or, PtSp)`.
    let (pt_sp, d_sp) = clib::circle_d1(spine, first);
    let dx = GpDir::from_vec(&GpVec::from_pnts(&pt_pl, &pt_sp)).map_err(|e| e.to_string())?;

    // `gp_Vec deru, derv; gp_Pnt PtCon; ElSLib::Parameters(Con, PtSp, u, v);
    // ElSLib::D1(u, v, Con, PtCon, deru, derv); gp_Dir Dcon(deru.Crossed
    // (derv)); if (Or2 == TopAbs_REVERSED) { Dcon.Reverse(); }` (`cxx:100-108`).
    let (u, v) =
        slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &pt_sp);
    let (_pt_con, deru, derv) = surface_eval::cone_d1(cone, u, v);
    let mut dcon = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dcon.reverse();
    }

    // `bool dedans = (Dx.Dot(Dcon) <= 0.); bool ouvert = (Dpl.Dot(Dcon) >= 0.);`
    // (`cxx:110-111`).
    let dedans = dx.dot(&dcon) <= 0.0;
    let ouvert = dpl.dot(&dcon) >= 0.0;

    // `double angCon = Con.SemiAngle(); double move; double ChamfRad,
    // SemiAngl; bool pointu = false; double dis; bool iscylinder = false; bool
    // isConPar = false;` (`cxx:114-121`). OCCT's local `dis` is kept apart from
    // the input `Dis` here as `dis_x`: it is the offset along the plane's X
    // direction used from `cxx:253` on, while the input `Dis` is only read in
    // the two `cxx:124-250` branches.
    let angcon = cone.semi_angle();
    let mut pointu = false;
    let mut iscylinder = false;
    let mut is_con_par = false;
    let move_;
    let dis_x;
    let chamf_rad;
    let mut semi_angl;

    if (plandab && dis_on_p) || (!plandab && !dis_on_p) {
        // `double tgang = std::tan(Angle), Dis11; double tgCon =
        // std::abs(std::tan(angCon)); if (ouvert) { move = Dis * tgang / (1. -
        // tgCon * tgang); Dis11 = move * tgCon; dis = Dis + Dis11; } else {
        // move = Dis * tgang / (1. + tgCon * tgang); Dis11 = move * tgCon; dis
        // = Dis - Dis11; }` (`cxx:124-139`).
        let tgang = angle.tan();
        let tgcon = angcon.tan().abs();
        if ouvert {
            move_ = dis * tgang / (1.0 - tgcon * tgang);
            dis_x = dis + move_ * tgcon;
        } else {
            move_ = dis * tgang / (1.0 + tgcon * tgang);
            dis_x = dis - move_ * tgcon;
        }

        // `if (dedans) { ChamfRad = Spine.Radius() - Dis; if (std::abs
        // (ChamfRad) < Precision::Confusion()) { pointu = true; } if (ChamfRad
        // < 0) { ... return false; } SemiAngl = M_PI / 2. - Angle; } else {
        // ChamfRad = Spine.Radius() + Dis; SemiAngl = Angle - M_PI / 2.; }`
        // (`cxx:140-161`). The `ChamfRad < 0` arm's `std::cout` is
        // `OCCT_DEBUG` only.
        if dedans {
            chamf_rad = spine.radius() - dis;
            if chamf_rad.abs() < CONFUSION {
                pointu = true;
            }
            if chamf_rad < 0.0 {
                return Ok(false);
            }
            semi_angl = std::f64::consts::FRAC_PI_2 - angle;
        } else {
            chamf_rad = spine.radius() + dis;
            semi_angl = angle - std::f64::consts::FRAC_PI_2;
        }

        // `if (ouvert) { if (std::abs(angCon) - std::abs(SemiAngl) >
        // -Precision::Confusion()) { ... return false; } }` (`cxx:163-171`).
        if ouvert && angcon.abs() - semi_angl.abs() > -CONFUSION {
            return Ok(false);
        }
    } else {
        // `double Dis1; move = Dis * std::cos(angCon); if (ouvert) { SemiAngl =
        // std::abs(angCon) + Angle; if ((M_PI / 2. - SemiAngl) <
        // Precision::Confusion()) { ... return false; } Dis1 = move *
        // std::tan(SemiAngl) - Dis * std::abs(std::sin(angCon)); if (!dedans) {
        // SemiAngl = -SemiAngl; } } else { SemiAngl = std::abs(angCon) -
        // Angle; if (std::abs(SemiAngl) < Precision::Confusion()) {
        // iscylinder = true; Dis1 = Dis * std::abs(std::sin(angCon)); } else {
        // Dis1 = Dis * std::abs(std::sin(angCon)) - move * std::tan(SemiAngl);
        // } if (SemiAngl > Precision::Confusion()) { isConPar = true; } if
        // (dedans) { SemiAngl = -SemiAngl; } }` (`cxx:172-216`).
        move_ = dis * angcon.cos();
        let dis1;
        if ouvert {
            semi_angl = angcon.abs() + angle;
            if std::f64::consts::FRAC_PI_2 - semi_angl < CONFUSION {
                return Ok(false);
            }
            dis1 = move_ * semi_angl.tan() - dis * angcon.sin().abs();
            if !dedans {
                semi_angl = -semi_angl;
            }
        } else {
            semi_angl = angcon.abs() - angle;
            if semi_angl.abs() < CONFUSION {
                iscylinder = true;
                dis1 = dis * angcon.sin().abs();
            } else {
                dis1 = dis * angcon.sin().abs() - move_ * semi_angl.tan();
            }
            if semi_angl > CONFUSION {
                is_con_par = true;
            }
            if dedans {
                semi_angl = -semi_angl;
            }
        }

        // `if (dedans) { ChamfRad = Spine.Radius() - Dis1; if (std::abs
        // (ChamfRad) < Precision::Confusion()) { pointu = true; } if (ChamfRad
        // < 0) { ... return false; } } else { ChamfRad = Spine.Radius() +
        // Dis1; }` (`cxx:222-240`).
        if dedans {
            chamf_rad = spine.radius() - dis1;
            if chamf_rad.abs() < CONFUSION {
                pointu = true;
            }
            if chamf_rad < 0.0 {
                return Ok(false);
            }
        } else {
            chamf_rad = spine.radius() + dis1;
        }

        // `if (ouvert) { dis = Dis1 + Dis * std::abs(std::sin(angCon)); } else
        // { dis = Dis1 - Dis * std::abs(std::sin(angCon)); }` (`cxx:242-250`).
        if ouvert {
            dis_x = dis1 + dis * angcon.sin().abs();
        } else {
            dis_x = dis1 - dis * angcon.sin().abs();
        }
    }

    // `Or.SetCoord(Or.X() + move * Dpl.X(), ...); gp_Pnt Pt(Or.X() + dis *
    // PosPl.XDirection().X(), ...); gp_Ax3 ChamfAx3(PtPl, Dpl, Dx);`
    // (`cxx:253-259`). `ElSLib::PlaneD0(u, v, PosPl, Or)` wrote the projected
    // point back into `Or` (`cxx:89-91`), so the offset starts from `PtPl`.
    // OCCT's `Pt` is overwritten before any use in both the cylindrical
    // (`cxx:314`) and the conical (`cxx:537`) arm, so the point is built there
    // instead of here.
    let or = offset(&pt_pl, move_, &dpl);
    let mut chamf_ax3 =
        GpAx3::from_point_dir_x(pt_pl, dpl, &dx).map_err(|e| e.to_string())?;

    if iscylinder {
        // `occ::handle<Geom_CylindricalSurface> gcyl = new
        // Geom_CylindricalSurface(ChamfAx3, ChamfRad);` (`cxx:262`).
        let mut gcyl = GpCylinder::new(chamf_ax3.clone(), chamf_rad).map_err(|e| e.to_string())?;

        // `if (!plandab) { gcyl->VReverse(); ChamfAx3 = gcyl->Position(); }`
        // (`cxx:266-270`). `Geom_ElementarySurface::VReverse` only reverses
        // the frame's Z axis (`Geom_ElementarySurface.cxx:51-54`).
        if !plandab {
            chamf_ax3.z_reverse();
            gcyl = GpCylinder::new(chamf_ax3.clone(), chamf_rad).map_err(|e| e.to_string())?;
        }

        // `if (ChamfAx3.YDirection().Dot(DSp) <= 0.) { ChamfAx3.YReverse();
        // gcyl->SetPosition(ChamfAx3); }` (`cxx:273-277`).
        if chamf_ax3
            .y_direction()
            .dot(&GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?)
            <= 0.0
        {
            chamf_ax3.y_reverse();
            gcyl = GpCylinder::new(chamf_ax3.clone(), chamf_rad).map_err(|e| e.to_string())?;
        }

        // `Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gcyl, DStr));`
        // (`cxx:280`).
        data.set_surf(ch_fi_kpart_index_surface_in_ds(
            Arc::new(GeomCylinder::new(gcyl.clone())),
            ds,
        ));

        // `bool torevcha = !ChamfAx3.Direct(); gp_Dir cylaxe =
        // (ChamfAx3.Axis()).Direction(); torevcha = ((torevcha && !plandab) ||
        // (!torevcha && plandab)); if (torevcha) { cylaxe.Reverse(); } bool
        // toreverse = (norf.Dot(cylaxe) < 0.);` (`cxx:282-290`).
        let torevcha = !chamf_ax3.is_direct();
        let mut cyl_axe = chamf_ax3.axis().direction().clone();
        let torevcha = (torevcha && !plandab) || (!torevcha && plandab);
        if torevcha {
            cyl_axe.reverse();
        }
        let toreverse = norf.dot(&cyl_axe) < 0.0;

        // `if ((toreverse && dedans) || (!toreverse && !dedans)) {
        // Data->ChangeOrientation() = TopAbs_REVERSED; } else { ...
        // TopAbs_FORWARD; }` (`cxx:292-299`).
        data.set_orientation(if (toreverse && dedans) || (!toreverse && !dedans) {
            Orientation::Reversed
        } else {
            Orientation::Forward
        });

        // Planar face (`cxx:302-364`).
        // `gp_Ax2 CirAx2 = ChamfAx3.Ax2(); CirAx2.SetLocation(PtPl);
        // Pt.SetCoord(PtPl.X() + ChamfRad * Dx.X(), ...); gp_Circ
        // CirPln(CirAx2, ChamfRad);` (`cxx:311-318`).
        let mut circ_ax2 = chamf_ax3.ax2();
        circ_ax2.set_location(pt_pl);
        let pt = offset(&pt_pl, chamf_rad, &dx);
        let circ_pln = GpCirc::new(circ_ax2, chamf_rad);

        // `ElSLib::PlaneParameters(PosPl, Pt, u, v); gp_Pnt2d p2dPln(u, v);
        // gp_Dir2d d2d(DSp.Dot(PosPl.XDirection()), DSp.Dot
        // (PosPl.YDirection())); gp_Ax22d ax2dPln(pt2dPln, gp_Dir2d
        // (gp_Vec2d(pt2dPln, p2dPln)), d2d); gp_Circ2d cir2dPln(ax2dPln,
        // ChamfRad);` (`cxx:321-327`).
        let (u, v) = slib::plane_parameters(&pos_pl, &pt);
        let p2d_pln = GpPnt2d::new(u, v);
        let d2d = GpDir2d::new(
            d_sp.dot(&vec_of(pos_pl.x_direction())),
            d_sp.dot(&vec_of(pos_pl.y_direction())),
        )
        .map_err(|e| e.to_string())?;
        let v2d = GpVec2d::new(p2d_pln.x() - pt2d_pln.x(), p2d_pln.y() - pt2d_pln.y());
        let vx = GpDir2d::new(v2d.x(), v2d.y()).map_err(|e| e.to_string())?;
        let ax2d_pln = GpAx22d {
            point: pt2d_pln,
            vxdir: vx,
            vydir: d2d,
        };
        let circ2d_pln = GpCirc2d::new(ax2d_pln, chamf_rad);

        // `gp_Pnt2d p2dch; p2dch.SetCoord(0., 0.); gp_Lin2d lin2dch(p2dch,
        // gp::DX2d());` (`cxx:329-334`).
        let mut p2dch = GpPnt2d::new(0.0, 0.0);
        let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
        let mut lin2dch = Geom2dLine::from_pnt_dir(p2dch, dx2d);

        // `TopAbs_Orientation trans; gp_Dir norpl = PosPl.XDirection().Crossed
        // (PosPl.YDirection()); toreverse = (norpl.Dot(cylaxe) < 0.);
        // toreverse = (toreverse && plandab) || (!toreverse && !plandab); if
        // ((toreverse && dedans) || (!toreverse && !dedans)) { trans =
        // TopAbs_FORWARD; } else { trans = TopAbs_REVERSED; }` (`cxx:337-349`).
        let norpl = pos_pl
            .x_direction()
            .crossed(pos_pl.y_direction())
            .map_err(|e| e.to_string())?;
        let toreverse = norpl.dot(&cyl_axe) < 0.0;
        let toreverse = (toreverse && plandab) || (!toreverse && !plandab);
        let trans = if (toreverse && dedans) || (!toreverse && !dedans) {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };

        // `if (plandab) { ...OnS1... } else { ...OnS2... }` (`cxx:351-364`).
        let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_pln));
        let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
        let on_face: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(circ2d_pln));
        let on_fil: Arc<dyn Curve2d> = Arc::new(lin2dch.clone());
        if plandab {
            data.interference_on_s1_mut()
                .set_interference(ic, trans, Some(on_face), Some(on_fil));
        } else {
            data.interference_on_s2_mut()
                .set_interference(ic, trans, Some(on_face), Some(on_fil));
        }

        // Conical face (`cxx:366-464`).
        // `double Rad; if (dedans) { Rad = ChamfRad + dis; } else { Rad =
        // ChamfRad - dis; } CirAx2.SetLocation(Or); gp_Circ CirCon(CirAx2,
        // Rad);` (`cxx:368-381`).
        let rad = if dedans { chamf_rad + dis_x } else { chamf_rad - dis_x };
        circ_ax2.set_location(or);
        let circ_con = GpCirc::new(circ_ax2, rad);

        // `if (plandab) { v = sqrt(dis * dis + move * move); } else { v =
        // -sqrt(...); } p2dch.SetCoord(0., v); ElSLib::CylinderD1(0., v,
        // ChamfAx3, ChamfRad, Pt, deru, derv); lin2dch.SetLocation(p2dch);`
        // (`cxx:384-395`).
        let chord = (dis_x * dis_x + move_ * move_).sqrt();
        p2dch.set_coord(0.0, if plandab { chord } else { -chord });
        // `ElSLib::CylinderD1(0., v, ChamfAx3, ChamfRad, Pt, deru, derv)`
        // (`cxx:392`) only refreshes `deru`/`derv`, which `ElSLib::D1(u, v,
        // Con, ...)` overwrites at `cxx:415` before they are read, so the call
        // has no observable effect and is not reproduced.
        lin2dch = Geom2dLine::from_pnt_dir(p2dch, dx2d);

        // `Pt.SetCoord(Or.X() + Rad * Dx.X(), ...); ElSLib::Parameters(Con, Pt,
        // u, v); double tol = Precision::PConfusion(); if (u >= 2 * M_PI - tol
        // && u <= 2 * M_PI) { u = 0.; } if (u >= fu - tol && u < fu) { u = fu;
        // } if (u <= lu + tol && u > lu) { u = lu; } if (u < fu || u > lu) { u
        // = ElCLib::InPeriod(u, fu, fu + 2 * M_PI); }` (`cxx:397-414`).
        let pt_con = offset(&or, rad, &dx);
        let (mut u, v) = slib::cone_parameters(
            &cone.position(),
            cone.radius(),
            cone.semi_angle(),
            &pt_con,
        );
        let tol = PCONFUSION;
        if u >= 2.0 * std::f64::consts::PI - tol && u <= 2.0 * std::f64::consts::PI {
            u = 0.0;
        }
        if u >= fu - tol && u < fu {
            u = fu;
        }
        if u <= lu + tol && u > lu {
            u = lu;
        }
        if u < fu || u > lu {
            u = clib::in_period(u, fu, fu + 2.0 * std::f64::consts::PI);
        }

        // `ElSLib::D1(u, v, Con, Pt, deru, derv); gp_Pnt2d p2dCon(u, v);
        // gp_Dir2d d2dCon; if (deru.Dot(DSp) <= 0.) { d2dCon = -gp::DX2d(); }
        // else { d2dCon = gp::DX2d(); }` (`cxx:415-427`).
        let (_pt, deru_con, derv_con) = surface_eval::cone_d1(cone, u, v);
        let mut d2d_con = dx2d;
        if deru_con.dot(&d_sp) <= 0.0
        {
            d2d_con.reverse();
        }
        let lin2d_con = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), d2d_con);

        // `gp_Dir norcon = deru.Crossed(derv); gp_Dir DirCon =
        // (Con.Axis()).Direction(); if (angCon > Precision::Confusion()) {
        // DirCon.Reverse(); } bool torevcon = (norcon.Dot(DirCon) < 0.); if
        // ((torevcon && dedans) || (!torevcon && !dedans)) { trans =
        // TopAbs_REVERSED; } else { trans = TopAbs_FORWARD; }` (`cxx:430-448`).
        let norcon = GpDir::from_vec(&deru_con.crossed(&derv_con)).map_err(|e| e.to_string())?;
        let mut dir_con = cone.axis().direction().clone();
        if angcon > CONFUSION {
            dir_con.reverse();
        }
        let torevcon = norcon.dot(&dir_con) < 0.0;
        let trans = if (torevcon && dedans) || (!torevcon && !dedans) {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };

        // `if (plandab) { ...OnS2... } else { ...OnS1... }` (`cxx:450-464`).
        let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_con));
        let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
        let on_face: Arc<dyn Curve2d> = Arc::new(lin2d_con);
        let on_fil: Arc<dyn Curve2d> = Arc::new(lin2dch);
        if plandab {
            data.interference_on_s2_mut()
                .set_interference(ic, trans, Some(on_face), Some(on_fil));
        } else {
            data.interference_on_s1_mut()
                .set_interference(ic, trans, Some(on_face), Some(on_fil));
        }
        return Ok(true);
    }

    // `occ::handle<Geom_ConicalSurface> gcon = new Geom_ConicalSurface
    // (ChamfAx3, SemiAngl, ChamfRad);` (`cxx:467`), then `if (!plandab) {
    // gcon->VReverse(); ChamfAx3 = gcon->Position(); SemiAngl =
    // gcon->SemiAngle(); }` (`cxx:471-475`).
    let mut gcon = build_cone(&chamf_ax3, chamf_rad, semi_angl)?;
    if !plandab {
        semi_angl = -semi_angl;
        chamf_ax3.z_reverse();
        gcon = build_cone(&chamf_ax3, chamf_rad, semi_angl)?;
    }

    // `if (ChamfAx3.YDirection().Dot(DSp) <= 0.) { ChamfAx3.YReverse();
    // gcon->SetPosition(ChamfAx3); }` (`cxx:479-483`).
    if chamf_ax3
        .y_direction()
        .dot(&GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?)
        <= 0.0
    {
        chamf_ax3.y_reverse();
        gcon = build_cone(&chamf_ax3, chamf_rad, semi_angl)?;
    }

    // `Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gcon, DStr));` (`cxx:485`).
    data.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::new(GeomCone::new(gcon.clone())),
        ds,
    ));

    // `gp_Pnt P; u = 0.; if (plandab) { v = sqrt(dis * dis + move * move); }
    // else { v = -sqrt(...); } ElSLib::ConeD1(u, v, ChamfAx3, ChamfRad,
    // SemiAngl, P, deru, derv); gp_Dir norchamf(deru.Crossed(derv));`
    // (`cxx:490-503`).
    let chord = (dis_x * dis_x + move_ * move_).sqrt();
    let (_, deru, derv) = surface_eval::cone_d1(&gcon, 0.0, if plandab { chord } else { -chord });
    let mut norchamf = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;

    // `bool toreverse = (norf.Dot(norchamf) < 0.); if (isConPar) { toreverse =
    // !toreverse; } if (toreverse) { Data->ChangeOrientation() =
    // TopAbs_REVERSED; } else { ...FORWARD; }` (`cxx:505-519`).
    let toreverse = norf.dot(&norchamf) < 0.0;
    let toreverse = if is_con_par { !toreverse } else { toreverse };
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // Planar face (`cxx:522-592`).
    // `gp_Ax2 CirAx2 = ChamfAx3.Ax2(); CirAx2.SetLocation(PtPl);` (`cxx:531-534`).
    let mut circ_ax2 = chamf_ax3.ax2();
    circ_ax2.set_location(pt_pl);
    let mut circ_pln: Option<GpCirc> = None;
    let mut circ2d_pln: Option<GpCirc2d> = None;
    if !pointu {
        // `Pt.SetCoord(PtPl.X() + ChamfRad * Dx.X(), ...); gp_Circ
        // CirPln(CirAx2, ChamfRad);` (`cxx:537-541`).
        let pt = offset(&pt_pl, chamf_rad, &dx);
        circ_pln = Some(GpCirc::new(circ_ax2, chamf_rad));

        // `ElSLib::PlaneParameters(PosPl, Pt, u, v); gp_Pnt2d p2dPln(u, v);
        // gp_Dir2d d2d(...); gp_Ax22d ax2dPln(pt2dPln, gp_Dir2d
        // (gp_Vec2d(pt2dPln, p2dPln)), d2d); gp_Circ2d cir2dPln(ax2dPln,
        // ChamfRad);` (`cxx:544-549`).
        let (u, v) = slib::plane_parameters(&pos_pl, &pt);
        let p2d_pln = GpPnt2d::new(u, v);
        let d2d = GpDir2d::new(
            d_sp.dot(&vec_of(pos_pl.x_direction())),
            d_sp.dot(&vec_of(pos_pl.y_direction())),
        )
        .map_err(|e| e.to_string())?;
        let v2d = GpVec2d::new(p2d_pln.x() - pt2d_pln.x(), p2d_pln.y() - pt2d_pln.y());
        let vx = GpDir2d::new(v2d.x(), v2d.y()).map_err(|e| e.to_string())?;
        let ax2d_pln = GpAx22d {
            point: pt2d_pln,
            vxdir: vx,
            vydir: d2d,
        };
        circ2d_pln = Some(GpCirc2d::new(ax2d_pln, chamf_rad));
    }

    // `gp_Pnt2d p2dch; p2dch.SetCoord(0., 0.); ElSLib::ConeD1(0., 0.,
    // ChamfAx3, ChamfRad, SemiAngl, Pt, deru, derv); gp_Lin2d lin2dch(p2dch,
    // gp::DX2d());` (`cxx:552-556`).
    let mut p2dch = GpPnt2d::new(0.0, 0.0);
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let (_, du, dv) = surface_eval::cone_d1(&gcon, 0.0, 0.0);
    let deru = du;
    let derv = dv;
    let mut lin2dch = Geom2dLine::from_pnt_dir(p2dch, dx2d);

    // `TopAbs_Orientation trans; gp_Dir norpl = PosPl.XDirection().Crossed
    // (PosPl.YDirection()); if (!pointu) { norchamf.SetXYZ(deru.Crossed(derv)
    // .XYZ()); } toreverse = (norchamf.Dot(norpl) <= 0.); if (isConPar) {
    // toreverse = !toreverse; } if ((toreverse && plandab) || (!toreverse &&
    // !plandab)) { trans = TopAbs_FORWARD; } else { trans = TopAbs_REVERSED;
    // }` (`cxx:560-581`).
    let norpl = pos_pl
        .x_direction()
        .crossed(pos_pl.y_direction())
        .map_err(|e| e.to_string())?;
    if !pointu {
        norchamf = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    }
    let toreverse = norchamf.dot(&norpl) <= 0.0;
    let toreverse = if is_con_par { !toreverse } else { toreverse };
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };

    // `if (plandab) { ...OnS1... } else { ...OnS2... }` (`cxx:583-594`).
    let circ_3d: Option<Arc<dyn GeomCurve>> =
        circ_pln.map(|c| Arc::new(GeomCircle::new(c)) as Arc<dyn GeomCurve>);
    let ic = ch_fi_kpart_index_curve_in_ds(circ_3d, ds);
    let on_face: Option<Arc<dyn Curve2d>> =
        circ2d_pln.map(|c| Arc::new(Geom2dCircle::new(c)) as Arc<dyn Curve2d>);
    let on_fil: Arc<dyn Curve2d> = Arc::new(lin2dch.clone());
    if plandab {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    } else {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, on_face, Some(on_fil));
    }

    // Conical face (`cxx:596-703`).
    // `double Rad; if (dedans) { Rad = ChamfRad + dis; } else { Rad = ChamfRad -
    // dis; } CirAx2.SetLocation(Or); gp_Circ CirCon(CirAx2, Rad);`
    // (`cxx:600-613`).
    let rad = if dedans { chamf_rad + dis_x } else { chamf_rad - dis_x };
    circ_ax2.set_location(or);
    let circ_con = GpCirc::new(circ_ax2, rad);

    // `if (plandab) { v = sqrt(...); } else { v = -sqrt(...); }
    // p2dch.SetCoord(0., v); ElSLib::ConeD1(0., v, ChamfAx3, ChamfRad,
    // SemiAngl, Pt, deru, derv); lin2dch.SetLocation(p2dch);`
    // (`cxx:616-627`).
    p2dch.set_coord(0.0, if plandab { chord } else { -chord });
    let (_, du, dv) = surface_eval::cone_d1(&gcon, 0.0, p2dch.y());
    let deru = du;
    let derv = dv;
    lin2dch = Geom2dLine::from_pnt_dir(p2dch, dx2d);

    // `norchamf.SetXYZ(deru.Crossed(derv).XYZ());` (`cxx:629`).
    norchamf = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;

    // `Pt.SetCoord(Or.X() + Rad * Dx.X(), ...); ElSLib::Parameters(Con, Pt, u,
    // v); double tol = Precision::PConfusion(); if (u >= 2 * M_PI - tol && u <=
    // 2 * M_PI) { u = 0.; } if (u >= fu - tol && u < fu) { u = fu; } if (u <= lu
    // + tol && u > lu) { u = lu; } if (u < fu || u > lu) { u =
    // ElCLib::InPeriod(u, fu, fu + 2 * M_PI); }` (`cxx:631-648`).
    let pt_con = offset(&or, rad, &dx);
    let (mut u, v) =
        slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &pt_con);
    let tol = PCONFUSION;
    if u >= 2.0 * std::f64::consts::PI - tol && u <= 2.0 * std::f64::consts::PI {
        u = 0.0;
    }
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    if u < fu || u > lu {
        u = clib::in_period(u, fu, fu + 2.0 * std::f64::consts::PI);
    }

    // `ElSLib::D1(u, v, Con, Pt, deru, derv); gp_Pnt2d p2dCon(u, v); gp_Dir2d
    // d2dCon; if (deru.Dot(DSp) <= 0.) { d2dCon = -gp::DX2d(); } else { d2dCon =
    // gp::DX2d(); }` (`cxx:650-660`).
    let (_pt, deru_con, derv_con) = surface_eval::cone_d1(cone, u, v);
    let mut d2d_con = dx2d;
    if deru_con.dot(&d_sp) <= 0.0 {
        d2d_con.reverse();
    }
    let lin2d_con = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), d2d_con);

    // `gp_Dir norcon = deru.Crossed(derv); gp_Dir DirCon = (Con.Axis())
    // .Direction(); gp_Dir DirChamf = (gcon->Axis()).Direction(); if (angCon >
    // Precision::Confusion()) { DirCon.Reverse(); } if (SemiAngl >
    // Precision::Confusion()) { DirChamf.Reverse(); } bool torevcon =
    // (norcon.Dot(DirCon) > 0.); bool torevcha = (norchamf.Dot(DirChamf) > 0.);
    // toreverse = ((torevcon && !torevcha) || (!torevcon && torevcha)); if
    // ((toreverse && plandab) || (!toreverse && !plandab)) { trans =
    // TopAbs_REVERSED; } else { trans = TopAbs_FORWARD; }` (`cxx:665-690`).
    let norcon = GpDir::from_vec(&deru_con.crossed(&derv_con)).map_err(|e| e.to_string())?;
    let mut dir_con = cone.axis().direction().clone();
    let mut dir_chamf = gcon.axis().direction().clone();
    if angcon > CONFUSION {
        dir_con.reverse();
    }
    if semi_angl > CONFUSION {
        dir_chamf.reverse();
    }
    let torevcon = norcon.dot(&dir_con) > 0.0;
    let torevcha = norchamf.dot(&dir_chamf) > 0.0;
    let toreverse = (torevcon && !torevcha) || (!torevcon && torevcha);
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };

    // `if (plandab) { ...OnS2... } else { ...OnS1... }` (`cxx:692-703`).
    let circ_3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ_con));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circ_3d), ds);
    let on_face: Arc<dyn Curve2d> = Arc::new(lin2d_con);
    let on_fil: Arc<dyn Curve2d> = Arc::new(lin2dch);
    if plandab {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    } else {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(on_face), Some(on_fil));
    }

    Ok(true)
}

/// `P.SetCoord(base.X() + d * dir.X(), ...)`, the pattern of `cxx:253-256`,
/// `:314-316`, `:398`, `:537-539`, `:631-633`.
fn offset(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// `new Geom_ConicalSurface(A3, Ang, R)`, whose constructor raises
/// `Standard_ConstructionError` outside the legal angle/radius range
/// (`Geom_ConicalSurface.cxx:52-62`).
fn build_cone(ax3: &GpAx3, radius: f64, semi_angle: f64) -> Result<GpCone, String> {
    GpCone::new(ax3.clone(), radius, semi_angle).map_err(|e| e.to_string())
}

/// `gp_Vec(Dir)` — the vector of a direction, used where OCCT dots a `gp_Vec`
/// with a `gp_Dir` (`cxx:324-325`, `:547-548`).
fn vec_of(d: &GpDir) -> GpVec {
    GpVec::from_xyz(d.xyz())
}
