//! `ChFiKPart_ComputeData_FilPlnCyl` — the fillet between a plane and a
//! cylinder, guideline a straight line.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_FilPlnCyl.cxx:52-276`
//! (OCCT 8.0.0). This is the `GeomAbs_Line` arm of
//! `ChFiKPart_ComputeData.cxx:107-175`; the `GeomAbs_Circle` arm is
//! [`super::fillet_plane_cylinder_circle`].
//!
//! The corner surface is the cylinder of radius `Radius` sitting on the
//! intersection of the plane offset by `Radius` along its oriented normal and
//! the cylinder offset radially by `Radius` (`cxx:78-128`).
//!
//! `plandab` is `true` when the plane is `S1` and the cylinder `S2`
//! (`ChFiKPart_ComputeData.cxx:122`) — it decides which of the two
//! `ChFiDS_FaceInterference` gets each of the two curves, and it flips the
//! `gp_Ax3`'s X/Y roles (`cxx:146-157`).

use std::sync::Arc;

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{GpAx3, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, ANGULAR};
use occt_geom::cylinder::GeomCylinder;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::intana::{quadric_quadric_plane_cylinder, QuadricIntersection};
use occt_geom::line::GeomLine;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChFiKPart_MakeFillet(DStr, Data, Pln, Cyl, fu, lu, Or1, Or2, Radius, Spine,
/// First, Ofpl, plandab)` (`cxx:52-276`), `Spine` a `gp_Lin`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_fillet_plane_cylinder_line(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cyl: &GpCylinder,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &GpLin,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> Result<bool, String> {
    // `gp_Ax3 AxPln = Pln.Position(); gp_Dir NorPln =
    // AxPln.XDirection().Crossed(AxPln.YDirection()); gp_Dir NorF(NorPln);
    // gp_Ax3 AxCyl = Cyl.Position();` (`cxx:69-72`).
    let ax_pln = pln.position();
    let nor_pln = ax_pln
        .x_direction()
        .crossed(ax_pln.y_direction())
        .map_err(|e| e.to_string())?;
    let mut norf = nor_pln;
    let ax_cyl = cyl.position();

    // `if (Or1 == TopAbs_REVERSED) { NorF.Reverse(); }` (`cxx:74-77`).
    if or1 == Orientation::Reversed {
        norf.reverse();
    }

    // `gp_Pln PlanOffset = Pln.Translated(Radius * gp_Vec(NorF));` (`cxx:78`).
    let plan_offset = pln.translated_vec(&GpVec::from_xyz(norf.xyz()).multiplied_scalar(radius));

    // `double ROff = Cyl.Radius(); if ((Or2 == TopAbs_FORWARD && Cyl.Direct())
    // || (Or2 == TopAbs_REVERSED && !Cyl.Direct())) { ROff += Radius; } else if
    // (Radius < ROff) { ROff -= Radius; } else { return false; }` (`cxx:81-97`).
    // The `else` arm's `std::cout` is `OCCT_DEBUG` only.
    let mut roff = cyl.radius();
    if (or2 == Orientation::Forward && cyl.is_direct())
        || (or2 == Orientation::Reversed && !cyl.is_direct())
    {
        roff += radius;
    } else if radius < roff {
        roff -= radius;
    } else {
        return Ok(false);
    }

    // `gp_Cylinder CylOffset(Cyl.Position(), ROff); IntAna_QuadQuadGeo
    // LInt(PlanOffset, CylOffset, Precision::Angular(),
    // Precision::Confusion());` (`cxx:99-100`). The `H` argument of the
    // `IntAna_QuadQuadGeo` plane/cylinder overload defaults to `0`
    // (`IntAna_QuadQuadGeo.hxx:94`), which the port passes explicitly.
    let cyl_offset = GpCylinder::new(ax_cyl.clone(), roff).map_err(|e| e.to_string())?;
    let lint = quadric_quadric_plane_cylinder(&plan_offset, &cyl_offset, ANGULAR, CONFUSION, 0.0);

    // `gp_Pnt OrSpine = ElCLib::Value(First, Spine);` (`cxx:101`).
    let or_spine = clib::line_value(spine, first);

    // `if (LInt.IsDone()) { DirFillet = LInt.Line(1).Direction(); gp_Pnt P1 =
    // ElCLib::Value(ElCLib::Parameter(LInt.Line(1), OrSpine), LInt.Line(1));
    // if (LInt.NbSolutions() == 2) { gp_Pnt P2 = ...; if
    // (P1.SquareDistance(OrSpine) < P2.SquareDistance(OrSpine)) { OrFillet = P1;
    // } else { OrFillet = P2; } } else { OrFillet = P1; } } else { return
    // false; }` (`cxx:104-128`).
    //
    // The non-parallel plane/cylinder answers are an ellipse or a
    // hyperbola/parabola, where OCCT's `LInt.Line(1)` raises
    // `Standard_OutOfRange`; the port answers `Err` for that case.
    let (dir_fillet, or_fillet) = match lint {
        QuadricIntersection::Line(l) => {
            let p1 = projected_on_line(&l, &or_spine);
            (l.direction(), p1)
        }
        QuadricIntersection::TwoLines(l1, l2) => {
            let p1 = projected_on_line(&l1, &or_spine);
            let p2 = projected_on_line(&l2, &or_spine);
            let closer = if p1.square_distance(&or_spine) < p2.square_distance(&or_spine) {
                p1
            } else {
                p2
            };
            (l1.direction(), closer)
        }
        QuadricIntersection::None => return Ok(false),
        _ => return Err("ChFiKPart_MakeFillet: LInt.Line(1)".into()),
    };

    // `if (DirFillet.Dot(Spine.Direction()) < 0.) { DirFillet.Reverse(); }`
    // (`cxx:131-134`).
    let mut dir_fillet = dir_fillet;
    if dir_fillet.dot(&spine.direction()) < 0.0 {
        dir_fillet.reverse();
    }

    // `double UOnCyl, VOnCyl, UOnPln, VOnPln; ElSLib::Parameters(Cyl, OrFillet,
    // UOnCyl, VOnCyl); constexpr double tesp = Precision::Confusion(); if
    // (UOnCyl < fu - tesp || UOnCyl > lu + tesp) { UOnCyl =
    // ElCLib::InPeriod(UOnCyl, fu, fu + 2 * M_PI); } ElSLib::Parameters(Pln,
    // OrFillet, UOnPln, VOnPln);` (`cxx:136-143`).
    let (mut u_on_cyl, v_on_cyl) = slib::cylinder_parameters(&ax_cyl, &or_fillet);
    let tesp = CONFUSION;
    if u_on_cyl < fu - tesp || u_on_cyl > lu + tesp {
        u_on_cyl = clib::in_period(u_on_cyl, fu, fu + 2.0 * std::f64::consts::PI);
    }
    let (u_on_pln, v_on_pln) = slib::plane_parameters(&ax_pln, &or_fillet);

    // `gp_Vec XDir, OtherDir; if (plandab) { XDir = NorF.Reversed(); OtherDir =
    // gp_Vec(OrFillet, ElSLib::Value(UOnCyl, VOnCyl, Cyl)); OtherDir.Normalize();
    // } else { OtherDir = NorF.Reversed(); XDir = gp_Vec(OrFillet,
    // ElSLib::Value(UOnCyl, VOnCyl, Cyl)); XDir.Normalize(); }` (`cxx:145-157`).
    let towards = GpVec::from_pnts(&or_fillet, &slib::cylinder_value(cyl, u_on_cyl, v_on_cyl));
    let norf_reversed = GpVec::from_xyz(norf.reversed().xyz());
    let (xdir, other_dir) = if plandab {
        let mut other = towards.clone();
        other.normalize();
        (norf_reversed, other)
    } else {
        let mut x = towards.clone();
        x.normalize();
        (x, norf_reversed)
    };

    // `gp_Ax3 AxFil(OrFillet, DirFillet, XDir); gp_Vec aProd =
    // XDir.Crossed(OtherDir); if (aProd.Dot(DirFillet) < 0.) { AxFil.YReverse();
    // }` (`cxx:159-164`).
    let x_dir = GpDir::from_vec(&xdir).map_err(|e| e.to_string())?;
    let mut ax_fil =
        GpAx3::from_point_dir_x(or_fillet, dir_fillet, &x_dir).map_err(|e| e.to_string())?;
    if xdir.crossed(&other_dir).dot(&GpVec::from_xyz(dir_fillet.xyz())) < 0.0 {
        ax_fil.y_reverse();
    }

    // `occ::handle<Geom_CylindricalSurface> Fillet = new
    // Geom_CylindricalSurface(AxFil, Radius); Data->ChangeSurf(
    // ChFiKPart_IndexSurfaceInDS(Fillet, DStr));` (`cxx:166-167`).
    let fillet_cyl = GpCylinder::new(ax_fil.clone(), radius).map_err(|e| e.to_string())?;
    data.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::new(GeomCylinder::new(fillet_cyl.clone())),
        ds,
    ));

    // Planar face (`cxx:169-225`).
    // `gp_Pnt2d PPln2d(UOnPln, VOnPln); gp_Dir2d VPln2d(DirFillet.Dot(
    // AxPln.XDirection()), DirFillet.Dot(AxPln.YDirection())); gp_Lin2d
    // Lin2dPln(PPln2d, VPln2d); gp_Pnt POnPln = ElSLib::Value(UOnPln, VOnPln,
    // Pln); gp_Lin C3d(POnPln, DirFillet);` (`cxx:171-175`).
    let p_pln2d = GpPnt2d::new(u_on_pln, v_on_pln);
    let v_pln2d = GpDir2d::new(
        dir_fillet.dot(ax_pln.x_direction()),
        dir_fillet.dot(ax_pln.y_direction()),
    )
    .map_err(|e| e.to_string())?;
    let lin2d_pln = Geom2dLine::from_pnt_dir(p_pln2d, v_pln2d);
    let p_on_pln = slib::plane_value(pln, u_on_pln, v_on_pln);

    // `double UOnFillet, V; ElSLib::CylinderParameters(AxFil, Radius, POnPln,
    // UOnFillet, V); if (UOnFillet > M_PI) { UOnFillet = 0.; } gp_Lin2d
    // LOnFillet(gp_Pnt2d(UOnFillet, V), gp::DY2d());` (`cxx:177-184`).
    let (mut u_on_fillet, v) = slib::cylinder_parameters(&ax_fil, &p_on_pln);
    if u_on_fillet > std::f64::consts::PI {
        u_on_fillet = 0.0;
    }
    let dy2d = GpDir2d::new(0.0, 1.0).map_err(|e| e.to_string())?;
    let l_on_fillet = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_on_fillet, v), dy2d);

    // `ElSLib::CylinderD1(UOnFillet, V, AxFil, Radius, P, deru, derv); gp_Dir
    // NorFil(deru.Crossed(derv)); bool toreverse = (NorFil.Dot(NorPln) <= 0.);
    // if (toreverse) { Data->ChangeOrientation() = TopAbs::Reverse(Ofpl); }
    // else { Data->ChangeOrientation() = Ofpl; }` (`cxx:190-201`).
    let (_, deru, derv) = surface_eval::cylinder_d1(&fillet_cyl, u_on_fillet, v);
    let nor_fil = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    let toreverse = nor_fil.dot(&nor_pln) <= 0.0;
    data.set_orientation(if toreverse { ofpl.reversed() } else { ofpl });

    // `if ((toreverse && plandab) || (!toreverse && !plandab)) { trans =
    // TopAbs_REVERSED; } else { trans = TopAbs_FORWARD; }` (`cxx:203-211`).
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };
    let l3d: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p_on_pln, dir_fillet)));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(l3d), ds);
    let l_fac: Arc<dyn Curve2d> = Arc::new(lin2d_pln);
    let l_fil: Arc<dyn Curve2d> = Arc::new(l_on_fillet);
    if plandab {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    } else {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    }

    // Cylindrical face (`cxx:227-274`).
    // `gp_Pnt2d PCyl2d(UOnCyl, VOnCyl); gp_Dir2d DPC = gp::DY2d(); if
    // (DirFillet.Dot(AxCyl.Direction()) < 0.) { DPC.Reverse(); } gp_Lin2d
    // Lin2dCyl(PCyl2d, DPC); gp_Pnt POnCyl = ElSLib::Value(UOnCyl, VOnCyl,
    // Cyl);` (`cxx:228-235`).
    let mut dpc = dy2d;
    if dir_fillet.dot(&ax_cyl.direction()) < 0.0 {
        dpc.reverse();
    }
    let p_on_cyl = slib::cylinder_value(cyl, u_on_cyl, v_on_cyl);

    // `ElSLib::CylinderParameters(AxFil, Radius, POnCyl, UOnFillet, V); if
    // (UOnFillet > M_PI) { UOnFillet = 0.; } LOnFillet = gp_Lin2d(
    // gp_Pnt2d(UOnFillet, V), gp::DY2d());` (`cxx:237-242`).
    let (mut u_on_fillet2, v2) = slib::cylinder_parameters(&ax_fil, &p_on_cyl);
    if u_on_fillet2 > std::f64::consts::PI {
        u_on_fillet2 = 0.0;
    }
    let l_on_fillet2 = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_on_fillet2, v2), dy2d);

    // `ElSLib::CylinderD1(UOnFillet, V, AxFil, Radius, P, deru, derv); NorFil =
    // gp_Dir(deru.Crossed(derv)); ElSLib::CylinderD1(UOnCyl, VOnCyl, AxCyl,
    // Cyl.Radius(), P, deru, derv); gp_Dir NorCyl(deru.Crossed(derv));` 
    // (`cxx:247-250`).
    let (_, deru, derv) = surface_eval::cylinder_d1(&fillet_cyl, u_on_fillet2, v2);
    let nor_fil = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    let (_, deru, derv) = surface_eval::cylinder_d1(cyl, u_on_cyl, v_on_cyl);
    let nor_cyl = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;

    // `toreverse = (NorFil.Dot(NorCyl) <= 0.); if ((toreverse && plandab) ||
    // (!toreverse && !plandab)) { trans = TopAbs_FORWARD; } else { trans =
    // TopAbs_REVERSED; }` (`cxx:252-260`).
    let toreverse = nor_fil.dot(&nor_cyl) <= 0.0;
    let trans = if (toreverse && plandab) || (!toreverse && !plandab) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let l3d: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p_on_cyl, dir_fillet)));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(l3d), ds);
    let l_fac: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(u_on_cyl, v_on_cyl),
        dpc,
    ));
    let l_fil: Arc<dyn Curve2d> = Arc::new(l_on_fillet2);
    if plandab {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    } else {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    }

    Ok(true)
}

/// `ElCLib::Value(ElCLib::Parameter(L, P), L)` (`cxx:107`, `:110`).
fn projected_on_line(l: &GpLin, p: &GpPnt) -> GpPnt {
    clib::line_value(l, clib::line_parameter(l, p))
}
