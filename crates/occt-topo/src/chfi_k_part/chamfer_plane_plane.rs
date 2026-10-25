//! `ChFiKPart_ComputeData_ChPlnPln` — the chamfer between two planes, guideline
//! a straight line, for the `ChFiDS_Sym` and `ChFiDS_TwoDist` kinds.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnPln.cxx:53-272`
//! (OCCT 8.0.0). The `ChFiDS_DistAngle` kind goes through
//! [`super::chamfer_asym_plane_plane`] instead (`ChFiKPart_ComputeData.cxx:501-513`).
//!
//! The corner surface is the plane through the guideline direction `Spine`
//! and the two points offset from the plane/plane intersection line by
//! `Dis1`/`Dis2` along `LinAx1 ^ D1` and `LinAx1 ^ D2` (`cxx:98-164`). `theMode`
//! only rewrites `Dis1`/`Dis2` (`cxx:111-130`): a constant throat halves the
//! section, the penetrating variant solves the right-angled section.
//!
//! Everything from `Data->ChangeSurf` on is shared with the asymmetric maker and
//! lives in [`super::chamfer_plane_plane_tail`].

use occt_core::elib::clib;
use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::intana::plane_plane_line;

use crate::abs::Orientation;
use crate::chfi_ds::{ChamfMode, SurfData};
use crate::topopebrep::ds::DataStructure;

use super::chamfer_plane_plane_tail::ch_fi_kpart_load_chamfer_plane_plane;

/// `ChFiKPart_MakeChamfer(DStr, Data, theMode, Pl1, Pl2, Or1, Or2, theDis1,
/// theDis2, Spine, First, Of1)` (`cxx:53-272`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_chamfer_plane_plane(
    ds: &mut DataStructure,
    data: &mut SurfData,
    the_mode: ChamfMode,
    pl1: &GpPln,
    pl2: &GpPln,
    or1: Orientation,
    or2: Orientation,
    the_dis1: f64,
    the_dis2: f64,
    spine: &GpLin,
    first: f64,
    of1: Orientation,
) -> Result<bool, String> {
    // `gp_Ax3 Pos1 = Pl1.Position(); gp_Dir D1 = Pos1.XDirection().Crossed(
    // Pos1.YDirection()); if (Or1 == TopAbs_REVERSED) D1.Reverse();`
    // (`cxx:70-75`).
    let pos1 = pl1.position();
    let mut d1 = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    if or1 == Orientation::Reversed {
        d1.reverse();
    }

    // `gp_Ax3 Pos2 = Pl2.Position(); ... D2 ...` (`cxx:76-81`).
    let pos2 = pl2.position();
    let mut d2 = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        d2.reverse();
    }

    // `IntAna_QuadQuadGeo LInt(Pl1, Pl2, Precision::Angular(),
    // Precision::Confusion()); double Fint; if (LInt.IsDone()) { Fint =
    // ElCLib::Parameter(LInt.Line(1), ElCLib::Value(First, Spine)); P =
    // ElCLib::Value(Fint, LInt.Line(1)); } else { return false; }`
    // (`cxx:84-96`).
    let (lint, fint, p) = match plane_plane_line(pl1, pl2, ANGULAR, CONFUSION) {
        Some(line) => {
            let fint = clib::line_parameter(&line, &clib::line_value(spine, first));
            (line, fint, clib::line_value(&line, fint))
        }
        None => return Ok(false),
    };

    // `gp_Dir LinAx1 = Spine.Direction(); gp_Dir VecTransl1 = LinAx1.Crossed
    // (D1); if (VecTransl1.Dot(D2) <= 0.) { VecTransl1.Reverse(); }`
    // (`cxx:98-103`). The `<=` here differs from the `<` of
    // `ChAsymPlnPln.cxx:101`.
    let lin_ax1 = spine.direction();
    let mut vt1 = lin_ax1.crossed(&d1).map_err(|e| e.to_string())?;
    if vt1.dot(&d2) <= 0.0 {
        vt1.reverse();
    }

    // `gp_Dir VecTransl2 = LinAx1.Crossed(D2); if (VecTransl2.Dot(D1) <= 0.) {
    // VecTransl2.Reverse(); }` (`cxx:105-109`).
    let mut vt2 = lin_ax1.crossed(&d2).map_err(|e| e.to_string())?;
    if vt2.dot(&d1) <= 0.0 {
        vt2.reverse();
    }

    // `double Dis1 = theDis1, Dis2 = theDis2; double Alpha =
    // VecTransl1.Angle(VecTransl2); double CosHalfAlpha = std::cos(Alpha / 2);
    // if (theMode == ChFiDS_ConstThroatChamfer) { Dis1 = Dis2 = theDis1 /
    // CosHalfAlpha; } else if (theMode ==
    // ChFiDS_ConstThroatWithPenetrationChamfer) { double aDis1 =
    // std::min(theDis1, theDis2); double aDis2 = std::max(theDis1, theDis2);
    // ... Dis1 = sqrt(dis2dis2 - dis1dis1) - aDis1 * CotanAlpha; double CosBeta
    // = sqrt(1 - dis1dis1 / dis2dis2) * CosAlpha + aDis1 / aDis2 * SinAlpha;
    // double FullDist1 = aDis2 / CosBeta; Dis2 = FullDist1 - aDis1 / SinAlpha;
    // }` (`cxx:111-130`).
    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    let alpha = vt1.angle(&vt2);
    let cos_half_alpha = (alpha / 2.0).cos();
    if the_mode == ChamfMode::ConstThroatChamfer {
        dis1 = the_dis1 / cos_half_alpha;
        dis2 = dis1;
    } else if the_mode == ChamfMode::ConstThroatWithPenetrationChamfer {
        let a_dis1 = the_dis1.min(the_dis2);
        let a_dis2 = the_dis1.max(the_dis2);
        let dis1dis1 = a_dis1 * a_dis1;
        let dis2dis2 = a_dis2 * a_dis2;
        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();
        let cotan_alpha = cos_alpha / sin_alpha;
        dis1 = (dis2dis2 - dis1dis1).sqrt() - a_dis1 * cotan_alpha;
        let cos_beta = (1.0 - dis1dis1 / dis2dis2).sqrt() * cos_alpha
            + a_dis1 / a_dis2 * sin_alpha;
        let full_dist1 = a_dis2 / cos_beta;
        dis2 = full_dist1 - a_dis1 / sin_alpha;
    }

    // `gp_Pnt P1(P.X() + Dis1 * VecTransl1.X(), ...); gp_Pnt P2(P.X() + Dis2 *
    // VecTransl2.X(), ...); gp_Pnt Po((P1.X() + P2.X()) / 2., ...);`
    // (`cxx:133-143`).
    let p1 = offset(&p, dis1, &vt1);
    let p2 = offset(&p, dis2, &vt2);
    let po = midpoint(&p1, &p2);

    // `gp_Pnt Pp = ElCLib::Value(Fint + 10., LInt.Line(1)); gp_Pnt P22(Pp.X() +
    // Dis2 * VecTransl2.X(), ...);` (`cxx:146-149`).
    let pp = clib::line_value(&lint, fint + 10.0);
    let p22 = offset(&pp, dis2, &vt2);

    // `gp_Dir V1(P2 - P1); gp_Dir V2(P22 - P1); gp_Dir AxisPlan =
    // V1.Crossed(V2);` (`cxx:152-154`).
    let v1 = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).map_err(|e| e.to_string())?;
    let v2 = GpDir::from_vec(&GpVec::from_pnts(&p1, &p22)).map_err(|e| e.to_string())?;
    let axis_plan = v1.crossed(&v2).map_err(|e| e.to_string())?;

    // `gp_Dir xdir = LinAx1; gp_Ax3 PlanAx3(Po, AxisPlan, xdir); if
    // (PlanAx3.YDirection().Dot(D2) >= 0.) { PlanAx3.YReverse(); }`
    // (`cxx:156-161`).
    let xdir = lin_ax1;
    let mut plan_ax3 =
        GpAx3::from_point_dir_x(po, axis_plan, &xdir).map_err(|e| e.to_string())?;
    if plan_ax3.y_direction().dot(&d2) >= 0.0 {
        plan_ax3.y_reverse();
    }

    ch_fi_kpart_load_chamfer_plane_plane(
        ds, data, &pos1, &pos2, &vt1, &vt2, &xdir, &p1, &p2, &plan_ax3, of1,
    )?;
    Ok(true)
}

/// The `P.X() + d * V.X()` pattern of `cxx:133-140`, `:147-149`.
fn offset(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// The midpoint of `cxx:143`.
fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        (a.x() + b.x()) / 2.0,
        (a.y() + b.y()) / 2.0,
        (a.z() + b.z()) / 2.0,
    )
}
