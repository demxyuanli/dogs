//! `ChFiKPart_ComputeData_ChAsymPlnPln` — the asymmetric chamfer between two
//! planes, guideline a straight line (`ChFiDS_DistAngle`).
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChAsymPlnPln.cxx:53-266`
//! (OCCT 8.0.0). `ChFiKPart_ComputeData.cxx:501-513` is the only caller: the
//! distance-angle kind passes `Dis`, `Angle` and `DisOnP = true`.
//!
//! Differs from the symmetric maker ([`super::chamfer_plane_plane`]) in two
//! places only: the `VecTransl1`/`VecTransl2` flip tests are strict `< 0.`
//! (`cxx:101`, `:107`) instead of `<= 0.`, and `Dis1`/`Dis2` come from
//! `cosP + sinP / tan(Angle)` (`cxx:112-125`) instead of `theMode`.
//! Everything from `Data->ChangeSurf` on is shared, see
//! [`super::chamfer_plane_plane_tail`].

use occt_core::elib::clib;
use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::intana::plane_plane_line;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::chamfer_plane_plane_tail::ch_fi_kpart_load_chamfer_plane_plane;

/// `ChFiKPart_MakeChAsym(DStr, Data, Pl1, Pl2, Or1, Or2, Dis, Angle, Spine,
/// First, Of1, DisOnP1)` (`cxx:53-266`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_ch_asym_plane_plane(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pl1: &GpPln,
    pl2: &GpPln,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &GpLin,
    first: f64,
    of1: Orientation,
    dis_on_p1: bool,
) -> Result<bool, String> {
    // `gp_Ax3 Pos1 = Pl1.Position(); gp_Dir D1 = ...; if (Or1 ==
    // TopAbs_REVERSED) D1.Reverse();` (`cxx:70-75`).
    let pos1 = pl1.position();
    let mut d1 = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    if or1 == Orientation::Reversed {
        d1.reverse();
    }

    // `gp_Ax3 Pos2 = Pl2.Position(); ... D2 ...` (`cxx:77-82`).
    let pos2 = pl2.position();
    let mut d2 = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        d2.reverse();
    }

    // `IntAna_QuadQuadGeo LInt(Pl1, Pl2, ...); ... Fint = ElCLib::Parameter(
    // LInt.Line(1), ElCLib::Value(First, Spine)); P = ElCLib::Value(Fint,
    // LInt.Line(1)); ... else { return false; }` (`cxx:85-97`).
    let (lint, fint, p) = match plane_plane_line(pl1, pl2, ANGULAR, CONFUSION) {
        Some(line) => {
            let fint = clib::line_parameter(&line, &clib::line_value(spine, first));
            (line, fint, clib::line_value(&line, fint))
        }
        None => return Ok(false),
    };

    // `gp_Dir LinAx1 = Spine.Direction(); gp_Dir VecTransl1 =
    // LinAx1.Crossed(D1); if (VecTransl1.Dot(D2) < 0.) { VecTransl1.Reverse();
    // }` (`cxx:99-104`).
    let lin_ax1 = spine.direction();
    let mut vt1 = lin_ax1.crossed(&d1).map_err(|e| e.to_string())?;
    if vt1.dot(&d2) < 0.0 {
        vt1.reverse();
    }

    // `gp_Dir VecTransl2 = LinAx1.Crossed(D2); if (VecTransl2.Dot(D1) < 0.) {
    // VecTransl2.Reverse(); }` (`cxx:106-110`).
    let mut vt2 = lin_ax1.crossed(&d2).map_err(|e| e.to_string())?;
    if vt2.dot(&d1) < 0.0 {
        vt2.reverse();
    }

    // `double cosP, sinP, dis1, dis2; cosP = VecTransl1.Dot(VecTransl2); sinP =
    // sqrt(1. - cosP * cosP); if (DisOnP1) { dis1 = Dis; dis2 = Dis / (cosP +
    // sinP / std::tan(Angle)); } else { dis1 = Dis / (cosP + sinP /
    // std::tan(Angle)); dis2 = Dis; }` (`cxx:112-125`).
    let cos_p = vt1.dot(&vt2);
    let sin_p = (1.0 - cos_p * cos_p).sqrt();
    let (dis1, dis2) = if dis_on_p1 {
        (dis, dis / (cos_p + sin_p / angle.tan()))
    } else {
        (dis / (cos_p + sin_p / angle.tan()), dis)
    };

    // `gp_Pnt P1(P.X() + dis1 * VecTransl1.X(), ...); gp_Pnt P2(...); gp_Pnt
    // Po((P1.X() + P2.X()) / 2., ...);` (`cxx:127-137`).
    let p1 = offset(&p, dis1, &vt1);
    let p2 = offset(&p, dis2, &vt2);
    let po = midpoint(&p1, &p2);

    // `gp_Pnt Pp = ElCLib::Value(Fint + 10., LInt.Line(1)); gp_Pnt P22(Pp.X() +
    // dis2 * VecTransl2.X(), ...);` (`cxx:140-143`).
    let pp = clib::line_value(&lint, fint + 10.0);
    let p22 = offset(&pp, dis2, &vt2);

    // `gp_Dir V1(P2 - P1); gp_Dir V2(P22 - P1); gp_Dir AxisPlan =
    // V1.Crossed(V2);` (`cxx:146-148`).
    let v1 = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).map_err(|e| e.to_string())?;
    let v2 = GpDir::from_vec(&GpVec::from_pnts(&p1, &p22)).map_err(|e| e.to_string())?;
    let axis_plan = v1.crossed(&v2).map_err(|e| e.to_string())?;

    // `gp_Dir xdir = LinAx1; gp_Ax3 PlanAx3(Po, AxisPlan, xdir); if
    // (PlanAx3.YDirection().Dot(D2) >= 0.) { PlanAx3.YReverse(); }`
    // (`cxx:150-155`).
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

/// The `P.X() + d * V.X()` pattern of `cxx:127-134`, `:141-143`.
fn offset(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}

/// The midpoint of `cxx:137`.
fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        (a.x() + b.x()) / 2.0,
        (a.y() + b.y()) / 2.0,
        (a.z() + b.z()) / 2.0,
    )
}
