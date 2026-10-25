//! `ChFiKPart_ComputeData_ChPlnCon` — the chamfer between a plane and a cone,
//! guideline a circle, for the `ChFiDS_Sym` and `ChFiDS_TwoDist` kinds.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnCon.cxx:49-166`
//! (OCCT 8.0.0). This is both cone arms of `ChFiKPart_ComputeData.cxx:216-345`
//! — like OCCT the dispatcher does not test `ctyp` for a cone, because the
//! plane/cone guideline is always a circle.
//!
//! OCCT does not build a surface here: it only turns the two distances into the
//! equivalent distance/angle pair and delegates to
//! [`super::chamfer_asym_plane_cone`] with `DisOnP == false` (`cxx:147-165`).
//! The angle follows from the section of the cone (`cxx:127-143`): with the
//! cone opening towards the material the section is a triangle, otherwise the
//! general case subtracts the cone's own slope.
//!
//! `ChPlnCon.cxx:164-405` is a commented-out copy of the surface build, kept in
//! OCCT as a reference and therefore not translated.

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{GpCirc, GpCone, GpDir, GpPln};
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::chfi_ds::{ChamfMode, SurfData};
use crate::topopebrep::ds::DataStructure;

use super::chamfer_asym_plane_cone::ch_fi_kpart_make_ch_asym_plane_cone;

/// `ChFiKPart_MakeChamfer(DStr, Data, theMode, Pln, Con, fu, lu, Or1, Or2,
/// theDis1, theDis2, Spine, First, Ofpl, plandab)` (`cxx:49-166`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_chamfer_plane_cone(
    ds: &mut DataStructure,
    data: &mut SurfData,
    the_mode: ChamfMode,
    pln: &GpPln,
    cone: &GpCone,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    the_dis1: f64,
    the_dis2: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> Result<bool, String> {
    // `double angcon = Con.SemiAngle(); double Dis1 = theDis1, Dis2 =
    // theDis2; double Alpha = M_PI / 2 - angcon; double CosHalfAlpha =
    // std::cos(Alpha / 2);` (`cxx:68-71`).
    let angcon = cone.semi_angle();
    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    let alpha = std::f64::consts::FRAC_PI_2 - angcon;
    let cos_half_alpha = (alpha / 2.0).cos();

    // `if (theMode == ChFiDS_ConstThroatChamfer) { Dis1 = Dis2 = theDis1 /
    // CosHalfAlpha; } else if (theMode ==
    // ChFiDS_ConstThroatWithPenetrationChamfer) { double aDis1 =
    // std::min(theDis1, theDis2); double aDis2 = std::max(theDis1, theDis2);
    // ... Dis1 = sqrt(dis2dis2 - dis1dis1) - aDis1 * CotanAlpha; double CosBeta
    // = sqrt(1 - dis1dis1 / dis2dis2) * CosAlpha + aDis1 / aDis2 * SinAlpha;
    // double FullDist1 = aDis2 / CosBeta; Dis2 = FullDist1 - aDis1 / SinAlpha;
    // }` (`cxx:72-85`).
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
        let cos_beta =
            (1.0 - dis1dis1 / dis2dis2).sqrt() * cos_alpha + a_dis1 / a_dis2 * sin_alpha;
        let full_dist1 = a_dis2 / cos_beta;
        dis2 = full_dist1 - a_dis1 / sin_alpha;
    }

    // `double sincon = std::abs(std::sin(angcon));` (`cxx:88`).
    let sincon = angcon.sin().abs();

    // `gp_Ax3 PosPl = Pln.Position(); gp_Dir Dpl =
    // PosPl.XDirection().Crossed(PosPl.YDirection()); if (Or1 ==
    // TopAbs_REVERSED) { Dpl.Reverse(); }` (`cxx:93-98`).
    let pos_pl = pln.position();
    let mut dpl = pos_pl
        .x_direction()
        .crossed(pos_pl.y_direction())
        .map_err(|e| e.to_string())?;
    if or1 == Orientation::Reversed {
        dpl.reverse();
    }

    // `gp_Pnt Or = Con.Location(); double u, v; ElSLib::PlaneParameters(PosPl,
    // Or, u, v); ElSLib::PlaneD0(u, v, PosPl, Or);` (`cxx:101-107`). The `u`,
    // `v` and the projected `Or` are only read back by the `OCCT_DEBUG` block
    // at `cxx:105-106`; `Or` itself stays the cone's own location.
    let or = cone.location();
    let (_u, _v) = slib::plane_parameters(&pos_pl, &or);

    // `gp_Pnt PtSp; gp_Vec DSp; ElCLib::D1(First, Spine, PtSp, DSp);`
    // (`cxx:109-111`).
    let (pt_sp, _d_sp) = clib::circle_d1(spine, first);

    // `gp_Vec deru, derv; gp_Pnt PtCon; ElSLib::Parameters(Con, PtSp, u, v);
    // ElSLib::D1(u, v, Con, PtCon, deru, derv); gp_Dir
    // Dcon(deru.Crossed(derv)); if (Or2 == TopAbs_REVERSED) { Dcon.Reverse();
    // }` (`cxx:114-124`).
    let (u, v) =
        slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &pt_sp);
    let (_pt_con, deru, derv) = surface_eval::cone_d1(cone, u, v);
    let mut dcon = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dcon.reverse();
    }

    // `bool ouvert = (Dpl.Dot(Dcon) >= 0.);` (`cxx:126`).
    let ouvert = dpl.dot(&dcon) >= 0.0;

    // `if (!ouvert) { if (std::abs(Dis1 - Dis2 * sincon) > Precision::Confusion())
    // { double abscos = std::abs(Dis2 - Dis1 * sincon); angle = std::atan((Dis1
    // * std::cos(angcon)) / abscos); } else { angle = angcon; } } else { angle =
    // std::atan((Dis1 * std::cos(angcon)) / (Dis2 + Dis1 * sincon)); }`
    // (`cxx:128-143`).
    let angle = if !ouvert {
        if (dis1 - dis2 * sincon).abs() > CONFUSION {
            let abscos = (dis2 - dis1 * sincon).abs();
            ((dis1 * angcon.cos()) / abscos).atan()
        } else {
            angcon
        }
    } else {
        ((dis1 * angcon.cos()) / (dis2 + dis1 * sincon)).atan()
    };

    // `bool DisOnP = false; IsResol = ChFiKPart_MakeChAsym(DStr, Data, Pln, Con,
    // fu, lu, Or1, Or2, Dis2, angle, Spine, First, Ofpl, plandab, DisOnP);`
    // (`cxx:145-165`).
    ch_fi_kpart_make_ch_asym_plane_cone(
        ds, data, pln, cone, fu, lu, or1, or2, dis2, angle, spine, first, ofpl, plandab, false,
    )
}
