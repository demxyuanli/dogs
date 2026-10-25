//! `ChFiKPart_ComputeData_ChPlnCyl` — the chamfer between a plane and a
//! cylinder, guideline a circle, for the `ChFiDS_Sym` and `ChFiDS_TwoDist`
//! kinds. The corner surface is a cone.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnCyl.cxx:64-363`
//! (OCCT 8.0.0). This is the `GeomAbs_Circle` arm of
//! `ChFiKPart_ComputeData.cxx:216-345`; the `GeomAbs_Line` arm is
//! [`super::chamfer_plane_cylinder_line`], and the `ChFiDS_DistAngle` kind goes
//! through [`super::chamfer_asym_plane_cylinder_circle`].
//!
//! `theMode` only rewrites `Dis1`/`Dis2` (`cxx:82-91`): a constant throat makes
//! both distances `theDis1 * sqrt(2)`, the penetrating variant solves the
//! section from the larger and smaller distance. The cone itself, its two
//! `FaceInterference` records and the `pointu` degenerate case are shared with
//! the asymmetric maker and live in
//! [`super::chamfer_cylinder_cone_tail`].

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{GpCirc, GpCylinder, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};

use crate::abs::Orientation;
use crate::chfi_ds::{ChamfMode, SurfData};
use crate::topopebrep::ds::DataStructure;

use super::chamfer_cylinder_cone_tail::ch_fi_kpart_load_chamfer_cylinder_cone;

/// `ChFiKPart_MakeChamfer(DStr, Data, theMode, Pln, Cyl, fu, lu, Or1, Or2,
/// theDis1, theDis2, Spine, First, Ofpl, plandab)` (`cxx:64-363`), `Spine` a
/// `gp_Circ`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_chamfer_plane_cylinder_circle(
    ds: &mut DataStructure,
    data: &mut SurfData,
    the_mode: ChamfMode,
    pln: &GpPln,
    cyl: &GpCylinder,
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
    // `double Dis1 = theDis1, Dis2 = theDis2; if (theMode ==
    // ChFiDS_ConstThroatChamfer) { Dis1 = Dis2 = theDis1 * sqrt(2.); } else if
    // (theMode == ChFiDS_ConstThroatWithPenetrationChamfer) { double aDis2 =
    // std::min(theDis1, theDis2); double aDis1 = std::max(theDis1, theDis2);
    // Dis2 = sqrt(aDis1 * aDis1 - aDis2 * aDis2); Dis1 = aDis1 * aDis1 / aDis2
    // - aDis2; }` (`cxx:82-91`).
    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    if the_mode == ChamfMode::ConstThroatChamfer {
        dis1 = the_dis1 * std::f64::consts::SQRT_2;
        dis2 = dis1;
    } else if the_mode == ChamfMode::ConstThroatWithPenetrationChamfer {
        let a_dis2 = the_dis1.min(the_dis2);
        let a_dis1 = the_dis1.max(the_dis2);
        dis2 = (a_dis1 * a_dis1 - a_dis2 * a_dis2).sqrt();
        dis1 = a_dis1 * a_dis1 / a_dis2 - a_dis2;
    }

    // `gp_Ax3 PosPl = Pln.Position(); gp_Dir Dpl = PosPl.XDirection().Crossed
    // (PosPl.YDirection()); gp_Dir norf = Dpl; if (Ofpl == TopAbs_REVERSED) {
    // norf.Reverse(); } if (Or1 == TopAbs_REVERSED) { Dpl.Reverse(); }`
    // (`cxx:94-108`).
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

    // `gp_Pnt Or = Cyl.Location(); double u, v; ElSLib::PlaneParameters(PosPl,
    // Or, u, v); gp_Pnt2d pt2dPln(u, v); ElSLib::PlaneD0(u, v, PosPl, Or);
    // gp_Pnt PtPl = Or;` (`cxx:110-115`). `PlaneD0` writes the projection back
    // into `Or`, so every later use of `Or` in the `.cxx` is the projected
    // point `PtPl`, not the cylinder origin.
    let or_cyl = cyl.location();
    let (u0, v0) = slib::plane_parameters(&pos_pl, &or_cyl);
    let pt2d_pln = GpPnt2d::new(u0, v0);
    let pt_pl = slib::plane_value(pln, u0, v0);

    // `gp_Pnt PtSp; gp_Vec DSp; ElCLib::D1(First, Spine, PtSp, DSp);
    // gp_Dir Dx(gp_Vec(Or, PtSp)); gp_Dir Dy(DSp);` (`cxx:117-122`), with `Or`
    // already the projection.
    let (pt_sp, d_sp) = clib::circle_d1(spine, first);
    let dx = GpDir::from_vec(&GpVec::from_pnts(&pt_pl, &pt_sp)).map_err(|e| e.to_string())?;
    let dy = GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?;

    // `ElSLib::Parameters(Cyl, PtSp, u, v); gp_Pnt PtCyl; gp_Vec Vu, Vv;
    // ElSLib::D1(u, v, Cyl, PtCyl, Vu, Vv); gp_Dir Dcyl(Vu.Crossed(Vv)); if
    // (Or2 == TopAbs_REVERSED) { Dcyl.Reverse(); } bool dedans =
    // (Dcyl.Dot(Dx) <= 0.);` (`cxx:123-131`).
    let (u, v) = slib::cylinder_parameters(&cyl.position(), &pt_sp);
    let (_pt_cyl, vu, vv) = surface_eval::cylinder_d1(cyl, u, v);
    let mut dcyl = GpDir::from_vec(&vu.crossed(&vv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dcyl.reverse();
    }
    let dedans = dcyl.dot(&dx) <= 0.0;

    // `Or.SetCoord(Or.X() + Dis2 * Dpl.X(), ...);` (`cxx:137`), from the
    // projected `Or`.
    let or_shifted = offset(&pt_pl, dis2, &dpl);

    // `gp_Dir Vec1(Or.X() - PtPl.X(), ...); gp_Pnt Pt(Or.X() + Dis1 *
    // PosPl.XDirection().X(), ...); gp_Dir Vec2(Pt.X() - PtPl.X(), ...);`
    // (`cxx:139-144`).
    let vec1 = GpDir::from_vec(&GpVec::from_pnts(&pt_pl, &or_shifted)).map_err(|e| e.to_string())?;
    let pt = offset(&or_shifted, dis1, pos_pl.x_direction());
    let vec2 = GpDir::from_vec(&GpVec::from_pnts(&pt_pl, &pt)).map_err(|e| e.to_string())?;

    ch_fi_kpart_load_chamfer_cylinder_cone(
        ds,
        data,
        &pos_pl,
        cyl,
        dpl,
        &norf,
        &dx,
        &dy,
        &d_sp,
        &or_shifted,
        &pt_pl,
        pt2d_pln,
        &vec1,
        &vec2,
        dedans,
        plandab,
        dis1,
        dis2,
        fu,
        lu,
    )
}

/// The `base.X() + d * dir.X()` pattern of `cxx:137`, `:139-141`.
fn offset(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}
