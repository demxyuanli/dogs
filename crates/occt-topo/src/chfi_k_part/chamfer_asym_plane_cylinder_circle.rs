//! `ChFiKPart_ComputeData_ChAsymPlnCyl` — the distance/angle chamfer between a
//! plane and a cylinder, guideline a circle. The corner surface is a cone.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChAsymPlnCyl.cxx:65-369`
//! (OCCT 8.0.0). This is the `GeomAbs_Circle` arm of
//! `ChFiKPart_ComputeData.cxx:501-513` / `:618-636`; the `GeomAbs_Line` arm is
//! [`super::chamfer_asym_plane_cylinder_line`], and the plane/plane distance
//! angle maker is [`super::chamfer_asym_plane_plane`].
//!
//! Unlike the `ChFiDS_Sym`/`TwoDist` maker, `Dis` is the given distance on one
//! of the two faces and the other distance is solved from the angle between the
//! plane normal and the cylinder normal at the guideline's start point
//! (`cxx:124-141`): `dis = Dis / (sin(N) / tan(Angle) - cos(N))`. `DisOnP`
//! says which of the two faces `Dis` belongs to, together with `plandab`.
//!
//! Everything from `// compute the parameters of the conical surface` on is
//! shared with the symmetric maker and lives in
//! [`super::chamfer_cylinder_cone_tail`].

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{GpCirc, GpCylinder, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::chamfer_cylinder_cone_tail::ch_fi_kpart_load_chamfer_cylinder_cone;

/// `ChFiKPart_MakeChAsym(DStr, Data, Pln, Cyl, fu, lu, Or1, Or2, Dis, Angle,
/// Spine, First, Ofpl, plandab, DisOnP)` (`cxx:65-369`), `Spine` a `gp_Circ`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_ch_asym_plane_cylinder_circle(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cyl: &GpCylinder,
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
    // (`cxx:85-95`).
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
    // gp_Pnt PtPl = Or;` (`cxx:97-102`). `PlaneD0` writes the projection back
    // into `Or`, so every later use of `Or` in the `.cxx` is the projected
    // point `PtPl`, not the cylinder origin.
    let or_cyl = cyl.location();
    let (u0, v0) = slib::plane_parameters(&pos_pl, &or_cyl);
    let pt2d_pln = GpPnt2d::new(u0, v0);
    let pt_pl = slib::plane_value(pln, u0, v0);

    // `gp_Pnt PtSp; gp_Vec DSp; ElCLib::D1(First, Spine, PtSp, DSp); gp_Dir
    // Dx(gp_Vec(Or, PtSp)); gp_Dir Dy(DSp);` (`cxx:105-109`), with `Or`
    // already the projection.
    let (pt_sp, d_sp) = clib::circle_d1(spine, first);
    let dx = GpDir::from_vec(&GpVec::from_pnts(&pt_pl, &pt_sp)).map_err(|e| e.to_string())?;
    let dy = GpDir::from_vec(&d_sp).map_err(|e| e.to_string())?;

    // `ElSLib::Parameters(Cyl, PtSp, u, v); gp_Pnt PtCyl; gp_Vec Vu, Vv;
    // ElSLib::D1(u, v, Cyl, PtCyl, Vu, Vv); gp_Dir Dcyl(Vu.Crossed(Vv)); if
    // (Or2 == TopAbs_REVERSED) { Dcyl.Reverse(); } bool dedans = (Dcyl.Dot(Dx)
    // <= 0.);` (`cxx:110-119`).
    let (u, v) = slib::cylinder_parameters(&cyl.position(), &pt_sp);
    let (_pt_cyl, vu, vv) = surface_eval::cylinder_d1(cyl, u, v);
    let mut dcyl = GpDir::from_vec(&vu.crossed(&vv)).map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        dcyl.reverse();
    }
    let dedans = dcyl.dot(&dx) <= 0.0;

    // `double dis1, dis2, cosNPCyl, sinNPCyl; if ((plandab && DisOnP) ||
    // (!plandab && !DisOnP)) { dis1 = Dis; cosNPCyl = Dpl.Dot(Dcyl); sinNPCyl =
    // sqrt(1. - cosNPCyl * cosNPCyl); dis2 = Dis / (sinNPCyl / tan(Angle) -
    // cosNPCyl); } else { dis2 = Dis; ... dis1 = Dis / (sinNPCyl / tan(Angle) -
    // cosNPCyl); }` (`cxx:124-141`).
    let cos_np_cyl = dpl.dot(&dcyl);
    let sin_np_cyl = (1.0 - cos_np_cyl * cos_np_cyl).sqrt();
    let other = dis / (sin_np_cyl / angle.tan() - cos_np_cyl);
    let (dis1, dis2) = if (plandab && dis_on_p) || (!plandab && !dis_on_p) {
        (dis, other)
    } else {
        (other, dis)
    };

    // `Or.SetCoord(Or.X() + dis2 * Dpl.X(), ...);` (`cxx:143`), from the
    // projected `Or`.
    let or_shifted = offset(&pt_pl, dis2, &dpl);

    // `gp_Dir Vec1(Or.X() - PtPl.X(), ...); gp_Pnt Pt(Or.X() + dis1 *
    // PosPl.XDirection().X(), ...); gp_Dir Vec2(Pt.X() - PtPl.X(), ...);`
    // (`cxx:145-150`).
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

/// The `base.X() + d * dir.X()` pattern of `cxx:143`, `:145-148`.
fn offset(base: &GpPnt, d: f64, dir: &GpDir) -> GpPnt {
    GpPnt::new(
        base.x() + d * dir.x(),
        base.y() + d * dir.y(),
        base.z() + d * dir.z(),
    )
}
