//! `ChFiKPart_ComputeData_ChAsymPlnCyl` — the distance/angle chamfer between a
//! plane and a cylinder, guideline a straight line. The corner surface is a
//! plane.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChAsymPlnCyl.cxx:376-655`
//! (OCCT 8.0.0). This is the `GeomAbs_Line` arm of
//! `ChFiKPart_ComputeData.cxx:501-513` / `:618-636`; the `GeomAbs_Circle` arm is
//! [`super::chamfer_asym_plane_cylinder_circle`].
//!
//! Here the second distance is solved twice: first `dis2` from the angle
//! between the tangent plane and the chord `OrSpine -> OrCyl` (`cxx:445-477`),
//! then, when `Dis` sits on the cylinder, `dis1` from the chord between the
//! resulting cylinder point and the guideline start (`cxx:501-528`).
//!
//! `OrCyl` — the point of the axis "in the plane normal to the axis containing
//! OrSpine" — is built at `cxx:411-417` from a component-wise product
//! (`gp_XYZ::Multiplied(const gp_XYZ&)` multiplies coordinate by coordinate,
//! `gp_XYZ.hxx:333-337`), not from the dot product that
//! `ChPlnCyl.cxx:415-417` uses. That is what OCCT computes, so the port keeps
//! it literally; only an axis along a coordinate direction makes the two agree.

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{GpAx1, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpTrsf, GpVec};
use occt_core::gp::trsf_ext::rotate_vector;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::chamfer_cylinder_plane_tail::ch_fi_kpart_load_chamfer_cylinder_plane;

/// `ChFiKPart_MakeChAsym(DStr, Data, Pln, Cyl, fu, lu, Or1, Or2, Dis, Angle,
/// Spine, First, Ofpl, plandab, DisOnP)` (`cxx:376-655`), `Spine` a `gp_Lin`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_ch_asym_plane_cylinder_line(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cyl: &GpCylinder,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &GpLin,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> Result<bool, String> {
    // `gp_Pnt OrSpine = ElCLib::Value(First, Spine); gp_Dir XDir =
    // Spine.Direction(); gp_Ax3 AxPln = Pln.Position(); gp_Dir NorPln =
    // AxPln.XDirection().Crossed(AxPln.YDirection()); gp_Dir NorF(NorPln); if
    // (Or1 == TopAbs_REVERSED) { NorF.Reverse(); }` (`cxx:398-408`).
    let or_spine = clib::line_value(spine, first);
    let x_dir = spine.direction();
    let ax_pln = pln.position();
    let nor_pln = ax_pln
        .x_direction()
        .crossed(ax_pln.y_direction())
        .map_err(|e| e.to_string())?;
    let mut norf = nor_pln;
    if or1 == Orientation::Reversed {
        norf.reverse();
    }

    // `gp_Ax3 AxCyl = Cyl.Position(); gp_Pnt Loc = AxCyl.Location(); gp_Vec
    // LocSp(Loc, OrSpine); gp_XYZ temp = AxCyl.Direction().XYZ(); temp =
    // temp.Multiplied(LocSp.XYZ().Multiplied(temp)); OrCyl.SetXYZ
    // ((Loc.XYZ()).Added(temp));` (`cxx:410-417`).
    let ax_cyl = cyl.position();
    let loc = ax_cyl.location();
    let ax_dir = ax_cyl.direction();
    let loc_sp = GpVec::from_pnts(&loc, &or_spine);
    let temp = GpVec::new(
        ax_dir.x() * (loc_sp.x() * ax_dir.x()),
        ax_dir.y() * (loc_sp.y() * ax_dir.y()),
        ax_dir.z() * (loc_sp.z() * ax_dir.z()),
    );
    let or_cyl = GpPnt::new(
        loc.x() + temp.x(),
        loc.y() + temp.y(),
        loc.z() + temp.z(),
    );

    // `gp_Vec VecTranslPln, tmp; tmp = gp_Vec(OrSpine, OrCyl); if ((Or2 ==
    // TopAbs_FORWARD && Cyl.Direct()) || (Or2 == TopAbs_REVERSED &&
    // !Cyl.Direct())) { tmp.Reverse(); } VecTranslPln =
    // gp_Vec(XDir.Crossed(NorPln)); if (VecTranslPln.Dot(tmp) <= 0.) {
    // VecTranslPln.Reverse(); }` (`cxx:419-434`).
    let mut tmp = GpVec::from_pnts(&or_spine, &or_cyl);
    if (or2 == Orientation::Forward && cyl.is_direct())
        || (or2 == Orientation::Reversed && !cyl.is_direct())
    {
        tmp.reverse();
    }
    let mut vec_transl_pln = GpVec::from_xyz(x_dir.crossed(&nor_pln)?.xyz());
    if vec_transl_pln.dot(&tmp) <= 0.0 {
        vec_transl_pln.reverse();
    }

    // `gp_Vec VecTranslCyl; VecTranslCyl = gp_Vec(OrSpine, OrCyl);` (`cxx:436-437`).
    let vec_transl_cyl = GpVec::from_pnts(&or_spine, &or_cyl);

    // `gp_Vec DirSOrC = VecTranslCyl.Normalized(); double cosA1 =
    // DirSOrC.Dot(VecTranslPln.Normalized()); double sinA1 = std::sqrt(1. -
    // cosA1 * cosA1); double dis1 = 0.; double dis2, ray = Cyl.Radius(); bool
    // IsDisOnP = ((plandab && DisOnP) || (!plandab && !DisOnP));`
    // (`cxx:439-445`).
    let dir_s_or_c = vec_transl_cyl.normalized();
    let cos_a1 = dir_s_or_c.dot(&vec_transl_pln.normalized());
    let sin_a1 = (1.0 - cos_a1 * cos_a1).sqrt();
    let mut dis1 = 0.0;
    let dis2;
    let ray = cyl.radius();
    let is_dis_on_p = (plandab && dis_on_p) || (!plandab && !dis_on_p);

    // `if (IsDisOnP) { dis1 = Dis; double sinAl = std::sin(Angle), cosAl =
    // std::cos(Angle); double h = dis1 * sinAl; double cosAhOC = cosA1 * sinAl
    // + sinA1 * cosAl; double sinAhOC = sinA1 * sinAl - cosA1 * cosAl; if
    // (cosA1 > 0) { sinAhOC = -sinAhOC; } double temp1 = h / ray, temp2 =
    // sinAhOC * sinAhOC + temp1 * cosAhOC; dis2 = temp2 + temp1 * (cosAhOC -
    // temp1); if (dis2 < -1.E-09) { ... return false; } else if (dis2 <
    // 1.E-09) { dis2 = ray * std::sqrt(2. * temp2); } else { dis2 = ray *
    // std::sqrt(2. * (temp2 - sinAhOC * std::sqrt(dis2))); } } else { dis2 =
    // Dis; }` (`cxx:447-478`). The `dis2 < -1.E-09` arm's `std::cout` is
    // `OCCT_DEBUG` only.
    if is_dis_on_p {
        dis1 = dis;
        let sin_al = angle.sin();
        let cos_al = angle.cos();
        let h = dis1 * sin_al;
        let cos_ah_oc = cos_a1 * sin_al + sin_a1 * cos_al;
        let mut sin_ah_oc = sin_a1 * sin_al - cos_a1 * cos_al;
        if cos_a1 > 0.0 {
            sin_ah_oc = -sin_ah_oc;
        }
        let temp1 = h / ray;
        let temp2 = sin_ah_oc * sin_ah_oc + temp1 * cos_ah_oc;
        let solved = temp2 + temp1 * (cos_ah_oc - temp1);
        if solved < -1.0e-9 {
            return Ok(false);
        } else if solved < 1.0e-9 {
            dis2 = ray * (2.0 * temp2).sqrt();
        } else {
            dis2 = ray * (2.0 * (temp2 - sin_ah_oc * solved.sqrt())).sqrt();
        }
    } else {
        dis2 = dis;
    }

    // `double alpha = (2 * std::asin(dis2 * 0.5 / ray)); gp_Vec VecTemp =
    // VecTranslCyl.Reversed(); if ((XDir.Crossed(gp_Dir(VecTranslCyl))).Dot
    // (NorF) < 0.) { VecTemp.Rotate(gp_Ax1(OrCyl, XDir), alpha); } else {
    // VecTemp.Rotate(gp_Ax1(OrCyl, XDir.Reversed()), alpha); }
    // POnCyl.SetXYZ(OrCyl.XYZ().Added(VecTemp.XYZ()));` (`cxx:481-495`).
    let alpha = 2.0 * (dis2 * 0.5 / ray).asin();
    let vec_temp = vec_transl_cyl.reversed();
    let side = x_dir
        .crossed(&GpDir::from_vec(&vec_transl_cyl).map_err(|e| e.to_string())?)?
        .dot(&norf);
    let axis = GpAx1::new(or_cyl, if side < 0.0 { x_dir } else { x_dir.reversed() });
    let vec_temp = rotate_about(&vec_temp, &axis, alpha)?;
    let p_on_cyl = GpPnt::new(
        or_cyl.x() + vec_temp.x(),
        or_cyl.y() + vec_temp.y(),
        or_cyl.z() + vec_temp.z(),
    );

    // `double UOnCyl, VOnCyl, UOnPln, VOnPln; gp_Vec DUOnCyl, DVOnCyl;
    // ElSLib::Parameters(Cyl, POnCyl, UOnCyl, VOnCyl); ElSLib::CylinderD1
    // (UOnCyl, VOnCyl, AxCyl, Cyl.Radius(), POnCyl, DUOnCyl, DVOnCyl);`
    // (`cxx:496-499`).
    let (u_on_cyl, v_on_cyl) = slib::cylinder_parameters(&ax_cyl, &p_on_cyl);
    let (_pc, du_on_cyl, dv_on_cyl) = surface_eval::cylinder_d1(cyl, u_on_cyl, v_on_cyl);

    // `if (!IsDisOnP) { gp_Vec Corde(POnCyl, OrSpine); gp_Vec TCyl =
    // DUOnCyl.Crossed(DVOnCyl); TCyl = XDir.Crossed(TCyl); TCyl.Normalize();
    // Corde.Normalize(); double cosCorTan = TCyl.Dot(Corde); double tgCorTan =
    // 1. / (cosCorTan * cosCorTan); tgCorTan = std::sqrt(tgCorTan - 1.); double
    // tgAng = tan(Angle); tgAng = (tgAng + tgCorTan) / (1. - tgAng * tgCorTan);
    // double cosA11 = dis2 / (2. * ray); double sinA11 = std::sqrt(1. -
    // cosA11 * cosA11); if (cosA1 > 0) { sinA11 = -sinA11; } dis1 = (sinA1 +
    // cosA1 * tgAng) * cosA11; dis1 -= (cosA1 - sinA1 * tgAng) * sinA11; dis1 =
    // (dis2 * tgAng) / dis1; }` (`cxx:501-528`).
    if !is_dis_on_p {
        let mut corde = GpVec::from_pnts(&p_on_cyl, &or_spine);
        let t_cyl = du_on_cyl.crossed(&dv_on_cyl);
        let t_cyl_dir = x_dir.crossed(&GpDir::from_vec(&t_cyl).map_err(|e| e.to_string())?)?;
        let mut t_cyl_vec = GpVec::from_xyz(t_cyl_dir.xyz());
        t_cyl_vec.normalize();
        corde.normalize();
        let cos_cor_tan = t_cyl_vec.dot(&corde);
        let mut tg_cor_tan = 1.0 / (cos_cor_tan * cos_cor_tan);
        tg_cor_tan = (tg_cor_tan - 1.0).sqrt();
        let mut tg_ang = angle.tan();
        tg_ang = (tg_ang + tg_cor_tan) / (1.0 - tg_ang * tg_cor_tan);
        let cos_a11 = dis2 / (2.0 * ray);
        let mut sin_a11 = (1.0 - cos_a11 * cos_a11).sqrt();
        if cos_a1 > 0.0 {
            sin_a11 = -sin_a11;
        }
        dis1 = (sin_a1 + cos_a1 * tg_ang) * cos_a11;
        dis1 -= (cos_a1 - sin_a1 * tg_ang) * sin_a11;
        dis1 = (dis2 * tg_ang) / dis1;
    }

    // `VecTranslPln.Multiply(dis1); POnPln.SetXYZ((OrSpine.XYZ()).Added
    // (VecTranslPln.XYZ()));` (`cxx:530-532`).
    let vec_transl_pln = vec_transl_pln.multiplied_scalar(dis1);
    let p_on_pln = GpPnt::new(
        or_spine.x() + vec_transl_pln.x(),
        or_spine.y() + vec_transl_pln.y(),
        or_spine.z() + vec_transl_pln.z(),
    );

    // `ElSLib::Parameters(Pln, POnPln, UOnPln, VOnPln); POnPln =
    // ElSLib::PlaneValue(UOnPln, VOnPln, AxPln);` (`cxx:534-535`).
    let (u_on_pln, v_on_pln) = slib::plane_parameters(&ax_pln, &p_on_pln);

    ch_fi_kpart_load_chamfer_cylinder_plane(
        ds,
        data,
        pln,
        cyl,
        &ax_pln,
        &ax_cyl,
        &nor_pln,
        &or_spine,
        &x_dir,
        u_on_cyl,
        v_on_cyl,
        u_on_pln,
        v_on_pln,
        ofpl,
        plandab,
        Some(is_dis_on_p),
    )
}

/// `gp_Vec::Rotate(A1, Ang)` (`cxx:486`, `:489`).
fn rotate_about(v: &GpVec, axis: &GpAx1, angle: f64) -> Result<GpVec, String> {
    let mut t = GpTrsf::identity();
    t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
    Ok(rotate_vector(&t, v))
}
