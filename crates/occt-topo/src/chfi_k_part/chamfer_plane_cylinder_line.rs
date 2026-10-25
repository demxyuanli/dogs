//! `ChFiKPart_ComputeData_ChPlnCyl` — the chamfer between a plane and a
//! cylinder, guideline a straight line, for the `ChFiDS_Sym` and
//! `ChFiDS_TwoDist` kinds. The corner surface is a plane.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnCyl.cxx:372-592`
//! (OCCT 8.0.0). This is the `GeomAbs_Line` arm of
//! `ChFiKPart_ComputeData.cxx:216-345`; the `GeomAbs_Circle` arm is
//! [`super::chamfer_plane_cylinder_circle`], and the `ChFiDS_DistAngle` kind
//! goes through [`super::chamfer_asym_plane_cylinder_line`].
//!
//! Only `ChFiDS_ClassicChamfer` is handled here (`cxx:392-395`): OCCT returns
//! `false` for the two throat modes rather than deriving their distances, which
//! keeps the `Dis1`/`Dis2` arguments below literal.
//!
//! The two touch points are built by walking from the guideline's start point:
//! `POnPln` is `Dis1` along `XDir ^ NorPln` (`cxx:424-436`), `POnCyl` is the
//! cylinder point rotated off the axis foot `OrCyl` by the half-angle
//! `2 * asin(Dis2 / 2R)` (`cxx:439-462`).

use occt_core::elib::{clib, slib};
use occt_core::gp::{
    GpAx1, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpTrsf, GpVec,
};
use occt_core::gp::trsf_ext::rotate_vector;

use crate::abs::Orientation;
use crate::chfi_ds::{ChamfMode, SurfData};
use crate::topopebrep::ds::DataStructure;

use super::chamfer_cylinder_plane_tail::ch_fi_kpart_load_chamfer_cylinder_plane;

/// `ChFiKPart_MakeChamfer(DStr, Data, theMode, Pln, Cyl, fu, lu, Or1, Or2,
/// dis1, dis2, Spine, First, Ofpl, plandab)` (`cxx:372-592`), `Spine` a
/// `gp_Lin`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_chamfer_plane_cylinder_line(
    ds: &mut DataStructure,
    data: &mut SurfData,
    the_mode: ChamfMode,
    pln: &GpPln,
    cyl: &GpCylinder,
    or1: Orientation,
    or2: Orientation,
    dis1: f64,
    dis2: f64,
    spine: &GpLin,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> Result<bool, String> {
    // `if (theMode != ChFiDS_ClassicChamfer) { return false; }`
    // (`cxx:392-395`).
    if the_mode != ChamfMode::ClassicChamfer {
        return Ok(false);
    }

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

    // `gp_Ax3 AxCyl = Cyl.Position(); gp_XYZ AxLoc = AxCyl.Location().XYZ();
    // gp_XYZ AxDir = AxCyl.Direction().XYZ(); double Parameter =
    // (OrSpine.XYZ() - AxLoc) * AxDir; OrCyl.SetXYZ(AxLoc + Parameter *
    // AxDir);` (`cxx:410-418`). `OrCyl` is the point of the axis in the plane
    // normal to it that contains `OrSpine`.
    let ax_cyl = cyl.position();
    let ax_loc = ax_cyl.location();
    let ax_dir = ax_cyl.direction();
    let loc_sp = GpVec::from_pnts(&ax_loc, &or_spine);
    let parameter = loc_sp.dot(&GpVec::from_xyz(ax_dir.xyz()));
    let or_cyl = GpPnt::new(
        ax_loc.x() + parameter * ax_dir.x(),
        ax_loc.y() + parameter * ax_dir.y(),
        ax_loc.z() + parameter * ax_dir.z(),
    );

    // `gp_Vec VecTranslPln, tmp; tmp = gp_Vec(OrSpine, OrCyl); if ((Or2 ==
    // TopAbs_FORWARD && Cyl.Direct()) || (Or2 == TopAbs_REVERSED &&
    // !Cyl.Direct())) { tmp.Reverse(); } VecTranslPln =
    // gp_Vec(XDir.Crossed(NorPln)); if (VecTranslPln.Dot(tmp) <= 0.) {
    // VecTranslPln.Reverse(); } VecTranslPln.Multiply(dis1);
    // POnPln.SetXYZ((OrSpine.XYZ()).Added(VecTranslPln.XYZ()));`
    // (`cxx:420-436`).
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
    vec_transl_pln = vec_transl_pln.multiplied_scalar(dis1);
    let p_on_pln = GpPnt::new(
        or_spine.x() + vec_transl_pln.x(),
        or_spine.y() + vec_transl_pln.y(),
        or_spine.z() + vec_transl_pln.z(),
    );

    // `double alpha = (2 * std::asin(dis2 * 0.5 / Cyl.Radius())); gp_Vec
    // VecCylTransl = gp_Vec(OrCyl, OrSpine); if ((XDir.Crossed
    // (gp_Dir(VecCylTransl))).Dot(NorF) > 0.) { VecCylTransl.Rotate
    // (gp_Ax1(OrCyl, XDir), alpha); } else { VecCylTransl.Rotate
    // (gp_Ax1(OrCyl, XDir.Reversed()), alpha); }
    // POnCyl.SetXYZ(OrCyl.XYZ().Added(VecCylTransl.XYZ()));` (`cxx:439-462`).
    let alpha = 2.0 * (dis2 * 0.5 / cyl.radius()).asin();
    let vec_cyl_transl = GpVec::from_pnts(&or_cyl, &or_spine);
    let turn = x_dir
        .crossed(&GpDir::from_vec(&vec_cyl_transl).map_err(|e| e.to_string())?)?
        .dot(&norf);
    let axis = GpAx1::new(or_cyl, if turn > 0.0 { x_dir } else { x_dir.reversed() });
    let vec_cyl_transl = rotate_about(&vec_cyl_transl, &axis, alpha)?;
    let p_on_cyl = GpPnt::new(
        or_cyl.x() + vec_cyl_transl.x(),
        or_cyl.y() + vec_cyl_transl.y(),
        or_cyl.z() + vec_cyl_transl.z(),
    );

    // `ElSLib::Parameters(Cyl, POnCyl, UOnCyl, VOnCyl); POnCyl =
    // ElSLib::CylinderValue(UOnCyl, VOnCyl, AxCyl, Cyl.Radius());
    // ElSLib::Parameters(Pln, POnPln, UOnPln, VOnPln); POnPln =
    // ElSLib::PlaneValue(UOnPln, VOnPln, AxPln);` (`cxx:465-468`).
    let (u_on_cyl, v_on_cyl) = slib::cylinder_parameters(&ax_cyl, &p_on_cyl);
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
        None,
    )
}

/// `gp_Vec::Rotate(A1, Ang)` (`cxx:456`, `:458`): `gp_Vec` rotates through
/// `gp_Trsf::SetRotation` and its vectorial part (`gp_Vec.cxx`), so the port
/// builds the same rotation and applies its linear part only.
fn rotate_about(v: &GpVec, axis: &GpAx1, angle: f64) -> Result<GpVec, String> {
    let mut t = GpTrsf::identity();
    t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
    Ok(rotate_vector(&t, v))
}
