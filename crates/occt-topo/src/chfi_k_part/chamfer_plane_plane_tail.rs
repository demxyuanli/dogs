//! The tail shared by the two plane/plane chamfer makers: everything from
//! `Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gpl, DStr))` to `return true`.
//!
//! `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnPln.cxx:163-271`
//! and `ChFiKPart_ComputeData_ChAsymPlnPln.cxx:157-265` (OCCT 8.0.0) are
//! textually identical in that range — the two functions differ only in how
//! `Dis1`/`Dis2` are derived (and in the `<= 0.` / `< 0.` sign tests of their
//! `VecTransl1`/`VecTransl2` flip). Keeping the common part here leaves each
//! maker's own head line-for-line comparable with its `.cxx`.
//!
//! Callers build `PlanAx3` (including its `YReverse`) and pass the values the
//! tail reads, so the argument list is the set of locals alive at `:163`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx3, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::line::GeomLine;
use occt_geom::plane::GeomPlane;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChPlnPln.cxx:163-271` / `ChAsymPlnPln.cxx:157-265`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_load_chamfer_plane_plane(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pos1: &GpAx3,
    pos2: &GpAx3,
    vt1: &GpDir,
    vt2: &GpDir,
    xdir: &GpDir,
    p1: &GpPnt,
    p2: &GpPnt,
    plan_ax3: &GpAx3,
    of1: Orientation,
) -> Result<(), String> {
    // `occ::handle<Geom_Plane> gpl = new Geom_Plane(PlanAx3);
    // Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gpl, DStr));` (`cxx:163-164`).
    let pln = GpPln::new(plan_ax3.clone());
    let surf: Arc<dyn GeomSurface> = Arc::new(GeomPlane::new(pln));
    data.set_surf(ch_fi_kpart_index_surface_in_ds(surf, ds));

    // `gp_Dir norpl = Pos1.XDirection().Crossed(Pos1.YDirection()); gp_Dir
    // norface1 = norpl; if (Of1 == TopAbs_REVERSED) norface1.Reverse();`
    // (`cxx:168-173`).
    let norpl = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    let mut norface1 = norpl;
    if of1 == Orientation::Reversed {
        norface1.reverse();
    }

    // `gp_Dir norplch = gpl->Pln().Position().XDirection().Crossed(...YDirection
    // ());` (`cxx:176`).
    let norplch = plan_ax3
        .x_direction()
        .crossed(plan_ax3.y_direction())
        .map_err(|e| e.to_string())?;

    // `gp_Dir DirCh12(gp_Vec(P1, P2)); bool toreverse = (norplch.Dot(norface1)
    // <= 0.); if (VecTransl1.Dot(DirCh12) > 0) { toreverse = !toreverse; }`
    // (`cxx:178-183`).
    let dir_ch12 = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).map_err(|e| e.to_string())?;
    let mut toreverse = norplch.dot(&norface1) <= 0.0;
    if vt1.dot(&dir_ch12) > 0.0 {
        toreverse = !toreverse;
    }
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // Face 1 (`cxx:197-231`).
    // `gp_Lin linPln(P1, xdir); occ::handle<Geom_Line> GLinPln1 = new
    // Geom_Line(linPln); ElSLib::PlaneParameters(Pos1, P1, u, v); gp_Pnt2d
    // p2dPln(u, v); gp_Dir2d dir2dPln(xdir.Dot(Pos1.XDirection()),
    // xdir.Dot(Pos1.YDirection()));` (`cxx:197-203`).
    let glin_pln1: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p1, *xdir)));
    let (u, v) = slib::plane_parameters(pos1, p1);
    let dir2d_pln = GpDir2d::new(xdir.dot(pos1.x_direction()), xdir.dot(pos1.y_direction()))
        .map_err(|e| e.to_string())?;
    let glin2d_pln1: Arc<dyn Curve2d> =
        Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), dir2d_pln));

    // `ElSLib::PlaneParameters(PlanAx3, P1, u, v); p2dPln.SetCoord(u, v);
    // lin2dPln.SetLocation(p2dPln); lin2dPln.SetDirection(gp::DX2d());`
    // (`cxx:207-210`).
    let (u, v) = slib::plane_parameters(plan_ax3, p1);
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let glin2d_plnch1: Arc<dyn Curve2d> =
        Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), dx2d));

    // `toreverse = (norplch.Dot(norpl) <= 0.); if (VecTransl1.Dot(DirCh12) > 0)
    // { toreverse = !toreverse; } if (toreverse) { trans = TopAbs_FORWARD; }
    // else { trans = TopAbs_REVERSED; }` (`cxx:214-226`).
    let mut toreverse = norplch.dot(&norpl) <= 0.0;
    if vt1.dot(&dir_ch12) > 0.0 {
        toreverse = !toreverse;
    }
    let trans = if toreverse {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let ic = ch_fi_kpart_index_curve_in_ds(Some(glin_pln1), ds);
    data.interference_on_s1_mut()
        .set_interference(ic, trans, Some(glin2d_pln1), Some(glin2d_plnch1));

    // Face 2 (`cxx:235-269`).
    // `linPln.SetLocation(P2); ... ElSLib::PlaneParameters(Pos2, P2, u, v);
    // p2dPln.SetCoord(u, v); dir2dPln.SetCoord(xdir.Dot(Pos2.XDirection()),
    // xdir.Dot(Pos2.YDirection()));` (`cxx:235-240`).
    let glin_pln2: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p2, *xdir)));
    let (u, v) = slib::plane_parameters(pos2, p2);
    let dir2d_pln = GpDir2d::new(xdir.dot(pos2.x_direction()), xdir.dot(pos2.y_direction()))
        .map_err(|e| e.to_string())?;
    let glin2d_pln2: Arc<dyn Curve2d> =
        Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), dir2d_pln));

    // `ElSLib::PlaneParameters(PlanAx3, P2, u, v); ... lin2dPln.SetDirection
    // (gp::DX2d());` (`cxx:245-248`).
    let (u, v) = slib::plane_parameters(plan_ax3, p2);
    let glin2d_plnch2: Arc<dyn Curve2d> =
        Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), dx2d));

    // `norpl = Pos2.XDirection().Crossed(Pos2.YDirection()); toreverse =
    // (norplch.Dot(norpl) <= 0.); if (VecTransl2.Dot(DirCh12) < 0) { toreverse
    // = !toreverse; }` (`cxx:251-256`).
    let norpl = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    let mut toreverse = norplch.dot(&norpl) <= 0.0;
    if vt2.dot(&dir_ch12) < 0.0 {
        toreverse = !toreverse;
    }
    let trans = if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };
    let ic = ch_fi_kpart_index_curve_in_ds(Some(glin_pln2), ds);
    data.interference_on_s2_mut()
        .set_interference(ic, trans, Some(glin2d_pln2), Some(glin2d_plnch2));

    Ok(())
}
