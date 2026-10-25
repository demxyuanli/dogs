//! The tail shared by the two plane/cylinder chamfer makers whose guideline is
//! a straight line: everything from `gp_Vec YDir(POnPln, POnCyl)` to
//! `return true`. The corner surface there is a plane.
//!
//! `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_ChPlnCyl.cxx:471-591`
//! and `ChFiKPart_ComputeData_ChAsymPlnCyl.cxx:536-653` (OCCT 8.0.0) are the
//! same statements in that range; the two makers differ only in how `POnPln`
//! and `POnCyl` are constructed, which stays in each caller's head, and in
//! whether `toreverse` is negated unconditionally or only when the distance
//! sits on the plane. That last difference is the `gate` argument: `None`
//! reproduces `ChPlnCyl.cxx:505-508` (`if (PosChamfPln) {
//! toreverse = !toreverse; }`) and `Some(is_dis_on_p)` reproduces
//! `ChAsymPlnCyl.cxx:569-572` / `:631-634`
//! (`if (!IsDisOnP && PosChamfPln)` / `if (IsDisOnP && PosChamfCyl)`).

use std::sync::Arc;

use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpAx3, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::line::GeomLine;
use occt_geom::plane::GeomPlane;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChPlnCyl.cxx:471-591` / `ChAsymPlnCyl.cxx:536-653`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_load_chamfer_cylinder_plane(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &GpPln,
    cyl: &GpCylinder,
    ax_pln: &GpAx3,
    ax_cyl: &GpAx3,
    nor_pln: &GpDir,
    or_spine: &GpPnt,
    x_dir: &GpDir,
    u_on_cyl: f64,
    v_on_cyl: f64,
    u_on_pln: f64,
    v_on_pln: f64,
    ofpl: Orientation,
    plandab: bool,
    gate: Option<bool>,
) -> Result<bool, String> {
    // `gp_Vec YDir(POnPln, POnCyl); if (!plandab) { YDir.Reverse(); } gp_Ax3
    // AxCh(POnPln, XDir.Crossed(YDir), XDir);` (`cxx:472-477`).
    let p_on_pln = slib::plane_value(pln, u_on_pln, v_on_pln);
    let p_on_cyl = slib::cylinder_value(cyl, u_on_cyl, v_on_cyl);
    let mut y_dir = GpVec::from_pnts(&p_on_pln, &p_on_cyl);
    if !plandab {
        y_dir.reverse();
    }
    let ax_ch = GpAx3::from_point_dir_x(
        p_on_pln,
        x_dir.crossed(&GpDir::from_vec(&y_dir).map_err(|e| e.to_string())?)?,
        x_dir,
    )
    .map_err(|e| e.to_string())?;

    // `occ::handle<Geom_Plane> Chamfer = new Geom_Plane(AxCh);
    // Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(Chamfer, DStr));`
    // (`cxx:479-480`).
    data.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::new(GeomPlane::new(GpPln::new(ax_ch.clone()))),
        ds,
    ));

    // `gp_Pnt2d PPln2d(UOnPln, VOnPln); gp_Dir2d VPln2d(XDir.Dot
    // (AxPln.XDirection()), XDir.Dot(AxPln.YDirection())); gp_Lin2d
    // Lin2dPln(PPln2d, VPln2d);` (`cxx:483-485`).
    let v_pln2d = GpDir2d::new(
        x_dir.dot(ax_pln.x_direction()),
        x_dir.dot(ax_pln.y_direction()),
    )
    .map_err(|e| e.to_string())?;
    let lin2d_pln = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_on_pln, v_on_pln), v_pln2d);

    // `POnPln = ElSLib::Value(UOnPln, VOnPln, Pln); gp_Lin C3d(POnPln, XDir);`
    // (`cxx:487-488`).
    let p_on_pln = slib::plane_value(pln, u_on_pln, v_on_pln);

    // `ElSLib::PlaneParameters(AxCh, POnPln, U, VOnChamfer); gp_Lin2d
    // LOnChamfer(gp_Pnt2d(U, VOnChamfer), gp::DX2d());` (`cxx:490-492`).
    let dx2d = GpDir2d::new(1.0, 0.0).map_err(|e| e.to_string())?;
    let (u, v_chamfer) = slib::plane_parameters(&ax_ch, &p_on_pln);
    let l_on_chamfer = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v_chamfer), dx2d);

    // `gp_Dir NorFil = AxCh.Direction(); bool toreverse = (NorFil.Dot(NorPln)
    // <= 0.); gp_Dir DirPlnCyl(gp_Vec(POnPln, POnCyl)); gp_Dir DirSPln
    // (gp_Vec(OrSpine, POnPln)); bool PosChamfPln = DirPlnCyl.Dot(DirSPln) > 0;`
    // (`cxx:498-503`).
    let nor_fil = ax_ch.direction();
    let toreverse = nor_fil.dot(nor_pln) <= 0.0;
    let dir_pln_cyl =
        GpDir::from_vec(&GpVec::from_pnts(&p_on_pln, &p_on_cyl)).map_err(|e| e.to_string())?;
    let dir_s_pln =
        GpDir::from_vec(&GpVec::from_pnts(or_spine, &p_on_pln)).map_err(|e| e.to_string())?;
    let pos_chamf_pln = dir_pln_cyl.dot(&dir_s_pln) > 0.0;

    // `if (PosChamfPln) { toreverse = !toreverse; }` (`cxx:505-508`) resp.
    // `if (!IsDisOnP && PosChamfPln) { ... }` (`cxx:569-572`).
    let flip = match gate {
        None => pos_chamf_pln,
        Some(is_dis_on_p) => !is_dis_on_p && pos_chamf_pln,
    };
    let toreverse = if flip { !toreverse } else { toreverse };

    // `if (toreverse) { Data->ChangeOrientation() = TopAbs::Reverse(Ofpl); }
    // else { Data->ChangeOrientation() = Ofpl; }` (`cxx:510-517`).
    data.set_orientation(if toreverse { ofpl.reversed() } else { ofpl });

    // `TopAbs_Orientation trans = TopAbs_FORWARD; if ((!plandab &&
    // toreverse) || (plandab && !toreverse)) { trans = TopAbs_REVERSED; }`
    // (`cxx:519-523`).
    let trans = if (!plandab && toreverse) || (plandab && !toreverse) {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };

    // `if (plandab) { Data->ChangeInterferenceOnS1().SetInterference(
    // IndexCurveInDS(L3d, DStr), trans, LFac, LFil); } else { ...OnS2... }`
    // (`cxx:527-539`).
    let l3d: Arc<dyn GeomCurve> =
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(p_on_pln, *x_dir)));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(l3d), ds);
    let l_fac: Arc<dyn Curve2d> = Arc::new(lin2d_pln);
    let l_fil: Arc<dyn Curve2d> = Arc::new(l_on_chamfer);
    if plandab {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    } else {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    }

    // Cylindrical face (`cxx:542-591`).
    // `gp_Pnt2d PCyl2d(UOnCyl, VOnCyl); gp_Dir2d VCyl2d = gp::DY2d(); if
    // (XDir.Dot(AxCyl.Direction()) < 0) { VCyl2d.Reverse(); } gp_Lin2d
    // Lin2dCyl(PCyl2d, VCyl2d);` (`cxx:543-549`).
    let mut v_cyl2d = GpDir2d::new(0.0, 1.0).map_err(|e| e.to_string())?;
    if x_dir.dot(&ax_cyl.direction()) < 0.0 {
        v_cyl2d.reverse();
    }
    let lin2d_cyl = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_on_cyl, v_on_cyl), v_cyl2d);

    // `POnCyl = ElSLib::Value(UOnCyl, VOnCyl, Cyl); C3d = gp_Lin(POnCyl,
    // XDir); ElSLib::PlaneParameters(AxCh, POnCyl, U, VOnChamfer); LOnChamfer
    // = gp_Lin2d(gp_Pnt2d(U, VOnChamfer), gp::DX2d());` (`cxx:551-556`).
    let p_on_cyl = slib::cylinder_value(cyl, u_on_cyl, v_on_cyl);
    let (u, v_chamfer) = slib::plane_parameters(&ax_ch, &p_on_cyl);
    let l_on_chamfer = Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v_chamfer), dx2d);

    // `ElSLib::CylinderD1(UOnCyl, VOnCyl, AxCyl, Cyl.Radius(), POnCyl, deru,
    // derv); gp_Dir NorCyl(deru.Crossed(derv));` (`cxx:561-563`).
    let (_, deru, derv) = surface_eval::cylinder_d1(cyl, u_on_cyl, v_on_cyl);
    let nor_cyl = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;

    // `toreverse = (NorFil.Dot(NorCyl) <= 0.); gp_Dir DirSCyl(gp_Vec(OrSpine,
    // POnCyl)); bool PosChamfCyl = DirPlnCyl.Dot(DirSCyl) < 0;`
    // (`cxx:565-568`, in the other order at `cxx:626-629`).
    let toreverse = nor_fil.dot(&nor_cyl) <= 0.0;
    let dir_s_cyl =
        GpDir::from_vec(&GpVec::from_pnts(or_spine, &p_on_cyl)).map_err(|e| e.to_string())?;
    let pos_chamf_cyl = dir_pln_cyl.dot(&dir_s_cyl) < 0.0;

    // `if (PosChamfCyl) { toreverse = !toreverse; }` (`cxx:569-572`) resp.
    // `if (IsDisOnP && PosChamfCyl) { ... }` (`cxx:631-634`).
    let flip = match gate {
        None => pos_chamf_cyl,
        Some(is_dis_on_p) => is_dis_on_p && pos_chamf_cyl,
    };
    let toreverse = if flip { !toreverse } else { toreverse };

    // `trans = TopAbs_REVERSED; if ((!plandab && toreverse) || (plandab &&
    // !toreverse)) { trans = TopAbs_FORWARD; }` (`cxx:574-578`).
    let trans = if (!plandab && toreverse) || (plandab && !toreverse) {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };

    // `if (plandab) { ...OnS2... } else { ...OnS1... }` (`cxx:580-592`).
    let l3d: Arc<dyn GeomCurve> =
        Arc::new(GeomLine::new(GpLin::from_pnt_dir(p_on_cyl, *x_dir)));
    let ic = ch_fi_kpart_index_curve_in_ds(Some(l3d), ds);
    let l_fac: Arc<dyn Curve2d> = Arc::new(lin2d_cyl);
    let l_fil: Arc<dyn Curve2d> = Arc::new(l_on_chamfer);
    if plandab {
        data.interference_on_s2_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    } else {
        data.interference_on_s1_mut()
            .set_interference(ic, trans, Some(l_fac), Some(l_fil));
    }

    Ok(true)
}
