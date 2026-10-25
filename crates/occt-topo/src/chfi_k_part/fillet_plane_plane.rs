//! `ChFiKPart_ComputeData_FilPlnPln` — the fillet between two planes.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_FilPlnPln.cxx:42-174`
//! (OCCT 8.0.0).
//!
//! The corner surface is the cylinder of radius `Radius` whose axis is the
//! guideline `Spine` and whose origin is the intersection of the two planes
//! (`cxx:68-94`). The two `ChFiDS_FaceInterference` it fills carry, for each
//! plane, the 3D line where the cylinder meets it and the two 2D lines that
//! line has on the plane (`v = 0`, `cxx:118-143`) and on the cylinder
//! (`u = angle`, `cxx:144-171`).

use std::sync::Arc;

use occt_core::elib::{slib, surface_eval};
use occt_core::elib::clib;
use occt_core::gp::{GpAx3, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::cylinder::GeomCylinder;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::line::GeomLine;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use occt_geom::intana::plane_plane_line;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `ChFiKPart_MakeFillet(DStr, Data, Pl1, Pl2, Or1, Or2, Radius, Spine, First,
/// Of1)` (`cxx:42-174`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_fillet_plane_plane(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pl1: &GpPln,
    pl2: &GpPln,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &GpLin,
    first: f64,
    of1: Orientation,
) -> Result<bool, String> {
    // `gp_Ax3 Pos1 = Pl1.Position(); gp_Dir D1 = Pos1.XDirection().Crossed(
    // Pos1.YDirection()); if (Or1 == TopAbs_REVERSED) D1.Reverse();`
    // (`cxx:55-61`).
    let pos1 = pl1.position();
    let mut d1 = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    if or1 == Orientation::Reversed {
        d1.reverse();
    }

    // `gp_Ax3 Pos2 = Pl2.Position(); ... D2 ...` (`cxx:62-66`).
    let pos2 = pl2.position();
    let mut d2 = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    if or2 == Orientation::Reversed {
        d2.reverse();
    }

    // `IntAna_QuadQuadGeo LInt(Pl1, Pl2, Precision::Angular(),
    // Precision::Confusion()); gp_Pnt Pv; if (LInt.IsDone()) { Pv =
    // ElCLib::Value(ElCLib::Parameter(LInt.Line(1), ElCLib::Value(First,
    // Spine)), LInt.Line(1)); } else { return false; }` (`cxx:68-78`).
    let pv = match plane_plane_line(pl1, pl2, ANGULAR, CONFUSION) {
        Some(line) => {
            let on_spine = clib::line_value(spine, first);
            let par = clib::line_parameter(&line, &on_spine);
            clib::line_value(&line, par)
        }
        None => return Ok(false),
    };

    // `gp_Dir AxisCylinder = Spine.Direction(); double Ang = D1.Angle(D2);
    // gp_Vec V = gp_Vec(D1) + gp_Vec(D2); gp_Dir S(V); gp_Pnt C; double Fac =
    // Radius / std::cos(Ang / 2.); C.SetCoord(Pv.X() + Fac * S.X(), ...);`
    // (`cxx:79-87`).
    let axis_cylinder = spine.direction();
    let ang = d1.angle(&d2);
    let v = GpVec::from_xyz(d1.xyz()).added(&GpVec::from_xyz(d2.xyz()));
    let s = GpDir::from_vec(&v).map_err(|e| e.to_string())?;
    let c = {
        let fac = radius / (ang / 2.0).cos();
        GpPnt::new(
            pv.x() + fac * s.x(),
            pv.y() + fac * s.y(),
            pv.z() + fac * s.z(),
        )
    };

    // `gp_Dir xdir = D1.Reversed(); gp_Ax3 CylAx3(C, AxisCylinder, xdir); if
    // (CylAx3.YDirection().Dot(D2) >= 0.) { CylAx3.YReverse(); }`
    // (`cxx:88-92`). `gp_Ax3(P, N, Vx)` projects `Vx` onto the plane
    // perpendicular to `N` (`gp_Ax3.hxx:84-91`), which is
    // `GpAx3::from_point_dir_x`.
    let xdir = d1.reversed();
    let mut cyl_ax3 =
        GpAx3::from_point_dir_x(c, axis_cylinder, &xdir).map_err(|e| e.to_string())?;
    if cyl_ax3.y_direction().dot(&d2) >= 0.0 {
        cyl_ax3.y_reverse();
    }

    // `occ::handle<Geom_CylindricalSurface> gcyl = new
    // Geom_CylindricalSurface(CylAx3, Radius);
    // Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gcyl, DStr));`
    // (`cxx:93-94`).
    let cyl_geom = GpCylinder::new(cyl_ax3.clone(), radius).map_err(|e| e.to_string())?;
    data.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::new(GeomCylinder::new(cyl_geom.clone())),
        ds,
    ));

    // `ElSLib::CylinderD1(0., 0., CylAx3, Radius, P, deru, derv); gp_Dir
    // norcyl(deru.Crossed(derv));` (`cxx:97-100`).
    let (p, deru, derv) = surface_eval::cylinder_d1(&cyl_geom, 0.0, 0.0);
    let norcyl = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;

    // `gp_Dir norpl = Pos1.XDirection().Crossed(Pos1.YDirection()); gp_Dir
    // norface = norpl; if (Of1 == TopAbs_REVERSED) norface.Reverse(); bool
    // toreverse = (norcyl.Dot(norface) <= 0.); if (toreverse) {
    // Data->ChangeOrientation() = TopAbs_REVERSED; } else { ... FORWARD; }`
    // (`cxx:101-115`).
    let norpl = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    let mut norface = norpl;
    if of1 == Orientation::Reversed {
        norface.reverse();
    }
    let mut toreverse = norcyl.dot(&norface) <= 0.0;
    data.set_orientation(if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // La face 1 (`cxx:118-143`).
    let (u, vv) = slib::plane_parameters(&pos1, &p);
    let p2d_pln = GpPnt2d::new(u, vv);
    let dir2d_pln = GpDir2d::new(
        axis_cylinder.dot(pos1.x_direction()),
        axis_cylinder.dot(pos1.y_direction()),
    )
    .map_err(|e| e.to_string())?;
    let glin2d_pln1: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d_pln, dir2d_pln));
    // `gp_Lin linPln(P, AxisCylinder);`
    let glin_pln1: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, axis_cylinder)));
    // `gp_Lin2d lin2dCyl(gp_Pnt2d(0., 0.), gp::DY2d());` (`cxx:128`).
    let dy2d = GpDir2d::new(0.0, 1.0).map_err(|e| e.to_string())?;
    let glin2d_cyl1: Arc<dyn Curve2d> =
        Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dy2d));

    // `toreverse = (norcyl.Dot(norpl) <= 0.); if (toreverse) { trans =
    // TopAbs_REVERSED; } else { trans = TopAbs_FORWARD; }` (`cxx:130-138`).
    toreverse = norcyl.dot(&norpl) <= 0.0;
    let trans = if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };
    let ic = ch_fi_kpart_index_curve_in_ds(glin_pln1, ds);
    data.interference_on_s1_mut()
        .set_interference(ic, trans, Some(glin2d_pln1), Some(glin2d_cyl1));

    // La face 2 (`cxx:145-171`).
    let (p, deru, derv) = surface_eval::cylinder_d1(&cyl_geom, ang, 0.0);
    let norcyl = GpDir::from_vec(&deru.crossed(&derv)).map_err(|e| e.to_string())?;
    let norpl = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    toreverse = norcyl.dot(&norpl) <= 0.0;

    let (u, vv) = slib::plane_parameters(&pos2, &p);
    let p2d_pln = GpPnt2d::new(u, vv);
    let dir2d_pln = GpDir2d::new(
        axis_cylinder.dot(pos2.x_direction()),
        axis_cylinder.dot(pos2.y_direction()),
    )
    .map_err(|e| e.to_string())?;
    let glin2d_pln2: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d_pln, dir2d_pln));
    let glin_pln2: Arc<dyn GeomCurve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, axis_cylinder)));
    // `lin2dCyl.SetLocation(gp_Pnt2d(Ang, 0.));` (`cxx:159-160`).
    let glin2d_cyl2: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(ang, 0.0),
        dy2d,
    ));

    // `if (toreverse) { trans = TopAbs_FORWARD; } else { trans =
    // TopAbs_REVERSED; }` (`cxx:161-168`) — the sense is inverted on face 2.
    let trans = if toreverse {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let ic = ch_fi_kpart_index_curve_in_ds(glin_pln2, ds);
    data.interference_on_s2_mut()
        .set_interference(ic, trans, Some(glin2d_pln2), Some(glin2d_cyl2));

    Ok(true)
}
