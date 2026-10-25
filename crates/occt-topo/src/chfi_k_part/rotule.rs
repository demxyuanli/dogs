//! `ChFiKPart_ComputeData_Rotule` — the "rotule" corner between three planes.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_Rotule.cxx:43-173`
//! (OCCT 8.0.0), `bool ChFiKPart_MakeRotule(DStr, Data, pl, pl1, pl2, opl,
//! opl1, opl2, r, ofpl)`.
//!
//! The surface is a torus of both radii `r` whose axis passes through the point
//! of the edge line `pl1 ∩ pl2` at distance `r` from `pl`, along the reversed
//! normal of `pl` (`cxx:83-101`). The `S1` p-curve is the circle of radius `r`
//! around that axis, projected onto `pl`; the pointed side (`S2`) is the
//! degenerate `v = PI` isoline, so it takes the null 3D curve and null p-curve
//! that `ChFiKPart_ComputeData.cxx:735-737` fills in.

use std::sync::Arc;

use occt_core::elib::{clib, slib, surface_eval};
use occt_core::gp::{
    GpAx2, GpAx3, GpAx22d, GpCirc, GpCirc2d, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d, GpTorus,
};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::intana::plane_plane_line;
use occt_geom::torus::GeomTorus;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;

use super::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};

/// `gp::DX2d()` (`gp.hxx`) — the 2D X axis direction.
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };

/// `ChFiKPart_MakeRotule(DStr, Data, pl, pl1, pl2, opl, opl1, opl2, r, ofpl)`
/// (`cxx:43-173`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_make_rotule(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pl: &GpPln,
    pl1: &GpPln,
    pl2: &GpPln,
    opl: Orientation,
    opl1: Orientation,
    opl2: Orientation,
    r: f64,
    ofpl: Orientation,
) -> Result<bool, String> {
    // `gp_Ax3 pos = pl.Position(); gp_Dir dpl = pos.XDirection().Crossed(
    // pos.YDirection()); gp_Dir dfpl = dpl; gp_Dir dplnat = dpl;` (`cxx:56-60`).
    let pos = pl.position().clone();
    let mut dpl = pos
        .x_direction()
        .crossed(pos.y_direction())
        .map_err(|e| e.to_string())?;
    let mut dfpl = dpl;
    let dplnat = dpl;

    // `if (opl == TopAbs_REVERSED) { dpl.Reverse(); } if (ofpl ==
    // TopAbs_REVERSED) { dfpl.Reverse(); }` (`cxx:61-68`).
    if opl == Orientation::Reversed {
        dpl.reverse();
    }
    if ofpl == Orientation::Reversed {
        dfpl.reverse();
    }

    // `pos = pl1.Position(); gp_Dir dpl1 = pos.XDirection().Crossed(
    // pos.YDirection()); if (opl1 == TopAbs_REVERSED) { dpl1.Reverse(); }`
    // (`cxx:69-74`).
    let pos1 = pl1.position().clone();
    let mut dpl1 = pos1
        .x_direction()
        .crossed(pos1.y_direction())
        .map_err(|e| e.to_string())?;
    if opl1 == Orientation::Reversed {
        dpl1.reverse();
    }

    // `pos = pl2.Position(); gp_Dir dpl2 = pos.XDirection().Crossed(
    // pos.YDirection()); if (opl2 == TopAbs_REVERSED) { dpl2.Reverse(); }`
    // (`cxx:75-80`).
    let pos2 = pl2.position().clone();
    let mut dpl2 = pos2
        .x_direction()
        .crossed(pos2.y_direction())
        .map_err(|e| e.to_string())?;
    if opl2 == Orientation::Reversed {
        dpl2.reverse();
    }

    // `double alpha = dpl1.Angle(dpl2);` (`cxx:82`).
    let alpha = dpl1.angle_tol(&dpl2);

    // `IntAna_QuadQuadGeo LInt(pl1, pl2, Precision::Angular(),
    // Precision::Confusion()); gp_Pnt ptor, pcirc; if (LInt.IsDone()) {
    // pcirc = ElCLib::Value(ElCLib::Parameter(LInt.Line(1), pl.Location()),
    // LInt.Line(1)); ptor.SetCoord(pcirc.X() + r * dpl.X(), ...); } else {
    // return false; }` (`cxx:84-96`). `plane_plane_line` is the same
    // constructor's `IsDone()` / `Line(1)` pair (`IntAna_QuadQuadGeo.cxx:389-511`).
    let Some(lin) = plane_plane_line(pl1, pl2, ANGULAR, CONFUSION) else {
        return Ok(false);
    };
    let pcirc = clib::line_value(&lin, clib::line_parameter(&lin, &pl.location()));
    let ptor = GpPnt::new(
        pcirc.x() + r * dpl.x(),
        pcirc.y() + r * dpl.y(),
        pcirc.z() + r * dpl.z(),
    );

    // `gp_Ax3 ppos(ptor, dpl.Reversed(), dpl1); if (ppos.YDirection().Dot(dpl2)
    // < 0.) { ppos.YReverse(); }` (`cxx:98-101`).
    let mut ppos = GpAx3::from_point_dir_x(ptor, dpl.reversed(), &dpl1).map_err(|e| e.to_string())?;
    if ppos.y_direction().dot(&dpl2) < 0.0 {
        ppos.y_reverse();
    }

    // `occ::handle<Geom_ToroidalSurface> gtor = new Geom_ToroidalSurface(ppos,
    // r, r); Data->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gtor, DStr));`
    // (`cxx:102-103`).
    let gtor = GpTorus::new(ppos.clone(), r, r).map_err(|e| e.to_string())?;
    let gtor_handle: Arc<dyn GeomSurface> = Arc::new(GeomTorus::new(gtor.clone()));
    data.set_surf(ch_fi_kpart_index_surface_in_ds(gtor_handle, ds));

    // `gp_Pnt pp; gp_Vec du, dv; ElSLib::TorusD1(0., M_PI / 2, ppos, r, r, pp,
    // du, dv); gp_Dir drot(du.Crossed(dv)); bool reversecur = (drot.Dot(dplnat)
    // <= 0.); bool reversefil = (drot.Dot(dfpl) <= 0.);` (`cxx:108-112`).
    let (pp, du, dv) = surface_eval::torus_d1(&gtor, 0.0, std::f64::consts::FRAC_PI_2);
    let drot = GpDir::from_vec(&du.crossed(&dv)).map_err(|e| e.to_string())?;
    let reversecur = drot.dot(&dplnat) <= 0.0;
    let reversefil = drot.dot(&dfpl) <= 0.0;

    // `if (reversefil) { Data->ChangeOrientation() = TopAbs_REVERSED; } else {
    // ... TopAbs_FORWARD; }` (`cxx:113-121`).
    data.set_orientation(if reversefil {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // `gp_Ax2 circAx2 = ppos.Ax2(); circAx2.SetLocation(pcirc);
    // occ::handle<Geom_Circle> GC = new Geom_Circle(circAx2, r);` (`cxx:127-129`).
    let mut circ_ax2: GpAx2 = ppos.ax2();
    circ_ax2.set_location(pcirc);
    let circle3d: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(GpCirc::new(circ_ax2, r)));

    // `double u, v; ElSLib::Parameters(pl, pcirc, u, v); gp_Pnt2d p2dcirc(u,
    // v);` (`cxx:130-132`).
    let (u, v) = slib::plane_parameters(&pl.position(), &pcirc);
    let p2dcirc = GpPnt2d::new(u, v);

    // `gp_Dir2d dx2d(dpl1.Dot(pl.Position().XDirection()),
    // dpl1.Dot(pl.Position().YDirection())); gp_Dir2d dy2d(
    // ppos.YDirection().Dot(pl.Position().XDirection()),
    // ppos.YDirection().Dot(pl.Position().YDirection()));` (`cxx:133-136`).
    let plpos = pl.position();
    let dx2d = GpDir2d::new(
        dpl1.dot(plpos.x_direction()),
        dpl1.dot(plpos.y_direction()),
    )
    .map_err(|e| e.to_string())?;
    let dy2d = GpDir2d::new(
        ppos.y_direction().dot(plpos.x_direction()),
        ppos.y_direction().dot(plpos.y_direction()),
    )
    .map_err(|e| e.to_string())?;

    // `gp_Ax22d circ2dax(p2dcirc, dx2d, dy2d); occ::handle<Geom2d_Circle> GC2d =
    // new Geom2d_Circle(circ2dax, r);` (`cxx:137-138`).
    let circ2dax = GpAx22d::new(p2dcirc, dx2d, dy2d).map_err(|e| e.to_string())?;
    let gc2d: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(GpCirc2d::new(circ2dax, r)));

    // `gp_Pnt2d p2dlin(0., M_PI / 2); occ::handle<Geom2d_Line> GL2d = new
    // Geom2d_Line(p2dlin, gp::DX2d());` (`cxx:139-140`).
    let gl2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(0.0, std::f64::consts::FRAC_PI_2),
        DX2D,
    ));

    // `TopAbs_Orientation trans = TopAbs_REVERSED; if (reversecur) { trans =
    // TopAbs_FORWARD; } Data->ChangeInterferenceOnS1().SetInterference(
    // ChFiKPart_IndexCurveInDS(GC, DStr), trans, GC2d, GL2d);` (`cxx:140-148`).
    let trans = if reversecur {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let ic = ch_fi_kpart_index_curve_in_ds(Some(circle3d), ds);
    data.interference_on_s1_mut()
        .set_interference(ic, trans, Some(gc2d), Some(gl2d));

    // `p2dlin.SetCoord(0., M_PI); occ::handle<Geom2d_Line> GL2dcoin = new
    // Geom2d_Line(p2dlin, gp::DX2d());
    // Data->ChangeInterferenceOnS2().SetInterference(
    // ChFiKPart_IndexCurveInDS(bid, DStr), trans, bid2d, GL2dcoin);`
    // (`cxx:155-159`) — `bid` / `bid2d` are null handles.
    let gl2dcoin: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(
        GpPnt2d::new(0.0, std::f64::consts::PI),
        DX2D,
    ));
    let ic = ch_fi_kpart_index_curve_in_ds(None, ds);
    data.interference_on_s2_mut()
        .set_interference(ic, trans, None, Some(gl2dcoin));

    // `Data->ChangeVertexFirstOnS1().SetPoint(pp);
    // ElSLib::TorusD0(alpha, M_PI / 2, ppos, r, r, pp);
    // Data->ChangeVertexLastOnS1().SetPoint(pp);` (`cxx:163-165`).
    data.vertex_first_on_s1_mut().set_point(&pp);
    let last = slib::torus_value(&gtor, alpha, std::f64::consts::FRAC_PI_2);
    data.vertex_last_on_s1_mut().set_point(&last);

    // `Data->ChangeInterferenceOnS1().SetFirstParameter(0.);
    // SetLastParameter(alpha); ...OnS2().SetFirstParameter(0.);
    // SetLastParameter(alpha);` (`cxx:166-169`).
    data.interference_on_s1_mut().set_first_parameter(0.0);
    data.interference_on_s1_mut().set_last_parameter(alpha);
    data.interference_on_s2_mut().set_first_parameter(0.0);
    data.interference_on_s2_mut().set_last_parameter(alpha);

    Ok(true)
}
