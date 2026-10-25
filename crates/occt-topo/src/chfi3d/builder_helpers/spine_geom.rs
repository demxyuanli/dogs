//! Corner guidelines built from two tangent frames: `ChFi3d_CircularSpine` and
//! `ChFi3d_Spine`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:846-903`.

use occt_core::elib::clib::{line_parameter, line_value};
use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpPln, GpPnt, GpVec};
use occt_core::precision;
use occt_geom::intana::plane_plane_line;
use occt_geom::{GeomBezierCurve, GeomCircle};

/// `ChFi3d_CircularSpine(WFirst, WLast, Pdeb, Vdeb, Pfin, Vfin, rad)`
/// (`cxx:846-880`): the arc of radius `rad` tangent to `Vdeb` at `Pdeb` and to
/// `Vfin` at `Pfin`.
///
/// The centre is the intersection of the two planes `(Pdeb, Vdeb)` and
/// `(Pfin, Vfin)`; `WFirst` / `WLast` come back as the angles of the two radius
/// vectors measured from the first one.
///
/// `None` is OCCT's null handle: `IntAna_QuadQuadGeo` not done (`cxx:878-880`)
/// or the two tangent frames crossing (`cxx:867-870`). `Err` carries the
/// `gp_Dir` / `gp_Ax2` construction errors OCCT lets escape — a null tangent or
/// a null cross product.
pub fn circular_spine(
    w_first: &mut f64,
    w_last: &mut f64,
    p_deb: &GpPnt,
    v_deb: &GpVec,
    p_fin: &GpPnt,
    v_fin: &GpVec,
    rad: f64,
) -> Result<Option<GeomCircle>, String> {
    let dir_deb = GpDir::from_vec(v_deb).map_err(|e| e.to_string())?;
    let dir_fin = GpDir::from_vec(v_fin).map_err(|e| e.to_string())?;
    // `gp_Pln Pl1(Pdeb, gp_Dir(Vdeb)), Pl2(Pfin, gp_Dir(Vfin))`.
    let pl1 = GpPln::new(GpAx3::from_point_dir(*p_deb, dir_deb).map_err(|e| e.to_string())?);
    let pl2 = GpPln::new(GpAx3::from_point_dir(*p_fin, dir_fin).map_err(|e| e.to_string())?);

    // `IntAna_QuadQuadGeo LInt(Pl1, Pl2, Precision::Angular(),
    // Precision::Confusion())`.
    let Some(li) = plane_plane_line(&pl1, &pl2, precision::ANGULAR, precision::CONFUSION) else {
        return Ok(None);
    };

    let cen_deb = line_value(&li, line_parameter(&li, p_deb));
    let cen_fin = line_value(&li, line_parameter(&li, p_fin));
    let vv_deb = GpVec::from_pnts(&cen_deb, p_deb);
    let vv_fin = GpVec::from_pnts(&cen_fin, p_fin);
    let dd_deb = GpDir::from_vec(&vv_deb).map_err(|e| e.to_string())?;
    let dd_fin = GpDir::from_vec(&vv_fin).map_err(|e| e.to_string())?;

    // `if (Vdeb.Crossed(vvdeb).Dot(Vfin.Crossed(vvfin)) > 0.) return null;`
    if v_deb.crossed(&vv_deb).dot(&v_fin.crossed(&vv_fin)) > 0.0 {
        return Ok(None);
    }

    // `gp_Ax2 circax2(cendeb, dddeb ^ ddfin, dddeb);`
    let normal = dd_deb.crossed(&dd_fin).map_err(|e| e.to_string())?;
    let ax2 = GpAx2::new(cen_deb, normal, dd_deb).map_err(|e| e.to_string())?;
    let circ = GpCirc::new(ax2, rad);

    *w_first = 0.0;
    *w_last = dd_deb.angle(&dd_fin);
    Ok(Some(GeomCircle::new(circ)))
}

/// `ChFi3d_Spine(pd, vd, pf, vf, R)` (`cxx:888-903`): the four-pole Bezier
/// guideline of a corner, poles at the two tangent points and at
/// `P + V * (0.5 * tan((pi - angle) / 2) * R)`.
///
/// `vd` and `vf` are modified in place in OCCT (`gp_Vec&`); the port takes them
/// by value so the caller's vectors are untouched.
pub fn spine(pd: &GpPnt, vd: &GpVec, pf: &GpPnt, vf: &GpVec, r: f64) -> GeomBezierCurve {
    let fac = 0.5 * ((std::f64::consts::PI - vd.angle(vf)) * 0.5).tan();
    let p1 = *pd;
    let vd = vd.multiply_scalar(fac * r);
    let p2 = GpPnt::new(pd.x() + vd.x(), pd.y() + vd.y(), pd.z() + vd.z());
    let p4 = *pf;
    let vf = vf.multiply_scalar(fac * r);
    let p3 = GpPnt::new(pf.x() + vf.x(), pf.y() + vf.y(), pf.z() + vf.z());
    // `pol(1) = pd; pol(2) = ...; pol(3) = ...; pol(4) = pf;` (the Bezier arm is
    // `pol(3)` = the `vf` one and `pol(4)` = `pf`, `cxx:899-901`).
    GeomBezierCurve::new(vec![p1, p2, p3, p4]).expect("four distinct poles")
}
