//! `ChFiKPart_ComputeData_CS` — the corner guideline of a torus-shaped joint.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData_CS.cxx:24-73`
//! (OCCT 8.0.0).
//!
//! `ChFiKPart_CornerSpine` builds the cylinder and the circle that the torus
//! section of a corner follows: the plane `S1` gives the two directions the
//! fillet axis is projected on, `S2` gives the point the radius is measured
//! from, and the sign of `V . Vorien` decides on which side the centre sits.
//!
//! The four out-parameters `cyl` / `circ` / `First` / `Last` become the return
//! tuple; `P2S1` is unused in OCCT too (`cxx:27` names it `/*P2S1*/`).

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpDir, GpPnt, GpPnt2d, GpVec};

use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFiKPart_CornerSpine(S1, S2, P1S1, P2S1, P1S2, P2S2, R, cyl, circ, First,
/// Last)` (`cxx:23-73`), returning `(cyl, circ, First, Last)`.
pub fn ch_fi_kpart_corner_spine(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    p1s1: GpPnt2d,
    _p2s1: GpPnt2d,
    p1s2: GpPnt2d,
    p2s2: GpPnt2d,
    r: f64,
) -> Result<(GpCylinder, GpCirc, f64, f64), String> {
    // `gp_Ax3 ax = S1->Plane().Position(); gp_Vec V1(ax.XDirection()); gp_Vec
    // V2(ax.YDirection());` (`cxx:38-39`).
    let ax = s1
        .plane()
        .ok_or_else(|| "ChFiKPart_CornerSpine: S1 is not a plane".to_string())?
        .position();
    let v1 = GpVec::from_xyz(ax.x_direction().xyz());
    let v2 = GpVec::from_xyz(ax.y_direction().xyz());

    // `S2->D1(P1S2.X(), P1S2.Y(), P, du, dv);` (`cxx:41-42`).
    let (p, du, dv) = s2
        .d1(p1s2.x(), p1s2.y())
        .ok_or_else(|| "ChFiKPart_CornerSpine: S2->D1".to_string())?;

    // `gp_Vec V(P, S1->Value(P1S1.X(), P1S1.Y())); V = V.Dot(V1) * V1 +
    // V.Dot(V2) * V2; V.Normalize();` (`cxx:43-45`).
    let to_s1 = s1
        .value(p1s1.x(), p1s1.y())
        .ok_or_else(|| "ChFiKPart_CornerSpine: S1->Value".to_string())?;
    let mut v = GpVec::from_pnts(&p, &to_s1);
    v = v1
        .multiplied_scalar(v.dot(&v1))
        .added(&v2.multiplied_scalar(v.dot(&v2)));
    v.normalize();

    // `gp_Pnt P2 = S2->Value(P2S2.X(), P2S2.Y()); gp_Vec Vorien(P, P2);
    // gp_Pnt cent; gp_Dir dx(V);` (`cxx:46-49`).
    let p2 = s2
        .value(p2s2.x(), p2s2.y())
        .ok_or_else(|| "ChFiKPart_CornerSpine: S2->Value".to_string())?;
    let vorien = GpVec::from_pnts(&p, &p2);
    let mut dx = GpDir::from_vec(&v).map_err(|e| e.to_string())?;
    let cent = if v.dot(&vorien) >= 0.0 {
        // `cxx:50-54`: the centre is on the `V` side and `dx` is reversed.
        let c = GpPnt::new(p.x() + r * v.x(), p.y() + r * v.y(), p.z() + r * v.z());
        dx.reverse();
        c
    } else {
        // `cxx:56-58`.
        GpPnt::new(p.x() - r * v.x(), p.y() - r * v.y(), p.z() - r * v.z())
    };

    // `gp_Dir dy(gp_Vec(cent, P2)); dy = (dx ^ dy) ^ dx;` (`cxx:59-60`).
    let dy = GpDir::from_vec(&GpVec::from_pnts(&cent, &p2)).map_err(|e| e.to_string())?;
    let dy = dx.cross_cross(&dy).map_err(|e| e.to_string())?;

    // `gp_Ax2 circax2(cent, dx ^ dy, dx); gp_Ax3 cylax3(circax2);`
    // (`cxx:61-62`).
    let zdir = dx.crossed(&dy).map_err(|e| e.to_string())?;
    let circax2 = GpAx2::new(cent, zdir, dx).map_err(|e| e.to_string())?;
    let mut cylax3 = GpAx3::from_ax2(&circax2);

    // `if ((du ^ dv).Dot(dx) < 0.) { cylax3.ZReverse(); }` (`cxx:63-66`).
    if du
        .crossed(&dv)
        .dot(&GpVec::from_xyz(dx.xyz()))
        < 0.0
    {
        cylax3.z_reverse();
    }

    // `First = 0.; Last = ElCLib::CircleParameter(circax2, P2);
    // circ.SetPosition(circax2); circ.SetRadius(R); cyl.SetPosition(cylax3);
    // cyl.SetRadius(R);` (`cxx:67-73`).
    let first = 0.0;
    let last = clib::circle_parameter(&circax2, &p2);
    let circ = GpCirc::new(circax2, r);
    let cyl = GpCylinder::new(cylax3, r).map_err(|e| e.to_string())?;
    Ok((cyl, circ, first, last))
}
