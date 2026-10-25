//! The three `ChFiKPart_ComputeData::ComputeCorner` overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData.cxx:641-762`
//! (OCCT 8.0.0). Rust has no overloading, so each C++ overload keeps its own
//! name, in declaration order:
//!
//! ```text
//! ch_fi_kpart_compute_corner_plane_cylinder   cxx:641-714
//! ch_fi_kpart_compute_corner_sphere           cxx:716-732  -> ChFiKPart_Sphere
//! ch_fi_kpart_compute_corner_rotule           cxx:734-762  -> ChFiKPart_MakeRotule
//! ```
//!
//! All three are reached from `ChFi3d_Builder::PerformSetOfKPart`
//! (`ChFi3d_Builder_2.cxx:3070-3072`) together with `Compute`; the first is the
//! torus-joint corner between a plane and any second surface whose section
//! spine is a circle.

use std::sync::Arc;

use occt_core::gp::{GpCylinder, GpPnt2d};

use crate::abs::Orientation;
use crate::chfi_ds::SurfData;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::corner_spine::ch_fi_kpart_corner_spine;
use super::fcts::ch_fi_kpart_p_curve;
use super::fillet_plane_cylinder_circle::ch_fi_kpart_make_fillet_plane_cylinder_circle;
use super::rotule::ch_fi_kpart_make_rotule;
use super::sphere_corner::ch_fi_kpart_sphere;

/// `ComputeCorner(DStr, Data, S1, S2, OrFace1, /*unnamed*/, Or1, Or2, minRad,
/// majRad, P1S1, P2S1, P1S2, P2S2)` (`cxx:641-714`).
///
/// The second orientation (`cxx:647`) is unnamed in OCCT and never read.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_corner_plane_cylinder(
    ds: &mut DataStructure,
    data: &mut SurfData,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    or_face1: Orientation,
    _or_face2: Orientation,
    or1: Orientation,
    or2: Orientation,
    min_rad: f64,
    maj_rad: f64,
    p1s1: GpPnt2d,
    p2s1: GpPnt2d,
    p1s2: GpPnt2d,
    p2s2: GpPnt2d,
) -> Result<bool, String> {
    // `GeomAbs_SurfaceType typ1 = S1->GetType(); GeomAbs_SurfaceType typ2 =
    // S2->GetType(); if (typ1 != GeomAbs_Plane) { throw
    // Standard_ConstructionError("la face du conge torique doit etre plane"); }`
    // (`cxx:656-662`).
    let typ1 = s1.get_type();
    let typ2 = s2.get_type();
    if typ1 != GeomAbsSurfaceType::Plane {
        return Err(
            "Standard_ConstructionError: la face du conge torique doit etre plane".into(),
        );
    }

    // `gp_Cylinder cyl; gp_Circ circ; double First, Last, fu, lu;
    // ChFiKPart_CornerSpine(S1, S2, P1S1, P2S1, P1S2, P2S2, majRad, cyl, circ,
    // First, Last);` (`cxx:665-669`).
    let (cyl0, circ, first, last) =
        ch_fi_kpart_corner_spine(s1, s2, p1s1, p2s1, p1s2, p2s2, maj_rad)?;

    // `if (typ2 == GeomAbs_Cylinder) { cyl = S2->Cylinder(); fu = P1S2.X(); lu =
    // P2S2.X(); } else { fu = First; lu = Last; }` (`cxx:670-682`).
    let (cyl, fu, lu): (GpCylinder, f64, f64) = if typ2 == GeomAbsSurfaceType::Cylinder {
        let cyl = s2
            .cylinder()
            .ok_or_else(|| "ComputeCorner: S2->Cylinder".to_string())?;
        (cyl, p1s2.x(), p2s2.x())
    } else {
        (cyl0, first, last)
    };

    // `surfok = ChFiKPart_MakeFillet(DStr, Data, S1->Plane(), cyl, fu, lu, Or1,
    // Or2, minRad, circ, First, OrFace1, true);` (`cxx:683-693`).
    let pl = s1
        .plane()
        .ok_or_else(|| "ComputeCorner: S1->Plane".to_string())?;
    let surfok = ch_fi_kpart_make_fillet_plane_cylinder_circle(
        ds, data, &pl, &cyl, fu, lu, or1, or2, min_rad, &circ, first, or_face1, true,
    )?;
    if !surfok {
        return Ok(false);
    }

    // `if (typ2 != GeomAbs_Cylinder) {
    // Data->ChangeInterferenceOnS2().ChangePCurveOnFace() =
    // ChFiKPart_PCurve(P1S2, P2S2, First, Last); }` (`cxx:697-700`).
    if typ2 != GeomAbsSurfaceType::Cylinder {
        let pc = ch_fi_kpart_p_curve(p1s2, p2s2, first, last)?;
        data.interference_on_s2_mut()
            .set_p_curve_on_face(Some(Arc::new(pc)));
    }

    // `Data->ChangeVertexFirstOnS1().SetPoint(S1->Value(P1S1.X(), P1S1.Y()));
    // ...LastOnS1(P2S1); ...FirstOnS2(P1S2); ...LastOnS2(P2S2);
    // ChangeInterferenceOnS1().SetFirstParameter(First);
    // SetLastParameter(Last); ChangeInterferenceOnS2().SetFirstParameter(First);
    // SetLastParameter(Last); return true;` (`cxx:701-709`).
    let v = s1
        .d0(p1s1.x(), p1s1.y())
        .ok_or_else(|| "ComputeCorner: S1->Value (P1S1)".to_string())?;
    data.vertex_first_on_s1_mut().set_point(&v);
    let v = s1
        .d0(p2s1.x(), p2s1.y())
        .ok_or_else(|| "ComputeCorner: S1->Value (P2S1)".to_string())?;
    data.vertex_last_on_s1_mut().set_point(&v);
    let v = s2
        .d0(p1s2.x(), p1s2.y())
        .ok_or_else(|| "ComputeCorner: S2->Value (P1S2)".to_string())?;
    data.vertex_first_on_s2_mut().set_point(&v);
    let v = s2
        .d0(p2s2.x(), p2s2.y())
        .ok_or_else(|| "ComputeCorner: S2->Value (P2S2)".to_string())?;
    data.vertex_last_on_s2_mut().set_point(&v);
    data.interference_on_s1_mut().set_first_parameter(first);
    data.interference_on_s1_mut().set_last_parameter(last);
    data.interference_on_s2_mut().set_first_parameter(first);
    data.interference_on_s2_mut().set_last_parameter(last);
    Ok(true)
}

/// `ComputeCorner(DStr, Data, S1, S2, OrFace1, OrFace2, Or1, Or2, Rad, PS1,
/// P1S2, P2S2)` (`cxx:716-732`): the spherical corner.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_corner_sphere(
    ds: &mut DataStructure,
    data: &mut SurfData,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    or_face1: Orientation,
    or_face2: Orientation,
    or1: Orientation,
    or2: Orientation,
    rad: f64,
    ps1: GpPnt2d,
    p1s2: GpPnt2d,
    p2s2: GpPnt2d,
) -> Result<bool, String> {
    // `return ChFiKPart_Sphere(DStr, Data, S1, S2, OrFace1, OrFace2, Or1, Or2,
    // Rad, PS1, P1S2, P2S2);` (`cxx:729`).
    ch_fi_kpart_sphere(
        ds, data, s1, s2, or_face1, or_face2, or1, or2, rad, ps1, p1s2, p2s2,
    )
}

/// `ComputeCorner(DStr, Data, S, S1, S2, OfS, OS, OS1, OS2, Radius)`
/// (`cxx:734-762`): the "rotule" corner between three planes.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_corner_rotule(
    ds: &mut DataStructure,
    data: &mut SurfData,
    s: &SurfaceAdaptor,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    of_s: Orientation,
    os: Orientation,
    os1: Orientation,
    os2: Orientation,
    radius: f64,
) -> Result<bool, String> {
    // `GeomAbs_SurfaceType typ = S->GetType(); ... if (typ != GeomAbs_Plane ||
    // typ1 != GeomAbs_Plane || typ2 != GeomAbs_Plane) { throw
    // Standard_ConstructionError("torus joint only between the planes"); }`
    // (`cxx:745-751`).
    if s.get_type() != GeomAbsSurfaceType::Plane
        || s1.get_type() != GeomAbsSurfaceType::Plane
        || s2.get_type() != GeomAbsSurfaceType::Plane
    {
        return Err("Standard_ConstructionError: torus joint only between the planes".into());
    }
    let pl = s
        .plane()
        .ok_or_else(|| "ComputeCorner: S->Plane".to_string())?;
    let pl1 = s1
        .plane()
        .ok_or_else(|| "ComputeCorner: S1->Plane".to_string())?;
    let pl2 = s2
        .plane()
        .ok_or_else(|| "ComputeCorner: S2->Plane".to_string())?;

    // `return ChFiKPart_MakeRotule(DStr, Data, S->Plane(), S1->Plane(),
    // S2->Plane(), OS, OS1, OS2, Radius, OfS);` (`cxx:752-761`).
    ch_fi_kpart_make_rotule(ds, data, &pl, &pl1, &pl2, os, os1, os2, radius, of_s)
}
