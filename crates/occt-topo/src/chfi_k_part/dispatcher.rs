//! `ChFiKPart_ComputeData::Compute` — the dispatch from a fillet / chamfer
//! guideline and its two faces to the analytic corner surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiKPart/ChFiKPart_ComputeData.cxx:52-637`
//! (OCCT 8.0.0).
//!
//! `Compute` picks the maker by the two `GeomAbs_SurfaceType`s and, for a
//! fillet, by the type of the elementary spine (`Spine->CurrentElementarySpine
//! (Iedge).GetType()`); for a chamfer it picks by `Spine->IsChamfer()`
//! (`ChFiDS_Sym` / `ChFiDS_TwoDist` / the distance-angle kind) and then by the
//! same surface pair (`cxx:216-637`).
//!
//! UNPORTED arms. Each of them answers `Ok(false)`, which is the value the
//! caller's `if (!Compute(...)) { }` empty arm handles
//! (`ChFi3d_Builder_2.cxx:3070-3072`), and names the `.cxx` it comes from:
//! - `ChFiKPart_ComputeData_ChPlnPln.cxx:53-269`,
//!   `_ChPlnCyl.cxx:64-368` / `:372-592`, `_ChPlnCon.cxx:49-406`,
//!   `_ChAsymPlnPln.cxx:53-263`, `_ChAsymPlnCyl.cxx:65-373` / `:376-653`,
//!   `_ChAsymPlnCon.cxx:57-705` — the chamfer arms (`cxx:216-636`);
//! - `ChFiKPart_ComputeData_Sphere.cxx:47-204` / `_Rotule.cxx:43-169` — the
//!   three `ComputeCorner` overloads (`cxx:641-714`, `:716-732`, `:734-762`),
//!   which come with the corner batch together with `ChFiKPart_ProjPC`.
//!
//! Any other surface pair is `Standard_NotImplemented("particular case not
//! written")` in OCCT (`cxx:210-212`, `:354-356`, `:492-494`, `:631-633`), so
//! the port answers `Err` there rather than `false`.

use crate::abs::Orientation;
use crate::chfi_ds::{ChamfMethod, ChamfMode, Spine, SpineClass, SpineHandle, SurfData};
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::chamfer_asym_plane_plane::ch_fi_kpart_make_ch_asym_plane_plane;
use super::chamfer_plane_plane::ch_fi_kpart_make_chamfer_plane_plane;
use super::fillet_plane_cone::ch_fi_kpart_make_fillet_plane_cone;
use super::fillet_plane_cylinder::ch_fi_kpart_make_fillet_plane_cylinder_line;
use super::fillet_plane_cylinder_circle::ch_fi_kpart_make_fillet_plane_cylinder_circle;
use super::fillet_plane_plane::ch_fi_kpart_make_fillet_plane_plane;

/// `ChFiKPart_ComputeData::Compute(DStr, Data, S1, S2, Or1, Or2, Sp, Iedge)`
/// (`cxx:50-638`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi_kpart_compute_data(
    ds: &mut DataStructure,
    data: &mut SurfData,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    or1: Orientation,
    or2: Orientation,
    spine: &SpineHandle,
    iedge: i32,
) -> Result<bool, String> {
    // `double Wref = 0.;` (`cxx:60`).
    let wref = 0.0;

    let mut sp = spine.write().unwrap_or_else(|e| e.into_inner());

    // `if (!Spine.IsNull()) { ctyp = Spine->CurrentElementarySpine(Iedge)
    // .GetType(); } else { ctyp = CSpine->CurrentElementarySpine(Iedge)
    // .GetType(); }` (`cxx:69-76`). The down-casts to `ChFiDS_FilSpine` /
    // `ChFiDS_ChamfSpine` are the two `SpineClass` kinds. `CurrentElementarySpine`
    // also caches the edge it is asked about (`ChFiDS_Spine.cxx:866-874`), so
    // the call happens whichever arm follows; the cylinder arms read the type
    // back (`cxx:109`, `:144`, `:160`).
    let is_fil = sp.class() == SpineClass::Fil;
    let is_chamf = sp.class() == SpineClass::Chamf;
    let ctyp = sp
        .current_elementary_spine(iedge)
        .ok_or_else(|| "ChFiKPart_ComputeData::Compute: CurrentElementarySpine".to_string())?
        .get_type();

    let typ1 = s1.get_type();
    let typ2 = s2.get_type();

    // `TopAbs_Orientation OrFace1 = TopAbs_FORWARD, OrFace2 = TopAbs_FORWARD;
    // occ::handle<BRepAdaptor_Surface> HS =
    // occ::down_cast<BRepAdaptor_Surface>(S1); if (!HS.IsNull()) { OrFace1 =
    // HS->Face().Orientation(); } ... S2 ...` (`cxx:79-89`).
    let or_face1 = match s1.face() {
        Some(f) => f.0.orientation(),
        None => Orientation::Forward,
    };
    // `OrFace2` (`cxx:79-89`): the Cylinder/Plane and Cone/Plane arms pass it
    // as `Ofpl` (`cxx:157`, `:205`).
    let or_face2 = match s2.face() {
        Some(f) => f.0.orientation(),
        None => Orientation::Forward,
    };

    if is_fil {
        // `double Radius = Spine->Radius(Iedge);` (`cxx:94`). OCCT raises
        // `Standard_DomainError` from `Radius(int)` when the radius is not
        // constant on the edge (`ChFiDS_FilSpine.cxx:340-343`).
        let radius = sp.fil_radius_of_edge_index(iedge).ok_or_else(|| {
            "ChFiKPart_ComputeData::Compute: Spine->Radius(Iedge) is not constant".to_string()
        })?;

        match (typ1, typ2) {
            // `if (typ1 == GeomAbs_Plane && typ2 == GeomAbs_Plane)`
            // (`cxx:94-106`).
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane) => {
                let pl1 = s1
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S1->Plane()".to_string())?;
                let pl2 = s2
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S2->Plane()".to_string())?;
                let line = sp
                    .line()
                    .ok_or_else(|| "ChFiKPart_ComputeData: Spine->Line()".to_string())?;
                ch_fi_kpart_make_fillet_plane_plane(
                    ds, data, &pl1, &pl2, or1, or2, radius, &line, wref, or_face1,
                )
            }
            // `else if (typ1 == GeomAbs_Plane && typ2 == GeomAbs_Cylinder)`
            // (`cxx:107-141`): `MakeFillet(Plane, Cylinder, S2's U bounds,
            // Or1, Or2, Radius, Spine->Line() / Circle(), Wref, OrFace1,
            // true)`.
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder) => {
                let pl = s1
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S1->Plane()".to_string())?;
                let cy = s2
                    .cylinder()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S2->Cylinder()".to_string())?;
                make_fillet_plane_cylinder(
                    ds,
                    data,
                    &pl,
                    &cy,
                    s2.first_u_parameter(),
                    s2.last_u_parameter(),
                    or1,
                    or2,
                    radius,
                    ctyp,
                    &mut sp,
                    wref,
                    or_face1,
                    true,
                )
            }
            // `else if (typ1 == GeomAbs_Cylinder && typ2 == GeomAbs_Plane)`
            // (`cxx:142-175`): the plane is passed first, `Or2`/`Or1` swapped
            // and `OrFace2` used as `Ofpl`, `plandab = false`.
            (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane) => {
                let pl = s2
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S2->Plane()".to_string())?;
                let cy = s1
                    .cylinder()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S1->Cylinder()".to_string())?;
                make_fillet_plane_cylinder(
                    ds,
                    data,
                    &pl,
                    &cy,
                    s1.first_u_parameter(),
                    s1.last_u_parameter(),
                    or2,
                    or1,
                    radius,
                    ctyp,
                    &mut sp,
                    wref,
                    or_face2,
                    false,
                )
            }
            // `else if (typ1 == GeomAbs_Plane && typ2 == GeomAbs_Cone)`
            // (`cxx:177-192`), a single arm — the cone makers take the
            // guideline as `Spine->Circle()` without consulting `ctyp`.
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone) => {
                let pl = s1
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S1->Plane()".to_string())?;
                let co = s2
                    .cone()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S2->Cone()".to_string())?;
                let circle = sp
                    .circle()
                    .ok_or_else(|| "ChFiKPart_ComputeData: Spine->Circle()".to_string())?;
                ch_fi_kpart_make_fillet_plane_cone(
                    ds,
                    data,
                    &pl,
                    &co,
                    s2.first_u_parameter(),
                    s2.last_u_parameter(),
                    or1,
                    or2,
                    radius,
                    &circle,
                    wref,
                    or_face1,
                    true,
                )
            }
            // `else if (typ1 == GeomAbs_Cone && typ2 == GeomAbs_Plane)`
            // (`cxx:193-208`).
            (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => {
                let pl = s2
                    .plane()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S2->Plane()".to_string())?;
                let co = s1
                    .cone()
                    .ok_or_else(|| "ChFiKPart_ComputeData: S1->Cone()".to_string())?;
                let circle = sp
                    .circle()
                    .ok_or_else(|| "ChFiKPart_ComputeData: Spine->Circle()".to_string())?;
                ch_fi_kpart_make_fillet_plane_cone(
                    ds,
                    data,
                    &pl,
                    &co,
                    s1.first_u_parameter(),
                    s1.last_u_parameter(),
                    or2,
                    or1,
                    radius,
                    &circle,
                    wref,
                    or_face2,
                    false,
                )
            }
            // `else { throw Standard_NotImplemented("particular case not
            // written"); }` (`cxx:210-212`).
            _ => Err("ChFiKPart_ComputeData::Compute: particular case not written".into()),
        }
    } else if is_chamf {
        // `ChFiDS_ChamfMode aMode = CSpine->Mode();` (`cxx:216`), then the
        // three `CSpine->IsChamfer()` kinds (`cxx:219`, `:358`, `:495`).
        let a_mode = sp.mode();
        // `CSpine->Line()` / `CSpine->Circle()` are read by the arm the pair of
        // types selects, so the guiders are fetched inside the arms.
        match sp.is_chamfer() {
            // `if (CSpine->IsChamfer() == ChFiDS_Sym)` (`cxx:219-356`).
            ChamfMethod::Sym => {
                // `double dis; CSpine->GetDist(dis);` (`cxx:221-222`). The
                // accessor raises `Standard_DomainError` for another kind
                // (`ChFiDS_ChamfSpine.cxx:47-52`).
                let dis = sp.get_dist().ok_or_else(|| {
                    "ChFiKPart_ComputeData::Compute: CSpine->GetDist".to_string()
                })?;
                // `ChFiKPart_MakeChamfer(..., dis, dis, CSpine->Line(), ...)`
                // (`cxx:223-355`).
                chamfer_symmetric_arm(
                    ds, data, a_mode, s1, s2, or1, or2, dis, dis, &mut sp, wref, or_face1,
                )
            }
            // `else if (CSpine->IsChamfer() == ChFiDS_TwoDist)`
            // (`cxx:358-493`).
            ChamfMethod::TwoDist => {
                // `double dis1, dis2; CSpine->Dists(dis1, dis2);`
                // (`cxx:360-361`).
                let (dis1, dis2) = sp.dists().ok_or_else(|| {
                    "ChFiKPart_ComputeData::Compute: CSpine->Dists".to_string()
                })?;
                // Note the `Cylinder/Plane` arms of this kind pass
                // `dis2, dis1` in that order (`cxx:429-430`), unlike the `Sym`
                // kind, which passes `dis, dis`; the swap has to survive.
                chamfer_symmetric_arm(
                    ds, data, a_mode, s1, s2, or1, or2, dis1, dis2, &mut sp, wref, or_face1,
                )
            }
            // `else { double dis, Angle; bool DisOnP = true;
            // CSpine->GetDistAngle(dis, Angle); ... }` (`cxx:496-512`).
            ChamfMethod::DistAngle => {
                let (dis, angle) = sp.get_dist_angle().ok_or_else(|| {
                    "ChFiKPart_ComputeData::Compute: CSpine->GetDistAngle".to_string()
                })?;
                match (typ1, typ2) {
                    // `ChFiKPart_MakeChAsym(DStr, Data, S1->Plane(),
                    // S2->Plane(), Or1, Or2, dis, Angle, CSpine->Line(),
                    // Wref, OrFace1, DisOnP)` (`cxx:502-513`).
                    (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane) => {
                        let pl1 = s1
                            .plane()
                            .ok_or_else(|| "ChFiKPart_ComputeData: S1->Plane()".to_string())?;
                        let pl2 = s2
                            .plane()
                            .ok_or_else(|| "ChFiKPart_ComputeData: S2->Plane()".to_string())?;
                        let line = sp
                            .line()
                            .ok_or_else(|| "ChFiKPart_ComputeData: CSpine->Line()".to_string())?;
                        ch_fi_kpart_make_ch_asym_plane_plane(
                            ds, data, &pl1, &pl2, or1, or2, dis, angle, &line, wref, or_face1,
                            true,
                        )
                    }
                    // UNPORTED `ChFiKPart_ComputeData_ChAsymPlnCyl.cxx:65-373`
                    // (`Plane/Cylinder`, `cxx:514-553`), `:376-653`
                    // (`Cylinder/Plane` with `dis` and `Angle` and the
                    // `Or2/Or1` swap, `cxx:554-591`), and
                    // `ChFiKPart_ComputeData_ChAsymPlnCon.cxx:57-705`
                    // (`Plane/Cone`, `cxx:592-611`; `Cone/Plane`,
                    // `cxx:612-628`). `surfok` stays `false` (`cxx:64`).
                    (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
                    | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane)
                    | (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
                    | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => Ok(false),
                    // `else { throw Standard_NotImplemented("particular case
                    // not written"); }` (`cxx:630-633`).
                    _ => Err(
                        "ChFiKPart_ComputeData::Compute: particular case not written".into(),
                    ),
                }
            }
        }
    } else {
        // `return surfok;` with `surfok = false` (`cxx:64`, `:636`).
        Ok(false)
    }
}

/// The `ChFiDS_Sym` arm (`cxx:223-355`) and the `ChFiDS_TwoDist` arm
/// (`cxx:362-492`). The two differ only in the distances they pass — `dis, dis`
/// against `dis1, dis2` — and in that the `Cylinder/Plane` arm of `TwoDist`
/// passes `dis2, dis1` in that order (`cxx:429-430`, against `cxx:308-309`
/// where both are `dis`). Folding them here keeps the `dis1, dis2` swap where
/// the `.cxx` writes it.
#[allow(clippy::too_many_arguments)]
fn chamfer_symmetric_arm(
    ds: &mut DataStructure,
    data: &mut SurfData,
    the_mode: ChamfMode,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    or1: Orientation,
    or2: Orientation,
    dis1: f64,
    dis2: f64,
    sp: &mut Spine,
    wref: f64,
    or_face1: Orientation,
) -> Result<bool, String> {
    let typ1 = s1.get_type();
    let typ2 = s2.get_type();
    match (typ1, typ2) {
        // `ChFiKPart_MakeChamfer(DStr, Data, aMode, S1->Plane(), S2->Plane(),
        // Or1, Or2, dis1, dis2, CSpine->Line(), Wref, OrFace1)`
        // (`cxx:225-237`, `:363-375`).
        (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane) => {
            let pl1 = s1
                .plane()
                .ok_or_else(|| "ChFiKPart_ComputeData: S1->Plane()".to_string())?;
            let pl2 = s2
                .plane()
                .ok_or_else(|| "ChFiKPart_ComputeData: S2->Plane()".to_string())?;
            let line = sp
                .line()
                .ok_or_else(|| "ChFiKPart_ComputeData: CSpine->Line()".to_string())?;
            ch_fi_kpart_make_chamfer_plane_plane(
                ds, data, the_mode, &pl1, &pl2, or1, or2, dis1, dis2, &line, wref, or_face1,
            )
        }
        // UNPORTED `ChFiKPart_ComputeData_ChPlnCyl.cxx:64-368` — the
        // `Plane/Cylinder` arm, guideline a line (`cxx:238-277` for `Sym`), and
        // `:372-592` for guideline a circle (`cxx:241-258`, `:380-395`). The
        // `ctyp` test picks between them (`cxx:241`, `:380`).
        //
        // UNPORTED `ChFiKPart_ComputeData_ChPlnCyl.cxx:372-592` — the
        // `Cylinder/Plane` arm, `S2->Plane()` / `S1->Cylinder()` with
        // `Or2, Or1`, `S1`'s U bounds, `OrFace2` as `Of1` and `plamdab = false`
        // (`cxx:278-315`, `:416-454`).
        //
        // UNPORTED `ChFiKPart_ComputeData_ChPlnCon.cxx:49-406` — the
        // `Plane/Cone` (`cxx:319-335`, `:455-471`) and `Cone/Plane`
        // (`cxx:336-352`, `:472-490`) arms, both guideline a circle.
        // `surfok` stays `false` (`cxx:64`).
        (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
        | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane)
        | (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
        | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => Ok(false),
        // `else { throw Standard_NotImplemented("particular case not
        // written"); }` (`cxx:353-356`, `:490-493`).
        _ => Err("ChFiKPart_ComputeData::Compute: particular case not written".into()),
    }
}

/// The `if (ctyp == GeomAbs_Line) { ... Spine->Line() ... } else { ...
/// Spine->Circle() ... }` pair shared by the Plane/Cylinder (`cxx:109-140`)
/// and Cylinder/Plane (`cxx:144-174`) arms.
#[allow(clippy::too_many_arguments)]
fn make_fillet_plane_cylinder(
    ds: &mut DataStructure,
    data: &mut SurfData,
    pln: &occt_core::gp::GpPln,
    cyl: &occt_core::gp::GpCylinder,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    ctyp: GeomAbsCurveType,
    sp: &mut Spine,
    wref: f64,
    ofpl: Orientation,
    plandab: bool,
) -> Result<bool, String> {
    if ctyp == GeomAbsCurveType::Line {
        let line = sp
            .line()
            .ok_or_else(|| "ChFiKPart_ComputeData: Spine->Line()".to_string())?;
        ch_fi_kpart_make_fillet_plane_cylinder_line(
            ds, data, pln, cyl, fu, lu, or1, or2, radius, &line, wref, ofpl, plandab,
        )
    } else {
        let circle = sp
            .circle()
            .ok_or_else(|| "ChFiKPart_ComputeData: Spine->Circle()".to_string())?;
        ch_fi_kpart_make_fillet_plane_cylinder_circle(
            ds, data, pln, cyl, fu, lu, or1, or2, radius, &circle, wref, ofpl, plandab,
        )
    }
}
