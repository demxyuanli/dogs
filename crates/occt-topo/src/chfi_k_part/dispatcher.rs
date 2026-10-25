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
//! - `ChFiKPart_ComputeData_FilPlnCyl.cxx:52-276` (Line spine) and `:283-596`
//!   (Circle spine) — Plane/Cylinder and Cylinder/Plane (`cxx:105-207`);
//! - `ChFiKPart_ComputeData_FilPlnCon.cxx:48-338` — Plane/Cone and Cone/Plane
//!   (`cxx:177-207`);
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
use crate::chfi_ds::{SpineClass, SpineHandle, SurfData};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

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
    // the call happens whichever arm follows; only the UNPORTED
    // cylinder / cone arms read the type back (`cxx:105-175`, `:177-206`).
    let is_fil = sp.class() == SpineClass::Fil;
    let is_chamf = sp.class() == SpineClass::Chamf;
    let _ctyp = sp
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
    // Read by the UNPORTED chamfer and corner makers (`cxx:295`, `:433`,
    // `:571`, `:626`, `:721`).
    let _or_face2 = match s2.face() {
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
            // UNPORTED `ChFiKPart_ComputeData_FilPlnCyl.cxx:52-276` (Line
            // spine) / `:283-596` (Circle spine) (`cxx:107-176`).
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
            | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane) => Ok(false),
            // UNPORTED `ChFiKPart_ComputeData_FilPlnCon.cxx:48-338`
            // (`cxx:177-207`).
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => Ok(false),
            // `else { throw Standard_NotImplemented("particular case not
            // written"); }` (`cxx:210-212`).
            _ => Err("ChFiKPart_ComputeData::Compute: particular case not written".into()),
        }
    } else if is_chamf {
        // UNPORTED `ChFiKPart_ComputeData_Ch*.cxx` — all three
        // `Spine->IsChamfer()` kinds are unported (`cxx:216-636`), so the
        // section kind does not steer the answer yet. The three throw arms
        // (`cxx:354-356`, `:492-494`, `:631-633`) still answer `Err`.
        match (typ1, typ2) {
            (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Plane)
            | (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cylinder)
            | (GeomAbsSurfaceType::Cylinder, GeomAbsSurfaceType::Plane)
            | (GeomAbsSurfaceType::Plane, GeomAbsSurfaceType::Cone)
            | (GeomAbsSurfaceType::Cone, GeomAbsSurfaceType::Plane) => Ok(false),
            _ => Err("ChFiKPart_ComputeData::Compute: particular case not written".into()),
        }
    } else {
        // `return surfok;` with `surfok = false` (`cxx:64`, `:636`).
        Ok(false)
    }
}
