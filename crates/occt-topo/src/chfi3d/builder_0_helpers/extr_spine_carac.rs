//! `ChFi3d_ExtrSpineCarac` — the point, the surface normal and the radius of a
//! fillet section at one parameter of a guideline.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:773-837`
//! (OCCT 8.0.0).
//!
//! The surface of the `SurfData` is inspected: an elementary cylinder answers
//! `Radius()` and the `v` derivative of `ElSLib::D1` (`cxx:794-801`), an
//! elementary torus its `MinorRadius()` and the `u` derivative
//! (`cxx:802-807`). Every other surface — including an approximated one,
//! which is why OCCT accepts the p-curve parameters as the elaboration ones —
//! takes the radius from the `ChFiDS_FilSpine` (constant or law) and the
//! direction from the guideline's `D1` (`cxx:808-832`).
//!
//! `V` comes back normalized and reversed when `sens == 1`
//! (`cxx:834-837`).
//!
//! Divergences that are not heuristics: `p` is the parameter of the p-curve on
//! the surface, while the law and the elementary guideline are evaluated at the
//! same `p` exactly as OCCT does; a null `DStr.Surface(...).Surface()`, a null
//! p-curve, a null `gp_Dir` normalization or a spine that answers no elementary
//! guideline is the `Err` of this port where OCCT dereferences.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::law::{LawFunction, LawHandle};

use crate::chfi_ds::stripe::Stripe;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_ExtrSpineCarac(DStr, cd, i, p, jf, sens, P, V, R)`
/// (`cxx:773-837`).
///
/// Answers `(P, V, R)`; `V` is the reversed surface normal when `sens == 1`.
pub fn ch_fi3d_extr_spine_carac(
    ds: &HDataStructure,
    cd: &Stripe,
    i: i32,
    p: f64,
    jf: i32,
    sens: i32,
) -> Result<(GpPnt, GpVec, f64), String> {
    // `const occ::handle<Geom_Surface>& fffil =
    // DStr.Surface(cd->SetOfSurfData()->Value(i)->Surf()).Surface();`
    // (`cxx:785-786`).
    let sd = i
        .checked_sub(1)
        .and_then(|k| usize::try_from(k).ok())
        .and_then(|k| cd.set_of_surf_data().get(k))
        .ok_or_else(|| format!("ChFi3d_ExtrSpineCarac: SurfData {i} is out of range"))?;
    let surface = ds
        .surface(sd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_ExtrSpineCarac: DStr.Surface(Surf()).Surface() is null".to_string())?;

    // `gp_Pnt2d pp = cd->SetOfSurfData()->Value(i)
    // ->Interference(jf).PCurveOnSurf()->Value(p);` (`cxx:788`).
    let pp = sd
        .interference(jf)
        .p_curve_on_surf()
        .ok_or_else(|| "ChFi3d_ExtrSpineCarac: the interference has no p-curve".to_string())?
        .value(p);

    // `GeomAdaptor_Surface gs(fffil); P = fffil->Value(pp.X(), pp.Y());`
    // (`cxx:789-790`).
    let adaptor = SurfaceAdaptor::of_surface(Arc::clone(&surface));
    let point = adaptor
        .d0(pp.x(), pp.y())
        .ok_or_else(|| "ChFi3d_ExtrSpineCarac: the surface has no value".to_string())?;

    // `switch (gs.GetType())` (`cxx:793-833`).
    let (mut v, r) = match adaptor.get_type() {
        // `case GeomAbs_Cylinder: { gp_Cylinder cyl = gs.Cylinder(); R =
        // cyl.Radius(); ElSLib::D1(pp.X(), pp.Y(), cyl, Pbid, Vbid, V); }`
        // (`cxx:794-801`).
        GeomAbsSurfaceType::Cylinder => {
            let cyl = adaptor
                .cylinder()
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: gs.Cylinder()".to_string())?;
            let r = cyl.radius();
            let (_, _, dv) = adaptor
                .d1(pp.x(), pp.y())
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ElSLib::CylinderD1".to_string())?;
            (dv, r)
        }
        // `case GeomAbs_Torus: { gp_Torus tor = gs.Torus(); R =
        // tor.MinorRadius(); ElSLib::D1(pp.X(), pp.Y(), tor, Pbid, V, Vbid); }`
        // (`cxx:802-807`).
        GeomAbsSurfaceType::Torus => {
            let tor = surface
                .gp_torus()
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: gs.Torus()".to_string())?;
            let r = tor.minor_radius();
            let (_, du, _) = adaptor
                .d1(pp.x(), pp.y())
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ElSLib::TorusD1".to_string())?;
            (du, r)
        }
        // `default:` (`cxx:808-832`).
        _ => {
            let spine = cd
                .spine()
                .cloned()
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: the stripe has no spine".to_string())?;
            let guard = spine
                .read()
                .map_err(|_| "ChFi3d_ExtrSpineCarac: poisoned spine".to_string())?;
            // `nbelspine = sp->NbEdges(); if (nbelspine == 1) hels =
            // sp->ElSpine(1); else hels = sp->ElSpine(p);` (`cxx:808-820`).
            let hels = if guard.nb_edges() == 1 {
                guard
                    .el_spines()
                    .first()
                    .cloned()
                    .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ChFiDS_Spine::ElSpine(1)".to_string())?
            } else {
                guard
                    .el_spine(p)
                    .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ChFiDS_Spine::ElSpine(p)".to_string())?
            };
            // `if (fsp->IsConstant()) R = fsp->Radius(); else R =
            // fsp->Law(hels)->Value(p);` (`cxx:821-829`).
            let r = if guard.is_fil_constant() {
                guard
                    .fil_radius()
                    .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ChFiDS_FilSpine::Radius".to_string())?
            } else {
                let law =
                    guard
                        .fil_law(&hels)
                        .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ChFiDS_FilSpine::Law".to_string())?;
                let law: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
                law.value(p)
            };
            // `hels->D1(p, Pbid, V);` (`cxx:830`).
            let (_, dv) = hels
                .eval_d1(p)
                .ok_or_else(|| "ChFi3d_ExtrSpineCarac: ChFiDS_ElSpine::D1".to_string())?;
            (dv, r)
        }
    };

    // `V.Normalize(); if (sens == 1) { V.Reverse(); }` (`cxx:834-837`).
    // `gp_Vec::Normalize` divides by the modulus and leaves the vector
    // untouched at zero; the port keeps that and only reverses.
    let _ = v.normalize();
    if sens == 1 {
        v.reverse();
    }
    Ok((point, v, r))
}
