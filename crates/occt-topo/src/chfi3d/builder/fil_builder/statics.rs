//! The two file statics of `ChFi3d_FilBuilder.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:77-143`
//! (OCCT 8.0.0).

use std::sync::Arc;

use occt_geom::law::LawFunction;

use crate::chfi_ds::{ElSpine, SpineHandle};
use crate::shape::TopoShape;

/// `MaxRad(fsp, E)` (`cxx:77-110`): the radius of `E` when the contour is
/// constant there, otherwise the maximum over `ParAndRad` and the laws.
///
/// The commented-out tail of the OCCT body (`cxx:93-109`, a six-sample scan
/// over a `Law_Composite`) is inside `/* */` and is not translated; the live
/// answer is `MaxRadFromSeqAndLaws()`.
///
/// `ChFiDS_FilSpine::IsConstant(IE)` / `Radius(IE)` raise when the index is out
/// of range or the radius is not constant; the port answers `Err`.
pub(crate) fn max_rad(spine: &SpineHandle, e: &TopoShape) -> Result<f64, String> {
    let guard = spine
        .read()
        .map_err(|_| "MaxRad: poisoned spine".to_string())?;
    let ie = guard.index_of_edge(e);
    // `if (fsp->IsConstant(IE)) { return (fsp->Radius(IE)); }` (`cxx:84-87`).
    let constant = guard
        .is_fil_constant_on_edge(ie)
        .ok_or_else(|| "MaxRad: ChFiDS_FilSpine::IsConstant raised Standard_OutOfRange".to_string())?;
    if constant {
        return guard
            .fil_radius_of_edge_index(ie)
            .ok_or_else(|| "MaxRad: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string());
    }
    // `else { return (fsp->MaxRadFromSeqAndLaws()); }` (`cxx:88-91`).
    Ok(guard.max_rad_from_seq_and_laws())
}

/// `SimulParams(HGuide, fsp, MaxStep, Fleche)` (`cxx:112-143`).
///
/// `MaxStep` is five percent of the guideline length; `Fleche` is five percent
/// of the largest section radius (six samples of the composite law when the
/// contour is not constant).
///
/// `ChFiDS_FilSpine::Law(HGuide)` raises `Standard_DomainError` when the
/// guideline is not one of the contour's elementary spines; the port answers
/// `Err` there instead of panicking.
pub(crate) fn simul_params(guide: &Arc<ElSpine>, spine: &SpineHandle) -> Result<(f64, f64), String> {
    let fi = guide.first_parameter();
    let la = guide.last_parameter();
    let longueur = la - fi;
    let max_step = longueur * 0.05;
    let guard = spine
        .read()
        .map_err(|_| "SimulParams: poisoned spine".to_string())?;
    // `if (fsp->IsConstant()) { radiussect = fsp->Radius(); }`
    // (`cxx:124-127`).
    let radiussect = if guard.is_fil_constant() {
        guard
            .fil_radius()
            .ok_or_else(|| "SimulParams: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string())?
    } else {
        // `occ::handle<Law_Composite> lc = fsp->Law(HGuide);`
        // `for (int i = 0; i <= 5; i++) { w = fi + i*longueur*0.2; ... }`
        // (`cxx:130-141`).
        let lc = guard
            .fil_law(guide)
            .ok_or_else(|| "SimulParams: ChFiDS_FilSpine::Law raised Standard_DomainError".to_string())?;
        let guard = lc.read().unwrap_or_else(|e| e.into_inner());
        let mut radiussect = 0.0f64;
        for i in 0..=5 {
            let w = fi + f64::from(i) * longueur * 0.2;
            let temp = guard.value(w);
            if temp > radiussect {
                radiussect = temp;
            }
        }
        radiussect
    };
    Ok((max_step, radiussect * 0.05))
}
