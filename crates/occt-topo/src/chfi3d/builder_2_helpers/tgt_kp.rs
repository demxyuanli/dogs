//! `TgtKP`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:373-398`
//! (OCCT 8.0.0).
//!
//! `PerformSetOfKPart` calls it twice per `SurfData` (`_2.cxx:3145-3146`) to
//! attach the tangents of the section to the extremities of the `ElSpine` it
//! is about to append.

use occt_core::gp::{GpPnt, GpVec};

use crate::chfi_ds::{SpineHandle, SurfData};

/// `TgtKP(CD, Spine, iedge, isfirst, ped, ded)` (`cxx:373-398`).
///
/// `ped` / `ded` are `PFirst` / `TFirst` (or `PLast` / `TLast`) of the caller:
/// the point and the *normalized* tangent of the support edge at the section
/// parameter of `CD` on `S1`.
///
/// Divergences that are not heuristics:
/// - `Spine->CurrentElementarySpine(iedge)` caches the edge adaptor on the
///   spine (`ChFiDS_Spine.cxx:866-874`), hence the write lock,
/// - `Edges(iedge)` out of range and a null `BRep_Tool::Curve` are the two
///   places OCCT indexes unchecked; the port answers `None`.
pub fn tgt_kp(
    cd: &SurfData,
    spine: &SpineHandle,
    iedge: i32,
    isfirst: bool,
) -> Option<(GpPnt, GpVec)> {
    // `double wtg = CD->InterferenceOnS1().Parameter(isfirst);`
    // (`cxx:380`).
    let wtg = cd.interference_on_s1().parameter(isfirst);

    let mut guard = spine.write().ok()?;
    let reversed = guard
        .edges(usize::try_from(iedge).ok()?)?
        .orientation()
        .is_reversed();
    let bc = guard.current_elementary_spine(iedge)?;

    // `if (Spine->Edges(iedge).Orientation() == TopAbs_FORWARD) bc.D1(wtg +
    // bc.FirstParameter(), ped, ded); else { bc.D1(-wtg + bc.LastParameter(),
    // ped, ded); ded.Reverse(); }` (`cxx:382-391`).
    let (ped, mut ded) = if !reversed {
        bc.d1(wtg + bc.first_parameter())?
    } else {
        bc.d1(-wtg + bc.last_parameter())?
    };
    if reversed {
        ded.reverse();
    }
    // `ded.Normalize();` (`cxx:393`).
    ded.normalize();
    Some((ped, ded))
}
