//! `ChFi3d_KParticular`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:604-683`
//! (OCCT 8.0.0).
//!
//! `PerformSetOfKPart` asks this predicate per edge of a spine: only the edges
//! that answer `true` are handed to `ChFiKPart_ComputeData::Compute`
//! (`ChFi3d_Builder_2.cxx:3055-3057`).
//!
//! Divergences that are not heuristics:
//! - `occ::down_cast<ChFiDS_FilSpine>` is [`SpineClass::Fil`]; OCCT's
//!   `IsConstant(IE)` is [`crate::chfi_ds::Spine::is_fil_constant_on_edge`],
//!   whose `None` (an `IE` outside `1..=NbEdges`) OCCT cannot produce,
//! - `Adaptor3d_Surface::Plane()` / `Cylinder()` / `Cone()` are
//!   [`occt_geom::Surface::gp_pln`] / `gp_cylinder` / `gp_cone`; a type that
//!   reports `GeomAbs_Plane` but carries no `gp_Pln` is the null downcast OCCT
//!   never guards, the port answers `false` instead of dereferencing,
//! - `Spine->CurrentElementarySpine(IE)` caches the edge adaptor on the spine
//!   (`ChFiDS_Spine.cxx:866-874`), hence the write lock.

use occt_core::gp::{GpCylinder, GpDir, GpPln};
use occt_core::precision::ANGULAR;

use crate::chfi_ds::{SpineClass, SpineHandle};
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_KParticular(Spine, IE, S1, S2)` (`cxx:604-683`).
///
/// `true` when the edge `IE` of `Spine` is a straight line or a circle and one
/// of the two support faces is a plane that meets the other (plane, cylinder,
/// cone) along it — the configurations `ChFiKPart_ComputeData` builds
/// analytically instead of walking a blend.
pub fn ch_fi3d_k_particular(
    spine: &SpineHandle,
    ie: i32,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
) -> bool {
    // `occ::handle<ChFiDS_FilSpine> fs = occ::down_cast<ChFiDS_FilSpine>(Spine);
    // if (!fs.IsNull() && !fs->IsConstant(IE)) { return !bRet; }`
    // (`cxx:612-617`).
    {
        let guard = match spine.read() {
            Ok(g) => g,
            Err(_) => return false,
        };
        if guard.class() == SpineClass::Fil
            && guard.is_fil_constant_on_edge(ie) != Some(true)
        {
            return false;
        }
    }

    // `aST1 = S1.GetType(); aST2 = S2.GetType(); bIsPlane1 = (aST1 ==
    // GeomAbs_Plane); ... if (!(bIsPlane1 || bIsPlane2)) return !bRet;`
    // (`cxx:626-630`).
    let a_st1 = s1.get_type();
    let a_st2 = s2.get_type();
    let is_plane1 = a_st1 == GeomAbsSurfaceType::Plane;
    let is_plane2 = a_st2 == GeomAbsSurfaceType::Plane;
    if !(is_plane1 || is_plane2) {
        return false;
    }

    // `const BRepAdaptor_Surface& aS1 = (bIsPlane1) ? S1 : S2;` (`cxx:633-636`).
    let (a_s1, a_s2) = if is_plane1 { (s1, s2) } else { (s2, s1) };
    let a_st2 = a_s2.get_type();

    // `if (aST2 != GeomAbs_Plane && aST2 != GeomAbs_Cylinder && aST2 !=
    // GeomAbs_Cone) return !bRet;` (`cxx:639-642`).
    if a_st2 != GeomAbsSurfaceType::Plane
        && a_st2 != GeomAbsSurfaceType::Cylinder
        && a_st2 != GeomAbsSurfaceType::Cone
    {
        return false;
    }

    // `const BRepAdaptor_Curve& bc = Spine->CurrentElementarySpine(IE); aCT =
    // bc.GetType(); if (aCT != GeomAbs_Line && aCT != GeomAbs_Circle) return
    // !bRet;` (`cxx:644-648`).
    let a_ct = {
        let mut guard = match spine.write() {
            Ok(g) => g,
            Err(_) => return false,
        };
        match guard.current_elementary_spine(ie) {
            Some(c) => c.get_type(),
            // `Edges(IE)` out of range: OCCT indexes the spine list unchecked.
            None => return false,
        }
    };
    if a_ct != GeomAbsCurveType::Line && a_ct != GeomAbsCurveType::Circle {
        return false;
    }

    // `aPA = Precision::Angular();` (`cxx:650`).
    let a_pa = ANGULAR;

    // `if (aST2 == GeomAbs_Plane) { if (aCT == GeomAbs_Line) return bRet; }`
    // (`cxx:652-658`).
    if a_st2 == GeomAbsSurfaceType::Plane {
        if a_ct == GeomAbsCurveType::Line {
            return true;
        }
    } else if a_st2 == GeomAbsSurfaceType::Cylinder {
        // `const gp_Dir aD1 = aS1.Plane().Axis().Direction(); const gp_Dir aD2
        // = aS2.Cylinder().Axis().Direction();` (`cxx:661-662`).
        let (Some(a_d1), Some(a_d2)) = (plane_axis(a_s1), cylinder_axis(a_s2)) else {
            return false;
        };
        // `cxx:664-671`.
        if a_ct == GeomAbsCurveType::Line && a_d1.is_normal_tol(&a_d2, a_pa) {
            return true;
        } else if a_ct == GeomAbsCurveType::Circle && a_d1.is_parallel_tol(&a_d2, a_pa) {
            return true;
        }
    } else if a_st2 == GeomAbsSurfaceType::Cone {
        // `const gp_Dir aD1 = aS1.Plane().Axis().Direction(); const gp_Dir aD2
        // = aS2.Cone().Axis().Direction(); if (aCT == GeomAbs_Circle &&
        // aD1.IsParallel(aD2, aPA)) return bRet;` (`cxx:674-680`).
        let (Some(a_d1), Some(a_d2)) = (plane_axis(a_s1), cone_axis(a_s2)) else {
            return false;
        };
        if a_ct == GeomAbsCurveType::Circle && a_d1.is_parallel_tol(&a_d2, a_pa) {
            return true;
        }
    }

    // `return !bRet;` (`cxx:683`).
    false
}

/// `aS1.Plane().Axis().Direction()`.
fn plane_axis(s: &SurfaceAdaptor) -> Option<GpDir> {
    let pln: GpPln = s.surface()?.gp_pln()?;
    Some(*pln.axis().direction())
}

/// `aS2.Cylinder().Axis().Direction()`.
fn cylinder_axis(s: &SurfaceAdaptor) -> Option<GpDir> {
    let cyl: GpCylinder = s.surface()?.gp_cylinder()?;
    Some(*cyl.axis().direction())
}

/// `aS2.Cone().Axis().Direction()`.
fn cone_axis(s: &SurfaceAdaptor) -> Option<GpDir> {
    let cone = s.surface()?.gp_cone()?;
    Some(*cone.axis().direction())
}
