//! `ToricRotule` (`ChFi3d_FilBuilder_C2.cxx:87-122`): the torus corner test.
//!
//! "Test if it is a particular case of torus routine. Three planes with two
//! constant incident fillets of the same radius and the third face
//! perpendicular to two others are required."

use crate::chfi3d::builder::fil_builder::perform_surf::fil_spine;
use crate::chfi_ds::StripeHandle;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `static bool ToricRotule(fac, s1, s2, c1, c2)` (`cxx:87-122`).
pub(super) fn toric_rotule(
    fac: &SurfaceAdaptor,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    c1: &StripeHandle,
    c2: &StripeHandle,
) -> Result<bool, String> {
    let tolesp = 1.0e-7;

    // `occ::handle<ChFiDS_FilSpine> sp1 =
    // occ::down_cast<ChFiDS_FilSpine>(c1->Spine()); if (sp1.IsNull() ||
    // sp2.IsNull()) return false;` (`cxx:95-101`).
    let sp1 = c1
        .read()
        .map_err(|_| "ToricRotule: poisoned stripe".to_string())?
        .spine()
        .cloned();
    let sp2 = c2
        .read()
        .map_err(|_| "ToricRotule: poisoned stripe".to_string())?
        .spine()
        .cloned();
    let (Some(sp1), Some(sp2)) = (fil_spine(sp1.as_ref()), fil_spine(sp2.as_ref())) else {
        return Ok(false);
    };

    // `if (!sp1->IsConstant() || !sp2->IsConstant()) return false;`
    // (`cxx:102-105`).
    let (r1, r2) = {
        let g1 = sp1
            .read()
            .map_err(|_| "ToricRotule: poisoned spine".to_string())?;
        let g2 = sp2
            .read()
            .map_err(|_| "ToricRotule: poisoned spine".to_string())?;
        if !g1.is_fil_constant() || !g2.is_fil_constant() {
            return Ok(false);
        }
        (
            g1.fil_radius()
                .ok_or_else(|| "ToricRotule: ChFiDS_FilSpine::Radius()".to_string())?,
            g2.fil_radius()
                .ok_or_else(|| "ToricRotule: ChFiDS_FilSpine::Radius()".to_string())?,
        )
    };

    // `if ((fac.GetType() != GeomAbs_Plane) || (s1.GetType() !=
    // GeomAbs_Plane) || (s2.GetType() != GeomAbs_Plane)) return false;`
    // (`cxx:106-111`).
    if fac.get_type() != GeomAbsSurfaceType::Plane
        || s1.get_type() != GeomAbsSurfaceType::Plane
        || s2.get_type() != GeomAbsSurfaceType::Plane
    {
        return Ok(false);
    }

    // `gp_Dir df = fac.Plane().Position().Direction();` and the two others
    // (`cxx:111-113`).
    let df = fac
        .plane()
        .ok_or_else(|| "ToricRotule: BRepAdaptor_Surface::Plane()".to_string())?
        .position()
        .direction();
    let ds1 = s1
        .plane()
        .ok_or_else(|| "ToricRotule: BRepAdaptor_Surface::Plane()".to_string())?
        .position()
        .direction();
    let ds2 = s2
        .plane()
        .ok_or_else(|| "ToricRotule: BRepAdaptor_Surface::Plane()".to_string())?
        .position()
        .direction();

    // `if (std::abs(df.Dot(ds1)) >= tolesp || std::abs(df.Dot(ds2)) >=
    // tolesp) return false;` (`cxx:114-117`).
    if df.dot(&ds1).abs() >= tolesp || df.dot(&ds2).abs() >= tolesp {
        return Ok(false);
    }

    // `double r1 = sp1->Radius(); double r2 = sp2->Radius(); return
    // std::abs(r1 - r2) < tolesp;` (`cxx:119-121`).
    Ok((r1 - r2).abs() < tolesp)
}
