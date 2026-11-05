//! `PutProfilAt` (`BRepFill_Evolved.cxx:2714-2768`): place the reference
//! profile, which lives in the `yOz` reference plane, at a point of an edge of
//! the spine so that its `x` axis follows the edge tangent.

use occt_core::gp::{GpAx3, GpDir, GpPnt, GpPnt2d, GpTrsf};
use occt_core::toploc::TopLocLocation;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::shape::{Edge, Face, TopoShape};

/// `PutProfilAt(ProfRef, AxeRef, E, F, AtStart)`
/// (`BRepFill_Evolved.cxx:2714-2768`).
///
/// `e` must be an edge of `f`: the tangent is read from the p-curve of `e` on
/// `f`. A missing p-curve is the reference's `Standard_ConstructionError`
/// (`cxx:2729-2731`).
pub fn put_profil_at(
    prof_ref: &TopoShape,
    axe_ref: &GpAx3,
    e: &TopoShape,
    f: &TopoShape,
    at_start: bool,
) -> Result<TopoShape, String> {
    // `cxx:2727-2731`.
    let (c2d, first, last) = curve_on_surface_range(&Edge(e.clone()), &Face(f.clone()))
        .ok_or_else(|| "ConstructionError in PutProfilAt".to_string())?;

    // `cxx:2733-2756`: point and first derivative at the end of the p-curve;
    // the derivative is reversed when the edge itself is reversed.
    let reversed = e.orientation() == Orientation::Reversed;
    let (p, d1): (GpPnt2d, _) = if reversed {
        let (p, d) = if !at_start {
            c2d.d1(first)
        } else {
            c2d.d1(last)
        };
        (p, d.reversed())
    } else if !at_start {
        c2d.d1(last)
    } else {
        c2d.d1(first)
    };

    // `cxx:2757-2758`.
    let p3d = GpPnt::new(p.x(), p.y(), 0.0);
    let v3d = GpDir::new(d1.x(), d1.y(), 0.0).map_err(|e| e.to_string())?;

    // `cxx:2759`: `gp_Ax3 Ax(P3d, gp::DZ(), V3d)`.
    let dz = GpDir::new(0.0, 0.0, 1.0).map_err(|e| e.to_string())?;
    let ax = GpAx3::from_point_dir_x(p3d, dz, &v3d).map_err(|e| e.to_string())?;

    // `cxx:2760-2762`: `gp_Trsf Trans; Trans.SetTransformation(Ax, AxeRef);`.
    let mut trans = GpTrsf::identity();
    trans.set_transformation_from_to(&ax, axe_ref);

    // `cxx:2763-2765`: `ProfRef.Moved(TopLoc_Location(Trans))`.
    let location = TopLocLocation::composed(&trans, prof_ref.location().clone());
    Ok(prof_ref.located(&location))
}
