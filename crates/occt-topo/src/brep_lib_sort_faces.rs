//! `BRepLib::SortFaces` (`BRepLib.cxx:2894-2951`) and its mirror
//! `BRepLib::ReverseSortFaces` (`cxx:2955-...`).
//!
//! `SortFaces` buckets the faces of a shape by surface type and appends the
//! buckets in a fixed order - planes, cylinders, cones, spheres, toruses,
//! everything else, and finally the faces without a surface ("triangulated") -
//! so an offset builds its analytic faces before the free-form ones.
//! `BRepOffset_MakeOffset::MakeOffsetFaces` (`BRepOffset_MakeOffset.cxx:1214`)
//! is the caller.
//!
//! Deviation, and the only one:
//!
//! * `BRep_Tool::Surface(F, l)` answers the surface *with* the face's location
//!   baked out into `l`; OCCT then hands `S` - still un-located - to
//!   `GeomAdaptor_Surface`, whose type only depends on the surface itself, so
//!   the location never affects the bucket. The port reads
//!   [`BRepTool::face_surface`] (the un-located surface, the same one) and only
//!   uses it for `rectangular_trimmed_basis` and `geometric kind`.
//! * `LF.Clear()` first (`cxx:2896`): the port returns a fresh `Vec`, so the
//!   caller's previous contents are dropped rather than cleared in place.

use std::sync::Arc;

use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::abs::ShapeType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// The surface type `GeomAdaptor_Surface(S)` reports, after stripping a
/// `Geom_RectangularTrimmedSurface` (`cxx:2908-2912`).
fn adaptor_type(s: &Arc<dyn Surface>) -> GeomAbsSurfaceType {
    let basis = match s.rectangular_trimmed_basis() {
        Some(b) => b,
        None => s.clone(),
    };
    SurfaceAdaptor::of_surface(basis).get_type()
}

/// `BRepLib::SortFaces(Sh, LF)` (`cxx:2894-2951`).
pub fn sort_faces(shape: &TopoShape) -> Vec<TopoShape> {
    let mut l_tri: Vec<TopoShape> = Vec::new();
    let mut l_plan: Vec<TopoShape> = Vec::new();
    let mut l_cyl: Vec<TopoShape> = Vec::new();
    let mut l_con: Vec<TopoShape> = Vec::new();
    let mut l_sphere: Vec<TopoShape> = Vec::new();
    let mut l_tor: Vec<TopoShape> = Vec::new();
    let mut l_other: Vec<TopoShape> = Vec::new();

    for f in explorer_shapes(shape, ShapeType::Face) {
        let face = Face(f);
        let Some(s) = BRepTool::face_surface(&face) else {
            // `cxx:2939-2941`: no surface at all -> the triangulated bucket.
            l_tri.push(face.0);
            continue;
        };
        match adaptor_type(&s) {
            GeomAbsSurfaceType::Plane => l_plan.push(face.0),
            GeomAbsSurfaceType::Cylinder => l_cyl.push(face.0),
            GeomAbsSurfaceType::Cone => l_con.push(face.0),
            GeomAbsSurfaceType::Sphere => l_sphere.push(face.0),
            GeomAbsSurfaceType::Torus => l_tor.push(face.0),
            _ => l_other.push(face.0),
        }
    }

    let mut lf: Vec<TopoShape> = Vec::new();
    lf.extend(l_plan);
    lf.extend(l_cyl);
    lf.extend(l_con);
    lf.extend(l_sphere);
    lf.extend(l_tor);
    lf.extend(l_other);
    lf.extend(l_tri);
    lf
}
