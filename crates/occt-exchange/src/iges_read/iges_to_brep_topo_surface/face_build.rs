use std::sync::Arc;


use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use occt_topo::builder::TopoBuilder;
use occt_topo::shape::Wire;


/// `BRep_Builder::MakeFace` + `B.UpdateFace(plane, surf, TopLoc_Location(),
/// Precision::Confusion())` (`IGESToBRep_TopoSurface.cxx:302-303`, `:1641`),
/// `BRepLib_MakeFace(surf, Precision::Confusion())` (`cxx:308`, `:1186`) and
/// `BRepBuilderAPI_MakeFace::Init(surf, ..., Precision::Confusion())`
/// (`cxx:811-815`, `:961-966`): every face this module creates carries
/// `Precision::Confusion()`, not the `BRep_TFace` default `RealEpsilon()`.
/// `BRepBndLib::Add` feeds that tolerance to `BndLib_AddSurface::Add`
/// (`BRepBndLib.cxx:109`), so it sets `Bnd_Box::Gap` and therefore
/// `Prs3d::GetDeflection`.
pub(super) fn make_face(surface: Arc<dyn Surface>, wires: &[Wire]) -> occt_topo::shape::Face {
    TopoBuilder::new().make_face_with_tolerance(surface, wires, CONFUSION)
}
