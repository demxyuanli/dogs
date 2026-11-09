//! Locals that `LocOpe_SplitDrafts::Perform` (`cxx:105-1436`) carries between
//! its stages.
//!
//! The reference declares all of these at the top of `Perform` or right before
//! the stage that first uses them; grouping them here keeps the stage modules
//! readable without changing the order in which they are filled.

use std::sync::Arc;

use occt_core::gp::GpAx1;
use occt_geom::Surface;

use crate::geom_int::IntSS;
use crate::locope::shape_maps::{ShapeIndexedMap, ShapeSet, ShapeValueMap};
use crate::shape::{TopoShape, Wire};
use crate::topo_tools::AncestorMap;

/// Scratch state of one `Perform` call.
pub(crate) struct PerformState {
    /// `NewSg` / `NewSd` (`cxx:172-173`): the two draft planes.
    pub new_sg: Arc<dyn Surface>,
    pub new_sd: Arc<dyn Surface>,
    /// `NormalFg` (`cxx:155`): the axis along which `Choice == 2` walks.
    pub normal_fg: GpAx1,
    /// `theLinePipe` (`cxx:174`): `Geom_Line(NormalFg)`.
    pub line_pipe: Arc<dyn occt_geom::Curve>,
    /// `i2s` (`cxx:175`): `GeomInt_IntSS(NewSg, NewSd, Confusion)`.
    pub i2s: IntSS,
    /// `theW` (`cxx:181`): the wire, possibly re-split by `splw`.
    pub the_w: Wire,
    /// `Res` (`cxx:374`): the `LocOpe_Spliter` result.
    pub res: Option<TopoShape>,
    /// `theLeft` (`cxx:375`): `Spls.DirectLeft()`.
    pub the_left: Vec<TopoShape>,
    /// `MapW` (`cxx:385`):
    /// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`; `None` is the
    /// reference's null shape, and `IsBound` still answers `true` for it.
    pub map_w: ShapeValueMap<Option<TopoShape>>,
    /// `MapEV` (`cxx:405`): `IndexedDataMap<TopoDS_Shape, TopoDS_Shape>`.
    pub map_ev: ShapeIndexedMap<TopoShape>,
    /// `theMapEF` (`cxx:401`): `TopExp::MapShapesAndAncestors(Res, EDGE, FACE)`.
    pub the_map_ef: AncestorMap,
    /// `MapSg` / `MapSd` (`cxx:407-408`): the images on the left / right.
    pub map_sg: ShapeValueMap<Vec<TopoShape>>,
    pub map_sd: ShapeValueMap<Vec<TopoShape>>,
    /// `MapSonS` (`cxx:471`).
    pub map_son_s: ShapeValueMap<TopoShape>,
    /// `theMap` (`cxx:177`): the `NCollection_Map` reused as a "seen" set.
    pub the_map: ShapeSet,
    /// `thefaces` (`cxx:1103`): the drafted faces.
    pub the_faces: Vec<TopoShape>,
    /// `theSubs` (`cxx:1337`): filled by [`super::draft_faces`], read again by
    /// [`super::descendants`].
    pub subs: crate::brep_tools_substitution::BRepToolsSubstitution,
}
