//! `IGESToBRep_TopoSurface::TransferTopoSurface` for entities 144, 108, and 120
//! (`IGESToBRep_TopoSurface.cxx:184-232`).

mod face_build;
mod face_transform;
mod perforate;
mod plane;
mod plane_location;
mod param_shift;
mod revolution;
mod revolution_axis;
mod transfer;
mod trimmed_surface;
mod uv_map;

pub(in crate::iges_read) use revolution::extract_curve_3d;
pub(in crate::iges_read) use transfer::transfer_topo_surface;
pub(in crate::iges_read) use uv_map::uv_map;

pub(in crate::iges_read::iges_to_brep_topo_surface) use face_build::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use face_transform::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use perforate::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use plane::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use plane_location::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use param_shift::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use revolution::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use revolution_axis::*;
pub(in crate::iges_read::iges_to_brep_topo_surface) use trimmed_surface::*;
