//! IGES read path from `IGESControl_Reader` through `IGESFile_Read` and
//! `IGESToBRep` transfer.
//!
//! `XSAlgo_ShapeProcessor::ProcessShape` (`IGESToBRep_Actor.cxx:204-218`)
//! runs the ported `FixShape` steps in `iges_process_shape`.

mod analiges;
mod basis_parameter_curve;
mod bounded_surface;
mod brep_entity;
mod conic_arc;
mod copious_data;
pub(crate) mod create_solids;
mod curve_on_surface;
mod elementary_bspline;
mod gtransform;
mod iges_boundary;
mod iges_control_reader;
mod iges_data_iges_model;
mod iges_data_param_reader;
mod iges_file_read;
mod iges_process_shape;
mod iges_select_work_library;
mod iges_surface_params;
mod iges_to_brep_actor;
mod iges_to_brep_basic_curve;
mod iges_to_brep_basic_surface;
mod iges_to_brep_curve_and_surface;
mod iges_to_brep_topo_curve;
mod iges_to_brep_topo_surface;
mod igesread;
mod igessolid;
mod liriges;
mod offset_curve;
mod offset_surface;
mod ruled_surface;
mod spline_curve;
mod spline_surface;
mod structiges;
mod tabulated_cylinder;
mod xscontrol_reader;

pub use iges_control_reader::read_iges_shape;
