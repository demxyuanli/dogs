//! IGES read path from `IGESControl_Reader` through `IGESFile_Read` and
//! `IGESToBRep` transfer.
//!
//! `XSAlgo_ShapeProcessor::ProcessShape` (`IGESToBRep_Actor.cxx:204-218`)
//! runs the ported `FixShape` steps in `iges_process_shape`.

mod iges_control_reader;
mod iges_process_shape;
mod create_solids;
mod iges_data_iges_model;
mod iges_surface_params;
mod iges_data_param_reader;
mod iges_file_read;
mod liriges;
mod analiges;
mod structiges;
mod igesread;
mod iges_select_work_library;
mod iges_to_brep_actor;
mod basis_parameter_curve;
mod iges_to_brep_basic_curve;
mod conic_arc;
mod copious_data;
mod spline_curve;
mod spline_surface;
mod iges_to_brep_basic_surface;
mod iges_to_brep_curve_and_surface;
mod iges_to_brep_topo_curve;
mod brep_entity;
mod offset_curve;
mod curve_on_surface;
mod iges_boundary;
mod igessolid;
mod iges_to_brep_topo_surface;
mod ruled_surface;
mod tabulated_cylinder;
mod offset_surface;
mod bounded_surface;
mod gtransform;
mod elementary_bspline;
mod xscontrol_reader;

pub use iges_control_reader::read_iges_shape;
