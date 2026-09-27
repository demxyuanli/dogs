//! OCCT topology (TKBRep) — boundary representation data structures.
//!
//! Maps OCCT's TopoDS shape hierarchy to Rust:
//! ```text
//! Shape (enum: Vertex/Edge/Wire/Face/Shell/Solid/Compound)
//!   ├── TShape: geometric + topological data
//!   ├── Location: optional transform
//!   └── Orientation: FORWARD/REVERSED/INTERNAL/EXTERNAL
//! ```

pub mod abs;
pub mod tshape;
pub mod shape;
pub mod builder;
pub mod iterator;
pub mod tools;
pub mod mesh;
pub mod transform;
pub mod validate;
pub mod primitives;
pub mod tgeometry;

pub use abs::{ShapeType, Orientation};
pub use tshape::{TShape, VertexShape, EdgeShape, WireShape, FaceShape, ShellShape, SolidShape, CompoundShape};
pub use shape::{TopoShape, Vertex, Edge, Wire, Face, Shell, Solid, Compound};
pub use builder::TopoBuilder;
pub use iterator::{ShapeIterator, cumulated_children};
pub mod model;
pub mod topexp;
pub mod brep_tool;
pub mod brep_tools;
pub mod brep_check;

// Phase 3 modules (filled in by task agents).
pub mod wireframe;
pub mod shape_mesh;
pub mod bbox_from_geometry;
pub mod face_face;
pub mod edge_split;
pub mod shell_check;
pub mod solid_union;
pub mod mesh_to_brep;
pub mod brep_exchange;
pub mod brep_scene;
pub mod brep_extrema;
pub mod brep_face_intersect;
pub mod shape_custom_surface;
pub mod shape_fix_compose_shell;
pub mod shape_analysis;
pub mod brep_gprop;
pub mod sweep;
pub mod brep_surface;
pub mod topo_tools_full;
pub mod boolean_ops;
pub mod brep_measure;
// Phase 4 modules.
pub mod brep_builder_api;
pub mod fillet;
pub mod loft;
pub mod step;
pub mod hlr;
pub mod vrml;
pub mod sweep_revolve;
pub mod shape_ops;
pub mod shape_naming;
pub mod brep_pipe;
pub mod brep_faces;
pub mod brep_sketch;
pub mod brep_compare;
pub mod brep_assembly;
pub mod xcaf;
pub mod shape_checks;
pub mod geometry_query;
pub mod brep_shell;
pub mod brep_connect;
pub mod brep_projection;
pub mod brep_pattern;
pub mod shape_metrics;
pub mod render_svg;
pub mod iges;
pub mod brepmesh;
mod bop_builder_core;
mod bop_builder_planar_geom;
mod bop_builder_planar_trace;
mod bop_builder_planar_weld;
mod bop_builder_planar;
mod bop_builder_dispatch;
mod bop_builder_repair;
mod bop_builder_report;
mod bop_builder_heal;
mod bop_builder_splitapi;
pub mod bop_builder;
pub mod inttools;
pub mod intpatch;
pub mod intwalk;
pub mod bop_curved;
pub mod fillet_edge;
pub mod brep_offset;
pub mod bincaf;
pub mod gltf;
pub mod brepfeat;
pub mod fillet_var;
pub mod fillet_curved;
pub mod shhealing;
pub mod xmlcaf;
pub mod rwmesh;
pub mod brep_builder_full;
pub mod brep_lib_make_face;
pub mod brep_lib_make_wire;
pub mod brep_lib_same_parameter;
pub mod viz_scene;
pub mod draw;
pub mod gprop_analytic;
pub mod brep_gprop_full;
pub mod meshing;
// Phase 15 modules (precise NURBS boolean — mechanical wave).
pub mod bopds_tools;
mod bopds_pave;
mod bopds_cb;
mod bopds_tree;
mod bopds_faceinfo;
mod bopds_ds_interf;
mod bopds_iter;
pub mod bopds;
pub mod bopds_ff;
pub mod builder_area;
pub mod builder_face;
pub mod shell_splitter;
mod shell_splitter_block;
pub mod wire_splitter;
mod wire_splitter_block;
pub mod connexity_block;
pub mod bopalgo_options;
pub mod bop_hist;
// Phase 16 modules (precise NURBS boolean — IntTools mechanical layer + pcurve).
pub mod inttools_data;
pub mod inttools_range;
pub mod pcurve;
// Shared 2D sample-count helpers (`Geom2dAdaptor_Curve::NbSamples` /
// `Geom2dInt_Geom2dCurveTool::NbSamples`).
mod curve_sampling_2d;
// Phase 16b — IntTools sample/localize + full pcurve.
pub mod inttools_sample;
pub mod pcurve_full;
// Phase 16c — IntTools static root/parameter helpers.
pub mod inttools_roots;
// Phase 17 modules (precise NURBS boolean — wave C1: curve-surface intersect + 2D classify).
pub mod intcurvesurface;
pub mod intcurvesurface_poly;
pub mod int_curves_face;
pub mod brep_class;
pub mod brep_class3d;
pub mod fclass2d;
// Phase 17b — wave C1 core: bean-face intersector + edge-face.
mod bean_face_range;
mod bean_face_sample;
mod bean_face_distance;
mod bean_face_grid;
mod bean_face_localize;
mod bean_face_kind;
mod bean_face_analytic;
mod bean_face_exact;
pub mod bean_face;
mod edge_face_kind;
mod edge_face_type;
pub mod edge_face;
// Phase 18 modules (precise NURBS boolean — wave C2a: edge-edge + face-face + algo tools).
pub mod edge_edge;
pub mod int_face_face;
pub mod geom_int;
pub mod int_tools_wline;
pub mod int_tools_lines;
pub mod algo_tools;
// Phase 18b — full IntTools_Context + AlgoTools2D.
pub mod int_tools_full;
pub mod boptools_2d;
// Phase 18c — IntTools_Curve data class.
pub mod int_curve;
// Phase 19 modules (precise NURBS boolean — wave C2b: PaveFiller core).
pub mod pave_filler;
pub mod pave_intersect;
pub mod pave_blocks;
pub mod pave_common;
pub mod pave_ff;
mod pave_ff_exist;
mod pave_ff_update;
mod pave_ff_paves;
mod pave_self;
pub mod bopalgo_tools;
mod bopalgo_tools_wires;
mod bopalgo_tools_class;
mod algo_tools3d;
mod algo_tools_face;
mod algo_tools_range;
mod pave_ef;
mod pave_ef_perform;
mod pave_ee_aux;
mod pave_ee_perform;
mod pave_ff_misc;
mod pave_new;
mod pave_split;
mod builder_face_occt;
pub(crate) mod geom_bnd_lib_elclib2d;
// The former `geom_bnd_lib_elclib2d_{d2,dn,param}` submodules were deleted:
// they were a *second* port of the 2D `ElCLib` derivative/parameter surface
// (`ElCLib::D1/D2/D3/DN/Parameter` for `gp_Lin2d`/`Circ2d`/`Elips2d`/`Hypr2d`/
// `Parab2d`) that nothing called, while `occt_core::elib::clib2d` already
// carries the faithful implementation the rest of the crate uses
// (`clib2d.rs:126-531`). They accounted for 45 `never used` warnings.
mod geom_bnd_lib_elclib_to3d;
mod geom_bnd_lib_inf2d;
mod geom_bnd_lib_line2d;
mod geom_bnd_lib_circle2d;
mod geom_bnd_lib_ellipse2d;
mod geom_bnd_lib_hyperbola2d;
mod geom_bnd_lib_parabola2d;
mod geom_bnd_lib_other2d;
mod geom_bnd_lib_sample2d;
mod geom_bnd_lib_bezier2d;
mod geom_bnd_lib_bspline2d;
mod geom_bnd_lib_offset2d;
mod geom_bnd_lib_curve2d;
mod geom_bnd_lib_spline_helpers;
mod geom_bnd_lib_inf3d;
mod geom_bnd_lib_circle3d;
mod geom_bnd_lib_curve3d;
mod geom_bnd_lib_analytic3d;
mod geom_bnd_lib_surface3d;
pub mod brep_bnd_lib;
mod bnd_lib_add2d;
mod brep_uv_bounds;
mod bop_box2d_tree;
mod bvh_box2d;
mod bvh_traverse2d;
mod bop_box_selector;
mod bop_pair_selector;
mod bnd_tools;
mod pave_force_ee;
mod pave_force_ef;
mod pave_vf;
mod pave_pcurves;
mod pave_de;
mod pave_vv;
mod pave_ve;
mod pave_shrunk;
mod pave_split_blocks;
mod pave_ff_se;
mod pave_update_sd;
mod int_tools_vertex_line;
mod int_tools_curve_box;
mod int_tools_segpln;
mod pave_ff_perform;
mod pave_ff_pave_put;
mod pave_ff_pave_stick;
mod pave_ff_pave_bound;
mod pave_ff_make;
mod pave_ff_rebuild;
mod pave_ff_shared;
mod pave_ff_face_info;
mod pave_ff_exist_ve;
mod pave_ff_unused;
mod pave_ff_is_exist;
mod pave_ff_exist_es;
mod pave_ff_exist_onin;
mod pave_ff_post;
mod pave_ff_ef_pnts;
// Phase 20 modules (precise NURBS boolean — wave C2b-2a: BOPAlgo_Builder rebuild).
pub mod bop_builder2;
pub mod bop_bop;
mod bop_build_bop;
pub mod brep_algo_api;
pub mod bop_cells;
pub mod bop_remove_features;
pub mod bop_build_faces;
pub mod bop_build_common;
pub mod bop_build_solids;
pub mod builder_solid;
mod bop_occt_util;
mod bop_split_seam;
mod bop_split_to_reverse;
mod bop_draft_face;
mod bop_tools_set;
mod bop_split_faces_occt;
mod bop_draft_solid_occt;
mod bop_classify_occt;
mod bop_fill_in3d;
mod bop_split_solids_occt;
mod bop_fill_internals_occt;
mod bop_aabb_faces;
mod bop_connexity_faces;
mod bop_geomlib_closed;
mod bop_images_solids;
mod bopalgo_tools_blocks;
mod bop_pair_sd;
mod bop_same_domain_faces;
mod bop_fill_internal_verts;
