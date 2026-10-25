//! `TopOpeBRepTool` — the utility layer of the classic boolean engine.
//!
//! Source package: `TKBool/TopOpeBRepTool` (OCCT 8.0.0, 20,036 lines).
//!
//! Ported so far:
//! - `TopOpeBRepTool_ShapeTool`'s orientation family (`ShapesSameOriented`,
//!   `FacesSameOriented`, `SurfacesSameOriented`, `EdgesSameOriented`,
//!   `CurvesSameOriented`) and `EdgeData` ([`shape_tool`]).
//! - The `BRepAdaptor_Surface` / `BRepAdaptor_Curve` slice those read
//!   ([`adaptor`]).
//! - `LProp_CurveUtils` ([`lprop`]), the template `BRepLProp_CLProps` is
//!   instantiated from.
//! - `FUN_tool_projPonS` ([`project`]).
//! - `TopOpeBRepTool_TOOL`'s closure / orientation family ([`closed`]):
//!   `OriinSor`, `OriinSorclosed`, `ClosedE`, `ClosedS`, `IsClosingE`,
//!   `Vertices`, `Vertex`, `ParE`, `OnBoundary`.
//! - `TopOpeBRepTool_TOOL`'s tangent family and `FUN_tool_curvesSO`
//!   ([`tangent`]).
//! - `FUN_tool_projPonC` / `FUN_tool_value` / `FUN_tool_parE` ([`project_curve`]).
//! - `FUN_tool_maxtol` / `FUN_tool_nbshapes` ([`measures`]).
//! - `TopOpeBRepTool_GEOMETRY`'s type predicates and `FUN_tool_ngS`
//!   ([`geometry`]).
//! - `TopOpeBRepTool_OutCurveType` / `TopOpeBRepTool_GeomTool`
//!   ([`out_curve_type`], [`geom_tool`]).
//! - `TopOpeBRepTool_TOOL::SplitE` with `FUN_tool_sortVonE` and
//!   `TrslUVModifE` ([`split_e`]).
//! - `TopOpeBRepTool_FuseEdges` ([`fuse_edges`]): the useless-vertex chain
//!   fusion `Builder::End` runs after a boolean operation.
//! - `TopOpeBRepTool_PurgeInternalEdges` ([`purge_internal_edges`]): drops the
//!   edges no face carries, rebuilding through `BRepTools_Substitution`.

pub mod adaptor;
pub mod box_sort;
pub mod c2d;
pub mod classi;
pub mod closed;
pub mod connexity;
pub mod corr_iso;
pub mod curvature;
pub mod curve_tool;
pub mod face;
pub mod fuse_edges;
pub mod geometry;
pub mod geom_tool;
pub mod h_box_tool;
pub mod iso;
pub mod lprop;
pub mod make_transition;
pub mod matter;
pub mod measures;
pub mod mk_shell;
pub mod mk_tondg_e;
pub mod ng_app;
pub mod out_curve_type;
pub mod project;
pub mod project_curve;
pub mod purge;
pub mod purge_internal_edges;
pub mod reguw;
pub mod regus;
pub mod regularize_s;
pub mod regularize_w;
pub mod sc;
pub mod shape_classifier;
pub mod solid_classifier;
pub mod shape_explorer;
pub mod shape_tool;
pub mod split_e;
pub mod tangent;
pub mod tol;
pub mod topology;
pub mod uv_app;
pub mod wire_to_faces;

pub use adaptor::{CurveAdaptor, SurfaceAdaptor};
pub use box_sort::BoxSort;
pub use classi::{
    fun_add_owlw, Classi, DIFF, ONE_IN_TWO, SAME, TWO_IN_ONE, UNKNOWN,
};
pub use connexity::Connexity;
pub use face::ToolFace;
pub use closed::{
    basis_surface, check_edge_length, closed_e, closed_s, is_closing_e, is_closing_e_in_wire,
    on_boundary, ori_code, ori_in_sor, ori_in_sor_closed, par_e, vertex, vertices, CLOSING,
    EXTERNAL, FORWARD, INTERNAL, REVERSED,
};
pub use corr_iso::{CorrIso, DECREASE, INCREASE, INIT, SPLITEDGE};
pub use curvature::{
    curv_e, curv_f, edge_on_face, get_geom_xx, get_stp_3d_f, getduv, is_quad_e, is_quad_f, min_duv,
    nt, stuv_f, xx,
};
pub use make_transition::MakeTransition;
pub use geometry::{
    fun_quad_ct, fun_tool_closed, fun_tool_closed_s, fun_tool_cylinder, fun_tool_is_closed_s,
    fun_tool_line, fun_tool_ng_s, fun_tool_onapex, fun_tool_plane, fun_tool_quad,
    fun_tool_quad_surface,
};
pub use iso::{
    basis_curve2d, iso_on_clo, iso_on_clo_c2df, par_iso, uv_iso_c2df, uv_iso_curve,
    uv_iso_edge_face, uvf, UvIso,
};
pub use geom_tool::GeomTool;
pub use h_box_tool::HBoxTool;
pub use shape_explorer::ShapeExplorer;
pub use curve_tool::CurveTool;
pub use lprop::{CLProps, LPropStatus};
pub use measures::{fun_tool_maxtol, fun_tool_maxtol_into, fun_tool_nbshapes};
pub use out_curve_type::OutCurveType;
pub use project::fun_tool_proj_pon_f;
pub use purge::{
    correct_on_uv_iso, fun_connex_x, fun_tool_closed_w, make_faces, purge_closing_edges,
    purge_closing_edges_faces,
};
pub use purge_internal_edges::PurgeInternalEdges;
pub use sc::{fsc_get_psc, fsc_state_eon_face, fsc_state_pon_face};
pub use split_e::{fun_tool_sort_von_e, split_e, trsl_uv_modif_e};
pub use shape_classifier::ShapeClassifier;
pub use solid_classifier::SolidClassifier;
pub use project_curve::{
    fun_tool_bounds, fun_tool_par_e, fun_tool_par_e_tol, fun_tool_proj_pon_c,
    fun_tool_proj_pon_c_range, fun_tool_proj_pon_c_tol, fun_tool_proj_pon_e,
    fun_tool_proj_pon_e_tol, fun_tool_value,
};
pub use shape_tool::{
    closed, curves_same_oriented, edge_data, edges_same_oriented, faces_same_oriented, pnt_of,
    resolution3d, resolution3d_face, resolution3d_u, resolution3d_v, shapes_same_oriented,
    surfaces_same_oriented, tolerance_of, EdgeDataResult,
};
pub use tangent::{
    curves_same_oriented_at, curves_same_oriented_point, fun_tool_dir_c, fun_tool_dir_c_adaptor,
    fun_tool_tg_geom_e, tg_geom_e, tg_geom_e_curve, tg_inside,
};
pub use tol::{
    ftol_face_tolerances, ftol_face_tolerances2d, ftol_face_tolerances3d_boxes,
    ftol_face_tolerances3d_faces,
};
pub use topology::{
    fun_tool_comparebndkole, fun_tool_ebound_f, fun_tool_eshared, fun_tool_etg_f,
    fun_tool_etg_oo_e, fun_tool_find_ancestor, fun_tool_find_p_in_bac,
    fun_tool_find_p_in_e,
    fun_tool_find_par_in_bac, fun_tool_find_par_in_e, fun_tool_in_s, fun_tool_ng_geom_f,
    fun_tool_ng_geom_f_par, fun_tool_ng_geom_f_par_tol, fun_tool_orient_ein_f,
    fun_tool_orient_ein_f_forward, fun_tool_orient_vin_e, fun_tool_par_f, fun_tool_par_f_tol,
    fun_tool_par_on_e_f, fun_tool_par_on_e_f_tol, fun_tool_par_von_e, fun_tool_pcurveon_f,
    fun_tool_pin_c, fun_tool_pin_c_range, fun_tool_sta_pin_e, fun_tool_sta_pin_e_tol,
    fun_tool_value_uv, FIRST, LAST,
};
pub use uv_app::{fun_nearest_iso, fun_tool_getdxx, fun_tool_isobounds, uv_app};
pub use mk_shell::{mk_shell, remove_from_list};
pub use fuse_edges::FuseEdges;
pub use regus::Regus;
pub use regularize_s::regularize_shells;
pub use ng_app::{n_c2d_inside_s, ng_app, out_uv_bounds, try_ng_app, try_ori_e_in_f};
