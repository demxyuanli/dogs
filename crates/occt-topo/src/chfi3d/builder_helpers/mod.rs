//! The free functions of `ChFi3d_Builder_0.cxx` that the `ChFi3d_Builder`
//! files call.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx`
//! (OCCT 8.0.0), declarations in `ChFi3d_Builder_0.hxx:26-197`.
//!
//! | OCCT | here |
//! |---|---|
//! | `ChFi3d_InPeriod` (`cxx:93-111`) | [`period::in_period`] |
//! | `ChFi3d_SearchPivot` (`cxx:4509-4544`) | [`period::search_pivot`] |
//! | `ChFi3d_Boite` 2 / 4 points (`cxx:242-284`) | [`box2d::boite`], [`box2d::boite_of_four`] |
//! | `ChFi3d_conexfaces` (`cxx:473-499`) | [`faces::conex_faces`] |
//! | `ChFi3d_EdgeState` (`cxx:506-580`) | [`faces::edge_state`] |
//! | `ChFi3d_evalconti` (`cxx:581-599`) | [`faces::eval_conti`] |
//! | `ChFi3d_cherche_face1` (`cxx:5334-5353`) | [`vertex_edges::cherche_face1`] |
//! | `ChFi3d_cherche_element` (`cxx:5361-5405`) | [`vertex_edges::cherche_element`] |
//! | `ChFi3d_cherche_edge` (`cxx:5410-5462`) | [`vertex_edges::cherche_edge`] |
//! | `ChFi3d_cherche_vertex` (`cxx:5716-5744`) | [`vertex_edges::cherche_vertex`] |
//! | `ChFi3d_nbface` (`cxx:5468-5490`) | [`vertex_edges::nbface`] |
//! | `ChFi3d_edge_common_faces` (`cxx:5498-5521`) | [`vertex_edges::edge_common_faces`] |
//! | `ChFi3d_AngleEdge` (`cxx:5527-5549`) | [`vertex_edges::angle_edge`] |
//! | `ChFi3d_ChercheBordsLibres` (`cxx:5557-5611`) | [`topology_stats::cherche_bords_libres`] |
//! | `ChFi3d_NbNotDegeneratedEdges` (`cxx:5619-5631`) | [`topology_stats::nb_not_degenerated_edges`] |
//! | `ChFi3d_NbSharpEdges` (`cxx:5639-5661`) | [`topology_stats::nb_sharp_edges`] |
//! | `ChFi3d_NumberOfEdges` (`cxx:5668-5685`) | [`topology_stats::number_of_edges`] |
//! | `ChFi3d_NumberOfSharpEdges` (`cxx:5691-5708`) | [`topology_stats::number_of_sharp_edges`] |
//! | `ChFi3d_Couture` (`cxx:5750-5769`) | [`seam::couture`] |
//! | `ChFi3d_CoutureOnVertex` (`cxx:5773-5801`) | [`seam::couture_on_vertex`] |
//! | `ChFi3d_IsPseudoSeam` (`cxx:5805-5832`) | [`seam::is_pseudo_seam`] |
//! | `ChFi3d_CircularSpine` (`cxx:846-880`) | [`spine_geom::circular_spine`] |
//! | `ChFi3d_Spine` (`cxx:888-903`) | [`spine_geom::spine`] |
//! | `ChFi3d_Parameters` (`cxx:4671-4703`) | [`surface_params::parameters`] |
//! | `ChFi3d_IsSmooth` (`cxx:5855-5942`) | [`smooth::is_smooth`] |
//! | `ChFi3d_CheckSameParameter` (`cxx:1565-1595`) | [`same_parameter::check_same_parameter`] |
//! | `ChFi3d_SameParameter` (`cxx:1602-1640`) | [`same_parameter::same_parameter`], [`same_parameter::same_parameter_range`] |
//! | `ChFi3d_TrsfTrans` (`cxx:2174-2186`) | [`fil_common_point::trsf_trans`] |
//! | `ChFi3d_FilCommonPoint` (`cxx:2193-2296`) | [`fil_common_point::fil_common_point`] |
//! | `CompParam` (`ChFi3d_Builder_6.cxx:209-311`) | [`comp_param::comp_param`] |
//!
//! UNPORTED in this batch, with their OCCT location:
//! - `ChFi3d_BoundSrfImpl` / `ChFi3d_BoundFac` / `ChFi3d_BoundSrf`
//!   (`cxx:158-240`, `:691-715`) — needs the `GeomAdaptor_Surface::Load`
//!   bound overload and `UResolution` / `VResolution`, which
//!   [`crate::topopebrep::tool::adaptor::SurfaceAdaptor`] does not carry.
//! - `ChFi3d_KParticular` (`cxx:604-683`) — needs
//!   `ChFiDS_Spine::CurrentElementarySpine(IE).GetType()` returning the
//!   *elementary* (per-interval) curve; the port's `CurveAdaptor` answers the
//!   whole edge curve.
//! - `ChFi3d_InterPlaneEdge` (`cxx:718-767`) — there is no ported
//!   `IntCurveSurface_HInter::Perform(curve, surface)` driver; only the
//!   per-arm pieces in [`crate::intcurvesurface`].
//! - `ChFi3d_ApproxByC2` (`cxx:5835-5851`) — needs
//!   `GeomAPI_PointsToBSpline(points, ChordLength, 3, 8, C2, 1.000001e-3)`,
//!   which the port has only as global cubic interpolation.

pub mod box2d;
pub mod comp_param;
pub mod faces;
pub mod fil_common_point;
pub mod period;
pub mod same_parameter;
pub mod seam;
pub mod smooth;
pub mod spine_geom;
pub mod surface_params;
pub mod topology_stats;
pub mod vertex_edges;

pub use box2d::{boite, boite_of_four};
pub use comp_param::comp_param;
pub use faces::{conex_faces, edge_state, eval_conti};
pub use fil_common_point::{fil_common_point, trsf_trans};
pub use period::{in_period, search_pivot};
pub use same_parameter::{check_same_parameter, same_parameter, same_parameter_range};
pub use seam::{couture, couture_on_vertex, is_pseudo_seam};
pub use smooth::is_smooth;
pub use spine_geom::{circular_spine, spine};
pub use surface_params::parameters;
pub use topology_stats::{
    cherche_bords_libres, nb_not_degenerated_edges, nb_sharp_edges, number_of_edges,
    number_of_sharp_edges,
};
pub use vertex_edges::{
    angle_edge, cherche_edge, cherche_element, cherche_face1, cherche_vertex, edge_common_faces,
    nbface,
};
