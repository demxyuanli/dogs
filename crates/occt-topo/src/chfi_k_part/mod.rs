//! `ChFiKPart` — the analytic corner surfaces of the fillet and chamfer
//! builder (OCCT `ModelingAlgorithms/TKFillet/ChFiKPart`, 14 `.cxx` files,
//! OCCT 8.0.0).
//!
//! `ChFi3d_Builder::PerformSetOfKPart` asks
//! [`ch_fi_kpart_compute_data`] for one `SurfData` per edge that
//! [`crate::chfi3d::builder_0_helpers::ch_fi3d_k_particular`] accepted
//! (`ChFi3d_Builder_2.cxx:3070-3072`). The function is the port of
//! `ChFiKPart_ComputeData::Compute`, which dispatches to one maker per surface
//! pair (`ChFiKPart_ComputeData.cxx:91-636`).
//!
//! Module layout:
//! ```text
//! fcts               ChFiKPart_InPeriod / _PCurve / _IndexCurveInDS /
//!                    _IndexSurfaceInDS (`ChFiKPart_ComputeData_Fcts.cxx`)
//! corner_spine       ChFiKPart_CornerSpine (`ChFiKPart_ComputeData_CS.cxx`)
//! fillet_plane_plane ChFiKPart_MakeFillet(Plane, Plane)
//!                    (`ChFiKPart_ComputeData_FilPlnPln.cxx`)
//! dispatcher         ChFiKPart_ComputeData::Compute
//!                    (`ChFiKPart_ComputeData.cxx`)
//! ```
//!
//! Not wired into [`crate::chfi3d::builder::perform_set_of_k_part`] yet: the
//! only OCCT consumer of `Compute`'s answer is the `else if (!SplitKPart(...))`
//! arm (`ChFi3d_Builder_2.cxx:3073-3076`), and `SplitKPart`
//! (`ChFi3d_Builder_SpKP.cxx:749-1160`) is UNPORTED — see
//! [`crate::chfi3d::builder::split_k_part`]. Answering `false` there keeps the
//! empty arm OCCT takes when no corner data was produced.

pub mod corner_spine;
pub mod dispatcher;
pub mod fcts;
pub mod fillet_plane_plane;

pub use corner_spine::ch_fi_kpart_corner_spine;
pub use dispatcher::ch_fi_kpart_compute_data;
pub use fcts::{
    ch_fi_kpart_in_period, ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds,
    ch_fi_kpart_p_curve,
};
pub use fillet_plane_plane::ch_fi_kpart_make_fillet_plane_plane;
