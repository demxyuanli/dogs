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
//!                    _IndexSurfaceInDS / _ProjPC
//!                    (`ChFiKPart_ComputeData_Fcts.cxx`)
//! corner_spine       ChFiKPart_CornerSpine (`ChFiKPart_ComputeData_CS.cxx`)
//! fillet_plane_plane ChFiKPart_MakeFillet(Plane, Plane)
//!                    (`ChFiKPart_ComputeData_FilPlnPln.cxx`)
//! fillet_plane_cylinder
//!                    ChFiKPart_MakeFillet(Plane, Cylinder), line spine
//!                    (`ChFiKPart_ComputeData_FilPlnCyl.cxx:52-276`)
//! fillet_plane_cylinder_circle
//!                    ChFiKPart_MakeFillet(Plane, Cylinder), circle spine
//!                    (`ChFiKPart_ComputeData_FilPlnCyl.cxx:283-599`)
//! fillet_plane_cone  ChFiKPart_MakeFillet(Plane, Cone)
//!                    (`ChFiKPart_ComputeData_FilPlnCon.cxx:48-341`)
//! chamfer_plane_plane  ChFiKPart_MakeChamfer(Plane, Plane), Sym / TwoDist
//!                    (`ChFiKPart_ComputeData_ChPlnPln.cxx:53-272`)
//! chamfer_asym_plane_plane
//!                    ChFiKPart_MakeChAsym(Plane, Plane), DistAngle
//!                    (`ChFiKPart_ComputeData_ChAsymPlnPln.cxx:53-266`)
//! chamfer_plane_plane_tail
//!                    the `ChangeSurf`..`return true` range shared by the two
//!                    (`ChPlnPln.cxx:163-271` = `ChAsymPlnPln.cxx:157-265`)
//! chamfer_plane_cylinder_circle
//!                    ChFiKPart_MakeChamfer(Plane, Cylinder), circle spine
//!                    (`ChFiKPart_ComputeData_ChPlnCyl.cxx:64-363`)
//! chamfer_plane_cylinder_line
//!                    ChFiKPart_MakeChamfer(Plane, Cylinder), line spine
//!                    (`ChFiKPart_ComputeData_ChPlnCyl.cxx:372-592`)
//! chamfer_cylinder_cone_tail
//!                    the cone tail shared by the two above
//!                    (`ChPlnCyl.cxx:145-362` = `ChAsymPlnCyl.cxx:151-367`)
//! chamfer_cylinder_plane_tail
//!                    the plane tail shared by the line overloads
//!                    (`ChPlnCyl.cxx:471-591` = `ChAsymPlnCyl.cxx:536-653`)
//! chamfer_plane_cone  ChFiKPart_MakeChamfer(Plane, Cone), Sym / TwoDist
//!                    (`ChFiKPart_ComputeData_ChPlnCon.cxx:49-166`)
//! chamfer_asym_plane_cylinder_circle / _line
//!                    ChFiKPart_MakeChAsym(Plane, Cylinder), DistAngle
//!                    (`ChAsymPlnCyl.cxx:65-369` / `:376-655`)
//! chamfer_asym_plane_cone
//!                    ChFiKPart_MakeChAsym(Plane, Cone), DistAngle
//!                    (`ChFiKPart_ComputeData_ChAsymPlnCon.cxx:57-705`)
//! dispatcher         ChFiKPart_ComputeData::Compute
//!                    (`ChFiKPart_ComputeData.cxx`)
//! sphere_corner      ChFiKPart_Sphere
//!                    (`ChFiKPart_ComputeData_Sphere.cxx:47-206`)
//! rotule             ChFiKPart_MakeRotule
//!                    (`ChFiKPart_ComputeData_Rotule.cxx:43-173`)
//! compute_corner     the three ComputeCorner overloads
//!                    (`ChFiKPart_ComputeData.cxx:641-762`)
//! ```
//!
//! Not wired into [`crate::chfi3d::builder::perform_set_of_k_part`] yet: the
//! only OCCT consumer of `Compute`'s answer is the `else if (!SplitKPart(...))`
//! arm (`ChFi3d_Builder_2.cxx:3073-3076`), and `SplitKPart`
//! (`ChFi3d_Builder_SpKP.cxx:749-1160`) is UNPORTED — see
//! [`crate::chfi3d::builder::split_k_part`]. Answering `false` there keeps the
//! empty arm OCCT takes when no corner data was produced.

pub mod chamfer_asym_plane_cone;
pub mod chamfer_asym_plane_cylinder_circle;
pub mod chamfer_asym_plane_cylinder_line;
pub mod chamfer_asym_plane_plane;
pub mod chamfer_cylinder_cone_tail;
pub mod chamfer_cylinder_plane_tail;
pub mod chamfer_plane_cone;
pub mod chamfer_plane_cylinder_circle;
pub mod chamfer_plane_cylinder_line;
pub mod chamfer_plane_plane;
pub mod chamfer_plane_plane_tail;
pub mod compute_corner;
pub mod corner_spine;
pub mod dispatcher;
pub mod fcts;
pub mod fillet_plane_cone;
pub mod fillet_plane_cylinder;
pub mod fillet_plane_cylinder_circle;
pub mod fillet_plane_plane;
pub mod rotule;
pub mod sphere_corner;

pub use chamfer_asym_plane_cone::ch_fi_kpart_make_ch_asym_plane_cone;
pub use chamfer_asym_plane_plane::ch_fi_kpart_make_ch_asym_plane_plane;
pub use chamfer_asym_plane_cylinder_circle::ch_fi_kpart_make_ch_asym_plane_cylinder_circle;
pub use chamfer_asym_plane_cylinder_line::ch_fi_kpart_make_ch_asym_plane_cylinder_line;
pub use chamfer_plane_cone::ch_fi_kpart_make_chamfer_plane_cone;
pub use chamfer_plane_cylinder_circle::ch_fi_kpart_make_chamfer_plane_cylinder_circle;
pub use chamfer_plane_cylinder_line::ch_fi_kpart_make_chamfer_plane_cylinder_line;
pub use chamfer_plane_plane::ch_fi_kpart_make_chamfer_plane_plane;
pub use corner_spine::ch_fi_kpart_corner_spine;
pub use compute_corner::{
    ch_fi_kpart_compute_corner_plane_cylinder, ch_fi_kpart_compute_corner_rotule,
    ch_fi_kpart_compute_corner_sphere,
};
pub use dispatcher::ch_fi_kpart_compute_data;
pub use fcts::{
    ch_fi_kpart_in_period, ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds,
    ch_fi_kpart_p_curve, ch_fi_kpart_proj_pc,
};
pub use fillet_plane_cone::ch_fi_kpart_make_fillet_plane_cone;
pub use fillet_plane_cylinder::ch_fi_kpart_make_fillet_plane_cylinder_line;
pub use fillet_plane_cylinder_circle::ch_fi_kpart_make_fillet_plane_cylinder_circle;
pub use fillet_plane_plane::ch_fi_kpart_make_fillet_plane_plane;
