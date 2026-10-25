//! `ChFi3d_Builder::SplitKPart` — **UNPORTED**.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx:749-1160`
//! (OCCT 8.0.0).
//!
//! `PerformSetOfKPart` hands every `SurfData` built by `ChFiKPart_ComputeData`
//! to this member (`ChFi3d_Builder_2.cxx:3073-3077`); it hatches the two
//! `PCurveOnFace` of the corner against the tangency lines of the two
//! `TopolTool`s (`Geom2dHatch_Hatcher` / `Geom2dHatch_Intersector`), splits the
//! `SurfData` per hatch domain with the file statics `Tri` / `FillSD` / `CpSD`
//! and updates `intf` / `intl` from `SearchFace`.
//!
//! Blockers that keep it out of this batch (each is a package of its own):
//! - `Geom2dHatch_Hatcher` / `Geom2dHatch_Intersector` / `HatchGen_Domain`
//!   (`ModelingAlgorithms/TKTopAlgo/Geom2dHatch`),
//! - the three `ChFi3d_Builder_SpKP.cxx` statics `Tri` (`:594-...`), `FillSD`
//!   and `CpSD`.
//!
//! The `Err` is only reachable when a `SurfData` was actually split, i.e. once
//! [`crate::chfi3d::chfi_k_part`] stops answering "no data".

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

/// `SplitKPart(Data, SetData, Spine, Iedge, S1, I1, S2, I2, intf, intl)`
/// (`ChFi3d_Builder_SpKP.cxx:749-760`) — **UNPORTED**, answers `Err`.
#[allow(clippy::too_many_arguments)]
pub fn split_k_part(
    _data: &Arc<SurfData>,
    _set_data: &mut Vec<Arc<SurfData>>,
    _spine: &SpineHandle,
    _iedge: i32,
    _s1: &SurfaceAdaptor,
    _i1: &Rc<RefCell<TopolTool>>,
    _s2: &SurfaceAdaptor,
    _i2: &Rc<RefCell<TopolTool>>,
    _intf: &mut bool,
    _intl: &mut bool,
) -> Result<bool, String> {
    Err("ChFi3d_Builder::SplitKPart: UNPORTED (ChFi3d_Builder_SpKP.cxx:749-1160)".into())
}
