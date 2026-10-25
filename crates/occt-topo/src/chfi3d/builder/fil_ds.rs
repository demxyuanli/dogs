//! `ChFi3d_FilDS` - fill the DS with one stripe's fillet surface data.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2565-3351`
//! (OCCT 8.0.0). `ChFi3d_Builder::Compute` calls it once per stripe after the
//! `ChFi3d_StripeEdgeInter` sweep (`ChFi3d_Builder.cxx:391-395`) and it throws
//! straight out of `Compute` when it fails (there is no `try` around that
//! call).
//!
//! UNPORTED. The body is 787 `.cxx` lines and needs the whole fillet-side DS
//! vocabulary, none of which the port has yet:
//!
//! - `ChFi3d_FilPointInDS` (`cxx:2341-2369`),
//! - `ChFi3d_FilVertexInDS` (`cxx:2371-2394`),
//! - `ChFi3d_FilCurveInDS` (`cxx:2156-...`),
//! - `ChFi3d_Contains` (`cxx:2445-2480`),
//! - `CutEdge` (`cxx:2488-2517`),
//! - `findIndexPoint` (`cxx:2523-2558`),
//! - both `ChFi3d_ComputePCurv` overloads (`cxx:1648-1700`, called at
//!   `cxx:2926` and `:2938`),
//! - `TopOpeBRepDS_SolidSurfaceInterference` / `SurfaceCurveInterference` /
//!   `CurvePointInterference` construction and `ChangeCurve().Nullify()` /
//!   `SetSCI` bookkeeping.
//!
//! The function is therefore an explicit hook instead of a partial body: a
//! silent no-op would leave `Compute` reporting success with an empty DS, and
//! an approximated body would be a heuristic OCCT does not have.

use crate::chfi_ds::StripeHandle;

use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `ChFi3d_FilDS(SolidIndex, CorDat, DStr, reglist, tol3d, tol2d)`
    /// (`ChFi3d_Builder_0.cxx:2565-3351`).
    pub fn fil_ds(
        &mut self,
        _solid_index: i32,
        _cor_dat: &StripeHandle,
        _tol3d: f64,
        _tol2d: f64,
    ) -> Result<(), String> {
        Err("ChFi3d_FilDS: UNPORTED (ChFi3d_Builder_0.cxx:2565-3351)".to_string())
    }
}
