//! The `SurfData` / DS accessors shared by the two `PerformThreeCorner` bodies.
//!
//! `ChFi3d_FilBuilder_C3.cxx` and `ChFi3d_ChBuilder_C3.cxx` both hold
//! `ChFiDS_SurfData` by handle and write it in place, so two names can stand
//! for one slot (`fddeb` and `fdpiv` when the two stripes are one). The port
//! reaches the slot through the stripe that owns it:
//!
//! - [`sd`] copies one slot out (a read),
//! - [`with_sd`] writes one slot in place under `Arc::make_mut`, which keeps the
//!   aliasing the C++ handles give,
//! - [`sd_len`] is `SetOfSurfData()->Length()`,
//! - [`pcurve_on_surf`] / [`ds_curve`] / [`raise_curve_tolerance`] /
//!   [`data_surface`] are the DS lookups of the update blocks.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::curve::Curve;

use crate::chfi_ds::{StripeHandle, SurfData};
use crate::topopebrep::ds::DataStructure;

use super::ChFi3dBuilder;

/// `Corner->SetOfSurfData()->Value(index)` as a copy, the 1-based `index`.
pub(crate) fn sd(stripe: &StripeHandle, index: i32) -> Result<SurfData, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .and_then(|i| guard.set_of_surf_data().get(i))
        .map(|sd| sd.as_ref().clone())
        .ok_or_else(|| format!("PerformThreeCorner: SurfData {index} is out of range"))
}

/// `Corner->SetOfSurfData()->ChangeValue(index)` written in place.
pub(crate) fn with_sd<R>(
    stripe: &StripeHandle,
    index: i32,
    f: impl FnOnce(&mut SurfData) -> R,
) -> Result<R, String> {
    let mut guard = stripe
        .write()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let k = index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .ok_or_else(|| format!("PerformThreeCorner: SurfData {index} is out of range"))?;
    let slot = guard
        .change_set_of_surf_data()
        .get_mut(k)
        .ok_or_else(|| format!("PerformThreeCorner: SurfData {index} is out of range"))?;
    Ok(f(Arc::make_mut(slot)))
}

/// `Corner->SetOfSurfData()->Length()`.
pub(crate) fn sd_len(stripe: &StripeHandle) -> Result<i32, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    i32::try_from(guard.set_of_surf_data().len())
        .map_err(|_| "PerformThreeCorner: the SurfData sequence is too long".to_string())
}

/// `Fd->InterferenceOnS1()/InterferenceOnS2().PCurveOnSurf()->Value(par)`.
pub(crate) fn pcurve_on_surf(
    data: &SurfData,
    on_s: i32,
    param: f64,
) -> Result<GpPnt2d, String> {
    data.interference(on_s)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the interference has no p-curve".to_string())
        .map(|c| c.value(param))
}

/// `DStr.Curve(ind).Curve()`.
pub(crate) fn ds_curve(base: &ChFi3dBuilder, ind: i32) -> Result<Arc<dyn Curve>, String> {
    base.my_ds
        .ds()
        .curve(ind)
        .curve()
        .cloned()
        .ok_or_else(|| format!("PerformThreeCorner: DStr.Curve({ind}).Curve() is null"))
}

/// `tc.Tolerance(std::max(tol, tc.Tolerance()))`.
pub(crate) fn raise_curve_tolerance(base: &mut ChFi3dBuilder, ind: i32, tol: f64) {
    let curve = base.my_ds.change_ds().change_curve(ind);
    let raised = curve.tolerance().max(tol);
    curve.set_tolerance(raised);
}

/// `DStr.Surface(Fd->Surf()).Surface()` for one `SurfData`.
pub(crate) fn data_surface(
    ds: &DataStructure,
    data: &SurfData,
) -> Result<Arc<dyn occt_geom::Surface>, String> {
    ds.surface(data.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformThreeCorner: DStr.Surface(Surf()).Surface() is null".to_string())
}
