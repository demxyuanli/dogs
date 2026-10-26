//! `ChFi3d_SelectStripe`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:558-588`
//! (OCCT 8.0.0).
//!
//! ```text
//! bool ChFi3d_SelectStripe(List<handle(ChFiDS_Stripe)>::Iterator& It,
//!                          const TopoDS_Vertex& Vtx,
//!                          const bool thePrepareOnSame)
//! {
//!   if (!thePrepareOnSame) return true;
//!   for (; It.More(); It.Next()) {
//!     int sens = 0; stripe = It.Value();
//!     ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
//!     stat = (sens == 1) ? FirstStatus() : LastStatus();
//!     if (stat == ChFiDS_OnSame) return true;
//!   }
//!   return false;
//! }
//! ```
//!
//! The iterator is an **in/out** argument: on `true` it stays on the matching
//! stripe (the caller continues from there), on `false` it has run off the end.
//! [`StripeIterator`] carries that state; it is the port of the
//! `NCollection_List<handle(ChFiDS_Stripe)>::Iterator` the C1 bodies walk.

use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi_ds::{StripeHandle, State};
use crate::shape::TopoShape;

/// `NCollection_List<handle(ChFiDS_Stripe)>::Iterator` over the stripe list of
/// one vertex: `More` / `Next` / `Value` / `Initialize`, 0-based internally.
#[derive(Debug, Clone)]
pub struct StripeIterator<'a> {
    /// The list being walked (`myVDataMap(Index)`).
    stripes: &'a [StripeHandle],
    /// `myCurrent`: the index of `Value()`, `stripes.len()` once exhausted.
    index: usize,
}

impl<'a> StripeIterator<'a> {
    /// `Iterator(list)`: positioned on the first stripe.
    pub fn new(stripes: &'a [StripeHandle]) -> Self {
        Self { stripes, index: 0 }
    }

    /// `More()`.
    pub fn more(&self) -> bool {
        self.index < self.stripes.len()
    }

    /// `Next()`; stepping past the last element leaves the iterator exhausted.
    pub fn next(&mut self) {
        if self.index < self.stripes.len() {
            self.index += 1;
        }
    }

    /// `Value()` / `ChangeValue()`; `None` once exhausted.
    pub fn value(&self) -> Option<&StripeHandle> {
        self.stripes.get(self.index)
    }

    /// `Initialize()`: back to the first stripe.
    pub fn initialize(&mut self) {
        self.index = 0;
    }
}

/// `ChFi3d_SelectStripe(It, Vtx, thePrepareOnSame)` (`cxx:558-588`).
///
/// With `prepare_on_same == false` it answers `true` without touching the
/// iterator (the caller walks every stripe). Otherwise it advances to the first
/// stripe whose end status at `vtx` is `ChFiDS_OnSame` and answers whether one
/// was found.
///
/// `ChFi3d_IndexOfSurfData` raises `Standard_ConstructionError` in its last arm;
/// that becomes the `Err` here (OCCT lets it escape the corner body).
pub fn ch_fi3d_select_stripe(
    it: &mut StripeIterator<'_>,
    vtx: &TopoShape,
    prepare_on_same: bool,
) -> Result<bool, String> {
    if !prepare_on_same {
        return Ok(true);
    }

    while it.more() {
        let Some(stripe) = it.value() else {
            break;
        };
        let mut sens = 0;
        ch_fi3d_index_of_surf_data(vtx, stripe, &mut sens)?;

        let spine = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_SelectStripe: poisoned stripe".to_string())?;
            guard.spine().cloned()
        }
        .ok_or_else(|| "ChFi3d_SelectStripe: the stripe has no spine".to_string())?;

        let stat = {
            let guard = spine
                .read()
                .map_err(|_| "ChFi3d_SelectStripe: poisoned spine".to_string())?;
            if sens == 1 {
                guard.first_status()
            } else {
                guard.last_status()
            }
        };
        if stat == State::OnSame {
            return Ok(true);
        }
        it.next();
    }

    Ok(false)
}
