//! The splice that closes one turn of `PerformSetOfSurfOnElSpine`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:2870-2928`
//! (OCCT 8.0.0).
//!
//! Once the section computation set the builder's `done` flag, every
//! `SurfData` of `SeqSD` is pushed into the stripe right after (forward walk)
//! or right before (backward walk) the previous one, and the slices that fell
//! into an extension zone are cleaned by `ChFi3d_Purge`:
//!
//! - `!Ok1 && !obstacleon1` purges the contact on `S1` with
//!   `ref->Vertex(!forward, 1)`,
//! - `!Ok2 && !obstacleon2` does the same for `S2`,
//! - `ref = refbis` afterwards, i.e. the reference of the next section is the
//!   last inserted one.
//!
//! Divergences that are not heuristics:
//! - `InsertAfter` / `InsertBefore` take a shared handle: the port's sequence
//!   stores `Arc<SurfData>` clones,
//! - `ch_fi3d_purge` answers the handle as it now lives in the stripe, so the
//!   reference the next section starts from sees the purged data even when
//!   `Arc::make_mut` had to replace the sequence slot. When the purge removed
//!   the section, the caller's handle is kept, exactly like the OCCT handle the
//!   purge still owns (`cxx:2067-2077`).

use std::sync::Arc;

use crate::chfi_ds::SurfData;

use super::super::stripe_surgery::{ch_fi3d_purge, insert_after, insert_before};
use super::Walk;

impl<'a> Walk<'a> {
    /// `refbis = ref; ...; ref = refbis;` (`cxx:2872-2927`).
    pub(super) fn splice(&mut self, seq_sd: &mut Vec<Arc<SurfData>>) -> Result<(), String> {
        // `refbis = ref;` (`cxx:2873`).
        let mut refbis = self.ref_sd.clone();

        if self.forward {
            // `for (ii = 1; ii <= SeqSD.Length(); ii++) { SD = SeqSD(ii);
            // ...; InsertAfter(Stripe, refbis, SD); refbis = SD; }`
            // (`cxx:2875-2891`).
            for sd in seq_sd.iter_mut() {
                self.set_indices(Arc::make_mut(sd))?;
                insert_after(self.stripe, refbis.as_ref(), sd)?;
                refbis = Some((*sd).clone());
            }
        } else {
            // `for (ii = SeqSD.Length(); ii >= 1; ii--) { ...;
            // InsertBefore(Stripe, refbis, SD); refbis = SD; }`
            // (`cxx:2894-2910`).
            for i in (0..seq_sd.len()).rev() {
                self.set_indices(Arc::make_mut(&mut seq_sd[i]))?;
                insert_before(self.stripe, refbis.as_ref(), &seq_sd[i])?;
                refbis = Some(seq_sd[i].clone());
            }
        }

        // `if (!Ok1 && !obstacleon1) ChFi3d_Purge(Stripe, SD,
        // ref->Vertex(!forward, 1), !forward, 1, intf, intl);`
        // (`cxx:2912-2918`).
        if !self.ok1 && !self.obstacleon1 {
            refbis = self.purge_contact(refbis, 1)?;
        }
        // `if (!Ok2 && !obstacleon2) ChFi3d_Purge(Stripe, SD,
        // ref->Vertex(!forward, 2), !forward, 2, intf, intl);`
        // (`cxx:2920-2926`).
        if !self.ok2 && !self.obstacleon2 {
            refbis = self.purge_contact(refbis, 2)?;
        }

        // `ref = refbis;` (`cxx:2925`).
        self.ref_sd = refbis;
        Ok(())
    }

    /// One `ChFi3d_Purge` of [`Self::splice`], with the `VRef` taken from the
    /// reference `SurfData` **before** the splice (OCCT reassigns `ref` only at
    /// `cxx:2925`).
    fn purge_contact(
        &mut self,
        refbis: Option<Arc<SurfData>>,
        ons: i32,
    ) -> Result<Option<Arc<SurfData>>, String> {
        let v_ref = self
            .ref_sd
            .as_ref()
            .map(|sd| sd.vertex(!self.forward, ons).clone())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: no reference SurfData".to_string())?;
        let target = refbis
            .clone()
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: no inserted SurfData".to_string())?;
        let current = ch_fi3d_purge(
            self.stripe,
            &target,
            &v_ref,
            !self.forward,
            ons,
            &mut self.intf,
            &mut self.intl,
        )?;
        // A removed section keeps the caller's handle alive in OCCT.
        Ok(current.or(Some(target)))
    }
}
