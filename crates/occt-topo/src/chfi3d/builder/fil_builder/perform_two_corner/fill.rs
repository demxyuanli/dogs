//! The second half of `ChFi3d_FilBuilder::PerformTwoCorner`: the corner
//! `ChFiDS_SurfData` that the two stripes do not already share.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder_C2.cxx:516-1152`
//! (OCCT 8.0.0):
//!
//! - `cxx:516-529` - the new `corner` stripe, its single `coin` `SurfData`
//!   and the `SameSide` test: this file,
//! - `cxx:528-859` - `SameSide`: the `ToricRotule` direct construction
//!   (`cxx:529-567`, body `cxx:87-122`) or the `GeomFill_ConstrainedFilling`
//!   construction (`cxx:568-682`), then the `resetcp1` / `resetcp2` common
//!   points and the DS update of `cxx:684-859`: [`same_side`],
//! - `cxx:860-1144` - `OnSame` / `OnDiff`: the filling between the two stripe
//!   boundaries and the common face, then the DS update: [`on_diff`],
//! - `cxx:1145-1152` - `myEVIMap` / `myListStripe`: this file.
//!
//! `CornerCtx` carries the locals of the prologue that the two halves share;
//! OCCT keeps them as locals of one function body.

mod on_diff;
mod same_side;

use std::sync::{Arc, RwLock};

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::chfi_ds::{CommonPoint, State, Stripe, StripeHandle, SurfData};
use crate::shape::{Face, TopoShape};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::super::ChFi3dBuilder;
use super::index_slot;

/// The `PerformTwoCorner` locals of `cxx:150-405` that `cxx:516-1152` reads.
pub(super) struct CornerCtx {
    /// `st1`
    pub st1: StripeHandle,
    /// `st2`
    pub st2: StripeHandle,
    /// `sd1` - a copy, see `super`'s note on [`CornerCtx::flush`]
    pub sd1: SurfData,
    /// `sd2`
    pub sd2: SurfData,
    /// `Isd1`
    pub isd1: i32,
    /// `Isd2`
    pub isd2: i32,
    /// `isfirst1`
    pub is_first1: bool,
    /// `isfirst2`
    pub is_first2: bool,
    /// `IFaCo1`
    pub i_fa_co1: i32,
    /// `IFaCo2`
    pub i_fa_co2: i32,
    /// `IFaArc1`
    pub i_fa_arc1: i32,
    /// `IFaArc2`
    pub i_fa_arc2: i32,
    /// `OFaCo`
    pub o_fa_co: crate::abs::Orientation,
    /// `FaCo`
    pub fa_co: Face,
    /// `FF1`
    pub ff1: Face,
    /// `FF2`
    pub ff2: Face,
    /// `pivot`
    pub pivot: TopoShape,
    /// `Hpivot` - null exactly when `yapiv` is false
    pub hpivot: Option<CurveAdaptor>,
    /// `yapiv`
    pub yapiv: bool,
    /// `sameparam`
    pub same_param: bool,
    /// `parCP1`
    pub par_cp1: f64,
    /// `parCP2`
    pub par_cp2: f64,
    /// `OkinterCC`
    pub okinter_cc: bool,
    /// `Stat1`
    pub stat1: State,
    /// `Stat2`
    pub stat2: State,
    /// `SameSide` (`cxx:250`) - which half of `cxx:516-1152` runs
    pub same_side: bool,
    /// `UIntPC1`
    pub u_int_pc1: f64,
    /// `UIntPC2`
    pub u_int_pc2: f64,
    /// `CP1` - the `sd1` common point on the pivot arc
    pub cp1: CommonPoint,
    /// `CP2`
    pub cp2: CommonPoint,
    /// `resetcp1` (`cxx:375`)
    pub reset_cp1: bool,
    /// `resetcp2` (`cxx:376`)
    pub reset_cp2: bool,
    /// `PCurveOnPiv` (`cxx:418`) is a prologue local in OCCT but is only
    /// written and read inside the `!sameparam` arm, which
    /// [`same_side::pivot_bound`] keeps as a local value.
    /// `pctrans` (`cxx:417`)
    pub pctrans: crate::abs::Orientation,
    /// `Sens1` (`cxx:245`), the `ChFi3d_mkbound` direction of `Bfac`
    /// `FaPiv` (`cxx:416`)
    pub fa_piv: Option<Face>,
    /// `Sens1` (`cxx:245`), the `ChFi3d_mkbound` direction of `Bfac`
    pub sens1: i32,
    /// `Sens2`
    pub sens2: i32,
    /// `Vtx`, the `myVDataMap` key of `cxx:151`
    pub vtx: TopoShape,
}

impl CornerCtx {
    /// Writes the two `SurfData` copies back into their stripes.
    ///
    /// `PerformTwoCornerbyInter` and the DS update of `cxx:684-859` reach the
    /// stripes both directly (through `st1` / `st2`) and through `sd1` / `sd2`;
    /// writing one slot of one sequence in one pass keeps the two writes
    /// coherent when the two stripes are the same handle (`st2 == st1`).
    pub(super) fn flush(&mut self) -> Result<(), String> {
        if Arc::ptr_eq(&self.st1, &self.st2) {
            let mut stripe = self
                .st1
                .write()
                .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
            let (isd1, isd2) = (self.isd1, self.isd2);
            let (sd1, sd2) = (self.sd1.clone(), self.sd2.clone());
            let seq = stripe.change_set_of_surf_data();
            if let Some(slot) = index_slot(seq, isd1) {
                *slot = Arc::new(sd1);
            }
            if let Some(slot) = index_slot(seq, isd2) {
                *slot = Arc::new(sd2);
            }
            return Ok(());
        }
        for (stripe, index, data) in [
            (&self.st1, self.isd1, &self.sd1),
            (&self.st2, self.isd2, &self.sd2),
        ] {
            let mut guard = stripe
                .write()
                .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
            if let Some(slot) = index_slot(guard.change_set_of_surf_data(), index) {
                *slot = Arc::new(data.clone());
            }
        }
        Ok(())
    }
}

/// `if (!c1biseau && !done)` (`cxx:516-1152`).
///
/// The `corner` stripe and its single `coin` `SurfData` are created here
/// (`cxx:520-527`); the two halves fill the surface and the DS. When either
/// half succeeded, `cxx:1146-1151` appends `coin->Surf()` to `myEVIMap(Vtx)`
/// and `corner` to `myListStripe`.
pub(super) fn two_corner_fill(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    ctx: &mut CornerCtx,
) -> Result<(), String> {
    // `occ::handle<ChFiDS_Stripe> corner = new ChFiDS_Stripe(); cornerset =
    // corner->ChangeSetOfSurfData(); cornerset = new
    // NCollection_HSequence<...>(); coin = new ChFiDS_SurfData();
    // cornerset->Append(coin);` (`cxx:520-526`).
    let corner_handle: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
    let mut coin = SurfData::new();

    let done = {
        let mut corner = corner_handle
            .write()
            .map_err(|_| "PerformTwoCorner: poisoned corner stripe".to_string())?;
        corner.change_set_of_surf_data().push(Arc::new(coin.clone()));

        let done = if ctx.same_side {
            same_side::same_side(base, shape, ctx, &mut corner, &mut coin)?
        } else {
            on_diff::on_diff(base, shape, ctx, &mut corner, &mut coin)?
        };

        // OCCT mutates the `SurfData` the sequence holds; the port works on a
        // copy, so the sequence catches up here.
        if let Some(slot) = corner.change_set_of_surf_data().first_mut() {
            *slot = Arc::new(coin.clone());
        }
        done
    };

    // `if (!myEVIMap.IsBound(Vtx)) { NCollection_List<int> li;
    // myEVIMap.Bind(Vtx, li); } myEVIMap.ChangeFind(Vtx).Append(coin->Surf());
    // myListStripe.Append(corner);` (`cxx:1146-1151`).
    if done {
        let surf = coin.surf();
        base.change_evi(&ctx.vtx).push(surf);
        base.my_list_stripe.push(corner_handle);
    }
    Ok(())
}
