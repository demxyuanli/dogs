//! `SameSide` half of `ChFi3d_FilBuilder::PerformTwoCorner`
//! (`ChFi3d_FilBuilder_C2.cxx:528-859`, OCCT 8.0.0).
//!
//! The body builds the corner surface from the two stripe boundaries and the
//! common face, either by the `ToricRotule` direct construction
//! (`cxx:529-567`, body `cxx:87-122`) or by `GeomFill_ConstrainedFilling`
//! (`cxx:568-682`), then rewrites the common points and the DS
//! (`cxx:684-859`):
//!
//! - [`toric_rotule`] - `cxx:87-122`, the three-plane test,
//! - [`pivot_bound`] - `cxx:604-665`, the fourth boundary and the pivot p-curve
//!   of the `!sameparam` filling,
//! - [`ds_update`] - `cxx:690-859`, the two `ChFi3d_ComputeArete` /
//!   `ChFi3d_ComputePCurv` pairs and the stripe / `SurfData` writes.

mod ds_update;
mod pivot_bound;
mod toric_rotule;

use std::sync::Arc;

use occt_geom::geom_fill::ConstrainedFilling;
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::chfi3d::builder_0_helpers::{mkbound_pfac2d, mkbound_segment};
use crate::chfi3d::side::{concave_side, next_side};
use crate::chfi_ds::{Stripe, SurfData};
use crate::chfi3d::builder::fil_builder::perform_surf::fil_spine;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi_k_part::ch_fi_kpart_compute_corner_rotule;
use crate::shape::Edge;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::CornerCtx;

/// `DStr.Surface(Fd->Surf()).Surface()`.
fn data_surface(
    ds: &crate::topopebrep::ds::DataStructure,
    data: &SurfData,
) -> Result<Arc<dyn Surface>, String> {
    ds.surface(data.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: DStr.Surface(Surf()).Surface() is null".to_string())
}

/// `occ::down_cast<ChFiDS_FilSpine>(st1->Spine())->Radius()` (`cxx:554`).
fn spine_radius(stripe: &crate::chfi_ds::StripeHandle) -> Result<f64, String> {
    let spine_handle = stripe
        .read()
        .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?
        .spine()
        .cloned();
    let spine = fil_spine(spine_handle.as_ref())
        .ok_or_else(|| "PerformTwoCorner: the spine is not a ChFiDS_FilSpine".to_string())?;
    let radius = {
        let guard = spine
            .read()
            .map_err(|_| "PerformTwoCorner: poisoned spine".to_string())?;
        guard
            .fil_radius()
            .ok_or_else(|| "PerformTwoCorner: ChFiDS_FilSpine::Radius()".to_string())?
    };
    Ok(radius)
}

/// `SameSide` arm of `cxx:528-859`.
pub(super) fn same_side(
    base: &mut ChFi3dBuilder,
    _shape: BlendFuncSectionShape,
    ctx: &mut CornerCtx,
    corner: &mut Stripe,
    coin: &mut SurfData,
) -> Result<bool, String> {
    let h_fa_co = SurfaceAdaptor::of_face(&ctx.fa_co);
    let h_brs1 = SurfaceAdaptor::of_face(&ctx.ff1);
    let h_brs2 = SurfaceAdaptor::of_face(&ctx.ff2);

    let done;

    if toric_rotule::toric_rotule(&h_fa_co, &h_brs1, &h_brs2, &ctx.st1, &ctx.st2)? {
        // Direct construction (`cxx:530-567`).
        let mut ori = ctx.o_fa_co;
        let (ori_s, ori_sff1) = {
            let guard = ctx
                .st1
                .read()
                .map_err(|_| "PerformTwoCorner: poisoned stripe".to_string())?;
            (
                guard.orientation_on_surface(ctx.i_fa_co1),
                guard.orientation_on_surface(ctx.i_fa_arc1),
            )
        };
        let mut off1 = ctx.ff1.0.orientation();
        let mut bid = 1;
        bid = next_side(&mut ori, &mut off1, ori_s, ori_sff1, bid);
        let mut op1 = Orientation::Forward;
        let mut op2 = Orientation::Forward;
        if ctx.yapiv {
            bid = concave_side(
                &h_brs1,
                &h_brs2,
                &Edge(ctx.pivot.clone()),
                &mut op1,
                &mut op2,
            )?;
        }
        let _ = bid;
        // `op1 = TopAbs::Reverse(op1); op2 = TopAbs::Reverse(op2);`
        // (`cxx:546-547`).
        let op1 = op1.reversed();
        let op2 = op2.reversed();
        let radius = spine_radius(&ctx.st1)?;
        // `done = ChFiKPart_ComputeData::ComputeCorner(DStr, coin, HFaCo,
        // HBRS1, HBRS2, OFaCo, ori, op1, op2, radius)` (`cxx:555-563`).
        done = ch_fi_kpart_compute_corner_rotule(
            base.my_ds.change_ds(),
            coin,
            &h_fa_co,
            &h_brs1,
            &h_brs2,
            ctx.o_fa_co,
            ori,
            op1,
            op2,
            radius,
        )?;
    } else {
        // Construction by filling (`cxx:568-682`).
        let i_arc1 = ctx.sd1.interference(ctx.i_fa_arc1);
        let i_co1 = ctx.sd1.interference(ctx.i_fa_co1);
        let i_arc2 = ctx.sd2.interference(ctx.i_fa_arc2);
        let i_co2 = ctx.sd2.interference(ctx.i_fa_co2);

        // `uPCArc1 = sd1->Interference(IFaArc1).Parameter(isfirst1);
        // p2da1 = ...PCurveOnSurf()->Value(uPCArc1); p2df1 =
        // sd1->Interference(IFaCo1).PCurveOnSurf()->Value(uPCArc1);
        // sd1->Interference(IFaCo1).PCurveOnFace()->D1(uPCArc1, p2dfac1,
        // v2dfac1);` (`cxx:571-577`).
        let u_pc_arc1 = i_arc1.parameter(ctx.is_first1);
        let p2da1 = i_arc1
            .p_curve_on_surf()
            .ok_or_else(|| "PerformTwoCorner: sd1/IFaArc1 has no p-curve".to_string())?
            .value(u_pc_arc1);
        let p2df1 = i_co1
            .p_curve_on_surf()
            .ok_or_else(|| "PerformTwoCorner: sd1/IFaCo1 has no p-curve".to_string())?
            .value(u_pc_arc1);
        let (p2dfac1, v2dfac1) = i_co1
            .p_curve_on_face()
            .ok_or_else(|| "PerformTwoCorner: sd1/IFaCo1 has no p-curve on face".to_string())?
            .d1(u_pc_arc1);
        let u_pc_arc2 = i_arc2.parameter(ctx.is_first2);
        let p2da2 = i_arc2
            .p_curve_on_surf()
            .ok_or_else(|| "PerformTwoCorner: sd2/IFaArc2 has no p-curve".to_string())?
            .value(u_pc_arc2);
        let p2df2 = i_co2
            .p_curve_on_surf()
            .ok_or_else(|| "PerformTwoCorner: sd2/IFaCo2 has no p-curve".to_string())?
            .value(u_pc_arc2);
        let (p2dfac2, v2dfac2) = i_co2
            .p_curve_on_face()
            .ok_or_else(|| "PerformTwoCorner: sd2/IFaCo2 has no p-curve on face".to_string())?
            .d1(u_pc_arc2);

        // `B1 = ChFi3d_mkbound(surf1, p2df1, p2da1, tolapp3d, 2.e-4); B2 =
        // ChFi3d_mkbound(surf2, p2df2, p2da2, tolapp3d, 2.e-4);` (`cxx:585-586`).
        let surf1 = data_surface(base.my_ds.ds(), &ctx.sd1)?;
        let surf2 = data_surface(base.my_ds.ds(), &ctx.sd2)?;
        let b1 = mkbound_segment(&surf1, &p2df1, &p2da1, base.tolapp3d, 2.0e-4, false)?;
        let b2 = mkbound_segment(&surf2, &p2df2, &p2da2, base.tolapp3d, 2.0e-4, false)?;

        // `Bfac = ChFi3d_mkbound(HFaCo, PCurveOnFace, Sens1, p2dfac1,
        // v2dfac1, Sens2, p2dfac2, v2dfac2, tolapp3d, 2.e-4);`
        // (`cxx:588-598`).
        let (pcurve_on_face, bfac) = mkbound_pfac2d(
            &h_fa_co,
            ctx.sens1,
            &p2dfac1,
            &v2dfac1,
            ctx.sens2,
            &p2dfac2,
            &v2dfac2,
            base.tolapp3d,
            2.0e-4,
        )?;

        // `GeomFill_ConstrainedFilling fil(8, 20);`
        // (`cxx:599`).
        let mut fil = ConstrainedFilling::new(8, 20);
        let mut h_fa_piv = SurfaceAdaptor::default();
        let mut pcurve_on_piv = None;
        if ctx.same_param {
            // `fil.Init(Bfac, B2, B1, true);` (`cxx:600-602`).
            fil.init_3(bfac, b2, b1, true)?;
        } else {
            let piv = pivot_bound::pivot_bound(base, ctx, &b1)?;
            h_fa_piv = piv.surface;
            pcurve_on_piv = piv.p_curve;
            // `fil.Init(Bfac, B2, Bpiv, B1, true);` (`cxx:608`).
            fil.init_4(bfac, b2, piv.boundary, b1, true)?;
        }

        // `occ::handle<Geom_Surface> Surfcoin = fil.Surface(); done =
        // CompleteData(coin, Surfcoin, HFaCo, PCurveOnFace, HFaPiv,
        // PCurveOnPiv, OFaCo, true, false, false, false, false);`
        // (`cxx:666-678`).
        let surfcoin = fil
            .surface()
            .cloned()
            .ok_or_else(|| "PerformTwoCorner: GeomFill_ConstrainedFilling::Surface()".to_string())?;
        let surfcoin: Arc<dyn Surface> = Arc::new(surfcoin);
        done = base.complete_data_geom(
            coin,
            &surfcoin,
            &h_fa_co,
            Some(&pcurve_on_face),
            &h_fa_piv,
            pcurve_on_piv.as_ref(),
            ctx.o_fa_co,
            true,
            false,
            false,
            false,
            false,
        )?;
    }

    if done {
        // `cxx:690-859`.
        ds_update::update_ds(base, ctx, corner, coin)?;
        ctx.flush()?;
    }

    Ok(done)
}
