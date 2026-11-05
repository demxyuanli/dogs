//! `OnSame` / `OnDiff` half of `ChFi3d_FilBuilder::PerformTwoCorner`
//! (`ChFi3d_FilBuilder_C2.cxx:860-1144`, OCCT 8.0.0).
//!
//! The body sorts the two stripes into `stsam` / `stdif` by their spine
//! status (`cxx:872-907`), checks that the pivot arc of the `OnDiff` stripe is
//! a loop of the `OnSame` face (`cxx:908-927`), fills the corner between the
//! two stripe boundaries and the bridge on that face (`cxx:930-1004`) and
//! rewrites the DS (`cxx:1006-1144`, [`ds_write`]).

mod ds_write;

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::geom_fill::ConstrainedFilling;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::chfi_ds::{CommonPoint, State, Stripe, StripeHandle, SurfData};
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi3d::builder_0_helpers::{
    build_p_curve_vec3d, mkbound_curve, mkbound_curve_segment,
};
use crate::shape::{Edge, Face};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

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

/// `sdsam->Interference(ifacosam).PCurveOnSurf()` and friends.
fn p_curve_on_surf(
    data: &SurfData,
    index: i32,
    u: f64,
    what: &str,
) -> Result<GpPnt2d, String> {
    Ok(data
        .interference(index)
        .p_curve_on_surf()
        .ok_or_else(|| format!("PerformTwoCorner: {what} has no p-curve"))?
        .value(u))
}

/// The two stripes sorted by spine status (`cxx:868-907`).
///
/// `sdsam` / `sddif` are working copies of `ctx.sd1` / `ctx.sd2`: OCCT holds
/// handles into the stripes' `SurfData` sequences, and the port flushes the
/// copies back through [`CornerCtx::flush`] once [`ds_write::write_ds`] is
/// done.
pub(super) struct SortedCtx {
    /// `stsam`
    pub st_sam: StripeHandle,
    /// `stdif`
    pub st_dif: StripeHandle,
    /// `sdsam`
    pub sd_sam: SurfData,
    /// `sddif`
    pub sd_dif: SurfData,
    /// `uintpcsam` - `UIntPC1` or `UIntPC2`
    pub u_int_pc_sam: f64,
    /// `uintpcdif`
    pub u_int_pc_dif: f64,
    /// `ifacosam`
    pub i_fa_co_sam: i32,
    /// `ifaopsam`
    pub i_fa_op_sam: i32,
    /// `ifacodif`
    pub i_fa_co_dif: i32,
    /// `ifaopdif`
    pub i_fa_op_dif: i32,
    /// `isfirstsam`
    pub is_first_sam: bool,
    /// `isfirstdif`
    pub is_first_dif: bool,
    /// `cpopdif` - `sddif->ChangeVertex(isfirstdif, ifaopdif)`
    pub cp_op_dif: CommonPoint,
    /// `pcFopsam`
    pub pc_fopsam: Option<Arc<dyn Curve2d>>,
    /// `pcsurfdif`
    pub pc_surf_dif: Option<Arc<dyn Curve2d>>,
    /// `HBRFopsam` - `BRepAdaptor_Surface(Fopsam)` (`cxx:975-977`)
    pub hbr_fopsam: SurfaceAdaptor,
}

/// `else` arm of `cxx:860-1144`.
pub(super) fn on_diff(
    base: &mut ChFi3dBuilder,
    _shape: BlendFuncSectionShape,
    ctx: &mut CornerCtx,
    corner: &mut Stripe,
    coin: &mut SurfData,
) -> Result<bool, String> {
    // `if (!OkinterCC) throw Standard_Failure("TwoCorner : No intersection pc
    // pc");` (`cxx:862-866`).
    if !ctx.okinter_cc {
        return Err("TwoCorner : No intersection pc pc".to_string());
    }

    // `if (Stat1 == OnSame && Stat2 == OnDiff) { ... } else if (Stat1 ==
    // OnDiff && Stat2 == OnSame) { ... } else throw
    // Standard_Failure("TwoCorner : Config unknown");` (`cxx:872-907`).
    let sam_is_1 = match (ctx.stat1, ctx.stat2) {
        (State::OnSame, State::OnDiff) => true,
        (State::OnDiff, State::OnSame) => false,
        _ => return Err("TwoCorner : Config unknown".to_string()),
    };

    let mut sorted = if sam_is_1 {
        SortedCtx {
            st_sam: Arc::clone(&ctx.st1),
            st_dif: Arc::clone(&ctx.st2),
            sd_sam: ctx.sd1.clone(),
            sd_dif: ctx.sd2.clone(),
            u_int_pc_sam: ctx.u_int_pc1,
            u_int_pc_dif: ctx.u_int_pc2,
            i_fa_co_sam: ctx.i_fa_co1,
            i_fa_op_sam: ctx.i_fa_arc1,
            i_fa_co_dif: ctx.i_fa_co2,
            i_fa_op_dif: ctx.i_fa_arc2,
            is_first_sam: ctx.is_first1,
            is_first_dif: ctx.is_first2,
            cp_op_dif: CommonPoint::new(),
            pc_fopsam: None,
            hbr_fopsam: SurfaceAdaptor::default(),
            pc_surf_dif: None,
        }
    } else {
        SortedCtx {
            st_sam: Arc::clone(&ctx.st2),
            st_dif: Arc::clone(&ctx.st1),
            sd_sam: ctx.sd2.clone(),
            sd_dif: ctx.sd1.clone(),
            u_int_pc_sam: ctx.u_int_pc2,
            u_int_pc_dif: ctx.u_int_pc1,
            i_fa_co_sam: ctx.i_fa_co2,
            i_fa_op_sam: ctx.i_fa_arc2,
            i_fa_co_dif: ctx.i_fa_co1,
            i_fa_op_dif: ctx.i_fa_arc1,
            is_first_sam: ctx.is_first2,
            is_first_dif: ctx.is_first1,
            cp_op_dif: CommonPoint::new(),
            pc_fopsam: None,
            hbr_fopsam: SurfaceAdaptor::default(),
            pc_surf_dif: None,
        }
    };

    // `ChFiDS_CommonPoint& cpopdif = sddif->ChangeVertex(isfirstdif,
    // ifaopdif); if (!cpopdif.IsOnArc()) throw Standard_Failure("TwoCorner :
    // No point on restriction on surface OnDiff");` (`cxx:910-916`).
    let cp_op_dif = sorted
        .sd_dif
        .vertex(sorted.is_first_dif, sorted.i_fa_op_dif)
        .clone();
    if !cp_op_dif.is_on_arc() {
        return Err("TwoCorner : No point on restriction on surface OnDiff".to_string());
    }
    sorted.cp_op_dif = cp_op_dif.clone();
    let arc_op_dif = Edge(
        cp_op_dif
            .arc()
            .ok_or_else(|| "TwoCorner : No point on restriction on surface OnDiff".to_string())?
            .clone(),
    );

    // `const TopoDS_Face& Fopsam =
    // TopoDS::Face(DStr.Shape(sdsam->Index(ifaopsam)));` (`cxx:917`).
    let i_fa_op_sam = sorted.sd_sam.index(sorted.i_fa_op_sam);
    let fopsam = Face(
        base.my_ds
            .ds()
            .shape(i_fa_op_sam, false)
            .ok_or_else(|| "TwoCorner : the OnSame arc face is not in the DS".to_string())?,
    );

    // `for (ex.Init(Fopsam, TopAbs_EDGE); ex.More(); ex.Next()) { if
    // (ex.Current().IsSame(Arcopdif)) break; else if (!ex.More()) throw
    // Standard_Failure("TwoCorner : No common face to loop the contour"); }`
    // (`cxx:918-927`). The `else if (!ex.More())` can never fire: `More()` is
    // the loop condition itself, so the loop only searches.
    for e in explorer_shapes(&fopsam.0, ShapeType::Edge) {
        if is_same(&e, &arc_op_dif.0) {
            break;
        }
    }

    // `gp_Pnt2d ppopsam = sdsam->Interference(ifaopsam).PCurveOnSurf()
    // ->Value(uintpcsam); gp_Pnt2d ppcosam =
    // sdsam->Interference(ifacosam).PCurveOnSurf()->Value(uintpcsam);`
    // (`cxx:932-933`).
    let pp_op_sam = p_curve_on_surf(
        &sorted.sd_sam,
        sorted.i_fa_op_sam,
        sorted.u_int_pc_sam,
        "sdsam/ifaopsam",
    )?;
    let pp_co_sam = p_curve_on_surf(
        &sorted.sd_sam,
        sorted.i_fa_co_sam,
        sorted.u_int_pc_sam,
        "sdsam/ifacosam",
    )?;

    // `occ::handle<Geom_Surface> surfsam =
    // DStr.Surface(sdsam->Surf()).Surface(); occ::handle<GeomAdaptor_Surface>
    // Hsurfsam = new GeomAdaptor_Surface(surfsam); occ::handle<Geom2d_Curve>
    // pcsurfsam; Bsam = ChFi3d_mkbound(Hsurfsam, pcsurfsam, ppopsam, ppcosam,
    // tolapp3d, 2.e-4);` (`cxx:934-937`).
    let surf_sam = data_surface(base.my_ds.ds(), &sorted.sd_sam)?;
    let h_surf_sam = SurfaceAdaptor::of_surface(Arc::clone(&surf_sam));
    let (pc_surf_sam, b_sam) = mkbound_curve_segment(
        &h_surf_sam,
        &pp_op_sam,
        &pp_co_sam,
        base.tolapp3d,
        2.0e-4,
        false,
    )?;

    // `double upcopdif = sddif->Interference(ifaopdif).Parameter(isfirstdif);
    // gp_Pnt2d ppopdif = sddif->Interference(ifaopdif).PCurveOnSurf()
    // ->Value(upcopdif); gp_Pnt2d ppcodif =
    // sddif->Interference(ifacodif).PCurveOnSurf()->Value(uintpcdif);`
    // (`cxx:939-941`).
    let up_co_p_dif = sorted
        .sd_dif
        .interference(sorted.i_fa_op_dif)
        .parameter(sorted.is_first_dif);
    let pp_op_dif = p_curve_on_surf(
        &sorted.sd_dif,
        sorted.i_fa_op_dif,
        up_co_p_dif,
        "sddif/ifaopdif",
    )?;
    let pp_co_dif = p_curve_on_surf(
        &sorted.sd_dif,
        sorted.i_fa_co_dif,
        sorted.u_int_pc_dif,
        "sddif/ifacodif",
    )?;

    // `surfdif = DStr.Surface(sddif->Surf()).Surface(); Hsurfdif = new
    // GeomAdaptor_Surface(surfdif); pcsurfdif; Bdif =
    // ChFi3d_mkbound(Hsurfdif, pcsurfdif, ppcodif, ppopdif, tolapp3d,
    // 2.e-4);` (`cxx:942-946`).
    let surf_dif = data_surface(base.my_ds.ds(), &sorted.sd_dif)?;
    let h_surf_dif = SurfaceAdaptor::of_surface(Arc::clone(&surf_dif));
    let (pc_surf_dif, b_dif) = mkbound_curve_segment(
        &h_surf_dif,
        &pp_co_dif,
        &pp_op_dif,
        base.tolapp3d,
        2.0e-4,
        false,
    )?;
    sorted.pc_surf_dif = Some(pc_surf_dif);

    // `sdsam->Interference(ifaopsam).PCurveOnFace()->D0(uintpcsam,
    // ppfacsam); const occ::handle<Geom_Curve>& curvopsam =
    // DStr.Curve(sdsam->Interference(ifaopsam).LineIndex()).Curve();
    // curvopsam->D1(uintpcsam, PPfacsam, VVfacsam);` (`cxx:949-952`).
    let pp_fac_sam = sorted
        .sd_sam
        .interference(sorted.i_fa_op_sam)
        .p_curve_on_face()
        .ok_or_else(|| "PerformTwoCorner: sdsam/ifaopsam has no p-curve on face".to_string())?
        .d0(sorted.u_int_pc_sam);
    let line_index_sam = sorted.sd_sam.interference(sorted.i_fa_op_sam).line_index();
    let cur_op_sam = base
        .my_ds
        .ds()
        .curve(line_index_sam)
        .curve()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: DStr.Curve(I).Curve() is null".to_string())?;
    let (_, vv_fac_sam) = cur_op_sam.d1(sorted.u_int_pc_sam);

    // `BRepAdaptor_Curve2d PCArcFac(Arcopdif, Fopsam);
    // PCArcFac.D0(cpopdif.ParameterOnArc(), ppfacdif);` (`cxx:953-954`).
    let par_on_arc = cp_op_dif
        .parameter_on_arc()
        .ok_or_else(|| "TwoCorner : No point on restriction on surface OnDiff".to_string())?;
    let pc_arc_fac = Curve2dAdaptor::of_edge_face(&arc_op_dif, &fopsam);
    let mut pp_fac_dif = pc_arc_fac
        .d0(par_on_arc)
        .ok_or_else(|| "PerformTwoCorner: null BRepAdaptor_Curve2d curve".to_string())?;

    // `BRepAdaptor_Surface SurFopsam(Fopsam); if (SurFopsam.IsUClosed()) {
    // double Uperiod = SurFopsam.LastUParameter() -
    // SurFopsam.FirstUParameter(); if (std::abs(ppfacsam.X() - ppfacdif.X()) >
    // Uperiod / 2) { ... } }` (`cxx:956-971`).
    let sur_fopsam = SurfaceAdaptor::of_face(&fopsam);
    let u_closed = sur_fopsam
        .surface()
        .map(|s| s.is_u_closed())
        .unwrap_or(false);
    if u_closed {
        let u_period = sur_fopsam.last_u_parameter() - sur_fopsam.first_u_parameter();
        if (pp_fac_sam.x() - pp_fac_dif.x()).abs() > u_period / 2.0 {
            if pp_fac_dif.x() < pp_fac_sam.x() {
                pp_fac_dif.set_x(pp_fac_dif.x() + u_period);
            } else {
                pp_fac_dif.set_x(pp_fac_dif.x() - u_period);
            }
        }
    }

    // `BRepAdaptor_Curve CArcFac(Arcopdif);
    // CArcFac.D1(cpopdif.ParameterOnArc(), PPfacdif, VVfacdif);` (`cxx:973-975`).
    let (_, vv_fac_dif) = CurveAdaptor::of_edge(&arc_op_dif)
        .d1(par_on_arc)
        .ok_or_else(|| "TwoCorner : the OnDiff arc has no curve".to_string())?;

    // `occ::handle<BRepAdaptor_Surface> HBRFopsam = new
    // BRepAdaptor_Surface(); BRFopsam.Initialize(Fopsam, false); occ::handle
    // <Geom2d_Curve> pcFopsam = ChFi3d_BuildPCurve(HBRFopsam, ppfacsam,
    // VVfacsam, ppfacdif, VVfacdif, true);` (`cxx:976-979`).
    //
    // `Initialize(F, false)` loads the surface with its own bounds, which is
    // the `BRepAdaptor_Surface(F, false)` constructor [`SurfaceAdaptor::of_face`]
    // mirrors.
    let hbr_fopsam = SurfaceAdaptor::of_face(&fopsam);
    sorted.pc_fopsam = Some(Arc::new(build_p_curve_vec3d(
        &hbr_fopsam,
        &pp_fac_sam,
        &vv_fac_sam,
        &pp_fac_dif,
        &vv_fac_dif,
        true,
    )?));
    sorted.hbr_fopsam = hbr_fopsam;

    // `Bfac = ChFi3d_mkbound(HBRFopsam, pcFopsam, tolapp3d, 2.e-4);`
    // (`cxx:980`).
    let bfac_surface = sorted
        .hbr_fopsam
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_mkbound: null Adaptor3d_Surface".to_string())?;
    let b_fac = mkbound_curve(
        &bfac_surface,
        sorted
            .pc_fopsam
            .as_ref()
            .ok_or_else(|| "PerformTwoCorner: null pcFopsam".to_string())?,
        base.tolapp3d,
        2.0e-4,
        false,
    );

    // `GeomFill_ConstrainedFilling fil(8, 20); fil.Init(Bsam, Bdif, Bfac,
    // true);` (`cxx:981-982`).
    let mut fil = ConstrainedFilling::new(8, 20);
    fil.init_3(b_sam, b_dif, b_fac, true)?;

    // `occ::handle<Geom_Surface> Surfcoin = fil.Surface(); TopAbs_Orientation
    // Osurfsam = sdsam->Orientation();` (`cxx:983-984`).
    let surfcoin = fil
        .surface()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: GeomFill_ConstrainedFilling::Surface()".to_string())?;
    let surfcoin: Arc<dyn Surface> = Arc::new(surfcoin);
    let o_surf_sam = sorted.sd_sam.orientation();

    // `done = CompleteData(coin, Surfcoin, Hsurfsam, pcsurfsam, HBRFopsam,
    // pcnul, Osurfsam, true, false, false, false, false);` (`cxx:986-997`).
    let pc_surf_sam = Some(pc_surf_sam);
    let done = base.complete_data_geom(
        coin,
        &surfcoin,
        &h_surf_sam,
        pc_surf_sam.as_ref(),
        &sorted.hbr_fopsam,
        None,
        o_surf_sam,
        true,
        false,
        false,
        false,
        false,
    )?;

    // `if (!done) throw Standard_Failure("concavites inverted : fail");`
    // (`cxx:1001-1004`).
    if !done {
        return Err("concavites inverted : fail".to_string());
    }

    // `cxx:1006-1144`.
    ds_write::write_ds(base, corner, coin, &mut sorted)?;

    // OCCT writes straight into the stripes' `SurfData`; the port carries
    // copies, so they go back to the halves they came from.
    if sam_is_1 {
        ctx.sd1 = sorted.sd_sam;
        ctx.sd2 = sorted.sd_dif;
    } else {
        ctx.sd2 = sorted.sd_sam;
        ctx.sd1 = sorted.sd_dif;
    }
    ctx.flush()?;
    Ok(true)
}
