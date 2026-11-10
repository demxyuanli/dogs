//! `ChFi3d_ChBuilder::PerformThreeCorner` - the `!c1plan` arm.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder_C3.cxx:798-886`
//! (OCCT 8.0.0). The corner is a `GeomFill_ConstrainedFilling` patch between
//! the pivot face `Fac`, the pivot surface `Surf` and the two stripe surfaces,
//! then `CompleteData` on the reversed patch.
//!
//! Translation notes:
//! - `issmooth` is the dead `false` of `cxx:207`; both of its arms are
//!   translated, and the tangents it feeds (`Tgpiv` / `Tg3`) are the null
//!   vectors the `.cxx` default-constructs,
//! - `PcurveOnFace` / `PcurveOnPiv` are out-parameters of the `ChFi3d_mkbound`
//!   calls, which is why the `mkbound` ports answer the p-curve they built,
//! - `fil.Init(Bpiv, Bfin, Bfac, Bdeb)` (`cxx:869`) takes the header default
//!   `NoCheck = false`,
//! - `Geom_Surface::VReverse()` of `cxx:877` is an in-place mutation of the
//!   patch `Surface()` returned; the port asks the `Surface` trait for the
//!   reversed copy, which is what `CompleteData` reads,
//! - `IntV` / `IntPiv` (`cxx:807-808`, `:813-814`) feed only the commented-out
//!   `cxx:812-813` / `:947-948` call, so they are kept as the `D1` lookups
//!   without the unused vector.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::geom_fill::boundary::Boundary;
use occt_geom::geom_fill::ConstrainedFilling;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::{
    mkbound_curve_segment, mkbound_pfac2d, mkbound_pfac3d, mkbound_segment,
};
use crate::chfi_ds::SurfData;

use super::super::super::corner_sd::{data_surface, ds_curve, sd};
use super::super::super::ChFi3dBuilder;
use super::{CornerProbe, ThreeChamfer};

/// `CD[k]->SetOfSurfData()->Value(i)->Interference(jf).PCurveOnSurf()->Value(p)`.
fn p_curve_on_surf_value(
    sd: &SurfData,
    jf: i32,
    param: f64,
) -> Result<GpPnt2d, String> {
    sd.interference(jf)
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the interference has no p-curve".to_string())
        .map(|c| c.value(param))
}

/// `else { // !c1plan ... done = CompleteData(...); }` (`cxx:797-886`).
/// Answers the resulting `done`.
pub(super) fn build_corner(
    base: &mut ChFi3dBuilder,
    c: &ThreeChamfer,
    probe: &CornerProbe,
    coin: &mut SurfData,
) -> Result<bool, String> {
    // `bool issmooth = false;` (`cxx:207`).
    let issmooth = false;

    // `occ::handle<Geom_Surface> Surfcoin; occ::handle<Geom2d_Curve>
    // PCurveOnFace, PCurveOnPiv;` (`cxx:801-802`).
    let p_curve_on_face: Option<Arc<dyn Curve2d>>;
    let p_curve_on_piv: Option<Arc<dyn Curve2d>>;

    // `int ind1 = fddeb->Interference(jf[deb][pivot]).LineIndex(); int ind2 =
    // fdfin->Interference(jf[fin][pivot]).LineIndex();` (`cxx:806-808`).
    let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
    let fdfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
    let fdpiv = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;

    // `DStr.Curve(ind1).Curve()->D1(p[deb][pivot], Pfin, vpfin);` and the same
    // for `fin` (`cxx:813-814`). Only the `issmooth` arm of `cxx:817-839` reads
    // the two vectors.
    let ind1 = fddeb.interference(c.jf[c.deb][c.pivot]).line_index();
    let (_, vpfin) = ds_curve(base, ind1)?.d1(c.p[c.deb][c.pivot]);
    let ind2 = fdfin.interference(c.jf[c.fin][c.pivot]).line_index();
    let (_, vpdeb) = ds_curve(base, ind2)?.d1(c.p[c.fin][c.pivot]);

    // `occ::handle<GeomFill_Boundary> Bdeb, Bfin, Bpiv, Bfac;`
    // (`cxx:805`).
    let (bfac, bpiv);
    if issmooth {
        // `if (issmooth)` (`cxx:817-839`).
        let (pc, bound) = mkbound_pfac2d(
            &probe.fac,
            c.sens[c.deb],
            &probe.p2d[c.pivot],
            &probe.tgpiv,
            c.sens[c.fin],
            &probe.p2d[3],
            &probe.tg3,
            base.tolapp3d,
            2.0e-4,
        )?;
        p_curve_on_face = Some(pc);
        bfac = bound;
        let (pc, bound) = mkbound_pfac3d(
            &probe.surf,
            c.sens[c.deb],
            &probe.p2d[c.fin],
            &vpfin,
            c.sens[c.fin],
            &probe.p2d[c.deb],
            &vpdeb,
            base.tolapp3d,
            2.0e-4,
        )?;
        p_curve_on_piv = Some(pc);
        bpiv = bound;
    } else {
        // `Bfac = ChFi3d_mkbound(Fac, PCurveOnFace, p2d[pivot], p2d[3],
        // tolapp3d, 2.e-4); Bpiv = ChFi3d_mkbound(Surf, PCurveOnPiv, p2d[fin],
        // p2d[deb], tolapp3d, 2.e-4);` (`cxx:845-846`).
        let (pc, bound) = mkbound_curve_segment(
            &probe.fac,
            &probe.p2d[c.pivot],
            &probe.p2d[3],
            base.tolapp3d,
            2.0e-4,
            false,
        )?;
        p_curve_on_face = Some(pc);
        bfac = bound;
        let (pc, bound) = mkbound_curve_segment(
            &probe.surf,
            &probe.p2d[c.fin],
            &probe.p2d[c.deb],
            base.tolapp3d,
            2.0e-4,
            false,
        )?;
        p_curve_on_piv = Some(pc);
        bpiv = bound;
    }

    // `gp_Pnt2d pdeb1 = fddeb->Interference(jf[deb][pivot]).PCurveOnSurf()
    // ->Value(p[deb][pivot]); gp_Pnt2d pdeb2 = fddeb->Interference(
    // jf[deb][fin]).PCurveOnSurf()->Value(p[deb][pivot]); ...` (`cxx:849-852`).
    let pdeb1 = p_curve_on_surf_value(&fddeb, c.jf[c.deb][c.pivot], c.p[c.deb][c.pivot])?;
    let pdeb2 = p_curve_on_surf_value(&fddeb, c.jf[c.deb][c.fin], c.p[c.deb][c.pivot])?;
    let pfin1 = p_curve_on_surf_value(&fdfin, c.jf[c.fin][c.pivot], c.p[c.fin][c.pivot])?;
    let pfin2 = p_curve_on_surf_value(&fdfin, c.jf[c.fin][c.deb], c.p[c.fin][c.pivot])?;

    // `DStr.Surface(fddeb->Surf()).Surface()` and the same for `fdfin`
    // (`cxx:855`, `:862-865`).
    let sdeb = data_surface(base.my_ds.ds(), &fddeb)?;
    let sfin = data_surface(base.my_ds.ds(), &fdfin)?;

    let (bdeb, bfin): (Arc<dyn Boundary>, Arc<dyn Boundary>) = if issmooth {
        // `Bdeb = ChFi3d_mkbound(DStr.Surface(fddeb->Surf()).Surface(), pdeb1,
        // pdeb2, tolapp3d, 2.e-4); ...` (`cxx:855-856`).
        (
            mkbound_segment(&sdeb, &pdeb1, &pdeb2, base.tolapp3d, 2.0e-4, false)?,
            mkbound_segment(&sfin, &pfin1, &pfin2, base.tolapp3d, 2.0e-4, false)?,
        )
    } else {
        // `Bdeb = ChFi3d_mkbound(..., true);` (`cxx:862-865`).
        (
            mkbound_segment(&sdeb, &pdeb1, &pdeb2, base.tolapp3d, 2.0e-4, true)?,
            mkbound_segment(&sfin, &pfin1, &pfin2, base.tolapp3d, 2.0e-4, true)?,
        )
    };

    // `GeomFill_ConstrainedFilling fil(8, 20); fil.Init(Bpiv, Bfin, Bfac,
    // Bdeb);` (`cxx:867-869`).
    let mut fil = ConstrainedFilling::new(8, 20);
    fil.init_4(bpiv, bfin, bfac, bdeb, false)?;

    // `Surfcoin = fil.Surface(); Surfcoin->VReverse();` (`cxx:871-877`).
    let surfcoin: Arc<dyn Surface> = fil
        .surface()
        .ok_or_else(|| "GeomFill_ConstrainedFilling::Surface".to_string())?
        .v_reversed()
        .ok_or_else(|| "Geom_Surface::VReverse".to_string())?;

    // `done = CompleteData(coin, Surfcoin, Fac, PCurveOnFace, Surf,
    // PCurveOnPiv, fdpiv->Orientation(), false, false, false, false, false);`
    // (`cxx:879-886`).
    let o_fdpiv = fdpiv.orientation();
    base.done = base.complete_data_geom(
        coin,
        &surfcoin,
        &probe.fac,
        p_curve_on_face.as_ref(),
        &probe.surf,
        p_curve_on_piv.as_ref(),
        o_fdpiv,
        false,
        false,
        false,
        false,
        false,
    )?;

    Ok(base.done)
}
