//! The three file statics of `ChFi3d_ChBuilder_C3.cxx` (OCCT 8.0.0).
//!
//! - [`co_planar`] - `CoPlanar` (`cxx:64-88`): are four points coplanar?
//! - [`bound_surf`] - `BoundSurf` (`cxx:96-114`): a `GeomAdaptor_Surface` over
//!   `S` trimmed to the box of two `(u, v)` points, enlarged by `20%` of the
//!   larger side and clamped to the surface's own bounds,
//! - [`compute_intersection`] - `ComputeIntersection` (`cxx:131-197`): the
//!   intersection of one of the three `SurfData` with the corner `SDCoin`,
//!   between two known extremities.
//!
//! Translation notes:
//! - `gp_Vec& derudeb` / `dervdeb` are out-parameters of `ComputeIntersection`
//!   (the caller keeps `deru` / `derv` for `norpiv`); they travel in
//!   [`Intersection`],
//! - `NCollection_Array1<double> Pardeb(1, 4)` is `[f64; 4]` counted from `0`,
//!   the order [`ch_fi3d_compute_curves`] takes,
//! - the 10-argument `ChFi3d_ComputeCurves` call of `cxx:194` leaves
//!   `wholeCurv` at its `false` default.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_boite, ch_fi3d_bound_surf, ch_fi3d_compute_curves, ch_fi3d_parameters, ComputeCurves,
};
use crate::chfi_ds::SurfData;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `static bool CoPlanar(const gp_Pnt& PntA, const gp_Pnt& PntB, const gp_Pnt&
/// PntC, const gp_Pnt& PntD)` (`cxx:64-88`).
///
/// The distance of `PntD` to the plane of `PntA` / `PntB` / `PntC`, relative to
/// `Alpha` (the squared area of that triangle times two), is compared with
/// `Precision::Confusion()`.
pub(super) fn co_planar(pnt_a: &GpPnt, pnt_b: &GpPnt, pnt_c: &GpPnt, pnt_d: &GpPnt) -> bool {
    let vec_ab = GpVec::from_pnts(pnt_a, pnt_b);
    let vec_ac = GpVec::from_pnts(pnt_a, pnt_c);
    let vec_ad = GpVec::from_pnts(pnt_a, pnt_d);

    // `double nor2AB = vecAB.SquareMagnitude(); ... double Alpha = nor2AB *
    // nor2AC - ProABAC * ProABAC;` (`cxx:70-75`).
    let nor2_ab = vec_ab.square_magnitude();
    let nor2_ac = vec_ac.square_magnitude();
    let pro_ab_ac = vec_ab.dot(&vec_ac);
    let alpha = nor2_ab * nor2_ac - pro_ab_ac * pro_ab_ac;

    // `if (Alpha < Precision::Confusion()) { return true; }` (`cxx:77-80`).
    if alpha < CONFUSION {
        return true;
    }

    // `double ProABAD = vecAB.Dot(vecAD); ... gp_Vec vecDABC = Alpha1 * vecAB +
    // Alpha2 * vecAC - Alpha * vecAD;` (`cxx:82-87`).
    let pro_ab_ad = vec_ab.dot(&vec_ad);
    let pro_ac_ad = vec_ac.dot(&vec_ad);
    let alpha1 = pro_ab_ad * nor2_ac - pro_ab_ac * pro_ac_ad;
    let alpha2 = pro_ac_ad * nor2_ab - pro_ab_ac * pro_ab_ad;
    let vec_dabc = vec_ab
        .multiplied_scalar(alpha1)
        .added(&vec_ac.multiplied_scalar(alpha2))
        .subtracted(&vec_ad.multiplied_scalar(alpha));

    // `return (vecDABC.Magnitude() / Alpha) < Precision::Confusion();`
    // (`cxx:88`).
    (vec_dabc.magnitude() / alpha) < CONFUSION
}

/// `static occ::handle<GeomAdaptor_Surface> BoundSurf(const occ::handle<Geom_
/// Surface>& S, const gp_Pnt2d& Pdeb, const gp_Pnt2d& Pfin)` (`cxx:96-114`).
pub(super) fn bound_surf(
    s: &Arc<dyn Surface>,
    pdeb: &GpPnt2d,
    pfin: &GpPnt2d,
) -> SurfaceAdaptor {
    // `GeomAdaptor_Surface& GAS = *HS; GAS.Load(S);` (`cxx:100-102`).
    // `S->Bounds(uuu1, uuu2, vvv1, vvv2);` (`cxx:107`).
    let (uuu1, uuu2) = s.u_range();
    let (vvv1, vvv2) = s.v_range();
    // `ChFi3d_Boite(Pdeb, Pfin, uu1, uu2, vv1, vv2);` (`cxx:108`).
    let b = ch_fi3d_boite(pdeb, pfin);
    let mut uu1 = b.u_min;
    let mut uu2 = b.u_max;
    let mut vv1 = b.v_min;
    let mut vv2 = b.v_max;
    // `double Step = std::max((uu2 - uu1), (vv2 - vv1)); Step *= 0.2;`
    // (`cxx:109-110`).
    let step = (uu2 - uu1).max(vv2 - vv1) * 0.2;
    // `uuu1 = std::max((uu1 - Step), uuu1); ...` (`cxx:111-114`).
    uu1 = (uu1 - step).max(uuu1);
    uu2 = (uu2 + step).min(uuu2);
    vv1 = (vv1 - step).max(vvv1);
    vv2 = (vv2 + step).min(vvv2);
    // `GAS.Load(S, uuu1, uuu2, vvv1, vvv2);` (`cxx:115`).
    SurfaceAdaptor::of_surface_with_range(Arc::clone(s), uu1, uu2, vv1, vv2)
}

/// The out-parameters of [`compute_intersection`].
pub(super) struct Intersection {
    /// `gc`.
    pub gc: Option<Arc<dyn Curve>>,
    /// `pc1`.
    pub pc1: Option<Arc<dyn Curve2d>>,
    /// `pc2`.
    pub pc2: Option<Arc<dyn Curve2d>>,
    /// `derudeb`.
    pub derudeb: GpVec,
    /// `dervdeb`.
    pub dervdeb: GpVec,
    /// `ptcoindeb` is only read by the `BoundSurf` call inside
    /// `ComputeIntersection`, so it is a local of that body and not a field.
    /// `tolreached`.
    pub tolreached: f64,
    /// The `bool` `ComputeIntersection` returns.
    pub done: bool,
}

/// `static bool ComputeIntersection(...)` (`cxx:131-197`).
///
/// `sd` is the `SurfData` of one of the three stripes, `sd_coin` the corner;
/// `pdeb` / `p2ddeb` and `pfin` / `p2dfin` are the two extremities, the first
/// `(u, v)` on `sd` and the second on `sd_coin` (`p2ddeb` on `sd`, `p2dfin` on
/// `sd` as well - see the `Pardeb` / `Parfin` fill of `cxx:159-183`).
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_intersection(
    ds: &mut DataStructure,
    sd: &SurfData,
    sd_coin: &SurfData,
    pdeb: &GpPnt,
    p2ddeb: &GpPnt2d,
    pfin: &GpPnt,
    p2dfin: &GpPnt2d,
    tol3d: f64,
    tol2d: f64,
) -> Result<Intersection, String> {
    // `HS1 = ChFi3d_BoundSurf(DStr, SD, 1, 2);` (`cxx:152-153`).
    let hs1 = ch_fi3d_bound_surf(ds, sd, 1, 2)?;

    // `const occ::handle<Geom_Surface>& gpl =
    // DStr.Surface(SDCoin->Surf()).Surface(); const ... gSD =
    // DStr.Surface(SD->Surf()).Surface();` (`cxx:155-156`).
    let gpl = ds
        .surface(sd_coin.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ComputeIntersection: DStr.Surface(SDCoin->Surf()) is null".to_string())?;
    let g_sd = ds
        .surface(sd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ComputeIntersection: DStr.Surface(SD->Surf()) is null".to_string())?;

    // `NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);` (`cxx:159`).
    let mut pardeb = [0.0f64; 4];
    let mut parfin = [0.0f64; 4];

    // `u = p2ddeb.X(); v = p2ddeb.Y(); gSD->D1(u, v, Pbidon, derudeb,
    // dervdeb); Pardeb(1) = u; Pardeb(2) = v;` (`cxx:161-167`).
    let (u, v) = (p2ddeb.x(), p2ddeb.y());
    let (_, derudeb, dervdeb) = g_sd.d1(u, v);
    pardeb[0] = u;
    pardeb[1] = v;

    // `ChFi3d_Parameters(gpl, pdeb, u, v); Pardeb(3) = u; Pardeb(4) = v;
    // ptcoindeb.SetCoord(u, v);` (`cxx:169-173`).
    let (u, v) = ch_fi3d_parameters(&gpl, pdeb)?;
    pardeb[2] = u;
    pardeb[3] = v;
    let ptcoindeb = GpPnt2d::new(u, v);

    // `u = p2dfin.X(); v = p2dfin.Y(); Parfin(1) = u; Parfin(2) = v;`
    // (`cxx:175-178`).
    parfin[0] = p2dfin.x();
    parfin[1] = p2dfin.y();

    // `ChFi3d_Parameters(gpl, pfin, u, v); Parfin(3) = u; Parfin(4) = v;
    // gp_Pnt2d cpf2(u, v);` (`cxx:181-185`).
    let (u, v) = ch_fi3d_parameters(&gpl, pfin)?;
    parfin[2] = u;
    parfin[3] = v;
    let cpf2 = GpPnt2d::new(u, v);

    // `HS2 = BoundSurf(gpl, ptcoindeb, cpf2);` (`cxx:189-190`).
    let hs2 = bound_surf(&gpl, &ptcoindeb, &cpf2);

    // `return ChFi3d_ComputeCurves(HS1, HS2, Pardeb, Parfin, gc, pc1, pc2,
    // tol3d, tol2d, tolreached);` (`cxx:193-194`).
    let mut out = ComputeCurves::default();
    let done = ch_fi3d_compute_curves(&hs1, &hs2, pardeb, parfin, &mut out, tol3d, tol2d, false);
    Ok(Intersection {
        gc: out.c3d,
        pc1: out.pc1,
        pc2: out.pc2,
        derudeb,
        dervdeb,
        tolreached: out.tolreached,
        done,
    })
}
