//! `ChFi3d_ChBuilder::PerformThreeCorner` - the `c1plan` arm.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder_C3.cxx:573-796`
//! (OCCT 8.0.0). The corner is a piece of the plane through the three points,
//! so the three edges against the pivot face and the two stripe surfaces are
//! `ComputeIntersection` results, and `coin` is filled with p-curves rather
//! than a `GeomFill_ConstrainedFilling` patch.
//!
//! Translation notes:
//! - `done = true` (`cxx:795`) is [`ChFi3dBuilder::done`](super::super::ChFi3dBuilder),
//!   set by the caller after this function answers,
//! - `ComputeIntersection` writes its `deru` / `derv` out-parameters
//!   (`cxx:619-621`); only the `coin-pivot` call is read afterwards, for
//!   `norpiv` (`cxx:631`), so the other two answers are dropped,
//! - `GeomInt_IntSS::Perform(S1, S2, Precision::Intersection())` (`cxx:709`)
//!   takes the header defaults `Approx = true`, `ApproxS1 = ApproxS2 = false`
//!   (`GeomInt_IntSS.hxx:54-59`),
//! - the `nbl > 1` arm of `cxx:714-718` only prints under `OCCT_DEBUG` and
//!   leaves `gcface` null, and `nbl == 0` leaves it null too; both then reach
//!   `gcface->FirstParameter()` of `cxx:777` on a null handle. The port states
//!   that as the `Err` of this function,
//! - `gcface->FirstParameter()` / `LastParameter()` are the `uf` / `ul`
//!   [`ch_fi3d_trim_curve_range`] answers; OCCT stores them as given.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt2d, GpVec};
use occt_core::precision::INTERSECTION;
use occt_geom::curve::Curve;
use occt_geom::plane::GeomPlane;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{ch_fi3d_project_pcurv, ch_fi3d_trim_curve_range};
use crate::chfi_ds::SurfData;
use crate::chfi_k_part::fcts::{ch_fi_kpart_index_curve_in_ds, ch_fi_kpart_index_surface_in_ds};
use crate::geom_int::IntSS;
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::super::corner_sd::{pcurve_on_surf, sd};
use super::super::super::ChFi3dBuilder;
use super::statics::{self, Intersection};
use super::{CornerProbe, PlanCorner, ThreeChamfer};

/// `gp_Dir::Dot(const gp_Vec&)`: the coordinate inner product, without the
/// normalisation of `gp_Dir` (a direction is already a unit vector).
fn dir_dot_vec(d: &GpDir, v: &GpVec) -> f64 {
    d.x() * v.x() + d.y() * v.y() + d.z() * v.z()
}

/// `if (!ComputeIntersection(...)) { throw StdFail_NotDone(message); }`
/// (`cxx:611-629`, `:647-663`, `:684-699`).
fn intersection_or_not_done(
    base: &mut ChFi3dBuilder,
    sd: &SurfData,
    sd_coin: &SurfData,
    pdeb: &occt_core::gp::GpPnt,
    p2ddeb: &GpPnt2d,
    pfin: &occt_core::gp::GpPnt,
    p2dfin: &GpPnt2d,
    message: &str,
) -> Result<Intersection, String> {
    let tol3d = base.tolapp3d;
    let tol2d = base.tol2d;
    let answer = statics::compute_intersection(
        base.my_ds.change_ds(),
        sd,
        sd_coin,
        pdeb,
        p2ddeb,
        pfin,
        p2dfin,
        tol3d,
        tol2d,
    )?;
    if !answer.done {
        return Err(format!("StdFail_NotDone: {message}"));
    }
    Ok(answer)
}

/// `if (c1plan) { ... done = true; }` (`cxx:573-796`).
pub(super) fn build_plan_corner(
    base: &mut ChFi3dBuilder,
    c: &ThreeChamfer,
    coin: &mut SurfData,
    probe: &CornerProbe,
    nor: &GpVec,
    v1: &GpVec,
) -> Result<PlanCorner, String> {
    // `gp_Dir ndir(nor);` (`cxx:577`): the `gp_Dir(gp_Vec)` constructor raises
    // `Standard_ConstructionError` on a null vector.
    let ndir = GpDir::from_vec(nor)
        .map_err(|_| "Standard_ConstructionError: gp_Dir(nor) is null".to_string())?;
    // `gp_Dir xdir = gp_Dir(gp_Vec(p3d[fin], p3d[deb]));` (`cxx:579`).
    let xvec = GpVec::from_pnts(&probe.p3d[c.fin], &probe.p3d[c.deb]);
    let xdir = GpDir::from_vec(&xvec)
        .map_err(|_| "Standard_ConstructionError: gp_Dir(p3d[fin], p3d[deb]) is null".to_string())?;

    // `gp_Ax3 planAx3(p3d[pivot], ndir, xdir);` (`cxx:580`).
    let mut plan_ax3 = GpAx3::new(probe.p3d[c.pivot], ndir, &xdir)
        .map_err(|e| format!("Standard_ConstructionError: {e}"))?;

    // `if (planAx3.YDirection().Dot(v1) <= 0.) { planAx3.YReverse(); }`
    // (`cxx:581-584`).
    if dir_dot_vec(plan_ax3.y_direction(), v1) <= 0.0 {
        plan_ax3.y_reverse();
    }

    // `occ::handle<Geom_Plane> gpl = new Geom_Plane(planAx3);
    // coin->ChangeSurf(ChFiKPart_IndexSurfaceInDS(gpl, DStr));` (`cxx:585-586`).
    let gpl: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(plan_ax3)));
    coin.set_surf(ch_fi_kpart_index_surface_in_ds(
        Arc::clone(&gpl),
        base.my_ds.change_ds(),
    ));

    // `gp_Vec norface = norpl; if (face[pivot].Orientation() == TopAbs_REVERSED)
    // { norface.Reverse(); }` (`cxx:589-593`).
    let mut norface = probe.norpl;
    let pivot_face = c.face[c.pivot]
        .clone()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no face".to_string())?;
    if pivot_face.orientation() == Orientation::Reversed {
        norface.reverse();
    }

    // `gp_Vec norcoin = gpl->Pln().Position().XDirection().Crossed(
    // gpl->Pln().Position().YDirection());` (`cxx:594-595`).
    let norcoin = plan_ax3
        .x_direction()
        .crossed(plan_ax3.y_direction())
        .map_err(|e| format!("Standard_ConstructionError: norcoin ({e})"))?;

    // `if (norcoin.Dot(norface) <= 0.) { coin->ChangeOrientation() =
    // TopAbs_REVERSED; } else { ... = TopAbs_FORWARD; }` (`cxx:596-601`).
    coin.set_orientation(if dir_dot_vec(&norcoin, &norface) <= 0.0 {
        Orientation::Reversed
    } else {
        Orientation::Forward
    });

    // `occ::handle<Geom_Curve> gcpiv, gcdeb, gcfin; occ::handle<Geom_TrimmedCurve>
    // gcface; occ::handle<Geom2d_Curve> pivpc1, pivpc2, debpc2, finpc2,
    // facepc1, facepc2;` (`cxx:604-607`).
    let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
    let fdfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
    let fdpiv = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;

    // `// intersection coin-pivot` (`cxx:610-630`).
    let mut piv = intersection_or_not_done(
        base,
        &fdpiv,
        coin,
        &probe.p3d[c.fin],
        &probe.p2d[c.fin],
        &probe.p3d[c.deb],
        &probe.p2d[c.deb],
        "echec calcul intersection coin-pivot",
    )?;
    // `gp_Vec norpiv = deru.Crossed(derv);` (`cxx:631`): `deru` / `derv` were
    // overwritten by the call above.
    let norpiv = piv.derudeb.crossed(&piv.dervdeb);
    let gcpiv = piv
        .gc
        .clone()
        .ok_or_else(|| "StdFail_NotDone: echec calcul intersection coin-pivot".to_string())?;
    let pivpc1 = piv.pc1.take();
    let pivpc2 = piv.pc2.take();

    // `// intersection coin-deb` (`cxx:633-663`).
    let p2d1 = if c.c1triangle {
        pcurve_on_surf(&fddeb, c.jf[c.deb][c.fin], c.p[c.deb][c.fin])?
    } else {
        pcurve_on_surf(&fddeb, c.jf[c.deb][c.fin], c.p[c.deb][c.pivot])?
    };
    let p2d2 = pcurve_on_surf(&fddeb, c.jf[c.deb][c.pivot], c.p[c.deb][c.pivot])?;
    let mut deb = intersection_or_not_done(
        base,
        &fddeb,
        coin,
        &probe.p3d[c.pivot],
        &p2d1,
        &probe.p3d[c.fin],
        &p2d2,
        "echec calcul intersection coin-deb",
    )?;
    let gcdeb = deb
        .gc
        .clone()
        .ok_or_else(|| "StdFail_NotDone: echec calcul intersection coin-deb".to_string())?;
    let debpc1 = deb.pc1.take();
    // `Icf = DStr.AddCurve(TopOpeBRepDS_Curve(gcdeb, tolrcoindeb));`
    // (`cxx:665`). `debpc2` is never read again, so it is not kept.
    let icf = base.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(gcdeb),
        deb.tolreached,
        false,
    ));

    // `// intersection coin-fin` (`cxx:667-701`).
    let (p3dface, p2d1) = if c.c1triangle {
        (
            probe.p3d[c.pivot],
            pcurve_on_surf(&fdfin, c.jf[c.fin][c.deb], c.p[c.fin][c.deb])?,
        )
    } else {
        (
            probe.p3d[3],
            pcurve_on_surf(&fdfin, c.jf[c.fin][c.deb], c.p[c.fin][c.pivot])?,
        )
    };
    let p2d2 = pcurve_on_surf(&fdfin, c.jf[c.fin][c.pivot], c.p[c.fin][c.pivot])?;
    let mut fin = intersection_or_not_done(
        base,
        &fdfin,
        coin,
        &p3dface,
        &p2d1,
        &probe.p3d[c.deb],
        &p2d2,
        "echec calcul intersection coin-face",
    )?;
    let gcfin = fin
        .gc
        .clone()
        .ok_or_else(|| "StdFail_NotDone: echec calcul intersection coin-face".to_string())?;
    let finpc1 = fin.pc1.take();
    // `Icl = DStr.AddCurve(TopOpeBRepDS_Curve(gcfin, tolrcoinfin));`
    // (`cxx:701`).
    let icl = base.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(gcfin),
        fin.tolreached,
        false,
    ));

    // `//! c1triangle: intersection coin-face[pivot]` (`cxx:704-737`).
    let mut gcface: Option<Arc<dyn Curve>> = None;
    let mut gcface_range = (0.0f64, 0.0f64);
    let mut facepc1: Option<Arc<dyn Curve2d>> = None;
    let mut facepc2: Option<Arc<dyn Curve2d>> = None;
    if !c.c1triangle {
        // `BRepAdaptor_Surface facebid(face[pivot]); occ::handle<Geom_Surface>
        // surfbid = facebid.GeomSurfaceTransformed();` (`cxx:707-708`).
        let fac_surf = probe
            .fac
            .surface()
            .cloned()
            .ok_or_else(|| "PerformThreeCorner: face[pivot] has no surface".to_string())?;

        // `GeomInt_IntSS inter; inter.Perform(gpl, surfbid,
        // Precision::Intersection());` (`cxx:705-709`).
        let mut inter = IntSS::new();
        inter.perform(
            Arc::clone(&gpl),
            Arc::clone(&fac_surf),
            INTERSECTION,
            true,
            false,
            false,
        );

        // `if (inter.IsDone()) { int nbl = inter.NbLines(); ... }`
        // (`cxx:710-737`).
        if inter.is_done() && inter.nb_lines() == 1 {
            let line = inter
                .line(1)
                .map(|l| Arc::clone(&l.curve))
                .ok_or_else(|| {
                    "StdFail_NotDone: GeomInt_IntSS::Line(1) is null".to_string()
                })?;
            // `ChFi3d_TrimCurve(inter.Line(1), p3d[pivot], p3dface, gcface);`
            // (`cxx:721`).
            let (trimmed, uf, ul) =
                ch_fi3d_trim_curve_range(&line, &probe.p3d[c.pivot], &p3dface)?;
            gcface_range = (uf, ul);
            let gcface_arc: Arc<dyn Curve> = Arc::new(trimmed);
            // `ChFi3d_ProjectPCurv(gac, gaf, facepc1, tolapp3d, tolr);` and
            // `ChFi3d_ProjectPCurv(gac, gas, facepc2, tolapp3d, tolr);`
            // (`cxx:731-732`).
            let tol = base.tolapp3d;
            ch_fi3d_project_pcurv(&gcface_arc, &fac_surf, &mut facepc1, tol)?;
            ch_fi3d_project_pcurv(&gcface_arc, &gpl, &mut facepc2, tol)?;
            gcface = Some(gcface_arc);
        } else if inter.is_done() {
            // `if (nbl > 1) { ... }` (`cxx:714-718`) and `nbl == 0` both leave
            // `gcface` null, which the fill of `cxx:777` dereferences.
            return Err(format!(
                "StdFail_NotDone: GeomInt_IntSS answered {} line(s) between the corner plane \
                 and face[pivot]",
                inter.nb_lines()
            ));
        }
    }

    // `// on remplit les donnees du coin oriente face-pivot` (`cxx:739-750`).
    coin.vertex_first_on_s1_mut().set_point(&probe.p3d[c.pivot]);
    coin.vertex_first_on_s2_mut().set_point(&probe.p3d[c.fin]);
    if c.c1triangle {
        coin.vertex_last_on_s1_mut().set_point(&probe.p3d[c.pivot]);
    } else {
        coin.vertex_last_on_s1_mut().set_point(&probe.p3d[3]);
    }
    coin.vertex_last_on_s2_mut().set_point(&probe.p3d[c.deb]);

    // `// sur face[pivot]` (`cxx:757-779`).
    let mut trans = if dir_dot_vec(&norcoin, &probe.norpl) <= 0.0 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    let bidpc: Option<Arc<dyn Curve2d>> = None;
    if c.c1triangle {
        // `fi1.SetInterference(0, trans, bidpc, bidpc);` (`cxx:769-772`).
        coin.interference_on_s1_mut()
            .set_interference(0, trans, bidpc.clone(), bidpc);
    } else {
        // `Igcface = ChFiKPart_IndexCurveInDS(gcface, DStr);
        // fi1.SetInterference(Igcface, trans, facepc1, facepc2);
        // fi1.SetFirstParameter(gcface->FirstParameter()); ...`
        // (`cxx:774-778`).
        let igcface = ch_fi_kpart_index_curve_in_ds(gcface.clone(), base.my_ds.change_ds());
        let fi1 = coin.interference_on_s1_mut();
        fi1.set_interference(igcface, trans, facepc1.clone(), facepc2.clone());
        fi1.set_first_parameter(gcface_range.0);
        fi1.set_last_parameter(gcface_range.1);
    }

    // `// sur le pivot` (`cxx:780-793`).
    if dir_dot_vec(&norcoin, &norpiv) <= 0.0 {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }
    let igcpiv = ch_fi_kpart_index_curve_in_ds(Some(Arc::clone(&gcpiv)), base.my_ds.change_ds());
    let fi2 = coin.interference_on_s2_mut();
    fi2.set_interference(igcpiv, trans, pivpc1, pivpc2);
    fi2.set_first_parameter(gcpiv.first_parameter());
    fi2.set_last_parameter(gcpiv.last_parameter());

    Ok(PlanCorner {
        icf,
        icl,
        debpc1,
        finpc1,
    })
}
