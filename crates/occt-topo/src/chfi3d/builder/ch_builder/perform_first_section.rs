//! `ChFi3d_ChBuilder::PerformFirstSection` — the initial section of a chamfer
//! band when no `KPart` is available.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:1319-1614`
//! (OCCT 8.0.0). The body is split in three by `chsp->IsChamfer()`
//! (`cxx:1340`, `:1419`, `:1535`):
//!
//! | arm | blend function | `Set` |
//! |---|---|---|
//! | `ChFiDS_Sym` | `Chamfer` / `ConstThroat` by `chsp->Mode()` (`cxx:1346-1355`) | `(dis, dis, Choix)` (`cxx:1356`) |
//! | `ChFiDS_TwoDist` | `Chamfer` (`cxx:1427-1429`) / `ConstThroatWithPenetration` on the offset guideline (`cxx:1450-1452`) | `(dis1, dis2, Choix)` / `(max, max, Choix)` |
//! | `ChFiDS_DistAngle` | `ChAsym` (`cxx:1541-1542`) | `(dis, angle, Choix)` (`cxx:1543`) |
//!
//! All three then run the same tail: `HGuide->D1(Par, ...)`, `Set(Par)`,
//! `Tangent(...)`, the `Choix % 2` / sign-derived reversals of the two
//! tangents (`cxx:1369-1395`, `:1467-1493`, `:1557-1583`), the two points one
//! distance away from the guideline point, two `Extrema_GenLocateExtPS`
//! projections that refresh `SolDep`, and `BRepBlend_Walking::PerformFirstSection`
//! for the actual solve. The tail is [`approximate_start`] plus
//! [`project_and_walk`]; only the two distances differ between the arms.
//!
//! `ChFiDS_TwoDist` + `ConstThroatWithPenetration` is the only arm that
//! substitutes the **offset** guideline for the blend function (`cxx:1431-1452`)
//! while the `Walking` still takes `HGuide` (`cxx:1453`). The distance
//! recomputation of `cxx:1497-1509` is commented out in OCCT (the `aDist1` /
//! `aDist2` locals are assigned `dis1` / `dis2` and the block is a comment), so
//! the port uses `dis1` / `dis2` verbatim.
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("PerformSurf : this is not the spine of
//!   a chamfer")` (`cxx:1332-1336`) is the `Err` of
//!   [`perform_first_section`]; the trait entry point answers `false`,
//! - the "Construction of offset guide failed!" print of `cxx:1445-1449` has
//!   its exception commented out and would dereference a null handle; the port
//!   answers `Err` with that message instead (the same choice as
//!   [`super::simul_surf`]),
//! - `Occ::handle` sharing: the port clones the offset guideline out of
//!   `ChFiDS_Spine::offset_el_spines`, so nothing is reframed in place here,
//! - `Extrema_GenLocateExtPS::Point().Parameter(U, V)` (`cxx:1410-1416`,
//!   `:1525-1531`, `:1606-1611`) is
//!   [`occt_geom::extrema_surf::GenLocateExtPs::point`] — read only after
//!   `IsDone`, which is what the `.cxx` guard does,
//! - `TopAbs_State& Pos1 / Pos2` is [`State`] at the trait boundary and
//!   [`FaceState`] at the `BRepBlend_Walking` boundary; [`state_face`] and
//!   [`state_of`] are the same rename the fillet subclass uses.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::extrema_surf::GenLocateExtPs;
use occt_geom::surface::Surface;
use occt_math::MathVector;

use crate::abs::State;
use crate::blend::{BlendFunction, ChAsym, Chamfer, ConstThroat, ConstThroatWithPenetration};
use crate::brep_blend::walking::Walking;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::fil_builder::{state_face, state_of};
use crate::chfi3d::builder::ChFi3dCtx;
use crate::chfi_ds::{ChamfMethod, ChamfMode, ElSpine, SpineClass, SpineHandle};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::simul_surf::guide_as_curve;

/// `ChFi3d_ChBuilder::PerformFirstSection(Spine, HGuide, Choix, S1, S2, I1, I2,
/// Par, SolDep, Pos1, Pos2)` (`cxx:1319-1614`).
///
/// `false` is OCCT's "section not found", which `StartSol` acts on.
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_first_section(
    spine: &SpineHandle,
    guide: &Arc<ElSpine>,
    choix: i32,
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    i2: &Rc<RefCell<TopolTool>>,
    par: f64,
    sol_dep: &mut MathVector,
    pos1: &mut State,
    pos2: &mut State,
    ctx: &ChFi3dCtx<'_>,
) -> Result<bool, String> {
    // `chsp = down_cast<ChFiDS_ChamfSpine>(Spine); if (chsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of a
    // chamfer");` (`cxx:1332-1336`).
    let is_chamf = {
        let guard = spine
            .read()
            .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?;
        guard.class() == SpineClass::Chamf
    };
    if !is_chamf {
        return Err("PerformSurf : this is not the spine of a chamfer".to_string());
    }

    // `double TolGuide = HGuide->Resolution(tolapp3d);` (`cxx:1338`).
    let tol_guide = guide
        .resolution(ctx.tolapp3d)
        .ok_or_else(|| "PerformFirstSection: HGuide has no curve".to_string())?;

    // The two adaptor surfaces, shared by the blend function and the walker.
    let s1_surf = s1
        .surface()
        .ok_or_else(|| "PerformFirstSection: S1 has no surface".to_string())?
        .clone();
    let s2_surf = s2
        .surface()
        .ok_or_else(|| "PerformFirstSection: S2 has no surface".to_string())?
        .clone();

    // `chsp->IsChamfer()`, `chsp->Mode()`, `GetDist` / `Dists` /
    // `GetDistAngle` (`cxx:1340-1344`, `:1421-1423`, `:1537-1538`). A spine of
    // the wrong kind would raise inside those calls; `IsChamfer()` already
    // told which one applies.
    let (method, mode, dis1, dis2, angle) = {
        let guard = spine
            .read()
            .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?;
        let method = guard.is_chamfer();
        let mode = guard.mode();
        let (d1, d2, a) = match method {
            ChamfMethod::Sym => (
                guard
                    .get_dist()
                    .ok_or_else(|| "PerformFirstSection: GetDist raised".to_string())?,
                0.0,
                0.0,
            ),
            ChamfMethod::TwoDist => {
                let (d1, d2) = guard
                    .dists()
                    .ok_or_else(|| "PerformFirstSection: Dists raised".to_string())?;
                (d1, d2, 0.0)
            }
            ChamfMethod::DistAngle => {
                let (d1, a) = guard
                    .get_dist_angle()
                    .ok_or_else(|| "PerformFirstSection: GetDistAngle raised".to_string())?;
                (d1, 0.0, a)
            }
        };
        (method, mode, d1, d2, a)
    };

    // The blend function of the arm. `BRepBlend_Walking TheWalk(S1, S2, I1, I2,
    // HGuide)` (`cxx:1357`, `:1453`, `:1544`) is built on `HGuide` in all three
    // arms.
    let mut func: Box<dyn BlendFunction> = match method {
        ChamfMethod::Sym => {
            // `if (chsp->Mode() == ChFiDS_ClassicChamfer) pFunc.reset(new
            // BRepBlend_Chamfer(S1, S2, HGuide)); else pFunc.reset(new
            // BRepBlend_ConstThroat(S1, S2, HGuide));` (`cxx:1346-1355`), then
            // `pFunc->Set(dis, dis, Choix);` (`cxx:1356`).
            let dis = dis1;
            if mode == ChamfMode::ClassicChamfer {
                let mut f = Chamfer::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                f.set_distances(dis, dis, choix);
                Box::new(f)
            } else {
                let mut f = ConstThroat::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                f.set_throat(dis, dis, choix);
                Box::new(f)
            }
        }
        ChamfMethod::TwoDist => {
            if mode == ChamfMode::ClassicChamfer {
                // `pFunc.reset(new BRepBlend_Chamfer(S1, S2, HGuide));
                // pFunc->Set(dis1, dis2, Choix);` (`cxx:1427-1429`).
                let mut f = Chamfer::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                f.set_distances(dis1, dis2, choix);
                Box::new(f)
            } else {
                // `OffsetHGuide` from the parallel walk of `ChangeElSpines()` /
                // `ChangeOffsetElSpines()` (`cxx:1431-1442`); a null handle
                // raises the "Construction of offset guide failed!" message
                // (`cxx:1445-1449`).
                let offset_arc: Arc<ElSpine> = {
                    let sg = spine
                        .read()
                        .map_err(|_| "PerformFirstSection: poisoned spine".to_string())?;
                    let offset = sg
                        .offset_el_spine(guide)
                        .ok_or_else(|| "Construction of offset guide failed!".to_string())?;
                    Arc::clone(offset)
                };
                // `pFunc.reset(new BRepBlend_ConstThroatWithPenetration(S1, S2,
                // OffsetHGuide)); double Throat = std::max(dis1, dis2);
                // pFunc->Set(Throat, Throat, Choix);` (`cxx:1450-1452`).
                let mut f = ConstThroatWithPenetration::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(&offset_arc),
                );
                let throat = dis1.max(dis2);
                f.set_throat(throat, throat, choix);
                Box::new(f)
            }
        }
        ChamfMethod::DistAngle => {
            // `BRepBlend_ChAsym Func(S1, S2, HGuide); Func.Set(dis1, angle,
            // Ch);` (`cxx:1541-1543`).
            let mut f = ChAsym::new(
                Arc::clone(&s1_surf),
                Arc::clone(&s2_surf),
                guide_as_curve(guide),
            );
            f.set_dist_angle(dis1, angle, choix);
            Box::new(f)
        }
    };

    // `HGuide->D1(Par, ptgui, d1gui); ... Func.Set(Par); Func.Tangent(...)`
    // plus the two reversals (`cxx:1362-1395`, `:1461-1493`, `:1551-1583`).
    let (ptgui, tgf, tgl) = approximate_start(&mut *func, guide, choix, par, sol_dep)?;

    // The two points of the approximate section, one distance away from
    // `ptgui` along each tangent; the distances are the only per-arm part
    // (`cxx:1397-1402`, `:1512-1517`, `:1585-1597`).
    let (pt1, pt2) = match method {
        ChamfMethod::Sym | ChamfMethod::TwoDist => {
            // `temp = TgF.XYZ().Multiplied(aDist1); pt1 = ptgui + temp;`
            // (`cxx:1397-1398`, `:1512-1513`) and the same with `aDist2` on
            // `TgL` (`cxx:1399-1400`, `:1515-1516`). The `TwoDist` distances
            // are `dis1` / `dis2` — `aDist1` / `aDist2` of `cxx:1494-1509`.
            let a_dist1 = dis1;
            let a_dist2 = dis2;
            (
                ptgui.translated_vec(&tgf.multiplied_scalar(a_dist1)),
                ptgui.translated_vec(&tgl.multiplied_scalar(a_dist2)),
            )
        }
        ChamfMethod::DistAngle => {
            // `temp = TgF.XYZ().Multiplied(dis1); pt1 = ptgui + temp;`
            // (`cxx:1585-1586`), then `tmpcos = TgF.Dot(TgL); tmpsin =
            // sqrt(1 - tmpcos*tmpcos); dis2 = dis1 / (tmpcos + tmpsin /
            // tan(angle)); temp = TgL.XYZ().Multiplied(dis2); pt2 = ptgui +
            // temp;` (`cxx:1589-1597`).
            let p1 = ptgui.translated_vec(&tgf.multiplied_scalar(dis1));
            let tmpcos = tgf.dot(&tgl);
            let tmpsin = (1.0 - tmpcos * tmpcos).sqrt();
            let d2 = dis1 / (tmpcos + tmpsin / angle.tan());
            let p2 = ptgui.translated_vec(&tgl.multiplied_scalar(d2));
            (p1, p2)
        }
    };

    project_and_walk(
        &mut *func, &s1_surf, &s2_surf, i1, i2, guide, par, sol_dep, pos1, pos2, ctx, pt1, pt2,
        tol_guide,
    )
}

/// The shared head of the three arms: `HGuide->D1(Par, ptgui, d1gui)`,
/// `Set(Par)`, `Tangent(SolDep(1), SolDep(2), SolDep(3), SolDep(4), ...)`, the
/// `Choix % 2` reversal and the sign-based swap of the two flags
/// (`cxx:1362-1395`, `:1461-1493`, `:1551-1583`).
///
/// Returns the guideline point and the two oriented tangents.
fn approximate_start(
    func: &mut dyn BlendFunction,
    guide: &Arc<ElSpine>,
    choix: i32,
    par: f64,
    sol_dep: &MathVector,
) -> Result<(GpPnt, GpVec, GpVec), String> {
    // `HGuide->D1(Par, ptgui, d1gui);` (`cxx:1364`, `:1463`, `:1553`). The
    // `.cxx` only uses `d1gui` for the sign test below. `ChFiDS_ElSpine::D1`
    // raises `Standard_NullObject` on a spine without a curve.
    let (ptgui, d1gui) = guide
        .eval_d1(par)
        .ok_or_else(|| "PerformFirstSection: ChFiDS_ElSpine::D1 raised".to_string())?;

    // `pFunc->Set(Par);` (`cxx:1367`, `:1466`, `:1556`).
    func.set(par);

    // `pFunc->Tangent(SolDep(1), SolDep(2), SolDep(3), SolDep(4), TgF, TgL,
    // tmp1, tmp2);` (`cxx:1368`, `:1467`, `:1557`). `tmp1` / `tmp2` are the
    // normals, unused by the `.cxx`.
    let (mut tgf, mut tgl) = (GpVec::zero(), GpVec::zero());
    let (mut tmp1, mut tmp2) = (GpVec::zero(), GpVec::zero());
    func.tangent(
        sol_dep.value(1),
        sol_dep.value(2),
        sol_dep.value(3),
        sol_dep.value(4),
        &mut tgf,
        &mut tgl,
        &mut tmp1,
        &mut tmp2,
    );

    // `bool rev1 = false; bool rev2 = false; double sign =
    // (TgF.Crossed(d1gui)).Dot(TgL); if (Choix % 2 == 1) rev1 = true; else
    // rev2 = true; if (sign < 0.) { rev1 = !rev1; rev2 = !rev2; }`
    // (`cxx:1369-1385`).
    let mut rev1 = false;
    let mut rev2 = false;
    let sign = tgf.crossed(&d1gui).dot(&tgl);
    if choix % 2 == 1 {
        rev1 = true;
    } else {
        rev2 = true;
    }
    if sign < 0.0 {
        rev1 = !rev1;
        rev2 = !rev2;
    }

    // `if (rev1) TgF.Reverse(); if (rev2) TgL.Reverse();`
    // (`cxx:1387-1393`).
    if rev1 {
        tgf.reverse();
    }
    if rev2 {
        tgl.reverse();
    }

    Ok((ptgui, tgf, tgl))
}

/// The shared tail of the three arms: the two `Extrema_GenLocateExtPS`
/// projections that refresh `SolDep` (`cxx:1404-1416`, `:1518-1531`,
/// `:1599-1611`) and `BRepBlend_Walking::PerformFirstSection` (`cxx:1418`,
/// `:1533`, `:1613`).
#[allow(clippy::too_many_arguments)]
fn project_and_walk(
    func: &mut dyn BlendFunction,
    s1_surf: &Arc<dyn Surface>,
    s2_surf: &Arc<dyn Surface>,
    i1: &Rc<RefCell<TopolTool>>,
    i2: &Rc<RefCell<TopolTool>>,
    guide: &Arc<ElSpine>,
    par: f64,
    sol_dep: &mut MathVector,
    pos1: &mut State,
    pos2: &mut State,
    ctx: &ChFi3dCtx<'_>,
    pt1: GpPnt,
    pt2: GpPnt,
    tol_guide: f64,
) -> Result<bool, String> {
    // `double tol = tolesp * 1.e2;` (`cxx:1403`, `:1518`, `:1599`).
    let tol = ctx.tolesp * 1.0e2;

    // `Extrema_GenLocateExtPS proj1(*S1, tol, tol); proj1.Perform(pt1,
    // SolDep(1), SolDep(2));` (`cxx:1404-1407`). The three-argument `Perform`
    // is the `isDistanceCriteria = false` default of `hxx:47-49`.
    let mut proj1 = GenLocateExtPs::new(s1_surf.as_ref(), tol, tol);
    proj1.perform_default(&pt1, sol_dep.value(1), sol_dep.value(2));
    let mut proj2 = GenLocateExtPs::new(s2_surf.as_ref(), tol, tol);
    proj2.perform_default(&pt2, sol_dep.value(3), sol_dep.value(4));

    // `if (proj1.IsDone()) (proj1.Point()).Parameter(SolDep(1), SolDep(2));`
    // (`cxx:1409-1416`).
    if let Some(point) = proj1.point() {
        let (u, v) = point.parameter();
        sol_dep.set_value(1, u);
        sol_dep.set_value(2, v);
    }
    if let Some(point) = proj2.point() {
        let (u, v) = point.parameter();
        sol_dep.set_value(3, u);
        sol_dep.set_value(4, v);
    }

    // `BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide); return
    // TheWalk.PerformFirstSection(*pFunc, Par, SolDep, tolapp3d, TolGuide,
    // Pos1, Pos2);` (`cxx:1357` + `:1418`, `:1453` + `:1533`, `:1544` +
    // `:1613`).
    let mut walk = Walking::new(
        Arc::clone(s1_surf),
        Arc::clone(s2_surf),
        Rc::clone(i1),
        Rc::clone(i2),
        Arc::clone(guide),
    );
    let (mut f1, mut f2) = (state_face(*pos1), state_face(*pos2));
    let done = walk.perform_first_section(
        func,
        par,
        sol_dep,
        ctx.tolapp3d,
        tol_guide,
        &mut f1,
        &mut f2,
    );
    *pos1 = state_of(f1);
    *pos2 = state_of(f2);
    Ok(done)
}
