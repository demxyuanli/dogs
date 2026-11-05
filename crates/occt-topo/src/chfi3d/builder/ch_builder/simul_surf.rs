//! `ChFi3d_ChBuilder::SimulSurf` — the face / face overload.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:717-1222`
//! (OCCT 8.0.0). It simulates the chamfer band between two faces along one
//! guideline and fills the `SurfData` with the lines of the section (`Simul`),
//! the four 2d points of the extremities and the four
//! `ChFiDS_CommonPoint`s.
//!
//! The body is split in three by `chsp->IsChamfer()` (`cxx:781`, `:921`,
//! `:1084`):
//!
//! | arm | blend function | `Set` |
//! |---|---|---|
//! | `ChFiDS_Sym` | `Chamfer` / `ChamfInv` (`cxx:791-792`) | `(dis, dis)` (`:800-801`) |
//! | | `ConstThroat` / `ConstThroatInv` (`cxx:797-798`) | `(dis, dis)` (`:800-801`) |
//! | `ChFiDS_TwoDist` | `Chamfer` / `ChamfInv` (`cxx:935-936`) | `(dis1, dis2)` (`:937-938`) |
//! | | `ConstThroatWithPenetration` / `...Inv` (`cxx:958-959`) | `(max, max)` (`:960-962`) |
//! | `ChFiDS_DistAngle` | `ChAsym` / `ChAsymInv` (`cxx:1093-1094`) | `(dis, angle)` (`:1096-1097`) |
//!
//! and by `chsp->Mode()` (`cxx:788`, `:932`): the classic mode uses the
//! `Chamfer` family on `HGuide`, the two constant-throat modes the throat
//! family. The `TwoDist` / non-classic arm is the only one that substitutes the
//! **offset** guideline for `HGuide` (`cxx:942-959`); every other arm hands
//! `HGuide` and a null `AdditionalHGuide` to `SimulData`.
//!
//! All three arms then share the same tail (`cxx:832-920`, `:994-1082`,
//! `:1129-1218`): the section loop over the blending line, the `SetSimul` /
//! `Set2dPoints` writes, the four `ChFi3d_FilCommonPoint` calls and the two
//! `SearchFace` retries. The port keeps one copy of that tail; the section loop
//! is [`section_from_line`].
//!
//! The `radiusspine` scan (`cxx:750-766`) is `max |D2|^2` over 21 samples of
//! the guideline, despite the name: OCCT compares `d2.SquareMagnitude()`
//! against it (`cxx:761-765`) and uses it as a length in the `max` of
//! `cxx:785`, `:925-927` and `:1088-1089` — kept verbatim.
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("SimulSurf : this is not the spine of a
//!   chamfer")` (`cxx:739-743`) is the `Err` of [`simul_surf_face_face`]; the
//!   trait entry point answers OCCT's `false`, which is what
//!   `PerformSetOfSurf` acts on,
//! - `handle` sharing: OCCT mutates `HGuide` (and the offset guideline of the
//!   `TwoDist` / non-classic arm) in place, so the blend functions, the
//!   `Walking` and the caller all see the reframed parameters. The port's
//!   `Arc<ElSpine>` can only be reframed through `Arc::make_mut`, which clones
//!   once a second reference is held — the same pre-existing divergence as
//!   [`crate::brep_blend::walking`]. The offset guideline is looked up and
//!   cloned out of `ChFiDS_Spine::offset_el_spines` (`cxx:942-951`), so its
//!   reframing stays in the local handle,
//! - the "Construction of offset guide failed!" print of `cxx:954-957` has its
//!   exception commented out and would dereference a null handle in
//!   `BRepBlend_ConstThroatWithPenetration`; the port answers `Err` with that
//!   message instead,
//! - `gp_Pnt2d pf1, pl1, pf2, pl2` are default-constructed in OCCT and only
//!   written on the first / last point (which always exists: `SimulData`
//!   answers `false` for an empty line), so the port starts them at the origin,
//! - `gp_Lin line; pFunc->Section(ww, u1, v1, u2, v2, p1, p2, line);
//!   isec.Set(line, p1, p2);` (`cxx:844-846`) is
//!   [`crate::blend::chamfer::gen_chamfer_section`] plus
//!   [`CircSection::set_lin`](crate::chfi_ds::CircSection::set_lin): the body
//!   reads only the two surfaces, so the four blend functions share it,
//! - `occ::handle<NCollection_HArray1<ChFiDS_CircSection>>` is
//!   [`CircSectionArray`].

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::MathVector;

use crate::blend::chamfer::gen_chamfer_section;
use crate::blend::{
    BlendFuncInv, BlendFunction, ChAsym, ChAsymInv, ChamfInv, Chamfer, ConstThroat, ConstThroatInv,
    ConstThroatWithPenetration, ConstThroatWithPenetrationInv,
};
use crate::brep_blend::line::Line;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::simul_data_face_face::simul_data_face_face;
use crate::chfi3d::builder::ChFi3dCtx;
use crate::chfi3d::builder_helpers::fil_common_point;
use crate::chfi_ds::{
    ChamfMethod, ChamfMode, CircSectionArray, ElSpine, SpineClass, SpineHandle, SurfData,
};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `HGuide` as the `GeomAdaptor_Curve` the blend functions take.
///
/// OCCT passes the one `handle(ChFiDS_ElSpine)` to the blend functions and to
/// `BRepBlend_Walking` alike (`cxx:791-798`, `:935-959`, `:1093-1094`); the
/// port hands out clones of the one `Arc<ElSpine>`.
pub(super) fn guide_as_curve(guide: &Arc<ElSpine>) -> Arc<dyn Curve> {
    let arc: Arc<ElSpine> = Arc::clone(guide);
    arc
}

/// `SimulSurf(Data, HGuide, Spine, Choix, S1, I1, S2, I2, TolGuide, First,
/// Last, Inside, Appro, Forward, RecOnS1, RecOnS2, Soldep, intf, intl)`
/// (`cxx:717-1222`).
///
/// `ctx` is the base members the body reads (`tolapp3d`, `myVEMap`,
/// `myEFMap`); OCCT reads them off `this` because the method is a member of the
/// same object.
#[allow(clippy::too_many_arguments)]
pub(crate) fn simul_surf_face_face(
    data: &mut SurfData,
    guide: &mut Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    choix: i32,
    s1: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    s2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    tol_guide: f64,
    first: &mut f64,
    last: &mut f64,
    inside: bool,
    approx: bool,
    forward: bool,
    rec_on_s1: bool,
    rec_on_s2: bool,
    soldep: &MathVector,
    intf: &mut i32,
    intl: &mut i32,
    ctx: &ChFi3dCtx<'_>,
) -> Result<bool, String> {
    // `chsp = down_cast<ChFiDS_ChamfSpine>(Spine); if (chsp.IsNull()) throw
    // Standard_ConstructionError("SimulSurf : this is not the spine of a
    // chamfer");` (`cxx:739-743`).
    let spine_handle = spine
        .ok_or_else(|| "SimulSurf : this is not the spine of a chamfer".to_string())?;
    let is_chamf = {
        let guard = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?;
        guard.class() == SpineClass::Chamf
    };
    if !is_chamf {
        return Err("SimulSurf : this is not the spine of a chamfer".to_string());
    }

    // `double la = HGuide->LastParameter(), fi = HGuide->FirstParameter();
    // double longueur = la - fi; double MaxStep = longueur * 0.05;`
    // (`cxx:747-750`).
    let fi = guide.first_parameter();
    let la = guide.last_parameter();
    let longueur = la - fi;
    let max_step = longueur * 0.05;

    // `for (i = 0; i <= 20; i++) { w = fi + i * MaxStep; HGuide->D2(w, ...);
    // double temp = d2.SquareMagnitude(); if (temp > radiusspine) radiusspine
    // = temp; }` (`cxx:755-766`).
    //
    // `ChFiDS_ElSpine::D2` raises `Standard_NullObject` on a spine without a
    // curve; the port answers that as `Err`.
    let mut radiusspine = 0.0f64;
    for i in 0..=20 {
        let w = fi + f64::from(i) * max_step;
        let (_p, _d1, d2) = guide.eval_d2(w).ok_or_else(|| {
            "SimulSurf: ChFiDS_ElSpine::D2 raised Standard_NullObject".to_string()
        })?;
        let temp = d2.square_magnitude();
        if temp > radiusspine {
            radiusspine = temp;
        }
    }

    // `double PFirst = First; if (intf) { First = chsp->FirstParameter(1); }
    // if (intl) { Last = chsp->LastParameter(chsp->NbEdges()); }`
    // (`cxx:769-778`).
    let mut lin: Option<Rc<RefCell<Line>>> = None;
    let p_first = *first;
    if *intf != 0 {
        *first = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .first_parameter_of(1)
            .ok_or_else(|| "SimulSurf: ChFiDS_Spine::FirstParameter(1) raised".to_string())?;
    }
    if *intl != 0 {
        let guard = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?;
        let nb_edges = guard.nb_edges() as i32;
        *last = guard
            .last_parameter_of(nb_edges)
            .ok_or_else(|| "SimulSurf: ChFiDS_Spine::LastParameter raised".to_string())?;
    }

    // The two adaptor surfaces, shared by the blend functions and by
    // `SimulData`; OCCT passes the `BRepAdaptor_Surface` handles themselves.
    let s1_surf = s1
        .surface()
        .ok_or_else(|| "SimulSurf: S1 has no surface".to_string())?
        .clone();
    let s2_surf = s2
        .surface()
        .ok_or_else(|| "SimulSurf: S2 has no surface".to_string())?
        .clone();
    let guide_curve = guide_as_curve(guide);

    // `if (chsp->IsChamfer() == ...)` (`cxx:781`) and `chsp->Mode()`
    // (`cxx:788`, `:932`, `:1084`): the two reads happen before the arms
    // because each of them borrows the spine for its own `SimulData`.
    let (method, mode, chamf_dis1, chamf_dis2, chamf_angle) = {
        let guard = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?;
        let method = guard.is_chamfer();
        let mode = guard.mode();
        // `chsp->GetDist(dis)` / `Dists(dis1, dis2)` / `GetDistAngle(dis,
        // angle)` raise `Standard_Failure` on a spine of another kind; the
        // `IsChamfer()` arm already told which one applies (`cxx:783-784`,
        // `:923-924`, `:1085-1086`).
        let (d1, d2, angle) = match method {
            ChamfMethod::Sym => (
                guard
                    .get_dist()
                    .ok_or_else(|| "SimulSurf: ChFiDS_ChamfSpine::GetDist raised".to_string())?,
                0.0,
                0.0,
            ),
            ChamfMethod::TwoDist => {
                let (d1, d2) = guard
                    .dists()
                    .ok_or_else(|| "SimulSurf: ChFiDS_ChamfSpine::Dists raised".to_string())?;
                (d1, d2, 0.0)
            }
            ChamfMethod::DistAngle => {
                let (d1, angle) = guard.get_dist_angle().ok_or_else(|| {
                    "SimulSurf: ChFiDS_ChamfSpine::GetDistAngle raised".to_string()
                })?;
                (d1, 0.0, angle)
            }
        };
        (method, mode, d1, d2, angle)
    };

    // The three arms of `cxx:781-1218`; each answers `done` and shares the
    // `lin` output.
    let done = match method {
        ChamfMethod::Sym => {
            // `chsp->GetDist(dis); radius = std::max(dis, radiusspine);
            // locfleche = radius * 1.e-2;` (`cxx:783-786`).
            let dis = chamf_dis1;
            let radius = dis.max(radiusspine);
            let locfleche = radius * 1.0e-2;
            // `if (chsp->Mode() == ChFiDS_ClassicChamfer) { Chamfer;
            // ChamfInv; } else { ConstThroat; ConstThroatInv; }`
            // (`cxx:789-799`), then `Set(dis, dis, Choix)` on both
            // (`cxx:800-801`).
            let (mut func, mut finv): (Box<dyn BlendFunction>, Box<dyn BlendFuncInv>) =
                if mode == ChamfMode::ClassicChamfer {
                    let mut f = Chamfer::new(
                        Arc::clone(&s1_surf),
                        Arc::clone(&s2_surf),
                        Arc::clone(&guide_curve),
                    );
                    let mut f_inv = ChamfInv::new(
                        Arc::clone(&s1_surf),
                        Arc::clone(&s2_surf),
                        guide_as_curve(guide),
                    );
                    f.set_distances(dis, dis, choix);
                    f_inv.set_distances(dis, dis, choix);
                    (Box::new(f), Box::new(f_inv))
                } else {
                    let mut f = ConstThroat::new(
                        Arc::clone(&s1_surf),
                        Arc::clone(&s2_surf),
                        Arc::clone(&guide_curve),
                    );
                    let mut f_inv = ConstThroatInv::new(
                        Arc::clone(&s1_surf),
                        Arc::clone(&s2_surf),
                        guide_as_curve(guide),
                    );
                    f.set_throat(dis, dis, choix);
                    f_inv.set_throat(dis, dis, choix);
                    (Box::new(f), Box::new(f_inv))
                };
            // `done = SimulData(Data, HGuide, OffsetHGuide, ...)` — the null
            // `OffsetHGuide` of this arm (`cxx:803-831`).
            simul_data_face_face(
                ctx.tolapp3d,
                guide,
                None,
                &mut lin,
                Arc::clone(&s1_surf),
                Rc::clone(i1),
                Arc::clone(&s2_surf),
                Rc::clone(i2),
                &mut *func,
                &mut *finv,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                approx,
                forward,
                soldep,
                4,
                rec_on_s1,
                rec_on_s2,
            )
        }
        ChamfMethod::TwoDist => {
            // `chsp->Dists(dis1, dis2); radius = max(dis1, dis2); radius =
            // max(radius, radiusspine); locfleche = radius * 1.e-2;`
            // (`cxx:923-928`).
            let (dis1, dis2) = (chamf_dis1, chamf_dis2);
            let radius = dis1.max(dis2).max(radiusspine);
            let locfleche = radius * 1.0e-2;
            if mode == ChamfMode::ClassicChamfer {
                // `Chamfer; ChamfInv; Set(dis1, dis2, Choix)` on both
                // (`cxx:933-938`).
                let mut f = Chamfer::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    Arc::clone(&guide_curve),
                );
                let mut f_inv = ChamfInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                f.set_distances(dis1, dis2, choix);
                f_inv.set_distances(dis1, dis2, choix);
                simul_data_face_face(
                    ctx.tolapp3d,
                    guide,
                    None,
                    &mut lin,
                    Arc::clone(&s1_surf),
                    Rc::clone(i1),
                    Arc::clone(&s2_surf),
                    Rc::clone(i2),
                    &mut f,
                    &mut f_inv,
                    p_first,
                    max_step,
                    locfleche,
                    tol_guide,
                    first,
                    last,
                    inside,
                    approx,
                    forward,
                    soldep,
                    4,
                    rec_on_s1,
                    rec_on_s2,
                )
            } else {
                // `NCollection_List<handle(ChFiDS_ElSpine)>& ll =
                // Spine->ChangeElSpines(); ... for (; ILES.More();
                // ILES.Next(), ILES_offset.Next()) { if (aHElSpine ==
                // HGuide) OffsetHGuide = ILES_offset.Value(); }`
                // (`cxx:942-951`).
                //
                // `if (OffsetHGuide.IsNull()) { std::cout << "Construction of
                // offset guide failed!"; }` (`cxx:953-957`): the exception is
                // commented out, and the null handle is dereferenced right
                // after; the port answers that message instead.
                let offset_arc: Arc<ElSpine> = {
                    let sg = spine_handle
                        .read()
                        .map_err(|_| "SimulSurf: poisoned spine".to_string())?;
                    let offset = sg.offset_el_spine(guide).ok_or_else(|| {
                        "Construction of offset guide failed!".to_string()
                    })?;
                    Arc::clone(offset)
                };
                let offset_curve = guide_as_curve(&offset_arc);
                // `ConstThroatWithPenetration(S1, S2, OffsetHGuide);
                // ConstThroatWithPenetrationInv(S1, S2, OffsetHGuide);
                // double Throat = std::max(dis1, dis2); Set(Throat, Throat,
                // Choix);` (`cxx:958-962`).
                let mut f = ConstThroatWithPenetration::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    offset_curve,
                );
                let mut f_inv = ConstThroatWithPenetrationInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    offset_arc.clone(),
                );
                let throat = dis1.max(dis2);
                f.set_throat(throat, throat, choix);
                f_inv.set_throat(throat, throat, choix);
                // `SimulData(Data, HGuide, OffsetHGuide, ...)`
                // (`cxx:964-991`). The reframed parameters stay in the local
                // handle: OCCT writes through the shared one — see the module
                // note.
                let mut offset = offset_arc;
                simul_data_face_face(
                    ctx.tolapp3d,
                    guide,
                    Some(&mut offset),
                    &mut lin,
                    Arc::clone(&s1_surf),
                    Rc::clone(i1),
                    Arc::clone(&s2_surf),
                    Rc::clone(i2),
                    &mut f,
                    &mut f_inv,
                    p_first,
                    max_step,
                    locfleche,
                    tol_guide,
                    first,
                    last,
                    inside,
                    approx,
                    forward,
                    soldep,
                    4,
                    rec_on_s1,
                    rec_on_s2,
                )
            }
        }
        ChamfMethod::DistAngle => {
            // `chsp->GetDistAngle(dis, angle); radius = max(dis, dis *
            // tan(angle)); radius = max(radius, radiusspine); locfleche =
            // radius * 1.e-2;` (`cxx:1085-1090`).
            let dis = chamf_dis1;
            let angle = chamf_angle;
            let radius = dis.max(dis * angle.tan()).max(radiusspine);
            let locfleche = radius * 1.0e-2;
            // `int Ch = Choix; BRepBlend_ChAsym Func(S1, S2, HGuide);
            // BRepBlend_ChAsymInv FInv(S1, S2, HGuide); Func.Set(dis, angle,
            // Ch); FInv.Set(dis, angle, Ch);` (`cxx:1092-1098`).
            let mut f = ChAsym::new(
                Arc::clone(&s1_surf),
                Arc::clone(&s2_surf),
                Arc::clone(&guide_curve),
            );
            let mut f_inv = ChAsymInv::new(
                Arc::clone(&s1_surf),
                Arc::clone(&s2_surf),
                guide_as_curve(guide),
            );
            f.set_dist_angle(dis, angle, choix);
            f_inv.set_dist_angle(dis, angle, choix);
            // `done = SimulData(Data, HGuide, OffsetHGuide, ...)`
            // (`cxx:1100-1127`).
            simul_data_face_face(
                ctx.tolapp3d,
                guide,
                None,
                &mut lin,
                Arc::clone(&s1_surf),
                Rc::clone(i1),
                Arc::clone(&s2_surf),
                Rc::clone(i2),
                &mut f,
                &mut f_inv,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                approx,
                forward,
                soldep,
                4,
                rec_on_s1,
                rec_on_s2,
            )
        }
    };

    // `if (!done) { return false; }` (`cxx:828-831`, `:989-992`,
    // `:1124-1127`).
    if !done {
        return Ok(false);
    }

    // `gp_Pnt2d pf1, pl1, pf2, pl2;` (`cxx:834`) — written by the section loop.
    let mut pf1 = GpPnt2d::new(0.0, 0.0);
    let mut pl1 = GpPnt2d::new(0.0, 0.0);
    let mut pf2 = GpPnt2d::new(0.0, 0.0);
    let mut pl2 = GpPnt2d::new(0.0, 0.0);
    // `int nbp = lin->NbPoints(); sec = new NCollection_HArray1<...>(1, nbp);
    // for (i = 1; i <= nbp; i++) { ... pFunc->Section(ww, u1, v1, u2, v2, p1,
    // p2, line); isec.Set(line, p1, p2); ... }` (`cxx:836-857`, `:996-1017`,
    // `:1131-1152`).
    let sec = section_from_line(&lin, &*s1_surf, &*s2_surf, &mut pf1, &mut pl1, &mut pf2, &mut pl2)?;

    // `Data->SetSimul(sec); Data->Set2dPoints(pf1, pl1, pf2, pl2);`
    // (`cxx:859-860`, `:1021-1022`, `:1156-1157`).
    data.set_simul(Arc::new(sec));
    data.set_2d_points(pf1, pl1, pf2, pl2);

    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let line = line.borrow();

    // The four `ChFi3d_FilCommonPoint(...)` calls (`cxx:861-880`).
    let trans1 = line.transition_on_s1()?;
    let trans2 = line.transition_on_s2()?;
    fil_common_point(
        line.start_point_on_first(),
        trans1,
        true,
        data.vertex_first_on_s1_mut(),
        ctx.tolapp3d,
    );
    fil_common_point(
        line.end_point_on_first(),
        trans1,
        false,
        data.vertex_last_on_s1_mut(),
        ctx.tolapp3d,
    );
    fil_common_point(
        line.start_point_on_second(),
        trans2,
        true,
        data.vertex_first_on_s2_mut(),
        ctx.tolapp3d,
    );
    fil_common_point(
        line.end_point_on_second(),
        trans2,
        false,
        data.vertex_last_on_s2_mut(),
        ctx.tolapp3d,
    );
    drop(line);

    // `bool reverse = (!Forward || Inside); if (intf && reverse) { ... }`
    // (`cxx:882-901`, `:1043-1062`, `:1179-1198`).
    let reverse = !forward || inside;
    if *intf != 0 && reverse {
        let mut ok = false;
        // `const ChFiDS_CommonPoint& cp1 = Data->VertexFirstOnS1(); if
        // (cp1.IsOnArc()) { TopoDS_Face F1 = S1->Face(); TopoDS_Face bid;
        // intf = !SearchFace(Spine, cp1, F1, bid); ok = intf != 0; }`
        // (`cxx:886-893`).
        if data.vertex_first_on_s1().is_on_arc() {
            let f1 = s1.face().cloned();
            let mut bid = None;
            let cp1 = data.vertex_first_on_s1().clone();
            let found = ctx.search_face(Some(spine_handle), &cp1, f1.as_ref(), &mut bid);
            *intf = i32::from(!found);
            ok = *intf != 0;
        }
        // `const ChFiDS_CommonPoint& cp2 = Data->VertexFirstOnS2(); if
        // (cp2.IsOnArc() && !ok) { ... }` (`cxx:894-900`).
        if data.vertex_first_on_s2().is_on_arc() && !ok {
            let f2 = s2.face().cloned();
            let mut bid = None;
            let cp2 = data.vertex_first_on_s2().clone();
            let found = ctx.search_face(Some(spine_handle), &cp2, f2.as_ref(), &mut bid);
            *intf = i32::from(!found);
        }
    }

    // `if (intl) { ... }` (`cxx:903-919`, `:1063-1081`, `:1200-1217`): the
    // same two arms on the last common points, with `intl` in place of `intf`
    // and no `reverse` guard.
    if *intl != 0 {
        let mut ok = false;
        if data.vertex_last_on_s1().is_on_arc() {
            let f1 = s1.face().cloned();
            let mut bid = None;
            let cp1 = data.vertex_last_on_s1().clone();
            let found = ctx.search_face(Some(spine_handle), &cp1, f1.as_ref(), &mut bid);
            *intl = i32::from(!found);
            ok = *intl != 0;
        }
        if data.vertex_last_on_s2().is_on_arc() && !ok {
            let f2 = s2.face().cloned();
            let mut bid = None;
            let cp2 = data.vertex_last_on_s2().clone();
            let found = ctx.search_face(Some(spine_handle), &cp2, f2.as_ref(), &mut bid);
            *intl = i32::from(!found);
        }
    }

    // `return true;` (`cxx:1219`).
    Ok(true)
}

/// The section loop of `cxx:836-857` (repeated at `:996-1017`, `:1131-1152`):
/// one `ChFiDS_CircSection` per blending point, plus the 2d parameters of the
/// first and last point.
///
/// `pFunc->Section(ww, u1, v1, u2, v2, p1, p2, line)` reads only the two
/// support surfaces, so the loop takes them instead of the four blend objects
/// (`cxx:844`, `:1005`, `:1141`).
fn section_from_line(
    lin: &Option<Rc<RefCell<Line>>>,
    s1: &dyn Surface,
    s2: &dyn Surface,
    pf1: &mut GpPnt2d,
    pl1: &mut GpPnt2d,
    pf2: &mut GpPnt2d,
    pl2: &mut GpPnt2d,
) -> Result<CircSectionArray, String> {
    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let line = line.borrow();
    // `int nbp = lin->NbPoints(); sec = new
    // NCollection_HArray1<ChFiDS_CircSection>(1, nbp);` (`cxx:836-838`).
    let nbp = line.nb_points();
    let mut sec = CircSectionArray::new(1, nbp as i32);
    for i in 1..=nbp {
        // `const Blend_Point& p = lin->Point(i);
        // p.ParametersOnS1(u1, v1); p.ParametersOnS2(u2, v2); ww =
        // p.Parameter();` (`cxx:841-844`).
        let p = line
            .point(i)
            .ok_or_else(|| format!("SimulSurf: Lin->Point({i}) out of range"))?;
        let (u1, v1) = p.parameters_on_s1()?;
        let (u2, v2) = p.parameters_on_s2()?;
        // `pFunc->Section(ww, u1, v1, u2, v2, p1, p2, line);
        // isec.Set(line, p1, p2);` (`cxx:844-846`).
        let (p1, p2, geom) = gen_chamfer_section(s1, s2, u1, v1, u2, v2)?;
        sec.change_value(i as i32)
            .ok_or_else(|| format!("SimulSurf: sec->ChangeValue({i}) out of range"))?
            .set_lin(&geom, p1, p2);
        write_extremity(i, nbp, u1, v1, u2, v2, pf1, pl1, pf2, pl2);
    }
    Ok(sec)
}

/// `if (i == 1) { pf1.SetCoord(u1, v1); pf2.SetCoord(u2, v2); }
/// if (i == nbp) { pl1.SetCoord(u1, v1); pl2.SetCoord(u2, v2); }`
/// (`cxx:847-856`, `:1008-1017`, `:1144-1153`).
#[allow(clippy::too_many_arguments)]
fn write_extremity(
    i: usize,
    nbp: usize,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    pf1: &mut GpPnt2d,
    pl1: &mut GpPnt2d,
    pf2: &mut GpPnt2d,
    pl2: &mut GpPnt2d,
) {
    if i == 1 {
        *pf1 = GpPnt2d::new(u1, v1);
        *pf2 = GpPnt2d::new(u2, v2);
    }
    if i == nbp {
        *pl1 = GpPnt2d::new(u1, v1);
        *pl2 = GpPnt2d::new(u2, v2);
    }
}
