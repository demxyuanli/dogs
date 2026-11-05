//! `ChFi3d_FilBuilder::SimulSurf` — the face / face overload.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:558-787`
//! (OCCT 8.0.0). It simulates the band between two faces along one guideline
//! and fills the `SurfData` with the circles of the section (`Simul`), the four
//! 2d points of the extremities and the four `ChFiDS_CommonPoint`s.
//!
//! The body is split in two by `fsp->IsConstant()` (`cxx:605`): the constant
//! radius arm runs `BRepBlend_ConstRad` / `BRepBlend_ConstRadInv`
//! (`BlendFunc_ConstRad.cxx:86-118`, `BlendFunc_ConstRadInv.cxx:60-...`), the
//! variable one `BRepBlend_EvolRad` / `BRepBlend_EvolRadInv` with
//! `fsp->Law(HGuide)` (`cxx:666-667`). Both arms then share the same section
//! loop (`cxx:637-661` / `:697-721`), the `SetSimul` / `Set2dPoints` writes,
//! the four `ChFi3d_FilCommonPoint` calls and the two `SearchFace` retries
//! (`cxx:722-783`).
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError` on a spine that is not a
//!   `ChFiDS_FilSpine` (`cxx:579-582`) is the `Err` of
//!   [`simul_surf_face_face`]; the trait entry point answers OCCT's `false`,
//! - `handle` sharing: OCCT mutates `HGuide` in place, so the blend function,
//!   the `Walking` and the caller all see the reframed parameters. The port's
//!   `Arc<ElSpine>` can only be reframed through `Arc::make_mut`, which clones
//!   once the blend function holds a second reference — the same pre-existing
//!   divergence as [`crate::brep_blend::walking`] (its `hguide` field is an
//!   `Arc` too),
//! - `gp_Pnt2d pf1, pl1, pf2, pl2` are default-constructed in OCCT and only
//!   written on the first / last point (which always exists: `SimulData`
//!   answers `false` for an empty line), so the port starts them at the origin,
//! - `occ::handle<NCollection_HArray1<ChFiDS_CircSection>>` is
//!   [`CircSectionArray`] and `isec.Set(ci, p1, p2)` is `set_circ`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom::law::{LawFunction, LawHandle};
use occt_math::MathVector;

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{ConstRad, ConstRadInv, EvolRad, EvolRadInv};
use crate::brep_blend::line::Line;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::simul_data_face_face::simul_data_face_face;
use crate::chfi3d::builder::ChFi3dCtx;
use crate::chfi3d::builder_helpers::fil_common_point;
use crate::chfi_ds::{CircSectionArray, ElSpine, SpineClass, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::statics::simul_params;

/// `HGuide` as the `GeomAdaptor_Curve` the blend functions take.
///
/// OCCT passes the one `handle(ChFiDS_ElSpine)` to `BRepBlend_ConstRad` /
/// `BRepBlend_EvolRad` and to `BRepBlend_Walking` alike (`cxx:607-608`,
/// `:666-667`); the port hands out clones of the one `Arc<ElSpine>`.
pub(crate) fn guide_as_curve(guide: &Arc<ElSpine>) -> Arc<dyn Curve> {
    let arc: Arc<ElSpine> = Arc::clone(guide);
    arc
}

/// `SimulSurf(Data, HGuide, Spine, Choix, S1, I1, S2, I2, TolGuide, First,
/// Last, Inside, Appro, Forward, RecOnS1, RecOnS2, Soldep, intf, intl)`
/// (`cxx:558-787`).
///
/// `shape` is the subclass member `myShape` (`hxx:377`), `ctx` the base members
/// the body reads (`tolapp3d`, `myVEMap`, `myEFMap`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn simul_surf_face_face(
    shape: BlendFuncSectionShape,
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
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("SimulSurf : this is not the spine of the
    // fillet");` (`cxx:578-582`).
    let spine_handle = spine
        .ok_or_else(|| "SimulSurf : this is not the spine of the fillet".to_string())?;
    let is_fil = {
        let guard = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?;
        guard.class() == SpineClass::Fil
    };
    if !is_fil {
        return Err("SimulSurf : this is not the spine of the fillet".to_string());
    }

    let mut lin: Option<Rc<RefCell<Line>>> = None;

    // `SimulParams(HGuide, fsp, MaxStep, locfleche);` (`cxx:588`).
    let (max_step, locfleche) = simul_params(guide, spine_handle)?;

    // `gp_Pnt2d pf1, pl1, pf2, pl2;` (`cxx:591`).
    let mut pf1 = GpPnt2d::new(0.0, 0.0);
    let mut pl1 = GpPnt2d::new(0.0, 0.0);
    let mut pf2 = GpPnt2d::new(0.0, 0.0);
    let mut pl2 = GpPnt2d::new(0.0, 0.0);

    // `double PFirst = First;` (`cxx:595`) — before the two reframings below.
    let p_first = *first;
    // `if (intf) { First = fsp->FirstParameter(1); }` (`cxx:596-599`).
    if *intf != 0 {
        *first = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .first_parameter_of(1)
            .ok_or_else(|| "SimulSurf: ChFiDS_Spine::FirstParameter(1) raised".to_string())?;
    }
    // `if (intl) { Last = fsp->LastParameter(fsp->NbEdges()); }` (`cxx:600-603`).
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
    // `SimulData`. OCCT passes the `BRepAdaptor_Surface` handles themselves.
    let s1_surf = s1
        .surface()
        .ok_or_else(|| "SimulSurf: S1 has no surface".to_string())?
        .clone();
    let s2_surf = s2
        .surface()
        .ok_or_else(|| "SimulSurf: S2 has no surface".to_string())?
        .clone();
    // `BRepBlend_ConstRad Func(S1, S2, HGuide)` / `BRepBlend_EvolRad ...` hold
    // the same `HGuide` handle as the walker: OCCT shares one `ChFiDS_ElSpine`,
    // the port shares the one `Arc<ElSpine>`.
    let guide_curve: Arc<dyn Curve> = guide_as_curve(guide);

    // `if (fsp->IsConstant())` (`cxx:605`).
    let constant = spine_handle
        .read()
        .map_err(|_| "SimulSurf: poisoned spine".to_string())?
        .is_fil_constant();
    let sec = if constant {
        // `BRepBlend_ConstRad Func(S1, S2, HGuide);
        // BRepBlend_ConstRadInv FInv(S1, S2, HGuide);` (`cxx:607-608`).
        let radius = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_radius()
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string()
            })?;
        let mut func = ConstRad::new(Arc::clone(&s1_surf), Arc::clone(&s2_surf), guide_curve);
        let mut finv = ConstRadInv::new(Arc::clone(&s1_surf), Arc::clone(&s2_surf), guide_as_curve(guide));
        // `Func.Set(fsp->Radius(), Choix); FInv.Set(fsp->Radius(), Choix);
        // Func.Set(myShape);` (`cxx:609-610`).
        func.set_radius(radius, choix);
        finv.set_radius(radius, choix);
        func.set_section_shape(shape);

        // `done = SimulData(Data, HGuide, EmptyHGuide, lin, S1, I1, S2, I2,
        // Func, FInv, PFirst, MaxStep, locfleche, TolGuide, First, Last,
        // Inside, Appro, Forward, Soldep, 4, RecOnS1, RecOnS2);`
        // (`cxx:611-635`).
        let done = simul_data_face_face(
            ctx.tolapp3d,
            guide,
            None,
            &mut lin,
            Arc::clone(&s1_surf),
            Rc::clone(i1),
            Arc::clone(&s2_surf),
            Rc::clone(i2),
            &mut func,
            &mut finv,
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
        );
        if !done {
            return Ok(false);
        }
        section_from_line(&lin, &mut func, &mut pf1, &mut pl1, &mut pf2, &mut pl2)?
    } else {
        // `BRepBlend_EvolRad Func(S1, S2, HGuide, fsp->Law(HGuide));
        // BRepBlend_EvolRadInv FInv(S1, S2, HGuide, fsp->Law(HGuide));`
        // (`cxx:666-667`).
        let law = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_law(guide)
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Law raised Standard_DomainError".to_string()
            })?;
        // `fsp->Law(HGuide)` is one `handle(Law_Composite)` shared by the two
        // blend functions and by the spine; [`LawHandle`] is that handle.
        let law_f: Arc<dyn LawFunction> = Arc::new(LawHandle::new(Arc::clone(&law)));
        let law_inv: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
        let mut func = EvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Arc::clone(&guide_curve),
            law_f,
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            guide_as_curve(guide),
            law_inv,
        );
        // `Func.Set(Choix); FInv.Set(Choix); Func.Set(myShape);`
        // (`cxx:668-669`).
        func.set_choix(choix);
        finv.set_choix(choix);
        func.set_section_shape(shape);

        // `done = SimulData(...);` (`cxx:670-694`), same argument list as the
        // constant arm.
        let done = simul_data_face_face(
            ctx.tolapp3d,
            guide,
            None,
            &mut lin,
            Arc::clone(&s1_surf),
            Rc::clone(i1),
            Arc::clone(&s2_surf),
            Rc::clone(i2),
            &mut func,
            &mut finv,
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
        );
        if !done {
            return Ok(false);
        }
        section_from_line_evol(&lin, &mut func, &mut pf1, &mut pl1, &mut pf2, &mut pl2)?
    };

    // `Data->SetSimul(sec);` (`cxx:722`).
    data.set_simul(Arc::new(sec));
    // `Data->Set2dPoints(pf1, pl1, pf2, pl2);` (`cxx:723`).
    data.set_2d_points(pf1, pl1, pf2, pl2);

    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let line = line.borrow();

    // The four `ChFi3d_FilCommonPoint(...)` calls (`cxx:724-741`).
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

    // `bool reverse = (!Forward || Inside);` (`cxx:742`).
    let reverse = !forward || inside;
    if *intf != 0 && reverse {
        // `bool ok = false; const ChFiDS_CommonPoint& cp1 =
        // Data->VertexFirstOnS1(); if (cp1.IsOnArc()) { TopoDS_Face F1 =
        // S1->Face(); TopoDS_Face bid; intf = !SearchFace(Spine, cp1, F1, bid);
        // ok = intf != 0; }` (`cxx:743-756`).
        let mut ok = false;
        let on_arc1 = data.vertex_first_on_s1().is_on_arc();
        if on_arc1 {
            let f1 = s1.face().cloned();
            let mut bid = None;
            let cp1 = data.vertex_first_on_s1().clone();
            let found = ctx.search_face(Some(spine_handle), &cp1, f1.as_ref(), &mut bid);
            *intf = i32::from(!found);
            ok = *intf != 0;
        }
        // The second arm has the extra `&& !ok` guard (`cxx:757-761`).
        if data.vertex_first_on_s2().is_on_arc() && !ok {
            let f2 = s2.face().cloned();
            let mut bid = None;
            let cp2 = data.vertex_first_on_s2().clone();
            let found = ctx.search_face(Some(spine_handle), &cp2, f2.as_ref(), &mut bid);
            *intf = i32::from(!found);
        }
    }

    // `if (intl) { ... }` (`cxx:763-782`): the same two arms on the last
    // common points, with `intl` in place of `intf`.
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

    Ok(true)
}

/// The section loop of the constant-radius arm (`cxx:637-661`): one
/// `ChFiDS_CircSection` per blending point, plus the 2d parameters of the first
/// and last point.
fn section_from_line(
    lin: &Option<Rc<RefCell<Line>>>,
    func: &mut ConstRad,
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
    // `int nbp = lin->NbPoints(); sec =
    // new NCollection_HArray1<ChFiDS_CircSection>(1, nbp);` (`cxx:637-639`).
    let nbp = line.nb_points();
    let mut sec = CircSectionArray::new(1, nbp as i32);
    for i in 1..=nbp {
        // `const Blend_Point& p = lin->Point(i);
        // p.ParametersOnS1(u1, v1); p.ParametersOnS2(u2, v2); w = p.Parameter();`
        // (`cxx:645-648`).
        let p = line
            .point(i)
            .ok_or_else(|| format!("SimulSurf: Lin->Point({i}) out of range"))?;
        let (u1, v1) = p.parameters_on_s1()?;
        let (u2, v2) = p.parameters_on_s2()?;
        let w = p.parameter();
        // `Func.Section(w, u1, v1, u2, v2, p1, p2, ci); isec.Set(ci, p1, p2);`
        // (`cxx:649-650`).
        let (ci, p1, p2) = func.section_circ(w, u1, v1, u2, v2)?;
        sec.change_value(i as i32)
            .ok_or_else(|| format!("SimulSurf: sec->ChangeValue({i}) out of range"))?
            .set_circ(&ci, p1, p2);
        write_extremity(i, nbp, u1, v1, u2, v2, pf1, pl1, pf2, pl2);
    }
    Ok(sec)
}

/// The section loop of the variable-radius arm (`cxx:697-721`): identical to
/// [`section_from_line`] with `BlendFunc_EvolRad::Section`.
fn section_from_line_evol(
    lin: &Option<Rc<RefCell<Line>>>,
    func: &mut EvolRad,
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
    let nbp = line.nb_points();
    let mut sec = CircSectionArray::new(1, nbp as i32);
    for i in 1..=nbp {
        let p = line
            .point(i)
            .ok_or_else(|| format!("SimulSurf: Lin->Point({i}) out of range"))?;
        let (u1, v1) = p.parameters_on_s1()?;
        let (u2, v2) = p.parameters_on_s2()?;
        let w = p.parameter();
        let (ci, p1, p2) = func.section_circ(w, u1, v1, u2, v2)?;
        sec.change_value(i as i32)
            .ok_or_else(|| format!("SimulSurf: sec->ChangeValue({i}) out of range"))?
            .set_circ(&ci, p1, p2);
        write_extremity(i, nbp, u1, v1, u2, v2, pf1, pl1, pf2, pl2);
    }
    Ok(sec)
}

/// `if (i == 1) { pf1.SetCoord(u1, v1); pf2.SetCoord(u2, v2); }
/// if (i == nbp) { pl1.SetCoord(u1, v1); pl2.SetCoord(u2, v2); }`
/// (`cxx:651-660`, `:710-719`).
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
