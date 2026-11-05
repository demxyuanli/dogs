//! `ChFi3d_FilBuilder::SimulSurf` — the two edge / face overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:787-1248`
//! (OCCT 8.0.0):
//!
//! - [`simul_surf_on_s1`] — `SimulSurf(Data, HGuide, Spine, Choix, HS1, I1, PC1,
//!   HSref1, PCref1, Decroch1, HS2, I2, Or2, Fleche, TolGuide, First, Last,
//!   Inside, Appro, Forward, RecP, RecS, RecRst, Soldep)` (`cxx:787-1017`), the
//!   band between a face (`HS2`) and the restriction `PC1` of `HS1`,
//! - [`simul_surf_on_s2`] — the mirror overload (`cxx:1020-1248`).
//!
//! Both walk with `BRepBlend_SurfRstLineBuilder` through
//! `ChFi3d_Builder::SimulData(SurfRst)` (`ChFi3d_Builder_6.cxx:1346-1501`) and
//! fill only `Data->SetSimul` / `Set2dPoints` / the four common points: these
//! two overloads are `void` in OCCT and throw on failure, which is the `Err`
//! here.
//!
//! Divergences that are not heuristics:
//! - `PetitChoix` starts at `1`, becomes `3` when the "other" orientation is
//!   `TopAbs_REVERSED`, then `++` when `Choix` is even (`cxx:837-845`,
//!   `:1071-1079`); `Choix % 2` of a negative `Choix` is negative in C++ as
//!   well, so the port keeps `choix % 2 == 0` verbatim,
//! - `Adaptor3d_CurveOnSurface HC; HC->Load(PC1, HS1);` (`cxx:831-832`) is
//!   [`CurveOnSurface`] over the same p-curve and surface,
//! - the `Set2dPoints` argument order and the S1 / S2 vertex pairing are
//!   **swapped** between the two overloads (`cxx:996-1015` vs `:1226-1247`);
//!   the port keeps each overload's own order,
//! - `throw Standard_Failure("SimulSurf : ...")` (`cxx:883`, `:965`, `:1114`,
//!   `:1196`) is the `Err` here,
//! - `const double /*Fleche*/` / `const double /*Arrow*/` (`cxx:801`,
//!   `:1034`) are unnamed in OCCT and unused; the trait keeps the slot.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpCirc, GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::law::{LawFunction, LawHandle};
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::func_inv::BlendFuncInv;
use crate::blend::{ConstRadInv, EvolRadInv};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::inv::{
    SurfCurvConstRadInv, SurfCurvEvolRadInv, SurfPointConstRadInv, SurfPointEvolRadInv,
};
use crate::brep_blend::{Line, SurfRstConstRad, SurfRstEvolRad};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::simul_data::simul_data_surf_rst;
use crate::chfi3d::builder::ChFi3dCtx;
use crate::chfi3d::builder_helpers::fil_common_point;
use crate::chfi_ds::{CircSectionArray, ElSpine, SpineClass, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::simul_surf::guide_as_curve;
use super::statics::simul_params;

/// `PetitChoix` of the two edge / face overloads (`cxx:837-845`,
/// `:1071-1079`) and of the two `PerformSurf` ones (`cxx:1768-1775`,
/// `:1832-1839`, `:1941-1950`, `:2006-2013`).
pub(super) fn petit_choix(choix: i32, other_or: Orientation) -> i32 {
    let mut petit = 1;
    // `if (Or2 == TopAbs_REVERSED) { petitchoix = 3; }`.
    if other_or == Orientation::Reversed {
        petit = 3;
    }
    // `if (Choix % 2 == 0) { petitchoix++; }`.
    if choix % 2 == 0 {
        petit += 1;
    }
    petit
}

/// `occ::down_cast<ChFiDS_FilSpine>(Spine)` followed by the null check of the
/// two bodies (`cxx:815-818`, `:1045-1049`).
fn fil_spine<'a>(spine: Option<&'a SpineHandle>, msg: &str) -> Result<&'a SpineHandle, String> {
    let handle = spine.ok_or_else(|| msg.to_string())?;
    let is_fil = handle
        .read()
        .map_err(|_| "SimulSurf: poisoned spine".to_string())?
        .class()
        == SpineClass::Fil;
    if !is_fil {
        return Err(msg.to_string());
    }
    Ok(handle)
}

/// `HS->Surface()` of one of the four `BRepAdaptor_Surface` handles.
fn adaptor_surface(
    adaptor: &SurfaceAdaptor,
    what: &str,
) -> Result<Arc<dyn occt_geom::Surface>, String> {
    adaptor
        .surface()
        .cloned()
        .ok_or_else(|| format!("SimulSurf: {what} has no surface"))
}

/// The four 2d coordinates of the section loop (`cxx:888-911`, `:969-991`,
/// `:1119-1142`, `:1200-1224`).
struct SectionPoints {
    /// `pf` / `pl`: `ParametersOnS` of the first / last point.
    pf: GpPnt2d,
    pl: GpPnt2d,
    /// `ppcf` / `ppcl`: `ParametersOnS2` of the first / last point.
    ppcf: GpPnt2d,
    ppcl: GpPnt2d,
}

/// The section loop shared by the four arms of the two overloads
/// (`cxx:886-911`, `:968-991`, `:1117-1142`, `:1199-1224`).
///
/// `section` is `func.Section(param, u, v, w, p1, p2, ci)`; the two `SurfRst`
/// functions share the signature, so the two loops are one here.
fn section_from_line(
    lin: &Option<Rc<RefCell<Line>>>,
    section: impl FnMut(f64, f64, f64, f64) -> (f64, f64, GpCirc),
) -> Result<(CircSectionArray, SectionPoints), String> {
    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let guard = line.borrow();
    // `int nbp = lin->NbPoints(); sec = new
    // NCollection_HArray1<ChFiDS_CircSection>(1, nbp);` (`cxx:886-887`).
    let nbp = guard.nb_points();
    let mut sec = CircSectionArray::new(1, nbp as i32);
    let mut section = section;
    let mut pts = SectionPoints {
        pf: GpPnt2d::new(0.0, 0.0),
        pl: GpPnt2d::new(0.0, 0.0),
        ppcf: GpPnt2d::new(0.0, 0.0),
        ppcl: GpPnt2d::new(0.0, 0.0),
    };
    for i in 1..=nbp {
        // `p.ParametersOnS(u, v); w = p.ParameterOnC(); param = p.Parameter();`
        let p = guard
            .point(i)
            .ok_or_else(|| format!("SimulSurf: Lin->Point({i}) out of range"))?;
        let (u, v) = p.parameters_on_s()?;
        let w = p.parameter_on_c()?;
        let param = p.parameter();
        // `func.Section(param, u, v, w, p1, p2, ci); isec.Set(ci, p1, p2);`
        let (p1, p2, ci) = section(param, u, v, w);
        sec.change_value(i as i32)
            .ok_or_else(|| format!("SimulSurf: sec->ChangeValue({i}) out of range"))?
            .set_circ(&ci, p1, p2);
        if i == 1 {
            pts.pf = GpPnt2d::new(u, v);
            let (u2, v2) = p.parameters_on_s2()?;
            pts.ppcf = GpPnt2d::new(u2, v2);
        }
        if i == nbp {
            pts.pl = GpPnt2d::new(u, v);
            let (u2, v2) = p.parameters_on_s2()?;
            pts.ppcl = GpPnt2d::new(u2, v2);
        }
    }
    Ok((sec, pts))
}

/// The four `ChFi3d_FilCommonPoint` calls that close both overloads
/// (`cxx:997-1015`, `:1228-1247`).
///
/// `on_s1` is the edge / face overload that walks *from* `S1`: it pairs the
/// walk's "first" extremity with the `S2` vertices and its "second" one with
/// the `S1` vertices (`cxx:997-1015`), the mirror of the `S2` overload
/// (`cxx:1228-1247`).
fn fill_common_points(
    data: &mut SurfData,
    line: &Line,
    tolapp3d: f64,
    on_s1: bool,
) -> Result<(), String> {
    let trans1 = line.transition_on_s1()?;
    let trans2 = line.transition_on_s2()?;
    if on_s1 {
        fil_common_point(
            line.start_point_on_first(),
            trans1,
            true,
            data.vertex_first_on_s2_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.end_point_on_first(),
            trans1,
            false,
            data.vertex_last_on_s2_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.start_point_on_second(),
            trans2,
            true,
            data.vertex_first_on_s1_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.end_point_on_second(),
            trans2,
            false,
            data.vertex_last_on_s1_mut(),
            tolapp3d,
        );
    } else {
        fil_common_point(
            line.start_point_on_first(),
            trans1,
            true,
            data.vertex_first_on_s1_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.end_point_on_first(),
            trans1,
            false,
            data.vertex_last_on_s1_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.start_point_on_second(),
            trans2,
            true,
            data.vertex_first_on_s2_mut(),
            tolapp3d,
        );
        fil_common_point(
            line.end_point_on_second(),
            trans2,
            false,
            data.vertex_last_on_s2_mut(),
            tolapp3d,
        );
    }
    Ok(())
}

/// `SimulSurf(...)` on `S1`: face `HS2` against the restriction `PC1` of `HS1`
/// (`cxx:787-1017`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn simul_surf_on_s1(
    shape: BlendFuncSectionShape,
    data: &mut SurfData,
    guide: &Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    choix: i32,
    hs1: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    pc1: Option<&Curve2dAdaptor>,
    hsref1: &SurfaceAdaptor,
    pcref1: &Curve2dAdaptor,
    decroch1: &mut bool,
    hs2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    or2: Orientation,
    tol_guide: f64,
    first: &mut f64,
    last: &mut f64,
    inside: bool,
    approx: bool,
    forward: bool,
    rec_p: bool,
    rec_s: bool,
    rec_rst: bool,
    soldep: &MathVector,
    ctx: &ChFi3dCtx<'_>,
) -> Result<(), String> {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of the
    // fillet");` (`cxx:815-818`).
    let spine_handle = fil_spine(spine, "PerformSurf : this is not the spine of the fillet")?;

    // `SimulParams(HGuide, fsp, MaxStep, locfleche);` (`cxx:821-822`).
    let (max_step, locfleche) = simul_params(guide, spine_handle)?;

    // `double PFirst = First;` (`cxx:826`).
    let p_first = *first;

    // `func(HS2, HS1, PC1, HGuide)` — the walk is on `HS2` against the
    // restriction `PC1` of `HS1` (`cxx:829`).
    let s2_surf = adaptor_surface(hs2, "HS2")?;
    let s1_surf = adaptor_surface(hs1, "HS1")?;
    let pc1_adaptor = pc1
        .cloned()
        .ok_or_else(|| "SimulSurf: PC1 is null".to_string())?;
    let pc1_curve = pc1_adaptor
        .curve()
        .cloned()
        .ok_or_else(|| "SimulSurf: PC1 has no p-curve".to_string())?;
    let ref_surf = adaptor_surface(hsref1, "HSref1")?;
    let ref_curve = pcref1
        .curve()
        .cloned()
        .ok_or_else(|| "SimulSurf: PCref1 has no p-curve".to_string())?;

    let constant = spine_handle
        .read()
        .map_err(|_| "SimulSurf: poisoned spine".to_string())?
        .is_fil_constant();

    let mut lin: Option<Rc<RefCell<Line>>> = None;
    let (sec, pts) = if constant {
        // `BRepBlend_SurfRstConstRad func(HS2, HS1, PC1, HGuide);
        // func.Set(HSref1, PCref1);` (`cxx:829-830`).
        let mut func = SurfRstConstRad::new(
            Arc::clone(&s2_surf),
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            guide_as_curve(guide),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        // `HC->Load(PC1, HS1); BRepBlend_SurfCurvConstRadInv finvc(HS2, HC,
        // HGuide); BRepBlend_SurfPointConstRadInv finvp(HS2, HGuide);
        // BRepBlend_ConstRadInv finv(HS2, HSref1, HGuide);
        // finv.Set(false, PCref1);` (`cxx:831-838`).
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc1_curve), Arc::clone(&s1_surf)));
        let mut finvc =
            SurfCurvConstRadInv::new(Arc::clone(&s2_surf), hc, guide_as_curve(guide));
        let mut finvp = SurfPointConstRadInv::new(Arc::clone(&s2_surf), guide_as_curve(guide));
        let mut finv = ConstRadInv::new(
            Arc::clone(&s2_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));

        let rad = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_radius()
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string()
            })?;
        let petit = petit_choix(choix, or2);
        // `finv.Set(rad, Choix); finvc.Set(rad, petitchoix);
        // finvp.Set(rad, petitchoix); func.Set(rad, petitchoix);
        // func.Set(myShape);` (`cxx:847-852`).
        finv.set_radius(rad, choix);
        finvc.set_radius(rad, petit);
        finvp.set_radius(rad, petit);
        func.set_radius(rad, petit);
        func.set_section_shape(shape);

        // `done = SimulData(Data, HGuide, lin, HS2, I2, HS1, PC1, I1,
        // Decroch1, func, finv, finvp, finvc, PFirst, MaxStep, locfleche,
        // TolGuide, First, Last, Soldep, 4, Inside, Appro, Forward, RecP, RecS,
        // RecRst);` (`cxx:854-880`).
        let done = simul_data_surf_rst(
            ctx.tolapp3d,
            ctx.tolapp2d,
            guide,
            &mut lin,
            Arc::clone(&s2_surf),
            Rc::clone(i2),
            Arc::clone(&s1_surf),
            pc1_adaptor,
            Rc::clone(i1),
            decroch1,
            &mut func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            locfleche,
            tol_guide,
            first,
            last,
            soldep,
            4,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
        );
        // `if (!done) throw Standard_Failure("SimulSurf : Failed process!");`
        // (`cxx:881-884`).
        if !done {
            return Err("SimulSurf : Failed process!".to_string());
        }
        section_from_line(&lin, |param, u, v, w| func.section_circle(param, u, v, w))?
    } else {
        // `BRepBlend_SurfRstEvolRad func(HS2, HS1, PC1, HGuide,
        // fsp->Law(HGuide));` (`cxx:914`).
        let law = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_law(guide)
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Law raised Standard_DomainError".to_string()
            })?;
        let law_func: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
        let mut func = SurfRstEvolRad::new(
            Arc::clone(&s2_surf),
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc1_curve), Arc::clone(&s1_surf)));
        let mut finvc = SurfCurvEvolRadInv::new(
            Arc::clone(&s2_surf),
            hc,
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        let mut finvp = SurfPointEvolRadInv::new(
            Arc::clone(&s2_surf),
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s2_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
            law_func,
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));
        let petit = petit_choix(choix, or2);
        // `finv.Set(Choix); finvc.Set(petitchoix); finvp.Set(petitchoix);
        // func.Set(petitchoix); func.Set(myShape);` (`cxx:930-935`).
        finv.set_choix(choix);
        finvc.set_choix(petit);
        finvp.set_choix(petit);
        func.set_choix(petit);
        func.set_section_shape(shape);

        let done = simul_data_surf_rst(
            ctx.tolapp3d,
            ctx.tolapp2d,
            guide,
            &mut lin,
            Arc::clone(&s2_surf),
            Rc::clone(i2),
            Arc::clone(&s1_surf),
            pc1_adaptor,
            Rc::clone(i1),
            decroch1,
            &mut func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            locfleche,
            tol_guide,
            first,
            last,
            soldep,
            4,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
        );
        // `if (!done) throw Standard_Failure("SimulSurf : Fail !");`
        // (`cxx:963-966`).
        if !done {
            return Err("SimulSurf : Fail !".to_string());
        }
        section_from_line(&lin, |param, u, v, w| func.section_circle(param, u, v, w))?
    };

    // `Data->SetSimul(sec); Data->Set2dPoints(ppcf, ppcl, pf, pl);`
    // (`cxx:993-996`): the two points on the surface side come first in this
    // overload.
    data.set_simul(Arc::new(sec));
    data.set_2d_points(pts.ppcf, pts.ppcl, pts.pf, pts.pl);
    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let guard = line.borrow();
    fill_common_points(data, &guard, ctx.tolapp3d, true)?;
    Ok(())
}

/// `SimulSurf(...)` on `S2`: face `HS1` against the restriction `PC2` of `HS2`
/// (`cxx:1020-1248`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn simul_surf_on_s2(
    shape: BlendFuncSectionShape,
    data: &mut SurfData,
    guide: &Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    choix: i32,
    hs1: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    or1: Orientation,
    hs2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    pc2: Option<&Curve2dAdaptor>,
    hsref2: &SurfaceAdaptor,
    pcref2: &Curve2dAdaptor,
    decroch2: &mut bool,
    tol_guide: f64,
    first: &mut f64,
    last: &mut f64,
    inside: bool,
    approx: bool,
    forward: bool,
    rec_p: bool,
    rec_s: bool,
    rec_rst: bool,
    soldep: &MathVector,
    ctx: &ChFi3dCtx<'_>,
) -> Result<(), String> {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : it is not the spine of a
    // fillet");` (`cxx:1045-1049`).
    let spine_handle = fil_spine(spine, "PerformSurf : it is not the spine of a fillet")?;

    let (max_step, locfleche) = simul_params(guide, spine_handle)?;
    let p_first = *first;

    let s1_surf = adaptor_surface(hs1, "HS1")?;
    let s2_surf = adaptor_surface(hs2, "HS2")?;
    let pc2_adaptor = pc2
        .cloned()
        .ok_or_else(|| "SimulSurf: PC2 is null".to_string())?;
    let pc2_curve = pc2_adaptor
        .curve()
        .cloned()
        .ok_or_else(|| "SimulSurf: PC2 has no p-curve".to_string())?;
    let ref_surf = adaptor_surface(hsref2, "HSref2")?;
    let ref_curve = pcref2
        .curve()
        .cloned()
        .ok_or_else(|| "SimulSurf: PCref2 has no p-curve".to_string())?;

    let constant = spine_handle
        .read()
        .map_err(|_| "SimulSurf: poisoned spine".to_string())?
        .is_fil_constant();

    let mut lin: Option<Rc<RefCell<Line>>> = None;
    let (sec, pts) = if constant {
        // `BRepBlend_SurfRstConstRad func(HS1, HS2, PC2, HGuide);
        // func.Set(HSref2, PCref2);` (`cxx:1061-1062`).
        let mut func = SurfRstConstRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Arc::clone(&pc2_curve),
            guide_as_curve(guide),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc2_curve), Arc::clone(&s2_surf)));
        let mut finvc =
            SurfCurvConstRadInv::new(Arc::clone(&s1_surf), hc, guide_as_curve(guide));
        let mut finvp = SurfPointConstRadInv::new(Arc::clone(&s1_surf), guide_as_curve(guide));
        let mut finv = ConstRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));

        let rad = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_radius()
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string()
            })?;
        let petit = petit_choix(choix, or1);
        finv.set_radius(rad, choix);
        finvc.set_radius(rad, petit);
        finvp.set_radius(rad, petit);
        func.set_radius(rad, petit);
        func.set_section_shape(shape);

        // `SimulData(Data, HGuide, lin, HS1, I1, HS2, PC2, I2, Decroch2, ...)`
        // (`cxx:1086-1112`).
        let done = simul_data_surf_rst(
            ctx.tolapp3d,
            ctx.tolapp2d,
            guide,
            &mut lin,
            Arc::clone(&s1_surf),
            Rc::clone(i1),
            Arc::clone(&s2_surf),
            pc2_adaptor,
            Rc::clone(i2),
            decroch2,
            &mut func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            locfleche,
            tol_guide,
            first,
            last,
            soldep,
            4,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
        );
        // `if (!done) throw Standard_Failure("SimulSurf : Failed Processing!");`
        // (`cxx:1113-1116`).
        if !done {
            return Err("SimulSurf : Failed Processing!".to_string());
        }
        section_from_line(&lin, |param, u, v, w| func.section_circle(param, u, v, w))?
    } else {
        // `BRepBlend_SurfRstEvolRad func(HS1, HS2, PC2, HGuide,
        // fsp->Law(HGuide));` (`cxx:1146`).
        let law = spine_handle
            .read()
            .map_err(|_| "SimulSurf: poisoned spine".to_string())?
            .fil_law(guide)
            .ok_or_else(|| {
                "SimulSurf: ChFiDS_FilSpine::Law raised Standard_DomainError".to_string()
            })?;
        let law_func: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
        let mut func = SurfRstEvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Arc::clone(&pc2_curve),
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc2_curve), Arc::clone(&s2_surf)));
        let mut finvc = SurfCurvEvolRadInv::new(
            Arc::clone(&s1_surf),
            hc,
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        let mut finvp = SurfPointEvolRadInv::new(
            Arc::clone(&s1_surf),
            guide_as_curve(guide),
            Arc::clone(&law_func),
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
            law_func,
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));
        let petit = petit_choix(choix, or1);
        finv.set_choix(choix);
        finvc.set_choix(petit);
        finvp.set_choix(petit);
        func.set_choix(petit);
        func.set_section_shape(shape);

        let done = simul_data_surf_rst(
            ctx.tolapp3d,
            ctx.tolapp2d,
            guide,
            &mut lin,
            Arc::clone(&s1_surf),
            Rc::clone(i1),
            Arc::clone(&s2_surf),
            pc2_adaptor,
            Rc::clone(i2),
            decroch2,
            &mut func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            locfleche,
            tol_guide,
            first,
            last,
            soldep,
            4,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
        );
        // `if (!done) throw Standard_Failure("SimulSurf : Fail !");`
        // (`cxx:1195-1198`).
        if !done {
            return Err("SimulSurf : Fail !".to_string());
        }
        section_from_line(&lin, |param, u, v, w| func.section_circle(param, u, v, w))?
    };

    // `Data->SetSimul(sec); Data->Set2dPoints(pf, pl, ppcf, ppcl);`
    // (`cxx:1226-1229`).
    data.set_simul(Arc::new(sec));
    data.set_2d_points(pts.pf, pts.pl, pts.ppcf, pts.ppcl);
    let line = lin
        .as_ref()
        .ok_or_else(|| "SimulSurf: Lin is null".to_string())?
        .clone();
    let guard = line.borrow();
    fill_common_points(data, &guard, ctx.tolapp3d, false)?;
    Ok(())
}
