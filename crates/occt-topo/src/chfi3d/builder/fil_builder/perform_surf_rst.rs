//! `ChFi3d_FilBuilder::PerformSurf` — the two edge / face overloads.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx` (OCCT
//! 8.0.0):
//!
//! - [`perform_surf_on_s1`] — `PerformSurf(SeqData, HGuide, Spine, Choix, HS1,
//!   I1, PC1, HSref1, PCref1, Decroch1, HS2, I2, Or2, MaxStep, Fleche,
//!   TolGuide, First, Last, Inside, Appro, Forward, RecP, RecS, RecRst,
//!   Soldep)` (`cxx:1721-1890`), the band between the face `HS2` and the
//!   restriction `PC1` of `HS1`,
//! - [`perform_surf_on_s2`] — the mirror overload (`cxx:1894-2063`).
//!
//! Both walk with `BRepBlend_SurfRstLineBuilder` through
//! `ChFi3d_Builder::ComputeData(SurfRst)` (`ChFi3d_Builder_6.cxx:1009-1171`),
//! approximate with `ChFi3d_Builder::CompleteData(Data, Func, lin, S1, S2, Or,
//! Reversed)` (`cxx:557-586`) and store one `ChFiDS_SurfData` of the result.
//!
//! Divergences that are not heuristics:
//! - the three `throw`s of each body — `Standard_ConstructionError("PerformSurf
//!   : this is not the spine of a fillet")` (`cxx:1748-1752`, `:1922-1926`),
//!   `Standard_Failure("PerformSurf : Failed processing!")` (`cxx:1810-1814`,
//!   `:1874-1878`, `:1982-1986`, `:2046-2050`) and
//!   `Standard_Failure("PerformSurf : Failed approximation!")`
//!   (`cxx:1817-1821`, `:1878-1882`, `:1988-1992`, `:2052-2056`) — are the
//!   trait's `false`; `Spine->SetErrorStatus(ChFiDS_WalkingFailure)` is written
//!   before the second one, exactly where OCCT writes it,
//! - `HS2->Face().Orientation()` (on `S1`) and `HS1->Face().Orientation()` (on
//!   `S2`) raise `Standard_NullObject` for an adaptor without a face; that arm
//!   answers `false` too,
//! - `Func.GetMinimalDistance()` answers a `Result` in the port; the `Err` arm
//!   is treated as "not singular", which is where `Blend_Function` leaves
//!   `distmin` until a section lowers it (same rule as the face / face
//!   overload),
//! - the two arms of each overload share [`run_arm_rst`]: the `ComputeData` /
//!   `CompleteData` / `SplitSurf` tail is byte-identical in OCCT
//!   (`cxx:1784-1824` vs `:1845-1885`, `:1956-1996` vs `:2020-2060`), so the
//!   blend function is passed as the generic `F` OCCT passes as
//!   `BRepBlend_SurfRstFunction&`,
//! - `const double /*TolGuide*/` is commented out in the `hxx` declaration of
//!   the two overloads but present in the `cxx` bodies; the trait keeps the
//!   slot and passes `tolesp`, as `CallPerformSurf` does.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::law::{LawFunction, LawHandle};
use occt_geom::Surface;
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::func_inv::BlendFuncInv;
use crate::blend::{
    BlendSurfCurvFuncInv, BlendSurfPointFuncInv, BlendSurfRstFunction, ConstRadInv, EvolRadInv,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::inv::{
    SurfCurvConstRadInv, SurfCurvEvolRadInv, SurfPointConstRadInv, SurfPointEvolRadInv,
};
use crate::brep_blend::{Line, SurfRstConstRad, SurfRstEvolRad};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi_ds::{ElSpine, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::perform_surf::{face_orientation, fil_spine};
use super::simul_surf::guide_as_curve;
use super::simul_surf_rst::petit_choix;

/// `Spine->SetErrorStatus(ChFiDS_WalkingFailure)` (`cxx:1812`, `:1876`,
/// `:1984`, `:2048`, `:2177`, `:2255`), followed by the `Standard_Failure`
/// that is the caller's `false`.
pub(super) fn walking_failure(spine: Option<&SpineHandle>) {
    if let Some(spine) = spine {
        if let Ok(mut guard) = spine.write() {
            guard.set_error_status(crate::chfi_ds::ErrorStatus::WalkingFailure);
        }
    }
}

/// `fsp->Radius()` (`cxx:1768`, `:1941`, `:2126`): a `Standard_DomainError` in
/// OCCT for a non-constant spine, which never reaches here because the caller
/// is on the `IsConstant()` arm.
pub(super) fn fil_radius(spine: &SpineHandle) -> Option<f64> {
    spine.read().ok()?.fil_radius()
}

/// `fsp->Law(HGuide)` (`cxx:1825`, `:1829`, `:2000`, `:2197`): a
/// `Standard_DomainError` in OCCT for a constant spine.
pub(super) fn fil_law(spine: &SpineHandle, guide: &Arc<ElSpine>) -> Option<Arc<dyn LawFunction>> {
    let law = spine.read().ok()?.fil_law(guide)?;
    Some(Arc::new(LawHandle::new(law)))
}

/// The `ComputeData` / `CompleteData` / `SplitSurf` tail shared by the four
/// arms of the two overloads.
///
/// `s1` / `i1` / `s2` / `pc2` / `i2` are the walk's own arguments, i.e. the
/// `HS2, I2, HS1, PC1, I1` of the `S1` overload (`cxx:1786-1792`) and the
/// `HS1, I1, HS2, PC2, I2` of the `S2` one (`cxx:1959-1965`); `hs1` / `hs2`
/// are the two adaptors `CompleteData` receives, which is `HS1, HS2` in *both*
/// overloads (`cxx:1817`, `:1882`, `:1988`, `:2054`).
#[allow(clippy::too_many_arguments)]
fn run_arm_rst<F: BlendSurfRstFunction + 'static>(
    base: &mut ChFi3dBuilder,
    seq_data: &mut Vec<Arc<SurfData>>,
    guide: &Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    s1: Arc<dyn Surface>,
    i1: &Rc<RefCell<TopolTool>>,
    s2: Arc<dyn Surface>,
    pc2: Curve2dAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    decroch: &mut bool,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    or: Option<Orientation>,
    reversed: bool,
    func: F,
    finv: &mut dyn BlendFuncInv,
    finv_p: &mut dyn BlendSurfPointFuncInv,
    finv_c: &mut dyn BlendSurfCurvFuncInv,
    p_first: f64,
    max_step: f64,
    fleche: f64,
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
) -> bool {
    if seq_data.is_empty() {
        return false;
    }
    // `Or` is read only after the walk succeeded (`cxx:1816`, `:1878`).
    let Some(or) = or else {
        return false;
    };

    // `BRepBlend_AppFuncRst` keeps the `Buff_Curve` of the line and the
    // function by reference in OCCT; the port shares it through
    // `Rc<RefCell<..>>`.
    let func = Rc::new(RefCell::new(func));
    let mut lin: Option<Rc<RefCell<Line>>> = None;

    // `done = ComputeData(Data, HGuide, lin, ..., Decroch, Func, FInv, FInvP,
    // FInvC, PFirst, MaxStep, Fleche, TolGuide, First, Last, Soldep, Inside,
    // Appro, Forward, RecP, RecS, RecRst);` (`cxx:1784-1809`, `:1845-1872`,
    // `:1956-1981`, `:2020-2045`).
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.compute_data_surf_rst(
            data,
            guide,
            &mut lin,
            s1,
            Rc::clone(i1),
            s2,
            pc2,
            Rc::clone(i2),
            decroch,
            &mut *func.borrow_mut(),
            finv,
            finv_p,
            finv_c,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            soldep,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
        )
    };
    // `if (!done) { Spine->SetErrorStatus(ChFiDS_WalkingFailure); throw
    // Standard_Failure("PerformSurf : Failed processing!"); }`.
    if !done {
        walking_failure(spine);
        return false;
    }

    // `done = CompleteData(Data, Func, lin, HS1, HS2, Or, true /*false*/);`
    let Some(line) = lin.clone() else {
        return false;
    };
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.complete_data_rst(
            data,
            Rc::clone(&func),
            Rc::clone(&line),
            hs1,
            Some(hs2),
            or,
            reversed,
        )
    };
    // `if (!done) throw Standard_Failure("PerformSurf : Failed
    // approximation!");` — see the module note.
    if !matches!(done, Ok(true)) {
        return false;
    }

    // `maybesingular = (Func.GetMinimalDistance() <= 100 * tolapp3d);`
    let maybe_singular = match func.borrow().get_minimal_distance() {
        Ok(distance) => distance <= 100.0 * base.tolapp3d,
        Err(_) => false,
    };

    // `if (maybesingular) { SplitSurf(SeqData, lin); }` (`cxx:1886-1889`,
    // `:2059-2062`).
    if maybe_singular {
        let guard = line.borrow();
        if base.split_surf(seq_data, &guard).is_err() {
            return false;
        }
    }
    true
}

/// `PerformSurf(...)` on `S1`: the face `HS2` against the restriction `PC1` of
/// `HS1` (`cxx:1721-1890`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_surf_on_s1(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    seq_data: &mut Vec<Arc<SurfData>>,
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
    max_step: f64,
    fleche: f64,
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
) -> bool {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of a
    // fillet");` (`cxx:1748-1752`).
    let Some(spine_handle) = fil_spine(spine) else {
        return false;
    };
    // `double PFirst = First;` (`cxx:1755`).
    let p_first = *first;

    // The two adaptor surfaces and the p-curves; a null adaptor is OCCT's null
    // handle, dereferenced below.
    let (Some(s1_surf), Some(s2_surf), Some(pc1_adaptor), Some(pc1_curve), Some(ref_surf), Some(ref_curve)) = (
        hs1.surface().cloned(),
        hs2.surface().cloned(),
        pc1.cloned(),
        pc1.and_then(|pc| pc.curve().cloned()),
        hsref1.surface().cloned(),
        pcref1.curve().cloned(),
    ) else {
        return false;
    };

    let constant = match spine_handle.read() {
        Ok(guard) => guard.is_fil_constant(),
        Err(_) => return false,
    };

    if constant {
        // `BRepBlend_SurfRstConstRad func(HS2, HS1, PC1, HGuide);
        // func.Set(HSref1, PCref1); HC->Load(PC1, HS1);
        // BRepBlend_SurfCurvConstRadInv finvc(HS2, HC, HGuide);
        // BRepBlend_SurfPointConstRadInv finvp(HS2, HGuide);
        // BRepBlend_ConstRadInv finv(HS2, HSref1, HGuide);
        // finv.Set(false, PCref1);` (`cxx:1759-1767`).
        let mut func = SurfRstConstRad::new(
            Arc::clone(&s2_surf),
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            guide_as_curve(guide),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
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

        // `double rad = fsp->Radius(); int petitchoix = ...` (`cxx:1768-1779`).
        let Some(rad) = fil_radius(spine_handle) else {
            return false;
        };
        let petit = petit_choix(choix, or2);
        // `finv.Set(rad, Choix); finvc.Set(rad, petitchoix);
        // finvp.Set(rad, petitchoix); func.Set(rad, petitchoix);
        // func.Set(myShape);` (`cxx:1776-1782`).
        finv.set_radius(rad, choix);
        finvc.set_radius(rad, petit);
        finvp.set_radius(rad, petit);
        func.set_radius(rad, petit);
        func.set_section_shape(shape);

        run_arm_rst(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s2_surf),
            i2,
            Arc::clone(&s1_surf),
            pc1_adaptor,
            i1,
            decroch1,
            hs1,
            hs2,
            // `Or = HS2->Face().Orientation();` (`cxx:1816`).
            face_orientation(hs2),
            true,
            func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        )
    } else {
        // `BRepBlend_SurfRstEvolRad func(HS2, HS1, PC1, HGuide,
        // fsp->Law(HGuide)); ...` (`cxx:1825-1834`).
        let Some(law) = fil_law(spine_handle, guide) else {
            return false;
        };
        let mut func = SurfRstEvolRad::new(
            Arc::clone(&s2_surf),
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc1_curve), Arc::clone(&s1_surf)));
        let mut finvc = SurfCurvEvolRadInv::new(
            Arc::clone(&s2_surf),
            hc,
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        let mut finvp = SurfPointEvolRadInv::new(
            Arc::clone(&s2_surf),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s2_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
            law,
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));

        let petit = petit_choix(choix, or2);
        // `finv.Set(Choix); finvc.Set(petitchoix); finvp.Set(petitchoix);
        // func.Set(petitchoix); func.Set(myShape);` (`cxx:1841-1846`).
        finv.set_choix(choix);
        finvc.set_choix(petit);
        finvp.set_choix(petit);
        func.set_choix(petit);
        func.set_section_shape(shape);

        run_arm_rst(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s2_surf),
            i2,
            Arc::clone(&s1_surf),
            pc1_adaptor,
            i1,
            decroch1,
            hs1,
            hs2,
            face_orientation(hs2),
            true,
            func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        )
    }
}

/// `PerformSurf(...)` on `S2`: the face `HS1` against the restriction `PC2` of
/// `HS2` (`cxx:1894-2063`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_surf_on_s2(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    seq_data: &mut Vec<Arc<SurfData>>,
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
    max_step: f64,
    fleche: f64,
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
) -> bool {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of a
    // fillet");` (`cxx:1922-1926`).
    let Some(spine_handle) = fil_spine(spine) else {
        return false;
    };
    let p_first = *first;

    let (Some(s1_surf), Some(s2_surf), Some(pc2_adaptor), Some(pc2_curve), Some(ref_surf), Some(ref_curve)) = (
        hs1.surface().cloned(),
        hs2.surface().cloned(),
        pc2.cloned(),
        pc2.and_then(|pc| pc.curve().cloned()),
        hsref2.surface().cloned(),
        pcref2.curve().cloned(),
    ) else {
        return false;
    };

    let constant = match spine_handle.read() {
        Ok(guard) => guard.is_fil_constant(),
        Err(_) => return false,
    };

    if constant {
        // `BRepBlend_SurfRstConstRad func(HS1, HS2, PC2, HGuide);
        // func.Set(HSref2, PCref2); HC->Load(PC2, HS2);
        // BRepBlend_SurfCurvConstRadInv finvc(HS1, HC, HGuide);
        // BRepBlend_SurfPointConstRadInv finvp(HS1, HGuide);
        // BRepBlend_ConstRadInv finv(HS1, HSref2, HGuide);
        // finv.Set(false, PCref2);` (`cxx:1932-1939`).
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

        // `double rad = fsp->Radius(); int petitchoix = ...` (`cxx:1941-1952`).
        let Some(rad) = fil_radius(spine_handle) else {
            return false;
        };
        let petit = petit_choix(choix, or1);
        // `finv.Set(rad, Choix); finvc.Set(rad, petitchoix); ...` (`cxx:1952`).
        finv.set_radius(rad, choix);
        finvc.set_radius(rad, petit);
        finvp.set_radius(rad, petit);
        func.set_radius(rad, petit);
        func.set_section_shape(shape);

        run_arm_rst(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s1_surf),
            i1,
            Arc::clone(&s2_surf),
            pc2_adaptor,
            i2,
            decroch2,
            hs1,
            hs2,
            // `Or = HS1->Face().Orientation();` (`cxx:1988`).
            face_orientation(hs1),
            false,
            func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        )
    } else {
        // `BRepBlend_SurfRstEvolRad func(HS1, HS2, PC2, HGuide,
        // fsp->Law(HGuide)); ...` (`cxx:1997-2007`).
        let Some(law) = fil_law(spine_handle, guide) else {
            return false;
        };
        let mut func = SurfRstEvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Arc::clone(&pc2_curve),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        func.set_refs(Arc::clone(&ref_surf), Arc::clone(&ref_curve));
        let hc: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc2_curve), Arc::clone(&s2_surf)));
        let mut finvc = SurfCurvEvolRadInv::new(
            Arc::clone(&s1_surf),
            hc,
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        let mut finvp = SurfPointEvolRadInv::new(
            Arc::clone(&s1_surf),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&ref_surf),
            guide_as_curve(guide),
            law,
        );
        finv.set(false, Some(Arc::clone(&ref_curve)));

        let petit = petit_choix(choix, or1);
        // `finv.Set(Choix); finvc.Set(petitchoix); ...` (`cxx:2015-2021`).
        finv.set_choix(choix);
        finvc.set_choix(petit);
        finvp.set_choix(petit);
        func.set_choix(petit);
        func.set_section_shape(shape);

        run_arm_rst(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s1_surf),
            i1,
            Arc::clone(&s2_surf),
            pc2_adaptor,
            i2,
            decroch2,
            hs1,
            hs2,
            face_orientation(hs1),
            false,
            func,
            &mut finv,
            &mut finvp,
            &mut finvc,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        )
    }
}
