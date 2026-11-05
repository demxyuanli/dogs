//! `ChFi3d_FilBuilder::PerformSurf` — the edge / edge overload.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:2065-2270`
//! (OCCT 8.0.0). It walks the band between the restriction `PC1` of `HS1` and
//! the restriction `PC2` of `HS2`; `fsp->IsConstant()` (`cxx:2110`) splits the
//! body between `BRepBlend_RstRstConstRad` and `BRepBlend_RstRstEvolRad`, and
//! both arms then run the same `ComputeData(RstRst)` / `CompleteData(Data,
//! Func, lin, HS1, HS2, Or)` / `SplitSurf` tail.
//!
//! Divergences that are not heuristics:
//! - the three `throw`s — `Standard_ConstructionError("PerformSurf : this is
//!   not the spine of a fillet")` (`cxx:2102-2106`),
//!   `Standard_Failure("PerformSurf : Failed processing!")`
//!   (`cxx:2175-2179`, `:2253-2257`) and
//!   `Standard_Failure("PerformSurf : Failed approximation!")`
//!   (`cxx:2181-2185`, `:2259-2263`) — are the trait's `Err`;
//!   `Spine->SetErrorStatus(ChFiDS_WalkingFailure)` is written before the
//!   second one, exactly where OCCT writes it,
//! - `HS1->Face().Orientation()` (`cxx:2179`, `:2258`) raises
//!   `Standard_NullObject` for an adaptor without a face; that arm answers
//!   `false` too,
//! - `Func.GetMinimalDistance()` answers a `Result` in the port; the `Err` arm
//!   is treated as "not singular", which is where `Blend_Function` leaves
//!   `distmin` until a section lowers it (same rule as the other three
//!   overloads),
//! - `CompleteData(Data, Func, lin, HS1, HS2, Or)` has no `Reversed` argument
//!   (`ChFi3d_Builder_6.cxx:594-624`), unlike the surf / rst overload, so
//!   [`ChFi3dBuilder::complete_data_rst_rst`] takes no such parameter either,
//! - the two arms share [`run_arm_contact`]: the tail is byte-identical in
//!   OCCT (`cxx:2142-2186` vs `:2222-2264`), so the blend function is passed
//!   as the generic `F` OCCT passes as `BRepBlend_RstRstFunction&`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::Surface;
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{
    BlendCurvPointFuncInv, BlendRstRstFunction, BlendSurfCurvFuncInv,
};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::inv::{CurvPointRadInv, SurfCurvConstRadInv, SurfCurvEvolRadInv};
use crate::brep_blend::{Line, SurfRstRstConstRad, SurfRstRstEvolRad};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi_ds::{ElSpine, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::perform_surf::{face_orientation, fil_spine};
use super::perform_surf_rst::{fil_law, fil_radius, walking_failure};
use super::simul_surf::guide_as_curve;

/// `ch1 = 1; if (Or1 == TopAbs_REVERSED) { ch1 = 3; }` (`cxx:2125-2132`,
/// `:2204-2213`).
fn contact_choix(or: Orientation) -> i32 {
    if or == Orientation::Reversed {
        3
    } else {
        1
    }
}

/// The `ComputeData(RstRst)` / `CompleteData` / `SplitSurf` tail shared by the
/// two arms of [`perform_surf_contact`] (`cxx:2142-2186`, `:2222-2264`).
#[allow(clippy::too_many_arguments)]
fn run_arm_contact<F: BlendRstRstFunction + 'static>(
    base: &mut ChFi3dBuilder,
    seq_data: &mut Vec<Arc<SurfData>>,
    guide: &Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    s1: Arc<dyn Surface>,
    pc1: Curve2dAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    decroch1: &mut bool,
    s2: Arc<dyn Surface>,
    pc2: Curve2dAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    decroch2: &mut bool,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    func: F,
    finv1: &mut dyn BlendSurfCurvFuncInv,
    finv_p1: &mut dyn BlendCurvPointFuncInv,
    finv2: &mut dyn BlendSurfCurvFuncInv,
    finv_p2: &mut dyn BlendCurvPointFuncInv,
    p_first: f64,
    max_step: f64,
    fleche: f64,
    tol_guide: f64,
    first: &mut f64,
    last: &mut f64,
    inside: bool,
    approx: bool,
    forward: bool,
    rec_p1: bool,
    rec_rst1: bool,
    rec_p2: bool,
    rec_rst2: bool,
    soldep: &MathVector,
) -> bool {
    if seq_data.is_empty() {
        return false;
    }
    // `TopAbs_Orientation Or = HS1->Face().Orientation();` — read only after
    // the walk succeeded (`cxx:2179`, `:2258`).
    let Some(or) = face_orientation(hs1) else {
        return false;
    };

    // `BRepBlend_AppFuncRstRst` keeps the `Buff_Curve` of the line and the
    // function by reference in OCCT; the port shares it through
    // `Rc<RefCell<..>>`.
    let func = Rc::new(RefCell::new(func));
    let mut lin: Option<Rc<RefCell<Line>>> = None;

    // `done = ComputeData(Data, HGuide, lin, HS1, PC1, I1, Decroch1, HS2,
    // PC2, I2, Decroch2, func, finv1, finvp1, finv2, finvp2, PFirst, MaxStep,
    // Fleche, TolGuide, First, Last, Soldep, Inside, Appro, Forward, RecP1,
    // RecRst1, RecP2, RecRst2);` (`cxx:2142-2174`, `:2222-2252`).
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.compute_data_rst_rst(
            data,
            guide,
            &mut lin,
            s1,
            pc1,
            Rc::clone(i1),
            decroch1,
            s2,
            pc2,
            Rc::clone(i2),
            decroch2,
            &mut *func.borrow_mut(),
            finv1,
            finv_p1,
            finv2,
            finv_p2,
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
            rec_p1,
            rec_rst1,
            rec_p2,
            rec_rst2,
        )
    };
    // `if (!done) { Spine->SetErrorStatus(ChFiDS_WalkingFailure); throw
    // Standard_Failure("PerformSurf : Failed processing!"); }`.
    if !done {
        walking_failure(spine);
        return false;
    }

    // `done = CompleteData(Data, func, lin, HS1, HS2, Or);`
    let Some(line) = lin.clone() else {
        return false;
    };
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.complete_data_rst_rst(
            data,
            Rc::clone(&func),
            Rc::clone(&line),
            hs1,
            Some(hs2),
            or,
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

    // `if (maybesingular) { SplitSurf(SeqData, lin); }` (`cxx:2266-2269`).
    if maybe_singular {
        let guard = line.borrow();
        if base.split_surf(seq_data, &guard).is_err() {
            return false;
        }
    }
    true
}

/// `PerformSurf(...)` edge / edge (`cxx:2065-2270`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_surf_contact(
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
    or1: Orientation,
    hs2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    pc2: Option<&Curve2dAdaptor>,
    hsref2: &SurfaceAdaptor,
    pcref2: &Curve2dAdaptor,
    decroch2: &mut bool,
    or2: Orientation,
    max_step: f64,
    fleche: f64,
    tol_guide: f64,
    first: &mut f64,
    last: &mut f64,
    inside: bool,
    approx: bool,
    forward: bool,
    rec_p1: bool,
    rec_rst1: bool,
    rec_p2: bool,
    rec_rst2: bool,
    soldep: &MathVector,
) -> bool {
    // `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw
    // Standard_ConstructionError("PerformSurf : this is not the spine of a
    // fillet");` (`cxx:2102-2106`).
    let Some(spine_handle) = fil_spine(spine) else {
        return false;
    };
    // `double PFirst = First;` (`cxx:2109`).
    let p_first = *first;

    let (Some(s1_surf), Some(s2_surf), Some(pc1_adaptor), Some(pc2_adaptor)) = (
        hs1.surface().cloned(),
        hs2.surface().cloned(),
        pc1.cloned(),
        pc2.cloned(),
    ) else {
        return false;
    };
    let (Some(pc1_curve), Some(pc2_curve)) = (
        pc1_adaptor.curve().cloned(),
        pc2_adaptor.curve().cloned(),
    ) else {
        return false;
    };
    let (Some(ref_surf1), Some(ref_curve1), Some(ref_surf2), Some(ref_curve2)) = (
        hsref1.surface().cloned(),
        pcref1.curve().cloned(),
        hsref2.surface().cloned(),
        pcref2.curve().cloned(),
    ) else {
        return false;
    };

    let constant = match spine_handle.read() {
        Ok(guard) => guard.is_fil_constant(),
        Err(_) => return false,
    };

    // `ch1` / `ch2` (`cxx:2125-2139`, `:2204-2217`) are fixed before the walk
    // and do not depend on the spine law.
    let ch1 = contact_choix(or1);
    let ch2 = contact_choix(or2);

    if constant {
        // `BRepBlend_RstRstConstRad func(HS1, PC1, HS2, PC2, HGuide);
        // func.Set(HSref1, PCref1, HSref2, PCref2); HC1->Load(PC1, HS1);
        // HC2->Load(PC2, HS2);` (`cxx:2112-2118`).
        let mut func = SurfRstRstConstRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            Arc::clone(&s2_surf),
            Arc::clone(&pc2_curve),
            guide_as_curve(guide),
        );
        func.set_refs(
            Arc::clone(&ref_surf1),
            Arc::clone(&ref_curve1),
            Arc::clone(&ref_surf2),
            Arc::clone(&ref_curve2),
        );
        let hc1: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc1_curve), Arc::clone(&s1_surf)));
        let hc2: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc2_curve), Arc::clone(&s2_surf)));
        // `BRepBlend_SurfCurvConstRadInv finv1(HSref1, HC2, HGuide);
        // BRepBlend_CurvPointRadInv finvp1(HGuide, HC2);
        // BRepBlend_SurfCurvConstRadInv finv2(HSref2, HC1, HGuide);
        // BRepBlend_CurvPointRadInv finvp2(HGuide, HC1);` (`cxx:2119-2122`).
        let mut finv1 = SurfCurvConstRadInv::new(
            Arc::clone(&ref_surf1),
            Arc::clone(&hc2),
            guide_as_curve(guide),
        );
        let mut finvp1 = CurvPointRadInv::new(guide_as_curve(guide), Arc::clone(&hc2));
        let mut finv2 = SurfCurvConstRadInv::new(
            Arc::clone(&ref_surf2),
            Arc::clone(&hc1),
            guide_as_curve(guide),
        );
        let mut finvp2 = CurvPointRadInv::new(guide_as_curve(guide), Arc::clone(&hc1));
        // `finv1.Set(PCref1); finv2.Set(PCref2);` (`cxx:2124-2125`).
        finv1.set(Some(Arc::clone(&ref_curve1)));
        finv2.set(Some(Arc::clone(&ref_curve2)));

        // `double rad = fsp->Radius();` (`cxx:2126`).
        let Some(rad) = fil_radius(spine_handle) else {
            return false;
        };
        // `finv1.Set(rad, ch1); finvp1.Set(Choix); finv2.Set(rad, ch2);
        // finvp2.Set(Choix); func.Set(rad, Choix); func.Set(myShape);`
        // (`cxx:2137-2141`).
        finv1.set_radius(rad, ch1);
        finvp1.set_choix(choix);
        finv2.set_radius(rad, ch2);
        finvp2.set_choix(choix);
        func.set_radius(rad, choix);
        func.set_section_shape(shape);

        run_arm_contact(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s1_surf),
            pc1_adaptor,
            i1,
            decroch1,
            Arc::clone(&s2_surf),
            pc2_adaptor,
            i2,
            decroch2,
            hs1,
            hs2,
            func,
            &mut finv1,
            &mut finvp1,
            &mut finv2,
            &mut finvp2,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p1,
            rec_rst1,
            rec_p2,
            rec_rst2,
            soldep,
        )
    } else {
        // `BRepBlend_RstRstEvolRad func(HS1, PC1, HS2, PC2, HGuide,
        // fsp->Law(HGuide));` (`cxx:2189`).
        let Some(law) = fil_law(spine_handle, guide) else {
            return false;
        };
        let mut func = SurfRstRstEvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&pc1_curve),
            Arc::clone(&s2_surf),
            Arc::clone(&pc2_curve),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        func.set_refs(
            Arc::clone(&ref_surf1),
            Arc::clone(&ref_curve1),
            Arc::clone(&ref_surf2),
            Arc::clone(&ref_curve2),
        );
        let hc1: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc1_curve), Arc::clone(&s1_surf)));
        let hc2: Arc<dyn Curve> =
            Arc::new(CurveOnSurface::new(Arc::clone(&pc2_curve), Arc::clone(&s2_surf)));
        // `BRepBlend_SurfCurvEvolRadInv finv1(HSref1, HC2, HGuide,
        // fsp->Law(HGuide)); BRepBlend_CurvPointRadInv finvp1(HGuide, HC2);
        // BRepBlend_SurfCurvEvolRadInv finv2(HSref2, HC1, HGuide,
        // fsp->Law(HGuide)); BRepBlend_CurvPointRadInv finvp2(HGuide, HC1);`
        // (`cxx:2197-2200`).
        let mut finv1 = SurfCurvEvolRadInv::new(
            Arc::clone(&ref_surf1),
            Arc::clone(&hc2),
            guide_as_curve(guide),
            Arc::clone(&law),
        );
        let mut finvp1 = CurvPointRadInv::new(guide_as_curve(guide), Arc::clone(&hc2));
        let mut finv2 = SurfCurvEvolRadInv::new(
            Arc::clone(&ref_surf2),
            Arc::clone(&hc1),
            guide_as_curve(guide),
            law,
        );
        let mut finvp2 = CurvPointRadInv::new(guide_as_curve(guide), Arc::clone(&hc1));
        // `finv1.Set(PCref1); finv2.Set(PCref2);` (`cxx:2202-2203`).
        finv1.set(Some(Arc::clone(&ref_curve1)));
        finv2.set(Some(Arc::clone(&ref_curve2)));

        // `finv1.Set(ch1); finvp1.Set(Choix); finv2.Set(ch2);
        // finvp2.Set(Choix); func.Set(Choix); func.Set(myShape);`
        // (`cxx:2214-2220`).
        finv1.set_choix(ch1);
        finvp1.set_choix(choix);
        finv2.set_choix(ch2);
        finvp2.set_choix(choix);
        func.set_choix(choix);
        func.set_section_shape(shape);

        run_arm_contact(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            Arc::clone(&s1_surf),
            pc1_adaptor,
            i1,
            decroch1,
            Arc::clone(&s2_surf),
            pc2_adaptor,
            i2,
            decroch2,
            hs1,
            hs2,
            func,
            &mut finv1,
            &mut finvp1,
            &mut finv2,
            &mut finvp2,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p1,
            rec_rst1,
            rec_p2,
            rec_rst2,
            soldep,
        )
    }
}
