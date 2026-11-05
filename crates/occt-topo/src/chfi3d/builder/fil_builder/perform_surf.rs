//! `ChFi3d_FilBuilder::PerformSurf` — the face / face overload.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:1538-1717`
//! (OCCT 8.0.0). It walks the band between two faces and stores one
//! `ChFiDS_SurfData`; `fsp->IsConstant()` (`cxx:1581`) splits the body between
//! `BRepBlend_ConstRad` / `BRepBlend_ConstRadInv` and `BRepBlend_EvolRad` /
//! `BRepBlend_EvolRadInv`. Both arms then run the same
//! `ComputeData` / `CompleteData` pair and the same trailing test
//! `Func.GetMinimalDistance() <= 100 * tolapp3d` that triggers `SplitSurf`.
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("PerformSurf : this is not the spine of
//!   a fillet")` (`cxx:1565-1568`) is the trait's `false`,
//! - `throw Standard_Failure("PerformSurf : Failed approximation!")`
//!   (`cxx:1641-1644`, `:1706-1709`) is also `false`: the trait answers "did
//!   not compute", which is what the caller's recovery path reacts to,
//! - `Func.GetMinimalDistance()` answers a `Result` in the port; the `Err` arm
//!   is treated as "not singular", which is where `Blend_Function` leaves
//!   `distmin` until a section lowers it,
//! - `Or = S1->Face().Orientation()` raises `Standard_NullObject` in OCCT for
//!   an adaptor without a face; the port answers `false` on the same arm,
//! - the two arms share [`run_arm`] (the `ComputeData` / `CompleteData` /
//!   `SplitSurf` tail is byte-identical in OCCT, `cxx:1593-1645` vs
//!   `:1659-1710`), so the blend function is passed as the generic `F` OCCT
//!   passes as `Blend_Function&`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::law::{LawFunction, LawHandle};
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{BlendFunction, BlendFuncInv, ConstRad, ConstRadInv, EvolRad, EvolRadInv};
use crate::brep_blend::Line;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi_ds::{ElSpine, SpineClass, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::simul_surf::guide_as_curve;

/// `HS->Face().Orientation()`.
pub(super) fn face_orientation(hs: &SurfaceAdaptor) -> Option<Orientation> {
    hs.face().map(|f| f.0.orientation())
}

/// `fsp = down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw`
/// (`cxx:1564-1568`).
pub(super) fn fil_spine(spine: Option<&SpineHandle>) -> Option<&SpineHandle> {
    let spine = spine?;
    let is_fil = spine.read().ok()?.class() == SpineClass::Fil;
    if is_fil {
        Some(spine)
    } else {
        None
    }
}

/// `PerformSurf(SeqData, HGuide, Spine, Choix, S1, I1, S2, I2, MaxStep, Fleche,
/// TolGuide, First, Last, Inside, Appro, Forward, RecOnS1, RecOnS2, Soldep,
/// intf, intl)` (`ChFi3d_FilBuilder.cxx:1538-1717`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_surf_face_face(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    seq_data: &mut Vec<Arc<SurfData>>,
    guide: &mut Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    choix: i32,
    s1: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    s2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    max_step: f64,
    fleche: f64,
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
) -> bool {
    // `occ::handle<ChFiDS_SurfData> Data = SeqData(1);` (`cxx:1563`).
    if seq_data.is_empty() {
        return false;
    }
    let Some(spine_handle) = fil_spine(spine) else {
        return false;
    };

    // `TopAbs_Orientation Or = S1->Face().Orientation();` (`cxx:1571`).
    let Some(or) = face_orientation(s1) else {
        return false;
    };
    // `double PFirst = First;` (`cxx:1572`), before the two reframings.
    let p_first = *first;
    // `if (intf) { First = fsp->FirstParameter(1); }` (`cxx:1573-1576`).
    if *intf != 0 {
        let param = match spine_handle.read() {
            Ok(guard) => guard.first_parameter_of(1),
            Err(_) => None,
        };
        match param {
            Some(param) => *first = param,
            None => return false,
        }
    }
    // `if (intl) { Last = fsp->LastParameter(fsp->NbEdges()); }` (`cxx:1577-1580`).
    if *intl != 0 {
        let param = match spine_handle.read() {
            Ok(guard) => guard.last_parameter_of(guard.nb_edges() as i32),
            Err(_) => None,
        };
        match param {
            Some(param) => *last = param,
            None => return false,
        }
    }

    let Some(s1_surf) = s1.surface().cloned() else {
        return false;
    };
    let Some(s2_surf) = s2.surface().cloned() else {
        return false;
    };
    let guide_curve: Arc<dyn Curve> = guide_as_curve(guide);

    let constant = match spine_handle.read() {
        Ok(guard) => guard.is_fil_constant(),
        Err(_) => return false,
    };

    if constant {
        // `BRepBlend_ConstRad Func(S1, S2, HGuide); BRepBlend_ConstRadInv
        // FInv(S1, S2, HGuide); Func.Set(fsp->Radius(), Choix); FInv.Set(...);
        // Func.Set(myShape);` (`cxx:1583-1587`).
        let radius = match spine_handle.read() {
            Ok(guard) => guard.fil_radius(),
            Err(_) => None,
        };
        let Some(radius) = radius else {
            return false;
        };
        let mut func = ConstRad::new(Arc::clone(&s1_surf), Arc::clone(&s2_surf), guide_curve);
        let mut finv = ConstRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            guide_as_curve(guide),
        );
        func.set_radius(radius, choix);
        finv.set_radius(radius, choix);
        func.set_section_shape(shape);
        run_arm(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            s1,
            i1,
            s2,
            i2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_on_s1,
            rec_on_s2,
            soldep,
            intf,
            intl,
            or,
            p_first,
            func,
            &mut finv,
        )
    } else {
        // `BRepBlend_EvolRad Func(S1, S2, HGuide, fsp->Law(HGuide));
        // BRepBlend_EvolRadInv FInv(...); Func.Set(Choix); FInv.Set(Choix);
        // Func.Set(myShape);` (`cxx:1649-1653`).
        let law = match spine_handle.read() {
            Ok(guard) => guard.fil_law(guide),
            Err(_) => None,
        };
        let Some(law) = law else {
            return false;
        };
        let law_f: Arc<dyn LawFunction> = Arc::new(LawHandle::new(Arc::clone(&law)));
        let law_inv: Arc<dyn LawFunction> = Arc::new(LawHandle::new(law));
        let mut func = EvolRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            guide_curve,
            law_f,
        );
        let mut finv = EvolRadInv::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            guide_as_curve(guide),
            law_inv,
        );
        func.set_choix(choix);
        finv.set_choix(choix);
        func.set_section_shape(shape);
        run_arm(
            base,
            seq_data,
            guide,
            Some(spine_handle),
            s1,
            i1,
            s2,
            i2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_on_s1,
            rec_on_s2,
            soldep,
            intf,
            intl,
            or,
            p_first,
            func,
            &mut finv,
        )
    }
}

/// The `ComputeData` / `CompleteData` / `SplitSurf` tail shared by the two arms
/// of [`perform_surf_face_face`] (`cxx:1593-1645` and `:1659-1710`).
#[allow(clippy::too_many_arguments)]
fn run_arm<F: BlendFunction + 'static>(
    base: &mut ChFi3dBuilder,
    seq_data: &mut Vec<Arc<SurfData>>,
    guide: &mut Arc<ElSpine>,
    spine: Option<&SpineHandle>,
    s1: &SurfaceAdaptor,
    i1: &Rc<RefCell<TopolTool>>,
    s2: &SurfaceAdaptor,
    i2: &Rc<RefCell<TopolTool>>,
    max_step: f64,
    fleche: f64,
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
    or: Orientation,
    p_first: f64,
    func: F,
    finv: &mut dyn BlendFuncInv,
) -> bool {
    let Some(s1_surf) = s1.surface().cloned() else {
        return false;
    };
    let Some(s2_surf) = s2.surface().cloned() else {
        return false;
    };
    // `BRepBlend_AppFunc` keeps the `Buff_Curve` of the line and the function
    // by reference in OCCT; the port shares it through `Rc<RefCell<..>>`.
    let func = Rc::new(RefCell::new(func));
    let mut lin: Option<Rc<RefCell<Line>>> = None;

    // The four `Gd` / `Gf` flags `ComputeData` answers (`cxx:1569`); OCCT
    // leaves them uninitialised before the call and `CompleteData` reads them
    // after it, so the initial value is never observed.
    let mut gd1 = false;
    let mut gd2 = false;
    let mut gf1 = false;
    let mut gf2 = false;

    // `done = ComputeData(Data, HGuide, Spine, lin, S1, I1, S2, I2, Func,
    // FInv, PFirst, MaxStep, Fleche, TolGuide, First, Last, Inside, Appro,
    // Forward, Soldep, intf, intl, gd1, gd2, gf1, gf2, RecOnS1, RecOnS2);`
    // (`cxx:1593-1620`).
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.compute_data_face_face(
            data,
            guide,
            spine,
            &mut lin,
            s1_surf,
            Rc::clone(i1),
            s1.face(),
            s2_surf,
            Rc::clone(i2),
            s2.face(),
            &mut *func.borrow_mut(),
            finv,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            soldep,
            intf,
            intl,
            &mut gd1,
            &mut gd2,
            &mut gf1,
            &mut gf2,
            rec_on_s1,
            rec_on_s2,
        )
    };
    // `if (!done) { return false; // recovery is possible PMN 14/05/1998 }`
    // (`cxx:1626-1629`).
    if !done {
        return false;
    }

    // `done = CompleteData(Data, Func, lin, S1, S2, Or, gd1, gd2, gf1, gf2);`
    // (`cxx:1635`). `Reversed` keeps its `false` default.
    let Some(lin) = lin else {
        return false;
    };
    let done = {
        let data = Arc::make_mut(&mut seq_data[0]);
        base.complete_data_app_func(
            data,
            Rc::clone(&func),
            Rc::clone(&lin),
            s1,
            Some(s2),
            or,
            gd1,
            gd2,
            gf1,
            gf2,
            false,
        )
    };
    // `if (!done) { throw Standard_Failure("PerformSurf : Failed
    // approximation!"); }` (`cxx:1641-1644`) — see the module note.
    if !matches!(done, Ok(true)) {
        return false;
    }

    // `maybesingular = (Func.GetMinimalDistance() <= 100 * tolapp3d);`
    // (`cxx:1645`).
    let maybe_singular = match func.borrow().get_minimal_distance() {
        Ok(distance) => distance <= 100.0 * base.tolapp3d,
        Err(_) => false,
    };

    // `if (maybesingular) { SplitSurf(SeqData, lin); }` (`cxx:1712-1715`).
    // `SplitSurf` only raises on a corrupt data structure, which OCCT lets
    // escape out of `PerformSurf`; the trait's `false` is that abort.
    if maybe_singular {
        let guard = lin.borrow();
        if base.split_surf(seq_data, &guard).is_err() {
            return false;
        }
    }
    true
}
