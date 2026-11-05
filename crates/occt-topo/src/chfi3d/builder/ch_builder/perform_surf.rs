//! `ChFi3d_ChBuilder::PerformSurf` — the face / face overload.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:1618-1853`
//! (OCCT 8.0.0). It walks the chamfer band between two faces and completes the
//! `ChFiDS_SurfData` of the stripe.
//!
//! `chsp->IsChamfer()` (`cxx:1663`, `:1721`, `:1799`) splits the body in the
//! same three arms as `SimulSurf` (see
//! [`super::simul_surf`]), with the same blend functions and the same `Set`
//! arguments; `chsp->Mode()` (`cxx:1666`, `:1727`) picks the classic or the
//! constant-throat family, and the `TwoDist` / non-classic arm again substitutes
//! the **offset** guideline for `HGuide` (`cxx:1736-1756`).
//!
//! All three arms then run the same `ComputeData` (`ChFi3d_Builder_2.cxx`) /
//! `CompleteData` pair — that tail is [`run_arm`]. Unlike the fillet version
//! (`ChFi3d_FilBuilder.cxx:1712-1715`) there is **no** `SplitSurf` here: the
//! chamfer `PerformSurf` ends after `CompleteData` (`cxx:1684-1718`,
//! `:1761-1798`, `:1811-1849`).
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("PerformSurf : this is not the spine of
//!   a chamfer")` (`cxx:1645-1648`) and `throw Standard_Failure("PerformSurf :
//!   Fail of approximation!")` (`cxx:1713-1717`, `:1793-1797`, `:1844-1848`)
//!   are the trait's `false`: the caller's recovery path reacts to "did not
//!   compute" exactly as it does to the exception,
//! - `Or = S1->Face().Orientation()` (`cxx:1651`) raises `Standard_NullObject`
//!   for an adaptor without a face; the port answers `false` on that arm,
//! - the four `Gd` / `Gf` flags (`cxx:1649`) are uninitialised in OCCT before
//!   `ComputeData` writes them; the port starts them at `false`,
//! - the offset guideline is cloned out of `ChFiDS_Spine::offset_el_spines`,
//!   so its reframing stays local — the same `Arc` / `handle` divergence as
//!   [`super::simul_surf`],
//! - `SeqData(1)` is `seq_data[0]`; an empty sequence is OCCT's
//!   `Standard_OutOfRange`, answered `false`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::blend::{
    BlendFuncInv, BlendFunction, ChAsym, ChAsymInv, ChamfInv, Chamfer, ConstThroat, ConstThroatInv,
    ConstThroatWithPenetration, ConstThroatWithPenetrationInv,
};
use crate::brep_blend::Line;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi_ds::{ChamfMethod, ChamfMode, ElSpine, SpineClass, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::simul_surf::guide_as_curve;

/// `chsp = down_cast<ChFiDS_ChamfSpine>(Spine); if (chsp.IsNull()) throw`
/// (`cxx:1645-1648`).
fn chamf_spine(spine: Option<&SpineHandle>) -> Option<&SpineHandle> {
    let spine = spine?;
    let is_chamf = spine.read().ok()?.class() == SpineClass::Chamf;
    if is_chamf {
        Some(spine)
    } else {
        None
    }
}

/// `PerformSurf(SeqData, HGuide, Spine, Choix, S1, I1, S2, I2, MaxStep, Fleche,
/// TolGuide, First, Last, Inside, Appro, Forward, RecOnS1, RecOnS2, Soldep,
/// intf, intl)` (`cxx:1618-1853`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_surf_face_face(
    base: &mut ChFi3dBuilder,
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
    // `occ::handle<ChFiDS_SurfData> Data = SeqData(1);` (`cxx:1641`).
    if seq_data.is_empty() {
        return false;
    }
    let Some(spine_handle) = chamf_spine(spine) else {
        return false;
    };

    // `TopAbs_Orientation Or = S1->Face().Orientation();` (`cxx:1651`).
    let Some(or) = s1.face().map(|f| f.0.orientation()) else {
        return false;
    };
    // `double PFirst = First;` (`cxx:1652`), before the two reframings.
    let p_first = *first;
    // `if (intf) { First = chsp->FirstParameter(1); }` (`cxx:1653-1656`).
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
    // `if (intl) { Last = chsp->LastParameter(chsp->NbEdges()); }`
    // (`cxx:1657-1660`).
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

    // `chsp->IsChamfer()` trials and the matching `chsp->Mode()` reads
    // (`cxx:1663`, `:1666`, `:1721`, `:1723-1724`, `:1727`, `:1799-1801`).
    let (method, mode, chamf_dis1, chamf_dis2, chamf_angle) = {
        let guard = match spine_handle.read() {
            Ok(guard) => guard,
            Err(_) => return false,
        };
        let method = guard.is_chamfer();
        let mode = guard.mode();
        let (d1, d2, angle) = match method {
            ChamfMethod::Sym => match guard.get_dist() {
                Some(dis) => (dis, 0.0, 0.0),
                None => return false,
            },
            ChamfMethod::TwoDist => match guard.dists() {
                Some((d1, d2)) => (d1, d2, 0.0),
                None => return false,
            },
            ChamfMethod::DistAngle => match guard.get_dist_angle() {
                Some((d1, angle)) => (d1, 0.0, angle),
                None => return false,
            },
        };
        (method, mode, d1, d2, angle)
    };

    match method {
        ChamfMethod::Sym => {
            // `if (chsp->Mode() == ChFiDS_ClassicChamfer) { Chamfer; ChamfInv;
            // } else { ConstThroat; ConstThroatInv; }`
            // (`cxx:1666-1678`), then `double dis; chsp->GetDist(dis);
            // Set(dis, dis, Choix)` on both (`cxx:1679-1680`).
            let dis = chamf_dis1;
            if mode == ChamfMode::ClassicChamfer {
                let mut func = Chamfer::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    Arc::clone(&guide_curve),
                );
                let mut finv = ChamfInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                func.set_distances(dis, dis, choix);
                finv.set_distances(dis, dis, choix);
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
                let mut func = ConstThroat::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    Arc::clone(&guide_curve),
                );
                let mut finv = ConstThroatInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                func.set_throat(dis, dis, choix);
                finv.set_throat(dis, dis, choix);
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
        ChamfMethod::TwoDist => {
            // `double d1, d2; chsp->Dists(d1, d2);` (`cxx:1723-1724`).
            let (d1, d2) = (chamf_dis1, chamf_dis2);
            if mode == ChamfMode::ClassicChamfer {
                // `Chamfer; ChamfInv; Set(d1, d2, Choix)` (`cxx:1727-1734`).
                let mut func = Chamfer::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    Arc::clone(&guide_curve),
                );
                let mut finv = ChamfInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    guide_as_curve(guide),
                );
                func.set_distances(d1, d2, choix);
                finv.set_distances(d1, d2, choix);
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
                // `OffsetHGuide` lookup over `ChangeElSpines()` /
                // `ChangeOffsetElSpines()` (`cxx:1736-1745`); the null arm is
                // the "Construction of offset guide failed!" print whose
                // exception is commented out (`cxx:1747-1751`), answered here
                // as `false`.
                let offset_arc: Arc<ElSpine> = {
                    let guard = match spine_handle.read() {
                        Ok(guard) => guard,
                        Err(_) => return false,
                    };
                    match guard.offset_el_spine(guide) {
                        Some(offset) => Arc::clone(offset),
                        None => return false,
                    }
                };
                // `ConstThroatWithPenetration(S1, S2, OffsetHGuide);
                // ConstThroatWithPenetrationInv(S1, S2, OffsetHGuide);
                // double Throat = std::max(d1, d2); Set(Throat, Throat,
                // Choix);` (`cxx:1752-1757`).
                let offset_curve = guide_as_curve(&offset_arc);
                let mut func = ConstThroatWithPenetration::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    offset_curve,
                );
                let mut finv = ConstThroatWithPenetrationInv::new(
                    Arc::clone(&s1_surf),
                    Arc::clone(&s2_surf),
                    offset_arc,
                );
                let throat = d1.max(d2);
                func.set_throat(throat, throat, choix);
                finv.set_throat(throat, throat, choix);
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
        ChamfMethod::DistAngle => {
            // `double d1, angle; chsp->GetDistAngle(d1, angle); int Ch =
            // Choix; BRepBlend_ChAsym Func(S1, S2, HGuide); ...
            // Func.Set(d1, angle, Ch); FInv.Set(d1, angle, Ch);`
            // (`cxx:1800-1808`).
            let (d1, angle) = (chamf_dis1, chamf_angle);
            let mut func = ChAsym::new(
                Arc::clone(&s1_surf),
                Arc::clone(&s2_surf),
                Arc::clone(&guide_curve),
            );
            let mut finv = ChAsymInv::new(
                Arc::clone(&s1_surf),
                Arc::clone(&s2_surf),
                guide_as_curve(guide),
            );
            func.set_dist_angle(d1, angle, choix);
            finv.set_dist_angle(d1, angle, choix);
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
}

/// The `ComputeData` / `CompleteData` tail shared by the three arms of
/// [`perform_surf_face_face`] (`cxx:1684-1718`, `:1761-1798`, `:1811-1849`).
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

    // `bool gd1, gd2, gf1, gf2;` (`cxx:1649`): `ComputeData` writes all four
    // before `CompleteData` reads them.
    let mut gd1 = false;
    let mut gd2 = false;
    let mut gf1 = false;
    let mut gf2 = false;

    // `done = ComputeData(Data, HGuide, Spine, lin, S1, I1, S2, I2, Func,
    // FInv, PFirst, MaxStep, Fleche, TolGuide, First, Last, Inside, Appro,
    // Forward, Soldep, intf, intl, gd1, gd2, gf1, gf2, RecOnS1, RecOnS2);`
    // (`cxx:1684-1709`).
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
    // (`cxx:1710-1713`).
    if !done {
        return false;
    }

    // `done = CompleteData(Data, *pFunc, lin, S1, S2, Or, gd1, gd2, gf1, gf2);`
    // (`cxx:1714`). `Reversed` keeps its `false` default.
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
    // `if (!done) { throw Standard_Failure("PerformSurf : Fail of
    // approximation!"); }` (`cxx:1715-1718`) — see the module note.
    matches!(done, Ok(true))
}
