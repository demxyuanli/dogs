//! `ChFi3d_FilBuilder::PerformThreeCorner` - the corner surface when the three
//! stripes leave neither an analytic torus nor a sphere.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder_C3.cxx:739-1015`
//! (OCCT 8.0.0):
//!
//! - `cxx:744-908` - the guideline `ChFi3d_CircularSpine` / `ChFi3d_Spine`, the
//!   `BRepBlend_ConstRad` or `BRepBlend_EvolRad` walk of `ComputeData` between
//!   the pivot face `Fac` and the pivot surface `Surf`, then `CompleteData`,
//! - `cxx:910-1015` - the `GeomFill_ConstrainedFilling` corner between the
//!   facet `Bfac` (absent when `c1pointu`), the pivot `Bpiv` and the two stripe
//!   surfaces `Bdeb` / `Bfin`, then `CompleteData` on the reversed patch.
//!
//! Divergences that are not heuristics:
//! - OCCT's `done` is [`ChFi3dBuilder::done`] and `filling` is
//!   [`ThreeCorner::filling`]; both blocks belong to one C++ body, so this
//!   function runs them in that order and answers the final `done`,
//! - `ChFi3d_CircularSpine` answers `Result<Option<Geom_Circle>>` in the port,
//!   and its `spinecoin.IsNull()` arm (`cxx:766-773`) is the `None`,
//! - `ChFiDS_Spine NullSpine` (`cxx:783`) is the null handle OCCT hands to
//!   `ComputeData`, so the port passes `None` for the spine,
//! - `ComputeData` / `CompleteData` write `coin` in place, and the `SurfData`
//!   handles OCCT keeps (`fddeb` / `fdfin` / `fdpiv`) are read through
//!   [`super::sd`], which answers the same values without the aliasing,
//! - the `DStr.Curve(kkk).Curve()->D1(...)` lookups of `cxx:937-946` feed only
//!   the commented-out `ChFi3d_mkbound` call of `cxx:947-948`; they are kept, so
//!   a missing interference or curve fails here as it would in OCCT,
//! - `Geom_Surface::VReverse()` of `cxx:997` is an in-place mutation of the
//!   patch `GeomFill_ConstrainedFilling::Surface()` returned; the port asks the
//!   `Surface` trait for the reversed copy, which is what the only later reader
//!   (`CompleteData`) sees.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_geom::curve::Curve;
use occt_geom::geom_fill::boundary::Boundary;
use occt_geom::geom_fill::ConstrainedFilling;
use occt_geom::law::{LawFunction, LawS};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::MathVector;

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{ConstRad, ConstRadInv, EvolRad, EvolRadInv};
use crate::brep_blend::line::Line;
use crate::chfi3d::builder_0_helpers::{mkbound_curve_segment, mkbound_pfac2d, mkbound_segment};
use crate::chfi3d::builder_helpers::{circular_spine, spine as bezier_spine};
use crate::chfi_ds::{ElSpine, StripeHandle, SurfData};

use super::super::simul_surf::guide_as_curve;
use super::super::ChFi3dBuilder;
use super::{sd, CornerFrame, ThreeCorner};

/// `DStr.Curve(kkk).Curve()->D1(param, P, V)` (`cxx:940`, `:944`).
fn curve_d1(base: &ChFi3dBuilder, index: i32, param: f64) -> Result<(GpPnt, GpVec), String> {
    let curve = base.my_ds.curve(index).curve().ok_or_else(|| {
        "PerformThreeCorner: DStr.Curve(kkk).Curve() is null".to_string()
    })?;
    Ok(curve.d1(param))
}

/// `CD[k]->SetOfSurfData()->Value(i)->Interference(jf).PCurveOnSurf()->Value(p)`.
fn p_curve_on_surf_value(
    stripe: &StripeHandle,
    index: i32,
    jf: i32,
    param: f64,
) -> Result<GpPnt2d, String> {
    let data = sd(stripe, index)?;
    let curve = data.interference(jf).p_curve_on_surf().ok_or_else(|| {
        "PerformThreeCorner: the interference has no p-curve on the surface".to_string()
    })?;
    Ok(curve.value(param))
}

/// `DStr.Surface(CD[k]->SetOfSurfData()->Value(i)->Surf()).Surface()`
/// (`cxx:979-983`).
fn ds_surface(
    base: &ChFi3dBuilder,
    stripe: &StripeHandle,
    index: i32,
) -> Result<Arc<dyn Surface>, String> {
    let surf = sd(stripe, index)?.surf();
    base.my_ds
        .ds()
        .surface(surf)
        .surface()
        .cloned()
        .ok_or_else(|| "PerformThreeCorner: DStr.Surface(Surf()).Surface() is null".to_string())
}

/// `if (!done) { if (!filling) { ... } if (filling) { ... } }`
/// (`cxx:741-1016`). Answers the resulting `done`.
pub(super) fn build_corner(
    base: &mut ChFi3dBuilder,
    shape: BlendFuncSectionShape,
    c: &mut ThreeCorner,
    frame: &CornerFrame,
    coin: &mut SurfData,
) -> Result<bool, String> {
    // `if (!filling)` (`cxx:744-908`).
    if !c.filling {
        // `double radpondere = (Rdp + Rfp) / 2.; double locfleche = fleche;`
        // (`cxx:759-760`).
        let radpondere = (frame.rdp + frame.rfp) / 2.0;
        let fleche = base.fleche();
        let mut locfleche = fleche;

        // `double WFirst, WLast; occ::handle<Geom_Curve> spinecoin =
        // ChFi3d_CircularSpine(WFirst, WLast, Pdeb, Vdeb, Pfin, Vfin,
        // radpondere);` (`cxx:762-764`).
        let mut w_first = 0.0;
        let mut w_last = 0.0;
        let circle = circular_spine(
            &mut w_first,
            &mut w_last,
            &frame.pdeb,
            &frame.vdeb,
            &frame.pfin,
            &frame.vfin,
            radpondere,
        )?;
        let spinecoin: Arc<dyn Curve> = match circle {
            Some(circle) => {
                // `locfleche = radpondere * (WLast - WFirst) * fleche;`
                // (`cxx:774-777`).
                locfleche = radpondere * (w_last - w_first) * fleche;
                Arc::new(circle)
            }
            None => {
                // `spinecoin = ChFi3d_Spine(Pdeb, Vdeb, Pfin, Vfin,
                // radpondere); WFirst = 0.; WLast = 1.;` (`cxx:769-772`).
                let bezier = bezier_spine(
                    &frame.pdeb,
                    &frame.vdeb,
                    &frame.pfin,
                    &frame.vfin,
                    radpondere,
                );
                w_first = 0.0;
                w_last = 1.0;
                Arc::new(bezier)
            }
        };

        // `double pasmax = (WLast - WFirst) * 0.05; cornerspine = new
        // ChFiDS_ElSpine(); cornerspine->SetCurve(spinecoin);
        // cornerspine->FirstParameter(WFirst - pasmax);
        // cornerspine->LastParameter(WLast + pasmax);` (`cxx:778-781`).
        let pasmax = (w_last - w_first) * 0.05;
        let mut cornerspine = ElSpine::new();
        cornerspine.set_curve(spinecoin);
        cornerspine.set_first_parameter(w_first - pasmax);
        cornerspine.set_last_parameter(w_last + pasmax);
        let mut cornerspine = Arc::new(cornerspine);

        // `math_Vector Soldep(1, 4); Soldep(1) = pfac1.X(); ... Soldep(4) =
        // psurf1.Y();` (`cxx:788-792`).
        let mut soldep = MathVector::new(1, 4);
        soldep.set_value(1, frame.pfac1.x());
        soldep.set_value(2, frame.pfac1.y());
        soldep.set_value(3, frame.psurf1.x());
        soldep.set_value(4, frame.psurf1.y());

        // `double ffi = WFirst, lla = WLast + pasmax;` (`cxx:797`).
        let mut ffi = w_first;
        let mut lla = w_last + pasmax;

        let s1 = frame
            .fac
            .surface()
            .cloned()
            .ok_or_else(|| "PerformThreeCorner: Fac has no surface".to_string())?;
        let s2 = frame
            .surf
            .surface()
            .cloned()
            .ok_or_else(|| "PerformThreeCorner: Surf has no surface".to_string())?;

        // `int intf = 3, intl = 3;` (`cxx:808`, `:861`).
        let mut intf = 3;
        let mut intl = 3;
        // `bool Gd1, Gd2, Gf1, Gf2;` (`cxx:794`): OCCT leaves them
        // uninitialised before `ComputeData` writes them.
        let mut gd1 = false;
        let mut gd2 = false;
        let mut gf1 = false;
        let mut gf2 = false;
        let mut lin: Option<Rc<RefCell<Line>>> = None;

        let done = if (frame.rdeb - frame.rfin).abs() <= base.tolapp3d {
            // `BRepBlend_ConstRad func(Fac, Surf, cornerspine); BRepBlend_
            // ConstRadInv finv(Fac, Surf, cornerspine); func.Set(Rdeb,
            // choix); func.Set(myShape); finv.Set(Rdeb, choix);`
            // (`cxx:801-806`).
            let mut func = ConstRad::new(
                Arc::clone(&s1),
                Arc::clone(&s2),
                guide_as_curve(&cornerspine),
            );
            let mut finv = ConstRadInv::new(
                Arc::clone(&s1),
                Arc::clone(&s2),
                guide_as_curve(&cornerspine),
            );
            func.set_radius(frame.rdeb, frame.choix);
            func.set_section_shape(shape);
            finv.set_radius(frame.rdeb, frame.choix);
            // `double TolGuide = cornerspine->Resolution(tolapp3d);`
            // (`cxx:807`).
            let tol_guide = cornerspine.resolution(base.tolapp3d).ok_or_else(|| {
                "PerformThreeCorner: ChFiDS_ElSpine::Resolution".to_string()
            })?;
            let func = Rc::new(RefCell::new(func));
            let mut done = base.compute_data_face_face(
                coin,
                &mut cornerspine,
                None,
                &mut lin,
                s1,
                Rc::clone(&frame.ifac),
                frame.fac.face(),
                s2,
                Rc::clone(&frame.isurf),
                frame.surf.face(),
                &mut *func.borrow_mut(),
                &mut finv,
                ffi,
                pasmax,
                locfleche,
                tol_guide,
                &mut ffi,
                &mut lla,
                false,
                false,
                true,
                &soldep,
                &mut intf,
                &mut intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                false,
                true,
            );
            // `if (done && Gf2) { done = CompleteData(coin, func, lin, Fac,
            // Surf, OFac, Gd1, false, Gf1, false); filling = !done; } else {
            // filling = true; }` (`cxx:843-850`).
            if done && gf2 {
                let line = lin.clone();
                done = match line {
                    Some(line) => base.complete_data_app_func(
                        coin,
                        Rc::clone(&func),
                        line,
                        &frame.fac,
                        Some(&frame.surf),
                        frame.o_fac,
                        gd1,
                        false,
                        gf1,
                        false,
                        false,
                    )?,
                    None => false,
                };
                c.filling = !done;
            } else {
                c.filling = true;
            }
            done
        } else {
            // `Law_S law = new Law_S(); law->Set(WFirst, Rdeb, WLast, Rfin);`
            // (`cxx:852-854`). OCCT hands the same handle to `func` and
            // `finv`.
            let mut law = LawS::new();
            law.set(w_first, frame.rdeb, w_last, frame.rfin)?;
            let law: Arc<dyn LawFunction> = Arc::new(law);
            // `BRepBlend_EvolRad func(Fac, Surf, cornerspine, law);
            // BRepBlend_EvolRadInv finv(Fac, Surf, cornerspine, law);
            // func.Set(choix); func.Set(myShape); finv.Set(choix);`
            // (`cxx:855-860`).
            let mut func = EvolRad::new(
                Arc::clone(&s1),
                Arc::clone(&s2),
                guide_as_curve(&cornerspine),
                Arc::clone(&law),
            );
            let mut finv = EvolRadInv::new(
                Arc::clone(&s1),
                Arc::clone(&s2),
                guide_as_curve(&cornerspine),
                law,
            );
            func.set_choix(frame.choix);
            func.set_section_shape(shape);
            finv.set_choix(frame.choix);
            // `double TolGuide = cornerspine->Resolution(tolapp3d);`
            // (`cxx:860`).
            let tol_guide = cornerspine.resolution(base.tolapp3d).ok_or_else(|| {
                "PerformThreeCorner: ChFiDS_ElSpine::Resolution".to_string()
            })?;
            let func = Rc::new(RefCell::new(func));
            let mut done = base.compute_data_face_face(
                coin,
                &mut cornerspine,
                None,
                &mut lin,
                s1,
                Rc::clone(&frame.ifac),
                frame.fac.face(),
                s2,
                Rc::clone(&frame.isurf),
                frame.surf.face(),
                &mut *func.borrow_mut(),
                &mut finv,
                ffi,
                pasmax,
                locfleche,
                tol_guide,
                &mut ffi,
                &mut lla,
                false,
                false,
                true,
                &soldep,
                &mut intf,
                &mut intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                false,
                true,
            );
            // `if (done && Gf2) { ... } else { filling = true; }`
            // (`cxx:894-901`).
            if done && gf2 {
                let line = lin.clone();
                done = match line {
                    Some(line) => base.complete_data_app_func(
                        coin,
                        Rc::clone(&func),
                        line,
                        &frame.fac,
                        Some(&frame.surf),
                        frame.o_fac,
                        gd1,
                        false,
                        gf1,
                        false,
                        false,
                    )?,
                    None => false,
                };
                c.filling = !done;
            } else {
                c.filling = true;
            }
            done
        };

        // `done` is the member OCCT keeps; the assignment mirrors it.
        base.done = done;
    }

    // `if (filling)` (`cxx:910-1016`).
    if c.filling {
        // `occ::handle<GeomFill_Boundary> Bdeb, Bfin, Bpiv, Bfac;
        // occ::handle<Geom2d_Curve> PCurveOnFace;` (`cxx:920-921`).
        let mut bfac: Option<Arc<dyn Boundary>> = None;
        let mut p_curve_on_face: Option<Arc<dyn Curve2d>> = None;
        if !c.c1pointu {
            // `Bfac = ChFi3d_mkbound(Fac, PCurveOnFace, sens[deb], pfac1,
            // vfac1, sens[fin], pfac2, vfac2, tolapp3d, 2.e-4);`
            // (`cxx:924-933`).
            let (pc, bound) = mkbound_pfac2d(
                &frame.fac,
                c.sens[c.deb],
                &frame.pfac1,
                &frame.vfac1,
                c.sens[c.fin],
                &frame.pfac2,
                &frame.vfac2,
                base.tolapp3d,
                2.0e-4,
            )?;
            p_curve_on_face = Some(pc);
            bfac = Some(bound);
        }

        // `kkk = CD[deb]->...->Interference(jf[deb][pivot]).LineIndex();
        // DStr.Curve(kkk).Curve()->D1(p[deb][pivot], ppbid, vp1);` and the same
        // for `fin` (`cxx:937-946`). Only the commented-out `mkbound` call of
        // `cxx:947-948` reads `vp1` / `vp2`.
        let kkk = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?
            .interference(c.jf[c.deb][c.pivot])
            .line_index();
        let _ = curve_d1(base, kkk, c.p[c.deb][c.pivot])?;
        let kkk = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?
            .interference(c.jf[c.fin][c.pivot])
            .line_index();
        let _ = curve_d1(base, kkk, c.p[c.fin][c.pivot])?;

        // `Bpiv = ChFi3d_mkbound(Surf, PCurveOnPiv, psurf1, psurf2, tolapp3d,
        // 2.e-4, false);` (`cxx:948`).
        let (p_curve_on_piv, bpiv) = mkbound_curve_segment(
            &frame.surf,
            &frame.psurf1,
            &frame.psurf2,
            base.tolapp3d,
            2.0e-4,
            false,
        )?;

        // `double pardeb2 = p[deb][pivot]; double parfin2 = p[fin][pivot];
        // if (c1pointu) { pardeb2 = p[deb][fin]; parfin2 = p[fin][deb]; }`
        // (`cxx:949-953`).
        let mut pardeb2 = c.p[c.deb][c.pivot];
        let mut parfin2 = c.p[c.fin][c.pivot];
        if c.c1pointu {
            pardeb2 = c.p[c.deb][c.fin];
            parfin2 = c.p[c.fin][c.deb];
        }

        // `pdeb1 = CD[deb]->...->Interference(jf[deb][pivot]).PCurveOnSurf()
        // ->Value(p[deb][pivot]); pdeb2 = ...Interference(jf[deb][fin])
        // .PCurveOnSurf()->Value(pardeb2); pfin1 / pfin2 likewise on `fin`
        // (`cxx:954-978`).
        let pdeb1 = p_curve_on_surf_value(
            &c.stripes[c.deb],
            c.i[c.deb][c.pivot],
            c.jf[c.deb][c.pivot],
            c.p[c.deb][c.pivot],
        )?;
        let pdeb2 = p_curve_on_surf_value(
            &c.stripes[c.deb],
            c.i[c.deb][c.pivot],
            c.jf[c.deb][c.fin],
            pardeb2,
        )?;
        let pfin1 = p_curve_on_surf_value(
            &c.stripes[c.fin],
            c.i[c.fin][c.pivot],
            c.jf[c.fin][c.pivot],
            c.p[c.fin][c.pivot],
        )?;
        let pfin2 = p_curve_on_surf_value(
            &c.stripes[c.fin],
            c.i[c.fin][c.pivot],
            c.jf[c.fin][c.deb],
            parfin2,
        )?;

        // `sdeb` / `sfin` (`cxx:979-983`).
        let sdeb = ds_surface(base, &c.stripes[c.deb], c.i[c.deb][c.pivot])?;
        let sfin = ds_surface(base, &c.stripes[c.fin], c.i[c.fin][c.pivot])?;

        // `Bdeb = ChFi3d_mkbound(sdeb, pdeb1, pdeb2, tolapp3d, 2.e-4); Bfin =
        // ChFi3d_mkbound(sfin, pfin1, pfin2, tolapp3d, 2.e-4);`
        // (`cxx:985-986`).
        let bdeb = mkbound_segment(&sdeb, &pdeb1, &pdeb2, base.tolapp3d, 2.0e-4, false)?;
        let bfin = mkbound_segment(&sfin, &pfin1, &pfin2, base.tolapp3d, 2.0e-4, false)?;

        // `GeomFill_ConstrainedFilling fil(11, 20); if (c1pointu) {
        // fil.Init(Bpiv, Bfin, Bdeb, true); } else { fil.Init(Bpiv, Bfin,
        // Bfac, Bdeb, true); }` (`cxx:988-995`).
        let mut fil = ConstrainedFilling::new(11, 20);
        if c.c1pointu {
            fil.init_3(Arc::clone(&bpiv), bfin, bdeb, true)?;
        } else {
            let bfac = bfac
                .clone()
                .ok_or_else(|| "PerformThreeCorner: Bfac is null".to_string())?;
            fil.init_4(bpiv, bfin, bfac, bdeb, true)?;
        }

        // `occ::handle<Geom_Surface> Surfcoin = fil.Surface();
        // Surfcoin->VReverse();` (`cxx:996-997`).
        let surfcoin = fil
            .surface()
            .ok_or_else(|| "GeomFill_ConstrainedFilling::Surface".to_string())?
            .v_reversed()
            .ok_or_else(|| "Geom_Surface::VReverse".to_string())?;

        // `fdpiv->Orientation()` (`cxx:1004`).
        let o_fdpiv = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?.orientation();

        // `done = CompleteData(coin, Surfcoin, Fac, PCurveOnFace, Surf,
        // PCurveOnPiv, fdpiv->Orientation(), false, false, false, false,
        // false);` (`cxx:998-1010`).
        base.done = base.complete_data_geom(
            coin,
            &surfcoin,
            &frame.fac,
            p_curve_on_face.as_ref(),
            &frame.surf,
            Some(&p_curve_on_piv),
            o_fdpiv,
            false,
            false,
            false,
            false,
            false,
        )?;
    }

    Ok(base.done)
}
