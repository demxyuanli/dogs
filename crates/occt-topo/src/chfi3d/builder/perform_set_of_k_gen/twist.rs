//! The twisted-section repair of `ChFi3d_Builder::PerformSetOfKGen`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3318-3648`
//! (OCCT 8.0.0): `It is attempted to reprocess the squares that bore.`
//!
//! A `SurfData` whose trace twist was detected has its two side traces replaced
//! by a `GeomFill_ConstrainedFilling` patch built on the traces of its two
//! neighbours; the neighbours' trace parameters are pulled onto the new patch.
//!
//! Translation notes:
//! - OCCT mutates `cursd` / `prevsd` / `nextsd` through the handles it read out
//!   of the sequence, so a neighbour list element is written in place. The port
//!   takes the three sections out as clones, mutates them, and writes them back
//!   at the end of the iteration; when the index arithmetic of `cxx:3329-3348`
//!   names one and the same section for two roles (`len == 1`, periodic) the
//!   neighbour write is merged into that single clone,
//! - `gp_Pnt2d p2d` is an out-parameter of `ChFi3d_IntTraces` and the port's
//!   [`ch_fi3d_int_traces`] takes its reference point as an input,
//! - `ChFi3d_mkbound` is six overloads in OCCT; here
//!   [`mkbound_pfac2d`] / [`mkbound_pfac3d`] / [`mkbound_curve`] /
//!   [`mkbound_segment`] are named after the arguments they take, and each
//!   answers the p-curve it built because OCCT writes it back through `curv`,
//! - `done = CompleteData(...)` (`cxx:3606`, `:3622`) is not read afterwards, so
//!   the two calls go to [`ChFi3dBuilder::complete_data_geom`] and their `bool`
//!   is dropped.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_geom::geom_fill::ConstrainedFilling;
use occt_geom::surface::Surface as GeomSurface;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_int_traces, ch_fi3d_reparam_pcurv, mkbound_curve, mkbound_pfac2d, mkbound_pfac3d,
    mkbound_segment,
};
use crate::chfi_ds::{SpineHandle, SurfData};
use crate::shape::Face;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;

/// `It is attempted to reprocess the squares that bore.`
/// (`cxx:3318-3648`).
pub(super) fn reprocess_twists(
    builder: &mut ChFi3dBuilder,
    spine: &SpineHandle,
    seq: &mut Vec<Arc<SurfData>>,
) -> Result<(), String> {
    let len = seq.len() as i32;
    if len == 0 {
        return Ok(());
    }
    let periodic = spine
        .read()
        .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
        .is_periodic();
    let tolapp3d = builder.tolapp3d;

    for i in 1..=len {
        let ci = (i - 1) as usize;
        let mut cur = (*seq[ci]).clone();
        // `bool tw1 = cursd->TwistOnS1(); bool tw2 = cursd->TwistOnS2();`
        // (`cxx:3325-3326`).
        let tw1 = cur.twist_on_s1();
        let tw2 = cur.twist_on_s2();

        // `int iprev = i - 1; if (iprev == 0) { if (periodic) iprev = len; }
        // int inext = i + 1; if (inext > len) { if (periodic) inext = 1; else
        // inext = 0; }` (`cxx:3329-3348`).
        let mut iprev = i - 1;
        if iprev == 0 && periodic {
            iprev = len;
        }
        let mut inext = i + 1;
        if inext > len {
            inext = if periodic { 1 } else { 0 };
        }

        // `if (!tw1 && !tw2) { continue; }` (`cxx:3352-3355`).
        if !tw1 && !tw2 {
            continue;
        }

        // `int cursurf1 = cursd->IndexOfS1(); int cursurf2 =
        // cursd->IndexOfS2();` (`cxx:3361-3362`).
        let cursurf1 = cur.index_of_s1();
        let cursurf2 = cur.index_of_s2();
        // `const gp_Pnt& pd1 = cpd1.Point(); ...` (`cxx:3366-3371`).
        let pd1 = *cur.vertex_first_on_s1().point();
        let pd2 = *cur.vertex_first_on_s2().point();
        let pf1 = *cur.vertex_last_on_s1().point();
        let pf2 = *cur.vertex_last_on_s2().point();
        let ddeb = pd1.distance(&pd2);
        let dfin = pf1.distance(&pf2);
        let don1 = pd1.distance(&pf1);
        let don2 = pd2.distance(&pf2);
        let possibleon1 = don1 < 2.0 * (ddeb + dfin);
        let possibleon2 = don2 < 2.0 * (ddeb + dfin);
        // `if ((tw1 && !possibleon1) || (tw2 && !possibleon2)) {
        // Spine->SetErrorStatus(ChFiDS_TwistedSurface); throw ... }`
        // (`cxx:3376-3381`).
        if (tw1 && !possibleon1) || (tw2 && !possibleon2) {
            set_twisted_status(spine)?;
            return Err("PerformSetOfKGen: adjustment by reprocessing the non-written points".into());
        }

        // `It is checked if there are presentable neighbors`
        // (`cxx:3384-3408`).
        let mut prev = if iprev != 0 {
            Some((*seq[(iprev - 1) as usize]).clone())
        } else {
            None
        };
        let mut next = if inext != 0 {
            Some((*seq[(inext - 1) as usize]).clone())
        } else {
            None
        };

        let mut yaprevon1 = false;
        let mut yaprevon2 = false;
        let mut samesurfon1 = false;
        let mut samesurfon2 = false;
        if let Some(p) = prev.as_ref() {
            yaprevon1 = !p.twist_on_s1();
            samesurfon1 = p.index_of_s1() == cursurf1;
            yaprevon2 = !p.twist_on_s2();
            samesurfon2 = p.index_of_s2() == cursurf2;
        }
        let mut yanexton1 = false;
        let mut yanexton2 = false;
        if let Some(n) = next.as_ref() {
            yanexton1 = !n.twist_on_s1();
            if samesurfon1 {
                samesurfon1 = n.index_of_s1() == cursurf1;
            }
            yanexton2 = !n.twist_on_s2();
            if samesurfon2 {
                samesurfon2 = n.index_of_s2() == cursurf2;
            }
        }

        // `PC1 = intf1.PCurveOnFace(); PC2 = intf2.PCurveOnFace();`
        // (`cxx:3411-3412`).
        let mut pc1 = cur.interference_on_s1().p_curve_on_face().cloned();
        let mut pc2 = cur.interference_on_s2().p_curve_on_face().cloned();

        // `S1->Initialize(F1); S2->Initialize(F2);` (`cxx:3413-3418`).
        let (f1, f2) = {
            let ds = builder.my_ds.change_ds();
            let f1 = ds
                .shape(cursurf1, false)
                .ok_or_else(|| "PerformSetOfKGen: null face of IndexOfS1".to_string())?;
            let f2 = ds
                .shape(cursurf2, false)
                .ok_or_else(|| "PerformSetOfKGen: null face of IndexOfS2".to_string())?;
            (f1, f2)
        };
        let s1 = SurfaceAdaptor::of_face(&Face(f1.clone()));
        let s2 = SurfaceAdaptor::of_face(&Face(f2.clone()));
        let s1_surface = s1
            .surface()
            .cloned()
            .ok_or_else(|| "PerformSetOfKGen: null Adaptor3d_Surface S1".to_string())?;
        let s2_surface = s2
            .surface()
            .cloned()
            .ok_or_else(|| "PerformSetOfKGen: null Adaptor3d_Surface S2".to_string())?;

        let mut bon1: Option<Arc<dyn occt_geom::geom_fill::boundary::Boundary>> = None;
        let mut bon2: Option<Arc<dyn occt_geom::geom_fill::boundary::Boundary>> = None;
        let mut pointuon1 = false;
        let mut pointuon2 = false;

        // `if (tw1) { if (!yaprevon1 || !yanexton1) { ... throw } ... }`
        // (`cxx:3420-3491`).
        if tw1 {
            if !yaprevon1 || !yanexton1 {
                set_twisted_status(spine)?;
                return Err(
                    "PerformSetOfKGen: adjustment by reprocessing the non-written points: no neighbor"
                        .into(),
                );
            }
            let prevsd = prev.as_mut().unwrap();
            let nextsd = next.as_mut().unwrap();
            // `double prevpar1 = previntf1.LastParameter(); double nextpar1 =
            // nextintf1.FirstParameter();` (`cxx:3429-3430`).
            let prevpar1 = prevsd.interference_on_s1().last_parameter();
            let nextpar1 = nextsd.interference_on_s1().first_parameter();
            if samesurfon1 {
                // `pcprev1->D1(prevpar1, pdeb1, vdeb1); pcnext1->D1(nextpar1,
                // pfin1, vfin1);` (`cxx:3463-3466`).
                let mut previntf1 = prevsd.interference_on_s1().clone();
                let mut nextintf1 = nextsd.interference_on_s1().clone();
                let hit = {
                    let zero = GpPnt2d::new(0.0, 0.0);
                    ch_fi3d_int_traces(
                        prevsd, prevpar1, 1, 1, nextsd, nextpar1, 1, -1, &zero, false, true,
                    )
                };
                if let Some((nprevpar1, nnextpar1)) = hit {
                    // `previntf1.SetLastParameter(nprevpar1);
                    // nextintf1.SetFirstParameter(nnextpar1); pointuon1 = true;
                    // PC1.Nullify();` (`cxx:3455-3461`).
                    previntf1.set_last_parameter(nprevpar1);
                    nextintf1.set_first_parameter(nnextpar1);
                    prevsd.set_interference_on_s1(previntf1);
                    nextsd.set_interference_on_s1(nextintf1);
                    pointuon1 = true;
                    pc1 = None;
                } else {
                    let pcprev1 = previntf1
                        .p_curve_on_face()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null pcurve of the previous trace".to_string())?;
                    let pcnext1 = nextintf1
                        .p_curve_on_face()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null pcurve of the next trace".to_string())?;
                    let (pdeb1, vdeb1) = d1_2d(&pcprev1, prevpar1)?;
                    let (pfin1, vfin1) = d1_2d(&pcnext1, nextpar1)?;
                    let (pc, b) = mkbound_pfac2d(
                        &s1, -1, &pdeb1, &vdeb1, 1, &pfin1, &vfin1, tolapp3d, 2.0e-4,
                    )?;
                    pc1 = Some(pc);
                    bon1 = Some(b);
                }
            } else {
                // `here the base is on 3D tangents of neighbors.`
                // (`cxx:3469-3487`).
                // `gp_Pnt Pdeb1, Pfin1; ... c3dprev1->D1(prevpar1, Pdeb1,
                // Vdeb1);` (`cxx:3473-3477`): OCCT computes the two 3D points
                // but only the vectors reach `mkbound`.
                let (_pdeb1, vdeb1) = {
                    let ds = builder.my_ds.change_ds();
                    let c3dprev1 = ds
                        .curve(prevsd.interference_on_s1().line_index())
                        .curve()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null 3D curve of the previous trace".to_string())?;
                    c3dprev1.d1(prevpar1)
                };
                let (_pfin1, vfin1) = {
                    let ds = builder.my_ds.change_ds();
                    let c3dnext1 = ds
                        .curve(nextsd.interference_on_s1().line_index())
                        .curve()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null 3D curve of the next trace".to_string())?;
                    c3dnext1.d1(nextpar1)
                };
                // `pdeb1 = PC1->Value(pardeb1); pfin1 = PC1->Value(parfin1);`
                // (`cxx:3481-3484`).
                let pc1_curve = pc1
                    .as_ref()
                    .ok_or_else(|| "PerformSetOfKGen: null PC1".to_string())?;
                let pardeb1 = cur.interference_on_s1().first_parameter();
                let parfin1 = cur.interference_on_s1().last_parameter();
                let pdeb1 = pc1_curve.d0(pardeb1);
                let pfin1 = pc1_curve.d0(parfin1);
                let (pc, b) =
                    mkbound_pfac3d(&s1, -1, &pdeb1, &vdeb1, 1, &pfin1, &vfin1, tolapp3d, 2.0e-4)?;
                pc1 = Some(pc);
                bon1 = Some(b);
            }
        } else {
            // `Bon1 = ChFi3d_mkbound(S1, PC1, tolapp3d, 2.e-4);` (`cxx:3490`).
            let curve = pc1
                .as_ref()
                .ok_or_else(|| "PerformSetOfKGen: null PC1".to_string())?;
            bon1 = Some(mkbound_curve(&s1_surface, curve, tolapp3d, 2.0e-4, false));
        }

        // `if (tw2) { ... } else { Bon2 = ChFi3d_mkbound(S2, PC2, tolapp3d,
        // 2.e-4); }` (`cxx:3492-3561`).
        if tw2 {
            if !yaprevon2 || !yanexton2 {
                return Err(
                    "PerformSetOfKGen: adjustment by reprocessing the non-written points: no neighbor"
                        .into(),
                );
            }
            let prevsd = prev.as_mut().unwrap();
            let nextsd = next.as_mut().unwrap();
            let prevpar2 = prevsd.interference_on_s2().last_parameter();
            let nextpar2 = nextsd.interference_on_s2().first_parameter();
            if samesurfon2 {
                let mut previntf2 = prevsd.interference_on_s2().clone();
                let mut nextintf2 = nextsd.interference_on_s2().clone();
                let hit = {
                    let zero = GpPnt2d::new(0.0, 0.0);
                    ch_fi3d_int_traces(
                        prevsd, prevpar2, 2, 1, nextsd, nextpar2, 2, -1, &zero, false, true,
                    )
                };
                if let Some((nprevpar2, nnextpar2)) = hit {
                    previntf2.set_last_parameter(nprevpar2);
                    nextintf2.set_first_parameter(nnextpar2);
                    prevsd.set_interference_on_s2(previntf2);
                    nextsd.set_interference_on_s2(nextintf2);
                    pointuon2 = true;
                    pc2 = None;
                } else {
                    let pcprev2 = previntf2
                        .p_curve_on_face()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null pcurve of the previous trace".to_string())?;
                    let pcnext2 = nextintf2
                        .p_curve_on_face()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null pcurve of the next trace".to_string())?;
                    let (pdeb2, vdeb2) = d1_2d(&pcprev2, prevpar2)?;
                    let (pfin2, vfin2) = d1_2d(&pcnext2, nextpar2)?;
                    let (pc, b) = mkbound_pfac2d(
                        &s2, -1, &pdeb2, &vdeb2, 1, &pfin2, &vfin2, tolapp3d, 2.0e-4,
                    )?;
                    pc2 = Some(pc);
                    bon2 = Some(b);
                }
            } else {
                let (_pdeb2, vdeb2) = {
                    let ds = builder.my_ds.change_ds();
                    let c3dprev2 = ds
                        .curve(prevsd.interference_on_s2().line_index())
                        .curve()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null 3D curve of the previous trace".to_string())?;
                    c3dprev2.d1(prevpar2)
                };
                let (_pfin2, vfin2) = {
                    let ds = builder.my_ds.change_ds();
                    let c3dnext2 = ds
                        .curve(nextsd.interference_on_s2().line_index())
                        .curve()
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: null 3D curve of the next trace".to_string())?;
                    c3dnext2.d1(nextpar2)
                };
                let pc2_curve = pc2
                    .as_ref()
                    .ok_or_else(|| "PerformSetOfKGen: null PC2".to_string())?;
                let pardeb2 = cur.interference_on_s2().first_parameter();
                let parfin2 = cur.interference_on_s2().last_parameter();
                let pdeb2 = pc2_curve.d0(pardeb2);
                let pfin2 = pc2_curve.d0(parfin2);
                let (pc, b) =
                    mkbound_pfac3d(&s2, -1, &pdeb2, &vdeb2, 1, &pfin2, &vfin2, tolapp3d, 2.0e-4)?;
                pc2 = Some(pc);
                bon2 = Some(b);
            }
        } else {
            let curve = pc2
                .as_ref()
                .ok_or_else(|| "PerformSetOfKGen: null PC2".to_string())?;
            bon2 = Some(mkbound_curve(&s2_surface, curve, tolapp3d, 2.0e-4, false));
        }

        // `The parameters of neighbor traces are updated, so straight lines uv
        // are pulled.` (`cxx:3562-3585`).
        let prevsd = prev.as_mut().unwrap();
        let nextsd = next.as_mut().unwrap();
        let (sprev, snext) = {
            let ds = builder.my_ds.change_ds();
            let sprev = ds
                .surface(prevsd.surf())
                .surface()
                .cloned()
                .ok_or_else(|| "PerformSetOfKGen: null surface of the previous section".to_string())?;
            let snext = ds
                .surface(nextsd.surf())
                .surface()
                .cloned()
                .ok_or_else(|| "PerformSetOfKGen: null surface of the next section".to_string())?;
            (sprev, snext)
        };
        let (prevpar1, nextpar1) = (
            prevsd.interference_on_s1().last_parameter(),
            nextsd.interference_on_s1().first_parameter(),
        );
        let (prevpar2, nextpar2) = (
            prevsd.interference_on_s2().last_parameter(),
            nextsd.interference_on_s2().first_parameter(),
        );
        let pdebs1 = d0_2d(
            prevsd
                .interference_on_s1()
                .p_curve_on_surf()
                .ok_or_else(|| "PerformSetOfKGen: null pcurve on surf".to_string())?,
            prevpar1,
        )?;
        let pdebs2 = d0_2d(
            prevsd
                .interference_on_s2()
                .p_curve_on_surf()
                .ok_or_else(|| "PerformSetOfKGen: null pcurve on surf".to_string())?,
            prevpar2,
        )?;
        let pfins1 = d0_2d(
            nextsd
                .interference_on_s1()
                .p_curve_on_surf()
                .ok_or_else(|| "PerformSetOfKGen: null pcurve on surf".to_string())?,
            nextpar1,
        )?;
        let pfins2 = d0_2d(
            nextsd
                .interference_on_s2()
                .p_curve_on_surf()
                .ok_or_else(|| "PerformSetOfKGen: null pcurve on surf".to_string())?,
            nextpar2,
        )?;
        // `Bdeb = ChFi3d_mkbound(sprev, pdebs1, pdebs2, tolapp3d, 2.e-4);`
        // `Bfin = ChFi3d_mkbound(snext, pfins1, pfins2, tolapp3d, 2.e-4);`
        // (`cxx:3583-3584`).
        let bdeb = mkbound_segment(&sprev, &pdebs1, &pdebs2, tolapp3d, 2.0e-4, false)?;
        let bfin = mkbound_segment(&snext, &pfins1, &pfins2, tolapp3d, 2.0e-4, false)?;

        // `GeomFill_ConstrainedFilling fil(11, 20); ... fil.Init(...)`
        // (`cxx:3586-3598`).
        let mut fil = ConstrainedFilling::new(11, 20);
        let bon_deb = bdeb;
        let bon_fin = bfin;
        if pointuon1 {
            fil.init_3(bon2.clone().unwrap(), bon_fin, bon_deb, true)?;
        } else if pointuon2 {
            fil.init_3(bon1.clone().unwrap(), bon_fin, bon_deb, true)?;
        } else {
            fil.init_4(
                bon1.clone().unwrap(),
                bon_fin,
                bon2.clone().unwrap(),
                bon_deb,
                true,
            )?;
        }

        // `ChFi3d_ReparamPcurv(0., 1., PC1); ChFi3d_ReparamPcurv(0., 1.,
        // PC2);` (`cxx:3600-3601`).
        ch_fi3d_reparam_pcurv(0.0, 1.0, &mut pc1)?;
        ch_fi3d_reparam_pcurv(0.0, 1.0, &mut pc2)?;

        // `Handle(Geom_Surface) newsurf = fil.Surface();` (`cxx:3602`).
        let mut new_surface = fil
            .surface()
            .cloned()
            .ok_or_else(|| "PerformSetOfKGen: GeomFill_ConstrainedFilling::Surface()".to_string())?;
        let new_surface: Arc<dyn GeomSurface> = if pointuon1 {
            // `newsurf->VReverse(); // we return to direction 1 from 2;`
            // (`cxx:3606`).
            new_surface.v_reverse();
            Arc::new(new_surface)
        } else {
            Arc::new(new_surface)
        };

        // `done = CompleteData(cursd, newsurf, S1, PC1, S2, PC2, ...)`
        // (`cxx:3607-3636`).
        let orientation = if pointuon1 {
            f2.orientation()
        } else {
            f1.orientation()
        };
        builder.complete_data_geom(
            &mut cur,
            &new_surface,
            &s1,
            pc1.as_ref(),
            &s2,
            pc2.as_ref(),
            orientation,
            !pointuon1,
            false,
            false,
            false,
            false,
        )?;
        if pointuon1 {
            // `cursd->ChangeIndexOfS1(0);` (`cxx:3618`).
            cur.set_index_of_s1(0);
        } else if pointuon2 {
            // `cursd->ChangeIndexOfS2(0);` (`cxx:3634`).
            cur.set_index_of_s2(0);
        }

        // `if (tw1) { prevsd->ChangeVertexLastOnS1().SetPoint(cpd1.Point());
        // nextsd->ChangeVertexFirstOnS1().SetPoint(cpf1.Point()); }`
        // (`cxx:3639-3643`).
        if tw1 {
            let mut cp = prevsd.vertex_last_on_s1().clone();
            cp.set_point(&pd1);
            prevsd.set_vertex_last_on_s1(cp);
            let mut cp = nextsd.vertex_first_on_s1().clone();
            cp.set_point(&pf1);
            nextsd.set_vertex_first_on_s1(cp);
        }
        // `if (tw2) { ... }` (`cxx:3644-3647`).
        if tw2 {
            let mut cp = prevsd.vertex_last_on_s2().clone();
            cp.set_point(&pd2);
            prevsd.set_vertex_last_on_s2(cp);
            let mut cp = nextsd.vertex_first_on_s2().clone();
            cp.set_point(&pf2);
            nextsd.set_vertex_first_on_s2(cp);
        }

        // The three handles OCCT mutated in place are written back. When the
        // index arithmetic names one section for several roles, the neighbour
        // writes are merged into that single clone.
        let prev_index = if iprev != 0 { Some((iprev - 1) as usize) } else { None };
        let next_index = if inext != 0 { Some((inext - 1) as usize) } else { None };
        if prev_index == next_index {
            if let (Some(idx), Some(p), Some(n)) = (prev_index, prev.as_ref(), next.as_ref()) {
                // Same object: the interference and vertex writes of both roles
                // land on one section, exactly as the two handles alias.
                let mut merged = p.clone();
                merged.set_interference_on_s1(n.interference_on_s1().clone());
                merged.set_interference_on_s2(n.interference_on_s2().clone());
                merged.set_vertex_first_on_s1(n.vertex_first_on_s1().clone());
                merged.set_vertex_first_on_s2(n.vertex_first_on_s2().clone());
                if idx == ci {
                    // The current section is that same object: keep the
                    // `CompleteData` products of `cur` and the neighbour writes.
                    let mut merged = cur.clone();
                    merged.set_interference_on_s1(p.interference_on_s1().clone());
                    merged.set_interference_on_s2(p.interference_on_s2().clone());
                    merged.set_vertex_first_on_s1(n.vertex_first_on_s1().clone());
                    merged.set_vertex_first_on_s2(n.vertex_first_on_s2().clone());
                    merged.set_vertex_last_on_s1(p.vertex_last_on_s1().clone());
                    merged.set_vertex_last_on_s2(p.vertex_last_on_s2().clone());
                    seq[ci] = Arc::new(merged);
                } else {
                    seq[idx] = Arc::new(merged);
                    seq[ci] = Arc::new(cur);
                }
            }
        } else {
            if let (Some(idx), Some(p)) = (prev_index, prev.as_ref()) {
                seq[idx] = Arc::new(p.clone());
            }
            if let Some(idx) = next_index {
                seq[idx] = Arc::new(next.as_ref().unwrap().clone());
            }
            seq[ci] = Arc::new(cur);
        }
    }
    Ok(())
}

/// `Spine->SetErrorStatus(ChFiDS_TwistedSurface)` (`cxx:3379`, `:3423`).
fn set_twisted_status(spine: &SpineHandle) -> Result<(), String> {
    spine
        .write()
        .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
        .set_error_status(crate::chfi_ds::ErrorStatus::TwistedSurface);
    Ok(())
}

/// `Geom2d_Curve::D1(par)` as `(gp_Pnt2d, gp_Vec2d)`.
fn d1_2d(
    c: &Arc<dyn occt_geom2d::curve::Curve2d>,
    par: f64,
) -> Result<(GpPnt2d, GpVec2d), String> {
    Ok(c.d1(par))
}

/// `Geom2d_Curve::Value(par)`.
fn d0_2d(c: &Arc<dyn occt_geom2d::curve::Curve2d>, par: f64) -> Result<GpPnt2d, String> {
    Ok(c.d0(par))
}
