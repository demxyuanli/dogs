//! The two closing blocks of `ChFi3d_Builder::PerformSetOfKGen`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx` (OCCT
//! 8.0.0):
//!
//! - `The tolerance of points is updated.` (`cxx:3650-3754`), the
//!   [`update_point_tolerances`] loop: the neighbouring sections share the
//!   common point that is on an arc, the transition is moved when the two
//!   sections are on the same face, and the box of both end points becomes the
//!   new tolerance,
//! - `The connections edge/new faces are updated.` (`cxx:3755-3878`), the
//!   [`update_edge_face_links`] loop: each `ElSpine` is cut into the edges of
//!   the spine it covers and every `SurfData` of that range is appended to
//!   `myEVIMap(edge)`.
//!
//! Translation notes:
//! - `ChFiDS_CommonPoint&` / `ChFiDS_FaceInterference&` are values here, so
//!   each section is taken out of the sequence as a clone, repaired, and put
//!   back; when OCCT's `j = i % len + 1` names the same section (`len == 1`)
//!   the two writes land on the same clone, in OCCT order,
//! - `myEFMap(Arc)` is [`crate::chfi_ds::ChFiDSMap::find_from_key`], the
//!   face list of an edge,
//! - `Spine->Index` and `Spine->Value` answer `Option` where OCCT reads an
//!   index unconditionally; the missing index is the `Err`,
//! - `Extrema_LocateExtPC` is [`occt_geom::extrema_pc::GenLocateExtPC`],
//!   whose `with_params` is the `(P, C, U0, TolU)` constructor.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::extrema_pc::GenLocateExtPC;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_enlarge_box_edge, ch_fi3d_enlarge_box_surf_data,
};
use crate::chfi3d::builder_2_helpers::{ch_fi3d_box_diag, change_transition};
use crate::chfi_ds::{SpineHandle, SurfData};

use super::super::ChFi3dBuilder;

/// `The tolerance of points is updated.` (`cxx:3650-3754`).
pub(super) fn update_point_tolerances(
    builder: &mut ChFi3dBuilder,
    seq: &mut [Arc<SurfData>],
    len: i32,
    last: i32,
) -> Result<(), String> {
    for i in 1..last {
        // `int j = i % len + 1;` (`cxx:3654`).
        let j = i % len + 1;
        let ci = (i - 1) as usize;
        let cj = (j - 1) as usize;

        // `ChFiDS_CommonPoint& curp1 = cursd->ChangeVertexLastOnS1(); ...`
        // (`cxx:3658-3659`).
        let mut cur = (*seq[ci]).clone();
        let mut nxt = (*seq[cj]).clone();

        let mut curp1 = cur.vertex_last_on_s1().clone();
        let mut nextp1 = nxt.vertex_first_on_s1().clone();
        let mut curp2 = cur.vertex_last_on_s2().clone();
        let mut nextp2 = nxt.vertex_first_on_s2().clone();

        // `if (cursd->IsOnCurve1()) curs1 = cursd->IndexOfC1(); else
        // curs1 = cursd->IndexOfS1();` (`cxx:3660-3677`).
        let curs1 = surface_or_curve_index1(&cur);
        let curs2 = surface_or_curve_index2(&cur);
        // `double tol1 = std::max(curp1.Tolerance(), nextp1.Tolerance());`
        let mut tol1 = curp1.tolerance().max(nextp1.tolerance());
        let mut tol2 = curp2.tolerance().max(nextp2.tolerance());
        let nexts1 = surface_or_curve_index1(&nxt);
        let nexts2 = surface_or_curve_index2(&nxt);

        // `if (!curp1.IsOnArc() && nextp1.IsOnArc()) { ... }`
        // (`cxx:3696-3730`).
        if !curp1.is_on_arc() && nextp1.is_on_arc() {
            curp1 = nextp1.clone();
            if curs1 == nexts1 && !nxt.is_on_curve1() {
                change_transition(&nextp1, &mut curp1, nexts1, builder.my_ds.change_ds())?;
            }
        } else if curp1.is_on_arc() && !nextp1.is_on_arc() {
            nextp1 = curp1.clone();
            if curs1 == nexts1 && !cur.is_on_curve1() {
                change_transition(&curp1, &mut nextp1, curs1, builder.my_ds.change_ds())?;
            }
        }
        if !curp2.is_on_arc() && nextp2.is_on_arc() {
            curp2 = nextp2.clone();
            if curs2 == nexts2 && !nxt.is_on_curve2() {
                change_transition(&nextp2, &mut curp2, nexts2, builder.my_ds.change_ds())?;
            }
        } else if curp2.is_on_arc() && !nextp2.is_on_arc() {
            nextp2 = curp2.clone();
            if curs2 == nexts2 && !cur.is_on_curve2() {
                change_transition(&curp2, &mut nextp2, curs2, builder.my_ds.change_ds())?;
            }
        }

        // `curp1.SetTolerance(tol1); nextp1.SetTolerance(tol1);
        // curp2.SetTolerance(tol2); nextp2.SetTolerance(tol2);`
        // (`cxx:3732-3735`).
        curp1.set_tolerance(tol1);
        nextp1.set_tolerance(tol1);
        curp2.set_tolerance(tol2);
        nextp2.set_tolerance(tol2);

        // `Bnd_Box b1, b2; if (curp1.IsOnArc()) ChFi3d_EnlargeBox(
        // curp1.Arc(), myEFMap(curp1.Arc()), curp1.ParameterOnArc(), b1);`
        // (`cxx:3737-3746`).
        let mut b1 = occt_core::bnd::BndBox::new();
        let mut b2 = occt_core::bnd::BndBox::new();
        if let (Some(arc), Some(param)) = (curp1.arc(), curp1.parameter_on_arc()) {
            let faces = builder.my_ef_map.find_from_key(arc);
            ch_fi3d_enlarge_box_edge(arc, faces, param, &mut b1);
        }
        if let (Some(arc), Some(param)) = (curp2.arc(), curp2.parameter_on_arc()) {
            let faces = builder.my_ef_map.find_from_key(arc);
            ch_fi3d_enlarge_box_edge(arc, faces, param, &mut b2);
        }

        // `Handle(ChFiDS_Stripe) bidst; ChFi3d_EnlargeBox(DStr, bidst, cursd,
        // b1, b2, false); ChFi3d_EnlargeBox(DStr, bidst, nextsd, b1, b2, true);`
        // (`cxx:3747-3748`).
        ch_fi3d_enlarge_box_surf_data(
            builder.my_ds.change_ds(),
            None,
            &cur,
            &mut b1,
            &mut b2,
            false,
        )?;
        ch_fi3d_enlarge_box_surf_data(
            builder.my_ds.change_ds(),
            None,
            &nxt,
            &mut b1,
            &mut b2,
            true,
        )?;

        // `tol1 = ChFi3d_BoxDiag(b1); tol2 = ChFi3d_BoxDiag(b2);
        // curp1.SetTolerance(tol1); ...` (`cxx:3749-3754`).
        tol1 = ch_fi3d_box_diag(&b1);
        tol2 = ch_fi3d_box_diag(&b2);
        curp1.set_tolerance(tol1);
        nextp1.set_tolerance(tol1);
        curp2.set_tolerance(tol2);
        nextp2.set_tolerance(tol2);

        // The four writes back into the sections (`cxx:3658-3659`'s
        // `Change*` references).
        cur.set_vertex_last_on_s1(curp1);
        cur.set_vertex_last_on_s2(curp2);
        nxt.set_vertex_first_on_s1(nextp1);
        nxt.set_vertex_first_on_s2(nextp2);
        if cj == ci {
            // OCCT wrote through one and the same `SurfData`; replay the two
            // first-end writes on the single clone in the same order.
            let (f1, f2) = (
                nxt.vertex_first_on_s1().clone(),
                nxt.vertex_first_on_s2().clone(),
            );
            cur.set_vertex_first_on_s1(f1);
            cur.set_vertex_first_on_s2(f2);
            seq[ci] = Arc::new(cur);
        } else {
            seq[ci] = Arc::new(cur);
            seq[cj] = Arc::new(nxt);
        }
    }
    Ok(())
}

/// `The connections edge/new faces are updated.` (`cxx:3755-3878`).
pub(super) fn update_edge_face_links(
    builder: &mut ChFi3dBuilder,
    spine: &SpineHandle,
    seq: &[Arc<SurfData>],
    len: i32,
    periodic: bool,
) -> Result<(), String> {
    let tolesp = builder.tolesp;
    let guides: Vec<Arc<crate::chfi_ds::ElSpine>> = spine
        .read()
        .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
        .el_spines()
        .to_vec();

    for curhels in &guides {
        // `double WF = curhels->FirstParameter(); double WL =
        // curhels->LastParameter(); int IF, IL; double nwf = WF, nwl = WL;
        // double period = 0.; int nbed = Spine->NbEdges();`
        // (`cxx:3760-3766`).
        let wf = curhels.first_parameter();
        let wl = curhels.last_parameter();
        let (if_, il_, _period, nbed) = {
            let guard = spine
                .write()
                .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?;
            let nbed = guard.nb_edges() as i32;
            let (nwf, nwl, period) = if periodic {
                let period = guard
                    .period()
                    .ok_or_else(|| "PerformSetOfKGen: Spine::Period()".to_string())?;
                (
                    clib::in_period(wf, -tolesp, period - tolesp),
                    clib::in_period(wl, tolesp, period + tolesp),
                    period,
                )
            } else {
                (wf, wl, 0.0)
            };
            // `IF = Spine->Index(nwf, true); ... IL = Spine->Index(nwl,
            // false); if (nwl < nwf + tolesp) IL += nbed;` (`cxx:3767-3779`).
            let if_ = guard
                .index(nwf, true)
                .ok_or_else(|| "PerformSetOfKGen: Spine::Index(WF)".to_string())?;
            let mut il_ = guard
                .index(nwl, false)
                .ok_or_else(|| "PerformSetOfKGen: Spine::Index(WL)".to_string())?;
            if periodic && nwl < nwf + tolesp {
                il_ += nbed;
            }
            (if_, il_, period, nbed)
        };

        if if_ == il_ {
            // `int IFloc = IF; if (periodic) IFloc = (IF - 1) % nbed + 1;`
            // (`cxx:3785-3790`).
            let mut ifloc = if_;
            if periodic {
                ifloc = (if_ - 1).rem_euclid(nbed) + 1;
            }
            let ej = spine
                .read()
                .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
                .edges(ifloc as usize)
                .cloned()
                .ok_or_else(|| "PerformSetOfKGen: Spine::Edges(IFloc)".to_string())?;
            for sd in seq.iter().take(len as usize) {
                // `double fp = cursd->FirstSpineParam(); double lp =
                // cursd->LastSpineParam(); if (lp < WF + tolesp || fp > WL -
                // tolesp) continue;` (`cxx:3793-3800`).
                let fp = sd.first_spine_param();
                let lp = sd.last_spine_param();
                if lp < wf + tolesp || fp > wl - tolesp {
                    continue;
                }
                builder.change_evi(&ej).push(sd.surf());
            }
        } else if if_ < il_ {
            // `NCollection_Array1<double> wv(IF, IL - 1); for (i = IF; i < IL;
            // i++) { ... }` (`cxx:3811-3832`).
            let mut wv: Vec<f64> = Vec::with_capacity((il_ - if_) as usize);
            for i in if_..il_ {
                let mut iloc = i;
                if periodic {
                    iloc = (i - 1).rem_euclid(nbed) + 1;
                }
                let (mut wi, pv) = {
                    let mut guard = spine
                        .write()
                        .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?;
                    let mut wi = guard
                        .last_parameter_of(iloc)
                        .ok_or_else(|| "PerformSetOfKGen: Spine::LastParameter".to_string())?;
                    if periodic {
                        wi = clib::in_period(wi, wf, wf + _period);
                    }
                    let pv = guard
                        .value(wi)
                        .ok_or_else(|| "PerformSetOfKGen: Spine::Value".to_string())?;
                    (wi, pv)
                };
                // `Extrema_LocateExtPC ext(pv, *curhels, wi, 1.e-8); wv(i) =
                // wi; if (ext.IsDone()) wv(i) = ext.Point().Parameter();`
                // (`cxx:3826-3831`).
                let ext = GenLocateExtPC::with_params(
                    &pv,
                    curhels.as_ref() as &dyn GeomCurve,
                    wi,
                    1.0e-8,
                );
                if ext.is_done() {
                    wi = ext.point(1).parameter();
                }
                wv.push(wi);
            }

            // `for (i = 1; i <= len; i++) { ... }` (`cxx:3833-3876`).
            for sd in seq.iter().take(len as usize) {
                let fp = sd.first_spine_param();
                let lp = sd.last_spine_param();
                if lp < wf + tolesp || fp > wl - tolesp {
                    continue;
                }
                let mut jf = 0;
                for j in if_..il_ {
                    jf = j;
                    if fp < wv[(j - if_) as usize] - tolesp {
                        break;
                    }
                }
                let mut jl = 0;
                for j in if_..il_ {
                    jl = j;
                    if lp < wv[(j - if_) as usize] + tolesp {
                        break;
                    }
                }
                for j in jf..=jl {
                    let mut jloc = j;
                    if periodic {
                        jloc = (j - 1).rem_euclid(nbed) + 1;
                    }
                    let ej = spine
                        .read()
                        .map_err(|_| "PerformSetOfKGen: poisoned spine".to_string())?
                        .edges(jloc as usize)
                        .cloned()
                        .ok_or_else(|| "PerformSetOfKGen: Spine::Edges(jloc)".to_string())?;
                    builder.change_evi(&ej).push(sd.surf());
                }
            }
        }
    }
    Ok(())
}

/// `cursd->IsOnCurve1() ? cursd->IndexOfC1() : cursd->IndexOfS1()`
/// (`cxx:3660-3665`).
fn surface_or_curve_index1(sd: &SurfData) -> i32 {
    if sd.is_on_curve1() {
        sd.index_of_c1().unwrap_or(0)
    } else {
        sd.index_of_s1()
    }
}

/// The `S2` twin (`cxx:3667-3674`).
fn surface_or_curve_index2(sd: &SurfData) -> i32 {
    if sd.is_on_curve2() {
        sd.index_of_c2().unwrap_or(0)
    } else {
        sd.index_of_s2()
    }
}
