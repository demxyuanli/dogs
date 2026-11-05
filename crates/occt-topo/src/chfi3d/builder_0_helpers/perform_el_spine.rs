//! `ChFi3d_PerformElSpine` — rebuild one `ChFiDS_ElSpine`'s 3D curve from the
//! edges of its spine.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4797-5327`
//! (OCCT 8.0.0), declaration `ChFi3d_Builder_0.hxx:539-544`.
//!
//! The body trims the first / last spine edge to the `ElSpine` range, converts
//! every piece to a B-spline, smooths the pieces that are neither `C2` nor
//! `ChFi3d_IsSmooth` with `ChFi3d_ApproxByC2`, concatenates them through
//! `GeomConvert_CompCurveToBSplineCurve(TC, Convert_QuasiAngular)`,
//! reparametrizes the knots onto `[Wrefdeb, Wreffin]`, builds the `C3` tangent
//! extensions at both extremities (`GeomLib::ExtendCurveToPoint` +
//! `GCPnts_AbscissaPoint` + the file static `GoodExt`, `cxx:4773-4791`), trims
//! to the spine domain, lowers the interior knot multiplicities and, for a
//! non-periodic `ElSpine`, shifts the end poles with
//! `GeomLib::AdjustExtremity`.
//!
//! `PerformSetOfKPart` calls it once per `ElSpine` of the spine
//! (`ChFi3d_Builder_2.cxx:3267-3270`) and once per offset `ElSpine` when
//! `Spine->Mode() == ChFiDS_ConstThroatWithPenetrationChamfer`
//! (`:3271-3278`); `PerformSetOfKGen` is the other caller.
//!
//! Divergences that are not heuristics:
//! - `GeomConvert::C0BSplineToC1BSplineCurve(BSpline, tol, 0.1)`
//!   (`cxx:5203`) comes from `occt_geom::geom_convert_concat`; the
//!   `Standard_ConstructionError` / `Standard_DimensionError` throws of that
//!   chain are the `Err` here, so a failed conversion ends `PerformElSpine`
//!   instead of leaving `BSpline` unconverted.
//! - `BRep_Tool::Curve(E, First, Last)` is [`BRepTool::edge_curve_world`]:
//!   OCCT hands back the curve in the edge's own frame, and every quantity it
//!   is compared against here (`Spine->D1`, `Spine->Parameter`, the pcurve
//!   adaptor behind `CurrentElementarySpine`) is world, so the port applies the
//!   location at the read.
//! - `BRepLib::BuildCurve3d(E)` is [`build_curve3d`] with the default arguments
//!   of `BRepLib.hxx:90-94` (`1.0e-5`, `GeomAbs_C1`, `14`, `0`).
//! - `Geom_BSplineCurve::NbKnots()` (`cxx:5117`, `:5118`, `:5223`) is
//!   `distinct_knots().len()`; the port's `GeomBSplineCurve::nb_knots()` is the
//!   flat knot vector length.
//! - `GeomLib::ExtendCurveToPoint` / `GeomLib::AdjustExtremity` take
//!   `&mut GeomBSplineCurve` (see `occt_geom::geom_lib_extend_curve`), so the
//!   `anExtCurve` handle copy of `cxx:5131` / `:5162` / `:5314` is a
//!   `bspline.clone()`.
//! - `gp_Pnt&` / `gp_Vec&` out-parameters are returned values; a null edge,
//!   curve, adaptor or spine query, and the `Standard_ConstructionError`
//!   throws of `cxx:5022` / `:5092`, are the `Err`.

use std::sync::Arc;

use occt_core::bspl::knots::reparametrize;
use occt_core::convert::ParameterisationType;
use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpDir, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::{curve_to_bspline_curve, CompCurveToBSplineCurve};
use occt_geom::curve::Curve;
use occt_geom::gcpnts::abscissa_point_with_guess;
use occt_geom::geom_convert_concat::c0_bspline_to_c1_bspline_curve;
use occt_geom::geom_lib_extend_curve::{adjust_extremity, extend_curve_to_point};
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_helpers::{approx_by_c2, is_smooth};
use crate::chfi_ds::{ElSpine, Spine, SpineHandle};
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::common_vertex;

use super::curve_cleaner::curve_cleaner;

/// `BRepLib::BuildCurve3d` default tolerance (`BRepLib.hxx:91`).
const BUILD_CURVE3D_TOL: f64 = 1.0e-5;
/// `BRepLib::BuildCurve3d` defaults for the degree and segment bounds
/// (`BRepLib.hxx:93-94`); `GeomAbs_C1` is the continuity default (`:92`).
const BUILD_CURVE3D_MAX_DEGREE: i32 = 14;
const BUILD_CURVE3D_MAX_SEGMENT: i32 = 0;

fn lock_err() -> String {
    "ChFi3d_PerformElSpine: poisoned spine".to_string()
}

fn spine_rd<T>(spine: &SpineHandle, f: impl FnOnce(&Spine) -> T) -> Result<T, String> {
    let guard = spine.read().map_err(|_| lock_err())?;
    Ok(f(&guard))
}

fn spine_wr<T>(spine: &SpineHandle, f: impl FnOnce(&mut Spine) -> T) -> Result<T, String> {
    let mut guard = spine.write().map_err(|_| lock_err())?;
    Ok(f(&mut guard))
}

/// `IsOffset ? Spine->OffsetEdges(Index) : Spine->Edges(Index)` (`cxx:4846`,
/// `:5006`), 1-based.
fn spine_edge(spine: &SpineHandle, index: i32, is_offset: bool) -> Result<TopoShape, String> {
    let i = usize::try_from(index)
        .map_err(|_| format!("ChFi3d_PerformElSpine: Spine edge index {index}"))?;
    spine_rd(spine, |s| {
        if is_offset {
            s.offset_edges(i).cloned()
        } else {
            s.edges(i).cloned()
        }
    })?
    .ok_or_else(|| format!("ChFi3d_PerformElSpine: Spine->Edges({index}) out of range"))
}

/// `Spine->Parameter(IE, W, U, Oriented = false)` (`ChFiDS_Spine.cxx:582-612`):
/// the absolute abscissa `w` on edge `index` as a curve parameter.
fn spine_parameter(spine: &SpineHandle, index: i32, w: f64) -> Result<f64, String> {
    spine_wr(spine, |s| s.parameter_of(index, w, false))?
        .ok_or_else(|| format!("ChFi3d_PerformElSpine: Spine->Parameter({index}, {w})"))
}

/// `gp_Ax1(PntForElSpine, DirForElSpine)` from `Cv->D1(ParForElSpine, ...)`
/// (`cxx:4858-4870`, `:5024-5031`). OCCT's `gp_Dir(gp_Vec)` raises on a
/// vanishing vector (`gp_Dir.cxx`), which is the `Err`.
fn axis_at(cv: &dyn Curve, par: f64) -> Result<GpAx1, String> {
    let (p, d) = cv.d1(par);
    let dir = GpDir::from_vec(&d)
        .map_err(|_| "ChFi3d_PerformElSpine: gp_Dir of a zero derivative".to_string())?;
    Ok(GpAx1::new(p, dir))
}

/// `GoodExt(C, V, f, l, a)` (`cxx:4773-4791`): the extension is accepted when
/// the angle between the tangent and `V` stays below `a * t + 0.002` at the six
/// samples `t = i * 0.2`, `i = 0..5`.
fn good_ext(c: &dyn Curve, v: &GpVec, f: f64, l: f64, a: f64) -> bool {
    for i in 0..6 {
        let t = i as f64 * 0.2;
        let (_, d1) = c.d1((1.0 - t) * f + t * l);
        let ang = d1.angle(v);
        let angref = a * t + 0.002;
        if ang > angref {
            return false;
        }
    }
    true
}

/// `bs->Multiplicity(ii)`, `ii` 1-based (`Geom_BSplineCurve_1.cxx:365-375`).
///
/// An index the curve no longer has answers `0`: OCCT would read its
/// multiplicity array out of range, which the `RemoveKnot` arms of
/// `cxx:5258-5269` make unreachable for the range that loop walks.
fn knot_multiplicity(bs: &GeomBSplineCurve, ii: i32) -> i32 {
    let index = match usize::try_from(ii - 1) {
        Ok(i) => i,
        Err(_) => return 0,
    };
    bs.multiplicities().get(index).copied().unwrap_or(0)
}

/// The common tail of both edge blocks: `TC = new Geom_TrimmedCurve(Cv, First,
/// Last); BS = GeomConvert::CurveToBSplineCurve(TC); CurveCleaner(BS,
/// |WL-WF|*1.e-4, 0);` followed by the `GeomAbs_C2` / `ChFi3d_IsSmooth` pair
/// (`cxx:4977-4993`, `:5064-5080`).
///
/// Answers the curve to concatenate (`TC`, which is the approximation when the
/// smoothing arm fires) and whether `ChFi3d_ApproxByC2` ran (`iToApproxByC2`).
/// `BS` is discarded when the arm does not fire, exactly as in OCCT: `TC` stays
/// the trimmed curve.
fn smoothed_piece(
    cv: Arc<dyn Curve>,
    first: f64,
    last: f64,
    clean_tol: f64,
) -> Result<(Arc<dyn Curve>, i32), String> {
    let tc: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(cv, first, last));
    let mut bs = curve_to_bspline_curve(tc.as_ref(), ParameterisationType::TgtThetaOver2)
        .map_err(|e| format!("ChFi3d_PerformElSpine: GeomConvert::CurveToBSplineCurve: {e:?}"))?;
    curve_cleaner(&mut bs, clean_tol, 0);
    let a_continuity = tc.continuity();
    let b_is_smooth = is_smooth(tc.as_ref())?;
    if a_continuity < Shape::C2 as u8 && !b_is_smooth {
        let approx = approx_by_c2(tc.as_ref())
            .ok_or_else(|| "ChFi3d_PerformElSpine: ChFi3d_ApproxByC2 failed".to_string())?;
        Ok((Arc::new(approx), 1))
    } else {
        Ok((tc, 0))
    }
}

/// `ChFi3d_PerformElSpine(HES, Spine, Continuity, Tol, IsOffset = false)`
/// (`cxx:4797-5327`).
pub fn ch_fi3d_perform_el_spine(
    hes: &mut Arc<ElSpine>,
    spine: &SpineHandle,
    continuity: Shape,
    tol: f64,
    is_offset: bool,
) -> Result<(), String> {
    // `ChFiDS_ElSpine& ES = *HES;` (`cxx:4819`). `Arc::make_mut` gives the
    // write-through the callers expect: they hand in the list element and store
    // the result back (`ChFi3d_Builder_2.cxx:3267-3278`).
    let es = Arc::make_mut(hes);

    // `WF = ES.FirstParameter(); WL = ES.LastParameter(); Wrefdeb = WF;
    // Wreffin = WL; nwf = WF; nwl = WL;` (`cxx:4820-4825`).
    let mut wf = es.first_parameter();
    let mut wl = es.last_parameter();
    let mut wrefdeb = wf;
    let mut wreffin = wl;
    let mut nwf = wf;
    let mut nwl = wl;
    let nbed = spine_rd(spine, |s| s.nb_edges())?;
    let periodic = spine_rd(spine, |s| s.is_periodic())?;

    // `if (periodic) { ... } else { ... }` (`cxx:4821-4843`).
    let (if_, il) = if periodic {
        let period = spine_rd(spine, |s| s.period())?.ok_or_else(|| {
            "ChFi3d_PerformElSpine: ChFiDS_Spine::Period() on a non-periodic spine".to_string()
        })?;
        nwf = clib::in_period(wf, -tol, period - tol);
        let if_ = spine_rd(spine, |s| s.index(nwf, true))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->Index(nwf, true)".to_string())?;
        nwl = clib::in_period(wl, tol, period + tol);
        let mut il = spine_rd(spine, |s| s.index(nwl, false))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->Index(nwl, false)".to_string())?;
        if nwl < nwf + tol {
            il += nbed as i32;
        }
        (if_, il)
    } else {
        let if_ = spine_rd(spine, |s| s.index(wf, true))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->Index(WF, true)".to_string())?;
        let il = spine_rd(spine, |s| s.index(wl, false))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->Index(WL, false)".to_string())?;
        wrefdeb = spine_rd(spine, |s| s.first_parameter_of(if_))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->FirstParameter(IF)".to_string())?
            .max(wf);
        wreffin = spine_rd(spine, |s| s.last_parameter_of(il))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->LastParameter(IL)".to_string())?
            .min(wl);
        (if_, il)
    };

    // `Spine->D1(WF, PDeb, VrefDeb); Spine->D1(WL, PFin, VrefFin);
    // VrefDeb.Normalize(); VrefFin.Normalize();` (`cxx:4845-4848`).
    let (p_deb, mut vref_deb) = spine_wr(spine, |s| s.d1(wf))?
        .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->D1(WF)".to_string())?;
    let (p_fin, mut vref_fin) = spine_wr(spine, |s| s.d1(wl))?
        .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->D1(WL)".to_string())?;
    vref_deb.normalize();
    vref_fin.normalize();

    // Traitment de la premiere arete (`cxx:4844-4951`).
    let mut cepadur = false;
    let mut e_old = Edge::wrap(spine_edge(spine, if_, is_offset)?)
        .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->Edges(IF) is not an edge".to_string())?;
    let edge = e_old.clone();
    // `Bof = BRepLib::BuildCurve3d(E);` (`cxx:4848`): the result is overwritten
    // before it is read.
    let _ = build_curve3d(
        &edge,
        BUILD_CURVE3D_TOL,
        Shape::C1,
        BUILD_CURVE3D_MAX_DEGREE,
        BUILD_CURVE3D_MAX_SEGMENT,
    );
    // `const BRepAdaptor_Curve& edc = Spine->CurrentElementarySpine(IF);
    // tolpared = edc.Resolution(tol);` (`cxx:4849-4850`).
    let tolpared = spine_wr(spine, |s| s.current_elementary_spine(if_))?
        .map(|c| c.resolution(tol))
        .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->CurrentElementarySpine(IF)".to_string())?;
    // `Cv = BRep_Tool::Curve(E, First, Last);` (`cxx:4851`).
    let (mut first, mut last) = BRepTool::edge_parameters(&edge);
    let mut cv = BRepTool::edge_curve_world(&edge)
        .ok_or_else(|| "ChFi3d_PerformElSpine: BRep_Tool::Curve(IF) is null".to_string())?;
    // `if (ES.IsPeriodic()) { ... ES.AddVertexWithTangent(...); }`
    // (`cxx:4853-4866`).
    if es.is_periodic() {
        let par_for_el_spine = if edge.orientation().is_forward() {
            first
        } else {
            last
        };
        let ax = axis_at(cv.as_ref(), par_for_el_spine)?;
        es.add_vertex_with_tangent(ax);
    }
    // `urefdeb = Spine->FirstParameter(IF); checkdeb = (nwf > urefdeb); if
    // (checkdeb) Spine->Parameter(IF, nwf, pared, false);` (`cxx:4868-4872`).
    let urefdeb = spine_rd(spine, |s| s.first_parameter_of(if_))?
        .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->FirstParameter(IF)".to_string())?;
    let checkdeb = nwf > urefdeb;
    let mut pared = if checkdeb {
        spine_parameter(spine, if_, nwf)?
    } else {
        0.0
    };
    if edge.orientation().is_reversed() {
        // `cxx:4874-4906`.
        let sov = first;
        // `ReversedParameter` is taken on the curve as it was read, before the
        // handle is replaced.
        first = cv.reversed_parameter(last);
        last = cv.reversed_parameter(sov);
        pared = if checkdeb {
            cv.reversed_parameter(pared)
        } else {
            first
        };
        if first < pared {
            first = pared;
        }
        if il == if_ {
            let ureffin = spine_rd(spine, |s| s.last_parameter_of(il))?
                .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->LastParameter(IL)".to_string())?;
            let checkfin = nwl < ureffin;
            pared = if checkfin {
                cv.reversed_parameter(spine_parameter(spine, il, nwl)?)
            } else {
                last
            };
            if pared < last {
                last = pared;
            }
        }
        let reversed = cv.reversed();
        cv = Arc::from(reversed);
    } else {
        // `cxx:4907-4929`.
        if !checkdeb {
            pared = first;
        }
        if first < pared {
            first = pared;
        }
        if il == if_ {
            let ureffin = spine_rd(spine, |s| s.last_parameter_of(il))?
                .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->LastParameter(IL)".to_string())?;
            let checkfin = nwl < ureffin;
            pared = if checkfin {
                spine_parameter(spine, il, nwl)?
            } else {
                last
            };
            if pared < last {
                last = pared;
            }
        }
    }
    // `if (std::abs(Last - First) < tolpared) cepadur = true;` (`cxx:4934-4937`).
    if (last - first).abs() < tolpared {
        cepadur = true;
    }
    // Petite veru pour les cas ou un KPart a bouffe l arete sans parvenir a
    // terminer. On tire une droite (`cxx:4941-4977`).
    if cepadur {
        let first_par = spine_rd(spine, |s| s.first_parameter_of(1))?
            .ok_or_else(|| "ChFi3d_PerformElSpine: Spine->FirstParameter(1)".to_string())?;
        let last_par = spine_rd(spine, |s| s.last_parameter_of(nbed as i32))?
            .ok_or_else(|| {
                format!("ChFi3d_PerformElSpine: Spine->LastParameter({nbed})")
            })?;
        if wl < first_par + tol {
            let (_, vtemp) = es.last_point_and_tgt();
            // `olin.SetLinearForm(-WL, d.XYZ(), PFin.XYZ())`.
            let dir = GpDir::from_vec(&vtemp)
                .map_err(|_| "ChFi3d_PerformElSpine: gp_Dir of a zero tangent".to_string())?;
            let origin = p_fin
                .translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(-wl));
            es.set_curve(Arc::new(GeomLine::from_pnt_dir(origin, dir)));
        } else if wf > last_par - tol {
            let (_, vtemp) = es.first_point_and_tgt();
            let dir = GpDir::from_vec(&vtemp)
                .map_err(|_| "ChFi3d_PerformElSpine: gp_Dir of a zero tangent".to_string())?;
            let origin = p_deb
                .translated_vec(&GpVec::from_xyz(dir.xyz()).multiplied_scalar(-wf));
            es.set_curve(Arc::new(GeomLine::from_pnt_dir(origin, dir)));
        }
        return Ok(());
    }

    // `TC = new GeomTrimmedCurve(Cv, First, Last); BS = ...; CurveCleaner(...);
    // if (aContinuity < GeomAbs_C2 && !bIsSmooth) { ... }` (`cxx:4977-4993`).
    let clean_tol = (wl - wf).abs() * 1.0e-4;
    let (mut tc, mut i_to_approx_by_c2) = smoothed_piece(cv, first, last, clean_tol)?;

    // `GeomConvert_CompCurveToBSplineCurve Concat(TC, Convert_QuasiAngular);`
    // (`cxx:4995`).
    let mut concat =
        CompCurveToBSplineCurve::from_curve(tc.as_ref(), ParameterisationType::QuasiAngular)
            .map_err(|e| {
                format!("ChFi3d_PerformElSpine: GeomConvert_CompCurveToBSplineCurve: {e:?}")
            })?;

    // Concatenation des aretes suivantes (`cxx:4997-5093`).
    for iedge in (if_ + 1)..=il {
        let iloc = if periodic {
            if nbed == 0 {
                return Err("ChFi3d_PerformElSpine: periodic spine without edges".to_string());
            }
            (iedge - 1) % nbed as i32 + 1
        } else {
            iedge
        };
        let e = spine_edge(spine, iloc, is_offset)?;
        let edge = Edge::wrap(e)
            .ok_or_else(|| format!("ChFi3d_PerformElSpine: Spine->Edges({iloc}) not an edge"))?;
        // `if (BRep_Tool::Degenerated(E)) continue;` (`cxx:5009-5012`).
        if BRepTool::is_degenerated(&edge) {
            continue;
        }
        // `epsV = tol; if (TopExp::CommonVertex(Eold, E, V)) epsV =
        // BRep_Tool::Tolerance(V);` (`cxx:5014-5019`).
        let mut eps_v = tol;
        if let Some(v) = common_vertex(&e_old, &edge).and_then(Vertex::wrap) {
            eps_v = BRepTool::vertex_tolerance(&v);
        }
        // `Bof = BRepLib::BuildCurve3d(E); if (!Bof) throw;` (`cxx:5021-5025`).
        let bof = build_curve3d(
            &edge,
            BUILD_CURVE3D_TOL,
            Shape::C1,
            BUILD_CURVE3D_MAX_DEGREE,
            BUILD_CURVE3D_MAX_SEGMENT,
        );
        if !bof {
            return Err("PerformElSpine : BuildCurve3d error".to_string());
        }
        // `Cv = BRep_Tool::Curve(E, First, Last);` (`cxx:5027`).
        let (mut first, mut last) = BRepTool::edge_parameters(&edge);
        let mut cv = BRepTool::edge_curve_world(&edge)
            .ok_or_else(|| format!("ChFi3d_PerformElSpine: BRep_Tool::Curve({iloc}) is null"))?;
        // `double ParForElSpine = (E.Orientation() == TopAbs_FORWARD) ? First :
        // Last; Cv->D1(...); ES.AddVertexWithTangent(...);` (`cxx:5028-5031`).
        let par_for_el_spine = if edge.orientation().is_forward() {
            first
        } else {
            last
        };
        let ax = axis_at(cv.as_ref(), par_for_el_spine)?;
        es.add_vertex_with_tangent(ax);
        // `if (IEdge == IL) { ... }` (`cxx:5034-5062`).
        if iedge == il {
            let ureffin = spine_rd(spine, |s| s.last_parameter_of(iloc))?
                .ok_or_else(|| format!("ChFi3d_PerformElSpine: Spine->LastParameter({iloc})"))?;
            let checkfin = nwl < ureffin;
            let mut pared = if checkfin {
                spine_parameter(spine, iloc, nwl)?
            } else {
                last
            };
            if edge.orientation().is_reversed() {
                let sov = first;
                first = cv.reversed_parameter(last);
                last = cv.reversed_parameter(sov);
                pared = if checkfin {
                    cv.reversed_parameter(pared)
                } else {
                    last
                };
                let reversed = cv.reversed();
                cv = Arc::from(reversed);
            }
            if pared < last {
                last = pared;
            }
        }
        // `TC = ...; BS = ...; CurveCleaner(...); if (aContinuity < GeomAbs_C2
        // && !bIsSmooth) { ++iToApproxByC2; BS = ChFi3d_ApproxByC2(TC); TC = BS;
        // }` (`cxx:5064-5080`).
        let (piece, inc) = smoothed_piece(cv, first, last, (wl - wf).abs() * 1.0e-4)?;
        tc = piece;
        i_to_approx_by_c2 += inc;
        // `tolrac = std::min(tol, epsV); Bof = Concat.Add(TC, 2.*tolrac, true);
        // if (!Bof) Bof = Concat.Add(TC, 2.*epsV, true); if (!Bof) { Bof =
        // Concat.Add(TC, 200.*epsV, true); if (!Bof) throw; }`
        // (`cxx:5082-5093`). `Add` keeps the `WithRatio = true` and `MinM = 0`
        // defaults of `GeomConvert_CompCurveToBSplineCurve.hxx:56-60`.
        let tolrac = tol.min(eps_v);
        let mut added = add_piece(&mut concat, tc.as_ref(), 2.0 * tolrac)?;
        if !added {
            added = add_piece(&mut concat, tc.as_ref(), 2.0 * eps_v)?;
        }
        if !added {
            added = add_piece(&mut concat, tc.as_ref(), 200.0 * eps_v)?;
            if !added {
                return Err("PerformElSpine: spine merged error".to_string());
            }
        }
        e_old = edge;
    }

    // `BSpline = Concat.BSplineCurve();` (`cxx:5098`).
    let mut bspline = concat
        .into_curve()
        .ok_or_else(|| "ChFi3d_PerformElSpine: null concatenated B-spline".to_string())?;

    // `NCollection_Array1<double> BSNoeuds(BSpline->Knots());
    // BSplCLib::Reparametrize(Wrefdeb, Wreffin, BSNoeuds);
    // BSpline->SetKnots(BSNoeuds);` (`cxx:5100-5105`).
    let mut knots = bspline.distinct_knots().to_vec();
    reparametrize(wrefdeb, wreffin, &mut knots);
    bspline
        .set_knots(&knots)
        .map_err(|e| format!("ChFi3d_PerformElSpine: Geom_BSplineCurve::SetKnots: {e}"))?;

    // Traitement des Extremites (`cxx:5107-5197`).
    let mut caredeb = 0i32;
    let mut carefin = 0i32;
    let angle = std::f64::consts::PI * 0.75;
    let mut local_wf = wf;
    let mut local_wl = wl;
    // `Prolongement C3 au debut` (`cxx:5115-5152`).
    if !es.is_periodic() && !p_deb.is_equal_with_tolerance(bspline.pole(0), tol) {
        if bspline.is_rational() {
            caredeb = 1;
        }
        let rabdist = wrefdeb - wf;
        let bout = p_deb.translated_vec(&vref_deb.multiplied_scalar(-20.0 * rabdist));
        let mut goodext = false;
        let mut newc = bspline.clone();
        let mut icont = 3;
        while icont >= 1 && !goodext {
            newc = bspline.clone();
            // `GeomLib::ExtendCurveToPoint(anExtCurve, Bout, icont, false);`
            // (`cxx:5133`).
            extend_curve_to_point(&mut newc, &bout, icont, false)?;
            // `GCPnts_AbscissaPoint GCP(gacurve, -rabdist, Wrefdeb, WF);`
            // (`cxx:5136`), the four-argument constructor.
            if let Some(param) = abscissa_point_with_guess(&newc, -rabdist, wrefdeb, wf) {
                wf = param;
                goodext = good_ext(&newc, &vref_deb, wrefdeb, wf, angle);
            }
            icont -= 1;
        }
        if caredeb != 0 {
            caredeb = newc.distinct_knots().len() as i32 - bspline.distinct_knots().len() as i32;
        }
        bspline = newc;
        local_wf = bspline.first_parameter();
    }
    // `Prolongement C3 en fin` (`cxx:5154-5183`).
    if !es.is_periodic()
        && !p_fin.is_equal_with_tolerance(bspline.pole(bspline.nb_poles() - 1), tol)
    {
        if bspline.is_rational() {
            carefin = 1;
        }
        let rabdist = wl - wreffin;
        let bout = p_fin.translated_vec(&vref_fin.multiplied_scalar(20.0 * rabdist));
        let mut goodext = false;
        let mut newc = bspline.clone();
        let mut icont = 3;
        while icont >= 1 && !goodext {
            newc = bspline.clone();
            // `GeomLib::ExtendCurveToPoint(anExtCurve, Bout, icont, true);`
            // (`cxx:5164`).
            extend_curve_to_point(&mut newc, &bout, icont, true)?;
            // `GCPnts_AbscissaPoint GCP(gacurve, rabdist, Wreffin, WL);`
            // (`cxx:5167`).
            if let Some(param) = abscissa_point_with_guess(&newc, rabdist, wreffin, wl) {
                wl = param;
                goodext = good_ext(&newc, &vref_fin, wreffin, wl, angle);
            }
            icont -= 1;
        }
        if carefin != 0 {
            carefin = newc.distinct_knots().len() as i32 - bspline.distinct_knots().len() as i32;
        }
        bspline = newc;
        local_wl = bspline.last_parameter();
    }
    // Reparametrisation et segmentation sur le domaine de la Spine
    // (`cxx:5186-5198`).
    if (bspline.first_parameter() - wf).abs() < tol {
        wf = bspline.first_parameter();
    }
    if (bspline.last_parameter() - wl).abs() < tol {
        wl = bspline.last_parameter();
    }
    if local_wf < wf || local_wl > wl {
        // `BSpline->Segment(WF, WL);` (`cxx:5194`), whose default tolerance is
        // `Precision::PConfusion()` (`Geom_BSplineCurve.hxx:322-324`).
        bspline
            .segment(wf, wl, occt_core::precision::PCONFUSION)
            .map_err(|e| format!("ChFi3d_PerformElSpine: Geom_BSplineCurve::Segment: {e}"))?;
        es.set_first_parameter(wf);
        es.set_last_parameter(wl);
    }

    // `if (BSpline->IsRational()) { C1 = Copy;
    // GeomConvert::C0BSplineToC1BSplineCurve(C1, tol, 0.1);
    // if (C1->FirstParameter() == BSpline->FirstParameter()) BSpline = C1; }`
    // (`cxx:5199-5213`).
    if bspline.is_rational() {
        let c1 = c0_bspline_to_c1_bspline_curve(&bspline, tol, 0.1)
            .map_err(|e| format!("ChFi3d_PerformElSpine: C0BSplineToC1BSplineCurve: {e}"))?;
        // Il faut s'assurer que l'origine n'a pas bouge (cts21158).
        if c1.first_parameter() == bspline.first_parameter() {
            bspline = c1;
        }
    }

    // Correction C2 or C3 (if possible) (`cxx:5216-5269`).
    // `if ((caredeb || carefin) && BSpline->Degree() < 8) BSpline->IncreaseDegree(8);`
    // (`cxx:5217-5222`): an extended end gained knots, so the degree is raised
    // before the multiplicity limits below are taken from it.
    if (caredeb != 0 || carefin != 0) && bspline.degree() < 8 {
        bspline
            .increase_degree(8)
            .map_err(|e| format!("ChFi3d_PerformElSpine: IncreaseDegree(8): {e}"))?;
    }
    let mut fk = 2i32;
    let mut lk = bspline.distinct_knots().len() as i32 - 1;
    if bspline.periodic() {
        fk = 1;
    }
    if caredeb != 0 {
        fk += caredeb;
    }
    if carefin != 0 {
        lk -= carefin;
    }
    let mult_max;
    if continuity == Shape::C3 {
        if bspline.degree() < 7 {
            bspline
                .increase_degree(7)
                .map_err(|e| format!("ChFi3d_PerformElSpine: IncreaseDegree(7): {e}"))?;
        }
        mult_max = bspline.degree() as i32 - 3;
    } else {
        if bspline.degree() < 5 {
            bspline
                .increase_degree(5)
                .map_err(|e| format!("ChFi3d_PerformElSpine: IncreaseDegree(5): {e}"))?;
        }
        mult_max = bspline.degree() as i32 - 2;
    }
    curve_cleaner(&mut bspline, (wl - wf).abs() * 1.0e-4, 1);
    curve_cleaner(&mut bspline, (wl - wf).abs() * 1.0e-2, mult_max);
    let mult_min = (bspline.degree() as i32 - 4).max(1);
    for ii in fk..=lk {
        if knot_multiplicity(&bspline, ii) > mult_max {
            let _ = bspline.remove_knot(ii, mult_max, (wl - wf).abs() / 10.0);
        }
        // See C4
        if knot_multiplicity(&bspline, ii) > mult_min {
            let _ = bspline.remove_knot(ii, mult_min, (wl - wf).abs() * 1.0e-4);
        }
    }

    // `elspine periodic => BSpline Periodic` (`cxx:5270-5284`).
    if es.is_periodic() {
        if !bspline.periodic() {
            bspline.set_periodic();
            if i_to_approx_by_c2 != 0 {
                let _ = bspline.remove_knot(1, mult_max, (wl - wf).abs() / 10.0);
            }
        }
    } else {
        // Otherwise is it necessary to move the poles to adapt them to new
        // tangents ? (`cxx:5285-5318`).
        let mut adjust = false;
        let (p1, v1) = bspline.d1(wf);
        let v1 = v1.normalized();
        let (p_deb, vref_deb) = es.first_point_and_tgt();
        let scaldeb = vref_deb.dot(&v1);
        let disdeb = p_deb.distance(&p1);
        if (wf - local_wf).abs() < 1.0e-12 && (scaldeb <= 0.9999999 || disdeb >= tol) {
            adjust = true;
        }
        let (p2, v2) = bspline.d1(wl);
        let v2 = v2.normalized();
        let (p_fin, vref_fin) = es.last_point_and_tgt();
        let scalfin = vref_fin.dot(&v2);
        let disfin = p_fin.distance(&p2);
        if (wl - local_wl).abs() < 1.0e-12 && (scalfin <= 0.9999999 || disfin >= tol) {
            adjust = true;
        }
        if adjust {
            // `GeomLib::AdjustExtremity(anExtCurve, PDeb, PFin, VrefDeb,
            // VrefFin); BSpline = down_cast(anExtCurve);` (`cxx:5314-5316`).
            adjust_extremity(&mut bspline, &p_deb, &p_fin, &vref_deb, &vref_fin)?;
        }
    }

    // Le Resultat (`cxx:5320-5321`).
    es.set_curve(Arc::new(bspline) as Arc<dyn Curve>);
    Ok(())
}

/// `Concat.Add(TC, Tolerance, true)` (`cxx:5083`, `:5088`, `:5091`).
fn add_piece(
    concat: &mut CompCurveToBSplineCurve,
    curve: &dyn Curve,
    tolerance: f64,
) -> Result<bool, String> {
    concat
        .add(curve, tolerance, true, true, 0)
        .map_err(|e| format!("ChFi3d_PerformElSpine: Concat.Add: {e:?}"))
}
