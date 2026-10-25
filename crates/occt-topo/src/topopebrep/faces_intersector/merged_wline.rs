//! `TopOpeBRep_FacesIntersector`'s arc search and walking-line merge helpers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesIntersector.cxx:1357-1752`.
//!
//! Translation notes:
//! - `Extrema_ExtPS` is [`ExtPs`]; `Extrema_ExtPC` (`cxx:1391`) is
//!   [`occt_geom::extrema_pc::point_curve_extrema_all`], which minimises the
//!   same point-to-curve distance over the curve's own range that
//!   `GeomAdaptor_Curve::Load` gives.
//! - `BRep_Tool::Curve(E, f, l)` is [`CurveAdaptor::of_edge`]'s
//!   `Curve()` / `FirstParameter()` / `LastParameter()`.
//! - `TopExp::Vertices(E, V1, V2)` (`cxx:1490`) is the edge's two direct vertex
//!   children, in order.
//! - The arc scan of step 2 (`cxx:1425-1577`) re-runs the whole domain
//!   iteration and skips every arc but `ArcNumber`; the port does the same.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtFlag};
use occt_geom::Surface;
use occt_geom2d::geom2d_int::find_parameter;

use crate::bop_occt_util::vertices_of_edge_direct;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{PatchPoint, PntOn2S, WLine};
use crate::intpatch_line::IntPatchLine;
use crate::intpatch_vparams::wline::compute_vertex_parameters;
use crate::intsurf_transition::IntSurfTypeTrans;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// `IsPointOK(Extrema_ExtPS, P, Tol)` (`cxx:1631-1652`): true when the closest
/// point of the surface is within `tol`.
pub fn is_point_ok(ext_ps: &mut ExtPs<'_>, test_pnt: &GpPnt, tol: f64) -> bool {
    ext_ps.perform(test_pnt);
    if ext_ps.is_done() && ext_ps.nb_ext() > 0 {
        let mut min_dist2 = 1.0e+200;
        for i in 0..ext_ps.nb_ext() {
            let d2 = ext_ps.square_distance(i);
            if d2 < min_dist2 {
                min_dist2 = d2;
            }
        }
        if min_dist2 <= tol * tol {
            return true;
        }
    }
    false
}

/// `GetPointOn2S(Extrema_ExtPS, P, Tol, ResultPoint)` (`cxx:1658-1685`): the
/// closest surface parameters, when within `tol`.
pub fn get_point_on_2s(ext_ps: &mut ExtPs<'_>, test_pnt: &GpPnt, tol: f64) -> Option<(f64, f64)> {
    ext_ps.perform(test_pnt);
    if ext_ps.is_done() && ext_ps.nb_ext() > 0 {
        let mut minext = 0usize;
        let mut min_dist2 = 1.0e+200;
        for i in 0..ext_ps.nb_ext() {
            let d2 = ext_ps.square_distance(i);
            if d2 < min_dist2 {
                minext = i;
                min_dist2 = d2;
            }
        }
        if min_dist2 <= tol * tol {
            let (u, v, _) = ext_ps.point(minext);
            return Some((u, v));
        }
    }
    None
}

/// `GetArc(Slin, RankS, SurfaceObj, DomainObj, SurfaceTool, TestPoint, VrtxTol,
/// LineOn2S, First, Last)` (`cxx:1357-1623`): the arc of `DomainObj` that every
/// walking segment sits on, plus the samples taken from it. Returns 0 when no
/// such arc exists.
#[allow(clippy::too_many_arguments)]
pub fn get_arc(
    slin: &[GeomIntLine],
    rank_s: i32,
    surface_obj: &dyn Surface,
    domain_obj: &mut TopolTool,
    surface_tool: &dyn Surface,
    test_point: &GpPnt,
    vrtx_tol: &mut f64,
    line_on2_s: &mut Vec<PntOn2S>,
    first: &mut f64,
    last: &mut f64,
) -> i32 {
    // 1. find the number of the arc (edge) the WLine segments sit on.
    let mut min_distance2 = 1.0e+200;
    let mut first_es1 = 0.0;
    let mut last_es1 = 0.0;
    let mut arc_number = 0;
    let mut cur_arc = 0;

    *first = 0.0;
    *last = 0.0;

    domain_obj.init();
    while domain_obj.more() {
        cur_arc += 1;
        if let Some(an_e) = domain_obj.edge().cloned() {
            // `cxx:1387-1390`: a null 3D curve means e.g. a degenerated edge.
            if let Some(a_c_edge) = BRepTool::edge_curve(&an_e) {
                let extrema = point_curve_extrema_all(a_c_edge.as_ref(), test_point);
                for pair in &extrema {
                    let d2 = pair.distance * pair.distance;
                    if d2 < min_distance2 {
                        min_distance2 = d2;
                        arc_number = cur_arc;
                    }
                }
            }
        }
        domain_obj.next();
    }

    if arc_number == 0 {
        return 0;
    }

    // 2. load the parameters of the found edge and its arc.
    cur_arc = 0;
    let mut points_from_arc: Vec<GpPnt> = Vec::new();
    let tol = 1.0e-7;
    let mut wl_vertex_parameters: Vec<f64> = Vec::new();
    let mut classify_ok = true;
    let mut check_tol = 1.0e-5;

    domain_obj.init();
    while domain_obj.more() {
        cur_arc += 1;
        if cur_arc != arc_number {
            domain_obj.next();
            continue;
        }

        let Some(arc) = domain_obj.value().cloned() else {
            domain_obj.next();
            continue;
        };
        let Some(arc_curve) = arc.curve().cloned() else {
            domain_obj.next();
            continue;
        };

        for line in slin.iter() {
            let Some(wline) = line.as_wline() else {
                continue;
            };
            let nbpnts = wline.nb_pnts();
            let p_on2s_f = wline.point(1);
            let p_on2s_l = wline.point(nbpnts);
            let (upf, vpf, upl, vpl) = if rank_s == 1 {
                let (a, b) = p_on2s_f.parameters_on_s1();
                let (c, d) = p_on2s_l.parameters_on_s1();
                (a, b, c, d)
            } else {
                let (a, b) = p_on2s_f.parameters_on_s2();
                let (c, d) = p_on2s_l.parameters_on_s2();
                (a, b, c, d)
            };
            let a_p_on_line_f = GpPnt2d::new(upf, vpf);
            let a_p_on_line_l = GpPnt2d::new(upl, vpl);
            let par_f = find_parameter(arc_curve.as_ref(), &a_p_on_line_f, tol);
            let par_l = find_parameter(arc_curve.as_ref(), &a_p_on_line_l, tol);
            wl_vertex_parameters.push(par_f);
            wl_vertex_parameters.push(par_l);
        }

        // `cxx:1469-1485`: the exchange pass of the source, verbatim (it is
        // not a library sort; it keeps comparing against the running slot).
        let n = wl_vertex_parameters.len();
        for i in 0..n {
            for j in i..n {
                if j == i {
                    continue;
                }
                if wl_vertex_parameters[i] > wl_vertex_parameters[j] {
                    wl_vertex_parameters.swap(i, j);
                }
            }
        }

        // `cxx:1485-1496`. OCCT dereferences the edge unconditionally here; the
        // guard only replaces the null dereference, the rest is unchanged.
        let Some(an_e) = domain_obj.edge().cloned() else {
            domain_obj.next();
            continue;
        };
        let verts = vertices_of_edge_direct(an_e.shape());
        let v1_tol = verts.first().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
        let v2_tol = verts.get(1).map(BRepTool::vertex_tolerance).unwrap_or(0.0);
        let max_vertex_tol = if v1_tol > v2_tol { v1_tol } else { v2_tol };
        // `cxx:1492-1493`: `theVrtxTol` is an out-parameter of the function.
        *vrtx_tol = max_vertex_tol;
        let edge_tol = BRepTool::edge_tolerance(&an_e);
        check_tol = if max_vertex_tol > edge_tol {
            max_vertex_tol
        } else {
            edge_tol
        };

        let curve_adaptor = CurveAdaptor::of_edge(&an_e);
        let Some(a_c_edge) = curve_adaptor.curve().cloned() else {
            domain_obj.next();
            continue;
        };
        first_es1 = curve_adaptor.first_parameter();
        last_es1 = curve_adaptor.last_parameter();

        // `cxx:1499-1507`: the projector used by the classification checks.
        {
            let (u0, u1) = surface_tool.u_range();
            let (v0, v1) = surface_tool.v_range();
            let mut an_ext_ps_tool = ExtPs::new();
            an_ext_ps_tool.initialize(surface_tool, u0, u1, v0, v1, check_tol, check_tol);
            an_ext_ps_tool.set_flag(ExtremaExtFlag::Min);

            // Empty only if the sequence holds no walking line; OCCT reads
            // `Value(1)` out of range there, the port reports "no arc".
            if wl_vertex_parameters.is_empty() {
                classify_ok = false;
            }

            // a. min - first (`cxx:1509-1520`)
            if classify_ok {
                let gap = (first_es1 - wl_vertex_parameters[0]).abs();
                if gap > arc.resolution(max_vertex_tol).unwrap_or(0.0) {
                    let param = (first_es1 + wl_vertex_parameters[0]) / 2.0;
                    let point = a_c_edge.d0(param);
                    if !is_point_ok(&mut an_ext_ps_tool, &point, check_tol) {
                        classify_ok = false;
                    }
                }
            }
            // b. max - last (`cxx:1521-1533`)
            if classify_ok {
                let last_par = *wl_vertex_parameters.last().unwrap();
                if (last_es1 - last_par).abs() > arc.resolution(max_vertex_tol).unwrap_or(0.0) {
                    let param = (last_es1 + last_par) / 2.0;
                    let point = a_c_edge.d0(param);
                    if !is_point_ok(&mut an_ext_ps_tool, &point, check_tol) {
                        classify_ok = false;
                    }
                }
            }
            // c. all middle gaps (`cxx:1534-1550`)
            if classify_ok {
                let nb_chk_pnts = wl_vertex_parameters.len() / 2 - 1;
                for i in 1..=nb_chk_pnts {
                    let a = wl_vertex_parameters[i * 2 - 1];
                    let b = wl_vertex_parameters[i * 2];
                    if (b - a).abs() > arc.resolution(max_vertex_tol).unwrap_or(0.0) {
                        let param = (a + b) / 2.0;
                        let point = a_c_edge.d0(param);
                        if !is_point_ok(&mut an_ext_ps_tool, &point, check_tol) {
                            classify_ok = false;
                            break;
                        }
                    }
                }
            }
        }

        if !classify_ok {
            break;
        }

        // `cxx:1557-1575`: the 100 samples of the arc become the merged line.
        let param_first = wl_vertex_parameters[0];
        let param_last = *wl_vertex_parameters.last().unwrap();
        let d_param = (param_last - param_first).abs() / 100.0;
        let mut c_param = param_first;
        for i in 0..100 {
            if i == 99 {
                c_param = param_last;
            }
            points_from_arc.push(a_c_edge.d0(c_param));
            c_param += d_param;
        }
        *first = param_first;
        *last = param_last;

        domain_obj.next();
    }

    if !classify_ok {
        return 0;
    }

    // `cxx:1584-1622`: project the arc samples onto both surfaces.
    let (o_u0, o_u1) = surface_obj.u_range();
    let (o_v0, o_v1) = surface_obj.v_range();
    let mut an_ext_ps_obj = ExtPs::new();
    an_ext_ps_obj.initialize(surface_obj, o_u0, o_u1, o_v0, o_v1, check_tol, check_tol);
    an_ext_ps_obj.set_flag(ExtremaExtFlag::Min);

    let (t_u0, t_u1) = surface_tool.u_range();
    let (t_v0, t_v1) = surface_tool.v_range();
    let mut an_ext_ps_tool2 = ExtPs::new();
    an_ext_ps_tool2.initialize(surface_tool, t_u0, t_u1, t_v0, t_v1, check_tol, check_tol);
    an_ext_ps_tool2.set_flag(ExtremaExtFlag::Min);

    for arcpoint in &points_from_arc {
        let is_on_s1 = get_point_on_2s(&mut an_ext_ps_obj, arcpoint, check_tol);
        let is_on_s2 = get_point_on_2s(&mut an_ext_ps_tool2, arcpoint, check_tol);
        if let (Some((u1, v1)), Some((u2, v2))) = (is_on_s1, is_on_s2) {
            line_on2_s.push(PntOn2S {
                p: *arcpoint,
                u1,
                v1,
                u2,
                v2,
            });
        }
    }

    arc_number
}

/// `GetMergedWLineOnRestriction(Slin, VrtxTol, LineOn2S)` (`cxx:1692-1752`):
/// one walking line over the merged samples, carrying only the two end
/// vertices.
pub fn get_merged_wline_on_restriction(
    slin: &[GeomIntLine],
    vrtx_tol: f64,
    line_on2_s: &[PntOn2S],
) -> Option<WLine> {
    if line_on2_s.is_empty() {
        return None;
    }

    let mut trans1 = IntSurfTypeTrans::Undecided;
    let mut trans2 = IntSurfTypeTrans::Undecided;

    for line in slin.iter() {
        let Some(wline) = line.as_wline() else {
            continue;
        };
        if wline.base.transition_on_s1() != IntSurfTypeTrans::Undecided
            && wline.base.transition_on_s1() != IntSurfTypeTrans::Touch
        {
            trans1 = wline.base.transition_on_s1();
        }
        if wline.base.transition_on_s2() != IntSurfTypeTrans::Undecided
            && wline.base.transition_on_s2() != IntSurfTypeTrans::Touch
        {
            trans2 = wline.base.transition_on_s2();
        }
    }

    // `new IntPatch_WLine(theLineOn2S, false, trans1, trans2)` (`cxx:1727`).
    let mut m_wline = WLine::new();
    m_wline.base = IntPatchLine::with_transitions(false, trans1, trans2);
    m_wline.points = line_on2_s.to_vec();

    let nb_pnts = m_wline.nb_pnts();
    if nb_pnts < 1 {
        return None;
    }

    let mut a_first_vertex = PatchPoint::new(m_wline.point(1).value(), 0.0, 0.0, 0.0, 0.0, 0.0);
    a_first_vertex.set_value(m_wline.point(1).value(), vrtx_tol, false);
    let mut a_last_vertex =
        PatchPoint::new(m_wline.point(nb_pnts).value(), 0.0, 0.0, 0.0, 0.0, 0.0);
    a_last_vertex.set_value(m_wline.point(nb_pnts).value(), vrtx_tol, false);

    let (u1, v1, u2, v2) = m_wline.point(1).parameters();
    a_first_vertex.set_parameters(u1, v1, u2, v2);

    let (u1, v1, u2, v2) = m_wline.point(nb_pnts).parameters();
    a_last_vertex.set_parameters(u1, v1, u2, v2);

    a_first_vertex.set_parameter(1.0);
    a_last_vertex.set_parameter(line_on2_s.len() as f64);

    m_wline.vertices.push(a_first_vertex);
    m_wline.vertices.push(a_last_vertex);

    compute_vertex_parameters(&mut m_wline, vrtx_tol);

    Some(m_wline)
}
