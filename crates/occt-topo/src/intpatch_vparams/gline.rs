//! `IntPatch_GLine::ComputeVertexParameters`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_GLine.cxx:414-1102`.
//!
//! Kept verbatim, including the two places where the source reuses the S1 flag
//! in the S2 pass (`cxx:570` tests `IsVertexOnS1` inside the OnS2 loop) and the
//! `SortAgain` / `ToBreak` seam fix-point of BUC61033 (`cxx:624-670`).
//!
//! The `fipt` / `lapt` / `indf` / `indl` members of `IntPatch_GLine.hxx:198-201`
//! are `has_first_point` / `has_last_point` / `first_index` / `last_index` of
//! [`GLine`]. The source adjusts them with strict `>` (and `== i` remapping)
//! here, unlike `IntPatch_WLine::ComputeVertexParameters` which uses `>=`.

use occt_core::precision::PCONFUSION;

use crate::geom_int::{GLine, GLineKind};
use crate::int_tools_wline::PatchPoint;

const REAL_LAST: f64 = f64::MAX;
const REAL_FIRST: f64 = -f64::MAX;

/// `IntPatch_GLine::ComputeVertexParameters(Tol)` (`IntPatch_GLine.cxx:421`).
pub fn compute_vertex_parameters_gline(g: &mut GLine) {
    let prec_pc = PCONFUSION * 1000.0;
    let param_min_on_line = if g.has_first_point {
        g.vertices[(g.first_index.max(1) as usize) - 1].parameter_on_line()
    } else {
        -100000.0
    };
    let param_max_on_line = if g.has_last_point {
        g.vertices[(g.last_index.max(1) as usize) - 1].parameter_on_line()
    } else {
        100000.0
    };

    drop_domainless_duplicates_on_dom(g, prec_pc);
    drop_duplicate_arc_params(g, prec_pc, true);
    drop_duplicate_arc_params(g, prec_pc, false);
    seam_fix_point(
        g,
        prec_pc,
        param_min_on_line,
        param_max_on_line,
    );
}

/// `fipt` / `lapt` fix-up for the `Remove(j)` sites of `IntPatch_GLine.cxx:460-475`,
/// `:505-543`, `:574-610`.
fn remove_vertex(g: &mut GLine, index: i32) {
    g.vertices.remove((index - 1) as usize);
    if g.has_last_point && g.last_index > index {
        g.last_index -= 1;
    }
    if g.has_first_point && g.first_index > index {
        g.first_index -= 1;
    }
}

/// `IntPatch_GLine.cxx:442-479`: a vertex on no restriction that shares the
/// parameter of a vertex on a restriction is dropped.
fn drop_domainless_duplicates_on_dom(g: &mut GLine, prec_pc: f64) {
    let mut nbvtx = g.nb_vertex();
    loop {
        let mut a_point_deleted = false;
        let mut i = 1;
        while i <= nbvtx && !a_point_deleted {
            let vtx_i = g.vertex(i).clone();
            if vtx_i.is_on_dom_s1() || vtx_i.is_on_dom_s2() {
                let mut j = 1;
                while j <= nbvtx && !a_point_deleted {
                    if i != j {
                        let vtx_j = g.vertex(j).clone();
                        if !vtx_j.is_on_dom_s1()
                            && !vtx_j.is_on_dom_s2()
                            && (vtx_i.parameter_on_line() - vtx_j.parameter_on_line()).abs()
                                <= prec_pc
                        {
                            remove_vertex(g, j);
                            nbvtx -= 1;
                            a_point_deleted = true;
                        }
                    }
                    j += 1;
                }
            }
            i += 1;
        }
        if !(a_point_deleted && nbvtx > 2) {
            break;
        }
    }
}

/// `IntPatch_GLine.cxx:482-548` (S1) and `:550-617` (S2): two vertices on the
/// same restriction and only on it with the same arc parameter collapse; the
/// non-vertex one is the one that goes.
fn drop_duplicate_arc_params(g: &mut GLine, prec_pc: f64, on_first: bool) {
    let mut nbvtx = g.nb_vertex();
    loop {
        let mut a_point_deleted = false;
        let mut i = 1;
        while i <= nbvtx && !a_point_deleted {
            let vtx_i = g.vertex(i).clone();
            if vtx_i.is_on_dom_s1() == on_first && vtx_i.is_on_dom_s2() != on_first {
                let mut j = 1;
                while j <= nbvtx && !a_point_deleted {
                    if i != j {
                        let vtx_j = g.vertex(j).clone();
                        if vtx_j.is_on_dom_s1() == on_first
                            && vtx_j.is_on_dom_s2() != on_first
                            && (vtx_i.param_on_arc(on_first) - vtx_j.param_on_arc(on_first)).abs()
                                <= prec_pc
                            && vtx_i.arc_id(on_first) == vtx_j.arc_id(on_first)
                        {
                            // `cxx:570` reads `IsVertexOnS1` in both passes.
                            if vtx_i.is_vertex_on_s1() {
                                remove_vertex(g, j);
                                nbvtx -= 1;
                            } else {
                                remove_vertex(g, i);
                                nbvtx -= 1;
                            }
                            a_point_deleted = true;
                        }
                    }
                    j += 1;
                }
            }
            i += 1;
        }
        if !a_point_deleted {
            break;
        }
    }
}

/// Sort by `ParameterOnLine`, remapping `indf` / `indl` across the swap
/// (`IntPatch_GLine.cxx:630-664`).
fn sort_vertices_remap(g: &mut GLine, nbvtx: i32) {
    loop {
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            if g.vertices[(i - 2) as usize].parameter_on_line()
                > g.vertices[(i - 1) as usize].parameter_on_line()
            {
                sort_is_ok = false;
                g.vertices.swap((i - 2) as usize, (i - 1) as usize);
                if g.has_first_point {
                    if g.first_index == i {
                        g.first_index = i - 1;
                    } else if g.first_index == i - 1 {
                        g.first_index = i;
                    }
                }
                if g.has_last_point {
                    if g.last_index == i {
                        g.last_index = i - 1;
                    } else if g.last_index == i - 1 {
                        g.last_index = i;
                    }
                }
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// The `SortAgain` / `ToBreak` fix point of `IntPatch_GLine.cxx:620-1102`,
/// including the seam-vertex parameter shift for circles and ellipses
/// (`cxx:852-1078`).
fn seam_fix_point(g: &mut GLine, prec_pc: f64, param_min_on_line: f64, param_max_on_line: f64) {
    let mut u1min = REAL_LAST;
    let mut u1max = REAL_FIRST;
    let mut u2min = REAL_LAST;
    let mut u2max = REAL_FIRST;
    let mut to_break = false;
    let mut sort_again = true;
    let pi_pi = std::f64::consts::PI + std::f64::consts::PI;
    let is_periodic_arc = matches!(g.kind, GLineKind::Circ(_) | GLineKind::Elips(_));

    loop {
        let mut nbvtx = g.nb_vertex();
        if sort_again {
            sort_vertices_remap(g, nbvtx);
        }
        if to_break {
            break;
        }
        sort_again = false;
        let mut sort_is_ok = true;
        let mut i = 2;
        while i <= nbvtx && sort_is_ok {
            let mut j = 1;
            while j <= nbvtx && sort_is_ok {
                if i != j {
                    let mut vtx = g.vertex(i).clone();
                    let mut vtxm1 = g.vertex(j).clone();
                    let mut kill = false;
                    let mut killm1 = false;
                    if (vtxm1.parameter_on_line() - vtx.parameter_on_line()).abs() < prec_pc {
                        if vtxm1.is_on_dom_s1() && vtx.is_on_dom_s1() {
                            if vtxm1.arc_id(true) == vtx.arc_id(true) {
                                if vtxm1.is_on_dom_s2() {
                                    if !vtx.is_on_dom_s2() {
                                        kill = true;
                                    } else if vtxm1.arc_id(false) == vtx.arc_id(false) {
                                        if vtxm1.is_vertex_on_s2() {
                                            kill = true;
                                        } else {
                                            killm1 = true;
                                        }
                                    }
                                } else if vtx.is_on_dom_s2() {
                                    killm1 = true;
                                }
                            }
                        } else if !vtxm1.is_on_dom_s2() && !vtx.is_on_dom_s2() {
                            if vtxm1.is_on_dom_s1() && !vtx.is_on_dom_s1() {
                                kill = true;
                            } else if vtx.is_on_dom_s1() && !vtxm1.is_on_dom_s1() {
                                killm1 = true;
                            }
                        }

                        if !(kill || killm1) {
                            if vtxm1.is_on_dom_s2() && vtx.is_on_dom_s2() {
                                if vtxm1.arc_id(false) == vtx.arc_id(false) {
                                    if vtxm1.is_on_dom_s1() {
                                        if !vtx.is_on_dom_s1() {
                                            kill = true;
                                        } else if vtxm1.arc_id(true) == vtx.arc_id(true) {
                                            if vtxm1.is_vertex_on_s1() {
                                                kill = true;
                                            } else {
                                                killm1 = true;
                                            }
                                        }
                                    } else if vtx.is_on_dom_s1() {
                                        killm1 = true;
                                    }
                                }
                            } else if !vtxm1.is_on_dom_s1() && !vtx.is_on_dom_s1() {
                                if vtxm1.is_on_dom_s2() && !vtx.is_on_dom_s2() {
                                    kill = true;
                                } else if vtx.is_on_dom_s2() && !vtxm1.is_on_dom_s2() {
                                    killm1 = true;
                                }
                            }
                        }

                        if kill {
                            sort_is_ok = false;
                            remove_kill(g, i, j, true);
                            nbvtx -= 1;
                        } else if killm1 {
                            sort_is_ok = false;
                            remove_kill(g, i, j, false);
                            nbvtx -= 1;
                        } else if is_periodic_arc {
                            if !seam_shift(
                                &mut vtx,
                                &mut vtxm1,
                                nbvtx,
                                pi_pi,
                                prec_pc,
                                param_min_on_line,
                                param_max_on_line,
                                &mut u1min,
                                &mut u1max,
                                &mut u2min,
                                &mut u2max,
                                &mut to_break,
                                &mut sort_again,
                                &mut sort_is_ok,
                            ) {
                                // `continue` of `cxx:871`: the pair is skipped
                                // without writing anything back.
                                j += 1;
                                continue;
                            }
                            g.vertex_mut(i).set_parameter(vtx.parameter_on_line());
                            g.vertex_mut(j).set_parameter(vtxm1.parameter_on_line());
                        }
                    }
                }
                j += 1;
            }
            i += 1;
        }
        if sort_is_ok {
            break;
        }
    }

    recalc_first_last(g);
}

/// `cxx:791-851`: the `kill` / `killm1` removals with their own index rules.
fn remove_kill(g: &mut GLine, i: i32, j: i32, kill: bool) {
    if kill {
        if g.has_last_point {
            if g.last_index > i {
                g.last_index -= 1;
            } else if g.last_index == i {
                g.last_index = j;
            }
        }
        if g.has_first_point {
            if g.first_index > i {
                g.first_index -= 1;
            } else if g.first_index == i {
                g.first_index = j;
            }
        }
        g.vertices.remove((i - 1) as usize);
    } else {
        if g.has_last_point {
            if g.last_index > j {
                g.last_index -= 1;
            } else if g.last_index == j {
                g.last_index = i - 1;
            }
        }
        if g.has_first_point {
            if g.first_index > j {
                g.first_index -= 1;
            } else if g.first_index == j {
                g.first_index = i - 1;
            }
        }
        g.vertices.remove((j - 1) as usize);
    }
}

/// `IntPatch_GLine.cxx:1084-1101`: final re-sort and `indf = 1`, `indl = nbvtx`.
fn recalc_first_last(g: &mut GLine) {
    let nbvtx = g.nb_vertex();
    if nbvtx != 0 {
        loop {
            let mut sort_is_ok = true;
            for i in 2..=nbvtx {
                if g.vertices[(i - 2) as usize].parameter_on_line()
                    > g.vertices[(i - 1) as usize].parameter_on_line()
                {
                    sort_is_ok = false;
                    g.vertices.swap((i - 2) as usize, (i - 1) as usize);
                }
            }
            if sort_is_ok {
                break;
            }
        }
        g.last_index = nbvtx;
        g.first_index = 1;
    }
}

/// The seam-edge block of `IntPatch_GLine.cxx:852-1078`. Returns `false` for
/// the `continue` of `cxx:871`.
#[allow(clippy::too_many_arguments)]
fn seam_shift(
    vtx: &mut PatchPoint,
    vtxm1: &mut PatchPoint,
    nbvtx: i32,
    pi_pi: f64,
    prec_pc: f64,
    param_min_on_line: f64,
    param_max_on_line: f64,
    u1min: &mut f64,
    u1max: &mut f64,
    u2min: &mut f64,
    u2max: &mut f64,
    to_break: &mut bool,
    sort_again: &mut bool,
    sort_is_ok: &mut bool,
) -> bool {
    let ponline = vtx.parameter_on_line();
    let is2pi = (ponline - pi_pi).abs() <= prec_pc;
    if nbvtx > 2 && !is2pi {
        return false;
    }
    let new_param = if is2pi {
        0.0
    } else if ponline.abs() <= prec_pc {
        pi_pi
    } else {
        ponline - pi_pi
    };

    let (u1a, v1a, u2a, v2a) = vtxm1.parameters();
    let (u1b, v1b, u2b, v2b) = vtx.parameters();
    let mut flag = 0u32;
    if (u1a - u1b).abs() <= prec_pc {
        flag |= 1;
    }
    if (v1a - v1b).abs() <= prec_pc {
        flag |= 2;
    }
    if (u2a - u2b).abs() <= prec_pc {
        flag |= 4;
    }
    if (v2a - v2b).abs() <= prec_pc {
        flag |= 8;
    }

    let mut test_on1 = false;
    let mut test_on2 = false;
    match flag {
        3 | 7 | 10 | 12 | 13 => {}
        11 => test_on2 = true,
        14 => test_on1 = true,
        _ => {}
    }

    if test_on1 {
        let u1_a = if u1a < u1b { u1a } else { u1b };
        let u1_b = if u1a < u1b { u1b } else { u1a };
        if *u1min == REAL_LAST {
            *u1min = u1_a;
            *u1max = u1_b;
        } else {
            if (u1_a - *u1min).abs() > prec_pc {
                *to_break = true;
            }
            if (u1_b - *u1max).abs() > prec_pc {
                *to_break = true;
            }
        }
        if new_param >= param_min_on_line && new_param <= param_max_on_line {
            *sort_again = true;
            *sort_is_ok = false;
            if new_param > ponline {
                if u1a < u1b {
                    vtx.set_parameter(new_param);
                } else {
                    vtxm1.set_parameter(new_param);
                }
            } else if u1a > u1b {
                vtx.set_parameter(new_param);
            } else {
                vtxm1.set_parameter(new_param);
            }
        }
    }

    if test_on2 {
        let u2_a = if u2a < u2b { u2a } else { u2b };
        let u2_b = if u2a < u2b { u2b } else { u2a };
        if *u2min == REAL_LAST {
            *u2min = u2_a;
            *u2max = u2_b;
        } else {
            if (u2_a - *u2min).abs() > prec_pc {
                *to_break = true;
            }
            if (u2_b - *u2max).abs() > prec_pc {
                *to_break = true;
            }
        }
        if new_param >= param_min_on_line && new_param <= param_max_on_line {
            *sort_again = true;
            *sort_is_ok = false;
            if new_param > ponline {
                if u2a < u2b {
                    vtx.set_parameter(new_param);
                } else {
                    vtxm1.set_parameter(new_param);
                }
            } else if u2a > u2b {
                vtx.set_parameter(new_param);
            } else {
                vtxm1.set_parameter(new_param);
            }
        }
    }
    true
}
