//! `IntPatch_ALine::ComputeVertexParameters`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_ALine.cxx:77-676`.
//!
//! `PCONFUSION` of the source is the file-local `#define PCONFUSION 0.00001`
//! (`IntPatch_ALine.cxx:70`), i.e. `Precision::Confusion() * 100`, not
//! `Precision::PConfusion()`; it is kept as a local constant.
//!
//! `ParamMinOnLine` / `ParamMaxOnLine` come from `FirstParameter(OpenFirst)` /
//! `LastParameter(OpenLast)` (`IntPatch_ALine.lxx:38-51`), which read
//! `curv.Domain` and `curv.IsFirstOpen` / `IsLastOpen`; the port reads the same
//! two flags through `Curve::is_first_open` / `Curve::is_last_open`.

use occt_core::gp::GpPnt;

use crate::geom_int::ALine;
use crate::int_tools_wline::PatchPoint;

/// `#define PCONFUSION 0.00001` (`IntPatch_ALine.cxx:70`).
const PCONFUSION: f64 = 0.00001;

/// `IntPatch_ALine::ComputeVertexParameters(Tol)` (`IntPatch_ALine.cxx:77`).
pub fn compute_vertex_parameters_aline(a: &mut ALine, tol: f64) {
    let param_min_on_line = a.curve.first_parameter();
    let param_max_on_line = a.curve.last_parameter();
    let open_first = a.curve.is_first_open();
    let open_last = a.curve.is_last_open();

    duplicate_period_shifted_vertices(a, tol, param_min_on_line, param_max_on_line);

    let nbvtx = a.vertices.len() as i32;
    if nbvtx <= 0 {
        return;
    }
    sort_vertices(a, nbvtx);
    drop_duplicate_arc_params(a, true, nbvtx);
    drop_duplicate_arc_params(a, false, nbvtx);
    merge_fix_point(a);
    close_periodic_domain(a, tol, open_first, open_last, param_min_on_line, param_max_on_line);
    drop_domainless_ends(a);
    identify_same_parameter_points(a);
}

/// `IntPatch_ALine.cxx:102-132`: a vertex near the seam gets its +/- 2*pi
/// representation when the line is evaluated there.
fn duplicate_period_shifted_vertices(
    a: &mut ALine,
    tol: f64,
    param_min_on_line: f64,
    param_max_on_line: f64,
) {
    let nbvtx = a.vertices.len() as i32;
    let pi_pi = std::f64::consts::PI + std::f64::consts::PI;
    for i in 1..=nbvtx {
        let vtx = a.vertex(i).clone();
        let p = vtx.parameter_on_line();
        let pmpimpi = p - pi_pi;
        if pmpimpi >= param_min_on_line {
            let p1: GpPnt = a.value(pmpimpi);
            if p1.distance(&vtx.p) < tol {
                let mut ovtx = vtx.clone();
                ovtx.set_parameter(pmpimpi);
                a.vertices.push(ovtx);
            }
        }
        let pmpimpi = p + pi_pi;
        if pmpimpi <= param_max_on_line {
            let p1: GpPnt = a.value(pmpimpi);
            if p1.distance(&vtx.p) < tol {
                let mut ovtx = vtx.clone();
                ovtx.set_parameter(pmpimpi);
                a.vertices.push(ovtx);
            }
        }
    }
}

/// `IntPatch_ALine.cxx:139-173`: the sort remaps `indf` / `indl`.
fn sort_vertices(a: &mut ALine, nbvtx: i32) {
    loop {
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            if a.vertices[(i - 2) as usize].parameter_on_line()
                > a.vertices[(i - 1) as usize].parameter_on_line()
            {
                sort_is_ok = false;
                a.vertices.swap((i - 2) as usize, (i - 1) as usize);
                if a.has_first_point {
                    if a.first_index == i {
                        a.first_index = i - 1;
                    } else if a.first_index == i - 1 {
                        a.first_index = i;
                    }
                }
                if a.has_last_point {
                    if a.last_index == i {
                        a.last_index = i - 1;
                    } else if a.last_index == i - 1 {
                        a.last_index = i;
                    }
                }
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// `IntPatch_ALine.cxx:175-268`: two vertices on the same restriction and only
/// on it with the same arc parameter and line parameter collapse.
fn drop_duplicate_arc_params(a: &mut ALine, on_first: bool, nbvtx: i32) {
    let mut nbvtx = nbvtx;
    loop {
        let mut a_point_deleted = false;
        let mut i = 1;
        while i <= nbvtx && !a_point_deleted {
            let vtx_i = a.vertex(i).clone();
            if vtx_i.is_on_dom_s1() == on_first && vtx_i.is_on_dom_s2() != on_first {
                let mut j = 1;
                while j <= nbvtx && !a_point_deleted {
                    if i != j {
                        let vtx_j = a.vertex(j).clone();
                        if vtx_j.is_on_dom_s1() == on_first
                            && vtx_j.is_on_dom_s2() != on_first
                            && (vtx_i.param_on_arc(on_first) - vtx_j.param_on_arc(on_first)).abs()
                                <= PCONFUSION
                            && (vtx_i.parameter_on_line() - vtx_j.parameter_on_line()).abs()
                                <= PCONFUSION
                            && vtx_i.arc_id(on_first) == vtx_j.arc_id(on_first)
                        {
                            remove_vertex(a, j);
                            nbvtx -= 1;
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

/// `IntPatch_ALine.cxx:271-471`: the `SortAgain` kill loop. The inner index goes
/// `j < i` here (unlike `IntPatch_GLine`, which scans every `j != i`).
fn merge_fix_point(a: &mut ALine) {
    let mut sort_again = true;
    let mut nbvtx: i32;
    loop {
        nbvtx = a.vertices.len() as i32;
        if sort_again {
            sort_vertices(a, nbvtx);
        }
        sort_again = false;
        let mut sort_is_ok = true;
        let mut i = 2;
        while i <= nbvtx && sort_is_ok {
            let mut j = 1;
            while j < i && sort_is_ok {
                let vtx = a.vertex(i).clone();
                let vtxm1 = a.vertex(j).clone();
                let mut kill = false;
                let mut killm1 = false;
                if (vtxm1.parameter_on_line() - vtx.parameter_on_line()).abs() < PCONFUSION {
                    if vtxm1.is_on_dom_s1() && vtx.is_on_dom_s1() {
                        if vtxm1.arc_id(true) == vtx.arc_id(true) {
                            if vtxm1.is_on_dom_s2() {
                                if !vtx.is_on_dom_s2() {
                                    kill = true;
                                } else if vtxm1.arc_id(false) == vtx.arc_id(false) {
                                    kill = true;
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
                                        kill = true;
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
                        remove_kill(a, i, j, true);
                        nbvtx -= 1;
                    } else if killm1 {
                        sort_is_ok = false;
                        remove_kill(a, i, j, false);
                        nbvtx -= 1;
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
}

/// `svtx.Remove(index)` of `IntPatch_ALine.cxx:197-213`, `:245-260`.
fn remove_vertex(a: &mut ALine, index: i32) {
    a.vertices.remove((index - 1) as usize);
    if a.has_last_point && a.last_index >= index {
        a.last_index -= 1;
    }
    if a.has_first_point && a.first_index >= index {
        a.first_index -= 1;
    }
}

/// `kill` / `killm1` removals of `IntPatch_ALine.cxx:411-466`.
fn remove_kill(a: &mut ALine, i: i32, j: i32, kill: bool) {
    if kill {
        if a.has_last_point {
            if a.last_index > i {
                a.last_index -= 1;
            } else if a.last_index == i {
                a.last_index = j;
            }
        }
        if a.has_first_point {
            if a.first_index > i {
                a.first_index -= 1;
            } else if a.first_index == i {
                a.first_index = j;
            }
        }
        a.vertices.remove((i - 1) as usize);
    } else {
        if a.has_last_point {
            if a.last_index > j {
                a.last_index -= 1;
            } else if a.last_index == j {
                a.last_index = i - 1;
            }
        }
        if a.has_first_point {
            if a.first_index > j {
                a.first_index -= 1;
            } else if a.first_index == j {
                a.first_index = i - 1;
            }
        }
        a.vertices.remove((j - 1) as usize);
    }
}

/// `IntPatch_ALine.cxx:476-513`: on a closed periodic line, add the missing
/// endpoint representation of the first vertex.
fn close_periodic_domain(
    a: &mut ALine,
    tol: f64,
    open_first: bool,
    open_last: bool,
    param_min_on_line: f64,
    param_max_on_line: f64,
) {
    if open_first || open_last {
        return;
    }
    let nbvtx = a.vertices.len() as i32;
    if nbvtx < 1 {
        return;
    }
    let vtx0 = a.vertices[0].clone();
    let vtxn = a.vertices[(nbvtx - 1) as usize].clone();
    if vtx0.parameter_on_line() == param_min_on_line {
        if vtxn.parameter_on_line() != param_max_on_line {
            let pn: GpPnt = a.value(param_max_on_line);
            if pn.distance(&vtx0.p) <= tol {
                let mut ovtx = vtx0;
                ovtx.set_parameter(param_max_on_line);
                a.vertices.push(ovtx);
            }
        } else if vtxn.parameter_on_line() == param_max_on_line
            && vtx0.parameter_on_line() != param_min_on_line
        {
            let p0: GpPnt = a.value(param_min_on_line);
            if p0.distance(&vtx0.p) <= tol {
                let mut ovtx = vtxn;
                ovtx.set_parameter(param_min_on_line);
                a.vertices.insert(0, ovtx);
            }
        }
    }
}

/// `IntPatch_ALine.cxx:516-543`: drop the first / last vertex when it lies on
/// no restriction.
fn drop_domainless_ends(a: &mut ALine) {
    let mut nbvtx = a.vertices.len() as i32;
    if nbvtx > 1 {
        let vtx0 = a.vertices[0].clone();
        if !vtx0.is_on_dom_s1() && !vtx0.is_on_dom_s2() {
            a.vertices.remove(0);
            nbvtx -= 1;
            if a.has_last_point {
                a.last_index -= 1;
            }
        }
    }
    if nbvtx > 1 {
        let vtx0 = a.vertices[(nbvtx - 1) as usize].clone();
        if !vtx0.is_on_dom_s1() && !vtx0.is_on_dom_s2() {
            a.vertices.remove((nbvtx - 1) as usize);
            if a.has_last_point {
                a.last_index -= 1;
            }
        }
    }
}

/// `IntPatch_ALine.cxx:545-676`: two vertices with the same parameter get the
/// arc / vertex representations merged and, when their 3D points differ, the
/// later one is moved onto the earlier one.
fn identify_same_parameter_points(a: &mut ALine) {
    loop {
        let nbvtx = a.vertices.len() as i32;
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            let mut vtx = a.vertex(i).clone();
            let mut vtxm1 = a.vertex(i - 1).clone();
            if (vtx.parameter_on_line() - vtxm1.parameter_on_line()).abs() < PCONFUSION {
                merge_arc_state(&mut vtx, &mut vtxm1);
                if vtx.p.square_distance(&vtxm1.p) > 1e-12 {
                    let copy_vtx = vtxm1.clone();
                    vtxm1.set_parameter(vtx.parameter_on_line());
                    vtxm1.p = vtx.p;
                    vtxm1.set_tolerance(vtx.tolerance);
                    let (u1, v1, u2, v2) = vtx.parameters();
                    vtxm1.set_parameters(u1, v1, u2, v2);
                    copy_arc_if_on_dom(&mut vtxm1, &copy_vtx, true);
                    copy_arc_if_on_dom(&mut vtxm1, &copy_vtx, false);
                    if copy_vtx.is_vertex_on_s1() {
                        vtxm1.set_vertex_flag(true);
                        copy_arc_side(&mut vtxm1, &copy_vtx, true);
                    }
                    if copy_vtx.is_vertex_on_s2() {
                        vtxm1.set_vertex_flag(false);
                        copy_arc_side(&mut vtxm1, &copy_vtx, false);
                    }
                    sort_is_ok = false;
                }
                a.vertices[(i - 1) as usize] = vtx;
                a.vertices[(i - 2) as usize] = vtxm1;
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// `IntPatch_ALine.cxx:555-624`: propagate the arc and vertex state between the
/// two same-parameter vertices, both directions.
fn merge_arc_state(vtx: &mut PatchPoint, vtxm1: &mut PatchPoint) {
    for on_first in [true, false] {
        let (later_has, earlier_has) = if on_first {
            (vtx.is_on_dom_s1(), vtxm1.is_on_dom_s1())
        } else {
            (vtx.is_on_dom_s2(), vtxm1.is_on_dom_s2())
        };
        if later_has && !earlier_has {
            copy_arc_side(vtxm1, vtx, on_first);
        } else if earlier_has && !later_has {
            copy_arc_side(vtx, vtxm1, on_first);
        }
        let (later_vtx, earlier_vtx) = if on_first {
            (vtx.is_vertex_on_s1(), vtxm1.is_vertex_on_s1())
        } else {
            (vtx.is_vertex_on_s2(), vtxm1.is_vertex_on_s2())
        };
        if later_vtx && !earlier_vtx {
            vtxm1.set_vertex_flag(on_first);
            copy_arc_side(vtxm1, vtx, on_first);
        } else if earlier_vtx && !later_vtx {
            vtx.set_vertex_flag(on_first);
            copy_arc_side(vtx, vtxm1, on_first);
        }
    }
}

/// `SetArc(OnFirst, Src.ArcOnS*, Src.ParameterOnArc*, Src.TransitionLineArc*,
/// Src.TransitionOnS*)` (`IntPatch_ALine.cxx:557-562`).
fn copy_arc_side(dst: &mut PatchPoint, src: &PatchPoint, on_first: bool) {
    dst.set_arc(
        on_first,
        src.arc_id(on_first),
        src.param_on_arc(on_first),
        if on_first { src.trans1 } else { src.trans2 },
    );
}

/// `if (CopyVtx.IsOnDomS1()) { VTXm1.SetArc(true, ...) }` (`cxx:632-648`).
fn copy_arc_if_on_dom(dst: &mut PatchPoint, src: &PatchPoint, on_first: bool) {
    let has = if on_first {
        src.is_on_dom_s1()
    } else {
        src.is_on_dom_s2()
    };
    if has {
        copy_arc_side(dst, src, on_first);
    }
}
