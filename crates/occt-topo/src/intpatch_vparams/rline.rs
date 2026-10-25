//! `IntPatch_RLine::ComputeVertexParameters`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_RLine.cxx:143-433`.
//!
//! The restriction line keeps vertices only; unlike the walking line the source
//! never touches its `curv` here, so no point insertion / removal happens.
//! `restrdiff` (the source flag that cancels the kill when the two vertices lie
//! on different restrictions) is kept as a local.

use crate::geom_int::RLine;

/// `IntPatch_RLine::ComputeVertexParameters` (`IntPatch_RLine.cxx:144`).
pub fn compute_vertex_parameters_rline(rline: &mut RLine) {
    drop_same_rst_duplicates(rline);
    sort_vertices(rline);
    drop_same_parameter_pairs(rline);
    drop_domainless_duplicates(rline);
    sort_vertices(rline);
    let nbvtx = rline.vertices.len() as i32;
    rline.has_first_point = true;
    rline.first_index = 1;
    rline.has_last_point = true;
    rline.last_index = nbvtx;
}

/// `svtx.Remove(j)` of `IntPatch_RLine.cxx:177-195`, `:225-243`.
fn remove_vertex(rline: &mut RLine, j: i32) {
    rline.vertices.remove((j - 1) as usize);
    if rline.has_last_point && rline.last_index >= j {
        rline.last_index -= 1;
    }
    if rline.has_first_point && rline.first_index >= j {
        rline.first_index -= 1;
    }
}

/// `IntPatch_RLine.cxx:253-268`.
fn sort_vertices(rline: &mut RLine) {
    let nbvtx = rline.vertices.len() as i32;
    loop {
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            if rline.vertices[(i - 2) as usize].parameter_on_line()
                > rline.vertices[(i - 1) as usize].parameter_on_line()
            {
                sort_is_ok = false;
                rline.vertices.swap((i - 2) as usize, (i - 1) as usize);
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// `IntPatch_RLine.cxx:155-251`: two vertices on the same restriction and only
/// on that one must not share the same parameter and arc parameter.
fn drop_same_rst_duplicates(rline: &mut RLine) {
    for on_first in [true, false] {
        loop {
            let nbvtx = rline.vertices.len() as i32;
            let mut a_point_deleted = false;
            let mut i = 1;
            while i <= nbvtx && !a_point_deleted {
                let vtx_i = rline.vertices[(i - 1) as usize].clone();
                if vtx_i.is_on_dom_s1() == on_first && vtx_i.is_on_dom_s2() != on_first {
                    let mut j = 1;
                    while j <= nbvtx && !a_point_deleted {
                        if i != j {
                            let vtx_j = rline.vertices[(j - 1) as usize].clone();
                            if vtx_j.is_on_dom_s1() == on_first
                                && vtx_j.is_on_dom_s2() != on_first
                                && vtx_i.parameter_on_line() == vtx_j.parameter_on_line()
                                && vtx_i.arc_id(on_first) == vtx_j.arc_id(on_first)
                                && vtx_i.param_on_arc(on_first) == vtx_j.param_on_arc(on_first)
                            {
                                remove_vertex(rline, j);
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
}

/// `IntPatch_RLine.cxx:270-386`.
fn drop_same_parameter_pairs(rline: &mut RLine) {
    loop {
        let nbvtx = rline.vertices.len() as i32;
        let mut a_point_deleted = false;
        let mut i = 1;
        while i <= nbvtx && !a_point_deleted {
            let mut j = 1;
            while j <= nbvtx && !a_point_deleted {
                if i != j {
                    let vtx = rline.vertices[(i - 1) as usize].clone();
                    let vtxm1 = rline.vertices[(j - 1) as usize].clone();
                    let mut kill = false;
                    let mut killm1 = false;
                    if vtxm1.parameter_on_line() == vtx.parameter_on_line() {
                        let mut restrdiff = false;
                        if vtxm1.is_on_dom_s1() && vtx.is_on_dom_s1() {
                            if vtxm1.arc_id(true) == vtx.arc_id(true) {
                                if vtx.param_on_arc(true) == vtxm1.param_on_arc(true) {
                                    if vtxm1.is_on_dom_s2() {
                                        if !vtx.is_on_dom_s2() {
                                            kill = true;
                                        } else if vtxm1.arc_id(false) == vtx.arc_id(false)
                                            && vtx.param_on_arc(false)
                                                == vtxm1.param_on_arc(false)
                                        {
                                            kill = true;
                                        }
                                    } else if vtx.is_on_dom_s2() {
                                        killm1 = true;
                                    }
                                }
                            } else {
                                restrdiff = true;
                            }
                        }
                        if !restrdiff && !(kill || killm1) {
                            if vtxm1.is_on_dom_s2() && vtx.is_on_dom_s2() {
                                if vtxm1.arc_id(false) == vtx.arc_id(false) {
                                    if vtx.param_on_arc(false) == vtxm1.param_on_arc(false) {
                                        if vtxm1.is_on_dom_s1() {
                                            if !vtx.is_on_dom_s1() {
                                                kill = true;
                                            } else if vtxm1.arc_id(true) == vtx.arc_id(true)
                                                && vtx.param_on_arc(true)
                                                    == vtxm1.param_on_arc(true)
                                            {
                                                kill = true;
                                            }
                                        } else if vtx.is_on_dom_s1() {
                                            killm1 = true;
                                        }
                                    }
                                } else {
                                    restrdiff = true;
                                }
                            }
                        }
                        if !restrdiff {
                            if kill {
                                a_point_deleted = true;
                                remove_vertex(rline, i);
                            } else if killm1 {
                                a_point_deleted = true;
                                remove_vertex(rline, j);
                            }
                        }
                    }
                }
                j += 1;
            }
            i += 1;
        }
        if !a_point_deleted {
            break;
        }
    }
}

/// `IntPatch_RLine.cxx:388-415`: drop the vertex of a same-parameter pair that
/// lies on no restriction.
fn drop_domainless_duplicates(rline: &mut RLine) {
    loop {
        let nbvtx = rline.vertices.len() as i32;
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            if !sort_is_ok {
                break;
            }
            let pim1 = rline.vertices[(i - 2) as usize].clone();
            let pii = rline.vertices[(i - 1) as usize].clone();
            if pim1.parameter_on_line() == pii.parameter_on_line() {
                if !pii.is_on_dom_s1() && !pii.is_on_dom_s2() {
                    sort_is_ok = false;
                    rline.vertices.remove((i - 1) as usize);
                } else if !pim1.is_on_dom_s1() && !pim1.is_on_dom_s2() {
                    sort_is_ok = false;
                    rline.vertices.remove((i - 2) as usize);
                }
            }
        }
        if sort_is_ok {
            break;
        }
    }
}
