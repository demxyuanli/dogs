//! `IntPatch_WLine::ComputeVertexParameters`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_WLine.cxx:323-1211`.
//!
//! The C++ body is one 890-line function; it is split here into the phases the
//! comments of the source delimit, each carrying the source line range. The
//! phase order, the loop bounds, the `break`-by-flag style (`!APointDeleted`,
//! `!SortIsOK`, `encoreunefois`) and every vertex index fix-up are kept
//! verbatim: `fipt`/`lapt`/`indf`/`indl` are the `WLine` members
//! `has_first_point` / `has_last_point` / `first_index` / `last_index`
//! (`IntPatch_WLine.hxx:205-208`).
//!
//! `myCreationWay == IntPatch_WLImpImp` branches on `WLineWay::ImpImp`.

use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};

use crate::int_tools_wline::{PatchPoint, PntOn2S, WLine, WLineWay};
use crate::intpatch_vparams::common::{
    compare_vertex_and_point, compare_vertices_on_s1, compare_vertices_on_s2,
    recadre_pnt_meme_periode, same_vtx_rst,
};

/// `IntPatch_WLine::ComputeVertexParameters(RTol)` (`IntPatch_WLine.cxx:323`).
pub fn compute_vertex_parameters(wl: &mut WLine, r_tol: f64) {
    let mut tol = r_tol;
    duplicate_points_on_different_rst(wl);
    tol = drop_same_rst_duplicates(wl, tol);
    sort_vertices(wl);
    drop_coincident_points(wl);
    map_vertices_on_points(wl);
    drop_same_parameter_pairs(wl);
    sort_and_drop_domainless_duplicates(wl);
    add_end_vertices(wl, tol);
    drop_intermediate_on_restriction(wl, true);
    drop_intermediate_on_restriction(wl, false);
    sort_vertices(wl);
}

/// `svtx.Remove(j)` followed by the `indl` / `indf` fix-up that every removal
/// site of `ComputeVertexParameters` performs.
fn remove_vertex(wl: &mut WLine, j: i32) {
    wl.vertices.remove((j - 1) as usize);
    if wl.has_last_point && wl.last_index >= j {
        wl.last_index -= 1;
    }
    if wl.has_first_point && wl.first_index >= j {
        wl.first_index -= 1;
    }
}

/// `svtx.Exchange(i - 1, i)`.
fn exchange_vertices(wl: &mut WLine, i: i32) {
    wl.vertices.swap((i - 1) as usize, i as usize);
}

/// Sort vertices by `ParameterOnLine` (`IntPatch_WLine.cxx:585-601`, `:1181-1195`).
fn sort_vertices(wl: &mut WLine) {
    let nbvtx = wl.nb_vertex();
    loop {
        let mut sort_is_ok = true;
        for i in 2..=nbvtx {
            if wl.vertex(i - 1).parameter_on_line() > wl.vertex(i).parameter_on_line() {
                sort_is_ok = false;
                exchange_vertices(wl, i);
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// Phase "Traitement des aretes de couture" (`IntPatch_WLine.cxx:358-487`):
/// duplicate the points that lie on different restrictions.
fn duplicate_points_on_different_rst(wl: &mut WLine) {
    loop {
        let nbvtx = wl.nb_vertex();
        let mut encoreunefois = false;
        let mut i = 1;
        while i <= nbvtx && !encoreunefois {
            let mut j = 1;
            while j <= nbvtx && !encoreunefois {
                if i != j
                    && wl.vertex(i).parameter_on_line() != wl.vertex(j).parameter_on_line()
                {
                    let vtx_i = wl.vertex(i).clone();
                    let vtx_j = wl.vertex(j).clone();
                    let d = vtx_i.p.distance(&vtx_j.p);
                    let toli = vtx_i.tolerance;
                    let tolj = vtx_j.tolerance;
                    let maxtol = if toli > tolj { toli } else { tolj };
                    if d < maxtol
                        || compare_vertices_on_s1(&vtx_i, &vtx_j)
                        || compare_vertices_on_s2(&vtx_i, &vtx_j)
                    {
                        let newtoli = if toli > tolj + d * 1.01 {
                            toli
                        } else {
                            tolj + d * 1.01
                        };
                        let newtolj = if tolj > toli + d * 1.01 {
                            tolj
                        } else {
                            toli + d * 1.01
                        };
                        for on_first in [true, false] {
                            let mut acreer = false;
                            if on_first {
                                if vtx_i.is_on_dom_s1() {
                                    if vtx_j.is_on_dom_s1() {
                                        if vtx_j.arc_id(true) != vtx_i.arc_id(true) {
                                            acreer = true;
                                        }
                                    } else {
                                        acreer = true;
                                    }
                                }
                            } else if vtx_i.is_on_dom_s2() {
                                if vtx_j.is_on_dom_s2() {
                                    if vtx_j.arc_id(false) != vtx_i.arc_id(false) {
                                        acreer = true;
                                    }
                                } else {
                                    acreer = true;
                                }
                            }
                            if !acreer {
                                continue;
                            }
                            let mut vtx = vtx_j.clone();
                            vtx.set_arc(
                                on_first,
                                vtx_i.arc_id(on_first),
                                vtx_i.param_on_arc(on_first),
                                if on_first { vtx_i.trans1 } else { vtx_i.trans2 },
                            );
                            let mut found = false;
                            for k in 1..=nbvtx {
                                if same_vtx_rst(wl.vertex(k), &vtx) {
                                    found = true;
                                    break;
                                }
                            }
                            if !found {
                                wl.vertex_mut(i).set_tolerance(newtoli);
                                wl.vertex_mut(j).set_tolerance(newtolj);
                                vtx.set_tolerance(newtolj);
                                wl.vertices.push(vtx);
                                encoreunefois = true;
                            }
                        }
                    }
                }
                j += 1;
            }
            i += 1;
        }
        if !encoreunefois {
            break;
        }
    }
}

/// Phases of `IntPatch_WLine.cxx:490-581`: two vertices on the same restriction
/// and only on that one must not share the same parameter.
fn drop_same_rst_duplicates(wl: &mut WLine, mut tol: f64) -> f64 {
    for on_first in [true, false] {
        loop {
            let nbvtx = wl.nb_vertex();
            let mut a_point_deleted = false;
            let mut i = 1;
            while i <= nbvtx && !a_point_deleted {
                let vtx_i = wl.vertex(i).clone();
                if vtx_i.tolerance > tol {
                    tol = vtx_i.tolerance;
                }
                if vtx_i.is_on_dom_s1() != on_first || vtx_i.is_on_dom_s2() == on_first {
                    i += 1;
                    continue;
                }
                let mut j = 1;
                while j <= nbvtx && !a_point_deleted {
                    if i != j {
                        let vtx_j = wl.vertex(j).clone();
                        if vtx_j.is_on_dom_s1() == on_first
                            && vtx_j.is_on_dom_s2() != on_first
                            && vtx_i.parameter_on_line() == vtx_j.parameter_on_line()
                            && vtx_i.arc_id(on_first) == vtx_j.arc_id(on_first)
                        {
                            remove_vertex(wl, j);
                            a_point_deleted = true;
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
    tol
}

/// `IntPatch_WLine.cxx:603-631`: drop the coincident points of the walking line
/// and shift the vertex parameters that sit after them.
fn drop_coincident_points(wl: &mut WLine) {
    let mut nbponline = wl.nb_pnts();
    let mut i = 2;
    while i <= nbponline && nbponline > 2 {
        let d = wl
            .point(i - 1)
            .value()
            .square_distance(&wl.point(i).value());
        if d < SQUARE_CONFUSION {
            wl.remove_point(i);
            nbponline -= 1;
            let nbvtx = wl.nb_vertex();
            for j in 1..=nbvtx {
                let indicevertex = wl.vertex(j).parameter_on_line();
                if indicevertex >= i as f64 {
                    wl.vertex_mut(j).set_parameter(indicevertex - 1.0);
                }
            }
            i -= 1;
        }
        i += 1;
    }
}

/// `IntPatch_WLine.cxx:632-832`: move every vertex onto the walking-line point
/// that carries it, inserting or substituting points as needed.
fn map_vertices_on_points(wl: &mut WLine) {
    let nbvtx = wl.nb_vertex();
    let mut nbponline = wl.nb_pnts();
    for i in 1..=nbvtx {
        let pp = wl.vertex(i).clone();
        let p = pp.p;
        let v_tol = pp.tolerance;
        let indicevertex = pp.parameter_on_line();
        let indicevertexonline = indicevertex as i32;
        let u1p = wl.u1_period;
        let v1p = wl.v1_period;
        let u2p = wl.u2_period;
        let v2p = wl.v2_period;
        if indicevertexonline < 1 {
            if compare_vertex_and_point(&p, &wl.point(1).value(), v_tol) {
                let mut pon2s: PntOn2S = wl.vertex(i).pnt_on_2s();
                recadre_pnt_meme_periode(&mut pon2s, &wl.point(1).clone(), u1p, v1p, u2p, v2p);
                if wl.creating_way == WLineWay::ImpImp {
                    wl.point_mut(1).set_point(pon2s.value());
                    let (mu1, mv1, mu2, mv2) = wl.point(1).parameters();
                    wl.vertex_mut(i).set_parameter(1.0);
                    wl.vertex_mut(i).set_parameters(mu1, mv1, mu2, mv2);
                } else {
                    *wl.point_mut(1) = pon2s;
                    wl.vertex_mut(i).set_parameter(1.0);
                }
            } else {
                let mut pon2s: PntOn2S = wl.vertex(i).pnt_on_2s();
                recadre_pnt_meme_periode(&mut pon2s, &wl.point(1).clone(), u1p, v1p, u2p, v2p);
                wl.insert_before(1, pon2s);
                wl.vertex_mut(i).set_parameter(1.0);
                nbponline += 1;
                for j in 1..=nbvtx {
                    if j != 1 {
                        let t = wl.vertex(j).parameter_on_line();
                        if t > 1.0 {
                            wl.vertex_mut(j).set_parameter(t + 1.0);
                        }
                    }
                }
            }
        } else {
            let (mut mu1, mut mv1, mut mu2, mut mv2) = wl.point(indicevertexonline).parameters();
            let mut substitution = false;
            let mut k = indicevertexonline + 1;
            while k >= indicevertexonline - 1 {
                if k > 0 && k <= nbponline && compare_vertex_and_point(&p, &wl.point(k).value(), v_tol)
                {
                    let mut pon2s: PntOn2S = wl.vertex(i).pnt_on_2s();
                    recadre_pnt_meme_periode(
                        &mut pon2s,
                        &wl.point(k).clone(),
                        u1p,
                        v1p,
                        u2p,
                        v2p,
                    );
                    if wl.creating_way == WLineWay::ImpImp {
                        wl.point_mut(k).set_point(pon2s.value());
                        wl.point_mut(k).set_uv_on(true, mu1, mv1);
                        wl.point_mut(k).set_uv_on(false, mu2, mv2);
                    } else {
                        *wl.point_mut(k) = pon2s;
                        (mu1, mv1, mu2, mv2) = pon2s.parameters();
                    }
                    wl.vertex_mut(i).set_parameter(k as f64);
                    wl.vertex_mut(i).set_parameters(mu1, mv1, mu2, mv2);
                    substitution = true;
                }
                k -= 1;
            }
            if substitution {
                let mut ind_point = 2;
                while ind_point <= nbponline && nbponline > 1 {
                    let d = wl
                        .point(ind_point - 1)
                        .value()
                        .square_distance(&wl.point(ind_point).value());
                    if d < SQUARE_CONFUSION {
                        wl.remove_point(ind_point);
                        nbponline -= 1;
                        for j in 1..=nbvtx {
                            let indicevertex = wl.vertex(j).parameter_on_line();
                            if indicevertex >= ind_point as f64 {
                                wl.vertex_mut(j).set_parameter(indicevertex - 1.0);
                            }
                        }
                        ind_point -= 1;
                    }
                    ind_point += 1;
                }
            } else {
                let mut pon2s: PntOn2S = wl.vertex(i).pnt_on_2s();
                if indicevertexonline >= nbponline {
                    recadre_pnt_meme_periode(
                        &mut pon2s,
                        &wl.point(nbponline).clone(),
                        u1p,
                        v1p,
                        u2p,
                        v2p,
                    );
                    wl.add(pon2s);
                } else {
                    recadre_pnt_meme_periode(
                        &mut pon2s,
                        &wl.point(indicevertexonline + 1).clone(),
                        u1p,
                        v1p,
                        u2p,
                        v2p,
                    );
                    wl.insert_before(indicevertexonline + 1, pon2s);
                }
                wl.vertex_mut(i).set_parameter(indicevertexonline as f64 + 1.0);
                nbponline += 1;
                for j in 1..=nbvtx {
                    if j != i {
                        let t = wl.vertex(j).parameter_on_line();
                        if t > indicevertexonline as f64 {
                            wl.vertex_mut(j).set_parameter(t + 1.0);
                        }
                    }
                }
            }
        }
    }
}

/// `IntPatch_WLine.cxx:834-923`: two vertices with the same parameter on the
/// same restriction(s) collapse into one.
fn drop_same_parameter_pairs(wl: &mut WLine) {
    loop {
        let mut a_point_deleted = false;
        let mut i = 1;
        let nbvtx = wl.nb_vertex();
        while i <= nbvtx && !a_point_deleted {
            let mut j = 1;
            while j <= nbvtx && !a_point_deleted {
                if i != j {
                    let vtx = wl.vertex(i).clone();
                    let vtxm1 = wl.vertex(j).clone();
                    let mut kill = false;
                    let mut killm1 = false;
                    if vtxm1.parameter_on_line() == vtx.parameter_on_line() {
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
                        }
                        if !(kill || killm1)
                            && vtxm1.is_on_dom_s2()
                            && vtx.is_on_dom_s2()
                            && vtxm1.arc_id(false) == vtx.arc_id(false)
                        {
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
                        if kill {
                            a_point_deleted = true;
                            remove_vertex(wl, i);
                        } else if killm1 {
                            a_point_deleted = true;
                            remove_vertex(wl, j);
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

/// `IntPatch_WLine.cxx:925-957`: sort, then drop the vertex of a same-parameter
/// pair that lies on no restriction.
fn sort_and_drop_domainless_duplicates(wl: &mut WLine) {
    loop {
        let mut sort_is_ok = true;
        let nbvtx = wl.nb_vertex();
        for i in 2..=nbvtx {
            if !sort_is_ok {
                break;
            }
            let pim1 = wl.vertex(i - 1).clone();
            let pii = wl.vertex(i).clone();
            if pim1.parameter_on_line() == pii.parameter_on_line() {
                if !pii.is_on_dom_s1() && !pii.is_on_dom_s2() {
                    sort_is_ok = false;
                    wl.vertices.remove((i - 1) as usize);
                } else if !pim1.is_on_dom_s1() && !pim1.is_on_dom_s2() {
                    sort_is_ok = false;
                    wl.vertices.remove((i - 2) as usize);
                }
            }
        }
        if sort_is_ok {
            break;
        }
    }
}

/// `IntPatch_WLine.cxx:959-1036`: make sure a vertex sits on the first and on
/// the last point of the walking line.
fn add_end_vertices(wl: &mut WLine, tol: f64) {
    let mut nbponline = wl.nb_pnts();
    let nbvtx = wl.nb_vertex();
    let mut b_first = false;
    let mut b_last = false;
    for i in 1..=nbvtx {
        let pol = wl.vertex(i).parameter_on_line();
        if pol == 1.0 {
            b_first = true;
            wl.has_first_point = true;
            wl.first_index = i;
        }
        if pol == nbponline as f64 {
            b_last = true;
            wl.has_last_point = true;
            wl.last_index = i;
        }
    }
    if !b_first {
        let (pu1, pv1, pu2, pv2) = wl.point(1).parameters();
        let mut vtxfound = false;
        let mut vtx = PatchPoint::new(wl.point(1).value(), 1.0, pu1, pv1, pu2, pv2);
        for i in 1..=wl.nb_vertex() {
            let v = wl.vertex(i).clone();
            if compare_vertex_and_point(&v.p, &wl.point(1).value(), CONFUSION) {
                vtx = v;
                vtx.set_parameters(pu1, pv1, pu2, pv2);
                vtxfound = true;
                break;
            }
        }
        if !vtxfound {
            vtx = PatchPoint::new(wl.point(1).value(), 1.0, pu1, pv1, pu2, pv2);
            vtx.set_tolerance(tol);
        }
        vtx.set_parameter(1.0);
        wl.vertices.insert(0, vtx);
        wl.has_first_point = true;
        wl.first_index = 1;
    }
    if !b_last {
        nbponline = wl.nb_pnts();
        let (pu1, pv1, pu2, pv2) = wl.point(nbponline).parameters();
        let mut vtxfound = false;
        let mut vtx = PatchPoint::new(wl.point(nbponline).value(), 1.0, pu1, pv1, pu2, pv2);
        for i in 1..=wl.nb_vertex() {
            let v = wl.vertex(i).clone();
            if compare_vertex_and_point(&v.p, &wl.point(nbponline).value(), CONFUSION) {
                vtx = v;
                vtx.set_parameters(pu1, pv1, pu2, pv2);
                vtxfound = true;
                break;
            }
        }
        if !vtxfound {
            vtx = PatchPoint::new(wl.point(nbponline).value(), 1.0, pu1, pv1, pu2, pv2);
            vtx.set_tolerance(tol);
        }
        vtx.set_parameter(nbponline as f64);
        wl.vertices.push(vtx);
        wl.has_last_point = true;
        wl.last_index = wl.nb_vertex();
    }
}

/// `IntPatch_WLine.cxx:1043-1107` (S1) and `:1112-1161` (S2): drop the vertex
/// of a pair of neighbours on the same restriction with the same transition.
fn drop_intermediate_on_restriction(wl: &mut WLine, on_first: bool) {
    loop {
        let nbvtx = wl.nb_vertex();
        let nb_pnts = wl.nb_pnts();
        let mut a_point_deleted = false;
        let mut i = 1;
        while i <= nbvtx && !a_point_deleted {
            let vtx_i = wl.vertex(i).clone();
            if vtx_i.is_on_dom_s1() == on_first && vtx_i.is_on_dom_s2() != on_first {
                let mut j = 1;
                while j <= nbvtx && !a_point_deleted {
                    if i != j {
                        let vtx_j = wl.vertex(j).clone();
                        if vtx_j.is_on_dom_s1() == on_first
                            && vtx_j.is_on_dom_s2() != on_first
                            && (vtx_i.parameter_on_line() == vtx_j.parameter_on_line() + 1.0
                                || vtx_i.parameter_on_line() == vtx_j.parameter_on_line() - 1.0)
                            && vtx_i.arc_id(on_first) == vtx_j.arc_id(on_first)
                        {
                            let t1 = if on_first { vtx_i.trans1 } else { vtx_i.trans2 };
                            let t2 = if on_first { vtx_j.trans1 } else { vtx_j.trans2 };
                            if t1 == t2
                                && vtx_j.parameter_on_line() != 1.0
                                && vtx_j.parameter_on_line() != nb_pnts as f64
                            {
                                remove_vertex(wl, j);
                                a_point_deleted = true;
                            }
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
