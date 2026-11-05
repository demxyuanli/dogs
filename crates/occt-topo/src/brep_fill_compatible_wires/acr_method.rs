//! `BRepFill_CompatibleWires::SameNumberByACR` (`cxx:1505-1785`).
//!
//! The open-wire branch when `report` is set or the sections disagree on their
//! edge count: every section is cut at the union of the reduced curvilinear
//! abscissae of all sections (`BRepFill::ComputeACR`) and re-cut with
//! `BRepFill::InsertACR`, then `myMap` is grown so every original edge maps to
//! the pieces that replaced it.
//!
//! Deviation, an out-of-bounds read / write in OCCT: the sorted-cut insertion
//! (`cxx:1580-1596`) reads `dec(k)` and, when it inserts, writes `dec(k)` for
//! `k = nbdec + 1` after the `k > nbdec` break. `NCollection_Array1` does not
//! bounds-check in a release build, so OCCT touches memory past the array. This
//! port drops such an insertion instead of touching memory out of bounds.

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, RESOLUTION};

use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::{edge_vertices_cum_ori, WireExplorer};
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::is_same;

/// `SameNumberByACR(report)` (`cxx:1505-1785`).
pub(crate) fn same_number_by_acr(
    cw: &mut CompatibleWires,
    report: bool,
) -> Result<(), String> {
    // `cxx:1507-1520`.
    let mut ideb = 1usize;
    let mut ifin = cw.my_work.len();
    if cw.my_degen1 {
        ideb += 1;
    }
    if cw.my_degen2 {
        ifin -= 1;
    }
    let v_closed = (!cw.my_degen1)
        && (!cw.my_degen2)
        && ideb <= ifin
        && is_same(&cw.my_work[ideb - 1], &cw.my_work[ifin - 1]);

    // `cxx:1522-1546`: the edge counts over every section.
    let nb_sects = cw.my_work.len();
    let mut nb_edges = vec![0usize; nb_sects + 1];
    let mut nbmax = 0usize;
    let mut nbmin = 0usize;
    for i in 1..=nb_sects {
        nb_edges[i] = wire_edge_count(&Wire(cw.my_work[i - 1].clone()));
        if i == 1 {
            nbmin = nb_edges[i];
        }
        if nbmax < nb_edges[i] {
            nbmax = nb_edges[i];
        }
        if nbmin > nb_edges[i] {
            nbmin = nb_edges[i];
        }
    }

    if nbmax > 1 && (report || nbmin < nbmax) {
        // `cxx:1552-1557`: `dec` is the sorted cut table, 1-based, with the
        // sentinel `1` at index 2.
        let mut nbdec = (nbmax - 1) * nb_sects + 1;
        let mut dec = vec![0.0f64; nbdec + 1];
        dec[2] = 1.0;

        // `cxx:1559-1597`: merge the ACR of every section into `dec`.
        let mut wire_len = vec![0.0f64; nb_sects + 1];
        for i in 1..=nb_sects {
            let wire1 = Wire(cw.my_work[i - 1].clone());
            let acr = crate::brep_fill::compute_acr(&wire1)?;
            wire_len[i] = acr[0];
            // `for (j = 1; j < ACR.Length() - 1; j++)`: every intermediate
            // abscissa, not the total and not the final `1`.
            for j in 1..(acr.len() - 1) {
                let value = acr[j];
                let mut k = 1usize;
                while k <= nbdec && dec[k] < value {
                    k += 1;
                }
                if k > nbdec {
                    // OCCT reads and writes past the array here; see the module
                    // note.
                    continue;
                }
                if dec[k - 1] < value && value < dec[k] {
                    for l in (k..nbdec).rev() {
                        dec[l + 1] = dec[l];
                    }
                    dec[k] = value;
                }
            }
        }

        // `cxx:1599-1615`: drop everything from the first `1` on.
        let mut k = 1usize;
        while k <= nbdec && dec[k] < 1.0 {
            k += 1;
        }
        nbdec = k - 1;
        let dec2: Vec<f64> = dec[1..=nbdec].to_vec();

        // `cxx:1617-1637`: drop the cuts whose every section piece is too short.
        let mut cuts_to_remove: Vec<usize> = Vec::new();
        for k in 1..=nbdec {
            let knot1 = dec2[k - 1];
            let knot2 = if k == nbdec { 1.0 } else { dec2[k] };
            let mut all_lengths_null = true;
            for i in 1..=nb_sects {
                let edge_len = (knot2 - knot1) * wire_len[i];
                if edge_len > CONFUSION {
                    all_lengths_null = false;
                    break;
                }
            }
            if all_lengths_null {
                cuts_to_remove.push(k);
            }
        }
        // `cxx:1638-1648`.
        let mut dec3: Vec<f64> = Vec::new();
        for k in 1..=nbdec {
            if !cuts_to_remove.contains(&k) {
                dec3.push(dec2[k - 1]);
            }
        }

        // `cxx:1650-1747`: insert the cuts into every section.
        for i in 1..=nb_sects {
            let old_wire = Wire(cw.my_work[i - 1].clone());
            let mut tol = CONFUSION;
            if wire_len[i] > RESOLUTION {
                tol /= wire_len[i];
            }
            let new_wire = crate::brep_fill::insert_acr(&old_wire, &dec3, tol)?;
            update_map_for_split(cw, &old_wire, &new_wire);
            cw.my_work[i - 1] = new_wire.0;
        }
    }

    // `cxx:1750-1753`.
    if v_closed {
        let last = cw.my_work.len();
        cw.my_work[last - 1] = cw.my_work[0].clone();
    }

    // `cxx:1755-1784`: the debug pass that also rejects an unequal cut count.
    nbmax = 0;
    nbmin = 0;
    for i in ideb..=ifin {
        nb_edges[i] = wire_edge_count(&Wire(cw.my_work[i - 1].clone()));
        if i == ideb {
            nbmin = nb_edges[i];
        }
        if nbmax < nb_edges[i] {
            nbmax = nb_edges[i];
        }
        if nbmin > nb_edges[i] {
            nbmin = nb_edges[i];
        }
    }
    if nbmax != nbmin {
        cw.my_status = ThruSectionErrorStatus::Failed;
        return Ok(());
    }

    Ok(())
}

/// `cxx:1659-1746`: walk `oldwire` and `newwire` together; wherever the new
/// wire replaced an old edge by several pieces, splice those pieces into
/// `myMap` in place of the old edge.
fn update_map_for_split(cw: &mut CompatibleWires, old_wire: &Wire, new_wire: &Wire) {
    let old_seq = wire_pairs(old_wire);
    let new_seq = wire_pairs(new_wire);
    let mut i2 = 0usize;
    for (v1, e_old) in &old_seq {
        // `anExp2.Current()` on an exhausted explorer is invalid in OCCT.
        if i2 >= new_seq.len() {
            break;
        }
        if is_same(&e_old.0, &new_seq[i2].1.0) {
            i2 += 1;
            continue;
        }

        // `cxx:1667-1711`: collect the pieces up to the point that lands on the
        // old edge's far vertex.
        let mut le: Vec<TopoShape> = Vec::new();
        let p1 = other_vertex_point(v1.as_ref(), e_old);
        let mut p2 = other_vertex_point(new_seq[i2].0.as_ref(), &new_seq[i2].1);
        while p1.distance(&p2) > 1.0e-3 {
            le.push(new_seq[i2].1.0.clone());
            i2 += 1;
            if i2 >= new_seq.len() {
                break;
            }
            p2 = other_vertex_point(new_seq[i2].0.as_ref(), &new_seq[i2].1);
            if p1.distance(&p2) <= 1.0e-3 {
                le.push(new_seq[i2].1.0.clone());
                i2 += 1;
                if i2 >= new_seq.len() {
                    break;
                }
            }
        }

        // `cxx:1713-1737`: the ancestor is the `myMap` key whose list holds the
        // old edge; the pieces go right before it and the old edge is removed.
        for list in cw.my_map.values_mut() {
            let Some(pos) = list.iter().position(|s| is_same(s, &e_old.0)) else {
                continue;
            };
            list.splice(pos..=pos, le.iter().cloned());
            break;
        }
    }
}

/// `TopExp::Vertices(E, VF, VR, true)` followed by the `V1/V2` matching of
/// `cxx:1671-1693`: the endpoint of `e` opposite `v`, or the origin when the
/// vertex is neither endpoint (OCCT leaves `P1` / `P2` at `gp_Pnt()`).
fn other_vertex_point(v: Option<&Vertex>, e: &Edge) -> GpPnt {
    let (vf, vr) = edge_vertices_cum_ori(e);
    let mut p = GpPnt::zero();
    if let Some(v) = v {
        if let Some(vf) = &vf {
            if is_same(&v.0, &vf.0) {
                if let Some(vr) = &vr {
                    p = BRepTool::vertex_point_world(vr);
                }
            }
        }
        if let Some(vr) = &vr {
            if is_same(&v.0, &vr.0) {
                if let Some(vf) = &vf {
                    p = BRepTool::vertex_point_world(vf);
                }
            }
        }
    }
    p
}

/// `anExp.CurrentVertex()` / `anExp.Current()` in step (`cxx:1660-1663`).
fn wire_pairs(w: &Wire) -> Vec<(Option<Vertex>, Edge)> {
    let mut out = Vec::new();
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    while exp.more() {
        let v = exp.current_vertex().cloned();
        if let Some(e) = exp.current_edge() {
            out.push((v, e));
        }
        exp.next();
    }
    out
}

/// `for (anExp.Init(w); anExp.More(); anExp.Next()) nbEdges++`.
fn wire_edge_count(w: &Wire) -> usize {
    let mut exp = WireExplorer::new();
    exp.init(w, None);
    let mut n = 0usize;
    while exp.more() {
        n += 1;
        exp.next();
    }
    n
}
