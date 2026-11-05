//! `BRepFill_CompatibleWires::Perform`.
//!
//! Source: `BRepFill_CompatibleWires.cxx:855-995`.
//!
//! The analysis part: the FORWARD normalisation of every section
//! (`cxx:870-874`), `WireContinuity` per wire (`cxx:876-881`), the edge counts
//! and `nbmax` / `nbmin` (`cxx:883-898`), the `report` decision (`cxx:902`),
//! the `myMap` initialisation (`cxx:905-914`), `myDegen1` / `myDegen2`
//! (`cxx:920-940`), the `wClosed` check that also accepts a wire whose first and
//! last vertices are the same (`cxx:942-958`) and the closed / open dispatch
//! (`cxx:960-994`).
//!
//! The dispatch targets live in their own modules:
//! [`same_number_by_polar_method`] (`cxx:1008-1503`), [`compute_origin`]
//! (`cxx:1787-2126`), [`search_origin`] (`cxx:2453-2629`) and
//! [`same_number_by_acr`] (`cxx:1505-1785`, with `BRepFill::InsertACR`).

use crate::abs::Orientation;
use crate::brep_fill_compatible_wires::acr_method::same_number_by_acr;
use crate::brep_fill_compatible_wires::compute_origin::compute_origin;
use crate::brep_fill_compatible_wires::continuity::wire_continuity;
use crate::brep_fill_compatible_wires::polar_method::same_number_by_polar_method;
use crate::brep_fill_compatible_wires::search_origin::search_origin;
use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::Wire;
use crate::topo_tools::{is_same, shape_map_key, vertices_of};

/// `Perform(WithRotation)` (`cxx:855-995`).
///
/// `WithRotation` is only consumed by `SameNumberByPolarMethod` (`cxx:967`) and
/// `SearchOrigin`'s runner-up arms; it is forwarded to both.
pub(crate) fn perform(cw: &mut CompatibleWires, with_rotation: bool) -> Result<(), String> {
    cw.my_status = ThruSectionErrorStatus::Done;

    let nb_sects = cw.my_work.len();
    // `NCollection_Array1<int> nbEdges(1, nbSects)` (`cxx:867`).
    let mut nb_edges = vec![0usize; nb_sects];
    let mut nbmax = 0usize;
    let mut nbmin = 0usize;
    // `GeomAbs_Shape contS = GeomAbs_CN` (`cxx:871`).
    let mut cont_s: u8 = crate::brep_lprop::CN;

    // `cxx:873-899`.
    for i in 1..=nb_sects {
        // `myWork(i) = TopoDS::Wire(myWork(i).Oriented(TopAbs_FORWARD))`
        // (`cxx:872-874`).
        let mut s = cw.my_work[i - 1].clone();
        s.set_orientation(Orientation::Forward);
        cw.my_work[i - 1] = s;
        let w = Wire(cw.my_work[i - 1].clone());

        // `WireContinuity(W, cont); if (cont < contS) contS = cont;`
        let cont = wire_continuity(&w);
        if cont < cont_s {
            cont_s = cont;
        }

        // `for (anExp.Init(W); anExp.More(); anExp.Next()) nbEdges(i)++`.
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        let mut count = 0usize;
        while exp.more() {
            count += 1;
            exp.next();
        }
        nb_edges[i - 1] = count;
        if i == 1 {
            nbmin = nb_edges[i - 1];
        }
        if nbmax < nb_edges[i - 1] {
            nbmax = nb_edges[i - 1];
        }
        if nbmin > nb_edges[i - 1] {
            nbmin = nb_edges[i - 1];
        }
    }

    // `cxx:902`: a vertex-to-vertex report is only used when the sections all
    // have the same edge count and none of them is at least C1.
    let report = nbmax != nbmin || cont_s >= crate::brep_lprop::C1;

    // `cxx:905-914`: `myMap.Bind(E, Empty); myMap(E).Append(E)`.
    for i in 1..=nb_sects {
        let w = Wire(cw.my_work[i - 1].clone());
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        while exp.more() {
            if let Some(e) = exp.current_edge() {
                cw.my_map.insert(shape_map_key(&e.0), vec![e.0.clone()]);
            }
            exp.next();
        }
    }

    // `cxx:920-931`: is the first section punctual?
    let mut ideb = 1usize;
    let mut ifin = nb_sects;
    cw.my_degen1 = true;
    if let Some(first) = cw.my_work.first() {
        let w = Wire(first.clone());
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        while exp.more() {
            if let Some(e) = exp.current_edge() {
                cw.my_degen1 = cw.my_degen1 && BRepTool::is_degenerated(&e);
            }
            exp.next();
        }
    }
    if cw.my_degen1 {
        ideb += 1;
    }
    // `cxx:932-940`: is the last section punctual?
    cw.my_degen2 = true;
    if let Some(last) = cw.my_work.last() {
        let w = Wire(last.clone());
        let mut exp = WireExplorer::new();
        exp.init(&w, None);
        while exp.more() {
            if let Some(e) = exp.current_edge() {
                cw.my_degen2 = cw.my_degen2 && BRepTool::is_degenerated(&e);
            }
            exp.next();
        }
    }
    if cw.my_degen2 && ifin > 0 {
        ifin -= 1;
    }

    // `cxx:942-958`.
    let mut all_closed = true;
    let mut all_open = true;
    if ideb <= ifin {
        for i in ideb..=ifin {
            let mut w_closed = cw.my_work[i - 1].is_closed();
            if !w_closed {
                // `TopExp::Vertices(TopoDS::Wire(myWork(i)), V1, V2)` followed
                // by `V1.IsSame(V2)`.
                let verts = vertices_of(&cw.my_work[i - 1]);
                if let (Some(v1), Some(v2)) = (verts.first(), verts.last()) {
                    if is_same(&v1.0, &v2.0) {
                        w_closed = true;
                    }
                }
            }
            all_closed = all_closed && w_closed;
            all_open = all_open && !w_closed;
        }
    }

    // `cxx:960-994`.
    if all_closed {
        if report {
            // `SameNumberByPolarMethod(WithRotation)` (`cxx:967`).
            return same_number_by_polar_method(cw, with_rotation);
        }
        // `ComputeOrigin(false)` (`cxx:972`).
        return compute_origin(cw, false);
    } else if all_open {
        // `SearchOrigin()` (`cxx:977`) then a `myStatus != Done` check.
        search_origin(cw)?;
        if cw.my_status != ThruSectionErrorStatus::Done {
            return Ok(());
        }
        // `cxx:982-985`: `SameNumberByACR(report)`.
        if report {
            return same_number_by_acr(cw, report);
        }
        return Ok(());
    }

    // `cxx:990-994`: open and closed sections mixed.
    cw.my_status = ThruSectionErrorStatus::NotSameTopology;
    Ok(())
}
