use occt_core::gp::GpTrsf2d;
use occt_core::precision::CONFUSION;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::shhealing::{first_vertex, fix_same_parameter, last_vertex, ShapeFixEdge};
use crate::tgeometry::GeometryRegistry;

use super::super::iges_boundary::{closure_failed, connect_next_wire};
use super::super::iges_to_brep_basic_curve::eps_geom;
use super::super::iges_to_brep_curve_and_surface::Session;

use super::curve_on_face::component_inconsistent;
use super::pcurve_transfer::{edges_from_2d, edges_of, transfer_pcurve};
use super::transfer_topo_curve;

pub(in crate::iges_read) fn boundary_max_tol(session: &Session<'_>) -> f64 {
    // `UpdateMinMaxTol` (`IGESToBRep_CurveAndSurface.cxx:133-138`).
    // `read.maxprecision.val` default is 1 (`DEIGES_Parameters.hxx:99`).
    let eps = eps_geom(session.model.resolution) * session.model.unit_factor;
    eps.max(1.0)
}

pub(in crate::iges_read) struct BoundaryPiece {
    pub out: Vec<Edge>,
    pub wire3: Vec<Edge>,
    pub wire2: Vec<Edge>,
    pub ok: bool,
    pub ok3: bool,
    pub ok2: bool,
}

pub(super) fn edges_when_both(
    session: &mut Session<'_>,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
    seq3d: &[usize],
    seq2d: &[usize],
    pref3d: bool,
    pref2d: bool,
    check_closure: bool,
) -> Result<Vec<Edge>, String> {
    let piece = boundary_both(
        session, face, trans, u_fact, seq3d, seq2d, pref3d, pref2d, false,
    )?;
    Ok(boundary_after_check(
        piece,
        boundary_max_tol(session),
        check_closure,
    ))
}

/// Both-sides loop of `IGESControl_IGESBoundary::Transfer` (`cxx:282-509`),
/// before `Check`. `reverse3d` is `Sense == 2` on entity 141
/// (`IGESToBRep_TopoCurve.cxx:1445`). Entity 142 passes false.
pub(in crate::iges_read) fn boundary_both(
    session: &mut Session<'_>,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
    seq3d: &[usize],
    seq2d: &[usize],
    pref3d: bool,
    pref2d: bool,
    reverse3d: bool,
) -> Result<BoundaryPiece, String> {
    let maxtol = boundary_max_tol(session);
    let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
    let mut wire3 = Vec::new();
    let mut wire2 = Vec::new();
    let mut out = Vec::new();
    let mut ok_curve = true;
    let mut ok3 = true;
    let mut ok2 = true;
    for (i, &child3) in seq3d.iter().enumerate() {
        let child2 = seq2d[i];
        let shape3 = transfer_topo_curve(session, child3)?;
        let mut edges3 = edges_of(&shape3)?;
        if reverse3d {
            edges3 = super::super::iges_boundary::reverse_curves_3d(&edges3, precision, maxtol);
        }
        let mut edges2 = edges_from_2d(session, child2, face, trans, u_fact).unwrap_or_default();
        // `cxx:315` and `cxx:333`: connect before the edge-count test, so a
        // segment that meets the previous one tail-to-tail is reversed first.
        // 142 passes `toreverse3d`/`toreverse2d` false (`IGESToBRep_IGESBoundary.cxx:87`).
        ok3 &= connect_next_wire(&mut wire3, &mut edges3, precision, maxtol);
        ok2 &= connect_next_wire(&mut wire2, &mut edges2, precision, maxtol);
        let mut use3 = !edges3.is_empty();
        let mut use2 = !edges2.is_empty();
        // `cxx:342-358`: edge-count mismatch, or a closed circular arc (`BadCase`).
        if use3
            && use2
            && component_inconsistent(session, child3, child2, edges3.len(), edges2.len())
        {
            use3 = pref3d;
            use2 = pref2d;
        }
        let mut chosen = if use3 && !use2 {
            edges3
        } else if !use3 && use2 {
            edges2
        } else if !use3 && !use2 {
            return Err("IGES 142 wire has no edges: IGESToBRep_TopoCurve.cxx:547-551".into());
        } else {
            // `cxx:369-434`: copy the pcurve, `FixReversed2d`, then
            // `FixSameParameter`. `maxdev` is the edge tolerance that call
            // writes (`cxx:404`). Above `maxtol`, `SetTolerance` forces
            // `Precision::Confusion` on the edge and its vertices
            // (`ShapeFix_ShapeTolerance.cxx:206-210`, default `TopAbs_SHAPE`)
            // and the pcurves are removed. `Preferred3d` false keeps the 2D edge.
            let chosen = edges3;
            let mut use_2d = false;
            let reg = GeometryRegistry::global();
            let face_fix = Face(face.clone());
            for (k, edge2) in edges2.iter().enumerate() {
                if k >= chosen.len() {
                    break;
                }
                if !transfer_pcurve(edge2, &chosen[k], face) {
                    continue;
                }
                let _ = ShapeFixEdge::new().fix_reversed_2d(&chosen[k], face);
                fix_same_parameter(&chosen[k], &face_fix);
                let maxdev = reg.edge_tolerance(&chosen[k].0);
                if maxdev > maxtol {
                    reg.set_edge_tolerance(&chosen[k].0, CONFUSION);
                    if let Some(v) = first_vertex(&chosen[k]) {
                        v.set_tolerance(CONFUSION);
                    }
                    if let Some(v) = last_vertex(&chosen[k]) {
                        v.set_tolerance(CONFUSION);
                    }
                    if pref3d {
                        for edge in chosen.iter().take(k + 1) {
                            reg.remove_pcurves_on_surface(&edge.0, face);
                        }
                    } else {
                        use_2d = true;
                    }
                    break;
                }
            }
            if use_2d { edges2 } else { chosen }
        };
        // `cxx:437-508`. Failure re-transfers every segment, runs
        // `ShapeFix_Wire::Perform` on each wire, and appends the one whose
        // max tolerance is smaller. `Gsewd2d` is then filled with the 3D
        // wire (`cxx:496`), and `okCurve` is set true so the loop continues.
        if !connect_next_wire(&mut out, &mut chosen, precision, maxtol) {
            let (healed3, keep) = rebuild_disconnected(
                session, face, trans, u_fact, seq3d, seq2d, precision,
            );
            wire3 = healed3.clone();
            wire2 = healed3;
            out.extend(keep);
            ok_curve = true;
            ok3 = true;
            ok2 = true;
        }
    }
    Ok(BoundaryPiece {
        out,
        wire3,
        wire2,
        ok: ok_curve,
        ok3,
        ok2,
    })
}

/// `IGESControl_IGESBoundary::Connect` (`cxx:115-142`): `FixConnected` on the
/// incoming wire when `connect_next`, then `ConnectNextWire`.
pub(super) fn connect_wire(
    acc: &mut Vec<Edge>,
    next: &mut Vec<Edge>,
    connect_next: bool,
    precision: f64,
    maxtol: f64,
) -> bool {
    if connect_next {
        let mut wire = TopoBuilder::new().make_wire(next);
        crate::shhealing::fix_connected_closed(&mut wire, precision, maxtol, false);
        *next = crate::topo_tools::edges_of_wire(&wire);
    }
    connect_next_wire(acc, next, precision, maxtol)
}

/// `ShapeAnalysis_ShapeTolerance::Tolerance(wire, 1)`: mode > 0 is the maximum
/// of every edge and vertex tolerance (`cxx:491-492`).
pub(super) fn rebuild_wire_max_tolerance(edges: &[Edge]) -> f64 {
    let reg = GeometryRegistry::global();
    let mut max_tol = 0.0_f64;
    for edge in edges {
        max_tol = max_tol.max(reg.edge_tolerance(&edge.0));
        if let Some(v) = first_vertex(edge) {
            max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
        }
        if let Some(v) = last_vertex(edge) {
            max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
        }
    }
    max_tol
}

/// `IGESControl_IGESBoundary::Check` (`cxx:60-109`). `checkclosure` is true
/// for `TransferBoundaryOnFace` and for 144's `TransferCurveOnFace`.
pub(in crate::iges_read) fn boundary_after_check(
    piece: BoundaryPiece,
    maxtol: f64,
    check_closure: bool,
) -> Vec<Edge> {
    let BoundaryPiece {
        out,
        wire3,
        wire2,
        mut ok,
        mut ok3,
        mut ok2,
    } = piece;
    // `Check` skips the closure test when `checkclosure` is false
    // (`IGESControl_IGESBoundary.cxx:69`). Free entity 142 passes
    // `!isCurveOnSurf` false (`TopoCurve.cxx:545`).
    if check_closure && ok && closure_failed(&out, maxtol) {
        if closure_failed(&wire3, maxtol) {
            ok3 = false;
        } else {
            ok2 = false;
        }
        ok = false;
    }
    if !ok {
        if ok3 && !wire3.is_empty() {
            return wire3;
        }
        if ok2 && !wire2.is_empty() {
            return wire2;
        }
        return Vec::new();
    }
    out
}

fn rebuild_disconnected(
    session: &mut Session<'_>,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
    seq3d: &[usize],
    seq2d: &[usize],
    precision: f64,
) -> (Vec<Edge>, Vec<Edge>) {
    // `cxx:446-488`. `toreverse3d` and `toreverse2d` are false for a 142
    // (`IGESToBRep_IGESBoundary.cxx:85-87`), so neither reverse runs.
    let mut edges3 = Vec::new();
    for &index in seq3d {
        if let Ok(shape) = transfer_topo_curve(session, index) {
            if let Ok(edges) = edges_of(&shape) {
                edges3.extend(edges);
            }
        }
    }
    let mut edges2 = Vec::new();
    for &index in seq2d {
        if let Ok(more) = edges_from_2d(session, index, face, trans, u_fact) {
            edges2.extend(more);
        }
    }
    let face_fix = Face(face.clone());
    let builder = TopoBuilder::new();
    let mut wire3 = builder.make_wire(&edges3);
    let mut wire2 = builder.make_wire(&edges2);
    crate::shhealing::perform_wire_constructor_defaults(&mut wire3, &face_fix, precision);
    crate::shhealing::perform_wire_constructor_defaults(&mut wire2, &face_fix, precision);
    let healed3 = crate::topo_tools::edges_of_wire(&wire3);
    let healed2 = crate::topo_tools::edges_of_wire(&wire2);
    // `ShapeAnalysis_ShapeTolerance::Tolerance(wire, 1)` (`cxx:491-492`):
    // mode > 0 is the maximum of every edge and vertex tolerance.
    let keep = if wire_max_tolerance(&healed3) < wire_max_tolerance(&healed2) {
        healed3.clone()
    } else {
        healed2
    };
    (healed3, keep)
}

fn wire_max_tolerance(edges: &[Edge]) -> f64 {
    let reg = GeometryRegistry::global();
    let mut max_tol = 0.0_f64;
    for edge in edges {
        max_tol = max_tol.max(reg.edge_tolerance(&edge.0));
        if let Some(v) = first_vertex(edge) {
            max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
        }
        if let Some(v) = last_vertex(edge) {
            max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
        }
    }
    max_tol
}
