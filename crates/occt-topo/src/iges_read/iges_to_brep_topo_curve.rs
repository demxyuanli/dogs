//! `IGESToBRep_TopoCurve` for entities 110, 100, 126, and 102, and
//! `TransferCurveOnFace` (`IGESToBRep_TopoCurve.cxx:512-578`).
//!
//! The 142 UV arm is `IGESControl_IGESBoundary::Transfer`
//! (`IGESControl_IGESBoundary.cxx:146-536`). `read.stdsameparameter.mode`
//! is 0 (`XSAlgo.cxx:46`), so the copied pcurve goes through `FixReversed2d`
//! then `FixSameParameter(edge)`. A failed `ConnectNextWire` rebuilds both
//! wires with `ShapeFix_Wire::Perform` (`cxx:441-508`). `Check`
//! (`cxx:60-109`) then keeps that wire, or the 3D or 2D wire when closure fails.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, TrsfForm};
use occt_core::toploc::TopLocLocation;
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::shhealing::{
    combine_vertex, copy_replace_vertices_with, first_vertex, fix_same_parameter, last_vertex,
    ShapeBuildEdge, ShapeFixEdge,
};

use super::iges_boundary::{closure_failed, connect_next_wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::is_same;

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{
    basic_curve_2d, basic_curve_3d, convert_location, edge_from_span, entity_location, eps_geom,
    Curve2dSpan, CurveSpan,
};
use super::iges_to_brep_curve_and_surface::{is_basic_curve, Session};

pub(in crate::iges_read) fn transfer_topo_curve(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `IGESToBRep_TopoCurve::TransferTopoCurve` (`cxx:859-905`).
    if let Some(hit) = session.cache.get(&index) {
        return Ok(hit.clone());
    }
    let ent = session.model.entity(index)?;
    let type_number = ent.type_number;
    let shape = match type_number {
        110 | 100 | 104 | 106 | 112 | 126 => {
            // `TransferTopoBasicCurve` sets `ModeTransfer(false)` (`cxx:967-968`),
            // so a directory 124 is applied inside `TransferBasicCurve`.
            let span = basic_curve_3d(session.model, index)?;
            transfer_topo_basic_curve(&span)?
        }
        102 => transfer_composite(session, index)?,
        142 => super::curve_on_surface::transfer_curve_on_surface(session, index)?,
        141 => super::bounded_surface::transfer_boundary(session, index)?,
        116 => transfer_point(session, index)?,
        130 => super::offset_curve::transfer_offset_curve(session, index)?,
        _ if is_basic_curve(type_number) => {
            return Err(format!(
                "UNPORTED basic curve {type_number}: IGESToBRep_TopoCurve.cxx:871-873"
            ));
        }
        _ => {
            return Err(format!(
                "UNPORTED topo curve {type_number}: IGESToBRep_TopoCurve.cxx:859"
            ));
        }
    };
    session.cache.insert(index, shape.clone());
    Ok(shape)
}

/// `read.iges.bspline.continuity` (`IGESData.cxx:66`) -- the
/// `IGESToBRep_CurveAndSurface::myContinuity` that `GetContinuity()` returns
/// (`IGESToBRep_CurveAndSurface.cxx:58-64`).
const CONTINUITY: i32 = 1;

/// `TransferTopoBasicCurve` (`IGESToBRep_TopoCurve.cxx:950-1123`). One edge,
/// except that a C0 B-spline is replaced by its C1 pieces (`cxx:1013-1038`);
/// the result is then that single edge (`cxx:1108-1110`) or the wire of the
/// pieces (`cxx:1099-1104`).
///
/// `TheCurves` holds exactly the one `TransferBasicCurve` result: the
/// `GetModeApprox()` arm (`cxx:979-991`) is off because `myModeApprox` is set
/// false by the constructor and reset false by `Init`
/// (`IGESToBRep_CurveAndSurface.cxx:64`, `:110`), so `ApproxBSplineCurve` is
/// unreachable and unported.
fn transfer_topo_basic_curve(span: &CurveSpan) -> Result<TopoShape, String> {
    // `cxx:1015-1018`: `IsKind(STANDARD_TYPE(Geom_BSplineCurve))` &&
    // `Continuity() < GeomAbs_C1` && `GetContinuity() >= 1`.
    let mut pieces = None;
    if CONTINUITY >= 1
        && occt_geom::convert_bspl::is_bspline(span.curve.as_ref())
        && span.curve.continuity() < crate::shape_upgrade::GEOMABS_C1
    {
        // `cxx:1021-1022`:
        // `ShapeAlgo::AlgoContainer()->C0BSplineToSequenceOfC1BSplineCurve`.
        pieces = crate::shape_upgrade::c0_bspline_to_sequence_of_c1_bspline_curve(
            span.curve.as_ref(),
        );
    }
    let Some(pieces) = pieces else {
        // `cxx:1039-1053`: `BRepBuilderAPI_MakeEdge ME(mycurve, a, b)`.
        return edge_from_span(span);
    };

    // `cxx:1023-1037`: `BRepBuilderAPI_MakeEdge ME(seqBS->Value(i))` ->
    // `BRepLib_MakeEdge::Init(C)` -> `Init(C, FirstParameter, LastParameter)`
    // (`BRepLib_MakeEdge.cxx:505-508`).
    let mut edges = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let pspan = CurveSpan {
            first: piece.first_parameter(),
            last: piece.last_parameter(),
            curve: piece,
        };
        edges.push(Edge(edge_from_span(&pspan)?));
    }
    if edges.is_empty() {
        return Err(
            "IGES_1156 TopoBasicCurve produced no curve: IGESToBRep_TopoCurve.cxx:1000-1008".into(),
        );
    }
    if edges.len() == 1 {
        // `cxx:1108-1110`: `myedge = TopoDS::Edge(TopoDS_Iterator(myshape).Value())`.
        return Ok(edges.pop().expect("one edge").0);
    }
    // `cxx:1059-1064`: `ShapeFix_Wire sfw; sfw.Load(sewd); sfw.FixConnected();
    // myshape = sewd->Wire()`. `myClosedMode` is true
    // (`ShapeFix_Wire.cxx:117-126`, `:167-172`), so every junction including
    // the last-to-first one is visited, and the default `prec = -1.` resolves
    // to `MaxTolerance() = Precision::Confusion()` (`ShapeFix_Wire.cxx:1485`,
    // `ShapeFix_Root.cxx:24-28`).
    let mut wire = TopoBuilder::new().make_wire(&edges);
    let _ = crate::shhealing::fix_connected_closed(&mut wire, CONFUSION, CONFUSION, true);
    // `cxx:1099-1104`: more than one edge -> the wire (plus the `IGES_1120`
    // warning).
    Ok(wire.0)
}

fn transfer_point(session: &Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferPoint` (`IGESToBRep_TopoCurve.cxx:129-171`).
    // `GetModeTransfer` is true (`IGESToBRep_Actor.cxx:167`), so the point is
    // `Value()`, then scaled, and a directory 124 is `Move`d (`cxx:157-168`).
    let xyz = match &session.model.entity(index)?.data {
        EntityData::Point { xyz } => *xyz,
        _ => return Err("IGES 116 has no point parameters".into()),
    };
    let unit = session.model.unit_factor;
    let point = GpPnt::new(xyz[0] * unit, xyz[1] * unit, xyz[2] * unit);
    let mut shape: TopoShape = TopoBuilder::new().make_vertex(point, CONFUSION).into();
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, unit) {
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    }
    Ok(shape)
}

fn transfer_composite(session: &mut Session<'_>, index: usize) -> Result<TopoShape, String> {
    // `TransferCompositeCurveGeneral` (`IGESToBRep_TopoCurve.cxx:210-347`), 3D.
    let curves = match &session.model.entity(index)?.data {
        EntityData::Composite { curves } => curves.clone(),
        _ => return Err("IGES 102 has no curve list".into()),
    };
    let mut edges = Vec::new();
    let mut prev: Option<usize> = None;
    for (n, child) in curves.into_iter().enumerate() {
        if prev == Some(child) {
            // Same component as the previous one (`cxx:240-250`).
            continue;
        }
        prev = Some(child);
        let shape = transfer_topo_curve(session, child).map_err(|e| {
            format!("IGES 102 component {} failed: {e}", n + 1)
        })?;
        edges.extend(edges_of(&shape)?);
    }
    if edges.is_empty() {
        return Err("IGES 102 produced no edges: IGESToBRep_TopoCurve.cxx:252-257".into());
    }
    let mut wire = TopoBuilder::new().make_wire(&edges);
    // `ClosedWireMode = false`, then `FixConnected(maxtol)` and
    // `FixConnected(1, precision)` (`IGESToBRep_TopoCurve.cxx:339-343`).
    fix_composite_connected(session, &mut wire);
    Ok(wire.0)
}

/// `ShapeFix_Wire::FixConnected` as called from `TransferCompositeCurveGeneral`
/// (`IGESToBRep_TopoCurve.cxx:339-343`). The wire tool is default-constructed,
/// so `myPrecision` is `Precision::Confusion()` (`ShapeFix_Root.cxx:26`,
/// `ShapeAnalysis_Wire.cxx:91`) and `myTopoMode` is false (`ClearModes`).
/// `ClosedWireMode` is false, so the prec overload stops before junction 1
/// (`ShapeFix_Wire.cxx:565-569`); the next call closes that junction at
/// `GetEpsGeom() * GetUnitFactor()`.
fn fix_composite_connected(session: &Session<'_>, wire: &mut Wire) {
    let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
    let maxtol = boundary_max_tol(session);
    let nb = match edges_of(&wire.0) {
        Ok(edges) => edges.len(),
        Err(_) => return,
    };
    if nb == 0 {
        return;
    }
    if nb >= 2 {
        for i in (2..=nb).rev() {
            fix_connected_junction(wire, CONFUSION, maxtol, i);
        }
    }
    fix_connected_junction(wire, CONFUSION, precision, 1);
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ConnectedCheck {
    None,
    /// `ShapeAnalysis_Wire.cxx:727`: distance within `gp::Resolution()`.
    Done1,
    /// `cxx:730`: distance within `myPrecision`.
    Done2,
    /// `cxx:733`: distance within the caller's `prec`.
    Done3,
}

/// `ShapeAnalysis_Wire::CheckConnected(num, prec)` (`cxx:693-762`).
fn check_connected(wire: &Wire, precision: f64, prec: f64, num: usize) -> ConnectedCheck {
    let Ok(edges) = edges_of(&wire.0) else {
        return ConnectedCheck::None;
    };
    let nb = edges.len();
    if nb < 1 {
        return ConnectedCheck::None;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    if n1 > nb || n2 > nb {
        return ConnectedCheck::None;
    }
    let (Some(v1), Some(v2)) = (
        last_vertex(&edges[n1 - 1]),
        first_vertex(&edges[n2 - 1]),
    ) else {
        return ConnectedCheck::None;
    };
    if is_same(&v1.0, &v2.0) {
        return ConnectedCheck::None;
    }
    let dist = BRepTool::vertex_point(&v1).distance(&BRepTool::vertex_point(&v2));
    if dist <= REAL_SMALL {
        return ConnectedCheck::Done1;
    }
    if dist <= precision {
        return ConnectedCheck::Done2;
    }
    if dist <= prec {
        return ConnectedCheck::Done3;
    }
    ConnectedCheck::None
}

/// `ShapeFix_Wire::FixConnected(num, prec, false)` (`ShapeFix_Wire.cxx:1476-1608`).
fn fix_connected_junction(wire: &mut Wire, precision: f64, prec: f64, num: usize) -> bool {
    let check = check_connected(wire, precision, prec, num);
    if check == ConnectedCheck::None {
        return false;
    }
    let Ok(edges) = edges_of(&wire.0) else {
        return false;
    };
    let nb = edges.len();
    if nb == 0 {
        return false;
    }
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let (Some(v1), Some(v2)) = (last_vertex(&e1), first_vertex(&e2)) else {
        return false;
    };
    let combined: Vertex = if check == ConnectedCheck::Done1 {
        let v2_is_last = last_vertex(&e2)
            .map(|l| is_same(&v2.0, &l.0))
            .unwrap_or(false);
        if v2_is_last { v2 } else { v1 }
    } else {
        combine_vertex(&v1, &v2, 1.0001)
    };
    if nb < 2 {
        let tmp = copy_replace_vertices_with(&e2, Some(&combined), Some(&combined));
        set_wire_edge(wire, n2, &tmp);
    } else {
        let tmp_e2 = copy_replace_vertices_with(&e2, Some(&combined), None);
        set_wire_edge(wire, n2, &tmp_e2);
        let first_is_last = first_vertex(&e2)
            .zip(last_vertex(&e2))
            .map(|(f, l)| is_same(&f.0, &l.0))
            .unwrap_or(false);
        if check != ConnectedCheck::Done1 || first_is_last {
            let tmp_e1 = copy_replace_vertices_with(&e1, None, Some(&combined));
            set_wire_edge(wire, n1, &tmp_e1);
        }
    }
    true
}

fn set_wire_edge(wire: &mut Wire, at: usize, edge: &Edge) {
    let pos = if at == 0 {
        wire.0
            .tshape
            .read()
            .expect("lock")
            .children
            .len()
            .saturating_sub(1)
    } else {
        at - 1
    };
    let mut ts = wire.0.tshape.write().expect("lock");
    if pos < ts.children.len() && ts.children[pos].is_edge() {
        ts.children[pos] = edge.0.clone();
        ts.set_modified(true);
    }
}

/// `TransferCurveOnFace` (`IGESToBRep_TopoCurve.cxx:512-578`) through
/// `IGESControl_IGESBoundary::Transfer` (`cxx:146-434`).
/// `read.surfacecurve.mode` stays 0 (`IGESToBRep_CurveAndSurface.cxx:67`).
pub(in crate::iges_read) fn transfer_curve_on_face(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
    check_closure: bool,
) -> Result<Wire, String> {
    let (curve_uv, curve_3d, preference) = match &session.model.entity(index)?.data {
        EntityData::CurveOnSurface {
            curve_uv,
            curve_3d,
            preference,
            ..
        } => (*curve_uv, *curve_3d, *preference),
        _ => return Err(format!("IGES entity {index} is not a 142")),
    };
    let seq3d = match curve_3d {
        Some(c) => curve_sequence(session.model, c)?,
        None => Vec::new(),
    };
    let seq2d = match curve_uv {
        Some(c) => curve_sequence(session.model, c)?,
        None => Vec::new(),
    };
    // `cxx:199-207`: `GetSurfaceCurve` is 0, so only an empty sequence drops a side.
    let mut use3d = !seq3d.is_empty();
    let mut use2d = !seq2d.is_empty();
    let (pref3d, pref2d) = file_preference(preference);
    if use3d && use2d && seq3d.len() != seq2d.len() {
        // `cxx:233-237`: lengths differ, keep the preferred side.
        use3d = pref3d;
        use2d = pref2d;
    }
    if use3d && !use2d {
        let shape = transfer_topo_curve(session, curve_3d.expect("3d"))?;
        let edges = edges_of(&shape)?;
        if edges.is_empty() {
            return Err("IGES 142 wire has no edges: IGESToBRep_TopoCurve.cxx:547-551".into());
        }
        return Ok(TopoBuilder::new().make_wire(&edges));
    }
    if !use3d && use2d {
        let mut edges = Vec::new();
        for child in seq2d {
            edges.extend(edges_from_2d(session, child, face, trans, u_fact)?);
        }
        if edges.is_empty() {
            return Err("IGES 142 UV wire has no edges: IGESControl_IGESBoundary.cxx:266-280".into());
        }
        return Ok(TopoBuilder::new().make_wire(&edges));
    }
    if !use3d && !use2d {
        return Err("IGES 142 has no 3D curve: IGESToBRep_TopoCurve.cxx:547-551".into());
    }
    // Both representations (`cxx:282-434`): 3D plus pcurve, or the 2D edge
    // when the pair is inconsistent or the same-parameter deviation exceeds
    // `GetMaxTol`.
    let edges = edges_when_both(
        session, face, trans, u_fact, &seq3d, &seq2d, pref3d, pref2d, check_closure,
    )?;
    if edges.is_empty() {
        return Err("IGES 142 wire has no edges: IGESToBRep_TopoCurve.cxx:547-551".into());
    }
    Ok(TopoBuilder::new().make_wire(&edges))
}

fn file_preference(preference: i32) -> (bool, bool) {
    // `TransferCurveOnFace` (`cxx:529-536`) then `cxx:209-231`.
    // `surfcurv` stays 0, so the file preference (else the final else) decides.
    let file = if preference == 1 {
        2
    } else if preference == 2 {
        3
    } else {
        0
    };
    let mut pref3d = true;
    let mut pref2d = true;
    if file == 2 {
        pref3d = false;
    } else if file == 3 {
        pref2d = false;
    } else {
        pref3d = false;
    }
    (pref3d, pref2d)
}

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

fn edges_when_both(
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
            edges3 = super::iges_boundary::reverse_curves_3d(&edges3, precision, maxtol);
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

/// An edge carrying only a pcurve on `face`, the shape
/// `ShapeBuild_Edge::MakeEdge(edge, C2d, face, first, last)`
/// (`ShapeBuild_Edge.cxx:666-702`) creates. The port's `EdgeGeom::curve` is
/// non-nullable, so the pcurve is wrapped in a `CurveOnSurface` 3D
/// representation, and `EdgeGeomCore::curve3d_null` records that this wrapper
/// is only `BRepAdaptor_Curve(E, F)`'s evaluator, not a `BRep_Curve3D` curve.
///
/// The representation the wrapper stands in for is OCCT's edge that holds a
/// `BRep_Curve3D` whose curve is null plus a `BRep_CurveOnSurface`. The null
/// `BRep_Curve3D` is what `BRep_Builder::Degenerated(E, true)` installs
/// (`BRep_Builder.cxx:1073-1085` `UpdateCurves(curves, nullCurve, loc)`), and
/// its range is the `BRep_GCurve` default `[-RealLast(), RealLast()]`
/// (`BRep_GCurve` ctor). The two are equivalent for every range query the port
/// makes: `BRep_Tool::Range` skips that null-curve representation and falls
/// through to the curve-on-surface one (`BRep_Tool.cxx:931-960`), and
/// `BRep_Tool::Curve` returns null with `[-RealLast(), RealLast()]`
/// (`BRep_Tool.cxx:172-196`), which is `curve3d_null` here. Verified
/// against the oracle: on the IGES sphere face 0 the port's two pcurve-only
/// edges carry the same `pcrange = [0, 2*PI]` and the same `BRep_Tool::Range`
/// value OCCT reports (`occt_probe --igespc`, `erange`), and on the cone face 0
/// edge 0 the oracle reports `null3d=0 c3=Geom_BSplineCurve` because
/// `ShapeFix_Wire::FixEdgeCurves`'s `FixAddCurve3d` arm
/// (`ShapeFix_Wire.cxx:881-918`) rebuilt it from that pcurve.
fn pcurve_edge(face: &TopoShape, pc: Arc<dyn Curve2d>, first: f64, last: f64) -> Edge {
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .expect("IGES face has no surface");
    // `ShapeBuild_Edge().MakeEdge(newedge, curve->Reversed(), face, ...)`
    // (`IGESToBRep_IGESBoundary.cxx:408-421`) -> `BRepBuilderAPI_MakeEdge(C2d,
    // S, p1, p2)` -> `BRepLib_MakeEdge::Init(C2d, S, null, null, p1, p2)`
    // (`BRepLib_MakeEdge.cxx:905-1094`). `ShapeBuild_Edge::MakeEdge` only
    // prints in `OCCT_DEBUG` and leaves `newedge` null when `!ME.IsDone()`
    // (`ShapeBuild_Edge.cxx:859-865`), and `ReverseCurves2d` then stores that
    // null edge; the caller here keeps the un-reversed pcurve range instead.
    let Some(build) = crate::brep_lib_make_edge::init_2d(&pc, &surf, first, last) else {
        let mut plain = TopoBuilder::new().make_edge(
            Arc::new(CurveOnSurface::new(pc.clone(), surf, first, last)),
            first,
            last,
        );
        GeometryRegistry::global().set_curve3d_null(&plain.0, true);
        attach_pcurve(&plain, face, pc, first, last);
        return plain;
    };
    let curve3d: Arc<dyn occt_geom::Curve> =
        Arc::new(CurveOnSurface::new(pc.clone(), surf, build.first, build.last));
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve3d, build.first, build.last);
    // `BRep_Tool::Curve` answers null on this edge: `ShapeBuild_Edge::MakeEdge`
    // leaves a `BRep_Curve3D` holding a null curve (`BRep_Builder::UpdateCurves`,
    // `BRep_Builder.cxx:86-92`). `curve` above is only what
    // `BRepAdaptor_Curve(E, F)` builds from the pcurve.
    GeometryRegistry::global().set_curve3d_null(&edge.0, true);
    crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    if build.reverse {
        edge.0.set_orientation(crate::abs::Orientation::Reversed);
    }
    attach_pcurve(&edge, face, pc, build.first, build.last);
    edge
}

/// `IGESToBRep_IGESBoundary::ReverseCurves2d` (`cxx:394-426`).
///
/// `ShapeExtend_WireData::Reverse(face)` (`ShapeExtend_WireData.cxx:545-572`)
/// reverses the edge order and every edge orientation before the pcurves are
/// reversed (`Reverse`, `ShapeExtend_WireData.cxx:483-506`), then runs
/// `ComputeSeams(true)` (`:163-222`) and `SwapSeam` on the seam occurrences
/// (`:511-543`). That seam half is not carried here. `SwapSeam` exchanges the
/// two PCurve representations an edge holds for the face (`:530-542`), while
/// every edge of this list is rebuilt by `pcurve_edge`, which gives it its own
/// TShape and a single pcurve registration for the face, so `ComputeSeams` (it
/// pairs a FORWARD and a REVERSED occurrence of the SAME TShape,
/// `ShapeExtend_WireData.cxx:190-215`) has nothing to bind. The port's seam
/// bookkeeping lives in `shape_fix_compose_shell::wire_data::swap_seam`.
///
/// Reachability: the `toreverse2d` arm is only fed by the 508 Loop path
/// (`IGESToBRep_BRepEntity.cxx:394-427`), and none of the fourteen
/// `_igesdump/oracle` models carries an entity 502, 504, 508, 510, 514 or 186.
fn reverse_curves_2d(edges: &[Edge], face: &TopoShape) -> Vec<Edge> {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let mut built = Vec::with_capacity(edges.len());
    for edge in edges.iter().rev() {
        let orientation = edge.0.orientation().reversed();
        let Some(pc) = reg.edge_pcurves(&edge.0, key).into_iter().next() else {
            let mut plain = edge.clone();
            plain.0.set_orientation(orientation);
            built.push(plain);
            continue;
        };
        let (p1, p2) = reg
            .pcurve_range(&edge.0, key)
            .unwrap_or_else(|| reg.edge_parameters(&edge.0));
        let (first, last) = if pc.is_periodic() {
            (pc.reversed_parameter(p2), pc.reversed_parameter(p1))
        } else {
            (
                pc.first_parameter().max(pc.reversed_parameter(p2)),
                pc.last_parameter().min(pc.reversed_parameter(p1)),
            )
        };
        let rev: Arc<dyn Curve2d> = Arc::from(pc.reversed());
        let mut neu = pcurve_edge(face, rev, first, last);
        neu.0.set_orientation(orientation);
        built.push(neu);
    }
    built
}

/// `IGESControl_IGESBoundary::Connect` (`cxx:115-142`): `FixConnected` on the
/// incoming wire when `connect_next`, then `ConnectNextWire`.
fn connect_wire(
    acc: &mut Vec<Edge>,
    next: &mut Vec<Edge>,
    connect_next: bool,
    precision: f64,
    maxtol: f64,
) -> bool {
    if connect_next {
        let mut wire = TopoBuilder::new().make_wire(next);
        crate::shhealing::fix_connected_closed(&mut wire, precision, maxtol, false);
        *next = crate::topo_tools_full::edges_of_wire(&wire);
    }
    connect_next_wire(acc, next, precision, maxtol)
}

/// `IGESToBRep_IGESBoundary` for the `usescurve == true` entry point
/// (`IGESToBRep_IGESBoundary.cxx:94-113`), whose only caller is
/// `IGESToBRep_BRepEntity::TransferLoop` (`IGESToBRep_BRepEntity.cxx:415-421`).
/// `mysewd`, `mysewd3d`, and `mysewd2d` become `out`, `wire3`, and `wire2`;
/// `myCS` is the `Session`; `ConnectNextWire`'s `distmin` / `revsewd` /
/// `revnextsewd` are never read back by OCCT and are not kept.
pub(in crate::iges_read) struct WireDataBoundary {
    pub precision: f64,
    pub maxtol: f64,
    /// `Init`'s `filepreference` (`cxx:52-63`). `TransferLoop` passes 3.
    pub file_preference: i32,
    pub out: Vec<Edge>,
    pub wire3: Vec<Edge>,
    pub wire2: Vec<Edge>,
    pub ok: bool,
    pub ok3: bool,
    pub ok2: bool,
}

impl WireDataBoundary {
    pub(in crate::iges_read) fn new(precision: f64, maxtol: f64, file_preference: i32) -> Self {
        Self {
            precision,
            maxtol,
            file_preference,
            out: Vec::new(),
            wire3: Vec::new(),
            wire2: Vec::new(),
            ok: true,
            ok3: true,
            ok2: true,
        }
    }

    /// One `Transfer` call (`cxx:146-536`). `usescurve` is true, so the 3D side
    /// is the single prebuilt edge `edge3d` and `toreverse3d` is false
    /// (`IGESToBRep_IGESBoundary.cxx:105-111`).
    ///
    /// Returns the `lsewd` of this call.
    pub(in crate::iges_read) fn transfer(
        &mut self,
        session: &mut Session<'_>,
        face: &TopoShape,
        trans: &GpTrsf2d,
        u_fact: f64,
        edge3d: &[Edge],
        param_curves: &[Option<usize>],
        toreverse2d: bool,
        number: usize,
    ) -> Result<Vec<Edge>, String> {
        let len3d = edge3d.len();
        // `cxx:190-197`: every entry of `curves2d` is flattened into `seq2d`.
        let mut seq2d: Vec<usize> = Vec::new();
        for child in param_curves.iter().flatten() {
            seq2d.extend(curve_sequence(session.model, *child)?);
        }
        let len2d = seq2d.len();

        // `cxx:199-207` with `surfcurv == 0` (`read.surfacecurve.mode`).
        let mut translate3d = true;
        let mut translate2d = true;
        if len3d == 0 {
            translate3d = false;
        } else if len2d == 0 {
            translate2d = false;
        }
        // `cxx:210-233` with `surfcurv == 0`: the file preference decides.
        let mut preferred3d = true;
        let mut preferred2d = true;
        if translate3d && translate2d {
            if self.file_preference == 2 {
                preferred3d = false;
            } else if self.file_preference == 3 {
                preferred2d = false;
            } else {
                preferred3d = false;
            }
        }
        // `cxx:235-239`.
        if translate3d && translate2d && len3d != len2d {
            translate3d = preferred3d;
            translate2d = preferred2d;
        }

        let mut lsewd3d = edge3d.to_vec();
        let mut lsewd2d = Vec::new();

        if translate3d && !translate2d {
            // `cxx:240-265`, `usescurve` arm: `Gsewd->Add(scurve3d->Wire())`.
        } else if !translate3d && translate2d {
            // `cxx:266-281`.
            for child in param_curves.iter().flatten() {
                if let Ok(edges) = edges_from_2d(session, *child, face, trans, u_fact) {
                    lsewd2d.extend(edges);
                }
            }
            if toreverse2d {
                lsewd2d = reverse_curves_2d(&lsewd2d, face);
            }
        } else {
            // `cxx:282-437`: both representations with `len3d == 1`.
            // `TC.SetBadCase(false)` / `BadCase()` (`cxx:291-311`): the 508
            // caller hands over a prebuilt edge, so nothing reports the
            // closed-arc bad case for the 3D side.
            let mut lsewd3d_c = std::mem::take(&mut lsewd3d);
            self.ok3 &= connect_next_wire(
                &mut self.wire3,
                &mut lsewd3d_c,
                self.precision,
                self.maxtol,
            );
            lsewd3d = lsewd3d_c;
            let mut bad2d = false;
            let mut local2d = false;
            if let Some(&child2) = seq2d.first() {
                bad2d = circular_bad_case(session, child2);
                match edges_from_2d(session, child2, face, trans, u_fact) {
                    Ok(edges) if !edges.is_empty() => {
                        lsewd2d = edges;
                        if toreverse2d {
                            lsewd2d = reverse_curves_2d(&lsewd2d, face);
                        }
                        self.ok2 &= connect_next_wire(
                            &mut self.wire2,
                            &mut lsewd2d,
                            self.precision,
                            self.maxtol,
                        );
                        local2d = true;
                    }
                    _ => {}
                }
            }
            // `cxx:342-358`. `usescurve` clears `isBSpline`.
            let mut local3d = true;
            if local3d && local2d && (lsewd3d.len() != 1 || lsewd2d.len() != 1 || bad2d) {
                local3d = preferred3d;
                local2d = preferred2d;
            }
            if local3d && local2d {
                // `cxx:369-434`: copy the pcurve, `FixReversed2d`, then
                // `FixSameParameter`. Above `maxtol`, `SetTolerance` forces
                // `Precision::Confusion` and the pcurves are removed;
                // `Preferred3d` false keeps the 2D edge.
                let reg = GeometryRegistry::global();
                let face_fix = Face(face.clone());
                let mut use_2d = false;
                for (k, edge2) in lsewd2d.iter().enumerate() {
                    if k >= lsewd3d.len() {
                        break;
                    }
                    if !transfer_pcurve(edge2, &lsewd3d[k], face) {
                        continue;
                    }
                    let _ = ShapeFixEdge::new().fix_reversed_2d(&lsewd3d[k], face);
                    fix_same_parameter(&lsewd3d[k], &face_fix);
                    let maxdev = reg.edge_tolerance(&lsewd3d[k].0);
                    if maxdev > self.maxtol {
                        reg.set_edge_tolerance(&lsewd3d[k].0, CONFUSION);
                        if let Some(v) = first_vertex(&lsewd3d[k]) {
                            v.set_tolerance(CONFUSION);
                        }
                        if let Some(v) = last_vertex(&lsewd3d[k]) {
                            v.set_tolerance(CONFUSION);
                        }
                        if preferred3d {
                            for edge in lsewd3d.iter().take(k + 1) {
                                reg.remove_pcurves_on_surface(&edge.0, face);
                            }
                        } else {
                            use_2d = true;
                        }
                        break;
                    }
                }
                if use_2d {
                    lsewd3d = lsewd2d.clone();
                }
            } else if local2d {
                lsewd3d = lsewd2d.clone();
            }
        }

        let mut lsewd = lsewd3d.clone();
        // `cxx:437-439`.
        if !connect_next_wire(&mut self.out, &mut lsewd, self.precision, self.maxtol) {
            // `cxx:446-508`: re-transfer every segment, run
            // `ShapeFix_Wire::Perform` on each wire, and keep the one whose
            // maximum tolerance is smaller. `Gsewd2d` is then filled with the
            // 3D wire (`cxx:496`, a 3D wire in both slots).
            let mut rebuilt3d = edge3d.to_vec();
            let mut rebuilt2d = Vec::new();
            for &child in &seq2d {
                if let Ok(edges) = edges_from_2d(session, child, face, trans, u_fact) {
                    rebuilt2d.extend(edges);
                }
            }
            let builder = TopoBuilder::new();
            let face_fix = Face(face.clone());
            let mut wire3 = builder.make_wire(&rebuilt3d);
            let mut wire2 = builder.make_wire(&rebuilt2d);
            crate::shhealing::perform_wire_constructor_defaults(
                &mut wire3,
                &face_fix,
                self.precision,
            );
            crate::shhealing::perform_wire_constructor_defaults(
                &mut wire2,
                &face_fix,
                self.precision,
            );
            rebuilt3d = crate::topo_tools_full::edges_of_wire(&wire3);
            rebuilt2d = crate::topo_tools_full::edges_of_wire(&wire2);
            let heal3 = rebuild_wire_max_tolerance(&rebuilt3d);
            let heal2 = rebuild_wire_max_tolerance(&rebuilt2d);
            let keep = if heal3 < heal2 {
                rebuilt3d.clone()
            } else {
                rebuilt2d
            };
            self.wire3 = rebuilt3d;
            self.wire2 = self.wire3.clone();
            self.out.extend(keep.clone());
            self.ok = true;
            self.ok3 = true;
            self.ok2 = true;
            lsewd = keep;
        }

        if number > 1 {
            // `cxx:512-529`.
            self.ok &= connect_wire(
                &mut self.out,
                &mut lsewd,
                len3d > 1 || len2d > 1,
                self.precision,
                self.maxtol,
            );
            let mut w3 = std::mem::take(&mut self.wire3);
            self.ok3 &= connect_wire(
                &mut w3,
                &mut lsewd3d,
                len3d > 1,
                self.precision,
                self.maxtol,
            );
            self.wire3 = w3;
            let mut w2 = std::mem::take(&mut self.wire2);
            self.ok2 &= connect_wire(
                &mut w2,
                &mut lsewd2d,
                len2d > 1,
                self.precision,
                self.maxtol,
            );
            self.wire2 = w2;
        } else {
            // `cxx:530-535`.
            self.out = lsewd.clone();
            self.wire3 = lsewd3d.clone();
            self.wire2 = lsewd2d;
        }
        Ok(lsewd)
    }
}

/// `ShapeAnalysis_ShapeTolerance::Tolerance(wire, 1)`: mode > 0 is the maximum
/// of every edge and vertex tolerance (`cxx:491-492`).
fn rebuild_wire_max_tolerance(edges: &[Edge]) -> f64 {
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
pub(in crate::iges_read) fn boundary_after_check(    piece: BoundaryPiece,
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
    let healed3 = crate::topo_tools_full::edges_of_wire(&wire3);
    let healed2 = crate::topo_tools_full::edges_of_wire(&wire2);
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

fn component_inconsistent(
    session: &Session<'_>,
    child3: usize,
    child2: usize,
    n3: usize,
    n2: usize,
) -> bool {
    let bspline = is_bspline_curve(session, child3) && is_bspline_curve(session, child2);
    let bad = circular_bad_case(session, child3) || circular_bad_case(session, child2);
    (bspline && n3 != n2) || (!bspline && (n3 != 1 || n2 != 1)) || bad
}

fn is_bspline_curve(session: &Session<'_>, index: usize) -> bool {
    session
        .model
        .entity(index)
        .map(|ent| ent.type_number == 126)
        .unwrap_or(false)
}

fn circular_bad_case(session: &Session<'_>, index: usize) -> bool {
    // `IGESGeom_CircularArc::IsClosed` (`cxx:148-151`) sets `TheBadCase`
    // (`IGESToBRep_TopoCurve.cxx:1116-1119` and `:1304-1307`).
    let Ok(ent) = session.model.entity(index) else {
        return false;
    };
    match &ent.data {
        EntityData::CircularArc { start, end, .. } => {
            (start[0] - end[0]).abs() < PCONFUSION && (start[1] - end[1]).abs() < PCONFUSION
        }
        _ => false,
    }
}

pub(in crate::iges_read) fn transfer_pcurve(from: &Edge, to: &Edge, face: &TopoShape) -> bool {
    // `IGESToBRep::TransferPCurve` (`IGESToBRep.cxx:400-476`) onto an edge
    // that does not already carry a pcurve on `face`.
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    let pcs = reg.edge_pcurves(&from.0, key);
    let Some(pc) = pcs.into_iter().next() else {
        return false;
    };
    let (first, last) = reg
        .pcurve_range(&from.0, key)
        .unwrap_or_else(|| reg.edge_parameters(&from.0));
    attach_pcurve(to, face, pc, first, last);
    true
}

pub(in crate::iges_read) fn curve_sequence(model: &super::iges_data_iges_model::IgesModel, index: usize) -> Result<Vec<usize>, String> {
    // `IGESToBRep::IGESCurveToSequenceOfIGESCurve`.
    match &model.entity(index)?.data {
        EntityData::Composite { curves } => {
            let mut out = Vec::new();
            for child in curves.clone() {
                out.extend(curve_sequence(model, child)?);
            }
            Ok(out)
        }
        _ => Ok(vec![index]),
    }
}

fn mapped_pcurve(
    session: &Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<(Arc<dyn Curve2d>, f64, f64), String> {
    // `Transfer2dTopoBasicCurve` (`cxx:1217-1230`).
    let span = basic_curve_2d(session.model, index)?;
    let (curve, mut first, mut last) = apply_uv_trsf(span, trans);
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .ok_or_else(|| "IGES face has no surface".to_string())?;
    let mut ntrsf = GpTrsf2d::identity();
    if surf.gp_pln().is_some() {
        ntrsf
            .set_scale(&GpPnt2d::new(0.0, 0.0), session.model.unit_factor)
            .map_err(|e| e.to_string())?;
    }
    let pc = ShapeBuildEdge.transform_pcurve(&curve, &ntrsf, u_fact, &mut first, &mut last);
    Ok((pc, first, last))
}

fn apply_uv_trsf(span: Curve2dSpan, trans: &GpTrsf2d) -> (Arc<dyn Curve2d>, f64, f64) {
    if trans.form() == TrsfForm::Identity {
        return (span.curve, span.first, span.last);
    }
    let nf = span.curve.transformed_parameter(span.first, trans);
    let nl = span.curve.transformed_parameter(span.last, trans);
    let mut cloned = span.curve.clone_dyn();
    cloned.transform(trans);
    let basis = match cloned.trimmed_basis() {
        Some(b) => Arc::from(b.clone_dyn()),
        None => Arc::from(cloned),
    };
    let (first, last) = if nf <= nl { (nf, nl) } else { (nl, nf) };
    (
        Arc::new(Geom2dTrimmedCurve::new(basis, first, last)),
        first,
        last,
    )
}

fn attach_pcurve(edge: &Edge, face: &TopoShape, pc: Arc<dyn Curve2d>, first: f64, last: f64) {
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(face);
    reg.set_edge_pcurve(&edge.0, key, pc);
    reg.set_pcurve_range(&edge.0, key, first, last);
    let (f3, l3) = reg.edge_parameters(&edge.0);
    if (f3 - first).abs() > occt_core::precision::PCONFUSION
        || (l3 - last).abs() > occt_core::precision::PCONFUSION
    {
        reg.set_same_parameter(&edge.0, false);
        reg.set_same_range(&edge.0, false);
    }
}

pub(in crate::iges_read) fn edges_from_2d(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<Vec<Edge>, String> {
    // `Transfer2dTopoCurve` (`cxx:910-946`). Entity 130 is
    // `Transfer2dOffsetCurve` (`cxx:941-943`), which may be a wire.
    if session.model.entity(index)?.type_number == 130 {
        let shape = super::offset_curve::transfer_2d_offset_curve(
            session, index, face, trans, u_fact,
        )?;
        return edges_of(&shape);
    }
    edges_from_2d_basic(session, index, face, trans, u_fact)
}

/// `Transfer2dTopoBasicCurve` (`IGESToBRep_TopoCurve.cxx:1127-1309`): the
/// single pcurve edge, except that a C0 B-spline pcurve becomes one edge per
/// C1 piece (`cxx:1232-1256`).
pub(in crate::iges_read) fn edges_from_2d_basic(
    session: &Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<Vec<Edge>, String> {
    let (pc, first, last) = mapped_pcurve(session, index, face, trans, u_fact)?;
    // `cxx:1232-1237`: `IsKind(STANDARD_TYPE(Geom2d_BSplineCurve))` &&
    // `Continuity() < GeomAbs_C1` && `GetContinuity() >= 1`.
    let mut pieces = None;
    if CONTINUITY >= 1
        && pc.is_bspline2d()
        && pc.continuity() < crate::shape_upgrade::GEOMABS_C1
    {
        // `cxx:1238`: `ShapeAlgo::AlgoContainer()->C0BSplineToSequenceOfC1BSplineCurve`.
        pieces = crate::shape_upgrade::c0_bspline2d_to_sequence_of_c1_bspline_curve(pc.as_ref());
    }
    let Some(pieces) = pieces else {
        // `cxx:1257-1274`: `ShapeBuild_Edge().MakeEdge(myedge, mycurve2d, face, a, b)`.
        return Ok(vec![edge_from_pcurve(pc, first, last, face)?]);
    };
    // `cxx:1239-1253`: `ShapeBuild_Edge().MakeEdge(myedge, seqBS->Value(i), face)`
    // -> `BRepLib_MakeEdge(C2d, S)` -> `Init(C2d, S, FirstParameter, LastParameter)`.
    let mut edges = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let f = piece.first_parameter();
        let l = piece.last_parameter();
        edges.push(edge_from_pcurve(piece, f, l, face)?);
    }
    Ok(edges)
}

/// The edge tail of `Transfer2dTopoBasicCurve` (`cxx:1259-1274`) for one
/// pcurve span, shared by the split and unsplit arms.
fn edge_from_pcurve(
    pc: Arc<dyn Curve2d>,
    first: f64,
    last: f64,
    face: &TopoShape,
) -> Result<Edge, String> {
    let surf = GeometryRegistry::global()
        .face_surface(face)
        .ok_or_else(|| "IGES face has no surface".to_string())?;
    // `ShapeBuild_Edge().MakeEdge(myedge, mycurve2d, face, a, b)`
    // (`IGESToBRep_TopoCurve.cxx:1263`) -> `BRepBuilderAPI_MakeEdge(C2d, S,
    // p1, p2)` -> `BRepLib_MakeEdge::Init(C2d, S, null, null, p1, p2)`
    // (`BRepLib_MakeEdge.cxx:905-1094`). `Err` is the `Transfer2dTopoBasicCurve`
    // failure arm for a null edge (`cxx:1265-1271`).
    let Some(build) = crate::brep_lib_make_edge::init_2d(&pc, &surf, first, last) else {
        return Err(
            "IGES_1005 edge construction error: IGESToBRep_TopoCurve.cxx:1265-1271".into(),
        );
    };
    let curve3d: Arc<dyn occt_geom::Curve> =
        Arc::new(CurveOnSurface::new(pc.clone(), surf, build.first, build.last));
    let b = TopoBuilder::new();
    let mut edge = b.make_edge(curve3d, build.first, build.last);
    // `ShapeBuild_Edge::MakeEdge(edge, C2d, S, L, p1, p2)`
    // (`ShapeBuild_Edge.cxx:851-870`): the edge carries the pcurve only, so
    // `BRep_Tool::Curve` is null. `ShapeFix_Edge::FixAddCurve3d`
    // (`ShapeFix_Edge.cxx:618-638`) rebuilds it from the pcurve.
    GeometryRegistry::global().set_curve3d_null(&edge.0, true);
    crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    // `cxx:1089-1092`: a swapped non-periodic range stores the edge REVERSED.
    if build.reverse {
        edge.0.set_orientation(crate::abs::Orientation::Reversed);
    }
    attach_pcurve(&edge, face, pc, build.first, build.last);
    Ok(edge)
}

pub(in crate::iges_read) fn edges_of(shape: &TopoShape) -> Result<Vec<Edge>, String> {
    if shape.is_edge() {
        return Ok(vec![Edge(shape.clone())]);
    }
    if shape.is_wire() {
        let kids = shape.tshape.read().expect("lock").children.clone();
        let mut edges = Vec::new();
        for kid in kids {
            if kid.is_edge() {
                edges.push(Edge(kid));
            } else {
                return Err("IGES wire child is not an edge".into());
            }
        }
        return Ok(edges);
    }
    Err("IGES curve transfer did not produce an edge or a wire".into())
}

