use occt_core::gp::GpTrsf2d;
use occt_core::precision::PCONFUSION;

use crate::builder::TopoBuilder;
use crate::shape::{TopoShape, Wire};

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;

use super::boundary::edges_when_both;
use super::pcurve_transfer::{curve_sequence, edges_from_2d, edges_of};
use super::transfer_topo_curve;

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

pub(super) fn component_inconsistent(
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

pub(super) fn circular_bad_case(session: &Session<'_>, index: usize) -> bool {
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
