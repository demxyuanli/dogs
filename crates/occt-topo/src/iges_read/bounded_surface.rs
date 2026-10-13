//! `IGESToBRep_TopoSurface::TransferBoundedSurface` (`cxx:1381-1489`).
//!
//! `ParamSurface` (`cxx:1790-2071`) transfers the basis, then
//! `TransferBoundaryOnFace` (`IGESToBRep_TopoCurve.cxx:1398-1483`) adds each
//! entity 141. The commented `ShapeFix_Face` block (`cxx:1480-1486`) is not
//! ported. `BRepTools::Update` (`cxx:1478`) does not change the face.

use occt_core::toploc::TopLocLocation;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Wire};

use super::iges_boundary::connect_next_wire;
use super::iges_data_iges_model::{BoundaryCurve, EntityData};
use super::iges_to_brep_basic_curve::{convert_location, entity_location, eps_geom};
use super::iges_to_brep_curve_and_surface::{is_topo_surface, Session};
use super::iges_to_brep_topo_curve::{
    boundary_after_check, boundary_both, boundary_max_tol, curve_sequence, edges_from_2d, edges_of,
    transfer_topo_curve, BoundaryPiece,
};
use super::iges_to_brep_topo_surface::transfer_topo_surface;

pub(in crate::iges_read) fn transfer_bounded_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    let (representation, surface, boundaries) = match &session.model.entity(index)?.data {
        EntityData::Bounded {
            representation,
            surface,
            boundaries,
        } => (*representation, *surface, boundaries.clone()),
        _ => return Err("IGES 143 has no bounded-surface parameters".into()),
    };
    // `RepresentationType == 0` warns `IGES_1275` and still transfers the basis
    // (`cxx:1395-1400`).
    if representation == 0 {
        // `IGES_1275` (`cxx:1395-1399`): the basis is still transferred.
    }
    if !is_topo_surface(session.model, surface) {
        return Err("XSTEP_166: IGESToBRep_TopoSurface.cxx:1405-1411".into());
    }
    let myshape = match transfer_topo_surface(session, surface) {
        Ok(shape) => shape,
        Err(_) => return Err("IGES_1005: IGESToBRep_TopoSurface.cxx:1803-1807".into()),
    };
    let face = match param_face(&myshape) {
        ParamFace::Shell(shape) => return Ok(shape),
        ParamFace::Face(face) => face,
        ParamFace::Fail => {
            return Err("IGES_1156: IGESToBRep_TopoSurface.cxx:1460-1467".into());
        }
    };
    let surf = BRepTool::face_surface(&face)
        .ok_or_else(|| "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1851".to_string())?;
    let (trans, u_fact) = super::iges_to_brep_topo_surface::uv_map(session, surface, &face.0)?;
    // `face.EmptyCopy` (`cxx:1472`) keeps the surface and the face location.
    let mut blank = TopoBuilder::new().make_face(surf, &[]);
    if !face.0.location().is_identity() {
        blank.0.move_location(&face.0.location());
    }
    for boundary in boundaries {
        if let Ok(wire) = transfer_boundary_on_face(session, boundary, &blank.0, &trans, u_fact) {
            TopoBuilder::new().add_wire(&mut blank, &wire);
        }
    }
    Ok(blank.0)
}

enum ParamFace {
    Face(Face),
    /// `ParamSurface` returns the shell when it has more than one face
    /// (`cxx:1828-1836`). `TransferBoundedSurface` does the same (`cxx:1447-1456`).
    Shell(TopoShape),
    Fail,
}

fn param_face(shape: &TopoShape) -> ParamFace {
    if shape.is_face() {
        return match Face::wrap(shape.clone()) {
            Some(face) => ParamFace::Face(face),
            None => ParamFace::Fail,
        };
    }
    if shape.is_shell() {
        let faces: Vec<TopoShape> = cumulated_children(shape)
            .into_iter()
            .filter(|child| child.is_face())
            .collect();
        if faces.len() != 1 {
            return ParamFace::Shell(shape.clone());
        }
        return match Face::wrap(faces.into_iter().next().expect("one face")) {
            Some(face) => ParamFace::Face(face),
            None => ParamFace::Fail,
        };
    }
    ParamFace::Fail
}

/// `TransferBoundaryOnFace` (`IGESToBRep_TopoCurve.cxx:1398-1483`).
pub(in crate::iges_read) fn transfer_boundary_on_face(
    session: &mut Session<'_>,
    index: usize,
    face: &TopoShape,
    trans: &occt_core::gp::GpTrsf2d,
    u_fact: f64,
) -> Result<Wire, String> {
    let (boundary_type, preference, curves) = match &session.model.entity(index)?.data {
        EntityData::Boundary {
            boundary_type,
            preference,
            curves,
            ..
        } => (*boundary_type, *preference, curves.clone()),
        _ => return Err(format!("IGES entity {index} is not a 141")),
    };
    let (pref3d, pref2d) = boundary_file_preference(preference);
    let mut acc: Option<BoundaryPiece> = None;
    for curve in &curves {
        if curve.parameters.is_empty() && boundary_type == 1 {
            // `IGES_1135` (`cxx:1431-1435`). The parameter list stays null.
        }
        let piece = transfer_one(session, face, trans, u_fact, curve, pref3d, pref2d)?;
        let len3d = curve_sequence(session.model, curve.model)?.len();
        let len2d: usize = curve
            .parameters
            .iter()
            .map(|p| curve_sequence(session.model, *p).map(|s| s.len()).unwrap_or(0))
            .sum();
        acc = Some(match acc {
            None => piece,
            Some(prev) => connect_pieces(session, prev, piece, len3d, len2d),
        });
    }
    let Some(piece) = acc else {
        return Err("IGES_1095: IGESToBRep_TopoCurve.cxx:1452-1456".into());
    };
    let edges = boundary_after_check(piece, boundary_max_tol(session), true);
    if edges.is_empty() {
        return Err("IGES_1095: IGESToBRep_TopoCurve.cxx:1452-1456".into());
    }
    let mut wire = TopoBuilder::new().make_wire(&edges);
    move_boundary_location(session, index, &mut wire)?;
    Ok(wire)
}

fn transfer_one(
    session: &mut Session<'_>,
    face: &TopoShape,
    trans: &occt_core::gp::GpTrsf2d,
    u_fact: f64,
    curve: &BoundaryCurve,
    pref3d: bool,
    pref2d: bool,
) -> Result<BoundaryPiece, String> {
    let seq3d = curve_sequence(session.model, curve.model)?;
    let mut seq2d = Vec::new();
    for parameter in &curve.parameters {
        seq2d.extend(curve_sequence(session.model, *parameter)?);
    }
    let mut use3d = !seq3d.is_empty();
    let mut use2d = !seq2d.is_empty();
    if use3d && use2d && seq3d.len() != seq2d.len() {
        use3d = pref3d;
        use2d = pref2d;
    }
    let reverse3d = curve.sense == 2;
    if use3d && !use2d {
        let shape = transfer_topo_curve(session, curve.model)?;
        let mut edges = edges_of(&shape)?;
        if reverse3d {
            let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
            edges = super::iges_boundary::reverse_curves_3d(
                &edges,
                precision,
                boundary_max_tol(session),
            );
        }
        if edges.is_empty() {
            return Err("IGES 141 wire has no edges: IGESControl_IGESBoundary.cxx:241-264".into());
        }
        return Ok(BoundaryPiece {
            out: edges.clone(),
            wire3: edges,
            wire2: Vec::new(),
            ok: true,
            ok3: true,
            ok2: true,
        });
    }
    if !use3d && use2d {
        let mut edges = Vec::new();
        for child in &seq2d {
            if let Ok(more) = edges_from_2d(session, *child, face, trans, u_fact) {
                edges.extend(more);
            }
        }
        if edges.is_empty() {
            return Err("IGES 141 UV wire has no edges: IGESControl_IGESBoundary.cxx:266-280".into());
        }
        return Ok(BoundaryPiece {
            out: edges.clone(),
            wire3: Vec::new(),
            wire2: edges,
            ok: true,
            ok3: true,
            ok2: true,
        });
    }
    if !use3d && !use2d {
        return Err("IGES 141 has no curve: IGESControl_IGESBoundary.cxx:199-207".into());
    }
    boundary_both(
        session, face, trans, u_fact, &seq3d, &seq2d, pref3d, pref2d, reverse3d,
    )
}

/// `Connect` (`IGESControl_IGESBoundary.cxx:123-141`) when `number > 1`
/// (`cxx:512-528`).
fn connect_pieces(
    session: &Session<'_>,
    mut acc: BoundaryPiece,
    mut next: BoundaryPiece,
    len3d: usize,
    len2d: usize,
) -> BoundaryPiece {
    let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
    let maxtol = boundary_max_tol(session);
    acc.ok &= next.ok;
    acc.ok3 &= next.ok3;
    acc.ok2 &= next.ok2;
    acc.ok &= connect_wire(
        &mut acc.out,
        &mut next.out,
        len3d > 1 || len2d > 1,
        precision,
        maxtol,
    );
    acc.ok3 &= connect_wire(
        &mut acc.wire3,
        &mut next.wire3,
        len3d > 1,
        precision,
        maxtol,
    );
    acc.ok2 &= connect_wire(
        &mut acc.wire2,
        &mut next.wire2,
        len2d > 1,
        precision,
        maxtol,
    );
    acc
}

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

/// `PreferenceType == 2` stores file preference 2, `== 1` stores 3
/// (`IGESToBRep_TopoCurve.cxx:1414-1421`), then `cxx:221-230`.
fn boundary_file_preference(preference: i32) -> (bool, bool) {
    let file = if preference == 2 {
        2
    } else if preference == 1 {
        3
    } else {
        0
    };
    if file == 2 {
        (false, true)
    } else if file == 3 {
        (true, false)
    } else {
        (false, true)
    }
}

pub(in crate::iges_read) fn transfer_boundary(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    // `TransferBoundary` (`IGESToBRep_TopoCurve.cxx:1314-1394`).
    let (boundary_type, surface) = match &session.model.entity(index)?.data {
        EntityData::Boundary {
            boundary_type,
            surface,
            ..
        } => (*boundary_type, *surface),
        _ => return Err("IGES 141 has no boundary parameters".into()),
    };
    if boundary_type == 0 {
        return Err("IGES_1125: IGESToBRep_TopoCurve.cxx:1325-1330".into());
    }
    let Some(surface) = surface else {
        return Err("XSTEP_124: IGESToBRep_TopoCurve.cxx:1337-1342".into());
    };
    if !is_topo_surface(session.model, surface) {
        return Err("XSTEP_124: IGESToBRep_TopoCurve.cxx:1337-1342".into());
    }
    let myshape = transfer_topo_surface(session, surface)?;
    let face = if myshape.is_face() {
        myshape
    } else if myshape.is_shell() {
        let kids = myshape.tshape.read().expect("lock").children.clone();
        let faces: Vec<TopoShape> = kids.iter().filter(|child| child.is_face()).cloned().collect();
        if kids.len() == 1 && faces.len() == 1 {
            faces.into_iter().next().expect("one face")
        } else {
            return Err("IGES_1061: IGESToBRep_TopoCurve.cxx:1370-1375".into());
        }
    } else {
        return Err("IGES_1156: IGESToBRep_TopoCurve.cxx:1379-1386".into());
    };
    let wrapped = Face(face.clone());
    let surf = BRepTool::face_surface(&wrapped)
        .ok_or_else(|| "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1851".to_string())?;
    let (trans, u_fact) = super::iges_to_brep_topo_surface::uv_map(session, surface, &face)?;
    // `face.EmptyCopy` (`cxx:1391`) keeps the surface and the location.
    let mut blank = TopoBuilder::new().make_face(surf, &[]);
    if !face.location().is_identity() {
        blank.0.move_location(face.location());
    }
    let wire = transfer_boundary_on_face(session, index, &blank.0, &trans, u_fact)?;
    Ok(wire.0)
}

fn move_boundary_location(
    session: &Session<'_>,
    index: usize,
    wire: &mut Wire,
) -> Result<(), String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if tr.form() != occt_core::gp::TrsfForm::Identity {
                wire.0.move_location(&TopLocLocation::composed(
                    &tr,
                    TopLocLocation::identity(),
                ));
            }
        }
    }
    Ok(())
}
