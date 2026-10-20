use occt_core::precision::{CONFUSION, REAL_SMALL};

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shhealing::{combine_vertex, copy_replace_vertices_with, first_vertex, last_vertex};
use crate::topo_tools::is_same;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_basic_curve::eps_geom;
use super::super::iges_to_brep_curve_and_surface::Session;

use super::boundary::boundary_max_tol;
use super::pcurve_transfer::edges_of;
use super::transfer_topo_curve;

pub(super) fn transfer_composite(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
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
