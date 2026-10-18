//! `BRepSweep_Revol` and the `BRepPrimAPI_MakeRevol` constructor IGES calls.
//!
//! Sources: `BRepSweep_Revol.cxx`, `BRepPrimAPI_MakeRevol.cxx`.
//! `check` (`BRepPrimAPI_MakeRevol.cxx:40-43`) runs `BRepLib::BuildCurves3d`
//! before the sweep. `UpdateInnerTolerances` (`cxx:105`) runs after the
//! degenerated-edge unshare.
//! History
//! (`Generated` / `IsDeleted`, `cxx:359-447`) is not used by the IGES caller,
//! which only takes `Shape()`.

use occt_core::gp::{GpAx1, GpLin, GpTrsf};
use occt_core::precision::{ANGULAR, CONFUSION, SQUARE_CONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::extrema_cc::curve_curve_extrema_all_range;
use occt_geom::GeomLine;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::is_same;

use super::num_shape::NumShape;
use super::rotation::Rotation;

const PI: f64 = std::f64::consts::PI;

/// `BRepPrimAPI_MakeRevol(S, A, D)` with `Copy = false`
/// (`BRepPrimAPI_MakeRevol.cxx:48-65`, `IGESToBRep_TopoSurface.cxx:841`).
///
/// `CheckValidity` failure nullifies `myShape` (`cxx:56-59`). That null result
/// is an error here: the caller cannot hold a null `TopoShape`.
pub fn make_revol(shape: &TopoShape, axis: &GpAx1, angle: f64) -> Result<TopoShape, String> {
    // `check` (`cxx:40-43`) and the second call on the `CheckValidity` argument
    // (`cxx:56`). The bool is ignored; a failed edge does not abort the sweep.
    let _ = crate::brep_lib_build_curves3d::build_curves3d(shape);
    if !check_validity(shape, axis) {
        return Err("BRepPrimAPI_MakeRevol::CheckValidity".into());
    }
    let swept = revol_angle(angle);
    if swept <= ANGULAR {
        return Err("BRepSweep_Revol::Constructor".into());
    }
    let mut axe = *axis;
    if angle < 0.0 {
        axe.reverse();
    }
    let mut num = NumShape::dummy();
    if (swept - 2.0 * PI).abs() <= ANGULAR {
        num.init(2, crate::abs::ShapeType::Edge, true, false, false);
    } else {
        num.init(2, crate::abs::ShapeType::Edge, false, false, false);
    }
    let mut gpt = GpTrsf::identity();
    gpt.set_rotation_ax1(&axe, swept)
        .map_err(|_| "BRepSweep_Revol::Location".to_string())?;
    let location = TopLocLocation::composed(&gpt, TopLocLocation::identity());
    let mut profile = shape.clone();
    profile.set_orientation(Orientation::Forward);
    let rotation = Rotation::new(profile, num, location, axe, swept, false)?;
    let mut result = rotation
        .shape()
        .ok_or_else(|| "BRepSweep_Revol::Shape".to_string())?;
    if let Some(err) = rotation.tables.error.borrow().clone() {
        return Err(err);
    }
    unshare_degenerated(&mut result);
    crate::brep_lib_inner_tolerances::update_inner_tolerances(&result);
    Ok(result)
}

/// `BRepSweep_Revol::Angle` (`BRepSweep_Revol.cxx:134-142`).
fn revol_angle(d: f64) -> f64 {
    let mut d = d.abs();
    while d > 2.0 * PI + ANGULAR {
        d -= 2.0 * PI;
    }
    d
}

/// `BRepPrimAPI_MakeRevol::CheckValidity` (`BRepPrimAPI_MakeRevol.cxx:290-341`).
fn check_validity(shape: &TopoShape, axis: &GpAx1) -> bool {
    let mut exp = Explorer::new(shape, crate::abs::ShapeType::Edge);
    while exp.more() {
        let edge = Edge(exp.current().clone());
        if !BRepTool::is_degenerated(&edge) {
            let Some(curve) = BRepTool::edge_curve_world(&edge) else {
                return false;
            };
            let (first, last) = BRepTool::edge_parameters(&edge);
            let mut ratio = 1.0;
            let lin = GpLin::from_pnt_dir(*axis.location(), *axis.direction());
            loop {
                let pp = curve.value(first + (last - first) / ratio);
                let dist = lin.distance(&pp);
                ratio += 1.0;
                if dist >= CONFUSION || ratio >= 100.0 {
                    break;
                }
            }
            if ratio < 100.0 && is_intersect(curve.as_ref(), first, last, axis) {
                return false;
            }
        }
        exp.next();
    }
    true
}

/// `IsIntersect` (`BRepPrimAPI_MakeRevol.cxx:246-286`).
fn is_intersect(curve: &dyn Curve, first: f64, last: f64, axe: &GpAx1) -> bool {
    let axis = GpLin::from_pnt_dir(*axe.location(), *axe.direction());
    if let Some(circ) = curve.gp_circ() {
        let mut r2 = circ.radius() - CONFUSION;
        r2 *= r2;
        if axis.square_distance(&circ.location()) > r2 {
            return false;
        }
    }
    let line = GeomLine::new(axis);
    let par_tol = curve.resolution(CONFUSION);
    let a = first + par_tol;
    let b = last - par_tol;
    let ext = curve_curve_extrema_all_range(curve, &line, first, last, -1.0e6, 1.0e6);
    for e in ext {
        if e.distance * e.distance > SQUARE_CONFUSION {
            continue;
        }
        if a < e.u1 && e.u1 < b {
            return true;
        }
    }
    false
}

/// Shared degenerated edges are copied once per extra face
/// (`BRepPrimAPI_MakeRevol.cxx:115-237`). `BRepTools_ReShape` history is not
/// kept; only the replaced shape is.
fn unshare_degenerated(shape: &mut TopoShape) {
    let mut exp = Explorer::new(shape, crate::abs::ShapeType::Edge);
    let mut groups: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    while exp.more() {
        let edge = exp.current().clone();
        if BRepTool::is_degenerated(&Edge(edge.clone())) {
            if let Some(group) = groups.iter_mut().find(|(key, _)| is_same(key, &edge)) {
                group.1.push(empty_copy_edge(&edge));
            } else {
                groups.push((edge.clone(), vec![edge]));
            }
        }
        exp.next();
    }
    for (_, copies) in &groups {
        if copies.len() < 2 {
            continue;
        }
        let mut pending = copies.clone();
        let mut faces = Explorer::new(shape, crate::abs::ShapeType::Face);
        while faces.more() {
            let face = faces.current().clone();
            let mut edges = Explorer::new(&face, crate::abs::ShapeType::Edge);
            while edges.more() {
                let edge = edges.current().clone();
                if BRepTool::is_degenerated(&Edge(edge.clone())) {
                    if let Some(pos) = pending.iter().position(|c| is_equal(c, &edge)) {
                        pending.remove(pos);
                    } else if let Some(pos) = pending
                        .iter()
                        .position(|c| c.orientation() == edge.orientation())
                    {
                        let found = pending.remove(pos);
                        replace_equal(&face, &edge, &found);
                    }
                }
                edges.next();
            }
            faces.next();
        }
    }
}

/// `TopoDS_Shape::IsEqual` (`TopoDS_Shape.hxx:276-280`) — the comparison
/// `BRepPrimAPI_MakeRevol.cxx:198` (`anIt.Value().IsEqual(anE)`) uses.
fn is_equal(a: &TopoShape, b: &TopoShape) -> bool {
    crate::topo_tools::is_equal(a, b)
}

fn empty_copy_edge(src: &TopoShape) -> TopoShape {
    let mut edge = TopoShape::new(crate::abs::ShapeType::Edge);
    edge.set_location(src.location());
    let (core, pcurves) = {
        let ts = src.tshape.read().expect("tshape");
        (ts.edge_core().cloned(), ts.edge_pcurves().cloned())
    };
    {
        let mut ts = edge.tshape.write().expect("tshape");
        if let Some(core) = core {
            *ts.edge_core_mut() = core;
        }
        if let Some(pc) = pcurves {
            ts.edge_pcurves = Some(pc);
        }
    }
    let stored = src.tshape.read().expect("tshape").children.clone();
    let builder = TopoBuilder::new();
    let mut host = edge.clone();
    host.set_orientation(Orientation::Forward);
    for child in stored {
        builder.add(&mut host, &child);
    }
    host.set_orientation(src.orientation());
    host
}

fn replace_equal(shape: &TopoShape, old: &TopoShape, new: &TopoShape) {
    let kids = shape.tshape.read().expect("tshape").children.clone();
    for (i, child) in kids.iter().enumerate() {
        if is_equal(child, old) {
            let mut ts = shape.tshape.write().expect("tshape");
            let mut placed = new.clone();
            placed.set_orientation(child.orientation());
            placed.set_location(child.location());
            ts.children[i] = placed;
        } else {
            replace_equal(child, old, new);
        }
    }
}
