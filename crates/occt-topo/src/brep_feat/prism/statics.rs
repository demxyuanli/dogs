//! The file-local helpers of `BRepFeat_MakePrism.cxx`.
//!
//! Source: `cxx:1032-1256` -- `HeightMax`, `SensOfPrism`, `MajMap`, `TestCurve`
//! and `ToFuse`, plus the `BRepAlgoAPI_*` result helper the `Perform`
//! overloads share.
//!
//! These are not the helpers of `BRepFeat_MakeDPrism.cxx`, which carries its own
//! copies with different semantics (`HeightMax` there answers the longest single
//! extent and `TestCurve` takes the plane normal of the profile), so the two
//! filesets stay apart.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision;
use occt_geom::{Curve, GeomLine};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::brep_feat::statics::parametric_barycenter;
use crate::brep_tool::BRepTool;
use crate::locope::cs_intersector::CSIntersector;
use crate::locope::sample_edges::sample_edges;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::Prism;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;

/// `HeightMax(theSbase, theSkface, theSFrom, theSUntil)` (`cxx:1032-1103`).
///
/// The box of the basis and the sketch face, plus the box of each limit that is
/// not "infinite": a limit is skipped when one of its edges has no vertex
/// (`cxx:1041-1078`). The answer is twice the longest single extent of the six
/// box coordinates (`cxx:1082-1101`).
pub(super) fn height_max(
    sbase: Option<&TopoShape>,
    skface: Option<&Face>,
    sfrom: Option<&TopoShape>,
    suntil: Option<&TopoShape>,
) -> f64 {
    let mut b = BndBox::new();
    // `cxx:1039-1040`: the two shapes that are always present.
    if let Some(sb) = sbase {
        brep_bnd_lib::add(sb, &mut b, true);
    }
    if let Some(sk) = skface {
        brep_bnd_lib::add(&sk.0, &mut b, true);
    }

    // `cxx:1041-1078`.
    for s in [sfrom, suntil].into_iter().flatten() {
        let mut infinite = false;
        for e in explorer_shapes(s, ShapeType::Edge) {
            if explorer_shapes(&e, ShapeType::Vertex).is_empty() {
                infinite = true;
                break;
            }
        }
        if !infinite {
            brep_bnd_lib::add(s, &mut b, true);
        }
    }

    // `cxx:1082-1101`. `Box.Get` raises on a void box in the reference; the port
    // answers `0.` for the empty extent, which is what the void box means.
    let Some((xmin, ymin, zmin, xmax, ymax, zmax)) = b.get() else {
        return 0.0;
    };
    let c = [xmin, zmin, ymin, xmax, zmax, ymax];
    let mut parmin = c[0];
    let mut parmax = c[0];
    for v in c {
        if v > parmax {
            parmax = v;
        }
        if v < parmin {
            parmin = v;
        }
    }
    (2.0 * (parmax - parmin)).abs()
}

/// `SensOfPrism(C, Until)` (`cxx:1109-1131`): `1` or `-1`, the direction the
/// prism has to grow to reach `Until` along `C`.
pub(super) fn sens_of_prism(c: &Arc<dyn Curve>, until: &TopoShape) -> i32 {
    let mut asi1 = CSIntersector::with_shape(until);
    let scur: Vec<Option<Arc<dyn Curve>>> = vec![Some(c.clone())];
    // `cxx:1112`: a failed intersection leaves `IsDone` false.
    let _ = asi1.perform_curves(&scur);
    let mut sens = 1;
    if asi1.is_done() && asi1.nb_points(1).map(|n| n >= 1).unwrap_or(false) {
        let nb = asi1.nb_points(1).unwrap_or(0);
        let prm1 = asi1.point(1, 1).map(|p| p.parameter()).unwrap_or(0.0);
        let prm2 = asi1.point(1, nb).map(|p| p.parameter()).unwrap_or(0.0);
        // `cxx:1117-1121`.
        if prm1 + precision::CONFUSION < 0.0 && prm2 + precision::CONFUSION < 0.0 {
            sens = -1;
        }
    } else if parametric_barycenter(until, c) < 0.0 {
        sens = -1;
    }
    sens
}

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:1135-1176`).
///
/// The first and the last wire of the prism become map keys listing the faces of
/// their shape; then every edge of the base without an entry gets the prism
/// descendants `theP.Shapes(edge)`.
pub(super) fn maj_map(
    b: &TopoShape,
    p: &Prism,
    map: &mut ShapeValueMap<Vec<TopoShape>>,
    fshape: &mut Option<TopoShape>,
    lshape: &mut Option<TopoShape>,
) {
    // `cxx:1144-1155`.
    if let Some(first) = p.first_shape().cloned() {
        if let Some(w) = explorer_shapes(&first, ShapeType::Wire).into_iter().next() {
            *fshape = Some(w.clone());
            map.bind(&w, Vec::new());
            let faces = explorer_shapes(&first, ShapeType::Face);
            if let Some(slot) = map.change(&w) {
                slot.extend(faces);
            }
        }
    }

    // `cxx:1157-1166`.
    if let Some(last) = p.last_shape().cloned() {
        if let Some(w) = explorer_shapes(&last, ShapeType::Wire).into_iter().next() {
            *lshape = Some(w.clone());
            map.bind(&w, Vec::new());
            let faces = explorer_shapes(&last, ShapeType::Face);
            if let Some(slot) = map.change(&w) {
                slot.extend(faces);
            }
        }
    }

    // `cxx:1168-1175`.
    for e in explorer_shapes(b, ShapeType::Edge) {
        if !map.is_bound(&e) {
            map.bind(&e, Vec::new());
            let desc = p.shapes(&e).map(|v| v.to_vec()).unwrap_or_default();
            map.bind(&e, desc);
        }
    }
}

/// `TestCurve(Base, V)` (`cxx:1180-1194`): the line through the barycentre of
/// the sampled edges of `Base`, along `V`.
///
/// The reference divides by `spt.Length()` without a guard, so an empty sample
/// is a division by zero; the port answers `None` there.
pub(super) fn test_curve(base: &TopoShape, v: &GpDir) -> Option<Arc<dyn Curve>> {
    let spt = sample_edges(base);
    if spt.is_empty() {
        return None;
    }
    let mut bar = GpPnt::new(0.0, 0.0, 0.0);
    for p in &spt {
        let c = bar.xyz();
        bar = GpPnt::new(c.x + p.x(), c.y + p.y(), c.z + p.z());
    }
    let n = spt.len() as f64;
    let c = bar.xyz();
    let bar = GpPnt::new(c.x / n, c.y / n, c.z / n);

    Some(Arc::new(GeomLine::from_pnt_dir(bar, *v)))
}

/// `ToFuse(F1, F2)` (`cxx:1198-1258`): coplanarity of the two faces.
///
/// A null face answers false; a rectangular trim is replaced by its basis; a
/// type mismatch answers false. Only two planes then test
/// `gp_Ax3::IsCoplanar` with the confusion / angular tolerances
/// (`cxx:1235-1255`); every other equal-type pair answers false.
pub(super) fn to_fuse(f1: Option<&Face>, f2: Option<&Face>) -> bool {
    let (f1, f2) = match (f1, f2) {
        (Some(a), Some(b)) => (a, b),
        _ => return false,
    };

    // `cxx:1211-1215`: `BRep_Tool::Surface` of both faces. The location is appled
    // right away (`cxx:1240-1248`), so the world-coordinate surface is used.
    let (s1, s2) = match (
        BRepTool::face_surface_world(f1),
        BRepTool::face_surface_world(f2),
    ) {
        (Some(a), Some(b)) => (a, b),
        _ => return false,
    };

    // `cxx:1217-1227`: the trimmed surfaces are replaced by their basis.
    let s1 = s1.rectangular_trimmed_basis().unwrap_or(s1);
    let s2 = s2.rectangular_trimmed_basis().unwrap_or(s2);

    // `cxx:1229-1232`: `DynamicType()` equality.
    if std::any::type_name_of_val(&*s1) != std::any::type_name_of_val(&*s2) {
        return false;
    }

    // `cxx:1234-1255`.
    let (Some(pl1), Some(pl2)) = (s1.gp_pln(), s2.gp_pln()) else {
        return false;
    };
    let pos1 = pl1.position();
    let pos2 = pl2.position();
    pos1.is_coplanar(&pos2, precision::CONFUSION, precision::ANGULAR)
}

/// `Build()` plus `Shape()` of a `BRepAlgoAPI_*` the reference constructs on the
/// stack: the port keeps the operation so `UpdateDescendants` can read its
/// history, and answers the `StdFail_NotDone` the reference's `Shape()` raises.
pub(super) fn boolean_result(
    mut bop: BooleanOperation,
) -> Result<(BooleanOperation, TopoShape), String> {
    bop.build()?;
    let shape = bop
        .shape()
        .cloned()
        .ok_or_else(|| "StdFail_NotDone".to_string())?;
    Ok((bop, shape))
}

/// `myDir` scaled into a `gp_Vec`: `gp_Vec(Dir) * Scalar`.
pub(super) fn scaled_dir(scalar: f64, dir: &GpDir) -> GpVec {
    GpVec::new(dir.x(), dir.y(), dir.z()).multiply_scalar(scalar)
}
