//! The file-local helpers of `BRepFeat_MakeDPrism.cxx`.
//!
//! Source: `cxx:1159-1313` -- `HeightMax`, `SensOfPrism`, `MajMap`, `TestCurve`.
//!
//! `BRepFeat_MakeDPrism.cxx` carries its own copies of these helpers, and they
//! are not the ones of `BRepFeat_MakePrism.cxx`: `HeightMax` here answers the
//! longest single extent without doubling it and without skipping infinite
//! faces (`cxx:1163-1198`), and `TestCurve` builds its axis from the plane
//! normal of the profile (`cxx:1285-1313`). They are therefore kept apart from
//! `prism/statics.rs` instead of being shared.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpLin, GpPnt};
use occt_geom::{Curve, GeomLine};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::brep_feat::statics::parametric_barycenter;
use crate::brep_tool::BRepTool;
use crate::locope::cs_intersector::CSIntersector;
use crate::locope::sample_edges::sample_edges;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::DPrism;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;

/// `HeightMax(theSbase, theSkface, theSFrom, theSUntil)` (`cxx:1163-1198`).
///
/// The bounding box of the basis, the sketch face and the two limits; the
/// answer is the longest of the three extents (`cxx:1189`).
pub(super) fn height_max(
    sbase: Option<&TopoShape>,
    skface: Option<&Face>,
    sfrom: Option<&TopoShape>,
    suntil: Option<&TopoShape>,
) -> f64 {
    let mut b = BndBox::new();
    // `cxx:1168-1177`: `BRepBndLib::Add` with the default `useTriangulation`.
    if let Some(sb) = sbase {
        brep_bnd_lib::add(sb, &mut b, true);
    }
    if let Some(sk) = skface {
        brep_bnd_lib::add(&sk.0, &mut b, true);
    }
    if let Some(s) = sfrom {
        brep_bnd_lib::add(s, &mut b, true);
    }
    if let Some(s) = suntil {
        brep_bnd_lib::add(s, &mut b, true);
    }

    // `cxx:1182-1189`. `Box.Get` raises on a void box in the reference; the port
    // answers `0.` for the empty extent, which is what the void box means.
    let Some((xmin, ymin, zmin, xmax, ymax, zmax)) = b.get() else {
        return 0.0;
    };
    (xmax - xmin)
        .abs()
        .max((ymax - ymin).abs())
        .max((zmax - zmin).abs())
}

/// `SensOfPrism(C, Until)` (`cxx:1204-1230`): `1` or `-1`, the direction the
/// prism has to grow to reach `Until` along `C`.
pub(super) fn sens_of_prism(c: &Arc<dyn Curve>, until: &TopoShape) -> i32 {
    let mut asi1 = CSIntersector::with_shape(until);
    let scur: Vec<Option<Arc<dyn Curve>>> = vec![Some(c.clone())];
    // `cxx:1208`: a failed intersection leaves `IsDone` false.
    let _ = asi1.perform_curves(&scur);
    let mut sens = 1;
    if asi1.is_done() && asi1.nb_points(1).map(|n| n >= 1).unwrap_or(false) {
        let nb = asi1.nb_points(1).unwrap_or(0);
        let prm1 = asi1.point(1, 1).map(|p| p.parameter()).unwrap_or(0.0);
        let prm2 = asi1.point(1, nb).map(|p| p.parameter()).unwrap_or(0.0);
        if prm1 < 0.0 && prm2 < 0.0 {
            sens = -1;
        }
    } else if parametric_barycenter(until, c) < 0.0 {
        sens = -1;
    }
    sens
}

/// `TestCurve(Base)` (`cxx:1287-1313`): the line through the barycentre of the
/// sampled edges of `Base`, along the normal of the plane the face lies on.
///
/// A face that is not planar, or a base without sampled points, answers the
/// null handle in the reference (`cxx:1303-1307` and the divide by
/// `spt.Length()`); the port answers `None`.
pub(super) fn test_curve(base: &Face) -> Option<Arc<dyn Curve>> {
    let spt = sample_edges(&base.0);
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

    // `cxx:1296-1307`: trimmed surfaces are replaced by their basis, then
    // only a `Geom_Plane` answers a curve.
    let s = BRepTool::face_surface(base)?;
    let s = s.rectangular_trimmed_basis().unwrap_or(s);
    let pln = s.gp_pln()?;

    // `cxx:1308-1311`: `gp_Dir(Position().XDirection() ^ Position().YDirection())`.
    let pos = pln.position();
    let normale = pos.x_direction().crossed(pos.y_direction()).ok()?;
    Some(Arc::new(GeomLine::new(GpLin::from_pnt_dir(bar, normale))))
}

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:1234-1283`).
///
/// The first and the last wire of the draft prism become map keys listing the
/// faces of their shape; then every edge of the base without an entry gets the
/// prism descendants `theP.Shapes(edge)`.
pub(super) fn maj_map(
    b: &TopoShape,
    p: &DPrism,
    map: &mut ShapeValueMap<Vec<TopoShape>>,
    fshape: &mut Option<TopoShape>,
    lshape: &mut Option<TopoShape>,
) {
    // `cxx:1242-1254`.
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

    // `cxx:1256-1271`.
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

    // `cxx:1273-1282`.
    for e in explorer_shapes(b, ShapeType::Edge) {
        if !map.is_bound(&e) {
            map.bind(&e, Vec::new());
            let desc = p.shapes(&e).unwrap_or_default();
            map.bind(&e, desc);
        }
    }
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

/// `const TopoDS_Shape& Base = theP.FirstShape(); TopExp_Explorer exp(Base,
/// TopAbs_FACE); exp.Current(); exp.Next(); exp.More()` (`cxx:230-238`,
/// `cxx:322-330`, `cxx:921-929`): `true` when the first shape of the prism
/// holds exactly one face.
///
/// `Err` is the reference's `exp.Current()` on an empty first shape
/// (`Standard_NoMoreObject`).
pub(super) fn first_shape_has_one_face(p: &DPrism) -> Result<bool, String> {
    let faces = p
        .first_shape()
        .map(|s| explorer_shapes(s, ShapeType::Face))
        .unwrap_or_default();
    let mut it = faces.into_iter();
    if it.next().is_none() {
        return Err("Standard_NoMoreObject".to_string());
    }
    Ok(it.next().is_none())
}
