//! The file-local helpers of `BRepFeat_MakeRevol.cxx`: `MajMap`
//! (`cxx:892-936`), `VerifGluedFaces` (`cxx:816-890`) and `ToFuse`
//! (`cxx:938-990`).

use std::sync::Arc;

use occt_core::precision::{ANGULAR, PCONFUSION};
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::brep_feat::statics::tool;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::{CSIntersector, Revol};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;

use super::MakeRevol;

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:892-936`).
pub(super) fn maj_map(
    b: &TopoShape,
    p: &Revol,
    map: &mut ShapeValueMap<Vec<TopoShape>>,
    fshape: &mut Option<TopoShape>,
    lshape: &mut Option<TopoShape>,
) {
    // `cxx:902-914`.
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

    // `cxx:916-926`.
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

    // `cxx:928-935`.
    for e in explorer_shapes(b, ShapeType::Edge) {
        if !map.is_bound(&e) {
            map.bind(&e, Vec::new());
            let desc = p.shapes(&e).map(|v| v.to_vec()).unwrap_or_default();
            map.bind(&e, desc);
        }
    }
}

/// `VerifGluedFaces(theSkface, thePbase, theBCurve, theCurves, theRevol,
/// theMap)` (`cxx:816-890`).
///
/// Cuts the revolution by the tool solid of the sketch face and drops the
/// glued-face map when no resulting solid keeps a face of `thePbase`.
pub(super) fn verif_glued_faces(
    this: &mut MakeRevol,
    the_skface: &Face,
    the_pbase: &TopoShape,
    the_revol: &Revol,
) -> Result<(), String> {
    // `cxx:829-838`.
    let vrairevol = the_revol
        .shape()
        .cloned()
        .ok_or_else(|| "StdFail_NotDone".to_string())?;
    this.take_revol_curves(the_revol);
    let bcurve = this.my_bcurve.clone();
    let scur: Vec<Option<Arc<dyn Curve>>> = match &bcurve {
        Some(c) => vec![Some(c.clone())],
        None => vec![None],
    };

    let mut asi = CSIntersector::with_shape(&the_skface.0);
    let _ = asi.perform_curves(&scur);
    let nb = if asi.is_done() {
        asi.nb_points(1).unwrap_or(0)
    } else {
        0
    };
    if nb >= 1 {
        // `cxx:840-850`.
        let or = asi
            .point(1, 1)
            .ok()
            .map(|p| p.orientation())
            .unwrap_or(crate::abs::Orientation::Forward);
        let mut tools: Vec<TopoShape> = Vec::new();
        if let Ok(p) = asi.point(1, 1) {
            if let Some(s) = tool(&the_skface.0, p.face(), or) {
                tools.push(s.0);
            }
        }
        let comp = TopoBuilder::new().make_compound_of(&tools);
        // `cxx:850-851`.
        let mut tr_p = crate::brep_feat::algo_api::cut(&vrairevol, &comp.0);
        if tr_p.build().is_err() {
            return Ok(());
        }
        let Some(cutsh) = tr_p.shape().cloned() else {
            return Ok(());
        };

        // `cxx:852-879`: every solid must keep a face of `thePbase`.
        let mut glued_faces = true;
        let pbase_faces = crate::topo_tools::faces_of(the_pbase);
        for solid in explorer_shapes(&cutsh, ShapeType::Solid) {
            let mut found = false;
            for fac1 in explorer_shapes(&solid, ShapeType::Face) {
                if pbase_faces
                    .iter()
                    .any(|fac2| crate::topo_tools::is_same(&fac1, &fac2.0))
                {
                    found = true;
                    break;
                }
            }
            if found {
                continue;
            }
            glued_faces = false;
            break;
        }
        // `cxx:880-889`.
        if !glued_faces {
            this.form.my_glued_f.clear();
        }
    }
    Ok(())
}

/// `ToFuse(F1, F2)` (`cxx:938-990`): `true` only when both faces lie on the
/// same plane, with the surfaces taken in world coordinates (`cxx:977-978`).
pub(super) fn to_fuse(f1: Option<&Face>, f2: Option<&Face>) -> bool {
    let (Some(f1), Some(f2)) = (f1, f2) else {
        return false;
    };
    let strip = |f: &Face| -> Option<Arc<dyn Surface>> {
        let s = BRepTool::face_surface_world(f)?;
        Some(s.rectangular_trimmed_basis().unwrap_or(s))
    };
    let (Some(s1), Some(s2)) = (strip(f1), strip(f2)) else {
        return false;
    };
    let (Some(pl1), Some(pl2)) = (s1.gp_pln(), s2.gp_pln()) else {
        return false;
    };
    pl1.position().is_coplanar(&pl2.position(), PCONFUSION, ANGULAR)
}

/// `Build()` plus `Shape()` of a `BRepAlgoAPI_*` built on the stack, shared with
/// the `MakePrism` port.
pub(super) fn boolean_result(
    mut bop: crate::brep_feat::algo_api::BooleanOperation,
) -> Result<(crate::brep_feat::algo_api::BooleanOperation, TopoShape), String> {
    bop.build()?;
    let shape = bop
        .shape()
        .cloned()
        .ok_or_else(|| "StdFail_NotDone".to_string())?;
    Ok((bop, shape))
}
