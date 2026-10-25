//! `ChFi3d::IsTangentFaces`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d.cxx:160-296`, declared in
//! `ChFi3d.hxx:44-48`.

use occt_core::gp::GpDir;
use occt_core::kernel::geomabs::Shape;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::local_analysis::{StatusErrorType, SurfaceContinuity};
use crate::shape::{Edge, Face};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::tool::curvature::nt;

/// `ChFi3d::IsTangentFaces(theEdge, theFace1, theFace2, theOrder)`
/// (`cxx:160-296`).
///
/// `Err` carries the `Standard_Failure` OCCT would let escape from a null
/// surface handle or from `TopOpeBRepTool_TOOL::Nt`; the boolean answer itself
/// never fails.
pub fn is_tangent_faces(
    edge: &Edge,
    face1: &Face,
    face2: &Face,
    order: Shape,
) -> Result<bool, String> {
    // `BRep_Tool::Continuity(theEdge, theFace1, theFace2) != GeomAbs_C0`
    // (`cxx:167-170`). `stored_continuity` is the encoded regularity table and
    // answers `C0` when nothing was encoded.
    if order == Shape::G1 && stored_continuity(edge, face1, face2) != 0 {
        return Ok(true);
    }

    // `TolC0 = max(0.001, 1.5 * BRep_Tool::Tolerance(theEdge))` (`cxx:170`).
    let tol_c0 = 0.001_f64.max(1.5 * BRepTool::edge_tolerance(edge));

    // The two `CurveOnSurface(..., aFirst, aLast)` calls share their output
    // parameters (`cxx:198-216`), so the second read is what `aFirst` / `aLast`
    // hold when the sampling loop starts.
    let (a_c2d1, a_c2d2, a_first, a_last);
    if !is_same(&face1.0, &face2.0)
        && BRepTool::is_closed_edge_face(edge, face1)
        && BRepTool::is_closed_edge_face(edge, face2)
    {
        // Find the edge in face 1: this edge will have the correct orientation
        // (`cxx:180-197`).
        let forward1 = face1.0.oriented(Orientation::Forward);
        let mut an_edge_in_face1 = None;
        for e in explorer_shapes(&forward1, ShapeType::Edge) {
            if is_same(&e, &edge.0) {
                an_edge_in_face1 = Edge::wrap(e);
                break;
            }
        }
        let an_edge_in_face1 = match an_edge_in_face1 {
            Some(e) => e,
            None => return Ok(false),
        };
        let r1 = curve_on_surface_range(&an_edge_in_face1, face1);
        // `TopoDS_Face aFace2 = theFace2; aFace2.Orientation(TopAbs_FORWARD);
        // anEdgeInFace1.Reverse();` (`cxx:200-203`).
        let forward2 = face2.0.oriented(Orientation::Forward);
        let face2_fwd = Face::wrap(forward2)
            .ok_or_else(|| "ChFi3d::IsTangentFaces: face2 is not a face".to_string())?;
        let mut reversed = an_edge_in_face1;
        reversed.0.reverse();
        let r2 = curve_on_surface_range(&reversed, &face2_fwd);
        match (r1, r2) {
            (Some((c1, _, _)), Some((c2, f2, l2))) => {
                a_c2d1 = c1;
                a_c2d2 = c2;
                a_first = f2;
                a_last = l2;
            }
            // `if (aC2d1.IsNull() || aC2d2.IsNull()) return false;`
            // (`cxx:218-221`).
            _ => return Ok(false),
        }
    } else {
        // Obtaining of pcurves of edge on two faces (`cxx:207-216`).
        let r1 = curve_on_surface_range(edge, face1);
        let mut ee = edge.clone();
        if is_same(&face1.0, &face2.0) {
            ee.0.reverse();
        }
        let r2 = curve_on_surface_range(&ee, face2);
        match (r1, r2) {
            (Some((c1, _, _)), Some((c2, f2, l2))) => {
                a_c2d1 = c1;
                a_c2d2 = c2;
                a_first = f2;
                a_last = l2;
            }
            _ => return Ok(false),
        }
    }

    // Obtaining of two surfaces from adjacent faces (`cxx:223-230`).
    let a_surf1 = BRepTool::face_surface_world(face1)
        .ok_or_else(|| "ChFi3d::IsTangentFaces: null surface theFace1".to_string())?;
    let a_surf2 = BRepTool::face_surface_world(face2)
        .ok_or_else(|| "ChFi3d::IsTangentFaces: null surface theFace2".to_string())?;

    // Computation of the number of samples on the edge (`cxx:233-241`).
    let mut a_tool1 = TopolTool::new();
    a_tool1.initialize(a_surf1.clone(), face1)?;
    let mut a_tool2 = TopolTool::new();
    a_tool2.initialize(a_surf2.clone(), face2)?;
    let a_nb_samples = a_tool1.nb_samples().max(a_tool2.nb_samples());

    // Computation of the continuity (`cxx:243-283`).
    let a_delta = (a_last - a_first) / (a_nb_samples as f64 - 1.0);
    let mut nb_not_done = 0;
    for i in 1..=a_nb_samples {
        let mut a_par = a_first + (i - 1) as f64 * a_delta;
        if i == a_nb_samples {
            a_par = a_last;
        }
        let a_cont = SurfaceContinuity::from_curves(
            a_c2d1.as_ref(),
            a_c2d2.as_ref(),
            a_par,
            a_surf1.as_ref(),
            a_surf2.as_ref(),
            order,
            0.001,
            tol_c0,
            0.1,
            0.1,
            0.1,
            0.01,
            10000.0,
        )?;
        if !a_cont.is_done() {
            // `if (theOrder == GeomAbs_C2 && aCont.StatusError() ==
            // LocalAnalysis_NullSecondDerivative) continue;` (`cxx:258-262`).
            if order == Shape::C2 && a_cont.status_error() == StatusErrorType::NullSecondDerivative {
                continue;
            }
            nb_not_done += 1;
            continue;
        }
        if order == Shape::G1 {
            if !a_cont.is_g1()? {
                return Ok(false);
            }
        } else if !a_cont.is_g2()? {
            return Ok(false);
        }
    }

    if nb_not_done == a_nb_samples {
        return Ok(false);
    }

    // Compare normals of tangent faces in the middle point (`cxx:286-295`).
    let mid_par = (a_first + a_last) / 2.0;
    let uv1 = a_c2d1.value(mid_par);
    let uv2 = a_c2d2.value(mid_par);
    // `gp_Dir normal1, normal2;` (`cxx:290`) starts uninitialised; a failed
    // `Nt` leaves it untouched. `GpDir::default_dir()` is the port's stand-in
    // for that indeterminate value.
    let normal1 = nt(&uv1, face1).unwrap_or_else(GpDir::default_dir);
    let normal2 = nt(&uv2, face2).unwrap_or_else(GpDir::default_dir);
    Ok(normal1.dot(&normal2) >= 0.0)
}
