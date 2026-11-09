//! `ShapeAnalysis_TransferParametersProj::CopyNMVertex(V, toFace, fromFace)`.

use super::prelude::*;

/// `ShapeAnalysis_TransferParametersProj::CopyNMVertex(V, toFace, fromFace)`
/// (`Proj.cxx:715-802`).
///
/// The copy is a `BRep_TVertex::EmptyCopy` (`BRep_TVertex.cxx:33-39`): the point
/// and the tolerance come along, the representation list does not, so the loop
/// at `cxx:742-781` rebuilds it entry by entry -
/// `BRep_PointOnCurveOnSurface` / `BRep_PointOnCurve` are copied verbatim
/// (`cxx:750-761`), a `BRep_PointOnSurface` whose surface is `fromSurf` is
/// folded into `(apar1, apar2)` instead of being copied (`cxx:766-771`), and
/// every other `BRep_PointOnSurface` is copied (`cxx:772-779`).
///
/// When nothing matched or the carrier changed (`cxx:783`), the point is
/// projected onto `toFace`'s surface with `ShapeAnalysis_Surface::ValueOfUV`
/// and the vertex tolerance is widened to `Gap()` (`cxx:785-797`). The closing
/// `BRep_Builder::UpdateVertex(anewV, apar1, apar2, toFace, aTol)` (`cxx:800`)
/// is [`crate::brep_builder_update_vertex::update_vertex_on_face`].
///
/// `TopLoc_Location` is unported (see the module note), so `fromSurf != toSurf`
/// reduces to the surface handle identity and `fromLoc != toLoc` never fires.
///
/// Returns `None` for a vertex that is neither INTERNAL nor EXTERNAL
/// (`cxx:720-723` returns a null `TopoDS_Vertex` there), and for a `toFace`
/// without a surface, where `cxx:785-786` would build a null
/// `ShapeAnalysis_Surface`.
pub fn copy_nm_vertex_face(v: &Vertex, to_face: &Face, from_face: &Face) -> Option<Vertex> {
    if !matches!(v.0.orientation(), Orientation::Internal | Orientation::External) {
        return None; // cxx:720-723
    }
    // `BRep_Tool::Surface(fromFace, fromLoc)` / `(toFace, toLoc)` (`cxx:727-729`).
    let from_surf = BRepTool::face_surface(from_face);
    let to_surf = BRepTool::face_surface(to_face);

    // `cxx:731-732`: `anewV = TopoDS::Vertex(theV.EmptyCopied())`, then
    // `gp_Pnt apv = BRep_Tool::Pnt(anewV)` - the empty copy kept the point.
    let anew = Vertex::wrap(v.0.empty_copied())?;
    let apv = BRepTool::vertex_point(v);
    let registry = GeometryRegistry::global();
    registry.set_vertex(
        &anew.0,
        crate::tgeometry::VertexGeom {
            point: apv,
            tolerance: BRepTool::vertex_tolerance(v),
        },
    );

    // `cxx:742-781`: walk the source representations, filling `anewV`'s list.
    let src_points = match v.0.tshape.read() {
        Ok(ts) => ts.vertex_core().map(|c| c.points.clone()).unwrap_or_default(),
        Err(_) => Vec::new(),
    };
    let mut copied: Vec<PointRepresentation> = Vec::with_capacity(src_points.len());
    let mut has_repr = false;
    let mut apar1 = 0.0;
    let mut apar2 = 0.0;
    for pr in &src_points {
        match pr {
            // `cxx:750-755`.
            PointRepresentation::OnCurveOnSurface { .. } => copied.push(pr.clone()),
            // `cxx:756-761`.
            PointRepresentation::OnCurve { .. } => copied.push(pr.clone()),
            // `cxx:762-779`: only the representation living on `fromSurf` is
            // consumed; the others are copied as they are.
            PointRepresentation::OnSurface { surface, parameter, parameter2 } => {
                if from_surf.as_ref().is_some_and(|s| Arc::ptr_eq(s, surface)) {
                    apar1 = *parameter;
                    apar2 = *parameter2;
                    has_repr = true;
                } else {
                    copied.push(pr.clone());
                }
            }
        }
    }
    {
        let mut ts = anew.0.tshape.write().expect("poisoned TShape lock");
        ts.vertex_core_mut().points = copied;
    }

    // `cxx:782`: `aTol = BRep_Tool::Tolerance(anewV)`.
    let mut a_tol = BRepTool::vertex_tolerance(&anew);
    // `cxx:783`: `!hasRepr || (fromSurf != toSurf || fromLoc != toLoc)`.
    let same_carrier = match (from_surf.as_ref(), to_surf.as_ref()) {
        (Some(a), Some(b)) => Arc::ptr_eq(a, b),
        (None, None) => true,
        _ => false,
    };
    if !has_repr || !same_carrier {
        let s = to_surf?; // `cxx:783`
        let (p2d, gap) = surface_value_of_uv_with_gap(s.as_ref(), &apv, CONFUSION);
        apar1 = p2d.x();
        apar2 = p2d.y();
        if a_tol < gap {
            a_tol = gap + 0.1 * CONFUSION; // `cxx:791-794`
        }
    }

    // `cxx:800`: `aB.UpdateVertex(anewV, apar1, apar2, toFace, aTol)`.
    let _ = crate::brep_builder_update_vertex::update_vertex_on_face(
        &anew.0, apar1, apar2, &to_face.0, a_tol,
    );
    Some(anew)
}

