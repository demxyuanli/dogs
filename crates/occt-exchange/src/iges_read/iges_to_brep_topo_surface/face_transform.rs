use std::sync::Arc;


use occt_core::gp::GpTrsf;
use occt_core::precision::CONFUSION;

use occt_topo::shape::TopoShape;
use occt_topo::tgeometry::GeometryRegistry;

use super::super::iges_to_brep_basic_curve::entity_location;
use super::super::iges_to_brep_curve_and_surface::Session;

pub(super) fn apply_trimmed_transform(
    session: &Session<'_>,
    index: usize,
    face: &mut TopoShape,
) -> Result<(), String> {
    // `IGESToBRep_TopoSurface.cxx:1321-1371`.
    let Some(g) = entity_location(session.model, index)? else {
        return Ok(());
    };
    // `s1`/`s2`/`s3` are squared column lengths (`cxx:1328-1334`), not moduli.
    let s1 = col_sq(&g, 1);
    let s2 = col_sq(&g, 2);
    let s3 = col_sq(&g, 3);
    if (s1 - s2).abs() > CONFUSION || (s1 - s3).abs() > CONFUSION {
        // `BRepBuilderAPI_GTransform::Perform(face, true)` (`cxx:1336-1343`).
        return super::super::gtransform::gtransform_face(face, &g);
    }
    let tx = g.value(1, 4);
    let ty = g.value(2, 4);
    let tz = g.value(3, 4);
    let tmp = (g.value(1, 1) - 1.0).abs()
        + g.value(1, 2).abs()
        + g.value(1, 3).abs()
        + g.value(2, 1).abs()
        + (g.value(2, 2) - 1.0).abs()
        + g.value(2, 3).abs()
        + g.value(3, 1).abs()
        + g.value(3, 2).abs()
        + (g.value(3, 3) - 1.0).abs()
        + (tx * tx + ty * ty + tz * tz).sqrt();
    if tmp <= CONFUSION {
        return Ok(());
    }
    let mut tr = GpTrsf::identity();
    tr.set_values(
        g.value(1, 1),
        g.value(1, 2),
        g.value(1, 3),
        tx,
        g.value(2, 1),
        g.value(2, 2),
        g.value(2, 3),
        ty,
        g.value(3, 1),
        g.value(3, 2),
        g.value(3, 3),
        tz,
    )
    .map_err(|e| e.to_string())?;
    move_face(face, &tr);
    Ok(())
}

pub(super) fn col_sq(g: &occt_core::gp::GpGTrsf, col: usize) -> f64 {
    let x = g.value(1, col);
    let y = g.value(2, col);
    let z = g.value(3, col);
    x * x + y * y + z * z
}

pub(super) fn move_face(face: &TopoShape, tr: &GpTrsf) {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(face);
    let old_key = reg.debug_repr_key(face_key);
    if let Some(surf) = reg.face_surface(face) {
        let mut moved = surf.clone_dyn();
        moved.transform(tr);
        reg.set_face(face, occt_topo::tgeometry::FaceGeom::new(Arc::from(moved)));
    }
    let new_key = reg.debug_repr_key(face_key);
    let wires = face.tshape.read().expect("lock").children.clone();
    for wire in wires {
        let edges = wire.tshape.read().expect("lock").children.clone();
        for edge in edges {
            if !edge.is_edge() {
                continue;
            }
            if old_key != new_key {
                let mut ts = edge.tshape.write().expect("lock");
                let slot = ts.edge_pcurves_mut();
                slot.curves.rename_key(old_key, new_key);
                slot.ranges.rename_key(old_key, new_key);
            }
            if let Some(curve) = reg.edge_curve(&edge) {
                let moved = Arc::from(curve.transformed(tr));
                edge.tshape.write().expect("lock").edge_core_mut().curve = Some(moved);
            }
            let verts = edge.tshape.read().expect("lock").children.clone();
            for vert in verts {
                if vert.is_vertex() {
                    let p = reg.vertex_point(&vert).transformed(tr);
                    occt_topo::shape::Vertex(vert).set_point(p);
                }
            }
        }
    }
}
