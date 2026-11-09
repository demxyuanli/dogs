//! The vertex-substitution tail shared by `Build` (`cxx:125-178`),
//! `BuildOnContext` (`cxx:405-450`) and `BuildFaces` (`cxx:613-660`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::topo_tools::{children_keep_orientation, explorer_shapes, shape_map_key, ShapeMapKey};

/// `BRep_TVertex` merge of `cxx:166-170` / `:437-441` / `:646-650`: raise the
/// substitute's tolerance to the old vertex's and append the old point
/// representations to it.
fn merge_vertex_geometry(v: &TopoShape, new_v: &TopoShape) {
    // `TV->Tolerance()` / `TV->ChangePoints()` are read off the old vertex
    // first, so the two `TShape` locks are never held at the same time (the
    // substitute can share `V`'s `TShape`).
    let old = match v.tshape.read() {
        Ok(t) => t.vertex_core().cloned(),
        Err(_) => None,
    };
    let Some(old) = old else {
        return;
    };
    let Ok(mut t) = new_v.tshape.write() else {
        return;
    };
    let new_core = t.vertex_core_mut();
    // `if (TV->Tolerance() > NewTV->Tolerance()) NewTV->Tolerance(...)`.
    if old.tolerance > new_core.tolerance {
        new_core.tolerance = old.tolerance;
    }
    // `NewTV->ChangePoints().Append(TV->ChangePoints())`.
    new_core.points.extend(old.points.iter().cloned());
}

/// The `for (j = 1; j <= ContextFaces/LF.Extent(); j++)` tail: for every image
/// face of `shapes`, walk its edges, and for every direct child vertex of those
/// edges that `map` binds, substitute the vertex, merge its geometry, rebind it
/// in `AsDes` and rewrite it on the edge.
///
/// `shapes` is `LF` for `Build`/`BuildFaces` and `LContext` for
/// `BuildOnContext`; the reference repeats this block verbatim in all three.
pub fn substitute_vertices(
    shapes: &[TopoShape],
    image: &BRepAlgoImage,
    as_des: &mut AsDes,
    map: &HashMap<ShapeMapKey, TopoShape>,
) {
    let builder = TopoBuilder::new();
    for f in shapes {
        // `Image.LastImage(F, LIF)` (`:141`, `:412`, `:624`).
        let mut lif: Vec<TopoShape> = Vec::new();
        image.last_image(f, &mut lif);
        for if_shape in &lif {
            // `TopExp_Explorer EdExp(IF, TopAbs_EDGE)`.
            for e_shape in explorer_shapes(if_shape, ShapeType::Edge) {
                // `TopoDS_Shape E = EdExp.Current();` - a copy sharing the
                // `TShape`, so the `Remove` / `Add` below rewrite the edge.
                let mut e = e_shape;
                // `TopoDS_Iterator VerExp(E)` (`CumOri = false`).
                let vlist = children_keep_orientation(&e);
                for v in &vlist {
                    // `if (myVerVerMap.IsBound(V))`.
                    let Some(new_v) = map.get(&shape_map_key(v)) else {
                        continue;
                    };
                    let mut new_v = new_v.clone();
                    // `E.Free(true)`.
                    e.set_free(true);
                    // `NewV.Orientation(V.Orientation())`.
                    new_v.set_orientation(v.orientation());
                    merge_vertex_geometry(v, &new_v);
                    // `AsDes->Replace(V, NewV)`.
                    as_des.replace(v, &new_v);
                    // `BB.Remove(E, V); BB.Add(E, NewV)`.
                    builder.remove(&mut e, v);
                    builder.add(&mut e, &new_v);
                }
            }
        }
    }
}
