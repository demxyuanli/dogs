//! `ComputeIntervals` (`BRepFill_Evolved.cxx:2866-2946`) and `VertexFromNode`
//! (`BRepFill_Evolved.cxx:3391-3443`).
//!
//! `ComputeIntervals` merges the two vertex sequences collected on the two
//! branches of a bisecting-locus edge into the valid intervals of the bisector:
//! an interval survives only when the midpoint of its bisector arc lies on the
//! trimmed face (`BRepFill_TrimSurfaceTool::IsOnFace`).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom2d::curve::Curve2d;

use crate::brep_fill_trim_surface_tool::TrimSurfaceTool;
use crate::mat::{MatNode, NodeId};
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::is_same;

use super::shape_map::ShapeMap;
use super::side::{distance_to_oz, pos_on_face};

/// `MapNodeVertex` (`BRepFill_Evolved.cxx:685-687`):
/// `NCollection_DataMap<handle<MAT_Node>, NCollection_DataMap<Shape, Shape>>`.
/// The node handle key is its arena slot, exactly as in
/// `brep_fill_offset_wire::vertices`.
pub type MapNodeVertex = HashMap<NodeId, ShapeMap<TopoShape>>;

/// The four output sequences of `ComputeIntervals` (`cxx:2878-2882`).
#[derive(Default)]
pub struct Intervals {
    pub first_par: Vec<f64>,
    pub last_par: Vec<f64>,
    pub first_v: Vec<TopoShape>,
    pub last_v: Vec<TopoShape>,
}

/// `ComputeIntervals(VOnF, VOnL, ParOnF, ParOnL, Trim, Bis, VS, VE, FirstPar,
/// LastPar, FirstV, LastV)` (`BRepFill_Evolved.cxx:2866-2946`).
///
/// `ParOnF` / `ParOnL` are `gp_Pnt` sequences and only the `X` component is
/// read (`cxx:2898`, `cxx:2906`). `vs` / `ve` are the optional `VS` / `VE`
/// parameters; a null `TopoDS_Vertex` maps to `None`.
pub fn compute_intervals(
    v_on_f: &[TopoShape],
    v_on_l: &[TopoShape],
    par_on_f: &[GpPnt],
    par_on_l: &[GpPnt],
    trim: &TrimSurfaceTool,
    bis: &Arc<dyn Curve2d>,
    vs: Option<&TopoShape>,
    ve: Option<&TopoShape>,
) -> Intervals {
    let mut out = Intervals::default();

    // `cxx:2883-2889`.
    let mut i_on_f = 0usize;
    let mut i_on_l = 0usize;
    let mut u1 = 0.0f64;
    let mut v1: Option<&TopoShape> = None;
    if let Some(vs) = vs {
        u1 = bis.first_parameter();
        v1 = Some(vs);
    }
    let mut u2;
    let mut v2: Option<&TopoShape>;

    // `cxx:2890-2928`.
    while i_on_f < v_on_f.len() || i_on_l < v_on_l.len() {
        if i_on_l >= v_on_l.len()
            || (i_on_f < v_on_f.len() && par_on_f[i_on_f].x() < par_on_l[i_on_l].x())
        {
            u2 = par_on_f[i_on_f].x();
            v2 = Some(&v_on_f[i_on_f]);
            i_on_f += 1;
        } else {
            u2 = par_on_l[i_on_l].x();
            v2 = Some(&v_on_l[i_on_l]);
            i_on_l += 1;
        }

        // `cxx:2914-2926`: `!V1.IsNull() && !V2.IsSame(V1)`.
        if let (Some(v1), Some(v2)) = (v1, v2) {
            if !is_same(v1, v2) {
                let p = bis.value((u2 + u1) * 0.5);
                if trim.is_on_face(&p) {
                    out.first_par.push(u1);
                    out.last_par.push(u2);
                    out.first_v.push(v1.clone());
                    out.last_v.push(v2.clone());
                }
            }
        }
        u1 = u2;
        v1 = v2;
    }

    // `cxx:2930-2945`: the closing interval up to `VE`.
    if let Some(ve) = ve {
        u2 = bis.last_parameter();
        v2 = Some(ve);
        if let Some(v1) = v1 {
            if !is_same(v1, v2.unwrap()) {
                let p = bis.value((u2 + u1) * 0.5);
                if trim.is_on_face(&p) {
                    out.first_par.push(u1);
                    out.last_par.push(u2);
                    out.first_v.push(v1.clone());
                    out.last_v.push(v2.unwrap().clone());
                }
            }
        }
    }

    out
}

/// `VertexFromNode(aNode, E, VF, VL, MapNodeVertex, VN)`
/// (`BRepFill_Evolved.cxx:3391-3443`).
///
/// `Status` is 1 on `VF`, 2 in the middle of `E` and 3 on `VL`; 0 means the node
/// was not usable, and `None` answers the null `VN` of that branch. The created
/// vertex is shared per node and per shape it sits on, like the two-level
/// `NCollection_DataMap` of the reference.
pub fn vertex_from_node(
    node: NodeId,
    mat_node: &MatNode,
    e: &TopoShape,
    vf: &Vertex,
    vl: &Vertex,
    map_node_vertex: &mut MapNodeVertex,
) -> (i32, Option<Vertex>) {
    // `cxx:3406-3409`.
    let mut status = 0;
    if !mat_node.infinite() {
        status = pos_on_face(mat_node.distance(), distance_to_oz(vf), distance_to_oz(vl));
    }

    // `cxx:3410-3421`.
    let shape_on_node: Option<&TopoShape> = match status {
        2 => Some(e),
        1 => Some(&vf.0),
        3 => Some(&vl.0),
        _ => None,
    };
    let Some(shape_on_node) = shape_on_node else {
        return (status, None);
    };

    // `cxx:3428-3431`.
    if let Some(map) = map_node_vertex.get(&node) {
        if let Some(found) = map.find(shape_on_node) {
            return (status, Some(Vertex(found.clone())));
        }
    }

    // `cxx:3432-3440`: `B.MakeVertex(VN)` then bind both levels.
    let vn = Vertex::new();
    map_node_vertex
        .entry(node)
        .or_insert_with(ShapeMap::new)
        .bind(shape_on_node, vn.0.clone());
    (status, Some(vn))
}
