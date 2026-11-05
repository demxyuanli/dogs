//! `BRepFill_OffsetWire::UpdateDetromp` (`cxx:1290-1402`): ranks the vertices
//! that limit the parts of the bisector between the spine and the offset, so
//! `TrimEdge` can drop the ones coming from tangency zones.

use std::collections::HashMap;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;

use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_fill_trim_edge_tool::TrimEdgeTool;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::{shape_map_key, ShapeMapKey};

impl BRepFillOffsetWire {
    /// `UpdateDetromp(Detromp, Shape1, Shape2, Vertices, Params, Bisec, SOnE,
    /// EOnE, Trim)` (`cxx:1302-1402`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_detromp(
        &self,
        detromp: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        shape1: &TopoShape,
        shape2: &TopoShape,
        vertices: &[Vertex],
        params: &[GpPnt],
        bisec: &occt_geom2d::bisector::BisectorBisec,
        son_e: bool,
        eon_e: bool,
        trim: &TrimEdgeTool,
    ) {
        let k1 = shape_map_key(shape1);
        let k2 = shape_map_key(shape2);
        let push = |v: &Vertex, detromp: &mut HashMap<ShapeMapKey, Vec<TopoShape>>| {
            detromp.entry(k1.clone()).or_default().push(v.0.clone());
            detromp.entry(k2.clone()).or_default().push(v.0.clone());
        };

        // `cxx:1315-1324`.
        if self.my_join_type == JoinType::Intersection {
            for v in vertices {
                push(v, detromp);
            }
            return;
        }

        // `cxx:1326-1338`: `ForceAdd` for a closed periodic bisector.
        let Some(bis) = bisec.value() else {
            return;
        };
        let bis: &dyn occt_geom2d::curve::Curve2d = bis;
        let mut force_add = false;
        if bis.is_periodic() {
            let p_f = bis.d0(bis.first_parameter());
            let p_l = bis.d0(bis.last_parameter());
            force_add = p_f.distance(&p_l) <= CONFUSION;
        }

        let mut u1 = bis.first_parameter();
        let mut ii = 0usize;
        let mut v1: Option<&Vertex> = None;

        // `cxx:1341-1345`.
        if son_e && !vertices.is_empty() {
            v1 = Some(&vertices[0]);
            ii = 1;
        }

        // `cxx:1347-1367`.
        while ii < vertices.len() && ii < params.len() {
            let u2 = params[ii].x();
            let v2 = &vertices[ii];
            let p = bis.d0((u2 + u1) * 0.5);
            if !trim.is_inside(&p) || force_add {
                if let Some(v) = v1 {
                    push(v, detromp);
                }
                push(v2, detromp);
            }
            u1 = u2;
            v1 = Some(v2);
            ii += 1;
        }

        // `cxx:1369-1395`.
        let u2 = bis.last_parameter();
        if !eon_e {
            if !u2.is_infinite() {
                let p = bis.d0((u2 + u1) * 0.5);
                if (!trim.is_inside(&p) || force_add) && v1.is_some() {
                    push(v1.unwrap(), detromp);
                }
            } else if let Some(v) = v1 {
                push(v, detromp);
            }
        }
    }
}
