//! BOPAlgo_WireSplitter::SplitBlock (BOPAlgo_WireSplitter_1.cxx).
//!
//! Irregular connexity blocks (a vertex with other than one IN and one OUT,
//! or TShape-coincident edges) are walked by Path: unused OUT edges, the
//! smallest ClockWiseAngle from the incoming IN angle, then MakeWire of
//! each closed buffer. Regular 1-in/1-out blocks without TShape coincidence
//! fall through to a single MakeWire.

mod prelude {
pub(super) use std::collections::HashMap;
pub(super) use std::f64::consts::PI;

pub(super) use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
pub(super) use occt_core::precision::{Precision, ANGULAR, PCONFUSION};
pub(super) use occt_geom2d::cl_props2d::GeomLPropCLProps2d;
pub(super) use occt_geom2d::line::Geom2dLine;

pub(super) use crate::abs::{Orientation, ShapeType};
pub(super) use crate::boptools_2d::make_2d;
pub(super) use crate::brep_tool::BRepTool;
pub(super) use crate::builder::TopoBuilder;
pub(super) use crate::shape::{Edge, Face, TopoShape, Vertex};
pub(super) use crate::tgeometry::GeometryRegistry;
pub(super) use crate::topo_tools::{is_equal, is_same};

pub(super) const TWO_PI: f64 = PI + PI;
}

use prelude::*;

mod angles;
mod path;
mod refine;
mod vertices;

use angles::*;
use path::*;
use refine::*;
use vertices::*;

/// `BOPAlgo_EdgeInfo`.
#[derive(Clone)]
pub(super) struct EdgeInfo {
    pub(super) edge: Edge,
    pub(super) passed: bool,
    pub(super) is_in: bool,
    pub(super) is_inside: bool,
    pub(super) angle: f64,
}

pub(super) struct SmartMap {
    pub(super) verts: Vec<TopoShape>,
    pub(super) infos: Vec<Vec<EdgeInfo>>,
    pub(super) index: HashMap<usize, usize>,
}

impl SmartMap {
    pub(super) fn new() -> Self {
        Self {
            verts: Vec::new(),
            infos: Vec::new(),
            index: HashMap::new(),
        }
    }

    pub(super) fn slot(&mut self, v: &TopoShape) -> usize {
        let k = GeometryRegistry::shape_key(v);
        if let Some(&i) = self.index.get(&k) {
            return i;
        }
        let i = self.verts.len();
        self.index.insert(k, i);
        self.verts.push(v.clone());
        self.infos.push(Vec::new());
        i
    }
}

/// `BOPAlgo_WireSplitter::SplitBlock`.
pub(crate) fn split_wire_block(face: &Face, edges: &[Edge]) -> Result<Vec<TopoShape>, String> {
    let mut map = SmartMap::new();
    let mut unique: HashMap<usize, usize> = HashMap::new();
    let mut vert_closed: HashMap<usize, bool> = HashMap::new();
    let fk = GeometryRegistry::shape_key(&face.0);

    for e in edges {
        if make_2d(e, face).is_err() {
            continue;
        }
        let ek = GeometryRegistry::shape_key(&e.0);
        let closed = edge_closed_on_face(e, face, fk);
        let n = unique.entry(ek).or_insert(0);
        *n += 1;
        if *n > 1 && !closed {
            unique.insert(ek, 0);
        }
        let kids = oriented_vertices(e);
        for (i, v) in kids.iter().enumerate() {
            let slot = map.slot(v);
            let is_in = v.orientation() == Orientation::Reversed;
            map.infos[slot].push(EdgeInfo {
                edge: e.clone(),
                passed: false,
                is_in,
                is_inside: false,
                angle: -1.0,
            });
            let vk = GeometryRegistry::shape_key(v);
            let cyc = closed || (kids.len() >= 2 && is_same(&kids[0], &kids[1]));
            let e0 = vert_closed.entry(vk).or_insert(false);
            if cyc || i > 0 && is_same(&kids[0], v) {
                *e0 = true;
            }
            if cyc {
                *e0 = true;
            }
        }
    }

    if map.verts.is_empty() {
        // `BOPAlgo_WireSplitter_1.cxx:144-148`: start elements without a
        // p-curve on the face are skipped, so nothing enters the smart map and
        // `bNothingToDo` (`cxx:199-283`) stays true -> a plain wire over the
        // whole block (`cxx:284-295`).
        return Ok(vec![make_wire(edges)]);
    }

    for infos in &mut map.infos {
        for ei in infos.iter_mut() {
            let ek = GeometryRegistry::shape_key(&ei.edge.0);
            ei.is_inside = unique.get(&ek).copied().unwrap_or(0) == 0;
        }
    }

    if nothing_to_do(&map, edges) {
        return Ok(vec![make_wire(edges)]);
    }

    for slot in 0..map.verts.len() {
        let v = map.verts[slot].clone();
        let n = map.infos[slot].len();
        for i in 0..n {
            let e = map.infos[slot][i].edge.clone();
            let is_in = map.infos[slot][i].is_in;
            let mut vv = v.clone();
            vv.set_orientation(if is_in {
                Orientation::Reversed
            } else {
                Orientation::Forward
            });
            map.infos[slot][i].angle = angle_2d(&vv, &e, face, is_in);
        }
    }
    refine_angles(face, &mut map);

    if std::env::var("ZZ_SB").is_ok() {
        dump_map(face, &map, &vert_closed);
    }

    let mut loops: Vec<TopoShape> = Vec::new();
    let nverts = map.verts.len();
    for slot in 0..nverts {
        let ninfo = map.infos[slot].len();
        for i in 0..ninfo {
            if map.infos[slot][i].is_in || map.infos[slot][i].passed {
                continue;
            }
            path(
                face,
                &vert_closed,
                slot,
                i,
                &mut map,
                &mut loops,
            );
        }
    }
    if loops.is_empty() {
        // `Path` (`BOPAlgo_WireSplitter_1.cxx:407-611`) can end without closing
        // a single loop on a block that is neither regular nor walkable; OCCT
        // only leaves those edges unrepresented (they go to `myShapesToAvoid`
        // in `BOPAlgo_BuilderFace::PerformLoops`, `BOPAlgo_BuilderFace.cxx:313-321`),
        // it is not an error.
        return Ok(Vec::new());
    }
    Ok(loops)
}

