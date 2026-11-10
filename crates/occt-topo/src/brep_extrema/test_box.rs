//! Shared unit-box fixture for the Phase 3 test modules.
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_geom::GeomLine;

use crate::builder::TopoBuilder;
use crate::abs::Orientation;
use crate::shape::{Edge, Face, Solid, Vertex};

pub(crate) struct UnitBox {
    pub solid: Solid,
    pub corners: [GpPnt; 8],
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
}

fn segment_edge(
    b: &TopoBuilder,
    p1: &GpPnt,
    p2: &GpPnt,
    v1: &Vertex,
    v2: &Vertex,
) -> Edge {
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).expect("distinct points");
    let lin = GpLin::from_pnt_dir(*p1, dir);
    let mut e = b.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
    b.add_edge_vertices(&mut e, v1, v2);
    e
}

/// Build a closed unit box solid spanning `[0,1]³`.
///
/// 8 shared vertices, 12 segment edges (each referencing two shared
/// vertices), 6 planar faces (each with a closed wire of 4 edges), one
/// shell, one solid. The face planes are oriented so their `(u,v)` grid
/// maps `[0,1]²` onto each face quad with a consistent (inward) normal,
/// which keeps the mesh signed-volume and centroid computations coherent.
pub(crate) fn unit_box() -> UnitBox {
    let b = TopoBuilder::new();
    let corners = [
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(1.0, 0.0, 0.0),
        GpPnt::new(1.0, 1.0, 0.0),
        GpPnt::new(0.0, 1.0, 0.0),
        GpPnt::new(0.0, 0.0, 1.0),
        GpPnt::new(1.0, 0.0, 1.0),
        GpPnt::new(1.0, 1.0, 1.0),
        GpPnt::new(0.0, 1.0, 1.0),
    ];
    let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();

    let edge_idx: [(usize, usize); 12] = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ];
    let mut edges = Vec::new();
    for &(i, j) in &edge_idx {
        edges.push(segment_edge(&b, &corners[i], &corners[j], &vertices[i], &vertices[j]));
    }

    // Each face's 4 edges (indices into `edges`), sharing edges across faces.
    let face_edge_sets: [[usize; 4]; 6] = [
        [0, 1, 2, 3],     // bottom (z = 0)
        [4, 5, 6, 7],     // top (z = 1)
        [0, 9, 4, 8],     // front (y = 0)
        [2, 10, 6, 11],   // back (y = 1)
        [3, 11, 7, 8],    // left (x = 0)
        [1, 10, 5, 9],    // right (x = 1)
    ];
    // (origin, outward normal, u-direction) per face.
    let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
        (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        (GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
        (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
        (GpPnt::new(0.0, 1.0, 0.0), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
        (GpPnt::new(0.0, 0.0, 0.0), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
        (GpPnt::new(1.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
    ];
    // Outward corner cycle of every face, in the same order as
    // `BRepPrim_GWedge` / `BRepPrimBox::make_box` (`primitives/brep_prim_box.rs`). An
    // edge that the face closes against the `edge_idx` chord is stored
    // Reversed, so two faces sharing an edge always see opposite
    // orientations — the invariant `BOPTools_AlgoTools::GetEdgeOff`
    // (`BOPTools_AlgoTools.cxx:1099-1127`) and
    // `BOPAlgo_ShellSplitter::SplitBlock` (`BOPAlgo_ShellSplitter.cxx:319`)
    // rely on.
    let face_cycles: [[usize; 4]; 6] = [
        [0, 1, 2, 3], // bottom (z = 0), reversed below to face -Z
        [4, 5, 6, 7], // top (z = 1)
        [0, 1, 5, 4], // front (y = 0)
        [3, 7, 6, 2], // back (y = 1)
        [0, 4, 7, 3], // left (x = 0)
        [1, 2, 6, 5], // right (x = 1)
    ];
    let mut faces = Vec::new();
    for fi in 0..6 {
        let (origin, normal, u_dir) = face_planes[fi];
        let ax3 = GpAx3::new(origin, normal, &u_dir).expect("perpendicular axes");
        let mut face = b.make_face_plane(&GpPln::new(ax3));
        let cycle = face_cycles[fi];
        let quad: Vec<Edge> = face_edge_sets[fi]
            .iter()
            .map(|&ei| {
                let (i, j) = edge_idx[ei];
                // Directed chord of the outward cycle closed by this edge.
                let k = (0..4)
                    .find(|&k| {
                        let a = cycle[k];
                        let b = cycle[(k + 1) % 4];
                        (a == i && b == j) || (a == j && b == i)
                    })
                    .expect("edge belongs to the face cycle");
                let (a, b) = (cycle[k], cycle[(k + 1) % 4]);
                let e = edges[ei].clone();
                if (a, b) == (i, j) {
                    e
                } else {
                    Edge(e.0.oriented(Orientation::Reversed))
                }
            })
            .collect();
        let mut wire = b.make_wire(&quad);
        // Bottom plane: origin (0,0,0), normal -Z, u = +Y so v = +X.
        // The +X then +Y edge loop is clockwise in that UV frame; reverse
        // the wire so the outer ring is CCW (material on the left), matching
        // BRepPrimAPI_MakeBox / IntTools_FClass2d::IsHole.
        if fi == 0 {
            wire.0.reverse();
        }
        b.add_wire(&mut face, &wire);
        faces.push(face);
    }
    let shell = b.make_shell(&faces);
    let solid = b.make_solid(&[shell]);
    UnitBox { solid, corners, vertices, edges, faces }
}
