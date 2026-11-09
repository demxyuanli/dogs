use super::prelude::*;


/// Wireframe edges of a box: pairs of corner indices into `box_corners`.
pub(super) const BOX_EDGES: [(usize, usize); 12] = [
    (0, 1), (1, 2), (2, 3), (3, 0), // bottom face
    (4, 5), (5, 6), (6, 7), (7, 4), // top face
    (0, 4), (1, 5), (2, 6), (3, 7), // vertical edges
];


/// Unit direction (panics on a degenerate axis — callers use valid axes).
pub(super) fn dir(x: f64, y: f64, z: f64) -> GpDir {
    GpDir::new(x, y, z).expect("primitives: degenerate direction")
}


/// Line curve from `p1` to `p2`, parameterized on [0, |p2 − p1|].
pub(super) fn line_curve(p1: &GpPnt, p2: &GpPnt) -> Arc<dyn Curve> {
    let d = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).unwrap_or_else(|_| dir(1.0, 0.0, 0.0));
    Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p1, d)))
}


/// Full circle curve in the plane `ax3` with the given radius (parameter 0..2π).
pub(super) fn circle_curve(ax3: &GpAx3, radius: f64) -> Arc<dyn Curve> {
    let ax2 = GpAx2::new(ax3.location(), ax3.direction(), *ax3.x_direction())
        .expect("primitives: circle axis is not orthonormal");
    Arc::new(GeomCircle::new(GpCirc::new(ax2, radius)))
}


/// Plane through `origin` with the given outward `normal`; the in-plane X axis
/// is auto-selected perpendicular to the normal.
pub(super) fn plane(origin: GpPnt, normal: GpDir) -> GpPln {
    let z = normal;
    let z_axis = dir(0.0, 0.0, 1.0);
    let x_dir = if z.is_normal(&z_axis) { z_axis } else { dir(1.0, 0.0, 0.0) };
    GpPln::new(GpAx3::new(origin, z, &x_dir).expect("primitives: non-orthogonal plane frame"))
}


/// Map a corner pair to the index of the box edge between them.
pub(super) fn edge_index(a: usize, b: usize) -> usize {
    BOX_EDGES
        .iter()
        .position(|&(i, j)| (i == a && j == b) || (i == b && j == a))
        .unwrap_or_else(|| panic!("primitives: box edge ({a},{b}) not in BOX_EDGES"))
}


/// The 8 corners of an axis-aligned box from the origin, index order:
/// 0=(0,0,0), 1=(dx,0,0), 2=(dx,dy,0), 3=(0,dy,0),
/// 4=(0,0,dz), 5=(dx,0,dz), 6=(dx,dy,dz), 7=(0,dy,dz).
pub(super) fn box_corners(dx: f64, dy: f64, dz: f64) -> [GpPnt; 8] {
    [
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(dx, 0.0, 0.0),
        GpPnt::new(dx, dy, 0.0),
        GpPnt::new(0.0, dy, 0.0),
        GpPnt::new(0.0, 0.0, dz),
        GpPnt::new(dx, 0.0, dz),
        GpPnt::new(dx, dy, dz),
        GpPnt::new(0.0, dy, dz),
    ]
}


/// Build a real-geometry box solid. Faces are planar, edges are lines, the
/// 8 corners are shared by all incident edges (OCCT's shared-TShape model).
pub(super) fn build_box(b: &TopoBuilder, corners: &[GpPnt; 8]) -> Solid {
    // `BRepPrim_Builder::MakeVertex` (`BRepPrim_Builder.cxx:129-132`) is
    // `BRep_Builder::MakeVertex(V, P, Precision::Confusion())`, so a primitive
    // vertex carries the 1e-7 tolerance. It is load-bearing:
    // `BRepTools_WireExplorer::Init` (`BRepTools_WireExplorer.cxx:139-165`)
    // returns an empty walk when every vertex *and* edge tolerance of the wire
    // is below `Precision::Confusion()`.
    let verts: [Vertex; 8] =
        core::array::from_fn(|i| b.make_vertex(corners[i], occt_core::precision::CONFUSION));

    // 12 edges, each carrying a line curve and the two shared corner vertices.
    let mut edges: Vec<Edge> = Vec::with_capacity(12);
    for &(i, j) in &BOX_EDGES {
        let p1 = &corners[i];
        let p2 = &corners[j];
        let mut e = b.make_edge(line_curve(p1, p2), 0.0, p1.distance(p2));
        b.add_edge_vertices(&mut e, &verts[i], &verts[j]);
        edges.push(e);
    }

    // 6 faces. Each face is a plane through 4 corners; its wire is the cyclic
    // list of the 4 boundary edges (each shared with a neighboring face).
    let faces_def: [(GpPnt, GpDir, [usize; 4]); 6] = [
        // (a corner point on the face, outward normal, corner cycle CCW-outward)
        (corners[0], dir(0.0, 0.0, -1.0), [0, 3, 2, 1]), // -Z
        (corners[4], dir(0.0, 0.0, 1.0), [4, 5, 6, 7]), // +Z
        (corners[0], dir(0.0, -1.0, 0.0), [0, 1, 5, 4]), // -Y
        (corners[3], dir(0.0, 1.0, 0.0), [3, 7, 6, 2]), // +Y
        (corners[0], dir(-1.0, 0.0, 0.0), [0, 4, 7, 3]), // -X
        (corners[1], dir(1.0, 0.0, 0.0), [1, 2, 6, 5]), // +X
    ];

    let mut faces: Vec<Face> = Vec::with_capacity(6);
    for (origin, normal, cycle) in faces_def {
        // Orient each shared edge to the face's CCW-outward traversal: the edge
        // TShape carries one fixed curve direction (BOX_EDGES), so a wire that
        // traverses it the other way stores it Reversed (the OCCT
        // `BRepPrim_GWedge` / `BRepBuilderAPI_MakeWire` edge-orientation model).
        let quad: Vec<Edge> = [0, 1, 2, 3]
            .iter()
            .map(|&k| {
                let a = cycle[k];
                let b = cycle[(k + 1) % 4];
                let idx = edge_index(a, b);
                let e = edges[idx].clone();
                if BOX_EDGES[idx] == (a, b) {
                    e
                } else {
                    Edge(e.0.oriented(Orientation::Reversed))
                }
            })
            .collect();
        let wire = b.make_wire(&quad);
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane(origin, normal)));
        faces.push(b.make_face(surface, &[wire]));
    }

    let shell = b.make_shell(&faces);
    b.make_solid(&[shell])
}
