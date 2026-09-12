use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// 3. Shell / solid offset (BRepOffsetAPI_MakeOffsetShape)
// ---------------------------------------------------------------------------

/// Offset every vertex of a convex planar-face polyhedron by intersecting the
/// offset planes of the incident faces, then rebuild the boundary tree with
/// shared edges so the result is a closed shell.
pub(super) fn offset_convex_polyhedron(shape: &TopoShape, faces: &[Face], distance: f64) -> Result<TopoShape, String> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return Err("offset_convex_polyhedron: no vertices".into());
    }
    let mut acc = GpXyz::zero();
    for v in &verts {
        acc = acc.added(&BRepTool::vertex_point(v).coord);
    }
    let centroid = GpPnt::from_xyz(&acc.divided(verts.len() as f64));

    pub(super) struct FaceOff {
        pub(super) pln: GpPln,
        pub(super) offset_pln: GpPln,
        pub(super) normal: GpVec,
    }
    let mut face_offs = Vec::with_capacity(faces.len());
    for f in faces {
        let pln = face_plane(f).ok_or("offset_convex_polyhedron: non-planar face")?;
        let n0 = *pln.axis().direction();
        let probe = GpVec::from_pnts(&centroid, &pln.location());
        let n = if n0.xyz().dot(&probe.xyz()) >= 0.0 { n0 } else { n0.reversed() };
        let nv = GpVec::from_xyz(n.xyz());
        for v in &verts {
            let p = BRepTool::vertex_point(v);
            let d = GpVec::from_pnts(&pln.location(), &p).xyz().dot(&n.xyz());
            if d > 1e-6 {
                return Err("offset_convex_polyhedron: shape is not convex".into());
            }
        }
        let offset_pln = pln.translated_vec(&nv.multiplied_scalar(distance));
        face_offs.push(FaceOff { pln, offset_pln, normal: nv });
    }

    let b = TopoBuilder::new();
    let mut new_verts: Vec<(GpPnt, Vertex, GpPnt)> = Vec::with_capacity(verts.len());
    for v in &verts {
        let p = BRepTool::vertex_point(v);
        let mut incident: Vec<&GpPln> = Vec::new();
        for fo in &face_offs {
            let d = GpVec::from_pnts(&fo.pln.location(), &p).xyz().dot(&fo.normal.xyz());
            if d.abs() < 1e-6 {
                incident.push(&fo.offset_pln);
            }
        }
        if incident.len() < 3 {
            return Err("offset_convex_polyhedron: vertex has fewer than 3 incident planes".into());
        }
        let np = plane_plane_plane_intersection(incident[0], incident[1], incident[2])
            .ok_or("offset_convex_polyhedron: parallel offset planes")?;
        new_verts.push((p, b.make_vertex(np, 0.0), np));
    }
    let find = |p: &GpPnt| new_verts.iter().find(|(orig, _, _)| orig.distance(p) < 1e-6);

    let edges = edges_of(shape);
    let mut emap: HashMap<usize, Edge> = HashMap::new();
    for e in &edges {
        let (a, z) = crate::topo_tools_full::edge_vertices(e);
        let (Some(va), Some(vb)) = (a, z) else {
            return Err("offset_convex_polyhedron: edge missing endpoint vertices".into());
        };
        let p1 = BRepTool::vertex_point(&va);
        let p2 = BRepTool::vertex_point(&vb);
        let (_, _, np1) = find(&p1).ok_or("offset_convex_polyhedron: endpoint not matched")?;
        let (_, _, np2) = find(&p2).ok_or("offset_convex_polyhedron: endpoint not matched")?;
        let seg = b.make_edge_segment(&np1, &np2);
        emap.insert(std::sync::Arc::as_ptr(&e.tshape) as usize, seg);
    }

    let mut new_faces = Vec::with_capacity(faces.len());
    for (f, fo) in faces.iter().zip(&face_offs) {
        let wires = wires_of_face(f);
        let mut new_wires = Vec::with_capacity(wires.len());
        for w in wires {
            let w_edges = edges_of_wire(&w);
            let mapped: Vec<Edge> = w_edges
                .iter()
                .map(|e| emap.get(&(std::sync::Arc::as_ptr(&e.tshape) as usize)).cloned())
                .collect::<Option<_>>()
                .ok_or("offset_convex_polyhedron: face edge not mapped")?;
            let nw = b.make_wire(&mapped);
            nw.set_closed(true);
            new_wires.push(nw);
        }
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(fo.offset_pln.clone()));
        new_faces.push(b.make_face(surface, &new_wires));
    }

    let shell = b.make_shell(&new_faces);
    if shape.shape_type() == ShapeType::Solid {
        Ok(b.make_solid(&[shell]).into())
    } else {
        Ok(shell.into())
    }
}

/// Offset a shell or solid. Convex planar-face shapes (boxes, prisms) are
/// rebuilt by intersecting the offset face planes; a single curved face (a
/// sphere) offsets in place. Other shapes return an error rather than a
/// broken result.
pub fn offset_shell(shape: &TopoShape, distance: f64) -> Result<TopoShape, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("offset_shell: no faces".into());
    }
    if faces.iter().all(face_is_planar) {
        if let Ok(s) = offset_convex_polyhedron(shape, &faces, distance) {
            return Ok(s);
        }
    }
    if faces.len() == 1 {
        let f = offset_face(&faces[0], distance)?;
        let b = TopoBuilder::new();
        let shell = b.make_shell(&[f]);
        if shape.shape_type() == ShapeType::Solid {
            return Ok(b.make_solid(&[shell]).into());
        }
        return Ok(shell.into());
    }
    Err("offset_shell: non-convex or curved offset unsupported".into())
}
