use super::prelude::*;


/// `BRepFill::Axe` (`BRepFill.cxx:672-878`).
///
/// Computes the coordinate system `AxeProf` of an evolved sweep: the normal of
/// the spine plane, and the point / tangent of the spine closest to the
/// profile. `ProfOnSpine` reports whether that closest point is within `tol`
/// (`cxx:874`), which the caller uses to reject an inconsistent
/// `ProfOnSpine` request.
///
/// Translation notes:
/// - `BRep_Tool::Surface(Face, L)` + `S->Transformed(L)` (`cxx:690-719`) is
///   [`BRepTool::face_surface`] plus an explicit composition of the
///   `BRepLib_FindSurface` location with the face location, the same way
///   [`crate::brep_fill_evolved::perform::FindLocation`] does it.
/// - `BRepExtrema_ExtPC::Perform(V)` + `IsMin(i)` (`cxx:845-857`) is
///   [`crate::brep_extrema::closest_point_on_edge`], which returns the
///   minimal solution of the point/curve extrema.
pub fn axe(spine: &TopoShape, profile: &Wire, tol: f64) -> Result<(GpAx3, bool), String> {
    // `cxx:684-709`: the face and its plane.
    let mut a_face = Face(TopoShape::new(ShapeType::Face));
    let mut surface: Option<Arc<dyn Surface>> = None;
    let mut location = TopLocLocation::identity();
    if spine.shape_type() == ShapeType::Face {
        let face = Face(spine.clone());
        location = face.0.location().clone();
        surface = BRepTool::face_surface(&face);
        // `if (!S->IsKind(STANDARD_TYPE(Geom_Plane)))` (`cxx:691`): the exact
        // type test, so a trimmed surface of a plane takes the
        // `BRepLib_FindSurface` branch.
        let is_plain_plane = surface
            .as_ref()
            .is_some_and(|s| s.rectangular_trimmed_basis().is_none() && s.gp_pln().is_some());
        if !is_plain_plane {
            // `cxx:693-702`.
            let fs = FindSurface::new(&face.0, -1.0, true, false);
            if !fs.is_found() {
                return Err(
                    "Standard_NoSuchObject: BRepFill_Evolved : The Face is not planar".into(),
                );
            }
            surface = fs.surface();
            location = TopLocLocation::composed(
                &fs.location().transformation(),
                face.0.location().clone(),
            );
        }
        a_face = face;
    } else if spine.shape_type() == ShapeType::Wire {
        // `cxx:705-709`.
        a_face = make_face_from_wire(&Wire(spine.clone()), true)
            .ok_or_else(|| "BRepLib_MakeFace: the wire is not planar".to_string())?;
        location = a_face.0.location().clone();
        surface = BRepTool::face_surface(&a_face);
    }

    // `cxx:711-714`.
    let mut surface = match surface {
        Some(s) => s,
        None => return Err("Standard_DomainError: BRepFill_Evolved::Axe".into()),
    };

    // `cxx:716-719`.
    if !location.is_identity() {
        surface = Arc::from(surface.transformed(&location.transformation()));
    }

    // `cxx:721`: `down_cast<Geom_Plane>(S)->Pln().Axis().Direction()`.
    let normal = surface
        .gp_pln()
        .map(|p| {
            let ax1 = p.axis();
            *ax1.direction()
        })
        .ok_or_else(|| "Standard_DomainError: BRepFill_Evolved::Axe".to_string())?;

    // `cxx:724-727`.
    let mut dist_min = INFINITE;
    // `double Tol2 = 1.e-10;` (`cxx:727`).
    let tol2 = 1.0e-10;

    // `cxx:734-760`: is there a vertex/vertex contact between the spine and the
    // profile? Note the loop does not stop at the first spine vertex in OCCT:
    // `SE.Next()` runs on every outer iteration, and the break only fires once
    // a contact is found (`cxx:739-760`).
    let forward_face = a_face.0.oriented(Orientation::Forward);
    let spine_vertices = explorer_shapes(&forward_face, ShapeType::Vertex);
    let profile_vertices = explorer_shapes(&profile.0, ShapeType::Vertex);

    let mut on_vertex: Option<Vertex> = None;
    for sv in &spine_vertices {
        let p1 = BRepTool::vertex_point_world(&Vertex(sv.clone()));
        let mut found = false;
        for pv in &profile_vertices {
            let p2 = BRepTool::vertex_point_world(&Vertex(pv.clone()));
            if p1.square_distance(&p2) <= tol2 {
                found = true;
                break;
            }
        }
        if found {
            on_vertex = Some(Vertex(sv.clone()));
            break;
        }
    }

    // `cxx:730`: `Par` only feeds the edge branch.
    let mut par = 0.0;
    let mut loc = GpPnt::new(0.0, 0.0, 0.0);
    let mut tang = GpVec::from_pnts(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 0.0, 0.0));

    if let Some(vonf) = on_vertex {
        // `cxx:762-833`.
        let map = map_shapes_and_ancestors(&forward_face, ShapeType::Vertex, ShapeType::Edge);
        let list = map
            .get(&shape_map_key(&vonf.0))
            .cloned()
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;
        let e1 = list
            .first()
            .map(|s| Edge(s.clone()))
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;
        let e2 = list
            .last()
            .map(|s| Edge(s.clone()))
            .ok_or_else(|| "Standard_NoSuchObject: TopExp_AncestorMap".to_string())?;

        // `CE1 = BRep_Tool::Curve(E1, L, f, l)` + `D1` (`cxx:778-789`).
        let (loc1, tang1) = edge_point_and_tangent(&e1, &vonf, &a_face)?;
        // `cxx:791-802`.
        let (loc2, tang2) = edge_point_and_tangent(&e2, &vonf, &a_face)?;

        // `cxx:804-820`: sum the |dot| of both tangents over the profile edges.
        let mut sca1 = 0.0;
        let mut sca2 = 0.0;
        for pe in explorer_shapes(&profile.0, ShapeType::Edge) {
            let e = Edge(pe);
            let (v1, v2) = edge_vertices(&e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                continue;
            };
            let p1 = BRepTool::vertex_point_world(&v1);
            let p2 = BRepTool::vertex_point_world(&v2);
            let vec = GpVec::from_pnts(&p1, &p2);
            sca1 += tang1.dot(&vec).abs();
            sca2 += tang2.dot(&vec).abs();
        }

        // `cxx:822-831`.
        if sca1.abs() < sca2.abs() {
            loc = loc1;
            tang = tang1;
        } else {
            loc = loc2;
            tang = tang2;
        }
        dist_min = 0.0;
    } else {
        // `cxx:834-872`.
        for se in explorer_shapes(&forward_face, ShapeType::Edge) {
            let e = Edge(se);
            for pv in &profile_vertices {
                let v = Vertex(pv.clone());
                let p2 = BRepTool::vertex_point_world(&v);
                // `BE.Perform(V)` (`cxx:845`): the minimal solution.
                let (u, _pt) = crate::brep_extrema::closest_point_on_edge(&e, &p2, 0);
                let dist = p2.distance(&evaluate_edge(&e, u)?);
                // `cxx:860-869`.
                if dist < dist_min {
                    dist_min = dist;
                    par = u;
                    // `BRepAdaptor_Curve BAC(E); BAC.D1(Par, Loc, Tang)`.
                    let curve = BRepTool::edge_curve_adaptor_world(&e)
                        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
                    let (p, mut t) = curve.d1(par);
                    if e.0.orientation() == Orientation::Reversed {
                        t.reverse();
                    }
                    loc = p;
                    tang = t;
                }
            }
        }
    }

    // `cxx:874`.
    let prof_on_spine = dist_min < tol;

    // `cxx:876-877`: `gp_Ax3 A3(Loc, Normal, Tang)`.
    let tang_dir = GpDir::from_vec(&tang)
        .map_err(|_| "Standard_ConstructionError: gp_Ax3".to_string())?;
    let axe_prof = GpAx3::new(loc, normal, &tang_dir)
        .map_err(|_| "Standard_ConstructionError: gp_Ax3".to_string())?;
    Ok((axe_prof, prof_on_spine))
}


/// Evaluates the 3D curve of `e` at `u` in world coordinates
/// (`BRepAdaptor_Curve::D0`).
fn evaluate_edge(e: &Edge, u: f64) -> Result<GpPnt, String> {
    let curve = BRepTool::edge_curve_world(e)
        .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
    Ok(curve.d0(u))
}


/// `CE1 = BRep_Tool::Curve(E, L, f, l); CE1->D1(Par, Loc, Tang)` with the
/// location applied and the tangent reversed for a `REVERSED` edge
/// (`BRepFill.cxx:778-789`, `:791-802`).
fn edge_point_and_tangent(
    e: &Edge,
    v: &Vertex,
    face: &Face,
) -> Result<(GpPnt, GpVec), String> {
    let curve = BRepTool::edge_curve_world(e)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let par = BRepTool::parameter_of_vertex_on_face(&v.0, e, face);
    let (p, mut t) = curve.d1(par);
    if e.0.orientation() == Orientation::Reversed {
        t.reverse();
    }
    // `Tang1.Normalize()` (`cxx:805-806`).
    t.normalize();
    Ok((p, t))
}

