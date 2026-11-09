//! Free-standing helpers of `BRepFill_TrimShellCorner` (see `mod.rs`).
//!
//! Moved here verbatim from the single-file port; the `BRepFill_TrimShellCorner.cxx`
//! branch structure is unchanged.

use super::*;

/// `myHistMap.Find(shape)` restricted to `lst.IsEmpty()` being false
/// (`cxx:479-497`, `cxx:556-566`): the first recorded image of `shape`.
pub(super) fn hist_first(
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_shape: &TopoShape,
) -> Option<TopoShape> {
    the_hist_map
        .get(&shape_map_key(the_shape))
        .and_then(|list| list.first().cloned())
}

/// `shape.Oriented(TopAbs_FORWARD)` (`TopoDS_Shape::Oriented`): a copy put
/// FORWARD, leaving the argument untouched.
pub(super) fn forward(the_shape: &TopoShape) -> TopoShape {
    let mut s = the_shape.clone();
    s.set_orientation(Orientation::Forward);
    s
}

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:185-198`): a straight
/// segment between two existing vertices. `None` is the reference's
/// `BRepLib_LineThroughIdenticPoints` (`cxx:191-195`); `cxx:196` then runs
/// `Init(GL, V1, V2, 0, l)`.
pub(super) fn make_edge_vv(v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let l = p1.distance(&p2);
    if l <= REAL_SMALL {
        return None;
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).ok()?;
    let line: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir)));
    let build = init_3d_vertices(&line, Some(v1), Some(v2), 0.0, l)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(line, build.first, build.last);
    add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `TopoDS_Iterator(shape)`'s components: `TopoDS_Iterator.cxx:32-46` walks
/// the shape's own child list, so the port reads it directly. The components
/// keep the orientation `TopoDS_Builder::Add` stored.
pub(super) fn shape_children(the_shape: &TopoShape) -> Vec<TopoShape> {
    let ts = the_shape.tshape.read().expect("poisoned TShape lock");
    ts.children.clone()
}

/// `BRep_Builder::Remove(theComp, theEdge)` on a compound held as an edge
/// list (`TopoDS_Builder::Remove`, `TopoDS_Builder.cxx:107-133`): every
/// component that `IsSame`s `the_edge` is dropped.
pub(super) fn comp_remove(the_comp: &mut Vec<TopoShape>, the_edge: &TopoShape) {
    the_comp.retain(|e| !is_same(e, the_edge));
}

/// `RemoveEdges` (`cxx:2520-2547`): the edges of `the_source_comp` that
/// `IsSame` none of `the_list_to_remove`.
pub(super) fn remove_edges(the_source_comp: &[TopoShape], the_list_to_remove: &[TopoShape]) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for an_edge in explorer_shapes_in(the_source_comp, ShapeType::Edge) {
        // `cxx:2531-2538`.
        let found = the_list_to_remove
            .iter()
            .any(|x| is_same(&an_edge, x));
        if !found {
            out.push(an_edge);
        }
    }
    out
}

/// `TopExp_Explorer(compound, TopAbs_EDGE)` over a compound given as an edge
/// list, in list order. The list only holds edges, so this is the identity
/// restricted to `Edge` components.
pub(super) fn explorer_shapes_in(the_comp: &[TopoShape], to_find: ShapeType) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for s in the_comp {
        if s.shape_type() == to_find {
            out.push(s.clone());
        }
    }
    out
}

/// `FilterSectionEdges` (`cxx:2549-2617`).
///
/// `MakeFacesSec` passes `aSecPlane` while it is still null (`cxx:800`,
/// `cxx:802`), so `theSecPlane.IsNull()` holds and the
/// `IntTools_BeanFaceIntersector` arm (`cxx:2576-2606`) never runs: every
/// pave block's section edge is added to the resulting compound
/// (`cxx:2608-2612`).
pub(super) fn filter_section_edges(the_b_curves: &[BopdsCurve], the_ds: &BopdsDS) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for a_b_curve in the_b_curves {
        for a_pb in a_b_curve.pave_blocks() {
            // `cxx:2571-2573`: `nSect = aPB->Edge(); aS = theDS->Shape(nSect)`.
            let Some(a_s) = the_ds.shape(a_pb.edge()) else {
                continue;
            };
            out.push(a_s.clone());
        }
    }
    out
}

/// `GetUEdges` (`cxx:1511-1575`): the two `u`-iso edges of the face column
/// `the_rank`, ordered so that the first one ends where the bound edge starts.
pub(super) fn get_u_edges(
    the_index: usize,
    the_rank: usize,
    the_u_edges: &ShapeTable,
    the_bound_edge: &TopoShape,
    the_face: &TopoShape,
) -> Option<(Edge, Edge)> {
    // `cxx:1519-1520`.
    let a_ue1 = the_u_edges.value(the_index, 1 + the_rank)?;
    let a_ue2 = the_u_edges.value(the_index + 1, 1 + the_rank)?;

    // `cxx:1522-1538`: the face's own occurrence of each.
    let mut a_face = the_face.clone();
    a_face.set_orientation(Orientation::Forward);
    let mut e1: Option<TopoShape> = None;
    let mut e2: Option<TopoShape> = None;
    for cur in explorer_shapes(&a_face, ShapeType::Edge) {
        if e1.is_none() && is_same(a_ue1, &cur) {
            e1 = Some(cur);
        } else if e2.is_none() && is_same(a_ue2, &cur) {
            e2 = Some(cur);
        }
    }
    // `cxx:1539-1542`.
    let (e1, e2) = (e1?, e2?);

    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face);
    // `cxx:1545-1551`: `BRep_Tool::CurveOnSurface(E1, aFace, f, l)`.
    let c1 = reg.edge_pcurve(&e1, key)?;
    let (c1_f, c1_l) = reg.pcurve_range(&e1, key)?;
    let pu1 = if the_rank == 0 {
        c1.d0(c1_l)
    } else {
        c1.d0(c1_f)
    };
    // `cxx:1552-1557`.
    let c2 = reg.edge_pcurve(the_bound_edge, key)?;
    // `cxx:1558-1559`: the 3D range, not the pcurve range.
    let (bf, _) = reg.edge_parameters(the_bound_edge);
    let pf = c2.d0(bf);

    // `cxx:1560-1561`.
    let a_v = if the_rank == 0 {
        last_vertex(&Edge(e1.clone()))
    } else {
        first_vertex(&Edge(e1.clone()))
    };
    let a_v = a_v?;
    let a_tolerance = BRepTool::vertex_tolerance(&a_v);

    // `cxx:1562-1568`: `BRepAdaptor_Surface::UResolution(aTolerance)`.
    let Some(surf) = reg.face_surface(&a_face) else {
        return Some((Edge(e1), Edge(e2)));
    };
    let u_tol = u_resolution(surf.as_ref(), a_tolerance);
    let pf_dist = ((pf.x() - pu1.x()).powi(2) + (pf.y() - pu1.y()).powi(2)).sqrt();
    if pf_dist > u_tol {
        Some((Edge(e2), Edge(e1)))
    } else {
        Some((Edge(e1), Edge(e2)))
    }
}

/// `BRepLib_MakeWire(edges)` (`BRepLib_MakeWire.cxx:123-453`) over an explicit
/// edge sequence. `None` is `!IsDone()`, i.e. the build failed and OCCT keeps
/// an empty `Shape()`.
pub(super) fn wire_from_edges(the_edges: &[Edge]) -> Option<Wire> {    let mut w = crate::brep_lib_make_wire::MakeWire::new();
    for e in the_edges {
        w.add(e).ok()?;
    }
    if w.is_done() {
        Some(w.wire())
    } else {
        None
    }
}

/// `BRepLib_MakeWire(W, E)` (`BRepLib_MakeWire.cxx:214-260`): the edges of
/// `the_wire` followed by `the_edge`.
pub(super) fn wire_with_edge(the_wire: &TopoShape, the_edge: &Edge) -> Option<Wire> {
    let mut edges = edges_of_wire(&Wire(the_wire.clone()));
    edges.push(the_edge.clone());
    wire_from_edges(&edges)
}

/// `BRepExtrema_ExtCC(E1, E2)` (`BRepExtrema_ExtCC.cxx`): the curve/curve
/// extrema of the two edges over their own parameter ranges, as
/// `(point on E1, point on E2, parameter on E1, parameter on E2)`.
pub(super) fn extrema_of_edges(the_e1: &Edge, the_e2: &Edge) -> Vec<(GpPnt, GpPnt, f64, f64)> {
    let reg = GeometryRegistry::global();
    let (Some(c1), Some(c2)) = (reg.edge_curve(&the_e1.0), reg.edge_curve(&the_e2.0)) else {
        return Vec::new();
    };
    let (a1, b1) = reg.edge_parameters(&the_e1.0);
    let (a2, b2) = reg.edge_parameters(&the_e2.0);
    occt_geom::extrema_cc::curve_curve_extrema_all_range(
        c1.as_ref(),
        c2.as_ref(),
        a1,
        b1,
        a2,
        b2,
    )
    .into_iter()
    .map(|p| (p.p1, p.p2, p.u1, p.u2))
    .collect()
}

// ---------------------------------------------------------------------------
// Section-walk statics (`FillGap` and its chain, `cxx:1577-2519`)
// ---------------------------------------------------------------------------

/// `shape.Oriented(theOrientation)` (`TopoDS_Shape::Oriented`): a copy of
/// `shape` carrying `the_orientation`.
pub(super) fn oriented_to(shape: &TopoShape, the_orientation: Orientation) -> TopoShape {
    let mut s = shape.clone();
    s.set_orientation(the_orientation);
    s
}

/// `TopExp::Vertices(E, Vfirst, Vlast, /*CumOri=*/true)` (`TopExp.cxx:214-252`).
///
/// The component orientations are composed with the edge's own orientation
/// (`TopoDS_Iterator(aShape, CumOri=true)`, `TopoDS_Iterator.cxx:32-46`) before
/// the FORWARD/REVERSED split, so a REVERSED edge reports the *stored* REVERSED
/// vertex first. `CheckAndOrientEdges` relies on exactly that (`cxx:1785`).
pub(super) fn edge_vertices_cum_ori(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let ori = edge.0.orientation();
    let stored = edge
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let mut first: Option<Vertex> = None;
    let mut last: Option<Vertex> = None;
    for child in stored {
        if child.shape_type() != ShapeType::Vertex {
            continue;
        }
        let cum = Orientation::compose(ori, child.orientation());
        if cum == Orientation::Forward && first.is_none() {
            first = Some(Vertex(child));
        } else if cum == Orientation::Reversed && last.is_none() {
            last = Some(Vertex(child));
        }
    }
    (first, last)
}

/// `BRepAdaptor_Surface(theFace).UResolution(tol)` / `VResolution(tol)`
/// (`BRepAdaptor_Surface.cxx:1410-1430` delegate to
/// `GeomAdaptor_Surface::UResolution`, `:1818-1896`): the reference takes the
/// larger of the two as the working 2D tolerance (`cxx:1729-1735`).
pub(super) fn surface_tolerance(the_face: &TopoShape, tol: f64) -> f64 {
    let reg = GeometryRegistry::global();
    let Some(surf) = reg.face_surface(the_face) else {
        return tol;
    };
    let u = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol);
    let v = occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol);
    if u > v {
        u
    } else {
        v
    }
}

/// `gp_Pnt2d::Distance` (`gp_Pnt2d.cxx`).
pub(super) fn pnt2d_distance(a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    (dx * dx + dy * dy).sqrt()
}

