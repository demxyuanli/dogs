use super::prelude::*;

/// `BRep_Tool::CurveOnSurface(edge, face, f, l)` (`BRep_Tool.cxx:327-361`): the
/// edge's pcurve on the face together with the **`BRep_CurveOnSurface`
/// representation's own** parameter range (`GC->Range(First, Last)`, cxx:353),
/// not the edge's 3D range. `ShapeAnalysis_Edge::PCurve` (cxx:192-208) hands
/// exactly these `f`/`l` to its callers (the 2D branch of
/// `ShapeAnalysis_Wire::CheckOrder`, cxx:627-649, reads `c2d->Value(f)` /
/// `c2d->Value(l)`). Falls back to the 3D range when the edge carries no stored
/// COS range for this face (what the port's earlier code always used).
///
/// The distinction only matters when the two differ: after a split of a
/// periodic edge the child keeps the parent's pcurve but its COS range is the
/// sliced window (`ShapeBuild_Edge::CopyRanges`, `ShapeBuild_Edge.cxx:206-334`),
/// so evaluating at the 3D range would move every 2D endpoint by one period.
/// `ShapeAnalysis_Wire::CheckOrder` (cxx:593-651) and `ShapeAnalysis::TotCross2D`
/// both read the pcurve *through* this accessor, so the port has to as well.
fn pcurve_and_range(e: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    if let Some(v) = crate::boptools_2d::curve_on_surface_oriented(e, face, false) {
        return Some(v);
    }
    let pc = make_pcurve_full(e, face).ok()?;
    let fwd = Edge(e.0.oriented(Orientation::Forward));
    let (a, b) = BRepTool::edge_parameters(&fwd);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    Some((pc, a, b))
}

/// Port of `BRepMesh_ShapeVisitor::addWire` (2D pcurve mode).
///
/// Orders the wire's edges into one connected chain
/// (`ShapeAnalysis_Wire::CheckOrder` → `ShapeAnalysis_WireOrder::Perform`) and,
/// for every non-EXTERNAL edge, registers its pcurve on the face and its
/// position/orientation in the wire chain. Returns `false` when the wire cannot
/// be ordered (a missing pcurve → `ShapeExtend_FAIL`); a wire whose edges needed
/// reversing sets the face `UNORIENTED_WIRE` status (`ShapeExtend_DONE3`).
fn add_wire(
    model: &mut MeshModel,
    face_index: usize,
    wire: &Wire,
    edge_index: &mut HashMap<usize, usize>,
) -> bool {
    let face = match model.face(face_index) {
        Ok(f) => f.face().clone(),
        Err(_) => return false,
    };
    // `ShapeExtend_WireData(theWire, chained=true, manifold=false)` keeps only
    // FORWARD/REVERSED edges in the main list; INTERNAL/EXTERNAL are non-manifold.
    // `ShapeExtend_WireData::Init`: FORWARD/REVERSED only; a REVERSED wire
    // prepends each edge so the list is the reverse of iterator order
    // (`ShapeExtend_WireData.cxx:115-122`).
    //
    // OCCT's mesh model builder applies NO orientation fixup here:
    // `BRepMesh_ShapeVisitor::addWire` (`BRepMesh_ShapeVisitor.cxx:93-134`) only
    // runs `ShapeExtend_WireData(wire, true, false)` + `ShapeAnalysis_Wire::
    // CheckOrder` and stores the resulting chain as-is. `ShapeFix_Face::
    // FixOrientation` is a `ShapeFix` pass and is not part of the mesher, so the
    // wire's UV winding is left untouched (the clockwise/ccw split is the
    // classifier's job, not this loop's).
    //
    // Reversing a single clockwise wire here flipped face 433/435 of
    // `data/occ/acs10.stp` from `ShapeAnalysis_WireOrder`'s reference chain
    // (`e0,e15,e14,...,e1`) to its mirror (`e15,e0,...,e14`) and starved the
    // Delaunay cavity: 49/50 triangles instead of OCCT's 73/74.
    let mut work = wire.clone();
    let face_fwd = Face(face.0.oriented(Orientation::Forward));
    let mut stored: Vec<Edge> = edges_of_wire(&work)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if work.0.orientation().is_reversed() {
        stored.reverse();
    }
    if stored.is_empty() {
        return false;
    }

    // `ShapeAnalysis_Wire::CheckOrder(..., isClosed=true, mode3d=false)`:
    // pcurve endpoints on `face.Oriented(FORWARD)` with `PCurve(..., orient=true)`
    // (reversed edges swap cf/cl).
    let mut order = WireOrder::new();
    for e in &stored {
        // Port of `ShapeAnalysis_Wire::CheckOrder` (`ShapeAnalysis_Wire.cxx:627-649`):
        //   ShapeAnalysis_Edge EA; EA.PCurve(E, face.Oriented(FORWARD), c2d, f, l);
        //   gp_Pnt2d (c2d->Value(f)), (c2d->Value(l));
        // `PCurve(..., orient=true)` (`ShapeAnalysis_Edge.cxx:192-208`) reads
        // `BRep_Tool::CurveOnSurface(E, F, f, l)` - the **COS representation's**
        // range - and swaps f/l on a REVERSED edge. Seam edges have two pcurves
        // (one per side) and each traversal must get its own side, so the pcurve
        // is selected by the edge's *actual* orientation.
        let (pc, a, b) = match pcurve_and_range(e, &face_fwd) {
            Some(v) => v,
            None => return false,
        };
        if !a.is_finite() || !b.is_finite() || (b - a).abs() < 1e-15 {
            return false;
        }
        // `ShapeAnalysis_Wire::CheckOrder` (`cxx:648-649`): `c2d->Value(f/l)`
        // after `PCurve(..., orient=true)` toggles f/l on a REVERSED edge.
        let (p_a, p_b) = (pc.d0(a), pc.d0(b));
        let (begin, end) = if e.0.orientation().is_reversed() {
            (p_b, p_a)
        } else {
            (p_a, p_b)
        };
        order.add_edge(begin, end);
    }
    order.perform();
    // `BRepMesh_ShapeVisitor::addWire` (`BRepMesh_ShapeVisitor.cxx:99-134`)
    // applies no chain-level reversal after `CheckOrder`: it only reads
    // `aOrderTool.Ordered(i)` and lets `ShapeExtend_WireData::Edge(signed)`
    // reverse each edge individually (`ShapeExtend_WireData.cxx:583-588`). The
    // wire's UV winding stays whatever the oriented edges give (hole wires stay
    // clockwise); the inside/outside split is done later by the classifier, not
    // by flipping the whole chain on the sign of a signed area.

    if order.status() == WireOrderStatus::Reversed {
        if let Ok(f) = model.face_mut(face_index) {
            f.set_status(MeshStatus::UNORIENTED_WIRE);
        }
    }
    if order.nb_edges() != stored.len() {
        return false;
    }

    // `ShapeExtend_WireData::Edge(signed)` (`cxx:583-588`): `Ordered < 0`
    // is `Edge.Reverse()`, then `AddPCurve` / `AddEdge` use that orientation
    // (`BRepMesh_ShapeVisitor.cxx:123-130`). On a seam `Reverse` selects
    // PCurve2 (`BRep_Tool.cxx:354-357`).
    let wire_index = model.add_wire(wire.clone());
    // TEMPORARY diagnostic (DELREGSKEL): the `ShapeExtend_WireData` list fed to
    // `ShapeAnalysis_WireOrder::Perform` and the `Ordered(i)` chain it returns, in
    // the `occt_probe --orderdbg` `ODBG` line format.
    if std::env::var("DELREGSKEL").is_ok() {
        {
            let raw = work.0.tshape.read().expect("lock").children.clone();
            let rawlist = raw
                .iter()
                .map(|c| format!("{:?}:{:?}", c.shape_type(), c.orientation()))
                .collect::<Vec<_>>()
                .join(" ");
            eprintln!(
                "WIRE face={} wire={} work_ori={:?} raw_children=[{}]",
                face_index, wire_index, work.0.orientation(), rawlist
            );
        }
        let list = stored
            .iter()
            .map(|e| format!("{:p}:{:?}", Arc::as_ptr(&e.0.tshape), e.0.orientation()))
            .collect::<Vec<_>>()
            .join(" ");
        let chain = (1..=stored.len())
            .map(|i| order.ordered(i).to_string())
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "WIRE face={} wire={} ori={:?} stored=[{}] ordered=[{}]",
            face_index,
            wire_index,
            work.0.orientation(),
            list,
            chain
        );
        for (i, e) in stored.iter().enumerate() {
            let n_pc = GeometryRegistry::global()
                .edge_pcurves(&e.0, GeometryRegistry::shape_key(&face_fwd.0))
                .len();
            match pcurve_and_range(e, &face_fwd) {
                Some((pc, a, b)) => {
                    let (pa, pb) = (pc.d0(a), pc.d0(b));
                    let (bu, eu) = if e.0.orientation().is_reversed() {
                        (pb, pa)
                    } else {
                        (pa, pb)
                    };
                    eprintln!(
                        "WIRE   i={} ptr={:p} ori={:?} npc={} pc2d=({:.9},{:.9})->({:.9},{:.9}) trav=({:.9},{:.9})->({:.9},{:.9})",
                        i + 1,
                        Arc::as_ptr(&e.0.tshape),
                        e.0.orientation(),
                        n_pc,
                        pa.x(),
                        pa.y(),
                        pb.x(),
                        pb.y(),
                        bu.x(),
                        bu.y(),
                        eu.x(),
                        eu.y()
                    );
                }
                None => {
                    eprintln!(
                        "WIRE   i={} ptr={:p} ori={:?} npc={} no-pcurve",
                        i + 1,
                        Arc::as_ptr(&e.0.tshape),
                        e.0.orientation(),
                        n_pc
                    );
                }
            }
        }
    }
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        let e = &stored[signed.unsigned_abs() as usize - 1];
        let mut orientation = e.0.orientation();
        if orientation == Orientation::External {
            continue;
        }        // `ShapeExtend_WireData::Edge(signed)` is `Edge.Reverse()`. On a plane
        // that is the same pcurve walked backwards. On a seam `Reverse` selects
        // PCurve2 (`BRep_Tool.cxx:354-357`), leaving the CheckOrder iso; walk
        // that measured chord backwards instead.
        let mut reverse_walk = false;
        if signed < 0 {
            let face_key = GeometryRegistry::shape_key(&face_fwd.0);
            let n_pc = GeometryRegistry::global().edge_pcurves(&e.0, face_key).len();
            if n_pc >= 2 {
                reverse_walk = true;
            } else {
                orientation = orientation.reversed();
            }
        }
        let key = Arc::as_ptr(&e.0.tshape) as usize;
        let eidx = match edge_index.get(&key) {
            Some(&i) => i,
            None => {
                let i = model.add_edge(e.clone());
                edge_index.insert(key, i);
                i
            }
        };
        if let Ok(edge) = model.edge_mut(eidx) {
            edge.add_pcurve(face_index, orientation);
        }
        if let Ok(w) = model.wire_mut(wire_index) {
            w.add_edge_ordered(eidx, orientation, reverse_walk);
        }
    }
    if let Ok(f) = model.face_mut(face_index) {
        f.add_wire(wire_index);
    }
    true
}

/// `ShapeAnalysis::OuterWire` over an already-oriented wire list.
pub(crate) fn outer_of_wires(wires: &[Wire], face: &Face) -> Option<Wire> {
    if wires.is_empty() {
        return None;
    }
    for (i, w) in wires.iter().enumerate() {
        if i == wires.len() - 1 {
            return Some(w.clone());
        }
        if wire_area_2d(w, face) >= 0.0 {
            return Some(w.clone());
        }
    }
    None
}

/// `ShapeExtend_WireData(wire, chained=true)` edge list: FORWARD/REVERSED
/// only; a REVERSED wire prepends so iterator order is reversed
/// (`ShapeExtend_WireData.cxx:115-122`).
fn wire_edges_sewd(wire: &Wire) -> Vec<Edge> {
    let mut stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if wire.0.orientation().is_reversed() {
        stored.reverse();
    }
    stored
}

fn edge_oriented_ends(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (a, b) = edge_vertices(e);
    if e.0.orientation().is_reversed() {
        (b, a)
    } else {
        (a, b)
    }
}

/// `BRepTools_WireExplorer` — the wire's edges in connection order on `face`.
///
/// Port of `BRepTools_WireExplorer::Init` (`BRepTools_WireExplorer.cxx:124-382`)
/// and `Next` (`cxx:394-705`). Init keys every edge by its cumulated first
/// vertex (`TopExp::Vertices(E, V1, V2, true)`, `cxx:244`) and starts from the
/// wire's first vertex (`cxx:303-373`). Next advances through
/// `TopExp::LastVertex(myEdge, true)` and takes the next edge from that vertex's
/// list. Storage order only decides the *first* edge: when several edges start
/// at the current vertex (`cxx:410-705`) OCCT runs the list through the
/// closed-edge-only pass, then the open-edge pass, and inside a pass picks the
/// edge whose pcurve start lies within `(myTolU, myTolV)` of the current 2D
/// point and whose direction turns least from the current edge's backward
/// tangent (`cxx:536-695`). Walking by vertex identity alone picks the wrong
/// occurrence of a seam edge on a closed surface (the cylinder lateral wire
/// would visit `seam(u = 2*pi)` instead of the bottom circle), which is what
/// makes the sampled UV ring self-intersect.
pub(crate) fn wire_edges_explorer(wire: &Wire, face: &Face) -> Vec<Edge> {
    // `BRepTools_WireExplorer::Init` leaves `myEdge` null in every path that
    // bails out here (`BRepTools_WireExplorer.cxx:160-164`, `:368-372`), so the
    // caller sees no edge at all.
    macro_rules! bail {
        ($why:expr) => {{
            let _ = $why;
            return Vec::new();
        }};
    }
    let edges: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o != Orientation::Internal && o != Orientation::External
        })
        .collect();
    if edges.is_empty() {
        bail!("no edges");
    }

    // `Init` (`cxx:139-217`): UV tolerances from the largest vertex tolerance of
    // the wire (falling back to the edge tolerances below Confusion).
    let mut df_vert_toler = 0.0f64;
    for e in &edges {
        let (a, b) = edge_oriented_ends(e);
        for v in [a, b].into_iter().flatten() {
            df_vert_toler = df_vert_toler.max(BRepTool::vertex_tolerance(&v));
        }
    }
    if df_vert_toler < CONFUSION {
        for e in &edges {
            df_vert_toler = df_vert_toler.max(BRepTool::edge_tolerance(e));
        }
        if df_vert_toler < CONFUSION {
            // empty wire (`cxx:160-164`)
            bail!("vert toler under confusion");
        }
    }
    let Some(surf) = BRepTool::face_surface_world(&face) else {
        bail!("no face surface");
    };
    let mut tol_u =
        2.0 * occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), df_vert_toler);
    let mut tol_v =
        2.0 * occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), df_vert_toler);
    {
        // The two-argument `BRepTools_WireExplorer(W, F)` entry point
        // (`cxx:100-116`) fills the UV bounds only for Cone / BSpline / Bezier
        // faces, from `BRepTools::UVBounds`; every other surface type keeps them
        // at 0.
        let (umin, umax, vmin, vmax) = match mesh_surface_type(surf.as_ref()) {
            MeshSurfaceType::Cone
            | MeshSurfaceType::BSplineSurface
            | MeshSurfaceType::BezierSurface => crate::brep_tools::uv_bounds(face),
            _ => (0.0, 0.0, 0.0, 0.0),
        };
        match mesh_surface_type(surf.as_ref()) {
            MeshSurfaceType::Cone => {
                // `cxx:170-200`.
                let maxtol = 0.0005 * (umax - umin);
                let a1 = surf.d1(umin, vmin).1.magnitude();
                let tol1 = if a1 <= CONFUSION {
                    maxtol
                } else {
                    maxtol.min(df_vert_toler / a1)
                };
                let a2 = surf.d1(umin, vmax).1.magnitude();
                let tol2 = if a2 <= CONFUSION {
                    maxtol
                } else {
                    maxtol.min(df_vert_toler / a2)
                };
                tol_u = 2.0 * tol1.max(tol2);
            }
            MeshSurfaceType::BSplineSurface | MeshSurfaceType::BezierSurface => {
                // `cxx:202-217`.
                let mut max_tol = tol_u.max(tol_v);
                let (_, du, dv) = surf.d1((umax - umin) / 2.0, (vmax - vmin) / 2.0);
                let modulus = (du.dot(&du) + dv.dot(&dv)).sqrt();
                if modulus > REAL_SMALL {
                    if modulus * max_tol / df_vert_toler < 1.5 {
                        max_tol = 1.5 * df_vert_toler / modulus;
                    }
                    tol_u = max_tol;
                    tol_v = max_tol;
                }
            }
            _ => {}
        }
    }

    // `myMap` / `vmap` (`cxx:226-302`): `myMap(V)` is the list of edges whose
    // cumulated first vertex is `V`, appended in wire order; `vmap` keeps the
    // vertices met an odd number of times (the open ends of the wire).
    let mut map: Vec<(usize, Vec<Edge>)> = Vec::new();
    let mut vmap: Vec<(usize, Orientation)> = Vec::new();
    for e in &edges {
        let (v1, v2) = edge_oriented_ends(e);
        if let Some(v1) = v1 {
            let key = vertex_key(&v1);
            match map.iter_mut().find(|(k, _)| *k == key) {
                Some((_, l)) => l.push(e.clone()),
                None => map.push((key, vec![e.clone()])),
            }
            vmap_toggle(&mut vmap, key, Orientation::Forward);
        }
        if let Some(v2) = v2 {
            vmap_toggle(&mut vmap, vertex_key(&v2), Orientation::Reversed);
        }
    }

    // First vertex (`cxx:303-373`): the wire's open end when it has one, else the
    // first vertex of the first edge.
    let first_key = if !vmap.is_empty() {
        vmap.iter()
            .find(|(_, o)| *o == Orientation::Forward)
            .map(|(k, _)| *k)
    } else {
        if let Some(e) = edges.iter().find(|e| {
            let (v1, v2) = edge_oriented_ends(e);
            v1.is_none() || v2.is_none()
        }) {
            // Infinite edge (`cxx:326-348`): the walk is that single edge.
            {
                let (f, l) = BRepTool::edge_parameters(e);
                let rev = e.0.orientation().is_reversed();
                let infinite = if rev {
                    l == f64::NEG_INFINITY || l == f64::INFINITY
                } else {
                    f == f64::NEG_INFINITY || f == f64::INFINITY
                };
                if infinite {
                    return vec![e.clone()];
                }
            }
        }
        edge_oriented_ends(&edges[0]).0.map(|v| vertex_key(&v))
    };
    let Some(first_key) = first_key else {
        bail!("no first vertex");
    };
    let Some(entry) = map.iter_mut().find(|(k, _)| *k == first_key) else {
        bail!("first vertex not in map");
    };
    let mut out = vec![entry.1.remove(0)];

    while let Some(lv) = edge_oriented_ends(out.last().unwrap()).1 {
        let Some(entry) = map.iter_mut().find(|(k, _)| *k == vertex_key(&lv)) else {
            break;
        };
        let list = &mut entry.1;
        if list.is_empty() {
            break;
        }
        if list.len() == 1 {
            // `cxx:414-482`: a single candidate, with the pcurve-end check for a
            // closed candidate edge.
            let candidate = list[0].clone();
            let (v1, v2) = edge_oriented_ends(&candidate);
            let Some(v1) = v1 else { break };
            if !same_vertex(&v1, &lv) {
                break;
            }
            // `cxx:429-475`: a closed candidate (`aV1.IsSame(aV2)`) is only
            // accepted when both pcurves exist and the current edge's 2D
            // endpoint is nearer the candidate's 2D start than its 2D end.
            if let Some(v2) = v2 {
                if same_vertex(&v1, &v2) {
                    let (Some((prev_pc, p11, p12)), Some((next_pc, p21, p22))) = (
                        pcurve_and_range(out.last().unwrap(), face),
                        pcurve_and_range(&candidate, face),
                    ) else {
                        break;
                    };
                    let prev_par = if out.last().unwrap().0.orientation().is_reversed() {
                        p11
                    } else {
                        p12
                    };
                    let (next_fpar, next_lpar) = if candidate.0.orientation().is_reversed() {
                        (p22, p21)
                    } else {
                        (p21, p22)
                    };
                    let p_prev = prev_pc.d0(prev_par);
                    let p_next_f = next_pc.d0(next_fpar);
                    let p_next_l = next_pc.d0(next_lpar);
                    if p_prev.square_distance(&p_next_f) > p_prev.square_distance(&p_next_l) {
                        break;
                    }
                }
            }
            list.clear();
            out.push(candidate);
        } else {
            // `cxx:521-705`: several candidates start at the current vertex.
            let current = out.last().unwrap().clone();
            let Some(idx) = wire_next_by_uv(&current, face, list, tol_u, tol_v) else {
                break;
            };
            out.push(list.remove(idx));
        }
    }
    out
}

/// `TopTools_ShapeMapHasher::IsSame` key for a vertex (`TopoDS_Shape`
/// identity = `TShape` + location).
fn vertex_key(v: &Vertex) -> usize {
    GeometryRegistry::shape_key(&v.0)
}
fn same_vertex(a: &Vertex, b: &Vertex) -> bool {
    vertex_key(a) == vertex_key(b)
}

/// `vmap.Add` / `RemoveKey` (`BRepTools_WireExplorer.cxx:252-268`): the map
/// compares shapes by `IsSame`, so a key added twice is removed again and only
/// the last orientation written survives.
fn vmap_toggle(vmap: &mut Vec<(usize, Orientation)>, key: usize, ori: Orientation) {
    if let Some(pos) = vmap.iter().position(|(k, _)| *k == key) {
        vmap.remove(pos);
    } else {
        vmap.push((key, ori));
    }
}

/// `BRepTools_WireExplorer::Next` multi-candidate selection
/// (`BRepTools_WireExplorer.cxx:521-705`): returns the index into `list` of the
/// next edge, or `None` when the walk stops (`myEdge = TopoDS_Edge()`).
fn wire_next_by_uv(
    my_edge: &Edge,
    face: &Face,
    list: &[Edge],
    tol_u: f64,
    tol_v: f64,
) -> Option<usize> {
    let (a_pcurve, df_f_par, df_l_par) = pcurve_and_range(my_edge, face)?;
    // The 2D point of the current vertex on the current edge (`cxx:548-547`):
    // MAX pcurve parameter for a FORWARD edge, MIN for a REVERSED one.
    let my_rev = my_edge.0.orientation().is_reversed();
    let p_ref = if my_rev {
        a_pcurve.d0(df_f_par)
    } else {
        a_pcurve.d0(df_l_par)
    };
    // The direction from the current vertex back into the current edge
    // (`cxx:551-556`), used as the angle reference.
    let df_m_par =
        get_next_param_on_pc(a_pcurve.as_ref(), &p_ref, df_f_par, df_l_par, tol_u, tol_v, my_rev);
    let p_refm = a_pcurve.d0(df_m_par);
    let an_e_ref_dir = GpVec2d::new(p_refm.x() - p_ref.x(), p_refm.y() - p_ref.y());
    if an_e_ref_dir.square_magnitude() < REAL_SMALL {
        return None;
    }

    let mut k_min = 0usize;
    let mut is_degenerated = true;
    let mut dmin = f64::MAX;
    let mut df_min_angle = 3.0 * PI;
    let mut df_cur_angle = 3.0 * PI;
    for _pass in 0..2 {
        let mut k = 1usize;
        for e in list {
            if Arc::ptr_eq(&e.0.tshape, &my_edge.0.tshape) {
                k += 1;
                continue;
            }
            let (a_v1, a_v2) = edge_oriented_ends(e);
            let (Some(a_v1), Some(a_v2)) = (a_v1, a_v2) else {
                k += 1;
                continue;
            };
            let Some((a_pc, a_f, a_l)) = pcurve_and_range(e, face) else {
                k += 1;
                continue;
            };
            // `cxx:602`: pass 1 keeps closed edges (`V1.IsSame(V2)`), pass 2 the
            // rest.
            if same_vertex(&a_v1, &a_v2) == is_degenerated {
                let a_peb = if e.0.orientation().is_reversed() {
                    a_pc.d0(a_l)
                } else {
                    a_pc.d0(a_f)
                };
                if (a_l - a_f).abs() > PCONFUSION {
                    let reverse = !e.0.orientation().is_reversed();
                    let a_epm =
                        get_next_param_on_pc(a_pc.as_ref(), &a_peb, a_f, a_l, tol_u, tol_v, reverse);
                    let mut a_pee = a_pc.d0(a_epm);
                    if a_peb.square_distance(&a_pee) <= REAL_SMALL {
                        // Very short curve: take its tangent (`cxx:626-643`).
                        let (_, d) = a_pc.d1(a_epm);
                        a_pee = if e.0.orientation().is_reversed() {
                            GpPnt2d::new(a_peb.x() - d.x(), a_peb.y() - d.y())
                        } else {
                            GpPnt2d::new(a_peb.x() + d.x(), a_peb.y() + d.y())
                        };
                        if a_peb.square_distance(&a_pee) <= REAL_SMALL {
                            k += 1;
                            continue;
                        }
                    }
                    let an_e_dir = GpVec2d::new(a_pee.x() - a_peb.x(), a_pee.y() - a_peb.y());
                    df_cur_angle = an_e_dir.angle(&an_e_ref_dir).abs();
                }
                if df_cur_angle <= df_min_angle {
                    let mut d = p_ref.square_distance(&a_peb);
                    if d <= PCONFUSION {
                        d = 0.0;
                    }
                    if (a_peb.x() - p_ref.x()).abs() < tol_u
                        && (a_peb.y() - p_ref.y()).abs() < tol_v
                        && d <= dmin
                    {
                        df_min_angle = df_cur_angle;
                        k_min = k;
                        dmin = d;
                    }
                }
            }
            k += 1;
        }
        if k_min == 0 {
            is_degenerated = false;
            dmin = f64::MAX;
        } else {
            break;
        }
    }
    if k_min == 0 { None } else { Some(k_min - 1) }
}

/// `GetNextParamOnPC` (`BRepTools_WireExplorer.cxx:798-873`): the parameter of
/// the first sample of `pc` more than `(tol_u, tol_v)` away from `p_ref`,
/// stepping from the `reverse` end of `[f, l]` in `|l - f| / 1000` increments.
fn get_next_param_on_pc(
    pc: &dyn Curve2d,
    p_ref: &GpPnt2d,
    f: f64,
    l: f64,
    tol_u: f64,
    tol_v: f64,
    reverse: bool,
) -> f64 {
    let d_p = (l - f).abs() / 1000.0;
    let on_ref = |p: &GpPnt2d| {
        (p_ref.x() - p.x()).abs() < tol_u && (p_ref.y() - p.y()).abs() < tol_v
    };
    if reverse {
        let mut result = f;
        let mut start_par = f;
        let mut found = false;
        while !found && start_par < l {
            start_par += d_p;
            if on_ref(&pc.d0(start_par)) {
                continue;
            }
            result = start_par;
            found = true;
        }
        if !found {
            result = l;
        }
        if result > l {
            result = l;
        }
        result
    } else {
        let mut result = l;
        let mut start_par = l;
        let mut found = false;
        while !found && start_par > f {
            start_par -= d_p;
            if on_ref(&pc.d0(start_par)) {
                continue;
            }
            result = start_par;
            found = true;
        }
        if !found {
            result = f;
        }
        if result < f {
            result = f;
        }
        result
    }
}

/// `ShapeAnalysis::TotCross2D` — signed 2D area of a wire's pcurves (trapezoid
/// rule over sampled pcurve points, sequence reversed per REVERSED edge).
fn wire_area_2d(wire: &Wire, face: &Face) -> f64 {
    let mut totcross = 0.0;
    let mut uv0: Option<GpPnt2d> = None;
    let mut fuv = GpPnt2d::new(0.0, 0.0);
    let mut nbc = 0usize;
    for e in wire_edges_sewd(wire) {
        // `BRep_Tool::CurveOnSurface(edge, face)` picks PCurve2 on a reversed
        // seam (`BRep_Tool.cxx:354-357`) and returns the COS range; then
        // TotCross2D reverses the samples for a REVERSED edge.
        let Some((pc, a, b)) = pcurve_and_range(&e, face) else {
            continue;
        };
        let mut pts = sample_pcurve(pc.as_ref(), a, b);
        if e.0.orientation().is_reversed() {
            pts.reverse();
        }
        nbc += 1;
        if nbc == 1 {
            fuv = pts[0];
            uv0 = Some(pts[0]);
        }
        for p in &pts {
            totcross += (fuv.x() - p.x()) * (fuv.y() + p.y()) / 2.0;
            fuv = *p;
        }
    }
    if let Some(u0) = uv0 {
        totcross += (fuv.x() - u0.x()) * (fuv.y() + u0.y()) / 2.0;
    }
    totcross
}

/// Uniform sample of a pcurve over `[a, b]` (endpoints + interior) with the
/// point count OCCT uses for both `ShapeAnalysis::TotCross2D`
/// (`ShapeAnalysis_Curve::GetSamplePoints`, `ShapeAnalysis_Curve.cxx:1317-1337`)
/// and the `BRepTopAdaptor_FClass2d` boundary polygon
/// (`BRepTopAdaptor_FClass2d.cxx:179-185`).
pub(crate) fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    let n = sample_count_2d(pc, a, b);
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = a + (b - a) * i as f64 / (n - 1) as f64;
        pts.push(pc.d0(t));
    }
    pts
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:73-91`) plus the `nbs *= 4` applied by the
/// two callers that share this count:
///
/// * `ShapeAnalysis_Curve.cxx:1325-1331` (`GetSamplePoints` for `TotCross2D`)
/// * `BRepTopAdaptor_FClass2d.cxx:182-184`
///
/// Both apply the `*4` at their own call site, so it stays here rather than in
/// [`crate::curve_sampling_2d::geom2d_nb_samples`].
pub(crate) fn sample_count_2d(pc: &dyn Curve2d, a: f64, b: f64) -> usize {
    let mut nbs = crate::curve_sampling_2d::geom2d_nb_samples(pc, a, b);
    if nbs > 2 {
        nbs *= 4;
    }
    // `nbPoints` never returns less than 2 (the `Geom2d_Line` arm).
    nbs.max(2)
}

/// Tool for building a discrete model from a topological shape.
/// Port of `BRepMesh_ModelBuilder` (+ the visiting logic of
/// `BRepMesh_ShapeVisitor`).
pub struct ModelBuilder;

impl ModelBuilder {
    /// Build the discrete model of `shape` under `params`.
    ///
    /// Mirrors `BRepMesh_ModelBuilder::performInternal` plus the
    /// `BRepMesh_ShapeVisitor` walk: every face gets its outer wire first
    /// (`ShapeAnalysis::OuterWire`), then its inner wires, each wire reordered
    /// into a connected chain (`ShapeAnalysis_Wire::CheckOrder`), with EXTERNAL
    /// edges skipped and shared edges deduplicated by TShape identity. Free edges
    /// (not bounding any face) are added as well.
    ///
    /// Returns `Err` when the shape is empty (void bounding box → `Message_Fail1`).
    pub fn build_model(shape: &TopoShape, params: &MeshParameters) -> Result<MeshModel, String> {
        let mut model = MeshModel::new(shape.clone());

        // Maximum size of the shape's bounding box (used by relative deflection).
        let bbox = shape_bbox(shape);
        if bbox.is_void() {
            return Err("BRepMesh_ModelBuilder::build_model: empty shape".to_string());
        }
        let max_size = if params.relative {
            box_max_dimension(&bbox).unwrap_or(0.0)
        } else {
            params.deflection.max(params.deflection_interior.max(0.0))
        };
        model.set_max_size(max_size);

        // Visit(Face): `IMeshTools_ShapeExplorer.cxx:101` stores only
        // `Oriented(FORWARD)` faces "to prevent inverse issue". A REVERSED
        // face would Compose hole wires to FORWARD and fill the hole.
        // Outer wire first; a failure on the outer wire fails the face, a
        // failure on an inner wire only marks it unoriented.
        let mut edge_index: HashMap<usize, usize> = HashMap::new();
        for f in faces_of(shape) {
            if std::env::var("SWEEP_DBG").is_ok() {
                let s = BRepTool::face_surface(&f);
                let sd = match &s {
                    Some(s) => {
                        let o = s.d0(0.0, 0.0);
                        let (_, du, dv) = s.d1(0.0, 0.0);
                        format!(
                            "o=({:.3},{:.3},{:.3}) du=({:.3},{:.3},{:.3}) dv=({:.3},{:.3},{:.3})",
                            o.x(), o.y(), o.z(), du.x(), du.y(), du.z(), dv.x(), dv.y(), dv.z()
                        )
                    }
                    None => "none".to_string(),
                };
                eprintln!(
                    "[model-dbg] face tsh={:?} loc={:?} surf={}",
                    Arc::as_ptr(&f.0.tshape) as usize,
                    f.0.location().transformation(),
                    sd
                );
            }
            let f = Face(f.0.oriented(Orientation::Forward));
            let face_index = model.add_face(f.clone());

            // `BRepMesh_ShapeVisitor::Visit` (`BRepMesh_ShapeVisitor.cxx:58-87`)
            // applies NO orientation fixup: it is `ShapeAnalysis::OuterWire` +
            // `TopExp_Explorer(face, TopAbs_WIRE)`, and `addWire` only runs
            // `ShapeAnalysis_Wire::CheckOrder`. `ShapeFix_Face::FixOrientation`
            // belongs to the Shape Healing pass the reader already ran
            // (`IGESControl_IGESBoundary` -> `ShapeFix_Face::Perform`,
            // port: `iges_read/iges_process_shape.rs`); repeating it inside the
            // mesher reversed the torus' `V in [-pi, pi]` wire and made the
            // FaceChecker report a self-intersection (`FaceChecker.cxx:150-158`)
            // so the face came out with zero triangles.
            let wires = wires_of_face(&f);
            let outer = outer_of_wires(&wires, &f);
            if let Some(outer_wire) = &outer {
                if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::FAILURE);
                    }
                    continue;
                }
            }

            for w in &wires {
                if let Some(outer_wire) = &outer {
                    if Arc::ptr_eq(&w.0.tshape, &outer_wire.0.tshape) {
                        continue;
                    }
                }
                if !add_wire(&mut model, face_index, w, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::UNORIENTED_WIRE);
                    }
                }
            }
        }

        // Free edges that do not bound any face.
        for e in edges_of(shape) {
            let key = Arc::as_ptr(&e.0.tshape) as usize;
            if !edge_index.contains_key(&key) {
                model.add_edge(e.clone());
            }
        }

        Ok(model)
    }
}
