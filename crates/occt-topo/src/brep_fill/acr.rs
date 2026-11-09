use super::prelude::*;


/// `BRepFill::ComputeACR(wire, ACR)` (`BRepFill.cxx:1027-1063`): the reduced
/// curvilinear abscissae of the wire's vertices.
///
/// `ACR` follows OCCT's `NCollection_Array1<double>(0, NbEdge)` bounds: entry
/// `0` holds the total length before the normalisation and is the divisor, so
/// the answer is returned as `ACR[0..=nb_edges]` with `ACR[0]` overwritten by
/// the normalised value `1.0` (or left at the total for a punctual wire).
///
/// `GCPnts_AbscissaPoint::Length(anEcur)` is the three-argument form: the
/// parameters are `BRepAdaptor_Curve`'s own first / last and there is no
/// tolerance argument.
pub fn compute_acr(wire: &Wire) -> Result<Vec<f64>, String> {
    let mut nb_edges = 0usize;
    let mut acr: Vec<f64> = vec![0.0];

    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(wire, None);
    while wexp.more() {
        nb_edges += 1;
        let current = wexp.current_edge();
        let Some(edge) = current else {
            wexp.next();
            continue;
        };
        // `ACR(nbEdges) = ACR(nbEdges - 1)` (`cxx:1039`).
        acr.push(acr[nb_edges - 1]);
        if !BRepTool::is_degenerated(&edge) {
            // `BRepAdaptor_Curve anEcur(Ecur); ACR(nbEdges) +=
            //  GCPnts_AbscissaPoint::Length(anEcur)` (`cxx:1042-1043`).
            let curve = BRepTool::edge_curve_adaptor_world(&edge)
                .ok_or_else(|| "Standard_NullObject: BRepAdaptor_Curve".to_string())?;
            acr[nb_edges] += occt_geom::gcpnts::curve_length_between(
                &*curve,
                curve.first_parameter(),
                curve.last_parameter(),
                None,
            );
        }
        wexp.next();
    }

    // `ACR(0) = ACR(nbEdges)` (`cxx:1048`).
    acr[0] = acr[nb_edges];

    // `cxx:1051-1062`.
    if acr[0] > CONFUSION {
        for i in 1..=nb_edges {
            acr[i] /= acr[0];
        }
    } else {
        // Punctual wire.
        acr[nb_edges] = 1.0;
    }
    Ok(acr)
}


/// `BRepFill::InsertACR(wire, ACRcuts, prec)` (`BRepFill.cxx:1067-1155`): cut
/// `wire` at every reduced abscissa of `acr_cuts` and rebuild it with
/// `BRepLib_MakeWire`.
///
/// `TrimEdge` is the `BRepFill.cxx:143-207` copy, i.e.
/// [`crate::brep_fill_compatible_wires::statics::trim_edge`] with
/// `reject_empty_piece = false`.
pub fn insert_acr(wire: &Wire, acr_cuts: &[f64], prec: f64) -> Result<Wire, String> {
    // `cxx:1072-1079`.
    let acr_wire = compute_acr(wire)?;

    let mut mw = crate::brep_lib_make_wire::MakeWire::new();
    let mut t1 = 0.0f64;
    let mut nb_edges = 0usize;

    // `cxx:1089-1147`: edge by edge.
    let mut wexp = crate::brep_tools_wire_explorer::WireExplorer::new();
    wexp.init(wire, None);
    while wexp.more() {
        nb_edges += 1;
        let t0 = t1;
        t1 = acr_wire[nb_edges];

        // `cxx:1097-1105`: the cuts strictly inside this edge's interval.
        let mut paradec: Vec<f64> = Vec::new();
        for &cut in acr_cuts {
            if t0 + prec < cut && cut < t1 - prec {
                paradec.push(cut);
            }
        }

        let e = match wexp.current_edge() {
            Some(e) => e,
            None => {
                wexp.next();
                continue;
            }
        };
        let v = wexp.current_vertex().cloned();

        if paradec.is_empty() || BRepTool::is_degenerated(&e) {
            // `cxx:1108-1112`: copy the edge.
            let _ = mw.add(&e);
        } else {
            // `cxx:1116-1140`.
            let so = match (&v, edge_vertices(&e).0) {
                (Some(v), Some(first)) => is_same(&v.0, &first.0),
                _ => false,
            };
            // "the wire is always FORWARD / it is necessary to modify the
            //  parameter of cut if the edge is REVERSED" (`cxx:1124-1138`).
            let sr: Vec<f64> = if e.orientation() == Orientation::Forward {
                paradec.clone()
            } else {
                paradec.iter().rev().map(|p| t0 + t1 - p).collect()
            };
            let se = crate::brep_fill_compatible_wires::statics::trim_edge(
                &e, &sr, t0, t1, so, false,
            );
            for piece in &se {
                let _ = mw.add(piece);
            }
        }
        wexp.next();
    }

    // `cxx:1148-1153`: the result keeps the input wire's orientation.
    let mut wres = mw.wire();
    wres.0.set_orientation(wire.0.orientation());
    Ok(wres)
}


/// `MakeWire(Edges, rangdeb, forward, newwire)` (`BRepFill.cxx:81-117`).
/// `rangdeb` is 1-based like OCCT's `NCollection_Array1<TopoDS_Shape>(1, n)`:
/// the `forward` arm appends `Edges(rangdeb..=n)` then `Edges(1..rangdeb)`, the
/// reverse arm prepends `Edges((1..=rangdeb).rev())` (already reversed edges)
/// then `Edges((rangdeb+1..=n).rev())`. Both arms finish with
/// `Orientation(FORWARD)` and `Closed(true)`.
pub(super) fn make_wire_from(
    edges: &[TopoShape],
    rangdeb: usize,
    forward: bool,
    newwire: &mut Wire,
) -> Result<(), String> {
    let nb_edges = edges.len();
    if rangdeb == 0 || rangdeb > nb_edges {
        return Err("Standard_OutOfRange: BRepFill::SearchOrigin rangdeb".to_string());
    }
    let b = TopoBuilder::new();
    let mut ordered: Vec<Edge> = Vec::with_capacity(nb_edges);
    let edge_at = |i: usize| -> Edge { Edge(edges[i - 1].clone()) };
    if forward {
        // `for (rang = rangdeb; rang <= nbEdges; rang++)`
        for rang in rangdeb..=nb_edges {
            ordered.push(edge_at(rang));
        }
        // `for (rang = 1; rang < rangdeb; rang++)`
        for rang in 1..rangdeb {
            ordered.push(edge_at(rang));
        }
    } else {
        // `for (rang = rangdeb; rang >= 1; rang--) E.Reversed()`
        for rang in (1..=rangdeb).rev() {
            let mut s = edge_at(rang).0;
            s.reverse();
            ordered.push(Edge(s));
        }
        // `for (rang = nbEdges; rang > rangdeb; rang--) E.Reversed()`
        for rang in ((rangdeb + 1)..=nb_edges).rev() {
            let mut s = edge_at(rang).0;
            s.reverse();
            ordered.push(Edge(s));
        }
    }
    *newwire = b.make_wire(&ordered);
    // `newwire.Orientation(TopAbs_FORWARD); newwire.Closed(true);`
    newwire.0.set_orientation(Orientation::Forward);
    newwire.0.set_closed(true);
    Ok(())
}


/// `CutEdge(CurrentEdge, Param, E1, E2, VRef)` (`BRepFill.cxx:119-142`).
///
/// Splits `current_edge` at `param`; `E1` is the piece that starts at `VRef`
/// (i.e. the piece the wire already owns and is parsed first).
pub(super) fn cut_edge_at(
    current_edge: &Edge,
    param: f64,
    v_ref: &Vertex,
) -> Result<(Edge, Edge), String> {
    let c = BRepTool::edge_curve_world(current_edge)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let (first, last) = BRepTool::edge_parameters(current_edge);
    // `B.MakeVertex(Vi, C->Value(Param), Precision::Confusion())` (`cxx:129`).
    let b = TopoBuilder::new();
    let vi = b.make_vertex(c.d0(param), CONFUSION);
    let (vf, vl) = edge_vertices(current_edge);
    let (Some(vf), Some(vl)) = (vf, vl) else {
        return Err("Standard_NullObject: TopExp::Vertices".to_string());
    };
    let piece = |v1: &Vertex, v2: &Vertex, p1: f64, p2: f64| -> Result<Edge, String> {
        let build = crate::brep_lib_make_edge::init_3d_vertices(&c, Some(v1), Some(v2), p1, p2)
            .ok_or_else(|| "BRepLib_EdgeDone".to_string())?;
        let mut e = b.make_edge(c.clone(), build.first, build.last);
        crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
        Ok(e)
    };
    if is_same(&v_ref.0, &vf.0) {
        // `E1 = MakeEdge(C, Vf, Vi, first, Param); E2 = MakeEdge(C, Vi, Vl, ...)`
        Ok((
            piece(&vf, &vi, first, param)?,
            piece(&vi, &vl, param, last)?,
        ))
    } else {
        Ok((
            piece(&vi, &vl, param, last)?,
            piece(&vf, &vi, first, param)?,
        ))
    }
}
