//! `StepToTopoDS_TranslateEdge`: `EDGE_CURVE` / `ORIENTED_EDGE` resolution and the edge curve sense.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let start_ref = parse_ref(&rec.args[1]).ok_or("EDGE_CURVE: bad start ref")?;
        let end_ref = parse_ref(&rec.args[2]).ok_or("EDGE_CURVE: bad end ref")?;
        let curve_ref = parse_ref(&rec.args[3]).ok_or("EDGE_CURVE: bad curve ref")?;
        // `StepToTopoDS_TranslateEdge.cxx:290-322`: `same_sense` picks which
        // VERTEX is the edge's geometric first / last. `.T.` keeps
        // `edge_start -> edge_end`; `.F.` swaps them, so the projected range
        // and the vertex children both follow the curve's own direction.
        let same_sense = parse_logical(rec.args.get(4).map(String::as_str), true);
        let (first_ref, last_ref) = if same_sense {
            (start_ref, end_ref)
        } else {
            (end_ref, start_ref)
        };
        let v1 = self.resolve_shape(first_ref)?;
        let v2 = self.resolve_shape(last_ref)?;
        if !v1.is_vertex() || !v2.is_vertex() {
            return Err("EDGE_CURVE: endpoints are not vertices".into());
        }
        let curve = self.resolve_curve(curve_ref)?;
        // `MakeFromCurve3D` (`TranslateEdge.cxx:433-441`): the projection points
        // are the two endpoint vertices; when the two translated vertices are the
        // same `TopoDS_Shape` (`V1.IsSame(V2)`), `GetCartesianPoints` re-reads
        // them from the EDGE_CURVE's own `edge_start` / `edge_end` vertex geometry
        // (STEP order, i.e. **not** the `same_sense`-swapped order).
        let mut p1 = GeometryRegistry::global().vertex_point(&v1);
        let mut p2 = GeometryRegistry::global().vertex_point(&v2);
        if crate::topo_tools_full::is_same(&v1, &v2) {
            let raw1 = self.resolve_shape(start_ref)?;
            let raw2 = self.resolve_shape(end_ref)?;
            p1 = GeometryRegistry::global().vertex_point(&raw1);
            p2 = GeometryRegistry::global().vertex_point(&raw2);
        }
        let (curve, first, last) = edge_from_curve3d(curve, &p1, &p2);
        // `MakeFromCurve3D` (`TranslateEdge.cxx:452-455`): distance at projected
        // params after `UpdateParam3d` (and any displaced-Line shift).
        let temp1 = curve.d0(first).distance(&p1);
        let temp2 = curve.d0(last).distance(&p2);
        let mut e = self.b.make_edge(curve, first, last);
        // `BRep_Builder` / `BRepLib_MakeEdge`: first vertex FORWARD, last REVERSED
        // (`TopoDS_Builder::Add`). Needed so `TopExp::LastVertex` exists on a
        // closed EDGE_CURVE (same TVertex stored twice) and
        // `ShapeAnalysis_Edge::FirstVertex` works on REVERSED seam uses.
        let v1 = Vertex(v1);
        let v2 = Vertex(v2);
        self.b.add_edge_vertices(&mut e, &v1, &v2);
        // `MakeFromCurve3D` (`TranslateEdge.cxx:478-479`): grow vertex tolerance
        // to cover Project residual (`BRep_Builder::UpdateVertex` only raises).
        // Enabled after ExtPC Project temps on Shape f7 seams are ~1e-7..1e-13
        // (well under preci=1e-3); invent sampler residuals previously densified.
        let t1 = crate::brep_tool::BRepTool::vertex_tolerance(&v1).max(1.000001 * temp1);
        let t2 = crate::brep_tool::BRepTool::vertex_tolerance(&v2).max(1.000001 * temp2);
        v1.set_tolerance(t1);
        v2.set_tolerance(t2);
        // Remember the edge's curve entity for the face-level pcurve association
        // (a SURFACE_CURVE's pcurve is matched to the face's surface by ref).
        self.edge_curve_ref
            .borrow_mut()
            .insert(Arc::as_ptr(&e.0.tshape) as usize, curve_ref);
        Ok(e.0)
    }

    pub(in crate::step) fn resolve_oriented_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let edge_ref = parse_ref(&rec.args[3]).ok_or("ORIENTED_EDGE: bad edge ref")?;
        let mut s = self.resolve_shape(edge_ref)?;
        // `StepToTopoDS_TranslateEdgeLoop.cxx:546-556`: the occurrence
        // orientation is the ORIENTED_EDGE boolean combined with the
        // EDGE_CURVE's `same_sense`:
        //   FORWARD  when (Orientation && SameSense) || (!Orientation && !SameSense)
        //   REVERSED otherwise.
        // The edge's stored forward direction already follows `same_sense`
        // (`resolve_edge`), so this is the XNOR of the two booleans. Without it
        // an `EDGE_CURVE(...,.F.)` occurrence comes back backwards.
        let ori = parse_logical(rec.args.get(4).map(String::as_str), true);
        let ec_same_sense = self.edge_curve_same_sense(edge_ref);
        s.set_orientation(if ori == ec_same_sense {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        Ok(s)
    }

    /// `same_sense` of the `EDGE_CURVE` an ORIENTED_EDGE references.
    /// `StepToTopoDS_TranslateEdgeLoop.cxx:306-308` follows a nested
    /// ORIENTED_EDGE down to its EDGE_CURVE (bug #29979). A chain that does not
    /// end at an EDGE_CURVE defaults to `.T.`.
    pub(super) fn edge_curve_same_sense(&self, mut id: usize) -> bool {
        for _ in 0..8 {
            let Ok(rec) = self.record(id) else {
                return true;
            };
            match rec.type_name.as_str() {
                "EDGE_CURVE" => {
                    return parse_logical(rec.args.get(4).map(String::as_str), true);
                }
                "ORIENTED_EDGE" => match parse_ref(&rec.args[3]) {
                    Some(next) => id = next,
                    None => return true,
                },
                _ => return true,
            }
        }
        true
    }

}
