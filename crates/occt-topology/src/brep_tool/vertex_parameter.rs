use std::sync::Arc;

use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::is_same;

use super::BRepTool;

impl BRepTool {
    /// `BRep_Tool::Parameter(V, E)` (`BRep_Tool.cxx:1340-1397`) ->
    /// `Parameter(V, E, theParam)` (`cxx:1334-1505`).
    ///
    /// The returned parameter is the edge's `First` or `Last` **chosen by the
    /// orientation the vertex carries inside the edge** (`cxx:1360-1383` walks
    /// the edge's vertex children and matches with `IsSame`), flipped when the
    /// edge handle itself is `REVERSED` (`cxx:1368`, `:1390`, `:1396`). The
    /// caller's own vertex orientation only matters for a degenerated edge that
    /// stores no vertex children (`cxx:1354-1357`).
    ///
    /// `None` is the `false` of [`Self::parameter_of_vertex_checked`], i.e. the
    /// `Standard_NoSuchObject` of the `Parameter(V, E)` wrapper (`cxx:1496-1503`).
    /// OCCT throws there; this entry point keeps the port's historical answer,
    /// the edge's `First` parameter - `parameter_of_vertex` is called from ~40
    /// places that predate the representation list, and every one of them
    /// expects a number.
    pub fn parameter_of_vertex(edge: &TopoShape, vertex: &TopoShape) -> f64 {
        match Self::parameter_of_vertex_checked(edge, vertex) {
            Some(p) => p,
            None => Self::edge_parameters(&Edge(edge.clone())).0,
        }
    }

    /// `BRep_Tool::Parameter(V, E, theParam)` (`BRep_Tool.cxx:1334-1505`):
    /// `Some` is the `true` return, `None` the `false` one.
    ///
    /// The four arms are the OCCT ones verbatim: `FORWARD` reads
    /// `(rev) ? l : f` (`cxx:1387-1392`), `REVERSED` reads `(rev) ? f : l`
    /// (`cxx:1393-1399`), and anything else [`Self::point_parameter_on_curve`] /
    /// [`Self::point_parameter_on_pcurve`] resolve through the vertex's point
    /// representations. The 3D-curve / p-curve choice is OCCT's: the p-curve
    /// arm runs only when `Curve(E)` is null **and** the edge is not degenerated
    /// (`cxx:1451`).
    pub fn parameter_of_vertex_checked(edge: &TopoShape, vertex: &TopoShape) -> Option<f64> {
        let key = GeometryRegistry::shape_key(vertex);
        let (first, last) = Self::edge_parameters(&Edge(edge.clone()));
        let mut found: Option<Orientation> = None;
        let mut rev = false;
        let mut has_vertex_children = false;
        if let Ok(ts) = edge.tshape.read() {
            if let Some(core) = ts.edge_core() {
                if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                    return Some(*p);
                }
            }
            for child in &ts.children {
                if child.shape_type() != ShapeType::Vertex {
                    continue;
                }
                has_vertex_children = true;
                if !is_same(child, vertex) {
                    continue;
                }
                match found {
                    None => found = Some(child.orientation()),
                    Some(cur) => {
                        rev = edge.orientation() == Orientation::Reversed;
                        if child.orientation() == vertex.orientation() {
                            found = Some(cur);
                        }
                    }
                }
            }
        }
        let orient = match found {
            Some(o) => o,
            None if !has_vertex_children && Self::is_degenerated(&Edge(edge.clone())) => {
                vertex.orientation()
            }
            None => Orientation::Internal,
        };
        match orient {
            Orientation::Forward => Some(if rev { last } else { first }),
            Orientation::Reversed => Some(if rev { first } else { last }),
            _ => {
                // `cxx:1437-1492`: the INTERNAL / EXTERNAL arm.
                let curve = Self::edge_curve(&Edge(edge.clone()));
                if curve.is_some() || Self::is_degenerated(&Edge(edge.clone())) {
                    let c = curve?;
                    Self::point_parameter_on_curve(vertex, &c, first, last)
                } else {
                    // `cxx:1458-1492`: no 3D curve, try the first p-curve.
                    let (surface, pc, _, _) = GeometryRegistry::global()
                        .edge_pcurve_reps(edge)
                        .into_iter()
                        .next()?;
                    Self::point_parameter_on_pcurve(vertex, &pc, &surface)
                }
            }
        }
    }

    /// `cxx:1425-1457` / `cxx:1615-1638`: the `IsPointOnCurve(C, L)` walk over
    /// the vertex's point representations, plus the closed-curve guard.
    ///
    /// `UNPORTED`: `TopLoc_Location`. OCCT compares `(myCurve == C) && (Location() == L)`
    /// and transforms `C->Value(f)` / `C->Value(l)` by `L`; the port stores no
    /// location on the representations and `Curve(E)` carries none, so the
    /// predicate is the curve handle identity and the two evaluations stay in
    /// the curve's own frame.
    fn point_parameter_on_curve(
        vertex: &TopoShape,
        curve: &Arc<dyn Curve>,
        first: f64,
        last: f64,
    ) -> Option<f64> {
        let reg = GeometryRegistry::global();
        let point = reg.vertex_point(vertex);
        let tol = reg.vertex_tolerance(vertex);
        let points = match vertex.tshape.read() {
            Ok(ts) => ts.vertex_core().map(|c| c.points.clone()).unwrap_or_default(),
            Err(_) => return None,
        };
        for pr in &points {
            if !pr.is_point_on_curve_of(curve) {
                continue;
            }
            let p = pr.parameter();
            let mut res = p;
            // Closed curves RLE 16 june 94.
            if occt_core::precision::Precision::is_negative_infinite(first) {
                return Some(res);
            }
            if occt_core::precision::Precision::is_positive_infinite(last) {
                return Some(res);
            }
            let pf = curve.d0(first);
            let pl = curve.d0(last);
            if pf.distance(&pl) < tol && pf.distance(&point) < tol {
                res = if vertex.orientation() == Orientation::Forward { first } else { last };
            }
            return Some(res);
        }
        None
    }

    /// `cxx:1458-1492`: `IsPointOnCurveOnSurface(PC, S, L)` plus the
    /// `PC->IsClosed()` guard, which pins a parameter that sits on either end of
    /// a closed p-curve to the end the vertex orientation names.
    ///
    /// `Geom2d_Curve::IsClosed()` (`Geom2d_Curve.cxx`) is
    /// `Value(FirstParameter()).Distance(Value(LastParameter())) <= Confusion`;
    /// the exact `p == FirstParameter() || p == LastParameter()` test after it is
    /// OCCT's own.
    fn point_parameter_on_pcurve(
        vertex: &TopoShape,
        pc: &Arc<dyn occt_geom2d::curve::Curve2d>,
        surface: &Arc<dyn Surface>,
    ) -> Option<f64> {
        let points = match vertex.tshape.read() {
            Ok(ts) => ts.vertex_core().map(|c| c.points.clone()).unwrap_or_default(),
            Err(_) => return None,
        };
        for pr in &points {
            if !pr.is_point_on_curve_on_surface_of(pc, surface) {
                continue;
            }
            let mut p = pr.parameter();
            let pc_first = pc.first_parameter();
            let pc_last = pc.last_parameter();
            let closed =
                pc.d0(pc_first).distance(&pc.d0(pc_last)) <= occt_core::precision::CONFUSION;
            if closed && (p == pc_first || p == pc_last) {
                p = if vertex.orientation() == Orientation::Forward { pc_first } else { pc_last };
            }
            return Some(p);
        }
        None
    }

    /// `BRep_Tool::Parameter(V, E, F)` (`BRep_Tool.cxx:1519-1524`) ->
    /// `Parameter(V, E, S, L)` (`cxx:1532-1620`): the parameter of `v` on the
    /// p-curve of `e` on `f`.
    ///
    /// The vertex-in-edge search and the `rev` latch are the `BRep_Tool::
    /// Parameter(V, E)` ones (`cxx:1537-1566`); the difference is the range the
    /// FORWARD / REVERSED arms read — `BRep_Tool::Range(E, S, L, f, l)`
    /// (`cxx:1572-1582`), i.e. the **p-curve** range, not the 3D one.
    ///
    /// OCCT's `Standard_NoSuchObject` (`cxx:1654`) surfaces as the `None` of
    /// [`Self::parameter_of_vertex_on_face_checked`]; this entry point keeps the
    /// port's historical answer, [`Self::parameter_of_vertex`].
    pub fn parameter_of_vertex_on_face(vertex: &TopoShape, edge: &Edge, face: &Face) -> f64 {
        match Self::parameter_of_vertex_on_face_checked(vertex, edge, face) {
            Some(p) => p,
            None => Self::parameter_of_vertex(&edge.0, vertex),
        }
    }

    /// `BRep_Tool::Parameter(V, E, S, L)` (`BRep_Tool.cxx:1532-1654`):
    /// `Some` is the `true` return, `None` the trailing `throw`
    /// (`cxx:1654`).
    ///
    /// The two representation walks are OCCT's and run in OCCT's order: the
    /// p-curve walk of the `INTERNAL` / `EXTERNAL` arm first
    /// (`cxx:1585-1597`, which answers the stored `Parameter()` with **no**
    /// closed-curve guard), then - after the arm falls through - the 3D-curve
    /// walk (`cxx:1615-1638`), gated on `!C.IsNull() || Degenerated(E)`
    /// (`cxx:1619`).
    pub fn parameter_of_vertex_on_face_checked(
        vertex: &TopoShape,
        edge: &Edge,
        face: &Face,
    ) -> Option<f64> {
        let (first, last) = Self::pcurve_range(edge, face);

        let key = GeometryRegistry::shape_key(vertex);
        let mut found: Option<Orientation> = None;
        let mut rev = false;
        if let Ok(ts) = edge.0.tshape.read() {
            if let Some(core) = ts.edge_core() {
                if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                    return Some(*p);
                }
            }
            for child in &ts.children {
                if child.shape_type() != ShapeType::Vertex || !is_same(child, vertex) {
                    continue;
                }
                match found {
                    None => found = Some(child.orientation()),
                    Some(cur) => {
                        rev = edge.0.orientation() == Orientation::Reversed;
                        if child.orientation() == vertex.orientation() {
                            found = Some(cur);
                        }
                    }
                }
            }
        }

        match found {
            Some(Orientation::Forward) => Some(if rev { last } else { first }),
            Some(Orientation::Reversed) => Some(if rev { first } else { last }),
            _ => {
                // `cxx:1583-1597`: `CurveOnSurface(E, S, L, f, l)` then the
                // `IsPointOnCurveOnSurface(PC, S, L)` walk. No closed-curve
                // guard here, unlike the `Parameter(V, E)` arm.
                let face_key = GeometryRegistry::shape_key(&face.0);
                if let (Some(pc), Some(surf)) = (
                    GeometryRegistry::global().edge_pcurve(edge, face_key),
                    Self::face_surface(face),
                ) {
                    let points = match vertex.tshape.read() {
                        Ok(ts) => {
                            ts.vertex_core().map(|c| c.points.clone()).unwrap_or_default()
                        }
                        Err(_) => Vec::new(),
                    };
                    for pr in &points {
                        if pr.is_point_on_curve_on_surface_of(&pc, &surf) {
                            return Some(pr.parameter());
                        }
                    }
                }
                // `cxx:1615-1638`: the 3D-curve walk, the same block the
                // `Parameter(V, E)` arm runs.
                let curve = Self::edge_curve(edge);
                if curve.is_some() || Self::is_degenerated(edge) {
                    let c = curve?;
                    Self::point_parameter_on_curve(vertex, &c, first, last)
                } else {
                    None
                }
            }
        }
    }
}
